!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! fluid_characteristic_step_mod
!
! Trace callback-defined advective characteristics through GORILLA's
! piecewise-linear tetrahedral geometry, using the shared event/step
! contract introduced by #80.
!
! A characteristic is a trajectory
!
!     xdot = u( x, t ; marker )
!
! supplied by the caller through a small data boundary (a velocity
! subroutine plus an optional marker state).  The tracer advances the
! position through the tetrahedral mesh, stopping at the earliest of
!
!   * the requested flight time,
!   * a tetrahedron-face crossing (CHAR_EVENT_FACE),
!   * a material (domain boundary / wall) intersection  (CHAR_EVENT_WALL),
!   * a caller-supplied stop event (CHAR_EVENT_COLLISION), e.g. a sampled
!     collision/reaction time.
!
! Deterministic face ownership and restart semantics across tetrahedra
! are preserved: a face is only an exit if the trajectory approaches it
! from the interior (inward normal speed negative); the crossing point is
! snapped onto the shared plane; and, when the caller restarts from a
! face event on the same plane (passing the exited tetrahedron/face back
! via restart_tet / restart_face), the point is owned by the entering
! neighbour tetrahedron rather than by the outgoing one.  A characteristic
! thus advances to the EARLIEST face, material-surface or caller-stop
! event; each call emits at most one event, so no missed, duplicated or
! zero-length face events are produced.
!
! Two stepping paths are provided per element:
!
!   * Fast path - when the caller certifies (exact_in_current = .true.)
!     that the characteristic is exactly integrable (uniform velocity) in
!     the current tetrahedron, the crossing time is computed analytically.
!     Constant-velocity trajectories therefore agree with the analytic
!     solution to roundoff.
!
!   * Controlled local integrator - otherwise an embedded Dormand-Prince
!     4(5) Runge-Kutta integrator with a relative/absolute error estimate
!     advances inside each tetrahedron, choosing steps so that the local
!     error meets the caller tolerance (defaults rtol = 1e-8, atol = 1e-10).
!     Affine velocity fields therefore meet the documented error target
!     (local error controlled to ~ rtol).
!
! The routine returns the segment tally needed by an external conservative
! transport application:
!
!     segment_integral(3) = int_0^t u( x(tau), tau ) dtau
!
! which, for the advective characteristic xdot = u, equals the net
! displacement x(t) - x(0) of the traced segment.  GORILLA returns this raw
! tally without performing mesh-field deposition itself.
!
! Geometry contract (matches the GORILLA mesh layout):
!   verts(3,nvert)          - vertex coordinates (grid coordinate system)
!   ind_knot(4,ntetr)       - vertex index (1-based) of the 4 corners per tetra
!   neighbour_tetr(4,ntetr) - tetrahedron index adjacent across each face;
!                             -1 : the face is a domain/wall boundary
!   neighbour_face(4,ntetr) - face index (1..4) on the neighbour side
!
! Face convention (consistent with tetra_physics_mod): face `iface` is the
! triangular facet opposite vertex `iface`; inward unit normals are used.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module fluid_characteristic_step_mod
!
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
!
    implicit none
!
    private
!
!------------------ Common characteristic result/event contract ---------------!
!
    ! Event kinds (shared vocabulary introduced by #80)
    integer, parameter, public :: CHAR_EVENT_NONE  = 0  ! flight finished within requested time
    integer, parameter, public :: CHAR_EVENT_FACE  = 1  ! tetrahedron face crossing
    integer, parameter, public :: CHAR_EVENT_WALL  = 2  ! material (wall) intersection
    integer, parameter, public :: CHAR_EVENT_COLLISION = 3 ! termination at caller stop event
!
    type, public :: characteristic_event_t
        integer  :: kind        = CHAR_EVENT_NONE ! event kind
        integer  :: tetrahedron = 0               ! tetrahedron index at the event
        integer  :: face        = 0               ! face index (1..4) for FACE/WALL events
        real(dp) :: time        = 0.0_dp          ! event time since start of this step
        real(dp) :: position(3) = 0.0_dp          ! event position (grid coordinates)
        real(dp) :: normal(3)   = 0.0_dp          ! unit normal; inward for FACE, outward for WALL
        real(dp) :: remaining   = 0.0_dp          ! remaining flight time after the event
    end type characteristic_event_t
!
    type, public :: characteristic_result_t
        real(dp) :: position(3) = 0.0_dp          ! final position (grid coordinates)
        real(dp) :: time        = 0.0_dp          ! total flight time used by this step
        type(characteristic_event_t) :: event     ! first event contract
        logical  :: finished     = .false.        ! reached requested time without material event
    end type characteristic_result_t
!
!------------------ Error codes -------------------------------------------------!
!
    integer, parameter, public :: CHAR_OK             = 0
    integer, parameter, public :: CHAR_ERR_INVALID    = 1  ! non-finite / non-positive / bad time
    integer, parameter, public :: CHAR_ERR_GRAZING    = 2  ! ambiguous grazing intersection
    integer, parameter, public :: CHAR_ERR_NOT_INSIDE = 3  ! starting point outside all tetrahedra
!
!------------------ Caller-supplied velocity field boundary --------------------!
!
    ! The caller supplies the local advective velocity.  The marker state is
    ! carried along the characteristic and may be read / updated by the
    ! callback (e.g. a carried density or colour field); it does not affect
    ! the geometric tracing except through the velocity it returns.
    abstract interface
        subroutine velocity_field_t(x, t, marker, v)
            import :: dp
            real(dp), intent(in)             :: x(3)
            real(dp), intent(in)             :: t
            real(dp), intent(inout)          :: marker(:)
            real(dp), intent(out)            :: v(3)
        end subroutine velocity_field_t
    end interface
!
    public :: trace_fluid_characteristic
!
!------------------ Tolerances --------------------------------------------------!
!
    real(dp), parameter :: eps_tol = 1.0e-10_dp     ! absolute geometric/time tolerance
    real(dp), parameter :: def_rtol = 1.0e-8_dp     ! default relative error (local integrator)
    real(dp), parameter :: def_atol = 1.0e-10_dp    ! default absolute error (local integrator)
!
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine trace_fluid_characteristic(x, velocity, marker, t_request, exact_in_current, &
                    & verts, ind_knot, neighbour_tetr, neighbour_face, ntetr, result, ierr, &
                    & segment_integral, t_stop, rtol, atol, restart_tet, restart_face)
!
        ! Trace one advective characteristic step of a caller-supplied velocity
        ! field through the tetrahedral mesh.
        !
        ! On entry:
        !   x                - position (grid coordinates), inout; on return the
        !                      final position reached by this step.
        !   velocity         - caller velocity subroutine  v = u(x,t;marker).
        !   marker           - optional marker state carried along the characteristic.
        !   t_request        - requested flight time (>= 0).
        !   exact_in_current - when .true. the caller certifies the characteristic
        !                      is exactly integrable (uniform velocity) in the
        !                      current element -> fast analytic path.  When .false.
        !                      a controlled local integrator with an error estimate
        !                      is used (optional rtol/atol tolerances).
        !   verts/ind_knot/neighbour_tetr/neighbour_face/ntetr - mesh contract.
        !   segment_integral - optional out; the segment tally  int_0^t u dtau
        !                      (= net displacement x(t)-x(0)).
        !   t_stop           - optional caller stop event time (>= 0) relative to the
        !                      start of this step; if reached first the flight
        !                      terminates there (CHAR_EVENT_COLLISION).
        !   rtol / atol      - optional local-integrator error tolerances.
        !   restart_tet / restart_face - optional determinism hints for restart
        !                      semantics: the (tetrahedron,face) just exited by a
        !                      previous face event.  When the starting point lies
        !                      exactly on that face the restart enters the
        !                      neighbour tetrahedron instead of re-exiting the same
        !                      face (deterministic face ownership); a zero-length
        !                      duplicate face event is thereby avoided.
        !
        ! On return `result` carries the common characteristic contract and `ierr`
        ! a CHAR_* error code.
        !
        implicit none
        !
        real(dp), dimension(3), intent(inout) :: x
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), intent(in)       :: t_request
        logical,  intent(in)       :: exact_in_current
        real(dp), dimension(3,*), intent(in)    :: verts
        integer,  dimension(4,*), intent(in)    :: ind_knot
        integer,  dimension(4,*), intent(in)    :: neighbour_tetr
        integer,  dimension(4,*), intent(in)    :: neighbour_face
        integer,  intent(in)                    :: ntetr
        type(characteristic_result_t), intent(out) :: result
        integer,  intent(out)                    :: ierr
        real(dp), dimension(3), intent(out), optional :: segment_integral
        real(dp), intent(in), optional           :: t_stop
        real(dp), intent(in), optional           :: rtol, atol
        integer,  intent(in), optional           :: restart_tet
        integer,  intent(in), optional           :: restart_face
!
        real(dp) :: t_remain, t_elapsed, t_stop_left, t_cross, t_elapsed_entry
        integer  :: ind_tetr, iface_out
        real(dp), dimension(3,4) :: anorm, anorm_t
        real(dp), dimension(4)   :: plane_d, plane_d_t
        real(dp), dimension(3)   :: v0, int_seg
        real(dp), dimension(1)   :: dummy_marker
        real(dp) :: rtol_local, atol_local
        logical  :: boole_stop
        real(dp) :: t_seg
        integer  :: if
        integer  :: n_tries
!
        ierr = CHAR_OK
        int_seg = 0.0_dp
        dummy_marker = 0.0_dp
!
        !-------------------------- State validation ---------------------------!
        if( .not. all( ieee_is_finite(x) ) ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
        if( .not. ieee_is_finite(t_request) .or. t_request.lt.0.0_dp ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
        if( present(marker) ) then
            if( .not. all( ieee_is_finite(marker) ) ) then
                ierr = CHAR_ERR_INVALID
                return
            endif
        endif
        if( present(t_stop) ) then
            boole_stop = .true.
            if( .not. ieee_is_finite(t_stop) .or. t_stop.lt.0.0_dp ) then
                ierr = CHAR_ERR_INVALID
                return
            endif
            t_stop_left = t_stop
        else
            boole_stop = .false.
            t_stop_left = huge(1.0_dp)
        endif
        rtol_local = def_rtol
        if( present(rtol) ) rtol_local = rtol
        atol_local = def_atol
        if( present(atol) ) atol_local = atol
!
        !------------------ Tolerance validation ---------------------------!
        ! Reject non-finite or invalid-sign tolerances before entering the
        ! controlled integrator: a NaN or negative value could make the
        ! step-size factor NaN or force pathological minimum-step iteration.
        if( .not. ieee_is_finite(rtol_local) .or. rtol_local.le.0.0_dp .or. &
            & .not. ieee_is_finite(atol_local) .or. atol_local.lt.0.0_dp ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
!
        ! Reset the contract
        result = characteristic_result_t()
        result%position(:) = x
        result%time        = 0.0_dp
!
        ! Zero requested time -> nothing to trace
        if( t_request.eq.0.0_dp ) then
            result%finished = .true.
            return
        endif
!
        !-------------------------- Locate starting tetrahedron -----------------!
        call locate_starting_tetrahedron(x, verts, ind_knot, ntetr, ind_tetr)
        if( ind_tetr.lt.1 ) then
            ierr = CHAR_ERR_NOT_INSIDE
            return
        endif
!
        ! Deterministic restart semantics: when resuming from a face event, a
        ! point lying exactly on the face just exited belongs to the entering
        ! (neighbour) tetrahedron, never to the outgoing one.  This keeps face
        ! ownership unique and prevents a zero-length re-crossing of the same
        ! face on the next step.
        if( present(restart_tet) .and. present(restart_face) ) then
            ! (Nested guards: Fortran does not guarantee short-circuit evaluation
            ! of .and., so an absent optional must never be referenced.)
            if( restart_tet.ge.1 .and. restart_face.ge.1 .and. restart_face.le.4 ) then
                if( point_on_face(x, restart_tet, restart_face, verts, ind_knot) ) then
                    if( neighbour_tetr(restart_face, restart_tet).ge.1 ) then
                        ind_tetr = neighbour_tetr(restart_face, restart_tet)
                    else
                        ! restart directly on a material surface boundary
                        result%finished    = .true.
                        result%event%kind  = CHAR_EVENT_WALL
                        result%event%tetrahedron = restart_tet
                        result%event%face  = restart_face
                        result%event%time  = 0.0_dp
                        result%event%position = x
                        call tet_geometry(restart_tet, verts, ind_knot, anorm_t, plane_d_t)
                        result%event%normal = -anorm_t(:,restart_face)
                        if( present(segment_integral) ) segment_integral = int_seg
                        return
                    endif
                endif
            endif
        endif
!
        !-------------------------- Trace through the mesh -----------------------!
        t_remain  = t_request
        t_elapsed = 0.0_dp
!
        tetra_loop: do n_tries = 1, ntetr + 1
!
            call tet_geometry(ind_tetr, verts, ind_knot, anorm, plane_d)
!
            if( exact_in_current ) then
                !================ Fast path: exactly integrable (uniform) =========!
                call eval_velocity(velocity, x, t_elapsed, marker, dummy_marker, v0)
                ndotv: block
                    real(dp) :: ndotv(4)
                    ndotv(1:4) = matmul( transpose(anorm), v0 )
                    iface_out = -1
                    t_cross   = huge(1.0_dp)
                    do if = 1,4
                        if( ndotv(if).lt.-eps_tol ) then
                            t_seg = -( dot_product(anorm(:,if), x) + plane_d(if) ) / ndotv(if)
                            if( t_seg.ge.-eps_tol .and. t_seg.lt.t_cross ) then
                                t_cross   = t_seg
                                iface_out = if
                            endif
                        endif
                    enddo
                end block ndotv
!
                ! Caller stop event may shorten the flight (exact in the fast path);
                ! a zero stop time is an immediate collision at the entry point.
                if( boole_stop .and. ( t_stop_left.eq.0.0_dp .or. &
                    & ( t_stop_left.lt.min( t_remain, t_cross ) .and. &
                    &   t_stop_left.ge.0.0_dp ) ) ) then
                    result%position(:) = x + v0*t_stop_left
                    result%time        = t_elapsed + t_stop_left
                    x                  = result%position
                    int_seg            = int_seg + v0*t_stop_left
                    result%event%kind  = CHAR_EVENT_COLLISION
                    result%event%tetrahedron = ind_tetr
                    result%event%face  = 0
                    result%event%time  = result%time
                    result%event%position = result%position
                    result%event%remaining = t_remain - t_stop_left
                    result%finished    = .true.
                    exit tetra_loop
                endif
!
                ! Flight ends inside this tetrahedron before any face
                if( iface_out.lt.1 .or. t_cross.ge.t_remain ) then
                    result%position(:) = x + v0*t_remain
                    result%time        = t_elapsed + t_remain
                    x                  = result%position
                    int_seg            = int_seg + v0*t_remain
                    result%finished    = .true.
                    result%event%kind  = CHAR_EVENT_NONE
                    result%event%tetrahedron = ind_tetr
                    exit tetra_loop
                endif
!
                !-------------------------- Cross face iface_out ---------------------!
                result%position(:) = x + v0*t_cross
                result%time        = t_elapsed + t_cross
                x                  = result%position
                int_seg            = int_seg + v0*t_cross
!
                result%event%kind        = CHAR_EVENT_FACE
                result%event%tetrahedron = ind_tetr
                result%event%face        = iface_out
                result%event%time        = result%time
                result%event%position    = x
                result%event%remaining   = t_remain - t_cross
                result%event%normal      = anorm(:,iface_out)   ! inward unit normal
!
                t_remain   = t_remain - t_cross
                t_elapsed  = t_elapsed + t_cross
                if( t_stop_left.gt.0.0_dp ) t_stop_left = t_stop_left - t_cross
!
                !------------------ Wall / material-surface event --------------------!
                ! The objective is to advance to the EARLIEST face, material-surface
                ! or caller-stop event and report it.  A tetrahedron-face crossing
                ! therefore stops the trace here (the caller restarts from the face);
                ! we must never silently hand over to a neighbour within one call,
                ! otherwise the first face event would be lost.
                if( neighbour_tetr(iface_out, ind_tetr).lt.1 ) then
                    ! No neighbour -> domain boundary: material (wall) event.
                    result%event%kind   = CHAR_EVENT_WALL
                    result%event%normal = -anorm(:,iface_out)   ! outward unit normal
                endif
                exit tetra_loop
!
            else
                !========= Controlled local integrator (RK45, error estimate) =======!
                t_elapsed_entry = t_elapsed
                call advance_in_tet(x, velocity, marker, dummy_marker, t_remain, t_elapsed, &
                     & anorm, plane_d, neighbour_tetr, ind_tetr, &
                     & t_stop_left, boole_stop, rtol_local, atol_local, &
                     & int_seg, result, ierr, iface_out)
                if( ierr.ne.CHAR_OK ) return
!
                select case( result%event%kind )
                case( CHAR_EVENT_FACE )
                    ! Advance to the earliest face event and report it; the caller
                    ! restarts from the face, so no hand-over continuation here.
                    exit tetra_loop
                case default
                    ! NONE (finished inside) / WALL / COLLISION -> done
                    exit tetra_loop
                end select
            endif
!
        enddo tetra_loop
!
        ! A consistent mesh must never run out of attempts (each face leads either
        ! to a wall or an unvisited tetrahedron).
        if( n_tries.gt.ntetr + 1 ) then
            ierr = CHAR_ERR_NOT_INSIDE
        endif
!
        if( present(segment_integral) ) segment_integral = int_seg
!
    end subroutine trace_fluid_characteristic
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine advance_in_tet(x, velocity, marker, dummy_marker, t_remain, t_elapsed, &
                    & anorm, plane_d, neighbour_tetr, ind_tetr, &
                    & t_stop_left, boole_stop, rtol, atol, &
                    & int_seg, result, ierr, iface_out)
!
        ! Controlled local-integrator leg inside one tetrahedron.
        !
        ! Advances the position with an embedded RK45 integrator (error estimate)
        ! until the earliest of
        !   * requested end time inside this tetrahedron (CHAR_EVENT_NONE),
        !   * a face crossing (CHAR_EVENT_FACE / CHAR_EVENT_WALL),
        !   * the caller stop event time (CHAR_EVENT_COLLISION).
        !
        ! The crossing point is snapped onto the shared plane; on a wall the
        ! outward normal is reported.  `t_stop_left` is the stop time remaining
        ! relative to the overall step start (monitored internally).
        !
        implicit none
        !
        real(dp), dimension(3), intent(inout) :: x
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), dimension(1), intent(inout) :: dummy_marker
        real(dp), intent(in)    :: t_remain, t_elapsed
        real(dp), dimension(3,4), intent(in) :: anorm
        real(dp), dimension(4), intent(in)   :: plane_d
        integer, dimension(4,*), intent(in)  :: neighbour_tetr
        integer,  intent(in)    :: ind_tetr
        real(dp), intent(in)    :: t_stop_left
        logical,  intent(in)    :: boole_stop
        real(dp), intent(in)    :: rtol, atol
        real(dp), dimension(3), intent(inout) :: int_seg
        type(characteristic_result_t), intent(out) :: result
        integer,  intent(out)   :: ierr
        integer,  intent(out)   :: iface_out
!
        real(dp) :: tau, h, hmin, avail, tol, err, fac, stop_due, tau_cross
        real(dp), dimension(3) :: xcur, xnew, x_enter, xcross
        real(dp), dimension(:), allocatable :: marker_new
        real(dp) :: d_prev, d_new
        integer  :: if, guard, if_tie
        logical  :: crossed, stop_limited
!
        ierr = CHAR_OK
        iface_out = -1
        x_enter = x
        xcur    = x
        tau     = 0.0_dp
        h       = min( t_remain, 1.0_dp )
        hmin    = 1.0e-12_dp
        guard   = 0
!
        ! Reset the result so callers can branch on the event kind
        result = characteristic_result_t()
!
        do
            stop_limited = .false.
            avail = t_remain - tau
            if( avail.le.eps_tol ) then
                ! requested end time reached inside this tetrahedron
                result%position(:) = xcur
                result%time        = t_elapsed + tau
                result%finished    = .true.
                result%event%kind  = CHAR_EVENT_NONE
                result%event%tetrahedron = ind_tetr
                x = xcur
                int_seg = int_seg + ( xcur - x_enter )
                return
            endif
!
            !------ caller stop event inside the next candidate step? -------!
            if( boole_stop ) then
                ! t_stop_left is the stop time remaining at entry (tau=0);
                ! the currently remaining portion is  stop_due = t_stop_left - tau
                stop_due = t_stop_left - tau
                if( stop_due.le.eps_tol .and. stop_due.gt.-eps_tol ) then
                    ! exactly at the stop time boundary -> report collision now
                    result%position(:) = xcur
                    result%time        = t_elapsed + tau
                    result%finished    = .true.
                    result%event%kind  = CHAR_EVENT_COLLISION
                    result%event%tetrahedron = ind_tetr
                    result%event%time  = result%time
                    result%event%position = xcur
                    result%event%remaining = t_remain - tau
                    x = xcur
                    int_seg = int_seg + ( xcur - x_enter )
                    return
                endif
                if( stop_due.lt.h .and. stop_due.gt.eps_tol ) then
                    ! stop time falls inside the coming step -> limit the candidate
                    ! step to the stop time and run the normal crossing logic, so a
                    ! face crossing within the interval is still reported as the
                    ! earliest event.
                    h = stop_due
                    stop_limited = .true.
                endif
            endif
!
            !--------------- one controlled RK45 step ------------------------!
            call rk45_step(xcur, marker, dummy_marker, t_elapsed + tau, h, velocity, xnew, err, &
                 & marker_out=marker_new)
            tol = atol + rtol*max( norm(xnew), norm(xcur) )
            if( err.gt.tol .and. h.gt.hmin .and. guard.lt.60 ) then
                fac   = 0.9_dp * ( tol / max( err, tiny(1.0_dp) ) )**0.2_dp
                h     = max( hmin, min( 5.0_dp*h, h*fac ) )
                guard = guard + 1
                cycle
            endif
            guard = 0
!
            !----------- detect a face crossing within this step ---------------!
            crossed = .false.
            do if = 1,4
                d_prev = dot_product( anorm(:,if), xcur ) + plane_d(if)
                d_new  = dot_product( anorm(:,if), xnew ) + plane_d(if)
                ! Detect a face crossing from the signed distance change.  A
                ! crossing is reported whenever the step ends clearly outside
                ! (d_new < -eps_tol) a face the leg did not start clearly outside
                ! (d_prev > -eps_tol), so trajectories that begin within eps_tol
                ! of a face and head outward are still caught.
                if( d_prev.gt.-eps_tol .and. d_new.lt.-eps_tol ) then
                    crossed = .true.
                    exit
                endif
            enddo
!
            if( crossed ) then
                call find_crossing(xcur, velocity, marker, dummy_marker, &
                     & t_elapsed + tau, h, anorm, plane_d, xnew, &
                     & iface_out, tau_cross, xcross, marker_out=marker_new)
                if( iface_out.ge.1 ) then
                    ! commit the carried state at the accepted crossing time
                    if( present(marker) ) marker = marker_new
                    ! snap onto the shared plane (deterministic ownership)
                    call project_onto_plane( xcross, anorm(:,iface_out), plane_d(iface_out), xcur )
                    int_seg = int_seg + ( xcur - x_enter )
                    tau     = tau + tau_cross
                    result%position(:) = xcur
                    result%time        = t_elapsed + tau
                    result%event%kind  = CHAR_EVENT_FACE
                    result%event%tetrahedron = ind_tetr
                    result%event%face        = iface_out
                    result%event%time        = result%time
                    result%event%position    = xcur
                    result%event%remaining   = max( 0.0_dp, t_remain - tau )
                    result%event%normal      = anorm(:,iface_out)   ! inward unit normal
                    x = xcur
                    if( neighbour_tetr(iface_out, ind_tetr).lt.1 ) then
                        ! material (wall) intersection
                        result%event%kind   = CHAR_EVENT_WALL
                        result%event%normal = -anorm(:,iface_out)   ! outward unit normal
                    endif
                    return
                endif
                ! no clean single-face crossing -> grazing; fall through
            endif
!
            !--- stop-limited leg ending exactly on a face plane (tie) ---------!
            ! When the caller stop time coincides with a face crossing, the
            ! fast path reports the face as the earliest event.  Match that
            ! tie-break here: a stop-limited step that ends on a plane the leg
            ! began strictly inside of is reported as that face crossing.
            if( .not. crossed .and. stop_limited ) then
                if_tie = 0
                do if = 1,4
                    d_prev = dot_product( anorm(:,if), xcur ) + plane_d(if)
                    d_new  = dot_product( anorm(:,if), xnew ) + plane_d(if)
                    if( d_prev.gt.eps_tol .and. abs( d_new ).le.eps_tol ) then
                        if_tie = if
                        exit
                    endif
                enddo
                if( if_tie.ge.1 ) then
                    if( present(marker) ) marker = marker_new
                    call project_onto_plane( xnew, anorm(:,if_tie), plane_d(if_tie), xcur )
                    int_seg = int_seg + ( xcur - x_enter )
                    tau     = tau + h
                    result%position(:) = xcur
                    result%time        = t_elapsed + tau
                    result%event%kind  = CHAR_EVENT_FACE
                    result%event%tetrahedron = ind_tetr
                    result%event%face        = if_tie
                    result%event%time        = result%time
                    result%event%position    = xcur
                    result%event%remaining   = max( 0.0_dp, t_remain - tau )
                    result%event%normal      = anorm(:,if_tie)   ! inward unit normal
                    x = xcur
                    if( neighbour_tetr(if_tie, ind_tetr).lt.1 ) then
                        result%event%kind   = CHAR_EVENT_WALL
                        result%event%normal = -anorm(:,if_tie)   ! outward unit normal
                    endif
                    return
                endif
            endif
!
            !------------------- accept the step -----------------------------!
            if( present(marker) ) marker = marker_new   ! commit accepted update
            xcur = xnew
            tau  = tau + h
            if( err.gt.0.0_dp ) then
                fac = 0.9_dp * ( tol / err )**0.2_dp
                h   = max( hmin, min( 5.0_dp*h, h*fac ) )
            endif
            h = min( h, t_remain - tau )
            if( h.le.hmin ) cycle
        enddo
!
    end subroutine advance_in_tet
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine find_crossing(xcur, velocity, marker, dummy_marker, t0, h, &
                    & anorm, plane_d, xnew, iface, tcross, xcross, marker_out)
!
        ! Locate the first time in (t0, t0+h] at which the trajectory crosses a
        ! tetrahedron face, and return the crossed face and crossing position.
        ! The trajectory is re-evaluated with the same local RK45 step so the
        ! crossing point is consistent with the accepted segment.  Only faces
        ! entered from the interior or from within eps_tol of the face
        ! (d(t0)>-eps_tol, d(t0+h)<-eps_tol) are considered, which keeps the
        ! ownership deterministic.
        ! `marker` is treated as committed state (preserved across evaluation);
        ! the carried state at the crossing time is returned through the
        ! optional `marker_out` for the caller to commit.
        !
        implicit none
        !
        real(dp), dimension(3), intent(in) :: xcur
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), dimension(1), intent(inout) :: dummy_marker
        real(dp), intent(in) :: t0, h
        real(dp), dimension(3,4), intent(in) :: anorm
        real(dp), dimension(4), intent(in)   :: plane_d
        real(dp), dimension(3), intent(in)   :: xnew
        integer, intent(out) :: iface
        real(dp), intent(out) :: tcross
        real(dp), dimension(3), intent(out) :: xcross
        real(dp), dimension(:), allocatable, intent(out), optional :: marker_out
!
        real(dp) :: lo, hi, mid, dlo, dhi, dmid
        real(dp), dimension(3) :: xmid, xsnap
        real(dp) :: errv
        integer :: f, iter
!
        iface   = -1
        tcross  = huge(1.0_dp)
        xcross  = xcur
!
        do f = 1,4
            dlo = dot_product( anorm(:,f), xcur ) + plane_d(f)
            dhi = dot_product( anorm(:,f), xnew ) + plane_d(f)
            if( .not.( dlo.gt.-eps_tol .and. dhi.lt.-eps_tol ) ) cycle
            ! refine with bisection on the plane distance
            lo = 0.0_dp
            hi = h
            do iter = 1, 70
                mid = 0.5_dp*( lo + hi )
                call rk45_step(xcur, marker, dummy_marker, t0, mid, velocity, xmid, errv)
                dmid = dot_product( anorm(:,f), xmid ) + plane_d(f)
                if( dmid.gt.0.0_dp ) then
                    lo = mid
                else
                    hi = mid
                endif
            enddo
            mid = 0.5_dp*( lo + hi )
            if( mid.lt.tcross ) then
                iface  = f
                tcross = mid
                call rk45_step(xcur, marker, dummy_marker, t0, tcross, velocity, xcross, errv, &
                     & marker_out=marker_out)
            endif
        enddo
!
        if( iface.ge.1 ) then
            ! snap the crossing point onto the plane for determinism
            call project_onto_plane( xcross, anorm(:,iface), plane_d(iface), xsnap )
            xcross = xsnap
        endif
!
    end subroutine find_crossing
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine rk45_step(x, marker, dummy_marker, t, h, velocity, xnew, err, marker_out)
!
        ! One embedded Dormand-Prince 4(5) step:
        !   5th-order solution  xnew  with embedded 4th-order error estimate  err
        !   (Euclidean norm of the difference of the two orders).
        ! The velocity callback may update the carried marker; the marker value
        ! at the end of the step reflects the carried state at the endpoint.
        ! The caller's `marker` is treated as committed state: it is preserved
        ! across evaluation (save/restore) so rejected trials, stage evaluations
        ! and crossing bisections can never corrupt it.  The candidate carried
        ! state at the endpoint is returned through the optional `marker_out`;
        ! the caller commits it to the caller marker only when a step is
        ! actually accepted.
        !
        implicit none
        !
        real(dp), dimension(3), intent(in)  :: x
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), dimension(1), intent(inout) :: dummy_marker
        real(dp), intent(in)  :: t, h
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(3), intent(out) :: xnew
        real(dp), intent(out) :: err
        real(dp), dimension(:), allocatable, intent(out), optional :: marker_out
!
        real(dp), dimension(3) :: k1, k2, k3, k4, k5, k6, k7, xish, e
        real(dp), dimension(:), allocatable :: marker_save
!
        ! Preserve committed caller state before the stage evaluations mutate
        ! the working marker copy across the RK stages.
        if( present(marker) ) then
            allocate( marker_save(size(marker)) )
            marker_save = marker
        endif
!
        call k_eval(velocity, x, t, marker, dummy_marker, k1)
        xish = x + h*( 1.0_dp/5.0_dp )*k1
        call k_eval(velocity, xish, t + (1.0_dp/5.0_dp)*h, marker, dummy_marker, k2)
        xish = x + h*( (3.0_dp/40.0_dp)*k1 + (9.0_dp/40.0_dp)*k2 )
        call k_eval(velocity, xish, t + (3.0_dp/10.0_dp)*h, marker, dummy_marker, k3)
        xish = x + h*( (44.0_dp/45.0_dp)*k1 - (56.0_dp/15.0_dp)*k2 + (32.0_dp/9.0_dp)*k3 )
        call k_eval(velocity, xish, t + (4.0_dp/5.0_dp)*h, marker, dummy_marker, k4)
        xish = x + h*( (19372.0_dp/6561.0_dp)*k1 - (25360.0_dp/2187.0_dp)*k2 &
                     & + (64448.0_dp/6561.0_dp)*k3 - (212.0_dp/729.0_dp)*k4 )
        call k_eval(velocity, xish, t + (8.0_dp/9.0_dp)*h, marker, dummy_marker, k5)
        xish = x + h*( (9017.0_dp/3168.0_dp)*k1 - (355.0_dp/33.0_dp)*k2 &
                     & + (46732.0_dp/5247.0_dp)*k3 + (49.0_dp/176.0_dp)*k4 &
                     & - (5103.0_dp/18656.0_dp)*k5 )
        call k_eval(velocity, xish, t + h, marker, dummy_marker, k6)
        xish = x + h*( (35.0_dp/384.0_dp)*k1 + (500.0_dp/1113.0_dp)*k3 &
                     & + (125.0_dp/192.0_dp)*k4 - (2187.0_dp/6784.0_dp)*k5 + (11.0_dp/84.0_dp)*k6 )
        call k_eval(velocity, xish, t + h, marker, dummy_marker, k7)
!
        xnew = x + h*( (35.0_dp/384.0_dp)*k1 + (500.0_dp/1113.0_dp)*k3 &
                     & + (125.0_dp/192.0_dp)*k4 - (2187.0_dp/6784.0_dp)*k5 + (11.0_dp/84.0_dp)*k6 )
!
        e = h*( (5179.0_dp/57600.0_dp - 35.0_dp/384.0_dp)*k1 &
            & + (7571.0_dp/16695.0_dp - 500.0_dp/1113.0_dp)*k3 &
            & + (393.0_dp/640.0_dp - 125.0_dp/192.0_dp)*k4 &
            & + ( -92097.0_dp/339200.0_dp + 2187.0_dp/6784.0_dp )*k5 &
            & + (187.0_dp/2100.0_dp - 11.0_dp/84.0_dp)*k6 &
            & - (1.0_dp/40.0_dp)*k7 )
        err = sqrt( dot_product(e,e) )
!
        ! Report the carried state at the step endpoint and restore the
        ! committed caller marker so evaluations never leak into caller state.
        if( present(marker_out) ) then
            if( present(marker) ) marker_out = marker
        endif
        if( present(marker) ) marker = marker_save
        if( allocated(marker_save) ) deallocate( marker_save )
!
    end subroutine rk45_step
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine k_eval(velocity, x, t, marker, dummy_marker, k)
        implicit none
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(3), intent(in)  :: x
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), dimension(1), intent(inout) :: dummy_marker
        real(dp), dimension(3), intent(out) :: k
        if( present(marker) ) then
            call velocity(x, t, marker, k)
        else
            call velocity(x, t, dummy_marker, k)
        endif
    end subroutine k_eval
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine eval_velocity(velocity, x, t, marker, dummy_marker, v)
        implicit none
        procedure(velocity_field_t) :: velocity
        real(dp), dimension(3), intent(in)  :: x
        real(dp), intent(in) :: t
        real(dp), dimension(:), intent(inout), optional :: marker
        real(dp), dimension(1), intent(inout) :: dummy_marker
        real(dp), dimension(3), intent(out) :: v
        if( present(marker) ) then
            call velocity(x, t, marker, v)
        else
            call velocity(x, t, dummy_marker, v)
        endif
    end subroutine eval_velocity
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    real(dp) function norm(v)
        implicit none
        real(dp), intent(in) :: v(3)
        norm = sqrt( dot_product(v,v) )
    end function norm
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine project_onto_plane(p, nrm, plane_d, pout)
        implicit none
        real(dp), dimension(3), intent(in)  :: p, nrm
        real(dp), intent(in) :: plane_d
        real(dp), dimension(3), intent(out) :: pout
        pout = p - nrm*( dot_product(nrm,p) + plane_d )
    end subroutine project_onto_plane
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    logical function point_on_face(x, itet, iface, verts, ind_knot) result(on_face)
!
        ! True if x lies (to within eps_tol) on the plane of face `iface` of
        ! tetrahedron `itet`.  Used to enforce deterministic restart ownership.
!
        implicit none
        !
        real(dp), dimension(3), intent(in)   :: x
        integer,  intent(in)                 :: itet, iface
        real(dp), dimension(3,*), intent(in) :: verts
        integer,  dimension(4,*), intent(in) :: ind_knot
        !
        real(dp), dimension(3,4) :: anorm
        real(dp), dimension(4)   :: plane_d
        !
        call tet_geometry(itet, verts, ind_knot, anorm, plane_d)
        on_face = abs( dot_product(anorm(:,iface), x) + plane_d(iface) ).le.eps_tol
!
    end function point_on_face
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine locate_starting_tetrahedron(x, verts, ind_knot, ntetr, ind_tetr)
!
        ! Find the tetrahedron that contains the starting point (boundary
        ! planes included).
!
        implicit none
        !
        real(dp), dimension(3), intent(in)   :: x
        real(dp), dimension(3,*), intent(in) :: verts
        integer, dimension(4,*), intent(in)  :: ind_knot
        integer, intent(in)                  :: ntetr
        integer, intent(out)                 :: ind_tetr
!
        integer :: i
        real(dp), dimension(3,4) :: anorm
        real(dp), dimension(4)   :: plane_d
!
        ind_tetr = -1
        do i = 1, ntetr
            call tet_geometry(i, verts, ind_knot, anorm, plane_d)
            if( all( ( matmul( transpose(anorm), x ) + plane_d ).ge.-eps_tol ) ) then
                ind_tetr = i
                return
            endif
        enddo
!
    end subroutine locate_starting_tetrahedron
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine tet_geometry(ind_tetr, verts, ind_knot, anorm, plane_d)
!
        ! Compute inward unit normals anorm(3,4) of the tetrahedron faces
        ! (face i opposes vertex i, consistent with tetra_physics_mod) and the
        ! signed plane constants such that a point p lies on the inward side of
        ! face i when  anorm(:,i).p + plane_d(i) >= 0  (equality on the plane).
!
        implicit none
        !
        integer, intent(in)                  :: ind_tetr
        real(dp), dimension(3,*), intent(in) :: verts
        integer, dimension(4,*), intent(in)  :: ind_knot
        real(dp), dimension(3,4), intent(out):: anorm
        real(dp), dimension(4),   intent(out):: plane_d
!
        real(dp), dimension(3,4) :: p
        real(dp), dimension(3)   :: u, w, n, opp
        integer :: i, j, l
        integer, dimension(3) :: corner
        real(dp) :: nrm
!
        do i = 1,4
            p(:,i) = verts(:, ind_knot(i, ind_tetr))
        enddo
!
        do i = 1,4
            l = 0
            do j = 1,4
                if( j.eq.i ) cycle
                l = l + 1
                corner(l) = j
            enddo
            u = p(:,corner(2)) - p(:,corner(1))
            w = p(:,corner(3)) - p(:,corner(1))
            n(1) = u(2)*w(3) - u(3)*w(2)
            n(2) = u(3)*w(1) - u(1)*w(3)
            n(3) = u(1)*w(2) - u(2)*w(1)
            nrm = sqrt( dot_product(n,n) )
            if( nrm.lt.eps_tol ) then
                n   = 0.0_dp
                nrm = 1.0_dp
            endif
            anorm(:,i) = n / nrm
!
            opp = p(:,i) - p(:,corner(1))
            if( dot_product(anorm(:,i), opp).lt.0.0_dp ) then
                anorm(:,i) = -anorm(:,i)
            endif
            plane_d(i) = -dot_product( anorm(:,i), p(:,corner(1)) )
        enddo
!
    end subroutine tet_geometry
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module fluid_characteristic_step_mod
