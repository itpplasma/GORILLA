!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Common characteristic-step / event API.
!
! This module defines the result/event contract shared by trajectory
! characteristic tracers in GORILLA and provides the neutral free-flight
! tracer requested in issue #81.
!
! The neutral tracer advances a zero-charge state ballistically (constant
! velocity in the grid coordinate system) through the GORILLA tetrahedral
! mesh until the requested time is reached, the first material (wall)
! event, or an application supplied collision time shortens the flight.
! It does NOT apply the guiding-centre equations and implements none of
! the reaction physics; a sampled collision/reaction time merely terminates
! the flight early.
!
! Geometry contract (matches the GORILLA mesh layout):
!   verts(3,nvert)        - vertex coordinates in the grid coordinate system
!   ind_knot(4,ntetr)     - vertex index (1-based) of the 4 corners per tetrahedron
!   neighbour_tetr(4,ntetr) - tetrahedron index adjacent across each face;
!                             -1 : the face is a domain/wall boundary
!   neighbour_face(4,ntetr) - face index (1..4) on the neighbour side
!
! Face convention (consistent with tetra_physics_mod): face `iface` is the
! triangular facet opposite vertex `iface`; inward unit normals are used.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module characteristic_step_mod
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
    ! Event kinds
    integer, parameter, public :: CHAR_EVENT_NONE  = 0  ! flight finished within requested time
    integer, parameter, public :: CHAR_EVENT_FACE  = 1  ! tetrahedron face crossing
    integer, parameter, public :: CHAR_EVENT_WALL  = 2  ! material (wall) intersection
    integer, parameter, public :: CHAR_EVENT_COLLISION = 3 ! termination at sampled collision time
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
    integer, parameter, public :: CHAR_ERR_INVALID    = 1  ! non-finite / non-positive mass / bad time
    integer, parameter, public :: CHAR_ERR_GRAZING    = 2  ! ambiguous grazing intersection
    integer, parameter, public :: CHAR_ERR_NOT_INSIDE = 3  ! starting point outside all tetrahedra
!
    public :: trace_neutral_free_flight
!
    real(dp), parameter :: eps_tol = 1.0e-10_dp     ! absolute geometric/time tolerance
!
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine trace_neutral_free_flight(x, v, mass, t_request, verts, ind_knot, &
                    & neighbour_tetr, neighbour_face, ntetr, result, ierr, t_collision)
!
        ! Trace one free-flight characteristic step of a neutral (zero-charge)
        ! particle through the tetrahedral mesh.
        !
        ! The position advances ballistically  x(t) = x + v*t  in the grid
        ! coordinate system (the same coordinates the mesh vertices use).
        !
        ! On return, `result` carries the common characteristic contract:
        !   - final position and flown time
        !   - the first event: either a completed flight (CHAR_EVENT_NONE,
        !     finished = .true.), a face crossing (CHAR_EVENT_FACE), a wall
        !     intersection (CHAR_EVENT_WALL) or an application-shortened
        !     collision (CHAR_EVENT_COLLISION).
        !
        ! Optional argument t_collision: application sampled collision/reaction
        ! time (>= 0, relative to the start of this step). If it is reached
        ! before any face/wall crossing, the flight terminates there. Core does
        ! not implement the reaction itself.
        !
        implicit none
        !
        real(dp), dimension(3), intent(inout) :: x          ! position (grid coordinates)
        real(dp), dimension(3), intent(in)    :: v          ! constant velocity (grid coords / time)
        real(dp), intent(in)                  :: mass       ! particle mass (> 0, validation only)
        real(dp), intent(in)                  :: t_request  ! requested flight time (>= 0)
        real(dp), dimension(3,*), intent(in)  :: verts      ! vertex coordinates
        integer, dimension(4,*), intent(in)   :: ind_knot   ! vertex indices per tetrahedron
        integer, dimension(4,*), intent(in)   :: neighbour_tetr  ! neighbour across each face (-1 = wall)
        integer, dimension(4,*), intent(in)   :: neighbour_face  ! entry face on the neighbour side
        integer, intent(in)                   :: ntetr      ! number of tetrahedra
        type(characteristic_result_t), intent(out) :: result
        integer, intent(out)                  :: ierr       ! CHAR_OK or an error code
        real(dp), intent(in), optional        :: t_collision
!
        real(dp) :: t_remain, t_coll, t_seg, t_cross, t_elapsed
        integer  :: ind_tetr, iface, iface_out, n_tries
        real(dp), dimension(4)   :: ndotv
        real(dp), dimension(3,4) :: anorm      ! inward unit normals of current tetrahedron
        real(dp), dimension(4)   :: plane_d    ! signed plane constant of each face
        logical :: boole_collision
!
        ierr = CHAR_OK
!
        !-------------------------- State validation ---------------------------!
        if( mass.le.0.0_dp ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
        if( .not. ( all( ieee_is_finite(x) ) .and. all( ieee_is_finite(v) ) ) ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
        if( .not. ieee_is_finite(t_request) .or. t_request.lt.0.0_dp ) then
            ierr = CHAR_ERR_INVALID
            return
        endif
        if( present(t_collision) ) then
            boole_collision = .true.
            t_coll = t_collision
            if( .not. ieee_is_finite(t_coll) .or. t_coll.lt.0.0_dp ) then
                ierr = CHAR_ERR_INVALID
                return
            endif
        else
            boole_collision = .false.
            t_coll = huge(1.0_dp)
        endif
!
        ! Reset contract
        result = characteristic_result_t()
        result%position(:) = x
        result%time        = 0.0_dp
!
        ! Zero requested time or zero speed -> nothing to trace
        if( t_request.eq.0.0_dp .or. all( v.eq.0.0_dp ) ) then
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
        !-------------------------- Trace through the mesh -----------------------!
        t_remain = t_request
        t_elapsed = 0.0_dp
!
        tetra_loop: do n_tries = 1, ntetr + 1
!
            ! Current tetrahedron geometry
            call tet_geometry(ind_tetr, verts, ind_knot, anorm, plane_d)
!
            ! Outward normal speeds through every face
            ndotv(1:4) = matmul( transpose(anorm), v )
!
            ! Earliest crossing time over candidate (outward) exit faces
            iface_out = -1
            t_cross   = huge(1.0_dp)
            do iface = 1,4
                if( ndotv(iface).lt.-eps_tol ) then
                    t_seg = -( dot_product(anorm(:,iface), x) + plane_d(iface) ) / ndotv(iface)
                    if( t_seg.lt.t_cross ) then
                        t_cross   = t_seg
                        iface_out = iface
                    endif
                endif
            enddo
!
            if( iface_out.lt.1 ) then
                ! No face is approached: the straight line stays inside this
                ! convex tetrahedron -> flight finishes inside.
                result%position(:) = x + v*t_remain
                result%time        = t_elapsed + t_remain
                result%finished    = .true.
                result%event%kind  = CHAR_EVENT_NONE
                x = result%position
                return
            endif
!
            ! Grazing guard: the crossing must be unambiguous. Two candidate
            ! faces crossed at (numerically) the same time mean the path hits an
            ! edge/vertex -> reject as an ambiguous grazing intersection.
            do iface = 1,4
                if( iface.eq.iface_out ) cycle
                if( ndotv(iface).lt.-eps_tol ) then
                    t_seg = -( dot_product(anorm(:,iface), x) + plane_d(iface) ) / ndotv(iface)
                    if( abs( t_seg - t_cross ).lt. eps_tol ) then
                        ierr = CHAR_ERR_GRAZING
                        result%position(:) = x
                        result%time        = t_elapsed
                        result%event%kind  = CHAR_EVENT_FACE
                        result%event%tetrahedron = ind_tetr
                        result%event%time        = t_elapsed
                        result%event%position    = x
                        result%event%remaining   = t_remain
                        return
                    endif
                endif
            enddo
!
            ! Application collision time may shorten the flight. t_collision is
            ! relative to the start of the whole step, so subtract the time
            ! already elapsed in earlier segments before comparing against the
            ! remaining segment time.
            if( boole_collision .and. (t_coll - t_elapsed).lt.min(t_remain,t_cross) ) then
                t_seg = t_coll - t_elapsed
                result%position(:) = x + v*t_seg
                result%time        = t_elapsed + t_seg
                result%event%kind        = CHAR_EVENT_COLLISION
                result%event%tetrahedron = ind_tetr
                result%event%time        = result%time
                result%event%position    = result%position
                result%event%remaining   = t_remain - t_seg
                result%finished          = .true.
                x = result%position
                return
            endif
!
            ! Flight ends inside the current tetrahedron before any face
            if( t_cross.ge.t_remain ) then
                result%position(:) = x + v*t_remain
                result%time        = t_elapsed + t_remain
                result%finished    = .true.
                result%event%kind  = CHAR_EVENT_NONE
                x = result%position
                return
            endif
!
            !-------------------------- Cross face iface_out ---------------------!
            result%position(:) = x + v*t_cross
            result%time        = t_elapsed + t_cross
            x                  = result%position
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
!
            !------------------ Wall / handover to neighbour ---------------------!
            if( neighbour_tetr(iface_out, ind_tetr).lt.1 ) then
                ! No neighbour -> domain boundary: material (wall) event.
                result%event%kind   = CHAR_EVENT_WALL
                result%event%normal = -anorm(:,iface_out)   ! outward unit normal
                return
            endif
!
            iface     = neighbour_face(iface_out, ind_tetr)
            ind_tetr  = neighbour_tetr(iface_out, ind_tetr)
!
        enddo tetra_loop
!
        ! A consistent mesh must never get here (each face leads either to a
        ! wall or to an unvisited tetrahedron).
        ierr = CHAR_ERR_NOT_INSIDE
        result%position(:) = x
        result%time        = t_elapsed
!
    end subroutine trace_neutral_free_flight
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
        integer :: i, j
        real(dp), dimension(3,4) :: anorm
        real(dp), dimension(4)   :: plane_d
        real(dp) :: dist
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
        ! signed plane constants plane_d such that a point p lies on the inward
        ! side of face i when  anorm(:,i).p + plane_d(i) >= 0  (equality on the
        ! face plane).
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
            ! Face i corners: the three vertices distinct from vertex i
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
                ! Degenerate tetrahedron: treat as a unit-empty face; keep an
                ! arbitrary normal so downstream checks remain safe.
                n = 0.0_dp
                nrm = 1.0_dp
            endif
            anorm(:,i) = n / nrm
!
            ! Orient inward: the opposing vertex i must be on the positive side
            opp = p(:,i) - p(:,corner(1))
            if( dot_product(anorm(:,i), opp).lt.0.0_dp ) then
                anorm(:,i) = -anorm(:,i)
            endif
!
            plane_d(i) = -dot_product( anorm(:,i), p(:,corner(1)) )
        enddo
!
    end subroutine tet_geometry
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module characteristic_step_mod
