module test_fluid_characteristic_step_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use fluid_characteristic_step_mod
   use funit
   implicit none
!
! Tests for trace_fluid_characteristic: callback-defined advective
! characteristic tracing through a tetrahedral (GORILLA) mesh, following the
! event/step contract introduced by #80.
!
! Mesh layout (analytic):
!   Tetrahedron A: (0,0,0), (2,0,0), (0,2,0), (0,0,2)
!   Tetrahedron B: (2,0,0), (0,2,0), (0,0,2), (0,0,4)
!   Shared face  : plane x+y+z=2 (A face 1 <-> B face 4)
! All remaining faces are domain boundaries (walls).
!
! The velocity is provided by an external callback subroutine. Two callbacks
! are used here:
!   * u_const: uniform velocity (fast exact-integration path certified), and

   integer, parameter :: ntetr = 2
   integer, parameter :: nvert = 5
   real(dp) :: verts(3, nvert)
   integer  :: ind_knot(4, ntetr), nb_tetr(4, ntetr), nb_face(4, ntetr)
   real(dp) :: uconst(3)          ! constant velocity chosen per test
   real(dp), dimension(2) :: marker_state ! optional marker state

contains

   subroutine setUp()
      verts(:, 1) = [0.d0, 0.d0, 0.d0]
      verts(:, 2) = [2.d0, 0.d0, 0.d0]
      verts(:, 3) = [0.d0, 2.d0, 0.d0]
      verts(:, 4) = [0.d0, 0.d0, 2.d0]
      verts(:, 5) = [0.d0, 0.d0, 4.d0]
      ind_knot(:, 1) = [1, 2, 3, 4]     ! A
      ind_knot(:, 2) = [2, 3, 4, 5]     ! B
      nb_tetr = -1
      nb_face = -1
      nb_tetr(1, 1) = 2               ! A face 1 -> tetrahedron B
      nb_face(1, 1) = 4               !     ... via B face 4
      uconst = [0.d0, 0.d0, 0.d0]
      marker_state = [1.d0, 0.d0]
   end subroutine setUp

   !--- Caller velocity callbacks ------------------------------------------!
   subroutine u_const(x, t, marker, v)
      real(dp), intent(in)   :: x(3), t
      real(dp), intent(inout):: marker(:)
      real(dp), intent(out)  :: v(3)
      v = uconst
   end subroutine u_const

   ! Exponential-decay field  u = (0,0,1+2z)  used for an exact analytic
   ! convergence reference in one tetrahedron.
   subroutine u_expz(x, t, marker, v)
      real(dp), intent(in)   :: x(3), t
      real(dp), intent(inout):: marker(:)
      real(dp), intent(out)  :: v(3)
      v = [0.d0, 0.d0, 1.d0 + 2.d0 * x(3)]
   end subroutine u_expz

   !---------------- Constant velocity (fast exact path) -----------------------!
   @test
   subroutine test_const_velocity_inside_fast()
      real(dp) :: x(3), seg(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      call trace_fluid_characteristic(x, u_const, marker_state, 0.5d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, segment_integral=seg)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(0.5d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 0.7d0], res%position, tolerance=1.d-12)
      @assertEqual(CHAR_EVENT_NONE, res%event%kind)
      ! segment tallies for an external conservative transport application
      @assertEqual([0.d0, 0.d0, 0.5d0], seg, tolerance=1.d-12)
   end subroutine test_const_velocity_inside_fast

   @test
   subroutine test_const_velocity_inside_integrator()
      ! Same trajectory, but forced through the controlled local integrator
      ! (exact_in_current = .false.).
      real(dp) :: x(3), seg(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      call trace_fluid_characteristic(x, u_const, marker_state, 0.5d0, .false., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, segment_integral=seg)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(0.5d0, res%time, tolerance=1.d-10)
      @assertEqual([0.2d0, 0.2d0, 0.7d0], res%position, tolerance=1.d-9)
      @assertEqual(CHAR_EVENT_NONE, res%event%kind)
      @assertEqual([0.d0, 0.d0, 0.5d0], seg, tolerance=1.d-9)
   end subroutine test_const_velocity_inside_integrator

   !---------------- Affine velocity (integrator, convergence) -----------------!
   @test
   subroutine test_affine_velocity_analytic_agreement()
      ! u = (0,0,1+2z): dz/dt = 1+2z,  z(t) = ((1+2 z0) e^{2t} - 1)/2.
      ! With z0 = 0.2 and t = 0.2 the flight stays inside tetrahedron A,
      ! giving an exact analytic comparison for the controlled integrator.
      real(dp) :: x(3), z_analytic
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.5d0, 0.5d0, 0.2d0]
      z_analytic = 0.5d0 * (1.4d0 * exp(0.4d0) - 1.d0)
      call trace_fluid_characteristic(x, u_expz, marker_state, 0.2d0, .false., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, rtol=1.d-8, atol=1.d-10)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      ! meet the documented error target: < 1e-6 relative (default rtol=1e-8)
      @assertEqual(z_analytic, res%position(3), tolerance=1.d-7)
   end subroutine test_affine_velocity_analytic_agreement

   !---------------- Boundary (wall) termination ------------------------------!
   @test
   subroutine test_const_wall_termination_fast()
      ! Vertical constant rise starting inside tetrahedron B (above the shared
      ! face) straight into the B face-3 wall.  No internal face is crossed, so
      ! the wall is the earliest (and only) event.
      real(dp) :: x(3), seg(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 1.7d0]
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, segment_integral=seg)
      @assertEqual(CHAR_OK, ierr)
      @assertFalse(res%finished)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(3, res%event%face)
      @assertEqual(1.5d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 3.2d0], res%position, tolerance=1.d-12)
      @assertEqual(28.5d0, res%event%remaining, tolerance=1.d-12)
      @assertEqual([2.d0/3.d0, 2.d0/3.d0, 1.d0/3.d0], res%event%normal, tolerance=1.d-12)
      ! tallies along the whole constant segment equal the net displacement
      @assertEqual([0.d0, 0.d0, 1.5d0], seg, tolerance=1.d-12)
   end subroutine test_const_wall_termination_fast

   @test
   subroutine test_const_wall_termination_integrator()
      ! Same wall on the integrator path, matching analytic to roundoff on time
      ! and agreeing with the reproducible boundary identifier (tet 2 face 3).
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 1.7d0]
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .false., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertFalse(res%finished)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(3, res%event%face)
      @assertEqual(1.5d0, res%time, tolerance=1.d-9)
      @assertEqual([2.d0/3.d0, 2.d0/3.d0, 1.d0/3.d0], res%event%normal, tolerance=1.d-9)
   end subroutine test_const_wall_termination_integrator

   !--------------- Tetrahedron crossings with deterministic restart ----------!
   @test
   subroutine test_multi_tetra_no_duplicate_face_events_fast()
      ! A characteristic that crosses A -> B -> wall must report exactly one
      ! face event (A face 1) with no zero-length duplicate or missed crossing,
      ! then a unique wall event (B face 3).  The restart hints give
      ! deterministic face ownership on the shared plane x+y+z=2.
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr, rt, rf, n_face

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      n_face = 0
      rt = -1; rf = -1
      ! step 1: cross from A into B
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, &
           restart_tet=rt, restart_face=rf)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_FACE, res%event%kind)
      @assertEqual(1, res%event%tetrahedron)
      @assertEqual(1, res%event%face)
      @assertEqual(1.4d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 1.6d0], res%position, tolerance=1.d-12)
      @assertTrue(res%time > 1.d-12)  ! not a zero-length event
      n_face = n_face + 1
      rt = res%event%tetrahedron; rf = res%event%face
      ! step 2: continue from the face into B, hitting the wall
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, &
           restart_tet=rt, restart_face=rf)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(3, res%event%face)
      @assertEqual(3.0d0, res%time + 1.4d0, tolerance=1.d-12)  ! total time 3.0
      if( res%event%kind.eq.CHAR_EVENT_FACE ) n_face = n_face + 1
      @assertEqual(1, n_face)  ! exactly one face crossing, no duplicate
   end subroutine test_multi_tetra_no_duplicate_face_events_fast

   @test
   subroutine test_multi_tetra_no_duplicate_face_events_integrator()
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr, rt, rf, n_face

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      n_face = 0
      rt = -1; rf = -1
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .false., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, &
           restart_tet=rt, restart_face=rf)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_FACE, res%event%kind)
      @assertEqual(1, res%event%tetrahedron)
      @assertEqual(1, res%event%face)
      @assertEqual(1.4d0, res%time, tolerance=1.d-9)
      @assertTrue(res%time > 1.d-9)
      n_face = n_face + 1
      rt = res%event%tetrahedron; rf = res%event%face
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .false., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, &
           restart_tet=rt, restart_face=rf)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(3, res%event%face)
      @assertEqual(3.0d0, res%time + 1.4d0, tolerance=1.d-8)
      if( res%event%kind.eq.CHAR_EVENT_FACE ) n_face = n_face + 1
      @assertEqual(1, n_face)
   end subroutine test_multi_tetra_no_duplicate_face_events_integrator

   !---------------- Caller stop (collision) event -----------------------------!
   @test
   subroutine test_call_stop_event()
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      call trace_fluid_characteristic(x, u_const, marker_state, 30.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr, t_stop=0.5d0)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_COLLISION, res%event%kind)
      @assertTrue(res%finished)
      @assertEqual(0.5d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 0.7d0], res%position, tolerance=1.d-12)
   end subroutine test_call_stop_event

   !---------------- Invalid-input rejection (mesh-free) -----------------------!
   @test
   subroutine test_reject_negative_time()
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      call trace_fluid_characteristic(x, u_const, marker_state, -1.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)
   end subroutine test_reject_negative_time

   @test
   subroutine test_reject_nonfinite_position()
      use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
      real(dp) :: x(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      uconst = [0.d0, 0.d0, 1.d0]
      x = [0.2d0, 0.2d0, 0.2d0]
      x(1) = ieee_value(1.0_dp, ieee_quiet_nan)
      call trace_fluid_characteristic(x, u_const, marker_state, 1.d0, .true., &
           verts, ind_knot, nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)
   end subroutine test_reject_nonfinite_position

end module test_fluid_characteristic_step_mod
