module test_characteristic_step_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
   use characteristic_step_mod
   use funit
   implicit none
!
! Tests for trace_neutral_free_flight: neutral free-flight characteristic
! tracing through a tetrahedral (GORILLA) mesh. An analytical two-tetrahedron
! assembly replaces the produced mesh so that exact results are known.
!
! Mesh layout:
!   Tetrahedron A: (0,0,0), (2,0,0), (0,2,0), (0,0,2)
!   Tetrahedron B: (2,0,0), (0,2,0), (0,0,2), (0,0,4)
!   Shared face  : plane x+y+z=2 (A face 1 <-> B face 4)
! All remaining faces are domain boundaries (walls).
!
   integer, parameter :: ntetr = 2
   integer, parameter :: nvert = 5
   real(dp) :: verts(3, nvert)
   integer  :: ind_knot(4, ntetr), nb_tetr(4, ntetr), nb_face(4, ntetr)
   real(dp) :: mass

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
      mass = 1.d0
   end subroutine setUp

   !----- Constant velocity path through one tetrahedron --------------------!
   @test
   subroutine test_constant_velocity_inside()
      real(dp) :: x(3), v(3), t_req
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      t_req = 0.5d0
      call trace_neutral_free_flight(x, v, mass, t_req, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(0.5d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 0.7d0], res%position, tolerance=1.d-12)
   end subroutine test_constant_velocity_inside

   !----- Exact face crossing and wall event after neighbour handover -------!
   @test
   subroutine test_handover_then_wall()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 30.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertFalse(res%finished)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(3.0d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 3.2d0], res%position, tolerance=1.d-12)
      @assertEqual(30.d0 - 3.0d0, res%event%remaining, tolerance=1.d-12)
      ! outward wall normal of B face 3
      @assertEqual([2.d0/3.d0, 2.d0/3.d0, 1.d0/3.d0], res%event%normal, tolerance=1.d-12)
   end subroutine test_handover_then_wall

   !----- Exact wall facet for a pure wall hit -------------------------------!
   @test
   subroutine test_wall_intersection_metadata()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, -1.d0, 0.d0]
      call trace_neutral_free_flight(x, v, mass, 5.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(3, res%event%face)                 ! face y=0
      @assertEqual([0.d0, -1.d0, 0.d0], res%event%normal, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.d0, 0.2d0], res%position, tolerance=1.d-12)
      @assertEqual(0.2d0, res%time, tolerance=1.d-12)
      @assertEqual(4.8d0, res%event%remaining, tolerance=1.d-12)
   end subroutine test_wall_intersection_metadata

   !----- Reversibility of event-free free flight ----------------------------!
   @test
   subroutine test_reversibility()
      real(dp) :: x(3), v(3), t_req, t0
      type(characteristic_result_t) :: res
      integer :: ierr

      ! forward: ends exactly on the internal shared face (A->B), finished
      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      t_req = 1.4d0
      call trace_neutral_free_flight(x, v, mass, t_req, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      t0 = res%time

      ! reverse at the end point
      x = res%position
      v = [0.d0, 0.d0, -1.d0]
      call trace_neutral_free_flight(x, v, mass, t_req, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(t0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 0.2d0], res%position, tolerance=1.d-12)
   end subroutine test_reversibility

   !----- Optional collision time shortens, reaction physics not needed ------!
   @test
   subroutine test_collision_shortens()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 5.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr, 1.1d0)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(CHAR_EVENT_COLLISION, res%event%kind)
      @assertEqual(1.1d0, res%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 1.3d0], res%position, tolerance=1.d-12)
   end subroutine test_collision_shortens

   !----- Collision scheduled after the first face crossing fires in the
   !----- second tetrahedron, at the step-relative absolute time -----------!
   @test
   subroutine test_collision_after_face_crossing()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      ! x=0.2,z=0.2 moving +z; crosses A->B at t=1.4 (plane x+y+z=2).
      ! Collision at t=1.6 lies after the face, inside tetrahedron B.
      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 5.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr, 1.6d0)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual(CHAR_EVENT_COLLISION, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)
      @assertEqual(1.6d0, res%time, tolerance=1.d-12)
      @assertEqual(1.6d0, res%event%time, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 1.8d0], res%position, tolerance=1.d-12)
      @assertEqual([0.2d0, 0.2d0, 1.8d0], x, tolerance=1.d-12)   ! intent(inout) updated
      @assertEqual(5.d0 - 1.6d0, res%event%remaining, tolerance=1.d-12)
   end subroutine test_collision_after_face_crossing

   !----- Flight completing inside a tetrahedron updates intent(inout) x ---!
   @test
   subroutine test_x_updated_when_finished_inside()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 0.5d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertTrue(res%finished)
      @assertEqual([0.2d0, 0.2d0, 0.7d0], x, tolerance=1.d-12)
   end subroutine test_x_updated_when_finished_inside

   !----- Collision later than the geometry event is ignored ----------------!
   @test
   subroutine test_collision_later_ignored()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.2d0, 0.2d0, 0.2d0]
      v = [0.d0, 0.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 5.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr, 9.d0)
      @assertEqual(CHAR_OK, ierr)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(3.0d0, res%time, tolerance=1.d-12)
   end subroutine test_collision_later_ignored

   !----- Rejections: mass, non-finite state, bad time, outside --------------!
   @test
   subroutine test_rejections()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      ! zero mass
      x = [0.2d0, 0.2d0, 0.2d0]; v = [1.d0, 0.d0, 0.d0]
      call trace_neutral_free_flight(x, v, 0.d0, 1.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)

      ! non-finite velocity
      x = [0.2d0, 0.2d0, 0.2d0]
      v = [ieee_value(1.d0, ieee_positive_inf), 0.d0, 0.d0]
      call trace_neutral_free_flight(x, v, mass, 1.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)

      ! non-finite position
      x = [0.2d0, ieee_value(1.d0, ieee_positive_inf), 0.2d0]
      v = [1.d0, 0.d0, 0.d0]
      call trace_neutral_free_flight(x, v, mass, 1.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)

      ! negative time
      x = [0.2d0, 0.2d0, 0.2d0]; v = [1.d0, 0.d0, 0.d0]
      call trace_neutral_free_flight(x, v, mass, -1.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_INVALID, ierr)

      ! starting point outside the mesh
      x = [5.d0, 5.d0, 5.d0]; v = [1.d0, 0.d0, 0.d0]
      call trace_neutral_free_flight(x, v, mass, 1.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_NOT_INSIDE, ierr)
   end subroutine test_rejections

   !----- Ambiguous grazing intersection rejected ----------------------------!
   @test
   subroutine test_grazing_rejected()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      ! hits the x=0,y=0 edge of tetrahedron A exactly (both faces at once)
      x = [0.2d0, 0.2d0, 0.2d0]
      v = [-1.d0, -1.d0, 1.d0]
      call trace_neutral_free_flight(x, v, mass, 2.d0, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_ERR_GRAZING, ierr)
   end subroutine test_grazing_rejected

   !----- Slow crossing: |ndotv| < 1e-10 must still find the exit face --------!
   ! Start in A at (0.5,0.5,0.5); the shared face x+y+z=2 is reached at t=2.5e11
   ! with a normal speed (2e-12)/sqrt(3) ~ 1.15e-12, below any velocity
   ! tolerance. The tracer must cross into B and later hit B's wall (face 3,
   ! outward normal (2/3,2/3,1/3)) at t=3.75e11, ending at (0.875,0.875,0.5),
   ! instead of silently finishing outside A like the previous velocity
   ! threshold did.
   @test
   subroutine test_slow_face_crossing()
      real(dp) :: x(3), v(3)
      type(characteristic_result_t) :: res
      integer :: ierr

      x = [0.5d0, 0.5d0, 0.5d0]
      v = [1.d-12, 1.d-12, 0.d0]
      call trace_neutral_free_flight(x, v, mass, 5.d11, verts, ind_knot, &
           nb_tetr, nb_face, ntetr, res, ierr)
      @assertEqual(CHAR_OK, ierr)
      @assertFalse(res%finished)
      @assertEqual(CHAR_EVENT_WALL, res%event%kind)
      @assertEqual(2, res%event%tetrahedron)       ! crossed A->B despite tiny speed
      @assertEqual(3, res%event%face)
      @assertEqual(3.75d11, res%time, tolerance=1.d-2)
      @assertEqual([0.875d0, 0.875d0, 0.5d0], res%position, tolerance=1.d-10)
      @assertEqual([2.d0/3.d0, 2.d0/3.d0, 1.d0/3.d0], res%event%normal, tolerance=1.d-12)
   end subroutine test_slow_face_crossing

end module test_characteristic_step_mod
