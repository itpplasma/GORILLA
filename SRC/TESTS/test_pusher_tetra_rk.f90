module test_pusher_tetra_rk
   use pusher_tetra_rk_mod, only: initialize_const_motion_rk !, perpinv, perpinv2
   use funit
   implicit none
!
! Background:
! GORILLA performs the integration piecewise, by iteratively solving the linear set of
! of equations of motion in the current tetrahedron. One way to solve the set is via an
! Runge-Kutta scheme trying to determing the exit point and exit time of the orbit inside
! the cell. The necessary routines are contained in the module pusher_tetra_rk.
!
! Georg Graßler (25.01.2023)
!       
contains
	
	@test
	@disable
    ! Description of test_initialize_const_motion_rk:
    ! For the pusher to work the value of the constant of motion J_perp (i.e. perpinv)
    ! has to be set at the beginning of integration. For that a writing routine is used to
    ! set this value from outside the actual module (initialize_const_motion_rk). The second
    ! variable perpinv2 should be just perpinv2 squared, but can be set independently.
    !
    ! However, in the current version the program structure does not support the reading
    ! out of the set values outside the module. Until further changes the test is incomplete
    ! therefore DISABLED.
    !
    ! Georg Graßler (25.01.2023)
    ! 
	subroutine test_initialize_const_motion_rk()
	
		double precision  :: perpinv_in,perpinv2_in
		perpinv_in = 3.0d0
		perpinv2_in = 2.5d0
		
		call initialize_const_motion_rk(perpinv_in, perpinv2_in)
		
		!perpinv and perpinv2 are not public, due to this they are not able to inherit from pusher_tetra_rk_mod, so they are not testable in this way
		!@assertEqual(perpinv, perpinv_in, tolerance = 1e-13, message = "Testing initialize_const_motion_rk perpinv")
		!@assertEqual(perpinv2, perpinv2_in, tolerance = 1e-13, message = "Testing initialize_const_motion_rk perpinv2")
		
	end subroutine test_initialize_const_motion_rk
	
end module
