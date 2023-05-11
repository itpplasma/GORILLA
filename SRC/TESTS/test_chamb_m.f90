module test_chamb_m
   use funit
   implicit none
!
! Background:
! Consider a point in 3D-space described by flux coordinates. One of its coordinates is
! the toroidal angle phi, while the remaining two (y(1),y(2)) describe the position in
! the designated poloidal crossection. The flux label y(1) may be normalized in respect to
! a specific flux surface, which forms the outer border of the compuational domain of interest.
! The function chamb_can(y, phi,ierr) checks now if a proposed guiding-center position lies either
! inside (ierr = 0) or outside (ierr = 1) this region (vacuum chamber). This function is currently
! not used in GORILLA itself, but in seperate diagnostic routines (GORILLA_APPLETS).
! Its function here is just as an example case for the pFUnit test framework.
!
! Georg Graßler (24.01.2023)
!   
contains

   @test
   ! Description of test_chamb_m_1:
   ! The tests provides a position inside the vacuum chamber (y(1)<1.d0). Only the flux label 
   ! coordinate y(1) plays a role when determining wether the position is inside/outside. The 
   ! coordinate y(2) as well as phi (not even initialized in the tests) do not change the result.
   !
   ! Georg Graßler (24.01.2023)
   !  
   subroutine test_chamb_m_1()
      
      integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y

      y = [0.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 0, message = "test chamb_m 1")

   end subroutine test_chamb_m_1
   
   @test
   ! Description of test_chamb_m_2:
   ! The tests provides a position exactly at the border (y(1) = 1.d0). In this edge case the 
   ! position should be classified as being already outside (ierr = 1) of the domain. 
   !
   ! Georg Graßler (24.01.2023)
   !  
   subroutine test_chamb_m_2()
   
	  	integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y
      
      y = [1.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 1, message = "test chamb_m 2")
      
   end subroutine test_chamb_m_2
   
   @test
   ! Description of test_chamb_m_3:
   ! The tests provides a position well outside the vacuum chamber (y(1)>1.d0). The expected 
   ! result is ierr = 1.
   !
   ! Georg Graßler (24.01.2023)
   !  
   subroutine test_chamb_m_3()
   
   	integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y
      
      y = [2.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 1, message = "test chamb_m 3")
      
   end subroutine test_chamb_m_3
   
   
end module test_chamb_m
