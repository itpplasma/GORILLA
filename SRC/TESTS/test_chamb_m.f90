module test_chamb_m
   use funit
   implicit none
   
contains

   @test
   subroutine test_chamb_m_1()
      
      integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y

      y = [0.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 0, message = "test chamb_m 1")

   end subroutine test_chamb_m_1
   
   @test
   subroutine test_chamb_m_2()
   
	  	integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y
      
      y = [1.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 1, message = "test chamb_m 2")
      
   end subroutine test_chamb_m_2
   
   @test
   subroutine test_chamb_m_3()
   
   	integer :: ierr
      double precision :: phi
      double precision, dimension(2) :: y
      
      y = [2.d0,1.d0]
      call chamb_can(y, phi,ierr)
      @assertEqual(ierr, 1, message = "test chamb_m 3")
      
   end subroutine test_chamb_m_3
   
   
end module test_chamb_m
