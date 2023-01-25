module test_various_functions_mod
   use various_functions_mod
   use funit
   implicit none
!
! Background:
! The module vaious_functions_mod (similar to supporting_functions_mod) is intended as a 
! collection of functions/subroutines are more general and not tied to a specifc module. 
!
! Georg Graßler (25.01.2023)
!       
contains

   @test
   ! Description of test_dmatinv3:
   ! The routine dmatinv3 calculates the inverse B of a (3 x 3) matrix A. In this test
   ! trial matrices are provided of which the inverses are already known.
   ! 
   ! Georg Graßler (25.01.2023)
   ! 
   subroutine test_dmatinv3()
      
      double precision,dimension(3,3) :: A, B, C, D
      double precision,dimension(3,3) :: D_known
      integer :: ierr
      
      A = reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/3, 3/))
      call dmatinv3(A, B, ierr)
      
      !Inverse of unit matrix is unit matrix itself -> A = B
      @assertEqual(A, B, tolerance=1e-13, message="test various_functions_mod with identity matrix")
      
      C = reshape((/0.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 0.d0/), (/3, 3/))
      call dmatinv3(C, D, ierr)
      
      !Analytically known inverse of C
      D_known = 1.d0/2.d0*reshape((/-1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 1.d0, 1.d0, 1.d0, -1.d0/), (/3, 3/))
      @assertEqual(D, D_known, tolerance=1e-13, message="test various_functions_mod with a symmetric matrix")

   end subroutine test_dmatinv3
   
   @test
   ! Description of test_dmatinv3_fail:
   ! The test uses a singular matrix (has no invers) intentionally to check, if the internal
   ! fail save of dmatinv3 can properly handle this case and if it throws the correct error 
   ! flag (ierr = 1).  
   ! 
   ! Georg Graßler (25.01.2023)
   ! 
   subroutine test_dmatinv3_fail()
      
      double precision,dimension(3,3) :: A, B, C, D
      integer :: ierr
      
      A = reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0/), (/3, 3/))
      call dmatinv3(A, B, ierr)
      
      @assertEqual(ierr, 1, message = "test various_functions_mod, matrix not invertible")


   end subroutine test_dmatinv3_fail
   
end module test_various_functions_mod
