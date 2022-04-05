module test_various_functions_mod
   use various_functions_mod
   use funit
   implicit none
   
contains

   @test
   subroutine test_dmatinv3()
      
      double precision,dimension(3,3) :: A, B, C, D
      double precision,dimension(3,3) :: D_known
      integer :: ierr
      
      A = reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 1.d0/), (/3, 3/))
      call dmatinv3(A, B, ierr)
      
      @assertEqual(A, B, tolerance=1e-13, message="test various_functions_mod with identity matrix")
      
      C = reshape((/0.d0, 1.d0, 1.d0, 1.d0, 0.d0, 1.d0, 1.d0, 1.d0, 0.d0/), (/3, 3/))
      call dmatinv3(C, D, ierr)
      
      D_known = 1.d0/2.d0*reshape((/-1.d0, 1.d0, 1.d0, 1.d0, -1.d0, 1.d0, 1.d0, 1.d0, -1.d0/), (/3, 3/))
      @assertEqual(D, D_known, tolerance=1e-13, message="test various_functions_mod with a symmetric matrix")

   end subroutine test_dmatinv3
   
   @test
   subroutine test_dmatinv3_fail()
      
      double precision,dimension(3,3) :: A, B, C, D
      integer :: ierr
      
      A = reshape((/1.d0, 0.d0, 0.d0, 0.d0, 1.d0, 0.d0, 0.d0, 0.d0, 0.d0/), (/3, 3/))
      call dmatinv3(A, B, ierr)
      
      @assertEqual(ierr, 1, message = "test various_functions_mod, matrix not invertible")


   end subroutine test_dmatinv3_fail
   
end module test_various_functions_mod
