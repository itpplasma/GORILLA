module test_various_functions_mod
   use various_functions_mod
   use funit
   implicit none
   
contains

   @test
   subroutine test_dmatinv3()
      
      double precision,dimension(3,3) :: A, B, C, D
      integer :: ierr
      logical:: resL
      
      A = reshape( (/1., 0., 0., 0., 1., 0., 0., 0., 1./), (/3, 3/) )
      call dmatinv3(A, B, ierr)
      
      resL = all(abs(A-B) <= 1e-13)
      @assertTrue(resL, message = "test various_functions_mod")
      
      C = reshape( (/0., 1., 1., 1., 0., 1., 1., 1., 0./), (/3, 3/) )
      call dmatinv3(C, D, ierr)
      
      resL = all(abs(D-1./2.*reshape( (/-1., 1., 1., 1., -1., 1., 1., 1., -1./), (/3, 3/) )) <=1e-13)
      @assertTrue(resL, message = "test various_functions_mod")

   end subroutine test_dmatinv3
   
   @test
   subroutine test_dmatinv3_fail()
      
      double precision,dimension(3,3) :: A, B, C, D
      integer :: ierr
      logical:: resL
      
      A = reshape( (/1., 0., 0., 0., 1., 0., 0., 0., 0./), (/3, 3/) )
      call dmatinv3(A, B, ierr)
      
      @assertEqual(ierr, 1, message = "test various_functions_mod, matrix not invertible")


   end subroutine test_dmatinv3_fail
   
   
end module test_various_functions_mod
