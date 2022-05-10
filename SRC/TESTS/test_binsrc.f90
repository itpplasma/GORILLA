module test_binsrc
   use funit
   implicit none
   
contains

   @test
   subroutine test_binsrc_1()
      
      double precision, dimension(10) :: p
      double precision :: xi = 5.
      integer :: i

      p = [1.5d0, 2.5d0, 3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5d0, 8.5d0, 9.5d0, 10.5d0]
      call binsrc(p, 1, 10, xi, i)
      
      @assertEqual(5, i, message = "test binsrc")

   end subroutine test_binsrc_1
   
end module test_binsrc
