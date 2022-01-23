module test_binsrc
   use funit
   implicit none
   
contains

   @test
   subroutine test_binsrc_1()
      
      double precision, dimension(10) :: p
      double precision :: xi = 5.
      integer :: i

      p = [1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5, 10.5]
      call binsrc(p, 1, 10, xi, i)
      
      @assertEqual(5, i, message = "test binsrc")

   end subroutine test_binsrc_1
   
   
end module test_binsrc
