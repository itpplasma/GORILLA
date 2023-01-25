module test_binsrc
   use funit
   implicit none
!
! Background:
! The function binsrc(p,nmin,nmax,xi,i) takes an 1D-array p and searches for the index i 
! that fullfils p(i-1) <  xi  <  p(i) using binary search in the index range [nmin,nmax].
! This is employed for example when searching for the discrete data points between which
! field data has to be interpolated.
!
! Georg Graßler (24.01.2023)
!
contains

   @test
   ! Description of test_binsrc_1:
   ! The test test_binsrc_1() checks via a trial array if binsrc gives back the correct index
   ! according to the above outlined condition.
   !
   ! Georg Graßler (24.01.2023)
   !
   subroutine test_binsrc_1()
      
      double precision, dimension(10) :: p
      double precision :: xi = 5.
      integer :: i

      p = [1.5d0, 2.5d0, 3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5d0, 8.5d0, 9.5d0, 10.5d0]
      call binsrc(p, 1, 10, xi, i)
      
      @assertEqual(5, i, message = "test binsrc")

   end subroutine test_binsrc_1
   
end module test_binsrc
