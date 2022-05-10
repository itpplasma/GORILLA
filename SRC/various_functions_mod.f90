!
  module various_functions_mod
!
  contains
!
    subroutine dmatinv3(A,B,ierr)
         !! Performs a direct calculation of the inverse of a 3×3 matrix.
         double precision,dimension(3,3), intent(in) :: A  !! Matrix
         double precision,dimension(3,3), intent(out) :: B(3,3)   !! Inverse matrix
         integer, intent(out) :: ierr
         double precision             :: detinv
!
         !Integer for error status (0 ... no error, 1 ... singular matrix)
         ierr = 0
!
         !Calculate the inverse determinant of the matrix
         detinv = (A(1,1)*A(2,2)*A(3,3) - A(1,1)*A(2,3)*A(3,2)&
                   - A(1,2)*A(2,1)*A(3,3) + A(1,2)*A(2,3)*A(3,1)&
                   + A(1,3)*A(2,1)*A(3,2) - A(1,3)*A(2,2)*A(3,1))
!
         if(detinv.eq.0.d0) then
           ierr = 1
           B = 0.d0
           print *,"Error in dmatinv3: Matrix is singular - Inverse matrix doesn't exist."
           return
         endif
!
         detinv = 1/detinv
!
         ! Calculate the inverse of the matrix
         B(1,1) = +detinv * (A(2,2)*A(3,3) - A(2,3)*A(3,2))
         B(2,1) = -detinv * (A(2,1)*A(3,3) - A(2,3)*A(3,1))
         B(3,1) = +detinv * (A(2,1)*A(3,2) - A(2,2)*A(3,1))
         B(1,2) = -detinv * (A(1,2)*A(3,3) - A(1,3)*A(3,2))
         B(2,2) = +detinv * (A(1,1)*A(3,3) - A(1,3)*A(3,1))
         B(3,2) = -detinv * (A(1,1)*A(3,2) - A(1,2)*A(3,1))
         B(1,3) = +detinv * (A(1,2)*A(2,3) - A(1,3)*A(2,2))
         B(2,3) = -detinv * (A(1,1)*A(2,3) - A(1,3)*A(2,1))
         B(3,3) = +detinv * (A(1,1)*A(2,2) - A(1,2)*A(2,1))
       end subroutine
!
  end module various_functions_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
