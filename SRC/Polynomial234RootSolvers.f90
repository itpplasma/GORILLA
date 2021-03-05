!! ----------------------------------------------------------------------------------
!! MODULE Polynomial234RootSolvers
!!
!!    This is a placeholder module. Please, include external library:
!!
!!    N. Flocke, “Algorithm 954: An Accurate and Efficient Cubic and Quartic Equation Solver for Physical Applications”
!!    <https://doi.org/10.1145/2699468>
!!
!!    * Download supplemental material `954.zip` from above webpage.
!!    * Copy `954/F90/Src/Polynomial234RootSolvers.f90` to `GORILLA/SRC/` and replace existing file with indentical name.
!!      (This file is a placeholder which is necessary for compilation.)
!!
!! ----------------------------------------------------------------------------------

module Polynomial234RootSolvers

use SetWorkingPrecision, ONLY : wp, wpformat

contains

subroutine cubicRoots (c2, c1, c0, nReal, root, printInfo)
!
  implicit none
!
  logical, optional  , intent (in)  :: printInfo
  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: c2, c1, c0
  real    (kind = wp), intent (out) :: root (1:3,1:2)
!
  print *, 'ERROR: Please, include external library Polynomial234RootSolvers. See README.md'
  stop
!
end subroutine cubicRoots
!
!
!
subroutine quadraticRoots (q1, q0, nReal, root)
!
  implicit none
!
  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: q1, q0
  real    (kind = wp), intent (out) :: root (1:2,1:2)
!
  print *, 'ERROR: Please, include external library Polynomial234RootSolvers. See README.md'
  stop
!
end subroutine quadraticRoots
!
!
!
subroutine quarticRoots (q3, q2, q1, q0, nReal, root, printInfo)
!
  implicit none
!
  logical, optional  , intent (in)  :: printInfo
  integer            , intent (out) :: nReal
  real    (kind = wp), intent (in)  :: q3, q2, q1, q0
  real    (kind = wp), intent (out) :: root (1:4,1:2)
!
  print *, 'ERROR: Please, include external library Polynomial234RootSolvers. See README.md'
  stop
!
end subroutine quarticRoots
!
end module Polynomial234RootSolvers
