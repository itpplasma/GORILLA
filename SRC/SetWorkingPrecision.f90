    MODULE SetWorkingPrecision
!
    !This module is required for uitilization of Polynomial234RootSolvers.f90
!
      INTEGER, PARAMETER            :: wp = KIND(0.0D0)
      CHARACTER(LEN=*), PARAMETER   :: wpformat = '(a, e23.16)'
!
    END MODULE SetWorkingPrecision

