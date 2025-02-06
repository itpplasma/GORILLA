module odeint_mod

contains

subroutine odeint_allroutines(y, nvar, x1, x2, eps, derivs)

    use rkf45_mod, only: r8_rkf45

    implicit none

    external :: derivs
    integer, intent(in) :: nvar
    double precision, intent(in) :: x1, x2, eps
    double precision, dimension(nvar) :: y
    double precision, dimension(nvar) :: yp

    double precision :: epsrel, epsabs
    double precision :: x1in, x2in
    integer :: flag

    flag = 1
    epsrel = eps
    epsabs = 1d-31

    ! Required so that x1 and x2 are not overwritten by r8_rkf45
    x1in = x1
    x2in = x2

    call r8_rkf45 ( derivs, nvar, y, yp, x1in, x2in, epsrel, epsabs, flag )

    if (flag == 6) then
        epsrel = 10*epsrel
        epsabs = 10*epsabs
        flag = 2
        call r8_rkf45 ( derivs, nvar, y, yp, x1in, x2in, epsrel, epsabs, flag )
    elseif (flag == 7) then
        flag = 2
        call r8_rkf45 ( derivs, nvar, y, yp, x1in, x2in, epsrel, epsabs, flag )
    endif

end subroutine odeint_allroutines

end module
