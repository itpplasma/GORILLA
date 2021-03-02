module runge_kutta_mod
!
private
!
public :: runge_kutta_allroutines

    contains
!
    subroutine runge_kutta_allroutines(rk_order,y,nvar,x1,x2,derivs)
!
        IMPLICIT NONE
!
        external        :: derivs
        integer         :: rk_order,nvar
        double precision :: x1,x2,h1
        double precision, dimension(nvar) :: y,dydx,yout
!
        h1=x2-x1
!
        select case(rk_order)
            case(4)
                call derivs(x1,y,dydx)
                call rk4(y,dydx,x1,h1,yout,derivs)
                y = yout
!
            case DEFAULT
                print *, 'Chosen Runge Kutta is not available'
                stop
!
        end select
!
    end subroutine runge_kutta_allroutines
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    SUBROUTINE rk4(y,dydx,x,h,yout,derivs)
!
        IMPLICIT NONE
!
        external        :: derivs
        double precision, DIMENSION(:), INTENT(IN) :: y,dydx
        double precision, INTENT(IN) :: x,h
        double precision, DIMENSION(:), INTENT(OUT) :: yout
!
!        INTERFACE
!                SUBROUTINE derivs(x,y,dydx)
!!
!                IMPLICIT NONE
!!
!                double precision, INTENT(IN) :: x
!                double precision, DIMENSION(:), INTENT(IN) :: y
!                double precision, DIMENSION(:), INTENT(OUT) :: dydx
!!
!                END SUBROUTINE derivs
!        END INTERFACE
!
        double precision :: h6,hh,xh
        double precision, DIMENSION(size(y)) :: dym,dyt,yt
!
        hh=h*0.5d0
        h6=h/6.0d0
        xh=x+hh
        yt=y+hh*dydx
        call derivs(xh,yt,dyt)
        yt=y+hh*dyt
        call derivs(xh,yt,dym)
        yt=y+h*dym
        dym=dyt+dym
        call derivs(x+h,yt,dyt)
        yout=y+h6*(dydx+dyt+2.0d0*dym)
!
    END SUBROUTINE rk4
!
end module
