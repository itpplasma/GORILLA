module strong_electric_field_mod
!
    implicit none
!
    private
!
    !Setting mode of getting the electric potential
    integer, parameter                           :: electric_potential_mode = 1
    !Setting mode of getting the electric field
    integer, parameter                           :: electric_field_mode = 1
!
    public :: get_electric_field
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_electric_field(R,phi,Z,E_x1,E_x2,E_x3)
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: E_x1, E_x2, E_x3
!
        double precision, dimension(6)           :: electric_potential
        double precision, dimension(3)           :: dx_array, dx
        integer                                  :: i
!
        select case(electric_field_mode)
!
            !Get electric field from scalar electric potential
            case(1)
                do i = 1,3
                    dx = 0.0d0
                    !Forward step
                    dx(i) = dx_array(i)
                    call get_electric_potential(R + dx(1), phi + dx(1), Z + dx(1), electric_potential(2*i-1))
                    !Backward step
                    dx = -dx
                    call get_electric_potential(R + dx(1), phi + dx(1), Z + dx(1), electric_potential(2*i))
                end do
                !Forming the central difference as approximation of derivative
                E_x1 = -(electric_potential(1)-electric_potential(2))/(2*dx_array(1))
                E_x2 = -(electric_potential(3)-electric_potential(4))/(2*dx_array(2))
                E_x3 = -(electric_potential(5)-electric_potential(6))/(2*dx_array(3))
!
            case default
                print *, 'ERROR: Invalid mode to get electric field. Change the parameters in module strong_electric_field_mod!'
                stop
        end select
!
    end subroutine get_electric_field
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_electric_potential(R,phi,Z,electric_potential)
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: electric_potential
!
        select case(electric_potential_mode)
!
            !Get scalar electric potential as a rescale of the poloidal magnetic flux (potential contours align with flux surfaces)
            case(1)
!
            case default
                print *, 'ERROR: Invalid mode to get electric potential. Change the parameters in module strong_electric_field_mod!'
                stop
        end select
!
    end subroutine get_electric_potential
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module strong_electric_field_mod