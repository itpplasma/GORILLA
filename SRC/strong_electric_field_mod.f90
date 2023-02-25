module strong_electric_field_mod
!
    implicit none
!
    private
!
    !Setting mode of getting the electric potential
    !2 ... potential as a rescaled poloidal flux (E-field orthogonal to magnetic surfaces)
    integer, parameter                           :: electric_potential_option = 2
    !Setting mode of getting the electric field
    !1 ... using central differences to estimate the gradient of a scalar potential
    integer, parameter                           :: electric_field_option = 1
!
    !Stepsize for central difference approach (only calculated at initial call)
    logical                                      :: boole_first_call = .true.
    double precision, dimension(3)               :: dx_array
!
    public :: get_electric_field, save_electric_field, get_v_E
!
contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_electric_field(R,phi,Z,E_x1,E_x2,E_x3,phi_elec)
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: E_x1, E_x2, E_x3,phi_elec
!
        select case(electric_field_option)
!
            !Get electric field from scalar electric potential via central difference
            case(1)
                call electric_field_from_central_difference(R,phi,Z,E_x1,E_x2,E_x3)
!
            case default
                print *, 'ERROR: Invalid option to get electric field. Change the parameters in module strong_electric_field_mod!'
                stop
        end select
!
        !Get scalar potential in way that is consistent with calculated E-field (hard-coded settings in modulo variables above)
        call get_electric_potential(R,phi,Z,phi_elec)
!
    end subroutine get_electric_field
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine electric_field_from_central_difference(R,phi,Z,E_x1,E_x2,E_x3)
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: E_x1, E_x2, E_x3
!
        double precision, dimension(6)           :: electric_potential
        double precision, dimension(3)           :: dx
        integer                                  :: i
        double precision, parameter              :: machine_epsilon_3 = 1d-6
!
        !Determine the stepsize
        if(boole_first_call) then
            !To balance roundoff-error and precission, an estimated h = (machine_epsilon)^(1/3)*characteristic_lengthscale
            !is used with characteristic_lengthscale ~ the average coordinate-distance between discretization points of grid
            !see also http://www.it.uom.gr/teaching/linearalgebra/NumericalRecipiesInC/c5-7.pdf
            call calc_characteristic_scales(dx_array)
            dx_array = dx_array*machine_epsilon_3
            boole_first_call = .false.
        endif
!
        do i = 1,3
            dx = 0.0d0
            !Forward step
            dx(i) = dx_array(i)
            call get_electric_potential(R + dx(1), phi + dx(2), Z + dx(3), electric_potential(2*i-1))
            !Backward step
            dx = -dx
            call get_electric_potential(R + dx(1), phi + dx(2), Z + dx(3), electric_potential(2*i))
        end do
        !Forming the central difference as approximation of derivative (note the minus sign as E = -grad(Phi))
        E_x1 = -(electric_potential(1)-electric_potential(2))/(2*dx_array(1))
        E_x2 = -(electric_potential(3)-electric_potential(4))/(2*dx_array(2))
        E_x3 = -(electric_potential(5)-electric_potential(6))/(2*dx_array(3))
!
    end subroutine electric_field_from_central_difference
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_electric_potential(R,phi,Z,electric_potential)
!
        use field_eq_mod,                 only : psif
        use gorilla_settings_mod,         only: eps_Phi
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: electric_potential
!
        double precision                         :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 !Dummy variables for field() to set psif
!
        select case(electric_potential_option)
!
            !Get scalar electric potential as a rescale of the poloidal magnetic flux (potential contours align with flux surfaces)
            case(2)
                call field(R,phi,Z,d1,d2,d3,d4,d5,d6  &
                          ,d7,d8,d9,d10,d11,d12)
                electric_potential = psif*eps_Phi
!
            case default
                print *, 'ERROR: Invalid way to get electric potential. Change the parameters in module strong_electric_field_mod!'
                stop
        end select
!
    end subroutine get_electric_potential
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine save_electric_field(vec_E_x1,vec_E_x2,vec_E_x3)
!
        use tetra_grid_mod,                 only : verts_rphiz, nvert
!
        implicit none
!
        double precision, dimension(:), intent(in)      :: vec_E_x1,vec_E_x2,vec_E_x3
!
        integer                                         :: iv
!
        open(123, file='./electric_field.dat')
        do iv = 1, nvert
            write(123,*) verts_rphiz(:,iv),vec_E_x1(iv),vec_E_x2(iv),vec_E_x3(iv)
        end do
        close(123)
!
        !Stop program after successfully saving the E-field
        !STOP
!
    end subroutine save_electric_field
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_v_E(R,E_x1,E_x2,E_x3,h_x1,h_x2,h_x3,Bmod,v_E_x1,v_E_x2,v_E_x3,v2_E_mod)
!
        use constants, only: clight
!
        implicit none
!
        double precision, intent(in)             :: R, E_x1, E_x2, E_x3, h_x1, h_x2, h_x3, Bmod
        double precision, intent(out)            :: v_E_x1, v_E_x2, v_E_x3, v2_E_mod
!
        !Writting out the cross-product for v_E (covariant components) using the cylindrical metric determinant g_lk = (1,R^2,1)
        !(v_E)_l = g_lk * epsilon^(kij) * E_i * h_j / (sqrt(g) * Bmod)
        !The factor c (speed of light) is left out for having unitless values (include later when determing energy and pusher)
        v_E_x1 = (E_x2 * h_x3 - E_x3 * h_x2) / (R * Bmod) * clight
        v_E_x2 = (E_x3 * h_x1 - E_x1 * h_x3) / (Bmod) * R * clight
        v_E_x3 = (E_x1 * h_x2 - E_x2 * h_x1) / (R * Bmod) * clight
!
        !Writing out the scalar-product of v_E with itself using the cylindrical metric determinant
        !v2_E_mod = (v_E)_l * (v_E)^l = (v_E)_l * g^{lk} * (v_E)_k
        v2_E_mod = v_E_x1 * v_E_x1 + v_E_x2 * 1/R**2 * v_E_x2 + v_E_x3 * v_E_x3
!
    end subroutine get_v_E
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccAUXILIARYFUNCTIONSccccccccccccccccccccccccccccc
!
    subroutine find_domain_limits(limits)
!
        use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi
        use gorilla_settings_mod, only: coord_system
!
        implicit none
!
        double precision, dimension(3,2), intent(out)          :: limits
!
        select case(coord_system)
            !coord_system = 1 -> cylindrical (R,phi,Z)
            case(1)
                limits(:,1) = minval(verts_rphiz,2)
                limits(:,2) = maxval(verts_rphiz,2)
!
            !coord_system = 1 -> symflux (s,theta,phi)
            case(2)
                limits(:,1) = minval(verts_sthetaphi,2)
                limits(:,2) = maxval(verts_sthetaphi,2)
        end select
!
    end subroutine find_domain_limits
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_characteristic_scales(characteristic_scales)
!
        use tetra_grid_settings_mod, only: n2
        use tetra_grid_mod, only: nvert
        use gorilla_settings_mod, only: coord_system
!
        implicit none
!
        double precision, dimension(3), intent(out)        :: characteristic_scales
!
        double precision, dimension(3,2)                   :: limits
        double precision                                   :: average_2D_n, N_x1, N_x2, N_x3
!
        call find_domain_limits(limits)
!
        !estimate of points per coordinate direction in 2D slice
        average_2D_n = sqrt(dble(nvert/n2))
        !Taking care of position of coordinates in triple
        select case(coord_system)
            !coord_system = 1 -> cylindrical (R,phi,Z)
            case(1)
                N_x1 = average_2D_n
                N_x2 = n2
                N_x3 = average_2D_n
!
            !coord_system = 1 -> symflux (s,theta,phi)
            case(2)
                N_x1 = average_2D_n
                N_x2 = average_2D_n
                N_x3 = n2
        end select
!
        characteristic_scales(1) = abs(limits(1,2)-limits(1,1))/N_x1
        characteristic_scales(2) = abs(limits(2,2)-limits(2,1))/N_x2
        characteristic_scales(3) = abs(limits(3,2)-limits(3,1))/N_x3
!
    end subroutine calc_characteristic_scales
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module strong_electric_field_mod