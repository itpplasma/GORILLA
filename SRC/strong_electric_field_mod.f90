module strong_electric_field_mod
!
    implicit none
!
    private
!
    !Electric potential: Phi_0(psi_pol) from er_profile_file if set,
    !otherwise the rescaled poloidal flux eps_Phi*psi_pol.
    !Setting mode of getting the electric field
    !1 ... using central differences to estimate the gradient of a scalar potential
    integer, parameter                           :: electric_field_option = 1
!
    !Stepsize for central difference approach (only calculated at initial call)
    logical                                      :: boole_first_call = .true.
    double precision, dimension(3)               :: dx_array
!
    !Phi_0(psi_pol) spline, loaded on first use
    logical                                      :: phi0_loaded = .false.
    integer                                      :: n_phi0 = 0
    double precision, allocatable                :: psi_pol_phi0(:), phi0_vals(:), phi0_dd(:)
!
    public :: get_electric_potential, get_electric_field, save_electric_field, get_v_E, save_v_E
    public :: load_phi0_from_er, eval_phi0
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
        use gorilla_settings_mod,         only: eps_Phi, er_profile_file, er_equil_mapping_file
        use field_divB0_mod, only: field
!
        implicit none
!
        double precision, intent(in)             :: R, phi, Z
        double precision, intent(out)            :: electric_potential
!
        double precision                         :: d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12 !Dummy variables for field() to set psif
!
        call field(R,phi,Z,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
!
        if (len_trim(er_profile_file) > 0) then
            if (.not. phi0_loaded) &
                call load_phi0_from_er(trim(er_profile_file), trim(er_equil_mapping_file))
            electric_potential = eval_phi0(psif) * eps_Phi
        else
            electric_potential = psif*eps_Phi
        end if
!
    end subroutine get_electric_potential
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    ! Phi_0(r_eff) = -int_0^r_eff Er dr' (trapezoidal, Er interpolated linearly
    ! onto the mapping grid), splined in psi_pol. Er file: r_eff [cm], Er
    ! [statV/cm]. Mapping file: R_beg, r_eff, q, psi_pol, psi_tor.
    subroutine load_phi0_from_er(filename_er, filename_mapping)
!
        character(len=*), intent(in) :: filename_er, filename_mapping
        double precision, allocatable :: er(:,:), eq(:,:), er_eq(:), phi0_eq(:)
        integer :: i, j, n_eq
!
        call read_table(filename_er, 2, er)
        call read_table(filename_mapping, 5, eq)
        n_eq = size(eq, 2)
        allocate(er_eq(n_eq), phi0_eq(n_eq))
        do j = 1, n_eq
            if (eq(2,j) <= er(1,1)) then
                er_eq(j) = er(2,1)
            else if (eq(2,j) >= er(1,size(er,2))) then
                er_eq(j) = er(2,size(er,2))
            else
                do i = 1, size(er,2) - 1
                    if (er(1,i+1) >= eq(2,j)) exit
                end do
                er_eq(j) = er(2,i) + (eq(2,j) - er(1,i)) / (er(1,i+1) - er(1,i)) &
                    * (er(2,i+1) - er(2,i))
            end if
        end do
        phi0_eq(1) = 0.0d0
        do j = 2, n_eq
            phi0_eq(j) = phi0_eq(j-1) - 0.5d0*(er_eq(j-1) + er_eq(j))*(eq(2,j) - eq(2,j-1))
        end do
!
        n_phi0 = n_eq
        if (eq(4,n_eq) >= eq(4,1)) then
            psi_pol_phi0 = eq(4,:)
            phi0_vals = phi0_eq
        else
            psi_pol_phi0 = eq(4,n_eq:1:-1)
            phi0_vals = phi0_eq(n_eq:1:-1)
        end if
        if (allocated(phi0_dd)) deallocate(phi0_dd)
        allocate(phi0_dd(n_phi0))
        call spline_natural(n_phi0, psi_pol_phi0, phi0_vals, phi0_dd)
        phi0_loaded = .true.
        print *, 'Loaded Er profile ', trim(filename_er), ': Phi_0 = ', &
            phi0_vals(1), ' ... ', phi0_vals(n_phi0), ' statV'
!
    end subroutine load_phi0_from_er
!
    subroutine read_table(filename, ncol, table)
!
        character(len=*), intent(in) :: filename
        integer, intent(in) :: ncol
        double precision, allocatable, intent(out) :: table(:,:)
        character(len=512) :: line
        integer :: iunit, ios, n, pass
!
        open(newunit=iunit, file=filename, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR (strong_electric_field_mod): Cannot open ', filename
            stop
        end if
        do pass = 1, 2
            n = 0
            do
                read(iunit, '(A)', iostat=ios) line
                if (ios /= 0) exit
                line = adjustl(line)
                if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
                n = n + 1
                if (pass == 2) read(line, *) table(:,n)
            end do
            if (pass == 1) allocate(table(ncol, n))
            rewind(iunit)
        end do
        close(iunit)
!
    end subroutine read_table
!
    double precision function eval_phi0(psi_pol_val) result(phi_out)
!
        double precision, intent(in) :: psi_pol_val
        integer :: i_lo, i_hi, i_mid
        double precision :: h, a, b
!
        if (psi_pol_val <= psi_pol_phi0(1)) then
            phi_out = phi0_vals(1)
            return
        else if (psi_pol_val >= psi_pol_phi0(n_phi0)) then
            phi_out = phi0_vals(n_phi0)
            return
        end if
        i_lo = 1
        i_hi = n_phi0
        do while (i_hi - i_lo > 1)
            i_mid = (i_lo + i_hi) / 2
            if (psi_pol_phi0(i_mid) > psi_pol_val) then
                i_hi = i_mid
            else
                i_lo = i_mid
            end if
        end do
        h = psi_pol_phi0(i_hi) - psi_pol_phi0(i_lo)
        a = (psi_pol_phi0(i_hi) - psi_pol_val) / h
        b = (psi_pol_val - psi_pol_phi0(i_lo)) / h
        phi_out = a * phi0_vals(i_lo) + b * phi0_vals(i_hi) &
                + ((a**3 - a) * phi0_dd(i_lo) + (b**3 - b) * phi0_dd(i_hi)) * h * h / 6.0d0
!
    end function eval_phi0
!
    ! Natural cubic spline second derivatives (Numerical Recipes).
    subroutine spline_natural(n, x, y, y2)
!
        integer, intent(in)               :: n
        double precision, intent(in)      :: x(n), y(n)
        double precision, intent(out)     :: y2(n)
        integer :: i, k
        double precision :: p, sig, u(n)
!
        y2(1) = 0.0d0
        u(1)  = 0.0d0
        do i = 2, n - 1
            sig = (x(i) - x(i-1)) / (x(i+1) - x(i-1))
            p   = sig * y2(i-1) + 2.0d0
            y2(i) = (sig - 1.0d0) / p
            u(i)  = (6.0d0 * ((y(i+1) - y(i))/(x(i+1) - x(i)) &
                            - (y(i) - y(i-1))/(x(i) - x(i-1))) / (x(i+1) - x(i-1)) - sig * u(i-1)) / p
        end do
        y2(n) = 0.0d0
        do k = n - 1, 1, -1
            y2(k) = y2(k) * y2(k+1) + u(k)
        end do
!
    end subroutine spline_natural
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine save_electric_field(vec_E_x1,vec_E_x2,vec_E_x3)
!
        use tetra_grid_mod,                 only : verts_rphiz, nvert
        use gorilla_settings_mod,           only : filename_electric_field
!
        implicit none
!
        double precision, dimension(:), intent(in)      :: vec_E_x1,vec_E_x2,vec_E_x3
!
        integer                                         :: iv, iunit
!
        open(newunit=iunit, file=filename_electric_field)
        do iv = 1, nvert
            write(iunit,*) verts_rphiz(:,iv),vec_E_x1(iv),vec_E_x2(iv),vec_E_x3(iv)
        end do
        close(iunit)
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
    subroutine save_v_E(vec_v_E_x1,vec_v_E_x2,vec_v_E_x3,vec_v2_E_mod)
!
        use tetra_grid_mod,                 only : verts_rphiz, nvert
        use gorilla_settings_mod,           only : filename_electric_drift
!
        implicit none
!
        double precision, dimension(:), intent(in)      :: vec_v_E_x1,vec_v_E_x2,vec_v_E_x3,vec_v2_E_mod
!
        integer                                         :: iv, iunit
!
        open(newunit=iunit, file=filename_electric_drift)
        do iv = 1, nvert
            write(iunit,*) verts_rphiz(:,iv),vec_v_E_x1(iv),vec_v_E_x2(iv),vec_v_E_x3(iv),vec_v2_E_mod(iv)
        end do
        close(iunit)
!
    end subroutine save_v_E
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