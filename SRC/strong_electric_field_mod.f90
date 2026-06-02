module strong_electric_field_mod
!
    implicit none
!
    private
!
    !Setting mode of getting the electric potential
    !1 ... loaded from Er(r_eff) profile + flux mapping, integrated to Phi_0(psi_pol)
    !2 ... potential as a rescaled poloidal flux (E-field orthogonal to magnetic surfaces)
    !
    !Selection is runtime: case(1) is chosen if er_profile_file (from
    !gorilla_settings_mod) is non-empty at first call, otherwise case(2).
    !Setting mode of getting the electric field
    !1 ... using central differences to estimate the gradient of a scalar potential
    integer, parameter                           :: electric_field_option = 1
!
    !Stepsize for central difference approach (only calculated at initial call)
    logical                                      :: boole_first_call = .true.
    double precision, dimension(3)               :: dx_array
!
    !State for the loaded Phi_0(psi_pol) spline (case 1).  Filled on first call.
    logical                                      :: phi0_loaded = .false.
    integer                                      :: n_phi0 = 0
    double precision, allocatable                :: psi_pol_phi0(:), phi0_vals(:), phi0_dd(:)
!
    public :: get_electric_field, save_electric_field, get_v_E, save_v_E
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
        ! Always evaluate the field to set psif at (R,phi,Z) -- needed by both branches.
        call field(R,phi,Z,d1,d2,d3,d4,d5,d6,d7,d8,d9,d10,d11,d12)
!
        if (len_trim(er_profile_file) > 0) then
            ! Case (1): scalar potential from integrated Er(r_eff), splined as Phi_0(psi_pol).
            ! On first call, lazily load and spline; subsequent calls just evaluate.
            if (.not. phi0_loaded) then
                call load_phi0_from_er(trim(er_profile_file), trim(er_equil_mapping_file))
            end if
            electric_potential = eval_phi0(psif) * eps_Phi
        else
            ! Case (2): synthetic potential as a rescaled poloidal flux (contours align with flux surfaces).
            electric_potential = psif*eps_Phi
        end if
!
    end subroutine get_electric_potential
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    ! ------------------------------------------------------------------
    ! Load Er.dat (two columns: r_eff [cm], Er [statV/cm]) and the
    ! equilibrium mapping file (columns: R_beg, r_eff, q, psi_pol, psi_tor),
    ! integrate Phi_0(r_eff) = -int_0^r_eff Er(r') dr' (trapezoidal),
    ! and build a natural cubic spline Phi_0(psi_pol) for evaluation.
    !
    ! Assumes Er is a flux-surface quantity, so Phi_0 has only a psi_pol
    ! dependence -- consistent with the user-specified simplification.
    ! ------------------------------------------------------------------
    subroutine load_phi0_from_er(filename_er, filename_mapping)
!
        implicit none
!
        character(len=*), intent(in) :: filename_er, filename_mapping
!
        integer :: ios, n_er, n_eq, i, j, iunit
        double precision, allocatable :: r_er(:), Er_arr(:)
        double precision, allocatable :: r_eq(:), psi_pol_eq(:), phi0_eq(:), Er_at_eq(:)
        double precision :: dummy, frac
        character(len=512) :: line
!
        iunit = 124
!
        ! ---- Read Er.dat (r_eff [cm], Er [statV/cm]) ----
        open(unit=iunit, file=filename_er, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR (strong_electric_field_mod): Cannot open Er profile: ', filename_er
            stop
        end if
        n_er = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            n_er = n_er + 1
        end do
        rewind(iunit)
        allocate(r_er(n_er), Er_arr(n_er))
        i = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            i = i + 1
            read(line, *) r_er(i), Er_arr(i)
        end do
        close(iunit)
!
        ! ---- Read equil-mapping file (R_beg, r_eff, q, psi_pol, psi_tor) ----
        open(unit=iunit, file=filename_mapping, status='old', action='read', iostat=ios)
        if (ios /= 0) then
            print *, 'ERROR (strong_electric_field_mod): Cannot open equil mapping: ', filename_mapping
            stop
        end if
        n_eq = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            n_eq = n_eq + 1
        end do
        rewind(iunit)
        allocate(r_eq(n_eq), psi_pol_eq(n_eq), phi0_eq(n_eq), Er_at_eq(n_eq))
        j = 0
        do
            read(iunit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            line = adjustl(line)
            if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
            j = j + 1
            ! columns: R_beg, r_eff, q, psi_pol, psi_tor
            read(line, *) dummy, r_eq(j), dummy, psi_pol_eq(j), dummy
        end do
        close(iunit)
!
        ! ---- Linearly interpolate Er onto the equil-mapping r_eff grid ----
        do j = 1, n_eq
            if (r_eq(j) <= r_er(1)) then
                Er_at_eq(j) = Er_arr(1)
            else if (r_eq(j) >= r_er(n_er)) then
                Er_at_eq(j) = Er_arr(n_er)
            else
                do i = 1, n_er - 1
                    if (r_er(i+1) >= r_eq(j)) exit
                end do
                frac = (r_eq(j) - r_er(i)) / (r_er(i+1) - r_er(i))
                Er_at_eq(j) = Er_arr(i) + frac * (Er_arr(i+1) - Er_arr(i))
            end if
        end do
!
        ! ---- Trapezoidal integration: Phi_0(r) = -int_0^r Er(r') dr' ----
        phi0_eq(1) = 0.0d0
        do j = 2, n_eq
            phi0_eq(j) = phi0_eq(j-1) - 0.5d0*(Er_at_eq(j-1) + Er_at_eq(j))*(r_eq(j) - r_eq(j-1))
        end do
!
        ! ---- Build spline of Phi_0 vs psi_pol; ensure psi_pol is ascending ----
        n_phi0 = n_eq
        if (allocated(psi_pol_phi0)) deallocate(psi_pol_phi0)
        if (allocated(phi0_vals))    deallocate(phi0_vals)
        if (allocated(phi0_dd))      deallocate(phi0_dd)
        allocate(psi_pol_phi0(n_phi0), phi0_vals(n_phi0), phi0_dd(n_phi0))
        if (psi_pol_eq(n_eq) >= psi_pol_eq(1)) then
            psi_pol_phi0 = psi_pol_eq
            phi0_vals    = phi0_eq
        else
            do j = 1, n_phi0
                psi_pol_phi0(j) = psi_pol_eq(n_eq - j + 1)
                phi0_vals(j)    = phi0_eq(n_eq - j + 1)
            end do
        end if
        call spline_natural(n_phi0, psi_pol_phi0, phi0_vals, phi0_dd)
!
        phi0_loaded = .true.
!
        print *, '  [strong_electric_field_mod] Loaded Er profile from: ', filename_er
        print *, '    Er n_points = ', n_er
        print *, '    Equil mapping n_points = ', n_eq
        print *, '    Phi_0 at psi_pol_min = ', phi0_vals(1), ' statV'
        print *, '    Phi_0 at psi_pol_max = ', phi0_vals(n_phi0), ' statV'
!
        deallocate(r_er, Er_arr, r_eq, psi_pol_eq, phi0_eq, Er_at_eq)
!
    end subroutine load_phi0_from_er
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function eval_phi0(psi_pol_val) result(phi_out)
!
        implicit none
!
        double precision, intent(in) :: psi_pol_val
        double precision             :: phi_out
        integer :: i_lo, i_hi, i_mid
        double precision :: h, a, b
!
        if (.not. phi0_loaded .or. n_phi0 < 2) then
            phi_out = 0.0d0
            return
        end if
!
        if (psi_pol_val <= psi_pol_phi0(1)) then
            phi_out = phi0_vals(1)
            return
        else if (psi_pol_val >= psi_pol_phi0(n_phi0)) then
            phi_out = phi0_vals(n_phi0)
            return
        end if
!
        ! Binary search bracketing [i_lo, i_hi].
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
!
        h = psi_pol_phi0(i_hi) - psi_pol_phi0(i_lo)
        a = (psi_pol_phi0(i_hi) - psi_pol_val) / h
        b = (psi_pol_val - psi_pol_phi0(i_lo)) / h
        phi_out = a * phi0_vals(i_lo) + b * phi0_vals(i_hi) &
                + ((a**3 - a) * phi0_dd(i_lo) + (b**3 - b) * phi0_dd(i_hi)) * h * h / 6.0d0
!
    end function eval_phi0
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    ! Natural cubic spline second-derivative computation (Numerical Recipes).
    subroutine spline_natural(n, x, y, y2)
!
        implicit none
!
        integer, intent(in)               :: n
        double precision, intent(in)      :: x(n), y(n)
        double precision, intent(out)     :: y2(n)
        integer :: i, k
        double precision :: p, sig
        double precision, allocatable :: u(:)
!
        allocate(u(n))
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
        deallocate(u)
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
        integer                                         :: iv
!
        open(123, file=filename_electric_field)
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
    subroutine save_v_E(vec_v_E_x1,vec_v_E_x2,vec_v_E_x3,vec_v2_E_mod)
!
        use tetra_grid_mod,                 only : verts_rphiz, nvert
        use gorilla_settings_mod,           only : filename_electric_drift
!
        implicit none
!
        double precision, dimension(:), intent(in)      :: vec_v_E_x1,vec_v_E_x2,vec_v_E_x3,vec_v2_E_mod
!
        integer                                         :: iv
!
        open(123, file=filename_electric_drift)
        do iv = 1, nvert
            write(123,*) verts_rphiz(:,iv),vec_v_E_x1(iv),vec_v_E_x2(iv),vec_v_E_x3(iv),vec_v2_E_mod(iv)
        end do
        close(123)
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