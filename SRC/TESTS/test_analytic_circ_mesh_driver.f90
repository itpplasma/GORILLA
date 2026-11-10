! grid_kind = 5 mesh of the analytic circular tokamak against its field, with a
! tabulated q that differs from the q0 + q1 (r/a)^2 parabola.
! Usage: test_analytic_circ_mesh_driver.x <coord_system> <sign of B0>
!   1: B^phi and B_pol of curl A at the cell centroids.
!   2: cell q = dpsi_tor/dpsi_pol against flux quadratures, |B| at the vertices,
!      and faces that coincide across the theta* seam after the periodic shift.
program test_analytic_circ_mesh_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, grid_size
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz, verts_sthetaphi
    use tetra_physics_mod, only: tetra_physics, metric_determinant
    use field_analytic_circ_mod, only: field_analytic_circ
    implicit none
    real(dp), parameter :: r0 = 181.092169_dp, minor = 46.0_dp
    real(dp), parameter :: s_min = 0.1_dp, s_max = 0.6_dp
    integer, parameter :: n1 = 16, n2 = 4, n3 = 32
    real(dp) :: pi, b0, s_edge, xv(3, 4), xc(3), rho, b(3), q_true, err(3)
    real(dp) :: s_node(0:n1), pt_node(0:n1), pp_node(0:n1)
    integer :: coord_system, ind, j, nb, nseam
    character(len=16) :: arg

    pi = acos(-1.0_dp)
    call get_command_argument(1, arg)
    read(arg, *) coord_system
    call get_command_argument(2, arg)
    read(arg, *) b0
    b0 = sign(19800.0_dp, b0)
    s_edge = r0 - sqrt(r0**2 - minor**2)
    call write_inputs()
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()

    do j = 0, n1  ! radial mesh nodes are uniform in rho
        rho = rho_of_s(s_min) + j*(rho_of_s(s_max) - rho_of_s(s_min))/n1
        s_node(j) = (r0 - sqrt(r0**2 - rho**2))/s_edge
        pt_node(j) = psi_tor(rho)
        pp_node(j) = psi_pol(rho)
    end do

    err = 0.0_dp
    nseam = 0
    do ind = 1, ntetr/grid_size(2)
        if (coord_system == 2) then
            xv = verts_sthetaphi(:, tetra_grid(ind)%ind_knot)
            q_true = (node_flux(maxval(xv(1, :)), pt_node) - node_flux(minval(xv(1, :)), pt_node)) &
                /(node_flux(maxval(xv(1, :)), pp_node) - node_flux(minval(xv(1, :)), pp_node))
            err(1) = max(err(1), abs(tetra_physics(ind)%curlA(3)/tetra_physics(ind)%curlA(2) &
                /q_true - 1.0_dp))
            xc = verts_rphiz(:, tetra_grid(ind)%ind_knot(1))
            b = field_at(xc(1), xc(3))
            err(2) = max(err(2), abs(tetra_physics(ind)%bmod1/norm2(b) - 1.0_dp))
            do j = 1, 4
                nb = tetra_grid(ind)%neighbour_tetr(j)
                if (nb < 1) cycle
                if (abs(face_theta(ind, j) - face_theta(nb, tetra_grid(ind)%neighbour_face(j)) &
                    - 2.0_dp*pi*tetra_grid(ind)%neighbour_perbou_theta(j)) > 1.0e-12_dp) &
                    nseam = nseam + 1
            end do
        else
            ! (R,phi,Z): sqrt(g) = R; linear A, so O(cell size) at the centroid
            xc = sum(verts_rphiz(:, tetra_grid(ind)%ind_knot), dim=2)/4.0_dp
            b = field_at(xc(1), xc(3))
            err(1) = max(err(1), abs(tetra_physics(ind)%curlA(2) &
                /metric_determinant(ind, xc)/(b(2)/xc(1)) - 1.0_dp))
            err(2) = max(err(2), abs((-tetra_physics(ind)%curlA(1)*xc(3) &
                + tetra_physics(ind)%curlA(3)*(xc(1) - r0))/xc(1) &
                /(-b(1)*xc(3) + b(3)*(xc(1) - r0)) - 1.0_dp))
        end if
    end do

    print '(a,i2,a,f8.0,a,2es10.3,a,i0)', 'coord_system =', coord_system, '  B0 =', b0, &
        '  errors', err(1:2), '  seam faces ', nseam
    if (coord_system == 2) then
        if (err(1) > 1.0e-5_dp) error stop 'cell q differs from dpsi_tor/dpsi_pol'
        if (err(2) > 1.0e-10_dp) error stop '|B| at vertices differs'
        if (nseam /= 0) error stop 'faces inconsistent across the theta seam'
    else
        if (err(1) > 2.0e-2_dp) error stop 'B^phi differs'
        if (err(2) > 1.0e-1_dp) error stop 'B_pol differs'
    end if
    print '(a)', 'PASS: grid_kind 5 mesh field agrees with the analytic field'

contains

    function field_at(r, z) result(bv)
        real(dp), intent(in) :: r, z
        real(dp) :: bv(3), d(9)
        call field_analytic_circ(r, 0.0_dp, z, bv(1), bv(2), bv(3), d(1), d(2), d(3), &
            d(4), d(5), d(6), d(7), d(8), d(9))
    end function field_at

    ! Mean theta* of face iface of tetrahedron it (seam vertices at 0 count as 2 pi)
    real(dp) function face_theta(it, iface) result(th)
        integer, intent(in) :: it, iface
        real(dp) :: tv(4)
        tv = verts_sthetaphi(2, tetra_grid(it)%ind_knot)
        if (any(tv >= pi)) where (tv == 0.0_dp) tv = 2.0_dp*pi
        th = (sum(tv) - tv(iface))/3.0_dp
    end function face_theta

    real(dp) function node_flux(s, f) result(v)
        real(dp), intent(in) :: s, f(0:n1)
        integer :: k
        k = minloc(abs(s_node - s), dim=1) - 1
        if (abs(s_node(k) - s) > 1.0e-10_dp) error stop 'vertex s is not a mesh node'
        v = f(k)
    end function node_flux

    real(dp) function q_table(r) result(q)
        real(dp), intent(in) :: r
        q = 1.05_dp + 0.3_dp*(r/minor) + 1.2_dp*(r/minor)**3
    end function q_table

    real(dp) function rho_of_s(s) result(r)
        real(dp), intent(in) :: s
        r = sqrt(r0**2 - (r0 - s*s_edge)**2)
    end function rho_of_s

    ! Toroidal flux per radian inside r: int B_phi dA / (2 pi), Simpson in rho
    real(dp) function psi_tor(r) result(acc)
        real(dp), intent(in) :: r
        integer, parameter :: nr = 200, nt = 256
        real(dp) :: rr, th, bv(3)
        integer :: k, l
        acc = 0.0_dp
        do k = 1, nr
            rr = k*r/nr
            do l = 1, nt
                th = 2.0_dp*pi*(l - 0.5_dp)/nt
                bv = field_at(r0 + rr*cos(th), rr*sin(th))
                acc = acc + merge(1, merge(4, 2, mod(k, 2) == 1), k == nr) &
                    *r/nr/3.0_dp*rr*bv(2)/nt
            end do
        end do
    end function psi_tor

    ! Poloidal flux per radian inside r: int_R0^{R0+r} R B_Z dR on the midplane
    real(dp) function psi_pol(r) result(acc)
        real(dp), intent(in) :: r
        integer, parameter :: nr = 2000
        real(dp) :: bv(3)
        integer :: k
        acc = 0.0_dp
        do k = 1, nr
            bv = field_at(r0 + k*r/nr, 0.0_dp)
            acc = acc + merge(1, merge(4, 2, mod(k, 2) == 1), k == nr) &
                *r/nr/3.0_dp*(r0 + k*r/nr)*bv(3)
        end do
    end function psi_pol

    subroutine write_inputs()
        integer :: iu, k
        real(dp) :: r
        open(newunit=iu, file='q_profile_sfl_test.dat', status='replace')
        do k = 0, 400
            r = minor*k/400.0_dp
            write(iu, '(2es24.16)') r, q_table(r)
        end do
        close(iu)
        open(newunit=iu, file='tetra_grid.inp', status='replace')
        write(iu, '(a)') '&TETRA_GRID_NML'
        write(iu, '(a,i0,a,i0,a,i0)') ' grid_kind = 5, n1 = ', n1, ', n2 = ', n2, &
            ', n3 = ', n3
        write(iu, '(a)') ' boole_n_field_periods = .false., n_field_periods_manual = 1,'
        write(iu, '(a)') ' boole_write_mesh_obj = .false.,'
        write(iu, '(a,f6.3,a,f6.3,a)') ' sfc_s_min = ', s_min, ', sfc_s_max = ', s_max, ','
        write(iu, '(a,f14.8,a,f10.4,a,f12.2,a)') ' R0_analytic_circ = ', r0, &
            ', a_analytic_circ = ', minor, ', B0_analytic_circ = ', b0, ','
        write(iu, '(a)') ' q0_analytic_circ = 1.0, q1_analytic_circ = 2.0,'
        write(iu, '(a)') " q_profile_file_analytic_circ = 'q_profile_sfl_test.dat',"
        write(iu, '(a)') " helical_pert_file_analytic_circ = '',"
        write(iu, '(a)') '/'
        close(iu)
        open(newunit=iu, file='gorilla.inp', status='replace')
        write(iu, '(a)') '&GORILLANML'
        write(iu, '(a,i0,a)') ' coord_system = ', coord_system, ', ispecies = 2,'
        write(iu, '(a)') ' eps_Phi = 0.0, boole_periodic_relocation = .true.,'
        write(iu, '(a)') ' ipusher = 2, boole_pusher_ode45 = .false., rel_err_ode45 = 1e-10,'
        write(iu, '(a)') ' boole_dt_dtau = .true., boole_newton_precalc = .false.,'
        write(iu, '(a)') ' poly_order = 2, i_precomp = 0, boole_guess = .true.,'
        write(iu, '(a)') ' i_time_tracing_option = 1, handover_processing_kind = 1,'
        write(iu, '(a)') ' boole_axi_noise_vector_pot = .false., axi_noise_eps_A = 0.0,'
        write(iu, '(a)') ' boole_axi_noise_elec_pot = .false., axi_noise_eps_Phi = 0.0,'
        write(iu, '(a)') ' boole_non_axi_noise_vector_pot = .false., non_axi_noise_eps_A = 0.0,'
        write(iu, '(a)') ' boole_pert_from_mephit = .false., boole_helical_pert = .false.,'
        write(iu, '(a)') ' helical_pert_eps_Aphi = 0.0, helical_pert_m_fourier = 0,'
        write(iu, '(a)') ' helical_pert_n_fourier = 0, boole_time_hamiltonian = .false.,'
        write(iu, '(a)') ' boole_gyrophase = .false., boole_vpar_int = .false.,'
        write(iu, '(a)') ' boole_vpar2_int = .false., boole_adaptive_time_steps = .false.,'
        write(iu, '(a)') ' desired_delta_energy = 1e-10, max_n_intermediate_steps = 10000,'
        write(iu, '(a)') ' boole_grid_for_find_tetra = .false., a_factor = 1, b_factor = 0,'
        write(iu, '(a)') ' c_factor = 0, boole_strong_electric_field = .false.,'
        write(iu, '(a)') " boole_save_electric = .false., er_profile_file = '',"
        write(iu, '(a)') " filename_electric_field = 'electric_field.dat',"
        write(iu, '(a)') " filename_electric_drift = 'electric_drift.dat',"
        write(iu, '(a)') " er_equil_mapping_file = './flux_functions.dat',"
        write(iu, '(a)') '/'
        close(iu)
        open(newunit=iu, file='field_divB0.inp', status='replace')
        write(iu, '(a)') '0'
        write(iu, '(a)') '1'
        write(iu, '(a)') '1.0'
        write(iu, '(a)') '1'
        write(iu, '(a)') '0.99'
        write(iu, '(a)') '4'
        write(iu, '(a)') "'none'"
        write(iu, '(a)') "''"
        write(iu, '(a)') "'none'"
        write(iu, '(a)') "''"
        write(iu, '(a)') '0'
        write(iu, '(a)') '0'
        write(iu, '(a)') '0'
        close(iu)
    end subroutine write_inputs

end program test_analytic_circ_mesh_driver
