! Field consistency of the analytic circular tokamak mesh (grid_kind = 5)
! against the analytic field, for coord_system = 1 (R,phi,Z) and 2 (s,theta*,phi).
!
! Usage: test_analytic_circular_sfl.x <coord_system> <sign of B0>
! The program writes its own input files into the working directory.
!
! Oracles are independent of the GORILLA mesh code:
!   - straight-field-line angle theta*(rho,theta) from numerical quadrature of the
!     field-line winding dphi/dtheta = B^phi/B^theta of field_analytic_circ;
!   - toroidal flux per radian from a surface quadrature of B_phi,
!     poloidal flux per radian from a midplane quadrature of R*B_Z.
! A tabulated q profile that differs from the q0 + q1 (r/a)^2 parabola is used,
! so a code path that ignores the table fails.
program test_analytic_circular_sfl
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
    integer, parameter :: n1 = 16, n2 = 4, n3 = 32, nquad = 400
    real(dp) :: pi, b0, s_edge, xv(3, 4), xc(3), rho, thg, ths, wtot
    real(dp) :: bphi_true, bth_true, bphi_g, bth_g, q_true, q_g, smin_t, smax_t
    real(dp) :: err_bphi, err_bth, err_q, err_b, err_h(3), h_true(3), b_true
    real(dp) :: br, bp, bz, d(9), xr(3)
    integer :: coord_system, ind, j, iv, nfail, nplane
    real(dp) :: s_node(0:n1), pt_node(0:n1), pp_node(0:n1), th_own, th_nb
    integer :: nb, nseam
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

    ! Oracle fluxes on the radial mesh nodes (uniform in rho, see tetra_grid_mod)
    do j = 0, n1
        rho = rho_of_s(s_min) + j*(rho_of_s(s_max) - rho_of_s(s_min))/n1
        s_node(j) = (r0 - sqrt(r0**2 - rho**2))/s_edge
        pt_node(j) = psi_tor(rho)
        pp_node(j) = psi_pol(rho)
    end do

    nplane = ntetr/grid_size(2)
    err_bphi = 0.0_dp; err_bth = 0.0_dp; err_q = 0.0_dp; err_b = 0.0_dp; err_h = 0.0_dp
    nfail = 0
    do ind = 1, nplane
        do j = 1, 4
            iv = tetra_grid(ind)%ind_knot(j)
            if (coord_system == 2) then
                xv(:, j) = verts_sthetaphi(:, iv)
            else
                xv(:, j) = verts_rphiz(:, iv)
            end if
        end do
        if (coord_system == 2) then
            ! theta* seam: a vertex stored at 0 belongs to 2*pi in this tetrahedron
            if (any(xv(2, :) >= pi)) where (xv(2, :) == 0.0_dp) xv(2, :) = 2.0_dp*pi
        end if
        xc = sum(xv, dim=2)/4.0_dp

        if (coord_system == 2) then
            ! Contravariant B of the mesh: curl A / sqrt(g)
            bth_g = tetra_physics(ind)%curlA(2)/metric_determinant(ind, xc)
            bphi_g = tetra_physics(ind)%curlA(3)/metric_determinant(ind, xc)
            rho = rho_of_s(xc(1))
            thg = theta_geom(rho, xc(2))
            call contravariant_true(rho, thg, bth_true, bphi_true)
            err_bphi = max(err_bphi, abs(bphi_g/bphi_true - 1.0_dp))
            err_bth = max(err_bth, abs(bth_g/bth_true - 1.0_dp))
            ! Mesh q is the cell secant dpsi_tor/dpsi_pol (sqrt(g) cancels)
            smin_t = minval(xv(1, :)); smax_t = maxval(xv(1, :))
            q_g = tetra_physics(ind)%curlA(3)/tetra_physics(ind)%curlA(2)
            q_true = (node_flux(smax_t, pt_node) - node_flux(smin_t, pt_node)) &
                     /(node_flux(smax_t, pp_node) - node_flux(smin_t, pp_node))
            err_q = max(err_q, abs(q_g/q_true - 1.0_dp))
            ! |B| and covariant h at the first vertex (physical point = verts_rphiz)
            xr = verts_rphiz(:, tetra_grid(ind)%ind_knot(1))
            call covariant_true(xr, h_true, b_true)
            err_b = max(err_b, abs(tetra_physics(ind)%bmod1/b_true - 1.0_dp))
            err_h(1) = max(err_h(1), abs(tetra_physics(ind)%h1_1 - h_true(1)))
            err_h(2) = max(err_h(2), abs(tetra_physics(ind)%h2_1 - h_true(2))/rho)
            err_h(3) = max(err_h(3), abs(tetra_physics(ind)%h3_1 - h_true(3))/r0)
        else
            ! (R,phi,Z): sqrt(g) = R, B^phi = curlA(2)/R, B^Z = curlA(3)/R
            call field_analytic_circ(xc(1), xc(2), xc(3), br, bp, bz, d(1), d(2), d(3), &
                d(4), d(5), d(6), d(7), d(8), d(9))
            bphi_g = tetra_physics(ind)%curlA(2)/metric_determinant(ind, xc)
            err_bphi = max(err_bphi, abs(bphi_g/(bp/xc(1)) - 1.0_dp))
            bth_g = (-tetra_physics(ind)%curlA(1)*xc(3) + tetra_physics(ind)%curlA(3) &
                     *(xc(1) - r0))/xc(1)
            bth_true = -br*xc(3) + bz*(xc(1) - r0)
            err_bth = max(err_bth, abs(bth_g/bth_true - 1.0_dp))
        end if
    end do

    ! theta seam: a face shared with a neighbour must coincide after the periodic
    ! shift of the pusher, theta_own - theta_neighbour = 2*pi*neighbour_perbou_theta
    nseam = 0
    if (coord_system == 2) then
        do ind = 1, nplane
            do j = 1, 4
                nb = tetra_grid(ind)%neighbour_tetr(j)
                if (nb < 1) cycle
                th_own = face_theta(ind, j)
                th_nb = face_theta(nb, tetra_grid(ind)%neighbour_face(j))
                if (abs(th_own - th_nb - 2.0_dp*pi*tetra_grid(ind)%neighbour_perbou_theta(j)) &
                    > 1.0e-12_dp) nseam = nseam + 1
            end do
        end do
    end if

    print '(a,i2,a,f8.0)', 'coord_system =', coord_system, '  B0 =', b0
    print '(a,es10.3)', 'max rel err B^phi   (centroid) ', err_bphi
    print '(a,es10.3)', 'max rel err B^theta (centroid) ', err_bth
    if (coord_system == 2) then
        ! A_theta is linear in s: only the vertex interpolation of sqrt(g) remains
        call check(err_bphi < 5.0e-3_dp, 'B^phi')
        print '(a,es10.3)', 'max rel err q (cell secant)    ', err_q
        print '(a,es10.3)', 'max rel err |B| (vertex)       ', err_b
        print '(a,3es10.3)', 'max abs err h_s,h_th/rho,h_ph/R0', err_h
        call check(err_bth < 2.0e-2_dp, 'B^theta*')
        call check(err_q < 1.0e-5_dp, 'q secant')
        call check(err_b < 1.0e-10_dp, '|B| at vertices')
        call check(maxval(err_h) < 1.0e-6_dp, 'covariant h at vertices')
        print '(a,i0)', 'faces inconsistent across the theta seam ', nseam
        call check(nseam == 0, 'theta periodic boundary')
    else
        ! Linear A on the tetrahedron: centroid vs secant point, O(cell size)
        call check(err_bphi < 2.0e-2_dp, 'B^phi')
        call check(err_bth < 1.0e-1_dp, 'B_pol')
    end if
    if (nfail > 0) error stop 1
    print '(a)', 'PASS: grid_kind 5 mesh field agrees with the analytic field'

contains

    subroutine check(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            nfail = nfail + 1
            print '(2a)', 'FAIL: ', what
        end if
    end subroutine check

    ! Mean theta* of face iface (opposite to local vertex iface) of tetrahedron it, in the
    ! frame of it (seam vertices stored at 0 count as 2*pi, as in make_tetra_physics)
    real(dp) function face_theta(it, iface) result(th)
        integer, intent(in) :: it, iface
        real(dp) :: tv(4)
        tv = verts_sthetaphi(2, tetra_grid(it)%ind_knot)
        if (any(tv >= pi)) where (tv == 0.0_dp) tv = 2.0_dp*pi
        th = (sum(tv) - tv(iface))/3.0_dp
    end function face_theta

    ! Oracle flux at a mesh node, identified by its s value
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

    ! Field-line winding dphi/dtheta_geom on the circle rho
    real(dp) function winding(r, th) result(w)
        real(dp), intent(in) :: r, th
        real(dp) :: rr, zz, br_, bp_, bz_, dd(9)
        rr = r0 + r*cos(th); zz = r*sin(th)
        call field_analytic_circ(rr, 0.0_dp, zz, br_, bp_, bz_, dd(1), dd(2), dd(3), &
            dd(4), dd(5), dd(6), dd(7), dd(8), dd(9))
        w = (bp_/rr)/((-br_*sin(th) + bz_*cos(th))/r)
    end function winding

    ! Composite Simpson integral of the winding from 0 to th
    real(dp) function wint(r, th) result(acc)
        real(dp), intent(in) :: r, th
        real(dp) :: hh
        integer :: k
        hh = th/nquad
        acc = winding(r, 0.0_dp) + winding(r, th)
        do k = 1, nquad - 1
            acc = acc + merge(4.0_dp, 2.0_dp, mod(k, 2) == 1)*winding(r, k*hh)
        end do
        acc = acc*hh/3.0_dp
    end function wint

    ! Geometric angle at straight-field-line angle ts (Newton on the quadrature)
    real(dp) function theta_geom(r, ts) result(th)
        real(dp), intent(in) :: r, ts
        integer :: it
        wtot = wint(r, 2.0_dp*pi)
        th = ts
        do it = 1, 30
            ths = 2.0_dp*pi*wint(r, th)/wtot
            th = th - (ths - ts)/(2.0_dp*pi*winding(r, th)/wtot)
            if (abs(ths - ts) < 1.0e-13_dp) exit
        end do
    end function theta_geom

    real(dp) function theta_star(r, th) result(ts)
        real(dp), intent(in) :: r, th
        ts = 2.0_dp*pi*wint(r, th)/wint(r, 2.0_dp*pi)
    end function theta_star

    ! True B^theta* and B^phi at (rho, theta_geom)
    subroutine contravariant_true(r, th, bth, bph)
        real(dp), intent(in) :: r, th
        real(dp), intent(out) :: bth, bph
        real(dp) :: rr, zz, br_, bp_, bz_, dd(9), dts_dth
        rr = r0 + r*cos(th); zz = r*sin(th)
        call field_analytic_circ(rr, 0.0_dp, zz, br_, bp_, bz_, dd(1), dd(2), dd(3), &
            dd(4), dd(5), dd(6), dd(7), dd(8), dd(9))
        dts_dth = 2.0_dp*pi*winding(r, th)/wint(r, 2.0_dp*pi)
        bph = bp_/rr
        bth = dts_dth*(-br_*sin(th) + bz_*cos(th))/r
    end subroutine contravariant_true

    ! True covariant h = B_k/|B| in (s, theta*, phi) at the physical point x = (R,phi,Z)
    subroutine covariant_true(x, h, bm)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: h(3), bm
        real(dp) :: r, th, s, ts, br_, bp_, bz_, dd(9), dx(2, 2), dr, eps_s, th_p, th_m
        r = hypot(x(1) - r0, x(3))
        th = atan2(x(3), x(1) - r0)
        s = (r0 - sqrt(r0**2 - r**2))/s_edge
        ts = theta_star(r, modulo(th, 2.0_dp*pi))
        call field_analytic_circ(x(1), x(2), x(3), br_, bp_, bz_, dd(1), dd(2), dd(3), &
            dd(4), dd(5), dd(6), dd(7), dd(8), dd(9))
        bm = sqrt(br_**2 + bp_**2 + bz_**2)
        ! d(R,Z)/ds and d(R,Z)/dtheta* by central differences of the oracle map
        eps_s = 1.0e-5_dp
        th_p = theta_geom(rho_of_s(s + eps_s), ts)
        th_m = theta_geom(rho_of_s(s - eps_s), ts)
        dr = 0.5_dp/eps_s
        dx(1, 1) = dr*(rho_of_s(s + eps_s)*cos(th_p) - rho_of_s(s - eps_s)*cos(th_m))
        dx(2, 1) = dr*(rho_of_s(s + eps_s)*sin(th_p) - rho_of_s(s - eps_s)*sin(th_m))
        dts_dth_cov: block
            real(dp) :: dth_dts
            dth_dts = wint(r, 2.0_dp*pi)/(2.0_dp*pi*winding(r, th))
            dx(1, 2) = -r*sin(th)*dth_dts
            dx(2, 2) = r*cos(th)*dth_dts
        end block dts_dth_cov
        h(1) = (br_*dx(1, 1) + bz_*dx(2, 1))/bm
        h(2) = (br_*dx(1, 2) + bz_*dx(2, 2))/bm
        h(3) = x(1)*bp_/bm
    end subroutine covariant_true

    ! Toroidal flux per radian inside rho: int B_phi dA / (2 pi)
    real(dp) function psi_tor(r) result(acc)
        real(dp), intent(in) :: r
        integer, parameter :: nr = 200, nt = 256
        real(dp) :: rr, th, br_, bp_, bz_, dd(9), hh, wr
        integer :: k, l
        hh = r/nr
        acc = 0.0_dp
        do k = 1, nr  ! the k = 0 node (axis) carries zero weight rho*B_phi
            rr = k*hh
            wr = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == 0 .or. k == nr)
            do l = 1, nt
                th = 2.0_dp*pi*(l - 0.5_dp)/nt
                call field_analytic_circ(r0 + rr*cos(th), 0.0_dp, rr*sin(th), br_, bp_, &
                    bz_, dd(1), dd(2), dd(3), dd(4), dd(5), dd(6), dd(7), dd(8), dd(9))
                acc = acc + wr*hh/3.0_dp*rr*bp_/nt
            end do
        end do
    end function psi_tor

    ! Poloidal flux per radian inside rho: int_R0^{R0+rho} R B_Z dR on the midplane
    real(dp) function psi_pol(r) result(acc)
        real(dp), intent(in) :: r
        integer, parameter :: nr = 2000
        real(dp) :: rr, br_, bp_, bz_, dd(9), hh, wr
        integer :: k
        hh = r/nr
        acc = 0.0_dp
        do k = 1, nr  ! B_Z vanishes on the axis node k = 0
            rr = r0 + k*hh
            wr = merge(1.0_dp, merge(4.0_dp, 2.0_dp, mod(k, 2) == 1), k == 0 .or. k == nr)
            call field_analytic_circ(rr, 0.0_dp, 0.0_dp, br_, bp_, bz_, dd(1), dd(2), &
                dd(3), dd(4), dd(5), dd(6), dd(7), dd(8), dd(9))
            acc = acc + wr*hh/3.0_dp*rr*bz_
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

end program test_analytic_circular_sfl
