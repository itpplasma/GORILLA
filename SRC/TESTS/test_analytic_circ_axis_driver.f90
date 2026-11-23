! grid_kind = 5 mesh with the magnetic axis inside (sfc_s_min = 0), coord_system = 1:
! the innermost ring is closed by a triangle fan around one axis vertex per toroidal plane.
!   1. exactly 3*n2*n3 flat tetrahedra are flagged (tetra_degenerate), all others have volume > 0
!   2. neighbour links are symmetric, never point to a flagged tetrahedron, and open faces lie on rho = a
!   3. linked faces have opposite normals (the property GORILLA's normal-vector check tests)
!   4. find_tetra finds points on and next to the axis in a valid tetrahedron that contains them
!   5. a guiding-centre orbit that circles the axis inside the fan cells is not lost and conserves energy
program test_analytic_circ_axis_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use tetra_grid_mod, only: ntetr, tetra_grid, verts_rphiz, tetra_degenerate
    use tetra_physics_mod, only: tetra_physics, isinside, particle_mass
    use find_tetra_mod, only: find_tetra
    use supporting_functions_mod, only: bmod_func
    implicit none
    real(dp), parameter :: r0 = 200.0_dp, minor = 46.0_dp, b0 = 20000.0_dp
    integer, parameter :: n1 = 6, n2 = 8, n3 = 24
    real(dp) :: pi, xv(3, 4), vol, vol_min, x(3), vpar, vperp, rho, th, e0, e1, de_max, d, hr, ev
    integer :: k, f, nb, nf, n_bad_link, n_open_inner, n_normal, n_found, n_steps, n_ring1, iface, ind, j
    integer :: nfp
    character(len=16) :: arg
    logical :: init
    logical, allocatable :: visited(:)

    pi = acos(-1.0_dp)
    ev = 1.602176634e-12_dp
    nfp = 1      ! argument: number of field periods (periodic phi seam inside the mesh for nfp > 1)
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *) nfp
    end if
    call write_inputs()
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()
    if (.not. allocated(tetra_degenerate)) error stop 'tetra_degenerate not allocated: axis not closed'

    ! 1. flagged = zero volume, all others clearly positive (relative to the mean cell volume)
    vol_min = huge(1.0_dp)
    d = 0.0_dp
    do k = 1, ntetr
        xv = verts_rphiz(:, tetra_grid(k)%ind_knot)
        vol = abs(det3(xv(:, 2) - xv(:, 1), xv(:, 3) - xv(:, 1), xv(:, 4) - xv(:, 1)))
        if (tetra_degenerate(k)) then
            d = max(d, vol)
        else
            vol_min = min(vol_min, vol)
        end if
    end do
    vol = 2.0_dp*pi/nfp*minor**2/ntetr      ! mean tetrahedron volume in (R, phi, Z) coordinate space
    print '(a,i0,a,i0,a,es10.3,a,es10.3)', 'flat tetrahedra ', count(tetra_degenerate), ' of ', ntetr, &
        ', largest flagged volume ', d, ', smallest other volume / mean ', vol_min/vol
    if (count(tetra_degenerate) == 0) error stop 'no flat tetrahedra flagged'
    if (d > 0.0_dp) error stop 'a flagged tetrahedron has nonzero volume'
    if (vol_min < 1.0e-3_dp*vol) error stop 'zero-volume tetrahedron left in the mesh'
    ! total physical volume (Jacobian R, exact for linear R) = torus with the n3-gon cross-section
    vol = 0.0_dp
    do k = 1, ntetr
        if (tetra_degenerate(k)) cycle
        xv = verts_rphiz(:, tetra_grid(k)%ind_knot)
        vol = vol + abs(det3(xv(:, 2) - xv(:, 1), xv(:, 3) - xv(:, 1), xv(:, 4) - xv(:, 1)))/6.0_dp &
            *sum(xv(1, :))/4.0_dp
    end do
    d = 2.0_dp*pi/nfp*r0*0.5_dp*n3*minor**2*sin(2.0_dp*pi/n3)
    print '(a,es16.8,a,es16.8,a,es10.3)', 'mesh volume ', vol, ' polygon torus ', d, ' rel. diff ', vol/d - 1.0_dp
    if (abs(vol/d - 1.0_dp) > 1.0e-12_dp) error stop 'mesh volume differs (gap or overlap)'

    ! 2. and 3. links, open faces, normals
    n_bad_link = 0
    n_open_inner = 0
    n_normal = 0
    do k = 1, ntetr
        if (tetra_degenerate(k)) cycle
        do f = 1, 4
            nb = tetra_grid(k)%neighbour_tetr(f)
            if (nb < 1) then
                ! open face: its three vertices must lie on rho = a
                do j = 1, 4
                    if (j == f) cycle
                    rho = hypot(verts_rphiz(1, tetra_grid(k)%ind_knot(j)) - r0, verts_rphiz(3, tetra_grid(k)%ind_knot(j)))
                    if (abs(rho - minor) > 1.0e-9_dp*minor) then
                        n_open_inner = n_open_inner + 1
                        exit
                    end if
                end do
                cycle
            end if
            nf = tetra_grid(k)%neighbour_face(f)
            if (tetra_degenerate(nb)) then
                n_bad_link = n_bad_link + 1
            else if (tetra_grid(nb)%neighbour_tetr(nf) /= k .or. tetra_grid(nb)%neighbour_face(nf) /= f) then
                n_bad_link = n_bad_link + 1
            else if (dot_product(tetra_physics(k)%anorm(:, f), tetra_physics(nb)%anorm(:, nf)) >= 0.0_dp) then
                n_normal = n_normal + 1
            end if
        end do
    end do
    print '(a,i0,a,i0,a,i0)', 'bad links ', n_bad_link, ', open faces inside the mesh ', n_open_inner, &
        ', linked faces without opposite normals ', n_normal
    if (n_bad_link /= 0) error stop 'neighbour links not symmetric or pointing to a flat tetrahedron'
    if (n_open_inner /= 0) error stop 'open faces inside the mesh (axis not closed)'
    if (n_normal /= 0) error stop 'linked faces without opposite normals'

    ! 4. find_tetra on and next to the axis
    hr = minor/n1
    n_found = 0
    do j = 0, 199
        rho = merge(0.0_dp, hr*real(mod(j*37, 100), dp)/100.0_dp, j < 20)
        th = 2.0_dp*pi*real(mod(j*61, 97), dp)/97.0_dp
        x = [r0 + rho*cos(th), 2.0_dp*pi/nfp*real(mod(j*13, 50), dp)/50.0_dp, rho*sin(th)]
        vpar = 1.0e7_dp
        vperp = 1.0e7_dp
        call find_tetra(x, vpar, vperp, ind, iface)
        if (ind < 1) cycle
        if (tetra_degenerate(ind)) cycle
        if (.not. isinside(ind, x)) cycle
        n_found = n_found + 1
    end do
    print '(a,i0,a)', 'find_tetra: ', n_found, ' of 200 points on or next to the axis found'
    if (n_found /= 200) error stop 'find_tetra failed near the axis'

    ! 5. guiding-centre orbit around the axis (deuteron, 1 keV, passing, starting 3 mm off axis)
    x = [r0 + 0.3_dp, 0.1_dp, 0.0_dp]
    d = sqrt(2.0_dp*1000.0_dp*ev/particle_mass)
    vpar = 0.8_dp*d
    vperp = 0.6_dp*d
    init = .false.
    ind = -1
    allocate(visited(ntetr))
    visited = .false.
    n_ring1 = 0
    de_max = 0.0_dp
    e0 = -1.0_dp
    do n_steps = 1, 2000
        call orbit_timestep_gorilla(x, vpar, vperp, 1.0e-7_dp, init, ind, iface)
        if (ind < 1) then
            print '(a,i0,a,3es14.6)', 'orbit lost after ', n_steps, ' steps at ', x
            error stop 'guiding-centre orbit lost near the axis'
        end if
        if (tetra_degenerate(ind)) error stop 'orbit entered a flat tetrahedron'
        rho = hypot(x(1) - r0, x(3))
        if (rho < hr) n_ring1 = n_ring1 + 1
        visited(ind) = .true.
        e1 = 0.5_dp*particle_mass*(vpar**2 + vperp**2)
        if (e0 < 0.0_dp) e0 = e1
        de_max = max(de_max, abs(e1/e0 - 1.0_dp))
    end do
    print '(a,i0,a,i0,a,es10.3)', 'orbit: ', n_ring1, ' of 2000 steps in the fan ring, ', count(visited), &
        ' tetrahedra visited, max |dE/E| ', de_max
    if (n_ring1 /= 2000) error stop 'orbit left the innermost ring (test setup)'
    if (count(visited) < 3*n3) error stop 'orbit did not circle the axis'
    if (de_max > 1.0e-8_dp) error stop 'energy not conserved'
    print '(a)', 'PASS: grid_kind 5 axis closed by a triangle fan'

contains

    pure function det3(a, b, c) result(dd)
        real(dp), intent(in) :: a(3), b(3), c(3)
        real(dp) :: dd
        dd = a(1)*(b(2)*c(3) - b(3)*c(2)) - a(2)*(b(1)*c(3) - b(3)*c(1)) + a(3)*(b(1)*c(2) - b(2)*c(1))
    end function det3

    subroutine write_inputs()
        integer :: iu
        open(newunit=iu, file='tetra_grid.inp', status='replace')
        write(iu, '(a)') '&TETRA_GRID_NML'
        write(iu, '(a,i0,a,i0,a,i0)') ' grid_kind = 5, n1 = ', n1, ', n2 = ', n2, ', n3 = ', n3
        write(iu, '(a,i0,a)') ' boole_n_field_periods = .false., n_field_periods_manual = ', nfp, ','
        write(iu, '(a)') ' boole_write_mesh_obj = .false.,'
        write(iu, '(a)') ' sfc_s_min = 0.0, sfc_s_max = 1.0,'
        write(iu, '(a,f10.4,a,f10.4,a,f12.2,a)') ' R0_analytic_circ = ', r0, ', a_analytic_circ = ', minor, &
            ', B0_analytic_circ = ', b0, ','
        write(iu, '(a)') ' q0_analytic_circ = 1.1, q1_analytic_circ = 2.0,'
        write(iu, '(a)') " q_profile_file_analytic_circ = '',"
        write(iu, '(a)') " helical_pert_file_analytic_circ = '',"
        write(iu, '(a)') '/'
        close(iu)
        open(newunit=iu, file='gorilla.inp', status='replace')
        write(iu, '(a)') '&GORILLANML'
        write(iu, '(a)') ' coord_system = 1, ispecies = 2,'
        write(iu, '(a)') ' eps_Phi = 0.0, boole_periodic_relocation = .true.,'
        write(iu, '(a)') ' ipusher = 2, boole_pusher_ode45 = .false., rel_err_ode45 = 1e-10,'
        write(iu, '(a)') ' boole_dt_dtau = .true., boole_newton_precalc = .false.,'
        write(iu, '(a)') ' poly_order = 4, i_precomp = 0, boole_guess = .true.,'
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

end program test_analytic_circ_axis_driver
