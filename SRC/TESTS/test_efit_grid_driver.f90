! grid_kind 2 on the AUG test equilibrium (MHD_EQUILIBRIA/g_file_for_test) in cylindrical coordinates
! with theta_geom_flux = 2 and theta_xpoint_clustering = 0.9:
!   - the mesh does not fold (GORILLA's overlap check would stop the run otherwise)
!   - the last flux ring is much finer at the X-point than on average
!   - the magnetic axis is taken from the g-file header
! Usage: test_efit_grid_driver.x <MHD_EQUILIBRIA directory>
program test_efit_grid_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, n2
    use tetra_grid_mod, only: nvert, verts_rphiz, verts_sthetaphi
    use tetra_physics_mod, only: mag_axis_R0, mag_axis_Z0
    implicit none
    character(len=512) :: eqdir
    integer, parameter :: n1 = 30, n3 = 60
    real(dp), parameter :: s_min = 0.1_dp

    call get_command_argument(1, eqdir)
    call copy_file(trim(eqdir)//'/g_file_for_test', 'g_file_for_test')
    call copy_file(trim(eqdir)//'/convex_wall_for_test.dat', 'convex_wall_for_test.dat')
    call write_inputs(1, 2, 0.9_dp)

    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()
    call check_clustering()

contains

    subroutine check_clustering()
        ! last ring of the first slice: spacing next to theta = 0 (X-point) against the mean spacing
        integer :: nv, i, i0, nring
        integer, allocatable :: idx(:)
        real(dp) :: smax, d_x, d_mean, pi
        pi = acos(-1.0_dp)
        nv = nvert/n2
        smax = maxval(verts_sthetaphi(1, 1:nv))
        idx = pack([(i, i = 1, nv)], abs(verts_sthetaphi(1, 1:nv) - smax) < 1.0e-12_dp)
        nring = size(idx)
        ! ring vertices are stored in theta order; i0 = vertex closest to theta = 0
        i0 = minloc(min(modulo(verts_sthetaphi(2, idx), 2*pi), 2*pi - modulo(verts_sthetaphi(2, idx), 2*pi)), 1)
        d_mean = 0.0_dp
        do i = 1, nring
            d_mean = d_mean + dist(idx(i), idx(modulo(i, nring) + 1))
        end do
        d_mean = d_mean/nring
        d_x = 0.5_dp*(dist(idx(i0), idx(modulo(i0, nring) + 1)) + dist(idx(i0), idx(modulo(i0 - 2, nring) + 1)))
        print '(a,f8.3,a,f8.3,a)', ' last ring: spacing at the X-point ', d_x, ' cm, mean ', d_mean, ' cm'
        if (abs(mag_axis_R0 - 170.8519220_dp) > 1.0e-6_dp .or. abs(mag_axis_Z0 - 8.922428638_dp) > 1.0e-6_dp) then
            print *, 'FAIL: magnetic axis', mag_axis_R0, mag_axis_Z0
            error stop 1
        end if
        if (d_x > 0.35_dp*d_mean) then
            print *, 'FAIL: clustering did not refine the X-point region'
            error stop 1
        end if
        print *, 'PASS'
    end subroutine check_clustering

    real(dp) function dist(i, j)
        integer, intent(in) :: i, j
        dist = hypot(verts_rphiz(1, i) - verts_rphiz(1, j), verts_rphiz(3, i) - verts_rphiz(3, j))
    end function dist

    subroutine copy_file(src, dst)
        character(*), intent(in) :: src, dst
        character(len=1024) :: line
        integer :: iin, iout, ios
        open(newunit=iin, file=src, status='old', action='read')
        open(newunit=iout, file=dst, status='replace', action='write')
        do
            read(iin, '(a)', iostat=ios) line
            if (ios /= 0) exit
            write(iout, '(a)') trim(line)
        end do
        close(iin)
        close(iout)
    end subroutine copy_file

    subroutine write_inputs(coord_system, theta_geom_flux, clustering)
        integer, intent(in) :: coord_system, theta_geom_flux
        real(dp), intent(in) :: clustering
        integer :: iu
        open(newunit=iu, file='tetra_grid.inp', status='replace')
        write(iu, '(a)') '&TETRA_GRID_NML'
        write(iu, '(a,i0,a,i0,a)') ' grid_kind = 2, n1 = ', n1, ', n2 = 4, n3 = ', n3, ','
        write(iu, '(a,es12.4,a)') ' sfc_s_min = ', s_min, ','
        write(iu, '(a,i0,a)') ' theta_geom_flux = ', theta_geom_flux, ', theta0_at_xpoint = .true.,'
        if (clustering > 0.0_dp) write(iu, '(a,f6.3,a)') ' theta_xpoint_clustering = ', clustering, &
            ', theta_xpoint_clustering_width = 0.1,'
        write(iu, '(a)') " g_file_filename = 'g_file_for_test',"
        write(iu, '(a)') " convex_wall_filename = 'convex_wall_for_test.dat',"
        write(iu, '(a)') " netcdf_filename = 'none.nc',"
        write(iu, '(a)') ' boole_n_field_periods = .true., n_field_periods_manual = 1,'
        write(iu, '(a)') " boole_write_mesh_obj = .false., filename_mesh_rphiz = 'm.obj', filename_mesh_sthetaphi = 'm2.obj',"
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
        write(iu, '(a)') '72'
        write(iu, '(a)') '0.99'
        write(iu, '(a)') '4'
        write(iu, '(a)') "'see_tetra_grid_inp'"
        write(iu, '(a)') "'none'"
        write(iu, '(a)') "'see_tetra_grid_inp'"
        write(iu, '(a)') "'none'"
        write(iu, '(a)') '0'
        write(iu, '(a)') '0'
        write(iu, '(a)') '0'
        close(iu)
        open(newunit=iu, file='preload_for_SYNCH.inp', status='replace')
        write(iu, '(a)') '3600      nstep'
        write(iu, '(a)') '500       nlabel'
        write(iu, '(a)') '500       ntheta'
        write(iu, '(a)') '10000     nsurfmax'
        close(iu)
    end subroutine write_inputs

end program test_efit_grid_driver
