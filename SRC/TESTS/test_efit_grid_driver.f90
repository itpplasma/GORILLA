! grid_kind 2 on the AUG test equilibrium (MHD_EQUILIBRIA/g_file_for_test):
!   fold_cs1      coord_system 1 + theta_geom_flux 1 folds near the separatrix -> run stops with a hint
!   fold_cs2      coord_system 2 + theta_geom_flux 2 folds                     -> run stops with a hint
!   clustering    theta_xpoint_clustering: no folds, much finer spacing of the last ring at the X-point;
!                 magnetic axis (coord_system 1) from the g-file
!   wall_missing  boole_extend_to_wall without a mesh                           -> run stops, prints the command
!   wall          boole_extend_to_wall with a mesh built by make_edge_mesh.sh: core flags, valid faces,
!                 magnetic axis from the g-file
!   wall_mismatch mesh built for other settings                                -> run stops, lists the difference
! The stopping cases are checked by ctest through the printed message (PASS_REGULAR_EXPRESSION).
! Usage: test_efit_grid_driver.x <case> <MHD_EQUILIBRIA directory>
program test_efit_grid_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use orbit_timestep_gorilla_mod, only: initialize_gorilla
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp, n2
    use tetra_grid_mod, only: ntetr, nvert, tetra_grid, verts_rphiz, verts_sthetaphi, &
                              verts_s_theta_core, tetra_in_core
    use tetra_physics_mod, only: mag_axis_R0, mag_axis_Z0
    implicit none
    character(len=32) :: tcase
    character(len=512) :: eqdir
    integer, parameter :: n1 = 30, n3 = 60
    real(dp), parameter :: s_min = 0.1_dp
    character(len=*), parameter :: wall_options = '--fill-axis --xref-radius 12'

    call get_command_argument(1, tcase)
    call get_command_argument(2, eqdir)
    call copy_file(trim(eqdir)//'/g_file_for_test', 'g_file_for_test')
    call copy_file(trim(eqdir)//'/convex_wall_for_test.dat', 'convex_wall_for_test.dat')

    select case (trim(tcase))
    case ('fold_cs1')
        call write_inputs(1, 1, 0.0_dp, .false., n3)
    case ('fold_cs2')
        call write_inputs(2, 2, 0.0_dp, .false., n3)
    case ('clustering')
        call write_inputs(1, 2, 0.9_dp, .false., n3)
    case ('wall_missing')
        call write_inputs(1, 2, 0.0_dp, .true., n3, 'no_such_mesh')
    case ('wall')
        call write_inputs(1, 2, 0.0_dp, .true., n3, 'mesh')
    case ('wall_mismatch')
        call write_inputs(1, 2, 0.0_dp, .true., 2*n3, 'mesh')
    case default
        print *, 'unknown case ', trim(tcase)
        error stop 1
    end select

    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()

    select case (trim(tcase))
    case ('clustering')
        call check_clustering()
    case ('wall')
        call check_wall()
    case default
        print *, 'FAIL: case ', trim(tcase), ' was expected to stop'
        error stop 1
    end select

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

    subroutine check_wall()
        integer :: iu, nv2, ntri, i, j, ncore_tri, nbad, region
        real(dp) :: st(2)
        ! expected core tetrahedra from the generator's region list
        open(newunit=iu, file='mesh/edge_mesh_core.dat', status='old', action='read')
        read(iu, *) nv2, ntri
        do i = 1, nv2
            read(iu, *) st
        end do
        ncore_tri = 0
        do i = 1, ntri
            read(iu, *) region
            if (region == 0) ncore_tri = ncore_tri + 1
        end do
        close(iu)
        if (count(tetra_in_core) /= 3*ncore_tri*n2) then
            print *, 'FAIL: core tetrahedra', count(tetra_in_core), ' expected', 3*ncore_tri*n2
            error stop 1
        end if
        do i = 1, ntetr
            if (tetra_in_core(i)) then
                if (any(verts_s_theta_core(1, tetra_grid(i)%ind_knot) < s_min - 1.0e-12_dp) .or. &
                    any(verts_s_theta_core(1, tetra_grid(i)%ind_knot) > 1.0_dp + 1.0e-12_dp)) then
                    print *, 'FAIL: core tetrahedron', i, ' has a vertex without s in [s_min, 1]'
                    error stop 1
                end if
            end if
        end do
        ! faces disabled by check_tetra_overlaps (overlapping neighbours)
        nbad = 0
        do i = 1, ntetr
            do j = 1, 4
                if (tetra_grid(i)%neighbour_tetr(j) > 0 .and. tetra_grid(i)%neighbour_face(j) == -1) nbad = nbad + 1
            end do
        end do
        if (nbad /= 0) then
            print *, 'FAIL: overlapping faces', nbad
            error stop 1
        end if
        if (abs(mag_axis_R0 - 170.8519220_dp) > 1.0e-6_dp .or. abs(mag_axis_Z0 - 8.922428638_dp) > 1.0e-6_dp) then
            print *, 'FAIL: magnetic axis', mag_axis_R0, mag_axis_Z0
            error stop 1
        end if
        print *, 'core tetrahedra', count(tetra_in_core), ' of', ntetr
        print *, 'PASS'
    end subroutine check_wall

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

    subroutine write_inputs(coord_system, theta_geom_flux, clustering, extend, n3_run, mesh_dir)
        integer, intent(in) :: coord_system, theta_geom_flux, n3_run
        real(dp), intent(in) :: clustering
        logical, intent(in) :: extend
        character(*), intent(in), optional :: mesh_dir
        integer :: iu
        open(newunit=iu, file='tetra_grid.inp', status='replace')
        write(iu, '(a)') '&TETRA_GRID_NML'
        write(iu, '(a,i0,a,i0,a)') ' grid_kind = 2, n1 = ', n1, ', n2 = 4, n3 = ', n3_run, ','
        write(iu, '(a,es12.4,a)') ' sfc_s_min = ', s_min, ','
        write(iu, '(a,i0,a)') ' theta_geom_flux = ', theta_geom_flux, ', theta0_at_xpoint = .true.,'
        if (clustering > 0.0_dp) write(iu, '(a,f6.3,a)') ' theta_xpoint_clustering = ', clustering, &
            ', theta_xpoint_clustering_width = 0.1,'
        if (extend) then
            write(iu, '(a)') ' boole_extend_to_wall = .true.,'
            write(iu, '(a)') " edge_mesh_dir = '"//mesh_dir//"',"
            write(iu, '(a)') " edge_mesh_options = '"//wall_options//"',"
        end if
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
