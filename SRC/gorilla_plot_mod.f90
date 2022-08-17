!
module gorilla_plot_mod
!
    implicit none
!
    private
!
    !Input file variables - Definitions can be found in Input File:
    double precision,public,protected   :: total_orbit_time, start_pos_x1_beg, start_pos_x1_end, start_pitch_parameter, &
                                        & energy_eV_start
    logical,public,protected            :: boole_poincare_phi_0, boole_poincare_vpar_0, boole_full_orbit, &
                                        & boole_e_tot, boole_p_phi, boole_J_par
    integer,public,protected            :: n_skip_phi_0, n_skip_vpar_0, n_skip_full_orbit, n_surfaces, &
                                        & i_orbit_options
    character(50),public,protected      :: filename_poincare_phi_0_rphiz, filename_poincare_phi_0_sthetaphi, &
                                        & filename_poincare_vpar_0_rphiz, filename_poincare_vpar_0_sthetaphi, &
                                        & filename_full_orbit_rphiz, filename_full_orbit_sthetaphi, &
                                        & filename_orbit_start_pos_rphiz, filename_orbit_start_pos_sthetaphi, &
                                        & filename_e_tot, filename_p_phi, filename_J_par
    double precision,public,protected   :: start_pos_x2, start_pos_x3
!
    !Namelist for GORILLA Plotting input
    NAMELIST /gorilla_plot_nml/ total_orbit_time, start_pos_x1_beg, start_pos_x1_end,boole_poincare_phi_0, boole_poincare_vpar_0, &
    & boole_full_orbit, boole_e_tot, boole_p_phi, boole_J_par, n_skip_phi_0, n_skip_vpar_0, n_skip_full_orbit, n_surfaces, &
    & i_orbit_options, filename_poincare_phi_0_rphiz, filename_poincare_phi_0_sthetaphi, filename_poincare_vpar_0_rphiz, &
    & filename_poincare_vpar_0_sthetaphi, filename_full_orbit_rphiz, filename_full_orbit_sthetaphi, &
    & filename_orbit_start_pos_rphiz, filename_orbit_start_pos_sthetaphi, start_pos_x2, start_pos_x3, &
    & filename_e_tot, filename_p_phi, filename_J_par, start_pitch_parameter,energy_eV_start
!
    public              :: phi_elec_func, &
                        & gorilla_plot
!    
    contains
!
        subroutine gorilla_plot()
!
            use tetra_physics_mod, only: coord_system,particle_mass
            use constants, only: ev2erg
            use supporting_functions_mod, only: sym_flux_in_cyl
!
            implicit none
!
            integer                                      :: file_id_phi_0,file_id_vpar_0,file_id_full_orbit,file_id_e_tot, &
                                                         & file_id_p_phi,file_id_J_par,file_id_read_start,counter_phi_0_mappings, &
                                                         & counter_vpar_0_mappings, counter_tetrahedron_passes, ind_tetr, i_surf
            integer                                      :: i,counter_lost_particles,n_start_pos,i_os
            double precision                             :: vmod,vpar,vperp
            double precision,dimension(3)                :: x
            double precision,dimension(:,:),allocatable  :: start_pos_pitch_mat
!
            !Load plot settings from input file
            call load_gorilla_plot_inp()
!
            !---------------------------------------------------------------------------------------------------------------------!
            !File Handling
!
            !Select file id for Poincaré mappings depending on coordinate system
            select case(coord_system)
                case(1) !Cylindrical coordinates
                    file_id_phi_0 = 91
                    file_id_vpar_0 = 93
                    file_id_full_orbit = 95
                case(2) !Symmetry Flux coordinates
                    file_id_phi_0 = 92
                    file_id_vpar_0 = 94
                    file_id_full_orbit = 96
            end select
!
            !Determine remaining file ids
            file_id_e_tot = 97
            file_id_p_phi = 98
            file_id_J_par = 99
            file_id_read_start = 100
!
            !Open files (optional - only if input options require)
            if(boole_poincare_phi_0) then
                open(unit=91, file=filename_poincare_phi_0_rphiz, status='unknown')
                open(unit=92, file=filename_poincare_phi_0_sthetaphi, status='unknown')
            endif
!
            if(boole_poincare_vpar_0) then
                open(unit=93, file=filename_poincare_vpar_0_rphiz, status='unknown')
                open(unit=94, file=filename_poincare_vpar_0_sthetaphi, status='unknown')
            endif
!
            if(boole_full_orbit) then
                open(unit=95, file=filename_full_orbit_rphiz, status='unknown')
                open(unit=96, file=filename_full_orbit_sthetaphi, status='unknown')
            endif
!
            if(boole_e_tot) then
                open(unit=file_id_e_tot, file=filename_e_tot, status='unknown')
            endif
!
            if(boole_p_phi) then
                open(unit=file_id_p_phi, file=filename_p_phi, status='unknown')
            endif
!
            if(boole_J_par) then
                open(unit=file_id_J_par, file=filename_J_par, status='unknown')
            endif
!
            !---------------------------------------------------------------------------------------------------------------------!
            !Load starting positions and starting pitch parameters from specific file (i_orbit_options = 1 & 3)
            !In the case (i_orbit_options = 2 & 4) the starting positions and starting pitch parameters are taken from input file.
!
            select case(i_orbit_options)
                case(1)
                    allocate(start_pos_pitch_mat(1,4))
!
                    !Open file with starting positions and starting pitch parameter
                    select case(coord_system)
                        case(1) !Cylindrical coordinates
                            open(unit=file_id_read_start, file=filename_orbit_start_pos_rphiz, iostat=i_os, status='old')
                        case(2) !Symmetry flux coordinates
                            open(unit=file_id_read_start, file=filename_orbit_start_pos_sthetaphi, iostat=i_os, status='old')
                    end select
!
                    !Error, if file does not exist.
                    if ( i_os /= 0 ) then
                        select case(coord_system)
                            case(1) !Cylindrical coordinates
                                print *, "Error opening file with starting positions and starting pitch parameter: ", &
                                & filename_orbit_start_pos_rphiz
                            case(2) !Symmetry flux coordinates
                                print *, "Error opening file with starting positions and starting pitch parameter: ", &
                                & filename_orbit_start_pos_sthetaphi
                        end select
                        stop
                    endif
!
                    read(file_id_read_start,*) start_pos_pitch_mat(1,:)
!
                    close(file_id_read_start)
!
                case(3)
                    !Open file with starting positions and starting pitch parameter
                    select case(coord_system)
                        case(1) !Cylindrical coordinates
                            open(unit=file_id_read_start, file=filename_orbit_start_pos_rphiz, iostat=i_os, status='old')
                        case(2) !Symmetry flux coordinates
                            open(unit=file_id_read_start, file=filename_orbit_start_pos_sthetaphi, iostat=i_os, status='old')
                    end select
!
                    !Error, if file does not exist.
                    if ( i_os /= 0 ) then
                        select case(coord_system)
                            case(1) !Cylindrical coordinates
                                print *, "Error opening file with starting positions and starting pitch parameter: ", &
                                & filename_orbit_start_pos_rphiz
                            case(2) !Symmetry flux coordinates
                                print *, "Error opening file with starting positions and starting pitch parameter: ", &
                                & filename_orbit_start_pos_sthetaphi
                        end select
                        stop
                    endif
!
                    !Count number of lines
                    n_start_pos = 0
!
                    do
                        read(file_id_read_start, '(A)', iostat=i_os)
                        if (i_os /= 0) exit
                        n_start_pos = n_start_pos + 1
                    end do
!
                    print*, "File with starting positions and starting pitch parameter contains ", n_start_pos, "starting values."
!
                    allocate(start_pos_pitch_mat(n_start_pos,4))
!
                    rewind(file_id_read_start)
!
                    do i = 1, n_start_pos
                        read(file_id_read_start,*) start_pos_pitch_mat(i,:)
                    end do
!
                    close(file_id_read_start)
!
            end select

            !---------------------------------------------------------------------------------------------------------------------!
            !Orbit integration and plotting
!
            !Compute velocity module from kinetic energy dependent on particle species
            vmod=sqrt(2.d0*energy_eV_start*ev2erg/particle_mass)
!
            !Initialize Counter for lost particles
            counter_lost_particles = 0
!
            !Select orbit plotting options (Single OR Multiple orbits)
            select case(i_orbit_options)
!
                !Single orbit - Starting positions for the orbit are taken from file [First Line]
                case(1)
!
                    !Starting velocities
                    vpar = start_pos_pitch_mat(1,4)*vmod
                    vperp = sqrt(vmod**2-vpar**2)
!
                    !Starting positions from Input File
                    x(1) = start_pos_pitch_mat(1,1)  ! Cylindrical coordinates: R - Symmetry Flux coordinates: s
                    x(2) = start_pos_pitch_mat(1,2)  ! Cylindrical coordinates: $\varphi$ - Symmetry Flux coordinates: $\vartheta$
                    x(3) = start_pos_pitch_mat(1,3)  ! Cylindrical coordinates: Z - Symmetry Flux coordinates: $\varphi$
!
                    call gorilla_plot_orbit_integration(x,vpar,vperp,total_orbit_time,ind_tetr,file_id_phi_0, &
                                                        & file_id_vpar_0, file_id_e_tot,file_id_p_phi, &
                                                        & file_id_J_par, file_id_full_orbit, counter_phi_0_mappings, &
                                                        & counter_vpar_0_mappings, counter_tetrahedron_passes)
!
                    if(ind_tetr.eq.-1) then
                        !$omp critical
                        counter_lost_particles = counter_lost_particles +1
                        !$omp end critical
                    endif
!
                    print *, 'i_surf:', 1, 'x1_start:', start_pos_pitch_mat(1,1)
                    if(boole_poincare_phi_0) print *,'Number of toroidal mappings:', counter_phi_0_mappings
                    if(boole_poincare_vpar_0) print *,'Number of vpar=0 banana bounces:', counter_vpar_0_mappings
                    print *,'Number of tetrahedron pushings:', counter_tetrahedron_passes
!
                !Single orbit - Starting positions for the orbit are taken from starting drift surfaces [Input File]
                case(2)
!
                    !Starting velocities
                    vpar = start_pitch_parameter*vmod
                    vperp = sqrt(vmod**2-vpar**2)
!
                    !Starting positions from Input File
                    x(1) = start_pos_x1_beg ! Cylindrical coordinates: R - Symmetry Flux coordinates: s
                    x(2) = start_pos_x2     ! Cylindrical coordinates: $\varphi$ - Symmetry Flux coordinates: $\vartheta$
                    x(3) = start_pos_x3     ! Cylindrical coordinates: Z - Symmetry Flux coordinates: $\varphi$
!
                    call gorilla_plot_orbit_integration(x,vpar,vperp,total_orbit_time,ind_tetr,file_id_phi_0, &
                                                        & file_id_vpar_0, file_id_e_tot,file_id_p_phi, &
                                                        & file_id_J_par, file_id_full_orbit, counter_phi_0_mappings, &
                                                        & counter_vpar_0_mappings, counter_tetrahedron_passes)
!
                    if(ind_tetr.eq.-1) then
                        !$omp critical
                        counter_lost_particles = counter_lost_particles +1
                        !$omp end critical
                    endif
!
                    print *, 'i_surf:', 1, 'x1_start:', start_pos_x1_beg
                    if(boole_poincare_phi_0) print *,'Number of toroidal mappings:', counter_phi_0_mappings
                    if(boole_poincare_vpar_0) print *,'Number of vpar=0 banana bounces:', counter_vpar_0_mappings
                    print *,'Number of tetrahedron pushings:', counter_tetrahedron_passes
!
                !Multiple orbits - Starting positions for orbits are taken from file [Every Line New Starting position]
                case(3)
!
                    !Loop over drift surfaces from input file
!
                    !$OMP PARALLEL &
                    !$OMP& DO DEFAULT(NONE) &
                    !$OMP& SHARED(n_start_pos,total_orbit_time,file_id_phi_0,file_id_vpar_0,file_id_e_tot,file_id_p_phi, &
                    !$OMP& file_id_J_par, file_id_full_orbit, counter_lost_particles,start_pos_x1_beg,start_pos_x1_end, &
                    !$OMP& boole_poincare_phi_0,boole_poincare_vpar_0,vmod,start_pos_pitch_mat) &
                    !$OMP& PRIVATE(i_surf,ind_tetr,counter_phi_0_mappings,counter_vpar_0_mappings, counter_tetrahedron_passes, &
                    !$OMP& x,vpar,vperp)
                    do i_surf = 1,n_start_pos
!
                        !Starting velocities
                        vpar = start_pos_pitch_mat(i_surf,4)*vmod
                        vperp = sqrt(vmod**2-vpar**2)
!
                        !Starting positions from Input File
                        x(1) = start_pos_pitch_mat(i_surf,1)  ! Cylindrical coordinates: R - Symmetry Flux coordinates: s
                        x(2) = start_pos_pitch_mat(i_surf,2)  ! Cylindrical coordinates: $\varphi$ - Symmetry Flux coordinates: $\vartheta$
                        x(3) = start_pos_pitch_mat(i_surf,3)  ! Cylindrical coordinates: Z - Symmetry Flux coordinates: $\varphi$
!
                        call gorilla_plot_orbit_integration(x,vpar,vperp,total_orbit_time,ind_tetr,file_id_phi_0, &
                                                            & file_id_vpar_0, file_id_e_tot,file_id_p_phi, &
                                                            & file_id_J_par, file_id_full_orbit, counter_phi_0_mappings, &
                                                            & counter_vpar_0_mappings, counter_tetrahedron_passes)
!
                        if(ind_tetr.eq.-1) then
                            !$omp critical
                            counter_lost_particles = counter_lost_particles +1
                            !$omp end critical
                        endif
!
                        print *, 'i_surf:', i_surf, 'x1_start:', start_pos_pitch_mat(i_surf,1)
                        if(boole_poincare_phi_0) print *,'Number of toroidal mappings:', counter_phi_0_mappings
                        if(boole_poincare_vpar_0) print *,'Number of vpar=0 banana bounces:', counter_vpar_0_mappings
                        print *,'Number of tetrahedron pushings:', counter_tetrahedron_passes
!
                    enddo !i_surf
                    !$OMP END PARALLEL DO
!
                !Multiple orbits - Starting positions for orbits are taken from drift surfaces with regular spacing [Input File]
                case(4)
!
                    !Loop over surfaces from input file
!
                    !$OMP PARALLEL &
                    !$OMP& DO DEFAULT(NONE) &
                    !$OMP& SHARED(n_surfaces,total_orbit_time,file_id_phi_0,file_id_vpar_0,file_id_e_tot,file_id_p_phi, &
                    !$OMP& file_id_J_par, file_id_full_orbit, counter_lost_particles,start_pos_x1_beg,start_pos_x1_end, &
                    !$OMP& boole_poincare_phi_0,boole_poincare_vpar_0,start_pitch_parameter,vmod,start_pos_x2,start_pos_x3) &
                    !$OMP& PRIVATE(x,vpar,vperp, i_surf,ind_tetr,counter_phi_0_mappings,counter_vpar_0_mappings, &
                    !$OMP& counter_tetrahedron_passes)
                    do i_surf = 1,n_surfaces
!
                        !Starting velocities
                        vpar = start_pitch_parameter*vmod
                        vperp = sqrt(vmod**2-vpar**2)
!
                        !Starting positions from Input File
                        ! Cylindrical coordinates: R - Symmetry Flux coordinates: s
                        x(1) = start_pos_x1_beg + (start_pos_x1_end - start_pos_x1_beg) * (dble(i_surf-1)/(n_surfaces-1))
                        x(2) = start_pos_x2     ! Cylindrical coordinates: $\varphi$ - Symmetry Flux coordinates: $\vartheta$
                        x(3) = start_pos_x3     ! Cylindrical coordinates: Z - Symmetry Flux coordinates: $\varphi$
!
                        call gorilla_plot_orbit_integration(x,vpar,vperp,total_orbit_time,ind_tetr,file_id_phi_0, &
                                                            & file_id_vpar_0, file_id_e_tot,file_id_p_phi, &
                                                            & file_id_J_par, file_id_full_orbit, counter_phi_0_mappings, &
                                                            & counter_vpar_0_mappings, counter_tetrahedron_passes)
!
                        if(ind_tetr.eq.-1) then
                            !$omp critical
                            counter_lost_particles = counter_lost_particles +1
                            !$omp end critical
                        endif
!
                        print *, 'i_surf:', i_surf, 'x1_start:', &
                        & start_pos_x1_beg + (start_pos_x1_end - start_pos_x1_beg) * (dble(i_surf-1)/(n_surfaces-1))
                        if(boole_poincare_phi_0) print *,'Number of toroidal mappings:', counter_phi_0_mappings
                        if(boole_poincare_vpar_0) print *,'Number of vpar=0 banana bounces:', counter_vpar_0_mappings
                        print *,'Number of tetrahedron pushings:', counter_tetrahedron_passes
!
                    enddo !i_surf
                    !$OMP END PARALLEL DO
!
            end select
!
            !Print lost particles for all starting positions
            print *, ''
            print *,'Number of lost particles for all starting positions:', counter_lost_particles
!
            !---------------------------------------------------------------------------------------------------------------------!
            !Transformation of coordinates and File handling
!
            if(boole_poincare_phi_0) then
                close(unit=91)
                close(unit=92)
            endif
!
            if(boole_poincare_vpar_0) then
                close(unit=93)
                close(unit=94)
            endif
!
            if(boole_full_orbit) then
                close(unit=95)
                close(unit=96)
            endif
!
            if(boole_e_tot) then
                close(unit=file_id_e_tot)
            endif
!
            if(boole_p_phi) then
                close(unit=file_id_p_phi)
            endif
!
            if(boole_J_par) then
                close(unit=file_id_J_par)
            endif
!
            !If symmetry flux coordinates are used, transform to cylindrical coordinates
            if(coord_system.eq.2) then
                if(boole_poincare_phi_0) then
                    call sym_flux_in_cyl(filename_poincare_phi_0_sthetaphi,filename_poincare_phi_0_rphiz,0)
                endif
!
                if(boole_poincare_vpar_0) then
                    call sym_flux_in_cyl(filename_poincare_vpar_0_sthetaphi,filename_poincare_vpar_0_rphiz,0)
                endif
!
                if(boole_full_orbit) then
                    call sym_flux_in_cyl(filename_full_orbit_sthetaphi,filename_full_orbit_rphiz,0)
                endif
            endif
!
            !---------------------------------------------------------------------------------------------------------------------!
            !Deallocations
!
            select case(i_orbit_options)
                case(1,3)
                    deallocate(start_pos_pitch_mat)
            end select
!
        end subroutine gorilla_plot
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine load_gorilla_plot_inp()
!
            open(unit=90, file='gorilla_plot.inp', status='unknown')
            read(90,nml=gorilla_plot_nml)
            close(90)
!
            print *,'GORILLA Plotting: Loaded input data from gorilla_plot.inp'
!
        end subroutine load_gorilla_plot_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine gorilla_plot_orbit_integration(x,vpar,vperp,t_step,ind_tetr,file_id_phi_0,file_id_vpar_0,file_id_e_tot, &
                                                 & file_id_p_phi, file_id_J_par, file_id_full_orbit, counter_phi_0_mappings, &
                                                 & counter_vpar_0_mappings, counter_tetrahedron_passes)
!
            use pusher_tetra_rk_mod, only: find_tetra,pusher_tetra_rk,initialize_const_motion_rk
            use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly,manage_intermediate_steps_arrays
            use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
                &alloc_precomp_poly_perpinv
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            !use tetra_grid_mod, only: ntetr
            use par_adiab_inv_poly_mod, only: par_adiab_inv_tetra_poly,counter_banana_mappings_poly => counter_banana_mappings
            use par_adiab_inv_rk_mod, only: par_adiab_inv_tetra_rk,counter_banana_mappings_rk => counter_banana_mappings
            use gorilla_settings_mod, only: ipusher, poly_order
            use orbit_timestep_gorilla_mod, only: check_coordinate_domain
!
            implicit none
!
            double precision, dimension(3),intent(inout)    :: x
            double precision, intent(inout)                 :: vpar,vperp
            double precision, intent(in)                    :: t_step
            integer, intent(in)                             :: file_id_phi_0,file_id_vpar_0,file_id_e_tot,file_id_p_phi, &
                                                            & file_id_J_par,file_id_full_orbit
            integer, intent(out)                            :: ind_tetr, counter_phi_0_mappings,counter_vpar_0_mappings, &
                                                            & counter_tetrahedron_passes
            integer                                         :: iface
            double precision, dimension(3)                  :: z_save
            double precision                                :: vperp2,t_remain,t_pass,vpar_save
            logical                                         :: boole_t_finished
            integer                                         :: ind_tetr_save,file_id,iper,k
            double precision                                :: perpinv,perpinv2
!
            !Check coordinate domain (optionally perform modulo operation)
            call check_coordinate_domain(x)

            !Find tetrahedron index and face index for position x
            call find_tetra(x,vpar,vperp,ind_tetr,iface)
!               
            !If particle doesn't lie inside any tetrahedron
            if(ind_tetr.eq.-1) then
                print *, 'Error in GORILLA Plotting: Particle does not lie inside any tetrahedron.'
                stop
            endif
!           
            !Exit the subroutine after initialization, if time step equals zero
            if(t_step.eq.0.d0) then
                print *, 'Error in GORILLA Plotting: Time step equals zero.'
                stop
            endif
!
            !Initialize mapping counter
            counter_phi_0_mappings = 0
            counter_vpar_0_mappings = 0
            counter_tetrahedron_passes = 0
!
            !Squared perpendicular velocity
            vperp2 = vperp**2
!                
            !Compute relative particle position
            z_save = x-tetra_physics(ind_tetr)%x1
!
            !Compute perpendicular invariant of particle
            perpinv=-0.5d0*vperp2/bmod_func(z_save,ind_tetr)
            perpinv2 = perpinv**2
!               
            !Initialize constants of motion in particle-private module
            select case(ipusher)
                case(1)
                    call initialize_const_motion_rk(perpinv,perpinv2)
                case(2)
                    call initialize_const_motion_poly(perpinv,perpinv2)
                    call manage_intermediate_steps_arrays()
            end select        
!
!             !NOT FULLY IMPLEMENTED YET: Precompute quatities dependent on perpinv
!             call alloc_precomp_poly_perpinv(1,ntetr)
!             call initialize_boole_precomp_poly_perpinv()
!             call make_precomp_poly_perpinv(perpinv,perpinv2)
!
            !Integrate particle orbit for given time step
            t_remain = t_step
!
            !Logical for handling time integration
            boole_t_finished = .false.
!
            !Loop for tetrahedron pushings until t_step is reached
            do
!
                !Domain Boundary
                if(ind_tetr.eq.-1) then
                    print *, 'WARNING: Particle lost at x(1) = ',x(1), 'x(2)',x(2),'x(3)',x(3)
                    exit
                endif
!                
                !Save the tetrahedron index for computation of vperp in the last step
                ind_tetr_save = ind_tetr
!
                !Save vpar for the computation of the parallel adiabatic invariant
                vpar_save = vpar
!  
                !t_remain (in) ... remaining time until t_step is finished
                !t_pass (out) ... time to pass the tetrahdron
!
                !Calculate trajectory
                select case(ipusher)
                    case(1)
                        call pusher_tetra_rk(ind_tetr,iface,x,vpar,z_save,t_remain,t_pass,boole_t_finished,iper)
                    case(2)
                        call pusher_tetra_poly(poly_order,ind_tetr,iface,x,vpar,z_save,t_remain,&
                                                         & t_pass,boole_t_finished,iper)
                end select
!
                t_remain = t_remain - t_pass
                counter_tetrahedron_passes = counter_tetrahedron_passes + 1
!print *, 'counter_tetrahedron_passes', counter_tetrahedron_passes
!
                !Optional plotting of full orbit (after each tetrahedron passing)
                if(boole_full_orbit) then
                    if(counter_tetrahedron_passes/n_skip_full_orbit*n_skip_full_orbit.eq.counter_tetrahedron_passes) then
                        !Write elapsed time and orbit position
                        !$omp critical
                            write(file_id_full_orbit,*) x
                        !$omp end critical
    !
                        !Compute canonical toroidal angular momentum
                        if(boole_p_phi) then
                            !Write elapsed time and canonical toroidal angular momentum
                            !$omp critical
                                write(file_id_p_phi,*) t_step - t_remain, p_phi_func(vpar,z_save,ind_tetr_save)
                            !$omp end critical
                        endif
    !
                        !Computation of total energy
                        if(boole_e_tot) then
                            !Write elapsed time and total energy
                            !$omp critical
                                write(file_id_e_tot,*) t_step - t_remain, E_tot_func(z_save,vpar,perpinv,ind_tetr_save)
                            !$omp end critical
                        endif !boole_e_tot
                    endif
                endif
!
                !Orbit stops within cell, because "flight"-time t_step has finished
                if(boole_t_finished) then              
                    exit
                endif
!
                !Optional computation of parallel adiabatic invariant $J_\parallel$
                if(boole_J_par.or.boole_poincare_vpar_0) then
                     select case(ipusher)
                         case(1)
                             call par_adiab_inv_tetra_rk(t_pass,vpar_save,vpar,file_id_vpar_0,file_id_J_par,n_skip_vpar_0, &
                                                        & boole_J_par,boole_poincare_vpar_0,boole_e_tot,file_id_e_tot)
                         case(2)
                             call par_adiab_inv_tetra_poly(poly_order,t_pass,vpar_save,vpar,file_id_vpar_0,file_id_J_par, &
                                                          & n_skip_vpar_0,boole_J_par,boole_poincare_vpar_0,boole_e_tot, &
                                                          & file_id_e_tot)
                     end select
                endif !boole_J_par
!
                !Toroidal mappings at $\varphi$ = 0
                !Option to plot Poincaré mappings
                !Option to compute canonical toroidal angular momentum
                if(iper.ne.0) then
!
                    !Toroidal mapping counter
                    counter_phi_0_mappings = counter_phi_0_mappings + iper
!
                    !Optional plot of poincare cuts
                    if(boole_poincare_phi_0) then
                        if(counter_phi_0_mappings/n_skip_phi_0*n_skip_phi_0.eq.counter_phi_0_mappings) then
!
                            !Write coordinates for Poincaré mappings at $\varphi$ = 0
                            !$omp critical
                                write(file_id_phi_0,*) x
                            !$omp end critical
                        endif
                    endif
!
                    !Compute canonical toroidal angular momentum
                    if(boole_p_phi.and.boole_poincare_phi_0) then
                        if(counter_phi_0_mappings/n_skip_phi_0*n_skip_phi_0.eq.counter_phi_0_mappings) then
!
                            !Write number of toroidal mappings and canonical toroidal angular momentum
                            !$omp critical
                                write(file_id_p_phi,*) counter_phi_0_mappings, p_phi_func(vpar,z_save,ind_tetr_save)
                            !$omp end critical
                        endif
                    endif
!
                    !Computation of total energy
                    if(boole_e_tot.and.boole_poincare_phi_0) then
                        if(counter_phi_0_mappings/n_skip_phi_0*n_skip_phi_0.eq.counter_phi_0_mappings) then
!
                            !Write number of toroidal mappings and total energy
                            !$omp critical
                                write(file_id_e_tot,*) counter_phi_0_mappings, E_tot_func(z_save,vpar,perpinv,ind_tetr_save)
                            !$omp end critical
                        endif
                    endif !boole_e_tot
!
                endif ! iper.ne.0
!
            enddo !Loop for tetrahedron pushings
!
            !Deallocate intermediate_steps_arrays if pusher_poly was used
            call manage_intermediate_steps_arrays()
!
            !Compute vperp from position
            vperp = vperp_func(z_save,perpinv,ind_tetr_save)
!
            !Read counter of banana mappings:
            if(boole_J_par.or.boole_poincare_vpar_0) then
                 select case(ipusher)
                     case(1)
                         counter_vpar_0_mappings = counter_banana_mappings_rk
                     case(2)
                         counter_vpar_0_mappings = counter_banana_mappings_poly
                 end select
            endif !boole_J_par
!            
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)
!         
        end subroutine gorilla_plot_orbit_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function bmod_func(z123,ind_tetr)
!
            use tetra_physics_mod, only: tetra_physics
!
            implicit none
!
            double precision :: bmod_func
            integer, intent(in) :: ind_tetr
            double precision, dimension(3),intent(in) :: z123
!
            bmod_func = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z123)
!        
        end function bmod_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function vperp_func(z123,perpinv,ind_tetr)
!        
            use tetra_physics_mod, only: tetra_physics
!    
            implicit none
!
            double precision :: vperp_func
            integer, intent(in) :: ind_tetr
            double precision, dimension(3),intent(in) :: z123
            double precision, intent(in) :: perpinv

                if(perpinv.ne.0.d0) then
                    vperp_func=sqrt(2.d0*abs(perpinv)*bmod_func(z123,ind_tetr))
                else
                    vperp_func = 0.d0
                endif
!
        end function vperp_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function E_tot_func(z123,vpar,perpinv,ind_tetr)
!
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
!
            implicit none
!
            double precision :: E_tot_func
            integer, intent(in) :: ind_tetr
            double precision, dimension(3),intent(in) :: z123
            double precision, intent(in) :: perpinv,vpar
            double precision :: vperp
!
                vperp = vperp_func(z123,perpinv,ind_tetr)
!
                E_tot_func = 0.5d0*particle_mass*(vpar**2+vperp**2) + &
                            & particle_charge*phi_elec_func(z123,ind_tetr)
!
        end function E_tot_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function phi_elec_func(z123,ind_tetr)
!
            use tetra_physics_mod, only: tetra_physics
!
            implicit none
!
            double precision :: phi_elec_func
            integer, intent(in) :: ind_tetr
            double precision, dimension(3),intent(in) :: z123
!
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z123)
!        
        end function phi_elec_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function pitchpar_func(vpar,z,ind_tetr,perpinv)
!
            use tetra_physics_mod, only: tetra_physics
!
            implicit none
!
            integer :: ind_tetr
            double precision :: vmod,vpar,vperp,perpinv,pitchpar_func
            double precision,dimension(3) :: z
!
            vperp = vperp_func(z,perpinv,ind_tetr)
!
            vmod = sqrt(vpar**2+vperp**2)
            pitchpar_func = vpar/vmod
!        
        end function pitchpar_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        function p_phi_func(vpar,z,ind_tetr)
!
            use tetra_physics_mod, only: tetra_physics,particle_mass,cm_over_e,coord_system
!
            implicit none
!
            integer :: ind_tetr
            double precision :: vpar,p_phi_func,hphi1
            double precision,dimension(3) :: z,ghphi
!
            select case(coord_system)
                case(1)
                    hphi1 = tetra_physics(ind_tetr)%h2_1
                    ghphi = tetra_physics(ind_tetr)%gh2
                case(2)
                    hphi1 = tetra_physics(ind_tetr)%h3_1
                    ghphi = tetra_physics(ind_tetr)%gh3
            end select    
!
            p_phi_func = particle_mass*vpar*(hphi1+sum(ghphi*z(1:3))) + &
                        &particle_mass/cm_over_e* &
                        & (tetra_physics(ind_tetr)%Aphi1+sum(tetra_physics(ind_tetr)%gAphi*z(1:3)))
!        
        end function p_phi_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module gorilla_plot_mod









