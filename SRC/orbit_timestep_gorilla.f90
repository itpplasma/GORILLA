!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module orbit_timestep_gorilla_mod
!
    implicit none
!
    private
!
    public :: orbit_timestep_gorilla,initialize_gorilla,check_coordinate_domain
!
    integer, dimension(:,:), allocatable, public, protected   :: equidistant_grid
    integer, dimension(:), allocatable, public, protected     :: entry_counter
    double precision, dimension(3,5), public, protected       :: dimension_parametres
    logical, public, protected                                :: boole_axi_symmetry
!    
    contains
!
        subroutine orbit_timestep_gorilla(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface,t_remain_out)
!
            use supporting_functions_mod, only: bmod_func, vperp_func
            use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
            use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
            use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
                & alloc_precomp_poly_perpinv
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            !use tetra_grid_mod, only: ntetr
            use gorilla_settings_mod, only: ipusher, poly_order
            use find_tetra_mod, only: find_tetra
!
            implicit none
!
            double precision, dimension(3),intent(inout)    :: x
            double precision, intent(inout)                 :: vpar,vperp
            double precision, intent(in)                    :: t_step
            logical, intent(inout)                          :: boole_initialized
            integer, intent(inout)                          :: ind_tetr,iface
            double precision, intent(out), optional         :: t_remain_out
            double precision, dimension(3)                  :: z_save
            double precision                                :: vperp2,t_remain,t_pass,vpar_save
            logical                                         :: boole_t_finished
            integer                                         :: ind_tetr_save,iper,k
            double precision                                :: perpinv,perpinv2
            integer                                         :: sign_t_step
!
            !If orbit_timestep is called for the first time without grid position
            if(.not.boole_initialized) then
!
                !Check coordinate domain (optionally perform modulo operation)
                call check_coordinate_domain(x)
!
                !Compute sign t_step
                sign_t_step = int(sign(1.d0,t_step))
!
                !Find tetrahedron index and face index for position x
                call find_tetra(x,vpar,vperp,ind_tetr,iface,sign_t_step)
!               
                !If particle doesn't lie inside any tetrahedron
                if(ind_tetr.eq.-1) then
                    return
                endif
!
                boole_initialized = .true.
            endif
!           
            !Exit the subroutine after initialization, if time step equals zero
            if(t_step.eq.0.d0) return
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
                    print *, 'WARNING: Particle lost.'
                    if( present(t_remain_out)) then
                        t_remain_out = t_remain
                    endif
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
!
                !Orbit stops within cell, because "flight"-time t_step has finished
                if(boole_t_finished) then
                    if( present(t_remain_out)) then
                        t_remain_out = t_remain
                    endif
                    exit
                endif
!
            enddo !Loop for tetrahedron pushings
!
            !Compute vperp from position
            vperp = vperp_func(z_save,perpinv,ind_tetr_save)
!            
!             !NOT FULLY IMPLEMENTED YET: Deallocate precomputed quantities dependent on perpinv
!             call alloc_precomp_poly_perpinv(2,ntetr)
!         
        end subroutine orbit_timestep_gorilla
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_gorilla(i_option,ipert_in,bmod_multiplier)!,boole_grid_for_find_tetra)
!
            use constants, only: echarge,ame,amp,clight
            use gorilla_settings_mod, only: eps_Phi,coord_system,ispecies,boole_grid_for_find_tetra
            use tetra_grid_mod, only: make_tetra_grid
            use tetra_physics_mod, only: make_tetra_physics,check_tetra_overlaps,cm_over_e,particle_charge,particle_mass
            use tetra_physics_poly_precomp_mod, only: make_precomp_poly
            use tetra_grid_settings_mod, only: grid_kind
            use find_tetra_mod, only: grid_for_find_tetra

!
            implicit none
!
            integer                                 :: iper,ipert
            !logical, intent(inout), optional        :: boole_grid_for_find_tetra
            integer,intent(in),optional             :: i_option,ipert_in
            logical                                 :: boole_grid,boole_physics,boole_bmod_multiplier
            double precision,intent(in),optional    :: bmod_multiplier
!
            if(present(i_option)) then
                select case(i_option)
                    case(1)
                        boole_grid = .true.
                        boole_physics = .false.
                    case(2)
                        boole_grid = .false.
                        boole_physics = .true.
                    case(3)
                        boole_grid = .true.
                        boole_physics = .true.
                end select
            else
                boole_grid = .true.
                boole_physics = .true.
            endif
!
            if(present(ipert_in)) then
                ipert = ipert_in
            else
                ipert = 0
            endif
!
            if(present(bmod_multiplier)) then
                boole_bmod_multiplier = .true.
            else
                boole_bmod_multiplier = .false.
            endif
!
            if(boole_grid) then
!  
                select case(ispecies)
                    case(1) !electron
                        !electric charge of particle
                        particle_charge = echarge
!                   
                        !mass of particle
                        particle_mass = ame
!                    
                        !mass charge ratio
                        cm_over_e = clight*ame/echarge
!      
                    case(2) !deuterium ion
                        !electric charge of particle
                        particle_charge = echarge
!                   
                        !mass of particle
                        particle_mass = 2.d0*amp
!                    
                        !mass charge ratio
                        cm_over_e=2.d0*clight*amp/echarge
!      
                    case(3) !alpha particle 
                        !electric charge of particle
                        particle_charge = 2.d0*echarge
!                   
                        !mass of particle
                        particle_mass = 4.d0*amp
!                    
                        !mass charge ratio
                        cm_over_e=2.d0*clight*amp/echarge
!
                    case(4) !ionised tungsten 
                        !electric charge of particle
                        particle_charge = 74.d0*echarge
!                   
                        !mass of particle
                        particle_mass = 184.d0*amp
!                    
                        !mass charge ratio
                        cm_over_e= 184.d0*clight*amp/(74.d0*echarge)
!
                    case default
                        print *, 'ERROR: Invalid ispecies option selected. Check gorilla.inp.'
                        stop
!
                end select  
!  
                print *, 'Start grid computation'
                call make_tetra_grid()
!
            endif !boole_grid
!                
            if(boole_physics) then   
!
                if(boole_bmod_multiplier) then
                    call make_tetra_physics(coord_system,ipert,bmod_multiplier)
                else
                    call make_tetra_physics(coord_system,ipert)
                endif
                print *, 'Physics calculation of mesh is finished'
!
                print *, 'Start check tetra_overlaps'
                call check_tetra_overlaps
!
                call make_precomp_poly()
!
                !if (present(boole_grid_for_find_tetra)) then
                    !if (grid_kind.eq.1) boole_grid_for_find_tetra = .false. !rectangular grid
                    if (boole_grid_for_find_tetra) call grid_for_find_tetra
                !endif
!
            endif !boole_physics
!
        end subroutine initialize_gorilla
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_coordinate_domain(x)
!
            !This subroutine evaluates, if the given particle position is inside the computation domain.
            !In the case of a periodic coordinate, a modulo operation can optionally be performed.
!
            use gorilla_settings_mod, only: boole_periodic_relocation
            use tetra_physics_mod, only: coord_system
            use tetra_grid_settings_mod, only: sfc_s_min,n_field_periods
            use constants, only: pi
!
            implicit none
!
            double precision, dimension(3),intent(inout) :: x
!
            select case(coord_system)
                case(1) !Cylindrical coordinates
!
                    !Computation domain of periodic coordinate $\varphi$
                    if(boole_periodic_relocation) then
                        x(2) = modulo(x(2), (2.d0*pi/n_field_periods) ) !$\varphi$
                    else
                        if( x(2).lt.0.d0 ) then
                            print *, 'Error: Particle coordinate $\varphi$ is smaller than domain minimum 0.d0.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                        if( x(2).gt.(2.d0*pi/n_field_periods) ) then
                            print *, 'Error: Particle coordinate $\varphi$ is larger than domain maximum 2.d0*pi/n_field_periods.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                    endif

                case(2) !Symmetry Flux coordinates
                    !Computation domain of s
                    if( x(1).lt.sfc_s_min ) then
                        print *, 'Error: Particle flux coordinate s is smaller than GRID annulus s_min = ',sfc_s_min
                        stop
                    endif
!
                    if( x(1).gt.1.d0 ) then
                        print *, 'Error: Particle flux coordinate s is larger than GRID maximum s_max = 1.d0'
                        stop
                    endif
!
                    !Computation domain of periodic coordinates $\vartheta$ and $\varphi$
                    if(boole_periodic_relocation) then
                        x(2) = modulo(x(2),2.d0*pi) !$\vartheta$
                        x(3) = modulo(x(3), (2.d0*pi/n_field_periods) ) !$\varphi$
                    else
                        if( x(2).lt.0.d0 ) then
                            print *, 'Error: Particle coordinate $\vartheta$ is smaller than domain minimum 0.d0.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                        if( x(2).gt.(2.d0*pi) ) then
                            print *, 'Error: Particle coordinate $\vartheta$ is larger than domain maximum 2.d0*pi.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                        if( x(3).lt.0.d0 ) then
                            print *, 'Error: Particle coordinate $\varphi$ is smaller than domain minimum 0.d0.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                        if( x(3).gt.(2.d0*pi/n_field_periods) ) then
                            print *, 'Error: Particle coordinate $\varphi$ is larger than domain maximum 2.d0*pi/n_field_periods.'
                            print *, 'Hint: Switch on boole_periodic_relocation!'
                            stop
                        endif
!
                    endif
!
            end select
!
        end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module





