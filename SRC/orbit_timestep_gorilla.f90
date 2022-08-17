!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module orbit_timestep_gorilla_mod
!
    implicit none
!
    private
!
    public :: orbit_timestep_gorilla,initialize_gorilla,phi_elec_func,check_coordinate_domain, bmod_func, vperp_func
!    
    contains
!
        subroutine orbit_timestep_gorilla(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface,t_remain_out)
!
            use pusher_tetra_rk_mod, only: find_tetra,pusher_tetra_rk,initialize_const_motion_rk
            use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly,manage_intermediate_steps_arrays
            use tetra_physics_poly_precomp_mod , only: make_precomp_poly_perpinv, initialize_boole_precomp_poly_perpinv, &
                & alloc_precomp_poly_perpinv
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            !use tetra_grid_mod, only: ntetr
            use gorilla_settings_mod, only: ipusher, poly_order
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
!
            !If orbit_timestep is called for the first time without grid position
            if(.not.boole_initialized) then
!
                !Check coordinate domain (optionally perform modulo operation)
                call check_coordinate_domain(x)
!
                !Find tetrahedron index and face index for position x
                call find_tetra(x,vpar,vperp,ind_tetr,iface)
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
            !Deallocate intermediate_steps_arrays if pusher_poly was used
            call manage_intermediate_steps_arrays()
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
        subroutine initialize_gorilla(i_option,ipert_in,bmod_multiplier)
!
            use constants, only: echarge,ame,amp,clight
            use gorilla_settings_mod, only: eps_Phi,coord_system,ispecies
            use tetra_grid_mod, only: make_tetra_grid
            use tetra_physics_mod, only: make_tetra_physics,check_tetra_overlaps,cm_over_e,particle_charge,particle_mass
            use tetra_physics_poly_precomp_mod, only: make_precomp_poly
!
            implicit none
!
            integer                                 :: iper,ipert
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
            double precision :: vpar,perpinv,p_phi_func,hphi1
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
end module





