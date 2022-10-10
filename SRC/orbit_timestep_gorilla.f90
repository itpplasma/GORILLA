!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module orbit_timestep_gorilla_mod
!
    implicit none
!
    private
!
    public :: orbit_timestep_gorilla,initialize_gorilla,check_coordinate_domain,find_tetra
!
    integer, dimension(:,:), allocatable, public, protected   :: equidistant_grid
    integer, dimension(:), allocatable, public, protected     :: entry_counter
    double precision, dimension(3,5), public, protected       :: dimension_parametres
    logical, public, protected                                :: boole_axi_symmetry
!    
    contains
!
        subroutine orbit_timestep_gorilla(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface,t_remain_out)!,boole_grid_for_find_tetra)
!
            use supporting_functions_mod, only: bmod_func, vperp_func
            use pusher_tetra_rk_mod, only: pusher_tetra_rk,initialize_const_motion_rk
            use pusher_tetra_poly_mod, only: pusher_tetra_poly,initialize_const_motion_poly
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
            !logical, intent(in), optional                   :: boole_grid_for_find_tetra
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
        subroutine initialize_gorilla(boole_grid_for_find_tetra,i_option,ipert_in,bmod_multiplier)
!
            use constants, only: echarge,ame,amp,clight
            use gorilla_settings_mod, only: eps_Phi,coord_system,ispecies
            use tetra_grid_mod, only: make_tetra_grid
            use tetra_physics_mod, only: make_tetra_physics,check_tetra_overlaps,cm_over_e,particle_charge,particle_mass
            use tetra_physics_poly_precomp_mod, only: make_precomp_poly
            use tetra_grid_settings_mod, only: grid_kind
!
            implicit none
!
            integer                                 :: iper,ipert
            logical, intent(inout), optional        :: boole_grid_for_find_tetra
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
                if (present(boole_grid_for_find_tetra)) then
                    if (grid_kind.eq.1) boole_grid_for_find_tetra = .false. !rectangular grid
                    if (boole_grid_for_find_tetra) call grid_for_find_tetra
                endif
!
            endif !boole_physics
!
        end subroutine initialize_gorilla
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine grid_for_find_tetra
!
            use constants, only: pi, eps
            use tetra_grid_mod, only: verts_rphiz, tetra_grid, ntetr
            use tetra_grid_settings_mod, only: grid_size, grid_kind
            use tetra_physics_mod, only: tetra_physics, coord_system
            !use pusher_tetra_poly_mod, only: normal_distance_func
!
            implicit none            
            double precision                                    :: amin, amax, cmin, cmax, delta_a, delta_c, delta_b, &
                                                                & tetr_amin, tetr_amax, tetr_cmin, tetr_cmax, half_diagonal
            integer                                             :: i,k,m,a,b,c, box_index, a_factor, b_factor, c_factor, &
                                                    & num_columns, num_hexahedra,na, nb, nc, maxtetr, maxb, index_change
            integer, dimension(2)                               :: steps_a, steps_b, steps_c
            logical                                             :: boole_reduce_entries
            integer, dimension(:,:), allocatable                :: box_centres
            double precision, dimension(4)                      :: normal_distances
            double precision, dimension(3)                      :: current_box_centre
!
            a_factor = 4
            b_factor = 2
            c_factor = 4
            boole_reduce_entries = .true.
            boole_axi_symmetry = .false.
!
            na = grid_size(1)*a_factor !(R in cylindrical coordinates, s in flux coordinates)
            nb = grid_size(2)*b_factor !(phi in cylindrical and flux coordinates)
            nc = grid_size(3)*c_factor !(z in cylindrical coordinates, theta in flux coordinates)
!
            if (coord_system.eq.1) then
                nb = grid_size(3)*c_factor
                nc = grid_size(2)*b_factor
            endif
            if (grid_kind.eq.4) then
                na = 200*a_factor
                nc = 200*c_factor
                boole_axi_symmetry = .true.
            endif
            !boole_axi_symmetry = .false.
            maxb = nb
            if (boole_axi_symmetry) maxb = b_factor
            num_hexahedra = na*maxb*nc
!PRINT*, 'number of hexahedra is', num_hexahedra
!
            allocate(entry_counter(num_hexahedra))
            allocate(box_centres(num_hexahedra,3))
!
            entry_counter = 0
            box_centres = 0
!
            amin = minval(verts_rphiz(1,:)) - 2*eps
            amax = maxval(verts_rphiz(1,:)) + 2*eps
            cmin = minval(verts_rphiz(3,:)) - 2*eps
            cmax = maxval(verts_rphiz(3,:)) + 2*eps
!
            delta_a = (amax - amin)/na
            delta_c = (cmax - cmin)/nc
            delta_b = 2*pi/nb
!
            dimension_parametres(1,:) = (/delta_a,dble(a_factor),dble(na),amin,amax/)
            dimension_parametres(2,:) = (/delta_b,dble(b_factor),dble(maxb),dble(0),2*pi/)
            dimension_parametres(3,:) = (/delta_c,dble(c_factor),dble(nc),cmin,cmax/)
            if (boole_axi_symmetry) dimension_parametres(2,4) = dimension_parametres(2,4)/grid_size(2)
!
            maxtetr = ntetr
            if (boole_axi_symmetry) maxtetr = ntetr/grid_size(2)
!
            do i = 1, maxtetr
                tetr_amin = minval(verts_rphiz(1,tetra_grid(i)%ind_knot([1,2,3,4]))) - eps
                tetr_amax = maxval(verts_rphiz(1,tetra_grid(i)%ind_knot([1,2,3,4]))) + eps
                tetr_cmin = minval(verts_rphiz(3,tetra_grid(i)%ind_knot([1,2,3,4]))) - eps
                tetr_cmax = maxval(verts_rphiz(3,tetra_grid(i)%ind_knot([1,2,3,4]))) + eps
!
!Print*, tetra_grid(i)%ind_knot([1,2,3,4])
                steps_a(1) = int((tetr_amin-amin)/delta_a) + 1
                steps_a(2) = int((tetr_amax-amin)/delta_a) + 1
                steps_b(1) = nint(minval(verts_rphiz(2,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b) + 1
                steps_b(2) = nint(maxval(verts_rphiz(2,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b)
                steps_c(1) = int((tetr_cmin-cmin)/delta_c) + 1
                steps_c(2) = int((tetr_cmax-cmin)/delta_c) + 1
!
                if (verts_rphiz(2,tetra_grid(i)%ind_knot(1)).gt.verts_rphiz(2,tetra_grid(i)%ind_knot(4))) then !correct the cases where max phi is actually 2*pi but is instead given as 0
                    steps_b(1) = steps_b(2) + 1
                    steps_b(2) = nb
                endif
!
! if ((i.eq.220853).or.(i.eq.320310)) then
!     PRINT*, 'steps_b(1) and steps_b(2) are :', steps_b(1), steps_b(2)
!     PRINT*, verts_rphiz(2,tetra_grid(i)%ind_knot([1,2,3,4]))/delta_b
! endif
!
                do b = steps_b(1), steps_b(2)
                    do c = steps_c(1), steps_c(2)
                        do a = steps_a(1), steps_a(2)
                            box_index = (b-1)*na*nc + & !go to the correct slice
                                    & (c-1)*na + & !go to correct height/tetha value
                                    & a !go to correct box
                                entry_counter(box_index) = entry_counter(box_index)+1
                        enddo
                    enddo
                enddo
            enddo
!
            num_columns = maxval(entry_counter)
            allocate(equidistant_grid(num_hexahedra,num_columns))
            equidistant_grid = -1
            entry_counter = 0
!
            do i = 1, maxtetr
                tetr_amin = minval(verts_rphiz(1,tetra_grid(i)%ind_knot([1,2,3,4])))
                tetr_amax = maxval(verts_rphiz(1,tetra_grid(i)%ind_knot([1,2,3,4])))
                tetr_cmin = minval(verts_rphiz(3,tetra_grid(i)%ind_knot([1,2,3,4])))
                tetr_cmax = maxval(verts_rphiz(3,tetra_grid(i)%ind_knot([1,2,3,4])))
!
!Print*, tetra_grid(i)%ind_knot([1,2,3,4])
                steps_a(1) = int((tetr_amin-amin)/delta_a) + 1
                steps_a(2) = int((tetr_amax-amin)/delta_a) + 1
                steps_b(1) = nint(minval(verts_rphiz(2,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b) + 1
                steps_b(2) = nint(maxval(verts_rphiz(2,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b)
                steps_c(1) = int((tetr_cmin-cmin)/delta_c) + 1
                steps_c(2) = int((tetr_cmax-cmin)/delta_c) + 1
!
                if (verts_rphiz(2,tetra_grid(i)%ind_knot(1)).gt.verts_rphiz(2,tetra_grid(i)%ind_knot(4))) then !correct the cases where max phi is actually 2*pi but is instead given as 0
                    ! if ((i.gt.ntetr/grid_size(2)).and.(steps_b(1).eq.1)) then !correct the cases where max phi is actually 2*pi but is instead given as 0
                    steps_b(1) = steps_b(2) + 1
                    steps_b(2) = nb
                endif
!
                do b = steps_b(1), steps_b(2)
                    do c = steps_c(1), steps_c(2)
                        do a = steps_a(1), steps_a(2)
                            box_index = (b-1)*na*nc + & !go to the correct slice
                                    & (c-1)*na + & !go to correct height/tetha value
                                    & a !go to correct box
                                entry_counter(box_index) = entry_counter(box_index)+1
                                equidistant_grid(box_index,entry_counter(box_index)) = i
                        enddo
                    enddo
                enddo
            enddo
!
! PRINT*, 'maximum number of entries is: ', maxval(entry_counter)
! PRINT*, 'location of maximum number of entries is: ', maxloc(entry_counter)
! PRINT*, 'average number of entries is: ', sum(entry_counter), 'divided by', num_hexahedra
! PRINT*, 'number of hexahedra is: ', num_hexahedra
! PRINT*, 'total number of entries is', sum(entry_counter)
! PRINT*, ntetr, maxtetr, num_columns
! PRINT*, delta_a, delta_b, delta_c
! PRINT*, amin, amax, cmin, cmax
! PRINT*, na, nb, nc
!
            if (boole_reduce_entries) then
                do a = 1,na!fill 1st column (R in cylindrical coordinates, s in flux coordinates)
                    do b = 1,maxb                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                    box_centres((b-1)*na*nc+a : b*na*nc : na ,1) = amin + (a-1)*delta_a + delta_a/2
                    enddo
                enddo
                do b = 1,maxb !fill 2nd column (phi in cylindrical and flux coordinates)
                    box_centres((b-1)*na*nc+1 : b*na*nc,2) = (b-1)*delta_b + delta_b/2
                enddo
                do c = 1,nc!fill 3rd column (z in cylindrical coordinates, theta in flux coordinates)
                    do b = 1,maxb
                    box_centres((b-1)*na*nc+(c-1)*na+1 : (b-1)*na*nc+c*na,3) = cmin + (c-1)*delta_c + delta_c/2
                    enddo
                enddo
!
                half_diagonal = sqrt(delta_a**2+delta_b**2+delta_c**2) + eps
!
                do i = 1,num_hexahedra
                    index_change = 0
                    do k = 1,entry_counter(i)
                        current_box_centre = box_centres(i,:) - tetra_physics(equidistant_grid(i,k+index_change))%x1
                        do m = 1,4
                            normal_distances(m) = sum(current_box_centre* & 
                                                  & tetra_physics(equidistant_grid(i,k+index_change))%anorm(:,m))
                        enddo
                        normal_distances(1) = normal_distances(1) + tetra_physics(equidistant_grid(i,k+index_change))%dist_ref
!PRINT*, tetra_physics(equidistant_grid(i,k))%dist_ref !why is this so often zero???
                        if (any(normal_distances.lt.-half_diagonal)) then
                            equidistant_grid(i,k+index_change:num_columns) = (/equidistant_grid(i,k+index_change+1:num_columns),-1/)
                            entry_counter(i) = entry_counter(i) - 1
                            index_change = index_change -1
                        elseif (all(normal_distances.gt.half_diagonal)) then
                            equidistant_grid(i,1) = equidistant_grid(i,k+index_change)
                            equidistant_grid(i,2:num_columns) = -1
                            entry_counter(i) = 1
                            exit
                        endif
                    enddo
                enddo
!
! PRINT*, 'number of entries is maximally: ', maxval(entry_counter), ' and minimally: ', minval(entry_counter)
! PRINT*, 'location of maximum number of entries is: ', maxloc(entry_counter)
! PRINT*, 'average number of entries is: ', sum(entry_counter)/num_hexahedra
! PRINT*, 'number of hexahedra is: ', num_hexahedra
! PRINT*, 'total number of entries is', sum(entry_counter)
!PRINT*, 'entries of hexahedra with maximum number of entries are :', equidistant_grid(maxloc(entry_counter),:)
            endif
!
            !deallocate(entry_counter, box_centres, equidistant_grid)
PRINT*, 'grid_for_find_tetra is finished'
        end subroutine grid_for_find_tetra
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine find_tetra(x,vpar,vperp,ind_tetr_out,iface,boole_grid_for_find_tetra)
!
        use tetra_grid_mod, only : Rmin,Rmax,Zmin,Zmax,ntetr,tetra_grid, verts_sthetaphi
        use tetra_grid_settings_mod, only: grid_kind, grid_size, n_field_periods
        use tetra_physics_mod, only: tetra_physics,cm_over_e,isinside,coord_system
        use pusher_tetra_func_mod, only: pusher_handover2neighbour
        use constants, only: pi,clight,eps
        use supporting_functions_mod, only: logical2integer
        use pusher_tetra_rk_mod, only: normal_velocity_func, normal_distances_func, perpinv, dist_min,initialize_const_motion_rk, &
                                    & initialize_pusher_tetra_rk_mod, rk4_step
        !use orbit_timestep_gorilla_mod, only: equidistant_grid, entry_counter, dimension_parametres, boole_axi_symmetry
!
        implicit none
!
        double precision, dimension(3), intent(inout) :: x
        double precision, intent(in) :: vpar,vperp
        logical, intent(in), optional :: boole_grid_for_find_tetra
!
        integer, intent(out) :: ind_tetr_out,iface
!
        integer :: ir,iphi,iz,ind_search_tetra, indtetr_start, ind_normdist, ind_tetr_save
        integer :: ind_plane_tetra_start, ind_phi_hexa, ntetr_in_plane, numerical_corr
        integer :: nr, nphi,nphi_tetr, nphi_hexa, nz
        integer :: iper_phi
        integer :: n_plane_conv, l, counter_vnorm_pos, iface_new, iface_new_save, i_tetra_try, i
        integer :: na, nb, nc, b_factor, a, b, c, hexahedron_index, index_phi, ntetr_searched
        double precision :: amin, bmin, cmin, delta_a, delta_b, delta_c
        double precision ::hr,hphi,hz,vnorm,vperp2
        double precision, dimension(4):: cur_dist_value, z
        double precision, dimension(3) :: x_save
        double precision, dimension(4)   :: dzdtau
        logical, dimension(4) :: boole_plane_conv,boole_plane_conv_temp
        logical :: use_grid
        integer, dimension(:), allocatable :: ind_tetr_tried
!
        ! Initialize numerical correction
        numerical_corr = 0
!
        !Calculation of search domain depending on grid_kind
        select case(grid_kind)
            case(1) !rectangular grid
!
                nr = grid_size(1)
                nphi = grid_size(2)
                nz = grid_size(3)
!
                hr=(Rmax-Rmin)/nr !calculating the discretization step sizes in these directions
                hphi=(2.d0 * pi)/nphi
                hz=(Zmax-Zmin)/nz
!
                ir=int((x(1)-Rmin)/hr) +1 !get the reference coordinate in the (indexwise: 1-based) grid, so for instance the particle is in box (ir,iphi,iz)=(4,2,7)
                iphi=int(x(2)/hphi) +1
                iz=int((x(3)-Zmin)/hz) +1
!
                !Diagnostic to check, if out of domain
                if(ir.lt.1.or.ir.gt.nr.or.iphi.lt.1.or.iphi.gt.nphi.or.iz.lt.1.or.iz.gt.nz) then
                    ind_tetr_out=-1
                    iface=-1
                    print *, 'Error in start_tetra_rect: Particle is outside of the domain!'
                    stop
                endif
!
                indtetr_start = int((dble(iz)-1.d0)*6.d0 + 6.d0*dble(nz)*(dble(ir)-1.d0) & 
                & + 6.d0*(dble(iphi)-1.d0)*dble(nr)*dble(nz) +1.d0)
                    !This calculates the starting tetrahedron index as a function of (ir,iphi,iz)
                    !+1, since formula is 0-based (like istarttetr in make_grid_rect)
                ntetr_searched = 6
                use_grid = .false.
!
            case(2,3,4) !EFIT field-aligned grid or VMEC field_aligned grid or SOLEDGE3X_EIRENE
!
                index_phi = 2
                if (coord_system.eq.2) index_phi = 3
                nphi_tetr = grid_size(index_phi)
                ntetr_in_plane = ntetr/nphi_tetr !number of tetrahedra in a slice which is delimited by two phi=const. planes (with phi2-phi1 = 1*delta_phi)
!
                if(present(boole_grid_for_find_tetra)) then
                    if(boole_grid_for_find_tetra) then
                        use_grid = .true.
                    else
                        use_grid = .false.
                    endif
                else
                    use_grid = .false.
                endif
                if (use_grid.eqv..true.) then
                    na = dimension_parametres(1,3)
                    nb = dimension_parametres(2,3)
                    nc = dimension_parametres(3,3)
                    amin = dimension_parametres(1,4)
                    bmin = dimension_parametres(2,4)
                    cmin = dimension_parametres(3,4)
                    b_factor = int(dimension_parametres(2,2))
                    delta_a = dimension_parametres(1,1)
                    delta_b = dimension_parametres(2,1)
                    delta_c = dimension_parametres(3,1)
!
                    nphi_hexa = nphi_tetr*b_factor
                    ind_phi_hexa = int(x(index_phi)*nphi_hexa/(2.d0*pi/n_field_periods))
                    if(abs(x(index_phi)*nphi_hexa/(2.d0*pi/n_field_periods) - &
                    &  dble(ind_phi_hexa)).gt.(1.d0-eps)) numerical_corr = 1
!
                    a = (x(1)-amin)/delta_a + 1
                    b = (x(index_phi)-bmin)/delta_b + 1
                    c = (x(5 - index_phi)-cmin)/delta_c + 1
                    hexahedron_index = (b-1)*na*nc + (c-1)*na + a
                    if (boole_axi_symmetry) hexahedron_index = (b-int((b-1)/b_factor)*b_factor-1)*na*nc + (c-1)*na + a
!
                    if ((b.eq.nb).or.((b.eq.b_factor).and.boole_axi_symmetry)) numerical_corr = 0 !num_correction is not carried out if we are in the last phi slice
! PRINT*, size(entry_counter), hexahedron_index, size(entry_counter) - hexahedron_index, numerical_corr
! PRINT*, int(b/b_factor), (b-int((b-1)/b_factor)*b_factor-1), a, b, c, b_factor
! PRINT*, amin, bmin, cmin
! PRINT*, entry_counter(hexahedron_index)
                    ntetr_searched = entry_counter(hexahedron_index) + & 
                                                    & numerical_corr*entry_counter(hexahedron_index + numerical_corr*na*nc)
                    if (ntetr_searched.eq.0) then
                        print*, 'starting position is out of computation domain'
                        ind_tetr_out = -1
                        iface = -1
                        return
                    endif
                else
                    ind_plane_tetra_start = int(x(index_phi)*nphi_tetr/(2.d0*pi/n_field_periods))
                    if(abs(x(index_phi)*nphi_tetr/(2.d0*pi/n_field_periods) - &
                    &  dble(ind_plane_tetra_start)).gt.(1.d0-eps)) numerical_corr = 1
!
                    indtetr_start = ind_plane_tetra_start*ntetr_in_plane +1
                    ntetr_searched = ntetr_in_plane*(1+numerical_corr)
                endif
        end select
!
!print*, count(entry_counter.eq.0),size(entry_counter)
        do i = 1,ntetr_searched
!
            if (use_grid.eqv..true.) then
                if (i.le.entry_counter(hexahedron_index)) then
                    ind_search_tetra = equidistant_grid(hexahedron_index,i)
                else
                    ind_search_tetra = equidistant_grid(hexahedron_index + na*nc,i - entry_counter(hexahedron_index))
                endif
                if (boole_axi_symmetry) then
                    ind_search_tetra = ind_search_tetra + int(((b-1)/b_factor))*ntetr_in_plane
                endif
            else
                ind_search_tetra = indtetr_start + i -1
            endif
!
            if (isinside(ind_search_tetra,x)) then !inside tetrahedron
!
            ind_tetr_out = ind_search_tetra
            iface = 0
!
            do ind_normdist = 1,4 !calculate distances
                if (ind_normdist .ne. 1) then
                cur_dist_value(ind_normdist) = &
                    & sum(tetra_physics(ind_search_tetra)%anorm(:,ind_normdist)*(x-tetra_physics(ind_search_tetra)%x1))
                else ! ind_normdist .eq. 1
                cur_dist_value(ind_normdist) = sum(tetra_physics(ind_search_tetra)%anorm(:,ind_normdist)*&
                &(x-tetra_physics(ind_search_tetra)%x1))+tetra_physics(ind_search_tetra)%dist_ref
                endif
            enddo
!
            !Check, if particle is in vicinity of a plane (converged on plane)
            boole_plane_conv = abs(cur_dist_value) .le. (eps*abs(tetra_physics(ind_search_tetra)%dist_ref))
            n_plane_conv = sum(logical2integer(boole_plane_conv),1)
!
!print *, 'n_plane_conv',n_plane_conv
!
            if ( n_plane_conv.gt.0 ) then !if it is inside and next to a plane, set iface to the index of where it is 0
!
                !Initialize starting values and 'working' constants in module for tetrahedron
                !Squared perpendicular velocity
                vperp2 = vperp**2
!
                !Parallel velocity and position
                z(1:3) = x-tetra_physics(ind_tetr_out)%x1
                z(4) = vpar
!
                !Temporary iface
                iface_new = minloc(abs(cur_dist_value),1)
!
                ! Allocate and initialize vector with tried indices
                allocate(ind_tetr_tried(2*n_plane_conv))
                ind_tetr_tried = 0
!
                !Compute perpendicular invariant of particle
                perpinv=-0.5d0*vperp2/(tetra_physics(ind_tetr_out)%bmod1+sum(tetra_physics(ind_tetr_out)%gb*z(1:3)))
                call initialize_const_motion_rk(perpinv,perpinv**2)
!
                !Loop over all possible tetrahedra options
                try_loop: do i_tetra_try = 1,(2*n_plane_conv)
!
!print *, 'i_tetra_try',i_tetra_try
!print *, 'ind_tetr',ind_tetr_out
!print *, 'x',x
!
                    !Write 'tried' tetrahedron in 'memory'-vector
                    ind_tetr_tried(i_tetra_try) = ind_tetr_out
!
                    z(1:3) = x-tetra_physics(ind_tetr_out)%x1
!
                    call initialize_pusher_tetra_rk_mod(ind_tetr_out,x,iface_new,vpar,-1.d0)
!
                    cur_dist_value = normal_distances_func(z(1:3))
!print *, 'norm',cur_dist_value
!
                    !Check, if particle is inside the tetrahedron
                    if (any(cur_dist_value.lt.(-dist_min) )) then
                        print *, 'Error in find_tetra: Particle is not inside the tetrahedron while searching.'
                        stop
                    endif
!
                    boole_plane_conv_temp = abs(cur_dist_value) .le. (eps*abs(tetra_physics(ind_tetr_out)%dist_ref))
!
                    !call rk4 at z with dtau = 0 to get the velocity
                    call rk4_step(z,0.d0,dzdtau)
!
                    !Now, ALL normal velocites, where particle is converged  on plane, must be positive
                    counter_vnorm_pos = 0
                    do l = 1,4
                        if(.not.boole_plane_conv_temp(l)) cycle
                        vnorm = normal_velocity_func(l,dzdtau,z)
!print *, 'l',l,'norm(l)', normal_distance_func(z(1:3),l),'vnorm(l)',vnorm, 'boole',boole_plane_conv_temp(l)
                        if (vnorm .gt. 0.d0) then
                            counter_vnorm_pos = counter_vnorm_pos + 1
                        endif
                    enddo
!
!print *, 'counter_vnorm_pos',counter_vnorm_pos

                    if (counter_vnorm_pos.eq.n_plane_conv) then
                        iface = iface_new
!print *, 'gotcha'
                        exit try_loop!this means the tetrahedron was found
                    else
                        ind_tetr_save = ind_tetr_out
                        x_save = x
                        iface_new_save = iface_new
!
                        do l = 1,4
                            if(.not.boole_plane_conv_temp(l)) cycle
                            vnorm = normal_velocity_func(l,dzdtau,z)
                            if (vnorm.gt.0.d0) cycle
                            iface_new = l
                            call pusher_handover2neighbour(ind_tetr_save,ind_tetr_out,iface_new,x,iper_phi)
!
                            if( any(ind_tetr_out.eq.ind_tetr_tried)) then
!print *, 'Denied ind_tetr',ind_tetr_out
!print *, 'Denied iface_new', l
                                x = x_save
                                iface_new = iface_new_save
                            else
!print *, 'Succesful iface_new',l
                                exit !New try
                            endif
                        enddo
                        
                    endif

                enddo try_loop !l = i_tetra_try,n_exactly_conv
!
                deallocate(ind_tetr_tried)
!
            endif ! n_plane_conv.gt.0
!
            exit !this means the tetrahedron was found
!
            else !not inside
                ind_tetr_out = -1
                iface = -1
            endif !isinside
        enddo
!
        if (ind_tetr_out .eq. -1) then
            print *, ' in start_tetra: Starting tetrahedron was not found.'
            return
        endif
!
    end subroutine find_tetra
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





