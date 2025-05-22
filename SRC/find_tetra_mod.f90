!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module find_tetra_mod

    use, intrinsic :: iso_fortran_env, only: dp => real64
!
    implicit none
!
    private
!
    public :: grid_for_find_tetra, find_tetra
!
    integer, dimension(:,:), allocatable           :: equidistant_grid
    integer, dimension(:), allocatable             :: entry_counter
    logical                                        :: boole_axi_symmetry
    integer                                        :: a,b,c,na, nb, nc, max_nb
    integer                                        :: ind_a, ind_b, ind_c
    real(dp)                                       :: amin, amax, cmin, cmax, bmin, bmax, delta_a, delta_b, delta_c
    real(dp), dimension(:,:), allocatable          :: box_centres, save_box_centres
    integer                                        :: num_hexahedra
!
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine grid_for_find_tetra !compare with chapter 4.1 of master thesis of Jonatan Schatzlmayr
!
        use constants, only: pi, eps
        use tetra_grid_mod, only: verts_rphiz, verts_sthetaphi, tetra_grid, ntetr
        use tetra_grid_settings_mod, only: grid_size, grid_kind, n_field_periods
        use tetra_physics_mod, only: tetra_physics, coord_system
        use gorilla_settings_mod, only: a_factor, b_factor, c_factor
        !use pusher_tetra_poly_mod, only: normal_distance_func
!
        implicit none            
        real(dp)                               :: tetr_amin, tetr_amax, tetr_cmin, tetr_cmax, half_diagonal
        integer                                :: i,k,m, box_index, num_columns, maxtetr, index_change
        integer, dimension(2)                  :: steps_a, steps_b, steps_c
        logical                                :: boole_reduce_entries
        real(dp), dimension(4)                 :: normal_distances, normalisations
        real(dp), dimension(3)                 :: current_box_centre
        real(dp), dimension(:,:), allocatable  :: verts_abc
!
print*, 'grid_for_find_tetra started'
        boole_reduce_entries = .true.
        boole_axi_symmetry = .false.
!
        if (coord_system.eq.1) allocate(verts_abc(size(verts_rphiz(:,1)),size(verts_rphiz(1,:))))
        if (coord_system.eq.2) allocate(verts_abc(size(verts_sthetaphi(:,1)),size(verts_sthetaphi(1,:))))
!
        if (coord_system.eq.1) verts_abc = verts_rphiz
        if (coord_system.eq.2) verts_abc = verts_sthetaphi
!
        ind_a = 1 !(R in cylindrical coordinates, s in flux coordinates)
        ind_b = 2 !(phi in cylindrical and flux coordinates)
        ind_c = 3 !(z in cylindrical coordinates, theta in flux coordinates)
!
        if (coord_system.eq.2) then
            ind_b = 3
            ind_c = 2
        endif
!
        amin = minval(verts_abc(ind_a,:)) - 2*eps
        amax = maxval(verts_abc(ind_a,:)) + 2*eps
        bmin = 0
        bmax = 2*pi
        cmin = minval(verts_abc(ind_c,:)) - 2*eps
        cmax = maxval(verts_abc(ind_c,:)) + 2*eps
        if (boole_axi_symmetry) bmax = 2*pi/grid_size(2)
        if (coord_system.eq.2) then
            cmin = -2*eps
            cmax = 2*pi + 2*eps
        endif
!
        na = grid_size(ind_a)*a_factor
        delta_a = (amax - amin)/na
!
        if (c_factor.eq.0) c_factor = maxval((/nint((cmax-cmin)/(grid_size(ind_c)*delta_a)),1/))
        nc = grid_size(ind_c)*c_factor
!
        if (grid_kind.eq.4) then
            na = 200*a_factor
            delta_a = (amax - amin)/na
            c_factor = maxval((/nint((cmax-cmin)/(200*delta_a)),1/))
            nc = 200*c_factor
            boole_axi_symmetry = .true.
        endif
!
        delta_c = (cmax - cmin)/nc
!
        if (b_factor.eq.0) b_factor = maxval((/nint(2*pi/(grid_size(2)*n_field_periods*sqrt((delta_a**2+delta_c**2)/2))),1/))
!
        nb = grid_size(2)*b_factor
        delta_b = 2*pi/nb
print*, 'abc_factor, delta_abc = ', a_factor, b_factor, c_factor, delta_a, delta_b, delta_c
!
        max_nb = nb
        if (boole_axi_symmetry) max_nb = b_factor
        num_hexahedra = na*max_nb*nc
!print*, 'num_hexahedra = ', num_hexahedra, na, max_nb, nc, na*nc, na*max_nb, max_nb*nc, na*max_nb*nc
!
        allocate(entry_counter(num_hexahedra))
        allocate(box_centres(num_hexahedra,3))
        allocate(save_box_centres(num_hexahedra,3))
!
        entry_counter = 0
        box_centres = 0
        save_box_centres = 0
!
        maxtetr = ntetr
        if (boole_axi_symmetry) maxtetr = ntetr/grid_size(2)
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(maxtetr,verts_abc,ind_a,ind_b,ind_c,amin,cmin,delta_a,delta_b,delta_c,tetra_grid,na,nb,nc,entry_counter, &
        !$OMP& grid_size,n_field_periods,coord_system) &
        !$OMP& PRIVATE(i,tetr_amin,tetr_amax,tetr_cmin,tetr_cmax,steps_a,steps_b,steps_c,box_index,a,b,c)
        !$OMP DO
        do i = 1, maxtetr
            tetr_amin = minval(verts_abc(ind_a,tetra_grid(i)%ind_knot([1,2,3,4]))) - eps
            tetr_amax = maxval(verts_abc(ind_a,tetra_grid(i)%ind_knot([1,2,3,4]))) + eps
            tetr_cmin = minval(verts_abc(ind_c,tetra_grid(i)%ind_knot([1,2,3,4]))) - eps
            tetr_cmax = maxval(verts_abc(ind_c,tetra_grid(i)%ind_knot([1,2,3,4]))) + eps
!
            if ((coord_system.eq.2).and.((tetr_cmax-tetr_cmin).gt.(2*pi/(grid_size(ind_c))+eps))) then !correct the case when theta is 0 but should be 2*pi
                tetr_cmin = tetr_cmax
                tetr_cmax = 2*pi
            endif
!
            steps_a(1) = int((tetr_amin-amin)/delta_a) + 1
            steps_a(2) = int((tetr_amax-amin)/delta_a) + 1
            steps_b(1) = nint(minval(verts_abc(ind_b,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b) + 1
            steps_b(2) = nint(maxval(verts_abc(ind_b,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b)
            steps_c(1) = int((tetr_cmin-cmin)/delta_c) + 1
            steps_c(2) = int((tetr_cmax-cmin)/delta_c) + 1
!
            if (verts_abc(ind_b,tetra_grid(i)%ind_knot(1)).gt.verts_abc(ind_b,tetra_grid(i)%ind_knot(4))) then !correct the cases where max phi is actually 2*pi but is instead given as 0
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
                            !$omp critical
                            entry_counter(box_index) = entry_counter(box_index)+1
                           !$omp end critical
                    enddo
                enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
!
        num_columns = maxval(entry_counter)
        allocate(equidistant_grid(num_hexahedra,num_columns))
        equidistant_grid = -1
        entry_counter = 0
!
        !$OMP PARALLEL DEFAULT(NONE) &
        !$OMP& SHARED(maxtetr,verts_abc,ind_a,ind_b,ind_c,amin,cmin,delta_a,delta_b,delta_c,tetra_grid,na,nb,nc,entry_counter, &
        !$OMP& equidistant_grid,grid_size,n_field_periods,coord_system,num_hexahedra) &
        !$OMP& PRIVATE(i,tetr_amin,tetr_amax,tetr_cmin,tetr_cmax,steps_a,steps_b,steps_c,box_index,a,b,c)
        !$OMP DO
        do i = 1, maxtetr
            tetr_amin = minval(verts_abc(ind_a,tetra_grid(i)%ind_knot([1,2,3,4])))
            tetr_amax = maxval(verts_abc(ind_a,tetra_grid(i)%ind_knot([1,2,3,4])))
            tetr_cmin = minval(verts_abc(ind_c,tetra_grid(i)%ind_knot([1,2,3,4])))
            tetr_cmax = maxval(verts_abc(ind_c,tetra_grid(i)%ind_knot([1,2,3,4])))
!
            if ((coord_system.eq.2).and.((tetr_cmax-tetr_cmin).gt.(2*pi/(grid_size(ind_c))+eps))) then !correct the case when theta is 0 but should be 2*pi
                tetr_cmin = tetr_cmax
                tetr_cmax = 2*pi
            endif
!
            steps_a(1) = int((tetr_amin-amin)/delta_a) + 1
            steps_a(2) = int((tetr_amax-amin)/delta_a) + 1
            steps_b(1) = nint(minval(verts_abc(ind_b,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b) + 1
            steps_b(2) = nint(maxval(verts_abc(ind_b,tetra_grid(i)%ind_knot([1,2,3,4])))/delta_b)
            steps_c(1) = int((tetr_cmin-cmin)/delta_c) + 1
            steps_c(2) = int((tetr_cmax-cmin)/delta_c) + 1
!
            if (verts_abc(ind_b,tetra_grid(i)%ind_knot(1)).gt.verts_abc(ind_b,tetra_grid(i)%ind_knot(4))) then !correct the cases where max phi is actually 2*pi but is instead given as 0
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
                           !$omp critical
                            entry_counter(box_index) = entry_counter(box_index)+1
                            equidistant_grid(box_index,entry_counter(box_index)) = i
                           !$omp end critical
                    enddo
                enddo
            enddo
        enddo
        !$OMP END DO
        !$OMP END PARALLEL
!
 PRINT*, 'average number of entries before reduction is: ', sum(entry_counter), 'divided by', num_hexahedra, ' = ', &
 dble(sum(entry_counter))/dble(num_hexahedra)
!
        if (boole_reduce_entries) then
            do a = 1,na!fill 1st column (R in cylindrical coordinates, s in flux coordinates)
                do b = 1,max_nb                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
                box_centres((b-1)*na*nc+a : b*na*nc : na ,ind_a) = amin + (a-1)*delta_a + delta_a/2
                enddo
            enddo
            do b = 1,max_nb !fill phi column (second in cylindrical and third in flux coordinates)
                box_centres((b-1)*na*nc+1 : b*na*nc,ind_b) = (b-1)*delta_b + delta_b/2
            enddo
            do c = 1,nc !fill z/theta column (third in cylindrical coordinates, second in flux coordinates)
                do b = 1,max_nb
                box_centres((b-1)*na*nc+(c-1)*na+1 : (b-1)*na*nc+c*na,ind_c) = cmin + (c-1)*delta_c + delta_c/2
                enddo
            enddo
!
            half_diagonal = sqrt(delta_a**2+delta_b**2+delta_c**2)/2 + eps
!
            !$OMP PARALLEL DEFAULT(NONE) &
            !$OMP& SHARED(num_hexahedra,entry_counter,box_centres,tetra_physics,equidistant_grid,half_diagonal,num_columns, &
            !$OMP& ind_a,ind_b,ind_c,verts_abc,tetra_grid) &
            !$OMP& PRIVATE(i,index_change,k,current_box_centre,m,normalisations,normal_distances)
            !$OMP DO
            do i = 1,num_hexahedra
                index_change = 0
                do k = 1,entry_counter(i)
                    current_box_centre = box_centres(i,:) - tetra_physics(equidistant_grid(i,k+index_change))%x1
                    do m = 1,4
                        normalisations(m) = sqrt(tetra_physics(equidistant_grid(i,k+index_change))%anorm(1,m)**2+ &
                                                 tetra_physics(equidistant_grid(i,k+index_change))%anorm(2,m)**2+ &
                                                 tetra_physics(equidistant_grid(i,k+index_change))%anorm(3,m)**2)
                    enddo
                    do m = 1,4
                        normal_distances(m) = sum(current_box_centre* & 
                                                & tetra_physics(equidistant_grid(i,k+index_change))%anorm(:,m)/normalisations(m))
                    enddo
                    normal_distances(1) = normal_distances(1) + &
                                          tetra_physics(equidistant_grid(i,k+index_change))%dist_ref/normalisations(1)
                    if (any(normal_distances.lt.-half_diagonal)) then
                        equidistant_grid(i,k+index_change:num_columns) = (/equidistant_grid(i,k+index_change+1:num_columns),-1/)
                        entry_counter(i) = entry_counter(i) - 1
                        index_change = index_change -1
! if (i.eq.14129) then
!     print*, 'hello', normal_distances, half_diagonal
!     print*, entry_counter(14129)
!     print*, equidistant_grid(14129,:)
! endif
                    elseif (all(normal_distances.gt.half_diagonal)) then
                        equidistant_grid(i,1) = equidistant_grid(i,k+index_change)
                        equidistant_grid(i,2:num_columns) = -1
                        entry_counter(i) = 1
                        exit
                    endif
                enddo
            enddo
            !$OMP END DO
            !$OMP END PARALLEL
!
PRINT*, 'average number of entries after reduction is: ', sum(entry_counter), 'divided by', num_hexahedra, ' = ', &
dble(sum(entry_counter))/dble(num_hexahedra)
!PRINT*, 'delta_a, delta_b, delta_c and b_factor are:  ', delta_a,delta_b,delta_c,b_factor
        endif
!
PRINT*, 'grid_for_find_tetra finished'
    end subroutine grid_for_find_tetra
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine find_tetra(x,vpar,vperp,ind_tetr_out,iface,sign_t_step_in)
!
        use gorilla_settings_mod, only: boole_grid_for_find_tetra, b_factor
        use tetra_grid_mod, only : Rmin,Rmax,Zmin,Zmax,ntetr,tetra_grid, verts_rphiz
        use tetra_grid_settings_mod, only: grid_kind, grid_size, n_field_periods
        use tetra_physics_mod, only: tetra_physics,cm_over_e,isinside, coord_system
        use pusher_tetra_func_mod, only: pusher_handover2neighbour
        use constants, only: pi,clight,eps
        use supporting_functions_mod, only: logical2integer
        use pusher_tetra_rk_mod, only: normal_velocity_func, normal_distances_func, dist_min,initialize_const_motion_rk, &
                                    & initialize_pusher_tetra_rk_mod, rk4_step
!
        implicit none
!
        real(dp), dimension(3), intent(inout) :: x
        real(dp), intent(in) :: vpar,vperp
!
        integer, intent(out) :: ind_tetr_out,iface
!
        integer, intent(in), optional :: sign_t_step_in
        !
        integer :: sign_t_step
        integer :: ir,iphi,iz,ind_search_tetra, indtetr_start, ind_normdist, ind_tetr_save
        integer :: ind_plane_tetra_start, ind_phi_hexa, ntetr_in_plane, numerical_corr_plus, numerical_corr_minus
        integer :: nr, nphi,nphi_tetr, nphi_hexa, nz
        integer :: iper_phi
        integer :: n_plane_conv, l, counter_vnorm_pos, iface_new, iface_new_save, i_tetra_try, i
        integer :: hexahedron_index, shifted_hexahedron_index, ntetr_searched
        real(dp) :: perpinv
        real(dp) ::hr,hphi,hz,vnorm,vperp2
        real(dp), dimension(4):: cur_dist_value, z
        real(dp), dimension(3) :: x_save
        real(dp), dimension(4)   :: dzdtau
        logical, dimension(4) :: boole_plane_conv,boole_plane_conv_temp
        logical :: use_grid
        integer, dimension(:), allocatable :: ind_tetr_tried
        integer :: index_tetr
!
        ! Initialize sign of the right hand side of ODE - ensures that tau is ALWAYS positive inside the algorithm
        if(present(sign_t_step_in)) then
            sign_t_step = sign_t_step_in
        else
            sign_t_step = +1
        endif

        ! Initialize numerical correction
        numerical_corr_plus = 0
        numerical_corr_minus = 0
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

                ! Once incorporating corrections when particle sits close to any hexahedron wall, consider corrections for all 6 walls
                ! if () numerical_corr_minus = 1
                ! if () numerical_corr_plus = 1
                ntetr_searched = 6
                use_grid = .false.
!
            case(2,3,4) !EFIT field-aligned grid or VMEC field_aligned grid or SOLEDGE3X_EIRENE
!
                ind_b = 2 !(phi in cylindrical and flux coordinates)
                if (coord_system.eq.2) ind_b = 3
!
                nphi_tetr = grid_size(2)
                nphi = nphi_tetr 
                ntetr_in_plane = ntetr/nphi_tetr !number of tetrahedra in a slice which is delimited by two phi=const. planes (with phi2-phi1 = 1*delta_phi)
!
                !if(present(boole_grid_for_find_tetra)) then
                    if(boole_grid_for_find_tetra) then
                        use_grid = .true.
                    else
                        use_grid = .false.
                    endif

                if (use_grid.eqv..true.) then
!
                    nphi_hexa = nphi_tetr*b_factor
                    ind_phi_hexa = int(x(ind_b)*nphi_hexa/(2.d0*pi/n_field_periods))
!
                    a = (x(ind_a)-amin)/delta_a + 1
                    b = (x(ind_b)-bmin)/delta_b + 1
                    c = (x(ind_c)-cmin)/delta_c + 1
                    if ((a.lt.0).or.(b.lt.0).or.(c.lt.0).or.(a.gt.na).or.(b.gt.nb).or.(c.gt.nc)) then
                        print*, 'starting position is out of computation domain'
                        ind_tetr_out = -1
                        iface = -1
                        return
                    endif
                    hexahedron_index = (b-1)*na*nc + (c-1)*na + a
                    if (boole_axi_symmetry) hexahedron_index = (b-int((b-1)/b_factor)*b_factor-1)*na*nc + (c-1)*na + a

!
                    if(abs(x(ind_b)*nphi_hexa/(2.0_dp*pi/n_field_periods) - dble(ind_phi_hexa)).gt.(1.0_dp-eps)) then
                        numerical_corr_plus = 1
                    endif
                    if((abs(x(ind_b)*nphi_hexa/(2.0_dp*pi/n_field_periods) - dble(ind_phi_hexa)).lt.eps)) then 
                        numerical_corr_minus = 1
                    endif
!
                    ntetr_searched = entry_counter(hexahedron_index) + & 
                        & numerical_corr_plus *entry_counter(mod(hexahedron_index+numerical_corr_plus *na*nc-1,num_hexahedra)+1) + &
                        & numerical_corr_minus* &
                        & entry_counter(mod(hexahedron_index-numerical_corr_minus*na*nc-1+num_hexahedra,num_hexahedra)+1)
                    if (ntetr_searched.eq.0) then
!print*, 'starting position is out of computation domain'
                        ind_tetr_out = -1
                        iface = -1
                        return
                    endif
                else
                    ind_plane_tetra_start = int(x(ind_b)*nphi_tetr/(2.d0*pi/n_field_periods))
                    if(abs(x(ind_b)*nphi_tetr/(2.d0*pi/n_field_periods) - dble(ind_plane_tetra_start)).gt.(1.0_dp-eps)) &
                        &  numerical_corr_plus = 1
                    if((abs(x(ind_b)*nphi_tetr/(2.d0*pi/n_field_periods) - dble(ind_plane_tetra_start)).lt.eps)) &
                        &  numerical_corr_minus = 1
!
                    indtetr_start = ind_plane_tetra_start*ntetr_in_plane +1
                    ntetr_searched = ntetr_in_plane*(1+numerical_corr_plus+numerical_corr_minus)
                endif
        end select
!
        !search tetrahedra one by one 
        do i = 1,ntetr_searched
!
!print*, 'tetra_physics(142)%dist_ref = ', tetra_physics(142)%dist_ref
!stop
            if (use_grid.eqv..true.) then
                if (i.le.entry_counter(hexahedron_index)) then
                    ind_search_tetra = equidistant_grid(hexahedron_index,i)
                else
                    shifted_hexahedron_index = hexahedron_index + numerical_corr_plus*na*nc - numerical_corr_minus*na*nc
                    if (shifted_hexahedron_index.gt.num_hexahedra) shifted_hexahedron_index = shifted_hexahedron_index-num_hexahedra
                    if (shifted_hexahedron_index.le.0) shifted_hexahedron_index = shifted_hexahedron_index+num_hexahedra
                    ind_search_tetra = equidistant_grid(shifted_hexahedron_index, i - entry_counter(hexahedron_index))
                endif
                if (boole_axi_symmetry) then
                    ind_search_tetra = ind_search_tetra + int(((b-1)/b_factor))*ntetr_in_plane
                endif
            else
                ind_search_tetra = indtetr_start + i - 1 - numerical_corr_minus*ntetr/nphi
            endif

            ! Once incorporating corrections when particle sits close to any hexahedron wall, consider corrections for all 6 walls
            ! if (grid_kind.eq.1) then
            !     if (i.le.6) then
            !         ind_search_tetra = indtetr_start + i - 1
            !     else 
            !         ind_search_tetra = indtetr_start + i - 7 - numerical_corr_minus*ntetr/nphi + numerical_corr_plus*ntetr/nphi
            !     endif
            ! endif
            if (ind_search_tetra.gt.ntetr) ind_search_tetra = ind_search_tetra - ntetr
            if (ind_search_tetra.le.0) ind_search_tetra = ind_search_tetra + ntetr

            if (isinside(ind_search_tetra,x)) then !inside tetrahedron
!
                ind_tetr_out = ind_search_tetra
                iface = 0

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
                        if (ind_tetr_out.eq.-1) exit
!
                        !Write 'tried' tetrahedron in 'memory'-vector
                        ind_tetr_tried(i_tetra_try) = ind_tetr_out
!
                        z(1:3) = x-tetra_physics(ind_tetr_out)%x1
!
                        call initialize_pusher_tetra_rk_mod(ind_tetr_out,x,iface_new,vpar,dble(sign_t_step))
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
                                if(any(ind_tetr_out.eq.ind_tetr_tried).or.(ind_tetr_out.eq.-1)) then
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

                if (ind_tetr_out.eq.-1) then
                    iface = -1
                else
                    exit !this means the tetrahedron was found
                endif
!
            else !not inside
                ind_tetr_out = -1
                iface = -1
            endif !isinside
        enddo
!
        if (ind_tetr_out .eq. -1) then
            if (.not.boole_grid_for_find_tetra) print *, ' in start_tetra: Starting tetrahedron was not found.'
            ! open(35, file = 'outliers.dat',position = 'append')
            !     write(35,'(2ES20.10E4)') x(1), x(3)
            ! close(35)
            ! print*, 'x = ', x
            return
        endif
!
    end subroutine find_tetra
end module find_tetra_mod