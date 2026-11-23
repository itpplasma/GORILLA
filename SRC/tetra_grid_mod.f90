!
  module tetra_grid_mod
!
  implicit none
!---------------------------------------------------------------------------------------------------------
    type tetrahedron_grid
      sequence
      integer, dimension(4) :: ind_knot              ! pointer from vertex index to knot index
      integer, dimension(4) :: neighbour_tetr        ! pointer from face index to the index of next tetrahedron
      integer, dimension(4) :: neighbour_face        ! index of next tetrahedron entry face from the exit face
      integer, dimension(4) :: neighbour_perbou_phi  ! 1 if the face is on periodic boundary $\phi=2\pi$
                                                     !-1 if on \phi=0 and 0 otherwise
      integer, dimension(4) :: neighbour_perbou_theta! 1 if the face is on periodic boundary $\theta=2\pi$
                                                     !-1 if on \theta=0 and 0 otherwise
    end type tetrahedron_grid
!---------------------------------------------------------------------------------------------------------
    type(tetrahedron_grid), dimension(:),   allocatable, public, protected :: tetra_grid
    double precision,  dimension(:,:), allocatable, public, protected :: verts_rphiz
    double precision,  dimension(:,:), allocatable, public, protected :: verts_xyz
    double precision,  dimension(:,:), allocatable, public, protected :: verts_sthetaphi
    double precision,  dimension(:), allocatable, public, protected :: verts_theta_vmec
!
    !grid_kind = 2 extended to the wall (boole_extend_to_wall): flux coordinates of the flux-aligned core.
    !verts_s_theta_core(:,i) = (s, theta) of vertex i, (-1,-1) for vertices outside the core (edge, axis fill);
    !tetra_in_core(i) = .true. if tetrahedron i belongs to the grid_kind 2 core. Inside such a tetrahedron, s
    !and theta can be interpolated linearly from its vertices (theta: mind the 2 pi seam at theta = 0).
    double precision,  dimension(:,:), allocatable, public, protected :: verts_s_theta_core
    logical,  dimension(:), allocatable, public, protected :: tetra_in_core
    !tetra_in_edge_region(i) = .true. if tetrahedron i lies outside the flux surface s_edge_region
    !(tetra_grid_settings_mod): all non-core tetrahedra (edge; axis fill only for s_edge_region = 0) and the core
    !tetrahedra whose mean vertex s exceeds s_edge_region. The region boundary inside the core is a ring of the
    !mesh, so it consists of tetrahedron faces.
    logical,  dimension(:), allocatable, public, protected :: tetra_in_edge_region
    ! grid_kind = 5 with the magnetic axis in the mesh (sfc_s_min = 0): flat tetrahedra of the collapsed
    ! innermost cells. They keep their index (structured numbering) but have no neighbours and are
    ! skipped by find_tetra; see close_axis_analytic_circ.
    logical,  dimension(:), allocatable, public, protected :: tetra_degenerate
    integer, public, protected :: ntetr,nvert
    double precision, public, protected :: Rmin,Rmax,Zmin,Zmax
!
    contains
!
      subroutine make_tetra_grid()
!
        use field_eq_mod, only : nrad,nzet,rad,zet
        use tetra_grid_settings_mod, only: boole_n_field_periods,n_field_periods_manual,grid_size,n1,n2,n3,grid_kind, &
                                         & n_field_periods,set_grid_size,set_n_field_periods, &
                                         & boole_write_mesh_obj,filename_mesh_rphiz,filename_mesh_sthetaphi, &
                                         & netcdf_filename, &
                                         & R0_analytic_circ, a_analytic_circ, B0_analytic_circ, &
                                         & q0_analytic_circ, q1_analytic_circ, &
                                         & q_exp_analytic_circ, q_profile_file_analytic_circ, sfc_s_min, sfc_s_max
        use tetra_grid_settings_mod, only: helical_pert_file_analytic_circ, &
            helical_pert_amp_analytic_circ, helical_pert_m_analytic_circ, &
            helical_pert_n_analytic_circ
        use new_vmec_stuff_mod, only: nper
        use spline_vmec_data_mod, only: spline_vmec_data
        use boozer_chartmap_mod, only: load_boozer_chartmap, chartmap_nfp
        use field_divB0_mod, only: field
        use field_mod, only: ianalytic_circ
        use field_analytic_circ_mod, only: set_field_analytic_circ, write_flux_functions_analytic_circ
        use field_analytic_circ_mod, only: load_helical_pert_circ
        !use make_grid_rect_mod, only: make_grid_rect
!
        implicit none
!
        integer :: efit_vmec,i,iunit
        double precision :: rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
        double precision :: rho_ac,s_edge_ac,twopi_ac,eps_ac,th_geom_ac
        ! Variables for the flux-aligned analytic-circ mesh (grid_kind=5)
        double precision :: rho_inner_ac, rho_outer_ac, rho_v, theta_v
        integer :: iv_loc, ind_tetr1_ac, ind_tetr2_ac, iface1_ac, iface2_ac
        integer :: ir_ac, iphi_ac, ii_ac, jj_ac, nr_ac, nphi_ac, nz_ac
!
        !Free any pre-existing grid arrays so a fresh grid can be built without external bookkeeping.
        call deallocate_tetra_grid()
!
        call set_grid_size([n1,n2,n3])      !Rectangular grid:                  [n1,n2,n3] = [nR,nphi,nZ]
                                            !Field-aligned grid:                [n1,n2,n3] = [ns,nphi,ntheta]
!        
        !Initialize efit_vmec
        efit_vmec = 0 !0 ... not defined, 1 ... efit, 2 ... vmec
!
        !Select number of field periods - automatically or manually
        if(boole_n_field_periods) then !Automatically
            !Select number of field periods depending on input equilibrium
            select case(grid_kind)
                case(1,2,4,5) !2D EFIT & WEST equilibria
                    call set_n_field_periods(1)
                case(3,6) !3D VMEC or Boozer chartmap equilibria
                    !The correct value is assigned below (before subroutine 'make_grid_aligned' is called.)
                    call set_n_field_periods(0)
            end select
!
        else !Manually from input file
            call set_n_field_periods(n_field_periods_manual)
        endif
!
        !Construction of selected grid
        select case(grid_kind)
          case(1) !rectangular grid
!
            !subroutine field must be called at start in order to read input file and provide necessary quantities in field_eq_mod
            rrr=1.d0
            ppp=0.d0
            zzz=0.d0
!
            !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &
                       & dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!                   
            !Settings for rectangular grid: Values are obtained from field_eq_mod
            Rmin=rad(1)
            Rmax=rad(nrad)
            Zmin=zet(1)
            Zmax=zet(nzet)
!
            !Allocate tetrahedron array
            ntetr = grid_size(1)*grid_size(2)*grid_size(3)*6                !Number of tetrahedra
            nvert = (grid_size(1)+1)*(grid_size(2)+1)*(grid_size(3)+1)      !Number of vertices
            allocate(tetra_grid(ntetr))
            allocate(verts_rphiz(3,nvert))
!
            !Construct rectangular grid in external subroutine
            call make_grid_rect(tetra_grid,verts_rphiz,grid_size,Rmin,Rmax,Zmin,Zmax)
! 
          case(2) !field-aligned grid EFIT
!           
            efit_vmec = 1
!
            !subroutine field must be called at start in order to read input file and provide necessary quantities in field_eq_mod
            rrr=1.d0
            ppp=0.d0
            zzz=0.d0
!
            !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &
                       & dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!         
            call make_grid_aligned(grid_size,efit_vmec,n_field_periods)
!
          case(3) !field-aligned grid VMEC
!
            efit_vmec = 2
!
            !Spline the VMEC equilibrium data
            call spline_vmec_data
            print *, 'Spline VMEC DATA finished'
!
            if(boole_n_field_periods) then
                !Assign value for number of field periods from VMEC equilibrium file
                !The correct value is assigned below (before subroutine 'make_grid_aligned' is called.)
                call set_n_field_periods(nper)
            endif
!
            call make_grid_aligned(grid_size,efit_vmec,n_field_periods)
!
          case(6) !field-aligned grid in direct Boozer chartmap coordinates
!
            efit_vmec = 3
            call load_boozer_chartmap(netcdf_filename)
            if(boole_n_field_periods) call set_n_field_periods(chartmap_nfp)
            call make_grid_aligned(grid_size,efit_vmec,n_field_periods)
!            
          case(4) !SOLEDGE3X_EIRENE-grid               
!
            !subroutine field must be called at start in order to read input file and provide necessary quantities in field_eq_mod
            rrr=1.d0
            ppp=0.d0
            zzz=0.d0
!
            !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &
                          & dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
            call make_grid_SOLEDGE3X_EIRENE(grid_size)

          case(5) !analytic circular tokamak — flux-surface-aligned mesh in (rho,theta,phi)
            ianalytic_circ = 1
            call set_field_analytic_circ(R0_analytic_circ, a_analytic_circ, &
                                         B0_analytic_circ, q0_analytic_circ, q1_analytic_circ, &
                                         q_exp_analytic_circ, q_profile_file_analytic_circ)
            ! grid_kind=5 is analytic: generate its own self-consistent equil-mapping
            ! (psi_pol == tetra_physics%Aphi1) so no external flux_functions.dat is needed.
            call write_flux_functions_analytic_circ()
            if (len_trim(helical_pert_file_analytic_circ) > 0) &
                call load_helical_pert_circ(trim(helical_pert_file_analytic_circ), &
                    helical_pert_amp_analytic_circ, helical_pert_m_analytic_circ, &
                    helical_pert_n_analytic_circ)
              twopi_ac  = 8.d0*atan(1.d0)
              s_edge_ac = R0_analytic_circ - sqrt(R0_analytic_circ**2 - a_analytic_circ**2)
              !
              ! Build the mesh in flux-surface-aligned (rho, theta, phi) coordinates:
              !   n1 cells in rho  from rho(sfc_s_min) to rho(sfc_s_max)
              !   n2 cells in phi  from 0 to 2*pi  (toroidal, unchanged)
              !   n3 cells in theta from 0 to 2*pi  (poloidal, passed as Zmin/Zmax=0/2pi)
              ! make_grid_rect generates the topology; vertices are post-processed to
              ! physical (R=R0+rho*cos(theta), phi, Z=rho*sin(theta)).
              ! The theta direction is made periodic by identifying the iz=n3 vertex
              ! layer with iz=0 (analogous to phi-periodicity in make_grid_rect).
              !
              rho_inner_ac = sqrt(R0_analytic_circ**2 - &
                                  (R0_analytic_circ - sfc_s_min*s_edge_ac)**2)
              rho_outer_ac = sqrt(R0_analytic_circ**2 - &
                                  (R0_analytic_circ - sfc_s_max*s_edge_ac)**2)
              Rmin = rho_inner_ac   ! "R" direction = rho
              Rmax = rho_outer_ac
              Zmin = 0.d0           ! "Z" direction = theta (0 to 2*pi)
              Zmax = twopi_ac
              ntetr = grid_size(1)*grid_size(2)*grid_size(3)*6
              nvert = (grid_size(1)+1)*(grid_size(2)+1)*(grid_size(3)+1)
              allocate(tetra_grid(ntetr))
              allocate(verts_rphiz(3, nvert))
              call make_grid_rect(tetra_grid, verts_rphiz, grid_size, Rmin, Rmax, Zmin, Zmax)
              !
              ! Post-process: verts_rphiz now holds (rho, phi, theta).
              ! Convert to actual cylindrical (R, phi, Z).
              do i = 1, nvert
                rho_v   = verts_rphiz(1, i)   ! rho  (= "R" from make_grid_rect)
                theta_v = verts_rphiz(3, i)   ! theta (= "Z" from make_grid_rect)
                verts_rphiz(1, i) = R0_analytic_circ + rho_v * cos(theta_v)  ! R
                verts_rphiz(3, i) = rho_v * sin(theta_v)                       ! Z
                ! Compress phi from the full torus [0,2*pi] built by make_grid_rect
                ! to a single toroidal field period [0, 2*pi/n_field_periods] for the
                ! straight-cylinder / large-aspect-ratio limit (n_field_periods = N).
                ! make_grid_rect's phi-periodic stitch is topological (index-based,
                ! see "Neighbours through the periodic boundary"), so it stays valid
                ! after this coordinate rescale; and the pusher already wraps phi by
                ! 2*pi/n_field_periods on a periodic crossing (pusher_tetra_func_mod).
                ! n_field_periods = 1 -> rescale by 1.0 -> exact no-op (full torus).
                verts_rphiz(2, i) = verts_rphiz(2, i) / dble(n_field_periods)  ! phi -> one field period
              end do
              !
              ! Theta-periodicity: the iz=n3 vertex layer (theta=2*pi) is physically
              ! identical to iz=0 (theta=0).  Replace iz=n3 vertex indices with iz=0.
              ! In make_grid_rect with nz=grid_size(3), a vertex iv is at iz=n3 when
              ! modulo(iv-1, nz+1) == nz; the corresponding iz=0 vertex is iv - nz.
              nz_ac = grid_size(3)
              do i = 1, ntetr
                do ii_ac = 1, 4
                  iv_loc = tetra_grid(i)%ind_knot(ii_ac)
                  if (modulo(iv_loc-1, nz_ac+1) .eq. nz_ac) then
                    tetra_grid(i)%ind_knot(ii_ac) = iv_loc - nz_ac
                  end if
                end do
              end do
              !
              ! Add neighbour connectivity across the theta=0/2*pi seam.
              ! After the index replacement above, tetrahedra at iz=n3 share vertices
              ! with iz=0; check_neighbour will find the matching faces.
              nr_ac   = grid_size(1)
              nphi_ac = grid_size(2)
              do iphi_ac = 1, nphi_ac
              do ir_ac   = 1, nr_ac
                ! Cell at iz=1  starts at tetra index 6*((iphi_ac-1)*nr_ac*nz_ac + (ir_ac-1)*nz_ac + 0) + 1
                ! Cell at iz=nz starts at tetra index 6*((iphi_ac-1)*nr_ac*nz_ac + (ir_ac-1)*nz_ac + nz_ac-1) + 1
                ! Check ALL 6x6 sub-tetra pairs between the iz=0 and iz=nz cells:
                ! the face shared across the theta=0/2pi seam is generally NOT at the
                ! same sub-tetra index in both cells, so a diagonal-only (ii,ii) check
                ! leaves most seam faces with neighbour=-1 and markers fall out there.
                do ii_ac = 1, 6
                do jj_ac = 1, 6
                  ind_tetr1_ac = 6*((iphi_ac-1)*nr_ac*nz_ac + (ir_ac-1)*nz_ac + 0)       + ii_ac
                  ind_tetr2_ac = 6*((iphi_ac-1)*nr_ac*nz_ac + (ir_ac-1)*nz_ac + nz_ac-1) + jj_ac
                  call check_neighbour(tetra_grid(ind_tetr1_ac)%ind_knot(:), &
                                       tetra_grid(ind_tetr2_ac)%ind_knot(:), &
                                       iface1_ac, iface2_ac)
                  if (iface1_ac .ne. -1) then
                    ! Seam faces are periodic boundaries in theta for coord_system=2
                    ! (pusher_tetra_func_mod): leaving the iz=nz cell (theta -> 2*pi)
                    ! subtracts 2*pi (+1), leaving the iz=1 cell (theta -> 0) adds it (-1).
                    if (tetra_grid(ind_tetr1_ac)%neighbour_tetr(iface1_ac) .eq. -1) then
                      tetra_grid(ind_tetr1_ac)%neighbour_tetr(iface1_ac) = ind_tetr2_ac
                      tetra_grid(ind_tetr1_ac)%neighbour_face(iface1_ac) = iface2_ac
                      tetra_grid(ind_tetr1_ac)%neighbour_perbou_theta(iface1_ac) = -1
                    end if
                    if (tetra_grid(ind_tetr2_ac)%neighbour_tetr(iface2_ac) .eq. -1) then
                      tetra_grid(ind_tetr2_ac)%neighbour_tetr(iface2_ac) = ind_tetr1_ac
                      tetra_grid(ind_tetr2_ac)%neighbour_face(iface2_ac) = iface1_ac
                      tetra_grid(ind_tetr2_ac)%neighbour_perbou_theta(iface2_ac) = 1
                    end if
                  end if
                end do
                end do
              end do
              end do
              !
              ! Magnetic axis inside the mesh: merge the coinciding axis vertices and close the
              ! innermost ring as a triangle fan (see close_axis_analytic_circ).
              if (rho_inner_ac .eq. 0.d0) call close_axis_analytic_circ(nr_ac, nphi_ac, nz_ac)
!
              ! Set verts_sthetaphi: for the flux-aligned mesh, the coordinates are exact.
              !   s     = (R0 - sqrt(R0^2 - rho^2)) / s_edge_ac
              !   theta = STRAIGHT-FIELD-LINE poloidal angle theta*  (NOT geometric)
              !   phi   = toroidal angle (pass-through)
              ! The geometric angle atan2(Z,R-R0) equals theta* only as eps=rho/R0 -> 0.
              !
              ! For the analytic field B_phi = B0*R0/R the field-line winding is
              !   dphi/dtheta_geom = q / (1 + eps*cos(theta_geom)),
              ! whose integral gives the exact SFL relation
              !   theta* = 2*atan( sqrt((1-eps)/(1+eps)) * tan(theta_geom/2) ).
              ! We use the atan2 form below so it is single-valued and robust for all
              ! theta_geom in [0, 2*pi).
              allocate(verts_sthetaphi(3, nvert))
              do i = 1, nvert
                rho_ac = sqrt((verts_rphiz(1,i) - R0_analytic_circ)**2 + verts_rphiz(3,i)**2)
                verts_sthetaphi(1,i) = (R0_analytic_circ - sqrt(R0_analytic_circ**2 - rho_ac**2)) / s_edge_ac
                eps_ac     = rho_ac / R0_analytic_circ
                th_geom_ac = atan2(verts_rphiz(3,i), verts_rphiz(1,i) - R0_analytic_circ)
                verts_sthetaphi(2,i) = modulo( 2.d0*atan2( sqrt(1.d0-eps_ac)*sin(0.5d0*th_geom_ac), &
                                                           sqrt(1.d0+eps_ac)*cos(0.5d0*th_geom_ac) ), twopi_ac)
                verts_sthetaphi(3,i) = verts_rphiz(2,i)
              end do
!
        end select
!
        !Compute vertices in Cartesian coordinates
        allocate(verts_xyz(3,nvert))
        do i = 1,nvert
          verts_xyz(1,i) = verts_rphiz(1,i)*cos(verts_rphiz(2,i))
          verts_xyz(2,i) = verts_rphiz(1,i)*sin(verts_rphiz(2,i))
          verts_xyz(3,i) = verts_rphiz(3,i)
        enddo
!
        !Write mesh to .obj file
        if(boole_write_mesh_obj) then
!
            ![R,phi,Z]: Write Mesh to File
            open(newunit=iunit, file=filename_mesh_rphiz)
            101 format(1000(e21.14,x))
            do i=1, (nvert / grid_size(2)) * 2
                write(iunit, '(A2)', advance="no") "v "
                write(iunit,101) verts_rphiz(1, i)*cos(verts_rphiz(2, i)), &
                  &verts_rphiz(1, i)*sin(verts_rphiz(2, i)), &
                  &verts_rphiz(3, i)
            end do
            do i=1, ntetr / grid_size(2)
                write(iunit, '(A2)', advance="no") "f "
                write(iunit, *) tetra_grid(i)%ind_knot([2, 3, 4])
                write(iunit, '(A2)', advance="no") "f "
                write(iunit, *) tetra_grid(i)%ind_knot([1, 4, 3])
                write(iunit, '(A2)', advance="no") "f "
                write(iunit, *) tetra_grid(i)%ind_knot([1, 2, 4])
                write(iunit, '(A2)', advance="no") "f "
                write(iunit, *) tetra_grid(i)%ind_knot([1, 3, 2])
            end do
            close(iunit)
    !
            ![s,theta,phi]: Write Mesh to File
            if ((grid_kind .eq. 2).or.(grid_kind .eq. 3).or.(grid_kind .eq. 6)) then
    !
                open(newunit=iunit, file=filename_mesh_sthetaphi)
                do i=1, (nvert / grid_size(2)) * 2
                    write(iunit, '(A2)', advance="no") "v "
                    write(iunit,101) verts_sthetaphi(1, i), &
                                verts_sthetaphi(2, i), &
                                verts_sthetaphi(3, i)
                end do
                do i=1, ntetr / grid_size(2)
                    write(iunit, '(A2)', advance="no") "f "
                    write(iunit, *) tetra_grid(i)%ind_knot([2, 3, 4])
                    write(iunit, '(A2)', advance="no") "f "
                    write(iunit, *) tetra_grid(i)%ind_knot([1, 4, 3])
                    write(iunit, '(A2)', advance="no") "f "
                    write(iunit, *) tetra_grid(i)%ind_knot([1, 2, 4])
                    write(iunit, '(A2)', advance="no") "f "
                    write(iunit, *) tetra_grid(i)%ind_knot([1, 3, 2])
                end do
                close(iunit)
            endif
!
        endif !boole_write_mesh_obj
!        
      end subroutine make_tetra_grid
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine make_grid_aligned(grid_size,efit_vmec,n_field_periods)
!
    use constants, only: pi
    use circular_mesh, only : calc_mesh, create_points, calc_points_circular, calc_n_tetras, calc_n_verts
    use scaling_r_theta, only: scaling_r, scaling_theta
    use preload_for_SYNCH_mod, only: preload_for_SYNCH
    use magdata_in_symfluxcoordinates_mod, only: load_magdata_in_symfluxcoord
!
    implicit none
!
    integer:: i,j,k,l,i_vertex, ind_tetr
    integer, intent(in) :: efit_vmec,n_field_periods
    integer, dimension(3), intent(in) :: grid_size
    integer, dimension(grid_size(1)) :: verts_per_ring
    integer, allocatable, dimension(:, :) :: verts, neighbours, neighbour_faces, perbou_phi, perbou_theta
!
    double precision, dimension(4):: cur_r,cur_phi,cur_z
    double precision, dimension(8)   :: davec_ds,davec_dp,davec_dz
    double precision, dimension(4,8) :: avec
    double precision, dimension(:),allocatable ::As,Ap,Az
    double precision, dimension(:),allocatable :: h_s,h_p,h_z,bmod,phi_elec
    double precision, dimension(:),         allocatable :: rnd_noise
    double precision :: scalprod
    double precision, dimension(3) :: curCoordi,curCoordk,curCoordj
    integer :: mesh_nr,nphi,mesh_ntheta
!
    mesh_nr = grid_size(1)
    nphi = grid_size(2)
    mesh_ntheta = grid_size(3)
!
    verts_per_ring = mesh_ntheta

    !Distinguish inbetween axisymmetric (EFIT) and non-axisymmetric (VMEC) data 
    select case(efit_vmec)
        case(1) !EFIT (axisymmetric)
            call preload_for_SYNCH
            call load_magdata_in_symfluxcoord
            call create_points(verts_per_ring, nphi, verts_rphiz, verts_sthetaphi, efit_vmec,n_field_periods,nvert, &
                       r_scaling_func = scaling_r, theta_scaling_func = scaling_theta, &
                       repeat_center_point = .true.)
        case(2) !VMEC (non-axisymmetric)
            call create_points(verts_per_ring, nphi, verts_rphiz, verts_sthetaphi, efit_vmec,n_field_periods,nvert, &
                       verts_theta_vmec,r_scaling_func = scaling_r, theta_scaling_func = scaling_theta, &
                       repeat_center_point = .true.)
        case(3) !direct Boozer chartmap
            call create_points(verts_per_ring, nphi, verts_rphiz, verts_sthetaphi, &
                               efit_vmec, n_field_periods, nvert, &
                               r_scaling_func=scaling_r, &
                               theta_scaling_func=scaling_theta, &
                               repeat_center_point=.true.)
    end select
!    
    call calc_mesh(verts_per_ring, nphi, verts_rphiz(:, :nvert / nphi), ntetr, &
                   verts, neighbours, neighbour_faces, perbou_phi, perbou_theta, &
                   repeat_center_point = .true.)
!
    allocate(tetra_grid(1:ntetr))
!
    do i=1,ntetr
      tetra_grid(i)%ind_knot = verts(:, i)
      tetra_grid(i)%neighbour_tetr = neighbours(:, i)
      tetra_grid(i)%neighbour_face = neighbour_faces(:, i)
      tetra_grid(i)%neighbour_perbou_phi(:) = perbou_phi(:, i)
      tetra_grid(i)%neighbour_perbou_theta(:) = perbou_theta(:, i)
    enddo
  !
  end subroutine make_grid_aligned
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
 subroutine make_grid_SOLEDGE3X_EIRENE(grid_size)
  !
    use constants, only: pi
    use circular_mesh_SOLEDGE3X_EIRENE, only : calc_mesh_SOLEDGE3X_EIRENE, create_points_SOLEDGE3X_EIRENE
    use tetra_grid_settings_mod, only: boole_extend_to_wall
!
    implicit none
!
    integer:: i, verts_per_slice
    integer, dimension(3), intent(in) :: grid_size
    integer, allocatable, dimension(:, :) :: verts, neighbours, neighbour_faces, perbou_phi
!
    integer :: mesh_nr,nphi,mesh_ntheta
!
    nphi = grid_size(2)
!

    call create_points_SOLEDGE3X_EIRENE(nphi, verts_rphiz,verts_per_slice)
!   
    nvert = verts_per_slice*nphi
    call calc_mesh_SOLEDGE3X_EIRENE(nphi, verts_rphiz, verts_per_slice, ntetr, verts, neighbours, neighbour_faces, perbou_phi)
!
    allocate(tetra_grid(1:ntetr))
!
    do i=1,ntetr
      tetra_grid(i)%ind_knot = verts(:, i)
      tetra_grid(i)%neighbour_tetr = neighbours(:, i)
      tetra_grid(i)%neighbour_face = neighbour_faces(:, i)
      tetra_grid(i)%neighbour_perbou_phi(:) = perbou_phi(:, i)
    enddo
!
    if (boole_extend_to_wall) call load_edge_mesh_core(verts_per_slice, nphi)
  !
  end subroutine make_grid_SOLEDGE3X_EIRENE
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_edge_mesh_core(verts_per_slice, nphi)
!
    !Reads edge_mesh_dir/edge_mesh_core.dat written by the edge generator:
    !  line 1: number of 2D vertices, number of triangles
    !  then one line "s theta" per 2D vertex (-1 -1 outside the grid_kind 2 core)
    !  then one region number per triangle (0 = core, 1 = edge, 2 = axis fill)
    !The external-mesh reader extrudes each triangle into 3 tetrahedra per slice, triangle by triangle, and
    !copies the vertices slice by slice, which gives the mapping below.
!
    use tetra_grid_settings_mod, only: edge_mesh_dir, s_edge_region
!
    implicit none
!
    integer, intent(in) :: verts_per_slice, nphi
    integer :: iu, nv2, ntri, i, k, itri, ios
    double precision, dimension(:,:), allocatable :: st
    integer, dimension(:), allocatable :: region
!
    open(newunit=iu, file=trim(edge_mesh_dir)//'/edge_mesh_core.dat', status='old', action='read', iostat=ios)
    if (ios .ne. 0) then
      print *, 'Error: cannot open ', trim(edge_mesh_dir)//'/edge_mesh_core.dat'
      error stop
    endif
    read(iu,*) nv2, ntri
    if (nv2 .ne. verts_per_slice .or. 3*ntri*nphi .ne. ntetr) then
      print *, 'Error: edge_mesh_core.dat does not match the mesh files (vertices, triangles):', nv2, ntri
      error stop
    endif
    allocate(st(2, nv2), region(ntri))
    do i = 1, nv2
      read(iu,*) st(:, i)
    enddo
    do i = 1, ntri
      read(iu,*) region(i)
    enddo
    close(iu)
!
    if (allocated(verts_s_theta_core)) deallocate(verts_s_theta_core)
    if (allocated(tetra_in_core)) deallocate(tetra_in_core)
    if (allocated(tetra_in_edge_region)) deallocate(tetra_in_edge_region)
    allocate(verts_s_theta_core(2, nvert), tetra_in_core(ntetr), tetra_in_edge_region(ntetr))
    do k = 0, nphi - 1
      verts_s_theta_core(:, k*nv2 + 1 : (k+1)*nv2) = st
    enddo
    do i = 1, ntetr
      itri = mod(i - 1, 3*ntri) / 3 + 1
      tetra_in_core(i) = (region(itri) .eq. 0)
      if (s_edge_region .le. 0.d0) then
        tetra_in_edge_region(i) = .true.
      elseif (tetra_in_core(i)) then
        tetra_in_edge_region(i) = (sum(verts_s_theta_core(1, tetra_grid(i)%ind_knot)) .gt. 4.d0*s_edge_region)
      else
        tetra_in_edge_region(i) = (region(itri) .eq. 1)
      endif
    enddo
    print *, 'Extended mesh: ', count(tetra_in_core), ' of ', ntetr, ' tetrahedra in the flux-aligned core'
    print *, 'Edge region (s_edge_region =', s_edge_region, '): ', count(tetra_in_edge_region), ' of ', ntetr, &
             ' tetrahedra'
!
  end subroutine load_edge_mesh_core
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine set_neighbour_face(ind_tetr,ind_face,set_value)
!
  !This subroutine is a setter for neighbour face. In case of overlapping tetrahedra, some faces need to be "deactivated".
  !The subroutine is needed in order to keep the property as public protected.
!
    implicit none
!
    integer :: ind_tetr, ind_face, set_value
!
    tetra_grid(ind_tetr)%neighbour_face(ind_face) = set_value
!
  end subroutine set_neighbour_face
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine set_verts_sthetaphi(i_vertex,i_coord,set_value)
!
  !This subroutine is a setter for verts_s_theta_phi(i_coord,i_vertes).
  !The subroutine is needed in order to keep the property as public protected.
!
    implicit none
!
    integer :: i_vertex, i_coord
    double precision :: set_value
!    
    verts_sthetaphi(i_coord,i_vertex) = set_value
!    
  end subroutine set_verts_sthetaphi  
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine make_grid_rect(tetra,verts_rphiz,grid_size,Rmin,Rmax,Zmin,Zmax)
!
    use constants, only: pi
!   
    implicit none
!
    type(tetrahedron_grid), dimension(ntetr), intent(out) :: tetra
    double precision,  dimension(3,nvert), intent(out) :: verts_rphiz
    integer, dimension(3), intent(in) :: grid_size
    double precision, intent(in) :: Rmin,Rmax,Zmin,Zmax
    integer :: ir,iphi,iz,iv,itetr,i,j,k,l
    integer :: ind_tetr,ind_tetr1,ind_tetr2,iface1,iface2
    integer :: nvertinner,nper
    integer, dimension(4) :: knots1,knots2  
    integer, dimension(3,4,3) :: ivert_prism
    double precision :: hr,hphi,hz,r,phi,z,x,y,s
    double precision :: R_c,Z_c
    double precision, dimension(4)   :: sp,pp,zp
    integer,          dimension(:,:,:),       allocatable :: inodes,itetrbeg
    integer :: nr,nphi,nz
!
    nr = grid_size(1)
    nphi = grid_size(2)
    nz = grid_size(3)
!
! vertices of the first tetrahedron:
! point 1:
    ivert_prism(1,1,1)=0
    ivert_prism(2,1,1)=0
    ivert_prism(3,1,1)=0
! point 2:
    ivert_prism(1,2,1)=1
    ivert_prism(2,2,1)=0
    ivert_prism(3,2,1)=0
! point 3:
    ivert_prism(1,3,1)=0
    ivert_prism(2,3,1)=1
    ivert_prism(3,3,1)=0
! point 4:
    ivert_prism(1,4,1)=0
    ivert_prism(2,4,1)=0
    ivert_prism(3,4,1)=1
!
! vertices of the second tetrahedron:
! point 1:
    ivert_prism(:,1,2)=ivert_prism(:,2,1)
! point 2:
    ivert_prism(:,2,2)=ivert_prism(:,3,1)
! point 3:
    ivert_prism(:,3,2)=ivert_prism(:,4,1)
! point 4:
    ivert_prism(1,4,2)=1
    ivert_prism(2,4,2)=0
    ivert_prism(3,4,2)=1
!
! vertices of the third tetrahedron:
! point 1:
    ivert_prism(:,1,3)=ivert_prism(:,3,1)
! point 2:
    ivert_prism(:,2,3)=ivert_prism(:,4,1)
! point 3:
    ivert_prism(:,3,3)=ivert_prism(:,4,2)
! point 4:
    ivert_prism(1,4,3)=0
    ivert_prism(2,4,3)=1
    ivert_prism(3,4,3)=1
!
    hr=(Rmax-Rmin)/nr
    hphi=2.d0*pi/nphi
    hz=(Zmax-Zmin)/nz
!
    R_c=0.5d0*(Rmax+Rmin)
    Z_c=0.5d0*(Zmax+Zmin)
    nper=0
!
!
    allocate(inodes(0:nr,0:nz,0:nphi))
!
    iv=0
!
    do iphi=0,nphi
    phi=hphi*iphi
    do ir=0,nr
        do iz=0,nz
        r=rmin+hr*ir
        z=zmin+hz*iz
!
        x=(r-R_c)
        y=(z-Z_c)
        r=R_c+x*cos(nper*phi)+y*sin(nper*phi)
        z=Z_c-x*sin(nper*phi)+y*cos(nper*phi)! +hz/3000
        s=r
!
        iv=iv+1
        inodes(ir,iz,iphi)=iv
!
        verts_rphiz(:,iv) = [s,phi,z]
!
        enddo
    enddo
    enddo
!
    allocate(itetrbeg(nr,nz,nphi))
!
    ind_tetr=0
!
    do iphi=1,nphi
    do ir=1,nr
        do iz=1,nz
!
        itetrbeg(ir,iz,iphi)=ind_tetr
!
        do itetr=1,3
            ind_tetr=ind_tetr+1
            do i=1,4
            iv=inodes(ir-1+ivert_prism(1,i,itetr),iz-1+ivert_prism(2,i,itetr),iphi-1+ivert_prism(3,i,itetr))
            tetra(ind_tetr)%ind_knot(i)=iv
            tetra(ind_tetr)%neighbour_tetr(i)=-1
            tetra(ind_tetr)%neighbour_face(i)=-1
            tetra(ind_tetr)%neighbour_perbou_phi(i)=0
            tetra(ind_tetr)%neighbour_perbou_theta(i)=0
            enddo
        enddo
!
        do itetr=1,3
            ind_tetr=ind_tetr+1
            do i=1,4
            iv=inodes(ir-ivert_prism(1,i,itetr),iz-ivert_prism(2,i,itetr),iphi-ivert_prism(3,i,itetr))
            tetra(ind_tetr)%ind_knot(i)=iv
            tetra(ind_tetr)%neighbour_tetr(i)=-1
            tetra(ind_tetr)%neighbour_face(i)=-1
            tetra(ind_tetr)%neighbour_perbou_phi(i)=0
            tetra(ind_tetr)%neighbour_perbou_theta(i)=0
            enddo
        enddo
!
        do i=1,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,6
            if(j.eq.i) cycle
            ind_tetr2=itetrbeg(ir,iz,iphi)+j
            call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
            if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            endif
            enddo
        enddo
        enddo
    enddo
    enddo
!
!
    do iphi=1,nphi
    do ir=1,nr
        do iz=1,nz
!
        if(ir.gt.1) then
            do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
                ind_tetr2=itetrbeg(ir-1,iz,iphi)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        if(ir.lt.nr) then
            do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
                ind_tetr2=itetrbeg(ir+1,iz,iphi)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        if(iz.gt.1) then
            do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
                ind_tetr2=itetrbeg(ir,iz-1,iphi)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        if(iz.lt.nz) then
            do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
                ind_tetr2=itetrbeg(ir,iz+1,iphi)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        if(iphi.gt.1) then
            do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
                ind_tetr2=itetrbeg(ir,iz,iphi-1)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
            do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
                ind_tetr2=itetrbeg(ir,iz,iphi-1)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        if(iphi.lt.nphi) then
            do i=1,3
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=1,3
                ind_tetr2=itetrbeg(ir,iz,iphi+1)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
            do i=4,6
            ind_tetr1=itetrbeg(ir,iz,iphi)+i
            do j=4,6
                ind_tetr2=itetrbeg(ir,iz,iphi+1)+j
                call check_neighbour(tetra(ind_tetr1)%ind_knot(:),              &
                                    tetra(ind_tetr2)%ind_knot(:),iface1,iface2)
                if(iface1.ne.-1) then
                tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
                tetra(ind_tetr1)%neighbour_face(iface1)=iface2
                endif
            enddo
            enddo
        endif
!
        enddo
    enddo
    enddo
!
! Neighbours through the periodic boundary
!
    nvertinner=(nr+1)*(nz+1)*nphi
!
    do ir=1,nr
    do iz=1,nz
!
        do i=1,3
        ind_tetr1=itetrbeg(ir,iz,1)+i
        knots1(:)=tetra(ind_tetr1)%ind_knot(:)
        do j=1,3
            ind_tetr2=itetrbeg(ir,iz,nphi)+j
            knots2(:)=tetra(ind_tetr2)%ind_knot(:)
            knots2(:)=modulo(knots2(:),nvertinner)
            call check_neighbour(knots1,knots2,iface1,iface2)
            if(iface1.ne.-1) then
            tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
            tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            tetra(ind_tetr1)%neighbour_perbou_phi(iface1)=-1
            tetra(ind_tetr2)%neighbour_tetr(iface2)=ind_tetr1
            tetra(ind_tetr2)%neighbour_face(iface2)=iface1
            tetra(ind_tetr2)%neighbour_perbou_phi(iface2)=1
            endif
        enddo
        enddo
!
        do i=4,6
        ind_tetr1=itetrbeg(ir,iz,1)+i
        knots1(:)=tetra(ind_tetr1)%ind_knot(:)
        do j=4,6
            ind_tetr2=itetrbeg(ir,iz,nphi)+j
            knots2(:)=tetra(ind_tetr2)%ind_knot(:)
            knots2(:)=modulo(knots2(:),nvertinner)
            call check_neighbour(knots1,knots2,iface1,iface2)
            if(iface1.ne.-1) then
            tetra(ind_tetr1)%neighbour_tetr(iface1)=ind_tetr2
            tetra(ind_tetr1)%neighbour_face(iface1)=iface2
            tetra(ind_tetr1)%neighbour_perbou_phi(iface1)=-1
            tetra(ind_tetr2)%neighbour_tetr(iface2)=ind_tetr1
            tetra(ind_tetr2)%neighbour_face(iface2)=iface1
            tetra(ind_tetr2)%neighbour_perbou_phi(iface2)=1
            endif
        enddo
        enddo
!
    enddo
    enddo
!
! Check the number of tetrahedrons at the XY boundary
    j=0
!
b: do ind_tetr=1,ntetr
    do i=1,4
        if(tetra(ind_tetr)%neighbour_tetr(i).lt.1) then
        j=j+1
        cycle b
        endif
    enddo
    enddo b
!
    deallocate(inodes,itetrbeg)
!
    end subroutine make_grid_rect
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine check_neighbour(knots1,knots2,iface1,iface2)
!
    implicit none
!
    logical :: result
    integer :: iface1,iface2,i,j,matches
    integer, dimension(4) :: knots1,knots2
    integer, dimension(3) :: match1,match2
    logical, dimension(4) :: unused
!
    iface1=-1
    iface2=-1
    matches=0
!
    a: do i=1,4
    do j=1,4
        if(knots1(i).eq.knots2(j)) then
        matches=matches+1
        match1(matches)=i
        match2(matches)=j
        if(matches.eq.3) exit a
        endif
    enddo
    if(i-matches.gt.1) return
    enddo a
!
    unused= .true.
    do i=1,3
    unused(match1(i))= .false.
    enddo
    do i=1,4
    if(unused(i)) then
        iface1=i
        exit
    endif
    enddo
!
    unused= .true.
    do i=1,3
    unused(match2(i))= .false.
    enddo
    do i=1,4
    if(unused(i)) then
        iface2=i
        exit
    endif
    enddo
!
    end subroutine check_neighbour
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine close_axis_analytic_circ(nr, nphi, nz)
!
! grid_kind = 5 with sfc_s_min = 0: the rho = 0 face of every innermost (rho, theta, phi) cell
! collapses to a line, so the cell is a wedge (triangle A, a, b in the plane phi times the slice).
! make_grid_rect still splits it into 6 tetrahedra with one axis vertex per theta, all at the same
! point; 3 of them are flat. Here the innermost ring is rebuilt as a triangle fan:
!   - the axis vertices of each toroidal plane are merged into one (A),
!   - every wedge is split into 3 tetrahedra; the diagonal of the side faces theta = const is
!     always A(phi)-a(phi'), so neighbouring wedges match, and the outer face takes the diagonal
!     of the adjacent ring-2 cell,
!   - the 3 new tetrahedra use 3 of the 6 index slots of the cell, the other 3 are flagged in
!     tetra_degenerate (no neighbours, skipped by find_tetra); the structured numbering is kept,
!   - the neighbours of ring 1 and the ring-2 faces towards it are matched again by their vertices.
! Cylindrical coordinates only (the axis is singular in symmetry flux coordinates).
!
    use gorilla_settings_mod, only: coord_system
    use tetra_grid_settings_mod, only: R0_analytic_circ, n_field_periods
    use constants, only: pi
!
    implicit none
!
    integer, intent(in) :: nr, nphi, nz
    integer :: i, j, k, k2, f, f1, f2, iphi, iz, ir, n_planes, n_link, n_open, n_v
    integer :: dphi_i, dz_i, jr, jz, jphi, base, base2, ra0, rb0, ra1, rb1, ax0, ax1
    integer :: cell_v(8), kn(4), kn2(4), new_tets(4,3)
    integer, dimension(:), allocatable :: rep, plane_rep, pid
    double precision, dimension(:), allocatable :: plane_phi
    logical, dimension(:), allocatable :: on_axis
    double precision :: tol, period, phi_bot, rz(2)
    logical :: found, diag_a_bp
!
    if (coord_system .ne. 1) then
        print *, 'Error: grid_kind = 5 with the magnetic axis in the mesh (sfc_s_min = 0) needs coord_system = 1;'
        print *, '       in symmetry flux coordinates the axis is singular: set sfc_s_min > 0.'
        error stop
    endif
    if (nphi .lt. 2) error stop 'grid_kind = 5 with the magnetic axis in the mesh needs n2 >= 2'
!
    tol = 1.d-12*R0_analytic_circ
    period = 2.d0*pi/dble(n_field_periods)
!
    ! one representative vertex per toroidal plane for all vertices on the axis
    allocate(rep(nvert), plane_rep(nphi+1), plane_phi(nphi+1), on_axis(nvert), pid(nvert))
    on_axis = .false.
    n_planes = 0
    do i = 1, nvert
        rep(i) = i
        if (hypot(verts_rphiz(1,i) - R0_analytic_circ, verts_rphiz(3,i)) .gt. tol) cycle
        on_axis(i) = .true.
        found = .false.
        do j = 1, n_planes
            if (abs(verts_rphiz(2,i) - plane_phi(j)) .lt. 1.d-6*period/nphi) then
                rep(i) = plane_rep(j)
                found = .true.
                exit
            endif
        enddo
        if (.not. found) then
            n_planes = n_planes + 1
            if (n_planes .gt. nphi + 1) error stop 'close_axis_analytic_circ: more axis planes than toroidal planes'
            plane_rep(n_planes) = i
            plane_phi(n_planes) = verts_rphiz(2,i)
        endif
    enddo
    do k = 1, ntetr
        tetra_grid(k)%ind_knot = rep(tetra_grid(k)%ind_knot)
    enddo
!
    ! periodic identity of the vertices: a vertex on the plane phi = period is the vertex at phi = 0
    ! with the same (R, Z) (make_grid_rect keeps a separate vertex layer there)
    do i = 1, nvert
        pid(i) = rep(i)
    enddo
    do i = 1, nvert
        if (rep(i) .ne. i) cycle
        if (abs(verts_rphiz(2,i) - period) .gt. 1.d-6*period/nphi) cycle
        do j = 1, nvert
            if (rep(j) .ne. j) cycle
            if (abs(verts_rphiz(2,j)) .gt. 1.d-6*period/nphi) cycle
            if (abs(verts_rphiz(1,j) - verts_rphiz(1,i)) + abs(verts_rphiz(3,j) - verts_rphiz(3,i)) .gt. tol) cycle
            pid(i) = j
            exit
        enddo
    enddo
!
    allocate(tetra_degenerate(ntetr))
    tetra_degenerate = .false.
!
    ! rebuild every wedge of ring 1
    do iphi = 1, nphi
        phi_bot = dble(iphi - 1)*period/dble(nphi)
        do iz = 1, nz
            base = 6*((iphi - 1)*nr*nz + iz - 1)
            ! the 6 distinct vertices of the wedge
            n_v = 0
            do k = base + 1, base + 6
                do i = 1, 4
                    if (any(cell_v(1:n_v) .eq. tetra_grid(k)%ind_knot(i))) cycle
                    n_v = n_v + 1
                    if (n_v .gt. 6) error stop 'close_axis_analytic_circ: wedge with more than 6 vertices'
                    cell_v(n_v) = tetra_grid(k)%ind_knot(i)
                enddo
            enddo
            if (n_v .ne. 6) error stop 'close_axis_analytic_circ: wedge with less than 6 vertices'
            ax0 = -1; ax1 = -1; ra0 = -1; rb0 = -1; ra1 = -1; rb1 = -1
            do i = 1, 6
                j = cell_v(i)
                if (abs(verts_rphiz(2,j) - phi_bot) .lt. 1.d-6*period/nphi) then   ! bottom plane
                    if (on_axis(j)) then
                        ax0 = j
                    elseif (ra0 .lt. 0) then
                        ra0 = j
                    else
                        rb0 = j
                    endif
                endif
            enddo
            do i = 1, 6
                j = cell_v(i)
                if (j .eq. ax0 .or. j .eq. ra0 .or. j .eq. rb0) cycle
                if (on_axis(j)) then
                    ax1 = j
                else
                    rz = verts_rphiz([1,3], j)
                    if (abs(rz(1) - verts_rphiz(1,ra0)) + abs(rz(2) - verts_rphiz(3,ra0)) .lt. tol) then
                        ra1 = j
                    else
                        rb1 = j
                    endif
                endif
            enddo
            if (min(ax0, ax1, ra0, rb0, ra1, rb1) .lt. 1) error stop 'close_axis_analytic_circ: cannot classify wedge vertices'
            ! diagonal of the outer quad face (a, b, a', b'): that of the adjacent ring-2 cell
            diag_a_bp = .true.
            if (nr .ge. 2) then
                base2 = 6*((iphi - 1)*nr*nz + nz + iz - 1)
                do k = base2 + 1, base2 + 6
                    kn = tetra_grid(k)%ind_knot
                    if (count(kn .eq. ra0 .or. kn .eq. rb0 .or. kn .eq. ra1 .or. kn .eq. rb1) .ne. 3) cycle
                    diag_a_bp = any(kn .eq. ra0) .and. any(kn .eq. rb1)
                    exit
                enddo
            endif
            ! split: side diagonals A-a', A-b'; outer diagonal a-b' or b-a'
            if (diag_a_bp) then
                new_tets(:,1) = [ax0, ra0, rb0, rb1]
                new_tets(:,2) = [ax0, ra0, ra1, rb1]
            else
                new_tets(:,1) = [ax0, ra0, rb0, ra1]
                new_tets(:,2) = [ax0, rb0, ra1, rb1]
            endif
            new_tets(:,3) = [ax0, ax1, ra1, rb1]
            do i = 1, 6
                k = base + i
                do f = 1, 4      ! drop links from ring 2 into this cell
                    k2 = tetra_grid(k)%neighbour_tetr(f)
                    if (k2 .gt. 0) then
                        f2 = tetra_grid(k)%neighbour_face(f)
                        if (tetra_grid(k2)%neighbour_tetr(f2) .eq. k) then
                            tetra_grid(k2)%neighbour_tetr(f2) = -1
                            tetra_grid(k2)%neighbour_face(f2) = -1
                            tetra_grid(k2)%neighbour_perbou_phi(f2) = 0
                            tetra_grid(k2)%neighbour_perbou_theta(f2) = 0
                        endif
                    endif
                enddo
                tetra_grid(k)%neighbour_tetr = -1
                tetra_grid(k)%neighbour_face = -1
                tetra_grid(k)%neighbour_perbou_phi = 0
                tetra_grid(k)%neighbour_perbou_theta = 0
                if (i .le. 3) then
                    tetra_grid(k)%ind_knot = new_tets(:, i)
                else
                    tetra_grid(k)%ind_knot = [ax0, ax0, ra0, rb0]   ! flat placeholder (zero volume), flagged
                    tetra_degenerate(k) = .true.
                endif
            enddo
        enddo
    enddo
!
    ! match the open faces of ring 1 and ring 2 (same or adjacent cells, periodic in phi and theta)
    n_link = 0
    do iphi = 1, nphi
    do ir = 1, min(2, nr)
    do iz = 1, nz
    do i = 1, 6
        k = 6*((iphi - 1)*nr*nz + (ir - 1)*nz + iz - 1) + i
        if (tetra_degenerate(k)) cycle
        if (all(tetra_grid(k)%neighbour_tetr .ne. -1)) cycle
        kn = pid(tetra_grid(k)%ind_knot)
        do dphi_i = -1, 1
        do jr = 1, min(2, nr)
        do dz_i = -1, 1
            jphi = modulo(iphi - 1 + dphi_i, nphi) + 1
            jz = modulo(iz - 1 + dz_i, nz) + 1
            do j = 1, 6
                k2 = 6*((jphi - 1)*nr*nz + (jr - 1)*nz + jz - 1) + j
                if (k2 .eq. k .or. tetra_degenerate(k2)) cycle
                kn2 = pid(tetra_grid(k2)%ind_knot)
                call check_neighbour(kn, kn2, f1, f2)
                if (f1 .eq. -1) cycle
                if (tetra_grid(k)%neighbour_tetr(f1) .ne. -1 .or. tetra_grid(k2)%neighbour_tetr(f2) .ne. -1) cycle
                tetra_grid(k)%neighbour_tetr(f1) = k2
                tetra_grid(k)%neighbour_face(f1) = f2
                tetra_grid(k2)%neighbour_tetr(f2) = k
                tetra_grid(k2)%neighbour_face(f2) = f1
                tetra_grid(k)%neighbour_perbou_phi(f1) = face_on_plane(k, f1, period) - face_on_plane(k, f1, 0.d0)
                tetra_grid(k2)%neighbour_perbou_phi(f2) = face_on_plane(k2, f2, period) - face_on_plane(k2, f2, 0.d0)
                tetra_grid(k)%neighbour_perbou_theta(f1) = 0
                tetra_grid(k2)%neighbour_perbou_theta(f2) = 0
                n_link = n_link + 1
            enddo
        enddo
        enddo
        enddo
    enddo
    enddo
    enddo
    enddo
!
    ! a valid tetrahedron of rings 1 and 2 may only have open faces on the outer boundary rho = a
    n_open = 0
    do iphi = 1, nphi
    do ir = 1, min(2, nr)
    do iz = 1, nz
    do i = 1, 6
        k = 6*((iphi - 1)*nr*nz + (ir - 1)*nz + iz - 1) + i
        if (tetra_degenerate(k)) cycle
        do f = 1, 4
            if (tetra_grid(k)%neighbour_tetr(f) .ne. -1) cycle
            if (ir .lt. nr) then
                n_open = n_open + 1
            elseif (any(on_axis(tetra_grid(k)%ind_knot(pack([1,2,3,4], [1,2,3,4] .ne. f))))) then
                n_open = n_open + 1
            endif
        enddo
    enddo
    enddo
    enddo
    enddo
    print '(a,i0,a,i0,a,i0,a)', ' grid_kind 5: magnetic axis closed by a triangle fan (', n_planes, &
        ' axis vertices, ', count(tetra_degenerate), ' tetrahedron slots unused, ', n_link, ' faces matched)'
    if (n_open .gt. 0) then
        print *, 'Error: close_axis_analytic_circ: open faces next to the axis: ', n_open
        error stop
    endif
!
    contains
!
        integer function face_on_plane(kk, ff, phi_plane)
            ! 1 if face ff of tetrahedron kk lies in the plane phi = phi_plane, else 0
            integer, intent(in) :: kk, ff
            double precision, intent(in) :: phi_plane
            integer :: m
            face_on_plane = 1
            do m = 1, 4
                if (m .eq. ff) cycle
                if (abs(verts_rphiz(2, tetra_grid(kk)%ind_knot(m)) - phi_plane) .gt. 1.d-6*period/nphi) face_on_plane = 0
            enddo
        end function face_on_plane
!
    end subroutine close_axis_analytic_circ
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine deallocate_tetra_grid
!
! Deallocates tetra_grid arrays to allow rebuilding with different parameters
!
    implicit none
!
    if (allocated(tetra_grid)) deallocate(tetra_grid)
    if (allocated(verts_rphiz)) deallocate(verts_rphiz)
    if (allocated(verts_xyz)) deallocate(verts_xyz)
    if (allocated(verts_sthetaphi)) deallocate(verts_sthetaphi)
    if (allocated(verts_theta_vmec)) deallocate(verts_theta_vmec)
    if (allocated(verts_s_theta_core)) deallocate(verts_s_theta_core)
    if (allocated(tetra_in_core)) deallocate(tetra_in_core)
    if (allocated(tetra_in_edge_region)) deallocate(tetra_in_edge_region)
    if (allocated(tetra_degenerate)) deallocate(tetra_degenerate)
!
    end subroutine deallocate_tetra_grid
!
end module tetra_grid_mod
