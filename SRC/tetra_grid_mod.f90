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
                                         & R0_analytic_circ, a_analytic_circ, B0_analytic_circ, q0_analytic_circ, q1_analytic_circ, &
                                         & sfc_s_min, sfc_s_max
        use new_vmec_stuff_mod, only: nper
        use spline_vmec_data_mod, only: spline_vmec_data
        use field_divB0_mod, only: field
        use field_mod, only: ianalytic_circ
        use field_analytic_circ_mod, only: set_field_analytic_circ
        !use make_grid_rect_mod, only: make_grid_rect
!
        implicit none
!
        integer :: efit_vmec,i
        double precision :: rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
        double precision :: rho_ac,s_edge_ac,twopi_ac,eps_ac,th_geom_ac
        ! Variables for the flux-aligned analytic-circ mesh (grid_kind=5)
        double precision :: rho_inner_ac, rho_outer_ac, rho_v, theta_v
        integer :: iv_loc, ind_tetr1_ac, ind_tetr2_ac, iface1_ac, iface2_ac
        integer :: ir_ac, iphi_ac, iz_ac, ii_ac, jj_ac, nr_ac, nphi_ac, nz_ac
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
                case(3) !3D VMEC equilibria
                    call set_n_field_periods(0) !The correct value is assigned below ( befor subroutine 'make_grid_aligned' is called.)
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
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &      !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
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
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &      !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
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
                call set_n_field_periods(nper) !The correct value is assigned below ( befor subroutine 'make_grid_aligned' is called.)
            endif
!
            call make_grid_aligned(grid_size,efit_vmec,n_field_periods)
!            
          case(4) !SOLEDGE3X_EIRENE-grid               
!
            !subroutine field must be called at start in order to read input file and provide necessary quantities in field_eq_mod
            rrr=1.d0
            ppp=0.d0
            zzz=0.d0
!
            call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &      !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
                          & dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
            call make_grid_SOLEDGE3X_EIRENE(grid_size)

          case(5) !analytic circular tokamak — flux-surface-aligned mesh in (rho,theta,phi)
            ianalytic_circ = 1
            call set_field_analytic_circ(R0_analytic_circ, a_analytic_circ, &
                                         B0_analytic_circ, q0_analytic_circ, q1_analytic_circ)
              twopi_ac  = 8.d0*atan(1.d0)
              s_edge_ac = R0_analytic_circ - sqrt(R0_analytic_circ**2 - a_analytic_circ**2)
              !
              ! Build the mesh in flux-surface-aligned (rho, theta, phi) coordinates.
              ! Replacing the previous rectangular (R,Z,phi) mesh: in that layout only
              ! ~1.6% of cells fell within the sfc_s annulus, diluting the delta-f signal
              ! by the same factor.  Here we map:
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
                ! verts_rphiz(2,i) = phi  (unchanged)
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
                    if (tetra_grid(ind_tetr1_ac)%neighbour_tetr(iface1_ac) .eq. -1) then
                      tetra_grid(ind_tetr1_ac)%neighbour_tetr(iface1_ac) = ind_tetr2_ac
                      tetra_grid(ind_tetr1_ac)%neighbour_face(iface1_ac) = iface2_ac
                    end if
                    if (tetra_grid(ind_tetr2_ac)%neighbour_tetr(iface2_ac) .eq. -1) then
                      tetra_grid(ind_tetr2_ac)%neighbour_tetr(iface2_ac) = ind_tetr1_ac
                      tetra_grid(ind_tetr2_ac)%neighbour_face(iface2_ac) = iface1_ac
                    end if
                  end if
                end do
                end do
              end do
              end do
              !
              ! Set verts_sthetaphi: for the flux-aligned mesh, the coordinates are exact.
              !   s     = (R0 - sqrt(R0^2 - rho^2)) / s_edge_ac
              !   theta = STRAIGHT-FIELD-LINE poloidal angle theta*  (NOT geometric)
              !   phi   = toroidal angle (pass-through)
              !
              ! The previous code stored the geometric angle atan2(Z,R-R0) under the
              ! false premise "= SFL angle for circular geometry".  That is only true
              ! as eps=rho/R0 -> 0.  Here eps reaches ~0.37, so geometric and SFL
              ! angles differ by up to ~19 deg at the q=3 surface, which de-coheres the
              ! delta-f helical phase factor exp(i(m*theta + n*phi)) along resonant
              ! field lines (m=-6 amplifies a theta error 6x) and corrupts both the
              ! amplitude and the phase of the screening current.
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
            open(123, file=filename_mesh_rphiz)
            101 format(1000(e21.14,x))
            do i=1, (nvert / grid_size(2)) * 2
                write(123, '(A2)', advance="no") "v "
                write(123,101) verts_rphiz(1, i)*cos(verts_rphiz(2, i)), &
                  &verts_rphiz(1, i)*sin(verts_rphiz(2, i)), &
                  &verts_rphiz(3, i)
            end do
            do i=1, ntetr / grid_size(2)
                write(123, '(A2)', advance="no") "f "
                write(123, *) tetra_grid(i)%ind_knot([2, 3, 4])
                write(123, '(A2)', advance="no") "f "
                write(123, *) tetra_grid(i)%ind_knot([1, 4, 3])
                write(123, '(A2)', advance="no") "f "
                write(123, *) tetra_grid(i)%ind_knot([1, 2, 4])
                write(123, '(A2)', advance="no") "f "
                write(123, *) tetra_grid(i)%ind_knot([1, 3, 2])
            end do
            close(123)
    !
            ![s,theta,phi]: Write Mesh to File
            if ((grid_kind .eq. 2).or.(grid_kind .eq. 3)) then
    !
                open(123, file=filename_mesh_sthetaphi)
                do i=1, (nvert / grid_size(2)) * 2
                    write(123, '(A2)', advance="no") "v "
                    write(123,101) verts_sthetaphi(1, i), &
                                verts_sthetaphi(2, i), &
                                verts_sthetaphi(3, i)
                end do
                do i=1, ntetr / grid_size(2)
                    write(123, '(A2)', advance="no") "f "
                    write(123, *) tetra_grid(i)%ind_knot([2, 3, 4])
                    write(123, '(A2)', advance="no") "f "
                    write(123, *) tetra_grid(i)%ind_knot([1, 4, 3])
                    write(123, '(A2)', advance="no") "f "
                    write(123, *) tetra_grid(i)%ind_knot([1, 2, 4])
                    write(123, '(A2)', advance="no") "f "
                    write(123, *) tetra_grid(i)%ind_knot([1, 3, 2])
                end do
                close(123)
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
  end subroutine make_grid_SOLEDGE3X_EIRENE
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
!
    end subroutine deallocate_tetra_grid
!
end module tetra_grid_mod

