!
  module tetra_physics_mod
!
    implicit none
!---------------------------------------------------------------------------------------------------------
    type tetrahedron_physics
      sequence
      !Geometric properties
      double precision, dimension(3)   :: x1     ! coordinates of the first vertex
      double precision :: dist_ref               ! normal distance of vertex1 to face1 * |n_1| <-- normal vector to face 1
      double precision, dimension(4)    :: dist_ref_vec
      double precision :: tetra_dist_ref         ! distance reference of a tetrahdron in unit length
      double precision, dimension(3,4) :: anorm  ! normals to faces pointing invards
      !Field properties
      double precision, dimension(3)   :: curlA  ! $\epsilon^{ijk}\difp{A_k}{x^j}$ = vector a_4^i ("curl A")
      double precision :: bmod1                  ! module of B in the first vertex
      double precision :: Aphi1                  ! phi-component of the vector potential A in the first vertex
      double precision :: h2_1                   ! 2nd component of the unit vector h in the first vertex
      double precision :: h3_1                   ! 3rd component of the unit vector h in the first vertex
      double precision :: Phi1                   ! electrostatic potential in the first vertex
      double precision :: R1                     ! major radius in the first vertex
      double precision :: Z1                     ! Z in the first vertex
      double precision :: Er_mod                 ! module of the electric field in r-direction
      double precision :: sqg1                   ! square root g (metric determinant) at the first vertex (ONLY for grid_kind = 3)
      double precision :: dt_dtau_const          ! Factor dt_dtau averaged of the four vertices
      double precision :: gBxcurlA               ! $\epsilon^{ijk}\difp{B}{x^i}\difp{A_k}{x^j}$ scalar
                                                 ! of "curl A" with gradient of B module
      double precision :: gPhixcurlA             ! $\epsilon^{ijk}\difp{Phi}{x^i}\difp{A_k}{x^j}$ scalar
                                                 ! of "curl A" with gradient of Phi
      double precision :: spalpmat               ! alpha matrix element(4,4)= trace of the real space
                                                 ! part of matrix alpha (see below)
      double precision :: spbetmat               ! beta matrix element(4,4)= trace of the real space
                                                 ! part of matrix beta (see below)
      double precision, dimension(3)   :: gBxh1  ! $\epsilon^{ijk} \difp{B}{x^j} h_k$ - vector product
                                                 ! of gradient B module with $\bh$ in the first vertex
      double precision, dimension(3)   :: gPhixh1  ! $\epsilon^{ijk} \difp{Phi}{x^j} h_k$ - vector product
                                                 ! of gradient Phi with $\bh$ in the first vertex
      double precision, dimension(3)   :: gB     ! gradient of B module
      double precision, dimension(3)   :: gPhi   ! gradient of Phi
      double precision, dimension(3)   :: gR     ! gradient of R (major radius)
      double precision, dimension(3)   :: gZ     ! gradient of Z
      double precision, dimension(3)   :: gsqg   ! gradient of square root g (metric determinant) (ONLY for grid_kind = 3)
      double precision, dimension(3)   :: gAphi  ! gradient of phi-component of the vector potential A
      double precision, dimension(3)   :: gh1    ! gradient of the 1st component of the unit vector h
      double precision, dimension(3)   :: gh2    ! gradient of the 2nd component of the unit vector h
      double precision, dimension(3)   :: gh3    ! gradient of the 3rd component of the unit vector h
      double precision, dimension(3)   :: curlh  ! $\epsilon^{ijk}\difp{h_k}{x^j}$ = "curl" of $\bh$
      !3DGeoInt properties
      double precision, dimension(3,3) :: alpmat ! real space part of matrix alpha (block from 1 to 3)
      double precision, dimension(3,3) :: betmat ! real space part of matrix beta (block from 1 to 3)
      double precision, dimension(4)     :: acoef_pre              !Factor of acoef for analytical quadratic approximation 
! 
    end type tetrahedron_physics
!
    type(tetrahedron_physics), dimension(:),   allocatable, public, protected :: tetra_physics
!
!---------------------------------------------------------------------------------------------------------
!
    type tetrahedron_skew_coord
      sequence
!
      !Position exchange in between tetrahedra via Cartesian coorinates
      double precision, dimension(3,3,4) :: skew_coord_x1x2x3      !Matrix of skew coordinates for all four faces in (coord_system)-coordinates
      double precision, dimension(3,3,4) :: skew_coord_xyz         !Matrix of skew coordinates for all four faces in Cartesian coordinates
      double precision, dimension(3,3,4) :: inv_skew_coord_x1x2x3  !Inverse Matrix of skew coordinates for all four faces in (coord_system)-coordinates
      double precision, dimension(3,3,4) :: inv_skew_coord_xyz     !Inverse Matrix of skew coordinates for all four faces in Cartesian coordinates
      double precision, dimension(3,4)   :: skew_ref_x1x2x3        !Reference vertex for skew coordintes in (coord_system)-coordinates
      double precision, dimension(3,4)   :: skew_ref_xyz           !Reference vertex for skew coordintes in Cartesian coordinates
    end type tetrahedron_skew_coord
!
    type(tetrahedron_skew_coord), dimension(:),   allocatable, public, protected :: tetra_skew_coord
!
!---------------------------------------------------------------------------------------------------------
!
    double precision,public :: cm_over_e,particle_charge,particle_mass,mag_axis_R0,mag_axis_Z0
    integer, public, protected :: coord_system
!   
  contains
!
    subroutine make_tetra_physics(coord_system_in,ipert_in,bmod_multiplier_in)
!
      use tetra_grid_mod, only: tetra_grid,verts_rphiz,verts_sthetaphi,verts_theta_vmec,ntetr,nvert, &
                                & set_verts_sthetaphi,verts_xyz
      use tetra_grid_settings_mod, only: grid_kind,grid_size,n_field_periods
      use constants, only: pi
      use field_mod, only: ipert
      use various_functions_mod, only: dmatinv3
      use gorilla_settings_mod, only: eps_Phi,handover_processing_kind, boole_axi_noise_vector_pot, &
            & boole_axi_noise_elec_pot, boole_non_axi_noise_vector_pot, axi_noise_eps_A, axi_noise_eps_Phi, &
            & non_axi_noise_eps_A
!
      integer, intent(in) :: ipert_in,coord_system_in
      double precision, intent(in),optional :: bmod_multiplier_in
      integer :: iv,i,j,k,l,ind_tetr,ierr
      integer :: nr,nphi,nz,navec,inp_label
      integer,dimension(4) :: vertex_indices
      double precision :: rnd_non_axi_noise_x1,rnd_non_axi_noise_x2,rnd_non_axi_noise_x3, bmod_multiplier,met_det
      double precision, dimension(4)   :: p_x1,p_x2,p_x3
      double precision, dimension(3) :: cur_comp1, cur_comp2, cur_comp3
      double precision, dimension(3) :: curCoordi
      double precision, dimension(:), allocatable   :: davec_dx1,davec_dx2,davec_dx3
      double precision, dimension(:,:),allocatable :: avec
      double precision, dimension(:),           allocatable :: A_x1,A_x2,A_x3
      double precision, dimension(:),           allocatable :: h_x1,h_x2,h_x3,bmod,phi_elec,sqg,dR_ds,dZ_ds
      double precision, dimension(:),         allocatable :: rnd_axi_noise
      double precision :: rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                          ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
      double precision :: s,theta,phi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                          R,Z,alam,dR_ds1,dR_dt,dR_dp,dZ_ds1,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp,q    
      double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_vartheta,Bcovar_varphi
      double precision :: Bcovar_r,B_vartheta,B_varphi
      double precision :: psi_pol,dq_ds,sqrtg,dbmod_dtheta,dR_dtheta
      double precision :: dZ_dtheta,bmod1    
      double precision :: r_minor
!
      !Allocation of quantities dependent on number of vertices
      allocate(tetra_physics(ntetr))
      allocate(A_x1(nvert),A_x2(nvert),A_x3(nvert))
      allocate(bmod(nvert),h_x1(nvert),h_x2(nvert),h_x3(nvert))
      allocate(phi_elec(nvert))  
!
      !Distinguish the Processing of particle handover to tetrahedron neighbour
      select case(handover_processing_kind)
          case(1) !Processing with special treatment of periodic boundaries and manipulation of periodic position values
              allocate(tetra_skew_coord(1)) !Dummy value for parallelization
          case(2) !Position exchange via Cartesian variables (skew_coordinates) - Necessary precomputation is included below
              allocate(tetra_skew_coord(ntetr))
      end select
!
      !Distinguish, in between EFIT and VMEC
      if(grid_kind.eq.3) then !VMEC
        allocate(sqg(nvert))
        navec = 11
      else !EFIT
        navec = 10
      endif
!
      allocate(avec(4,navec))
      allocate(davec_dx1(navec),davec_dx2(navec),davec_dx3(navec))
!
      !Set ipert value in module according to ipert_in
      ipert = ipert_in
!
      !Set coord_system in module according to coord_system_in
      if( ( (grid_kind.eq.1).or.(grid_kind.eq.4) ) .and.(coord_system_in.ne.1)) then
        print *, 'Error: Wrong coordinate system - Only RPhiZ-coordinates are allowed for rectangular or SOLEDGE3X_EIRENE grid.'
        stop
      elseif((grid_kind.eq.3).and.(coord_system_in.ne.2)) then
        print *, 'Error: Wrong coordinate system - Only (s,theta,phi)-coordinates are allowed for field aligned VMEC grid.'
        stop
      else
        coord_system = coord_system_in
      endif 
!
      if(coord_system.eq.2) then
        allocate(dR_ds(nvert))
        allocate(dZ_ds(nvert))
      endif
!
      !Manipulation of the axisymmetric electromagnetic field with noise
      ! - Check that noise is only added in the case of axisymmetric electromagnetic field
      ! - Precompute noise in the case of axisymmetric noise
      if(boole_axi_noise_vector_pot.or.boole_axi_noise_elec_pot.or.boole_non_axi_noise_vector_pot) then
            select case(grid_kind)
                case(1,2)
                    if(boole_axi_noise_vector_pot.or.boole_axi_noise_elec_pot) then
                        allocate(rnd_axi_noise(nvert/grid_size(2)))
                        call random_number(rnd_axi_noise)
                    endif
                case default
                    print *, 'Error: Manipulation of electromagnetic field with noise ONLY for axisymmetric configurations.'
                    stop
            end select
      endif
!
      !Pre-processing for EFIT
!
      if( (grid_kind.eq.1).or.(grid_kind.eq.2).or.(grid_kind.eq.4) ) then
        !subroutine field must be called in order to realize perturbation according to ipert_in
        rrr=1.d0
        ppp=0.d0
        zzz=0.d0
!
        call field(rrr,ppp,zzz,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ,  &      !field must be called in order to obtain Rmin/Rmax and Zmin/Zmax values
                   & dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
      endif
                   
      !Optional input argument bmod_multiplier: Option to calculate magnetic field line (0th order drift orbits)
      if(present(bmod_multiplier_in)) then
        bmod_multiplier = bmod_multiplier_in
      else
        bmod_multiplier = 1.d0  !For normal conditions (drift orbits)
      endif
! 
      !Store cylindrical coordinates and bmod of magnetic axis
      select case(coord_system)
        case(1) !RphiZ --> Cylindrical coordinate system
            print *, 'Magnetic axis for cylindrical coordinates are NOT computed. Magnetic axis is HARD CODED!!!'
            mag_axis_R0 = 170.80689536175018d0
            mag_axis_Z0 = 8.9514398817462517d0
!
        case(2) !sthetaphi --> Symmetry flux coordinate system
            select case(grid_kind)
                case(2) !EFIT (axisymmetric data)
                    !inp_label - input switch: 1 for s and 2 for psi
                    inp_label = 1
                    call magdata_in_symfluxcoord_ext(inp_label,0.d0,psi_pol,2.706d0,q,dq_ds, &
                                         sqrtg,bmod1,dbmod_dtheta,mag_axis_R0,dR_ds1,dR_dtheta,       &
                                         mag_axis_Z0,dZ_ds1,dZ_dtheta)
                case(3) !VMEC (non-axisymmetric)
                    call splint_vmec_data(1.d-16,0.1d0,0.1d0,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota, &
                                            mag_axis_R0,mag_axis_Z0,alam,dR_ds1,dR_dt,dR_dp,dZ_ds1,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
            end select
      end select
!
      !Get field quantities at vertices
      do iv = 1,nvert
!
        select case(coord_system)
          case(1) !RphiZ --> Cylindrical coordinate system
!
            call vector_potential_rphiz(verts_rphiz(1,iv),verts_rphiz(2,iv),verts_rphiz(3,iv),ipert_in,bmod_multiplier, &
                              & A_x1(iv),A_x2(iv),A_x3(iv),h_x1(iv),h_x2(iv),h_x3(iv),Bmod(iv))
!
          case(2) !sthetaphi --> Symmetry flux coordinate system
            select case(grid_kind)
              case(2) !EFIT (axisymmetric data)
                call vector_potential_sthetaphi(verts_sthetaphi(1,iv),verts_sthetaphi(2,iv),verts_sthetaphi(3,iv),&
                                & ipert_in,bmod_multiplier,A_x1(iv),A_x2(iv),A_x3(iv),h_x1(iv),h_x2(iv),h_x3(iv),&
                                & Bmod(iv),q,dR_ds(iv),dZ_ds(iv))
!
              case(3) !VMEC (non-axisymmetric)
                !Read VMEC
                call vector_potential_sthetaphi_vmec(verts_sthetaphi(1,iv),verts_theta_vmec(iv),verts_sthetaphi(3,iv), &
                                  & ipert_in,bmod_multiplier,A_x1(iv),A_x2(iv),A_x3(iv),h_x1(iv),h_x2(iv),h_x3(iv),Bmod(iv), &
                                  & sqg(iv),dR_ds(iv),dZ_ds(iv))  
!               
            end select  
!            
        end select
!
!
!
        !Optionally add axisymmetric noise to vector potential
        if(boole_axi_noise_vector_pot) then
            A_x1(iv) = A_x1(iv)+A_x1(iv)*axi_noise_eps_A*rnd_axi_noise(modulo(iv-1,(nvert/grid_size(2))) +1)
            A_x2(iv) = A_x2(iv)+A_x2(iv)*axi_noise_eps_A*rnd_axi_noise(modulo(iv-1,(nvert/grid_size(2))) +1)
            A_x3(iv) = A_x3(iv)+A_x3(iv)*axi_noise_eps_A*rnd_axi_noise(modulo(iv-1,(nvert/grid_size(2))) +1)
        endif

        !Optionally add non-axisymmetric noise to vector potential
        if(boole_non_axi_noise_vector_pot) then
            call random_number(rnd_non_axi_noise_x1)
            call random_number(rnd_non_axi_noise_x2)
            call random_number(rnd_non_axi_noise_x3)
            A_x1(iv) = A_x1(iv)+A_x1(iv)*non_axi_noise_eps_A*rnd_non_axi_noise_x1
            A_x2(iv) = A_x2(iv)+A_x2(iv)*non_axi_noise_eps_A*rnd_non_axi_noise_x2
            A_x3(iv) = A_x3(iv)+A_x3(iv)*non_axi_noise_eps_A*rnd_non_axi_noise_x3
        endif
!
        !Unit vector in the direction of the magnetic field
        h_x1(iv)=h_x1(iv)/Bmod(iv)
        h_x2(iv)=h_x2(iv)/Bmod(iv)
        h_x3(iv)=h_x3(iv)/Bmod(iv)
!
        !Electrostatic potential as a product of co-variant component of vector potential and a factor eps_Phi (gorilla.inp)
        phi_elec(iv) = A_x2(iv)*eps_Phi
!
        !Optionally add axisymmetric noise to electrostatic potential
        if(boole_axi_noise_elec_pot) then
            phi_elec(iv) = phi_elec(iv)+phi_elec(iv)*axi_noise_eps_Phi*rnd_axi_noise(modulo(iv-1,(nvert/grid_size(2))) +1)
        endif
!
      enddo !iv (index vertex)
!
  !$OMP PARALLEL &
  !$OMP& DO DEFAULT(NONE) &
  !$OMP& SHARED(ntetr,tetra_grid,coord_system,verts_rphiz,grid_kind,verts_sthetaphi, &
  !$OMP& grid_size,A_x1,A_x2,A_x3,h_x1,h_x2,h_x3,Bmod,phi_elec,handover_processing_kind, &
  !$OMP& sqg,tetra_physics,tetra_skew_coord,navec,verts_xyz,mag_axis_R0,mag_axis_Z0,dR_ds,dZ_ds,n_field_periods) &
  !$OMP& PRIVATE(ind_tetr,i,j,k,l,iv,p_x1,p_x2,p_x3,A_phi,A_theta,dA_phi_ds, &
  !$OMP& dA_theta_ds,aiota,R,Z,alam,dR_ds1,dR_dt,dR_dp,dZ_ds1,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp, &
  !$OMP& avec,cur_comp1,cur_comp2,cur_comp3,curCoordi,davec_dx1,davec_dx2,davec_dx3,met_det, &
  !$OMP& vertex_indices,ierr,r_minor)
!
  do ind_tetr=1,ntetr
!
    do i=1,4
      iv=tetra_grid(ind_tetr)%ind_knot(i)
!
      select case(coord_system)
        case(1) !R,Phi,Z --> Cylindrical coordinate system
          p_x1(i)=verts_rphiz(1,iv)
          p_x2(i)=verts_rphiz(2,iv)
          p_x3(i)=verts_rphiz(3,iv)
        case(2) !s,theta,phi --> Symmetry flux coordinate system
          select case(grid_kind)
            case(2) !EFIT (axisymmetric data)
              p_x1(i)=verts_sthetaphi(1,iv)
              p_x2(i)=verts_sthetaphi(2,iv)
              p_x3(i)=verts_sthetaphi(3,iv)
!              
            case(3) !VMEC (non-axisymmetric)
              p_x1(i)=verts_sthetaphi(1,iv)
              p_x2(i)=verts_sthetaphi(2,iv)
              p_x3(i)=verts_sthetaphi(3,iv)
!
          end select    
      end select
!
      !Consider periodic boundary conditions
      select case(grid_kind)
        case(1) !rectangular grid
          !nothing to change
        case(2) !field-aligned grid EFIT
          select case(coord_system)
            case(1) !R,Phi,Z --> Cylindrical coordinate system
              if ((ind_tetr .ge. ntetr-ntetr/grid_size(2)+1) .and. (verts_rphiz(2,iv) .eq. 0.d0)) then
                p_x2(i) = 2.d0*pi
              endif
            case(2) !s,theta,phi --> Symmetry flux coordinate system
              if ((ind_tetr .ge. ntetr-ntetr/grid_size(2)+1) .and. (verts_sthetaphi(3,iv) .eq. 0.d0)) then !if(tetra in last phi_slice & phi .eq. 0)
                p_x3(i) = 2.d0*pi/n_field_periods
              endif
          end select
        case(3) !VMEC field-aligned grid in Symmetry flux coordinate system
          if ((ind_tetr .ge. ntetr-ntetr/grid_size(2)+1) .and. (verts_sthetaphi(3,iv) .eq. 0.d0)) then !if(tetra in last phi_slice & phi .eq. 0)
            p_x3(i) = 2.d0*pi/n_field_periods
          endif
        case(4) !SOLEDGE3X_EIRENE grid
                !R,Phi,Z --> Cylindrical coordinate system
                if ((ind_tetr .ge. ntetr-ntetr/grid_size(2)+1) .and. (verts_rphiz(2,iv) .eq. 0.d0)) then
                p_x2(i) = 2.d0*pi
                endif  
!
      end select
!
      avec(i,1)=A_x1(iv)
      avec(i,2)=A_x2(iv)
      avec(i,3)=A_x3(iv)
      avec(i,4)=h_x1(iv)
      avec(i,5)=h_x2(iv)
      avec(i,6)=h_x3(iv)
      avec(i,7)=Bmod(iv)
      avec(i,8)=phi_elec(iv)
      avec(i,9)=verts_rphiz(1,iv) !major radius R
      avec(i,10)=verts_rphiz(3,iv) !Z
!      
      if(grid_kind.eq.3) avec(i,11)=sqg(iv)
    enddo
!
    !At least one vertex lies on theta = 2pi, but reference point in verts-matrix is 0
    if(coord_system.eq.2) then
        do j = 1,4
            if( (p_x2(j).eq.0.d0) .and. any(p_x2.ge.pi)) then
                p_x2(j) = 2.d0*pi
            endif
        enddo
    endif
!
    ! Calculate geometric properties (x1,anorm,dist_ref)
    
       ! (1) Coordinates of first vertex within each tetrahedron
        tetra_physics(ind_tetr)%x1(:) = [p_x1(1),p_x2(1),p_x3(1)]
!
        ! (2.1) Major radius R of first vertex within each tetrahedron
        tetra_physics(ind_tetr)%R1 = verts_rphiz(1,tetra_grid(ind_tetr)%ind_knot(1))
!
        ! (2.2) Z of first vertex within each tetrahedron
        tetra_physics(ind_tetr)%Z1 = verts_rphiz(3,tetra_grid(ind_tetr)%ind_knot(1))
!
        ! (3) Normal vectors anorm (pointing inwards) and dist_ref
        do i=1,4 !Iterating over all 4 faces for which i want to calculate the normal vectors
!
          k = 0
          do j = 1,4
            if(j.eq.i) cycle
            k = k+1
            cur_comp1(k) = p_x1(j)
            cur_comp2(k) = p_x2(j)
            cur_comp3(k) = p_x3(j)
          enddo
!
          !calculate the vectors in the plane
          cur_comp1(1:2)=cur_comp1(1:2)-cur_comp1(3)
          cur_comp2(1:2)=cur_comp2(1:2)-cur_comp2(3)
          cur_comp3(1:2)=cur_comp3(1:2)-cur_comp3(3)
!
          !calculate the cross product of these 2 vectors
          tetra_physics(ind_tetr)%anorm(1,i)=cur_comp2(1)*cur_comp3(2)-cur_comp2(2)*cur_comp3(1)
          tetra_physics(ind_tetr)%anorm(2,i)=cur_comp3(1)*cur_comp1(2)-cur_comp3(2)*cur_comp1(1)
          tetra_physics(ind_tetr)%anorm(3,i)=cur_comp1(1)*cur_comp2(2)-cur_comp1(2)*cur_comp2(1)
!
          !Coordinates of the current i-th vertex opposing the current face i
          curCoordi = [p_x1(i),p_x2(i),p_x3(i)]
!
          tetra_physics(ind_tetr)%dist_ref_vec(i)                           &
                  =tetra_physics(ind_tetr)%anorm(1,i)*(curCoordi(1)-cur_comp1(3)) &
                  +tetra_physics(ind_tetr)%anorm(2,i)*(curCoordi(2)-cur_comp2(3)) &
                  +tetra_physics(ind_tetr)%anorm(3,i)*(curCoordi(3)-cur_comp3(3))
          if(tetra_physics(ind_tetr)%dist_ref_vec(i).lt.0.d0) then
            tetra_physics(ind_tetr)%dist_ref_vec(i)=-tetra_physics(ind_tetr)%dist_ref_vec(i)
            tetra_physics(ind_tetr)%anorm(:,i)=-tetra_physics(ind_tetr)%anorm(:,i)
          endif
!
        enddo
        !Hard coded 4th plane
        tetra_physics(ind_tetr)%dist_ref = tetra_physics(ind_tetr)%dist_ref_vec(4)
!
    !module of B in the first vertex:
    tetra_physics(ind_tetr)%bmod1 = avec(1,7)
!
    !phi-component of the vector potential A in the first vertex
    select case(coord_system)
        case(1)
            tetra_physics(ind_tetr)%Aphi1 = avec(1,2)
        case(2)
            tetra_physics(ind_tetr)%Aphi1 = avec(1,3)
    end select
!
! 2nd component of the unit vector h in the first vertex
    tetra_physics(ind_tetr)%h2_1 = avec(1,5)
!
! 3rd component of the unit vector h in the first vertex
    tetra_physics(ind_tetr)%h3_1 = avec(1,6)
!
! Electrostatic Potential Phi in the first vertex:
    tetra_physics(ind_tetr)%Phi1 = avec(1,8)
!
! Square root g in the first vertex:
    if(grid_kind.eq.3) tetra_physics(ind_tetr)%sqg1 = avec(1,11)
! 
! derivatives:
    call differentiate(p_x1,p_x2,p_x3,navec,avec,davec_dx1,davec_dx2,davec_dx3)
!
! "curl(A)" -  $\epsilon^{ijk}\difp{A_k}{x^j}$ :
    tetra_physics(ind_tetr)%curlA(1) = davec_dx2(3)-davec_dx3(2)
    tetra_physics(ind_tetr)%curlA(2) = davec_dx3(1)-davec_dx1(3)
    tetra_physics(ind_tetr)%curlA(3) = davec_dx1(2)-davec_dx2(1)
!
! gradient of B times "curl(A)" - $\epsilon^{ijk}difp{B}{x^i}\difp{A_k}{x^j} :
    tetra_physics(ind_tetr)%gBxcurlA = davec_dx1(7)*tetra_physics(ind_tetr)%curlA(1) &
                             + davec_dx2(7)*tetra_physics(ind_tetr)%curlA(2) &
                             + davec_dx3(7)*tetra_physics(ind_tetr)%curlA(3)
!
! gradient of Phi times "curl(A)" - $\epsilon^{ijk}difp{Phi}{x^i}\difp{A_k}{x^j} :
    tetra_physics(ind_tetr)%gPhixcurlA = davec_dx1(8)*tetra_physics(ind_tetr)%curlA(1) &
                             + davec_dx2(8)*tetra_physics(ind_tetr)%curlA(2) &
                             + davec_dx3(8)*tetra_physics(ind_tetr)%curlA(3)
!
! gradient of B module - $\difp{B}{x^i} :
    tetra_physics(ind_tetr)%gB(1) = davec_dx1(7)
    tetra_physics(ind_tetr)%gB(2) = davec_dx2(7)
    tetra_physics(ind_tetr)%gB(3) = davec_dx3(7)
!
! gradient of Phi - $\difp{Phi}{x^i} :
    tetra_physics(ind_tetr)%gPhi(1) = davec_dx1(8)
    tetra_physics(ind_tetr)%gPhi(2) = davec_dx2(8)
    tetra_physics(ind_tetr)%gPhi(3) = davec_dx3(8)
!
! gradient of R - $\difp{R}{x^i} :
    tetra_physics(ind_tetr)%gR(1) = davec_dx1(9)
    tetra_physics(ind_tetr)%gR(2) = davec_dx2(9)
    tetra_physics(ind_tetr)%gR(3) = davec_dx3(9)
!
! gradient of Z - $\difp{Z}{x^i} :
    tetra_physics(ind_tetr)%gZ(1) = davec_dx1(10)
    tetra_physics(ind_tetr)%gZ(2) = davec_dx2(10)
    tetra_physics(ind_tetr)%gZ(3) = davec_dx3(10)
!
! gradient of sqrtg - $\difp{\sqrt(g)}{x^i} :
    if(grid_kind.eq.3) then
      tetra_physics(ind_tetr)%gsqg(1) = davec_dx1(11)
      tetra_physics(ind_tetr)%gsqg(2) = davec_dx2(11)
      tetra_physics(ind_tetr)%gsqg(3) = davec_dx3(11)
    endif  
!
! gradient of Aphi - $\difp{Aphi}{x^i} :
    select case(coord_system)
        case(1)
            tetra_physics(ind_tetr)%gAphi(1) = davec_dx1(2)
            tetra_physics(ind_tetr)%gAphi(2) = davec_dx2(2)
            tetra_physics(ind_tetr)%gAphi(3) = davec_dx3(2)
        case(2)
            tetra_physics(ind_tetr)%gAphi(1) = davec_dx1(3)
            tetra_physics(ind_tetr)%gAphi(2) = davec_dx2(3)
            tetra_physics(ind_tetr)%gAphi(3) = davec_dx3(3)
    end select
!
! gradient of h1 - $\difp{h1}{x^i} :
    tetra_physics(ind_tetr)%gh1(1) = davec_dx1(4)
    tetra_physics(ind_tetr)%gh1(2) = davec_dx2(4)
    tetra_physics(ind_tetr)%gh1(3) = davec_dx3(4)
!
! gradient of h2 - $\difp{h2}{x^i} :
    tetra_physics(ind_tetr)%gh2(1) = davec_dx1(5)
    tetra_physics(ind_tetr)%gh2(2) = davec_dx2(5)
    tetra_physics(ind_tetr)%gh2(3) = davec_dx3(5)
!
! gradient of h3 - $\difp{h3}{x^i} :
    tetra_physics(ind_tetr)%gh3(1) = davec_dx1(6)
    tetra_physics(ind_tetr)%gh3(2) = davec_dx2(6)
    tetra_physics(ind_tetr)%gh3(3) = davec_dx3(6)
!
! "curl(h)" -  $\epsilon^{ijk}\difp{h_k}{x^j}$ :
    tetra_physics(ind_tetr)%curlh(1) = davec_dx2(6)-davec_dx3(5)
    tetra_physics(ind_tetr)%curlh(2) = davec_dx3(4)-davec_dx1(6)
    tetra_physics(ind_tetr)%curlh(3) = davec_dx1(5)-davec_dx2(4)
!
! "vector product" of grad B times vector h in the first vertex:
    tetra_physics(ind_tetr)%gBxh1(1) = davec_dx2(7)*avec(1,6)-davec_dx3(7)*avec(1,5)
    tetra_physics(ind_tetr)%gBxh1(2) = davec_dx3(7)*avec(1,4)-davec_dx1(7)*avec(1,6)
    tetra_physics(ind_tetr)%gBxh1(3) = davec_dx1(7)*avec(1,5)-davec_dx2(7)*avec(1,4)
!
! "vector product" of grad Phi times vector h in the first vertex:
    tetra_physics(ind_tetr)%gPhixh1(1) = davec_dx2(8)*avec(1,6)-davec_dx3(8)*avec(1,5)
    tetra_physics(ind_tetr)%gPhixh1(2) = davec_dx3(8)*avec(1,4)-davec_dx1(8)*avec(1,6)
    tetra_physics(ind_tetr)%gPhixh1(3) = davec_dx1(8)*avec(1,5)-davec_dx2(8)*avec(1,4)
!
! real space part of matrix alpha:
    tetra_physics(ind_tetr)%alpmat(1,1) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gB(1) &
                                + davec_dx2(7)*davec_dx1(6)-davec_dx3(7)*davec_dx1(5)
    tetra_physics(ind_tetr)%alpmat(1,2) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gB(2) &
                                + davec_dx2(7)*davec_dx2(6)-davec_dx3(7)*davec_dx2(5)
    tetra_physics(ind_tetr)%alpmat(1,3) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gB(3) &
                                + davec_dx2(7)*davec_dx3(6)-davec_dx3(7)*davec_dx3(5)
!
    tetra_physics(ind_tetr)%alpmat(2,1) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gB(1) &
                                + davec_dx3(7)*davec_dx1(4)-davec_dx1(7)*davec_dx1(6)
    tetra_physics(ind_tetr)%alpmat(2,2) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gB(2) &
                                + davec_dx3(7)*davec_dx2(4)-davec_dx1(7)*davec_dx2(6)
    tetra_physics(ind_tetr)%alpmat(2,3) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gB(3) &
                                + davec_dx3(7)*davec_dx3(4)-davec_dx1(7)*davec_dx3(6)
!
    tetra_physics(ind_tetr)%alpmat(3,1) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gB(1) &
                                + davec_dx1(7)*davec_dx1(5)-davec_dx2(7)*davec_dx1(4)
    tetra_physics(ind_tetr)%alpmat(3,2) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gB(2) &
                                + davec_dx1(7)*davec_dx2(5)-davec_dx2(7)*davec_dx2(4)
    tetra_physics(ind_tetr)%alpmat(3,3) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gB(3) &
                                + davec_dx1(7)*davec_dx3(5)-davec_dx2(7)*davec_dx3(4)
!
! real space part of matrix beta:
    tetra_physics(ind_tetr)%betmat(1,1) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gPhi(1) &
                                + davec_dx2(8)*davec_dx1(6)-davec_dx3(8)*davec_dx1(5)
    tetra_physics(ind_tetr)%betmat(1,2) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gPhi(2) &
                                + davec_dx2(8)*davec_dx2(6)-davec_dx3(8)*davec_dx2(5)
    tetra_physics(ind_tetr)%betmat(1,3) = 2.d0*tetra_physics(ind_tetr)%curlh(1)*tetra_physics(ind_tetr)%gPhi(3) &
                                + davec_dx2(8)*davec_dx3(6)-davec_dx3(8)*davec_dx3(5)
!
    tetra_physics(ind_tetr)%betmat(2,1) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gPhi(1) &
                                + davec_dx3(8)*davec_dx1(4)-davec_dx1(8)*davec_dx1(6)
    tetra_physics(ind_tetr)%betmat(2,2) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gPhi(2) &
                                + davec_dx3(8)*davec_dx2(4)-davec_dx1(8)*davec_dx2(6)
    tetra_physics(ind_tetr)%betmat(2,3) = 2.d0*tetra_physics(ind_tetr)%curlh(2)*tetra_physics(ind_tetr)%gPhi(3) &
                                + davec_dx3(8)*davec_dx3(4)-davec_dx1(8)*davec_dx3(6)
!
    tetra_physics(ind_tetr)%betmat(3,1) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gPhi(1) &
                                + davec_dx1(8)*davec_dx1(5)-davec_dx2(8)*davec_dx1(4)
    tetra_physics(ind_tetr)%betmat(3,2) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gPhi(2) &
                                + davec_dx1(8)*davec_dx2(5)-davec_dx2(8)*davec_dx2(4)
    tetra_physics(ind_tetr)%betmat(3,3) = 2.d0*tetra_physics(ind_tetr)%curlh(3)*tetra_physics(ind_tetr)%gPhi(3) &
                                + davec_dx1(8)*davec_dx3(5)-davec_dx2(8)*davec_dx3(4)
!
! trace of real space part of matrix alpha (element $\alpha^4_4$):
    tetra_physics(ind_tetr)%spalpmat = tetra_physics(ind_tetr)%alpmat(1,1) &
                             + tetra_physics(ind_tetr)%alpmat(2,2) &
                             + tetra_physics(ind_tetr)%alpmat(3,3)
!
! trace of real space part of matrix beta (element $\beta^4_4$):
    tetra_physics(ind_tetr)%spbetmat = tetra_physics(ind_tetr)%betmat(1,1) &
                             + tetra_physics(ind_tetr)%betmat(2,2) &
                             + tetra_physics(ind_tetr)%betmat(3,3)   
!
! precomputation of acef_pre for analytical quadratic approximation
    tetra_physics(ind_tetr)%acoef_pre = matmul(tetra_physics(ind_tetr)%curlA,tetra_physics(ind_tetr)%anorm)
!
! Constant dt_dtau averaged over all 4 vertices of tetrahedron
    tetra_physics(ind_tetr)%dt_dtau_const = 0.d0
    do j = 1,4
        ! Calculate metric_determinant dependening on coordinate system and vertex j
        select case(coord_system)
          case(1) !R,Phi,Z --> Cylindrical coordinate system
            met_det = p_x1(j)
          case(2) !s,theta,phi --> Symmetry flux coordinate system 
!
            select case(grid_kind)
              case(2) !EFIT -> Use linearized quantity
                met_det = metric_determinant(ind_tetr,[p_x1(j),p_x2(j),p_x3(j)])
!                            
              case(3) !VMEC
                met_det = avec(j,11)
!                                             
            end select  
        end select
!    
        !dt_dtau = metric_determinant * bmod 
        tetra_physics(ind_tetr)%dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const + met_det*  avec(j,7)
!
    enddo
    !Average over all 4 vertices
    tetra_physics(ind_tetr)%dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const/4.d0
!
    !Constant electric field modulus in r-direction averaged over all four vertices
    tetra_physics(ind_tetr)%Er_mod = 0.d0
    select case(coord_system)
        case(1) !R,Phi,Z --> Cylindrical coordinate system
            do j = 1,4
                iv=tetra_grid(ind_tetr)%ind_knot(j)
!                
                !Minor radius at the position of the vertex Sqrt((R-R0)^2 + (Z-Z0)^2)
                r_minor = sqrt((avec(j,9)-mag_axis_R0)**2+(avec(j,10)-mag_axis_Z0)**2)
!                
                tetra_physics(ind_tetr)%Er_mod = tetra_physics(ind_tetr)%Er_mod + & !dPhi_dR*dR_dr + dPhi_dZ*dZ_dr
                                & (tetra_physics(ind_tetr)%gPhi(1) * r_minor/(avec(j,9)-mag_axis_R0) + &
                                &  tetra_physics(ind_tetr)%gPhi(3) * r_minor/(avec(j,10)-mag_axis_Z0))
            enddo
        case(2) !s,theta,phi --> Symmetry flux coordinate system 
            do j = 1,4
                iv=tetra_grid(ind_tetr)%ind_knot(j)
                tetra_physics(ind_tetr)%Er_mod = tetra_physics(ind_tetr)%Er_mod + & !dPhi_ds * ds_dr
                                & tetra_physics(ind_tetr)%gPhi(1) * &
                                & sqrt((avec(j,9)-mag_axis_R0)**2+(avec(j,10)-mag_axis_Z0)**2) &
                                & /((avec(j,9)-mag_axis_R0)*dR_ds(iv) + (avec(j,10)-mag_axis_Z0)*dZ_ds(iv))  
            enddo
    end select
    tetra_physics(ind_tetr)%Er_mod = abs(tetra_physics(ind_tetr)%Er_mod/4.d0)

! Tetrahedron reference distance in the first vertex
    tetra_physics(ind_tetr)%tetra_dist_ref = 2.d0*pi/n_field_periods*tetra_physics(ind_tetr)%R1/grid_size(2)
!
!
        !-------------- Precomputation of matrices for position exchange in between tetrahedra via Cartesian coorinates ---------!
!
        !Precomputation is only performed, if position exchange via Cartesian variables (skew_coordinates) is selected
        if(handover_processing_kind.eq.2) then
!
            !Loop over all four 'ifaces'
            do k = 1,4
!
              !The indices of the vertices must be chosen in such a way, that iface is not the reference point.
              !Meaning, that two of the skew coord-vectors MUST lie in the plane of the iface.
              !vertex_indices(1) is ALWAYS the reference vertex and it takes never the value of iface
              do l = 1,4
                vertex_indices(l) = modulo(k+l-1,4) + 1
              enddo
!
              !Reference vertex for skew coordinates in (coord_system)-coordinates
              tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(1,k) = p_x1(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(2,k) = p_x2(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_ref_x1x2x3(3,k) = p_x3(vertex_indices(1))
!
              !Matrix of skew coordinates for all four faces in (coord_system)-coordinates
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(1,1,k) = p_x1(vertex_indices(2))-p_x1(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(1,2,k) = p_x1(vertex_indices(3))-p_x1(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(1,3,k) = p_x1(vertex_indices(4))-p_x1(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(2,1,k) = p_x2(vertex_indices(2))-p_x2(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(2,2,k) = p_x2(vertex_indices(3))-p_x2(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(2,3,k) = p_x2(vertex_indices(4))-p_x2(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(3,1,k) = p_x3(vertex_indices(2))-p_x3(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(3,2,k) = p_x3(vertex_indices(3))-p_x3(vertex_indices(1))
              tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(3,3,k) = p_x3(vertex_indices(4))-p_x3(vertex_indices(1))
!
              !Compute inverse of matrix of skew coordinates for all four faces in (coord_system)-coordinates
              call dmatinv3(tetra_skew_coord(ind_tetr)%skew_coord_x1x2x3(:,:,k), &
                            & tetra_skew_coord(ind_tetr)%inv_skew_coord_x1x2x3(:,:,k),ierr)
!
              !Reference vertex for skew coordinates in Cartesian coordinates
              tetra_skew_coord(ind_tetr)%skew_ref_xyz(1,k) = verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_ref_xyz(2,k) = verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_ref_xyz(3,k) = verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
!
              !Matrix of skew coordinates for all four faces in Cartesian coordinates
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(1,1,k) = verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(2)))-&
                                          & verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(1,2,k) = verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(3)))-&
                                          & verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(1,3,k) = verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(4)))-&
                                          & verts_xyz(1,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(2,1,k) = verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(2)))-&
                                          & verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(2,2,k) = verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(3)))-&
                                          & verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(2,3,k) = verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(4)))-&
                                          & verts_xyz(2,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(3,1,k) = verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(2)))-&
                                          & verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(3,2,k) = verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(3)))-&
                                          & verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
              tetra_skew_coord(ind_tetr)%skew_coord_xyz(3,3,k) = verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(4)))-&
                                          & verts_xyz(3,tetra_grid(ind_tetr)%ind_knot(vertex_indices(1)))
!
              !Compute inverse of matrix of skew coordinates for all four faces in Cartesian coordinates
              call dmatinv3(tetra_skew_coord(ind_tetr)%skew_coord_xyz(:,:,k), &
                            & tetra_skew_coord(ind_tetr)%inv_skew_coord_xyz(:,:,k),ierr)
!
            enddo !Loop over all four 'ifaces'
!
        endif !handover_processing_kind.eq.2
!
enddo
!$OMP END PARALLEL DO
!
        !Deallocate quantities
        deallocate(davec_dx1,davec_dx2,davec_dx3)
        deallocate(avec)
        deallocate(A_x1,A_x2,A_x3)
        deallocate(h_x1,h_x2,h_x3,bmod,phi_elec)
        if(coord_system.eq.2) deallocate(dR_ds,dZ_ds)
        if(grid_kind.eq.3) deallocate(sqg)
        if(boole_axi_noise_vector_pot.or.boole_axi_noise_elec_pot) deallocate(rnd_axi_noise)
!
      end subroutine make_tetra_physics
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      function isinside(ind_tetr,current_x) result(bool_isinside)
!
        use constants, only: eps
!
        implicit none
!
        integer, intent(in) :: ind_tetr
        double precision,dimension(3), intent(in):: current_x
!
        logical :: bool_isinside
!
        integer:: ind_normdist, i
        double precision :: dist_min
        double precision,dimension(4) :: cur_dist_value
!
        !It is needed to understand, if vector is inside the tetrahedron
        dist_min=eps*abs(tetra_physics(ind_tetr)%dist_ref)
        bool_isinside = .false.
!
        do ind_normdist = 1,4
          if (ind_normdist .ne. 1) then
            cur_dist_value(ind_normdist) = sum(tetra_physics(ind_tetr)%anorm(:,ind_normdist)*(current_x-tetra_physics(ind_tetr)%x1))
          else ! ind_normdist .eq. 1
            cur_dist_value(ind_normdist) = sum(tetra_physics(ind_tetr)%anorm(:,ind_normdist)*&
            &(current_x-tetra_physics(ind_tetr)%x1))+tetra_physics(ind_tetr)%dist_ref
          endif
        enddo
!
        if (all(cur_dist_value .ge. 0.d0)) then
          bool_isinside  = .true.
        else
          bool_isinside = .false.
        endif

      end function isinside
! 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine vector_potential_rphiz(r,phi,z,ipert,bmod_multiplier,A_r,A_phi,A_z,B_r,B_phi,B_z,bmod)
!
        use field_eq_mod,                 only : rtf,btf,psif
        use field_mod,                    only : ampl
        use getout_vector_potentials_mod, only : ar,az

!
        implicit none
!
        integer, intent(in) :: ipert
        double precision, intent(in) :: bmod_multiplier
        double precision :: r,phi,z,A_r,A_phi,A_z,B_phi,bmod
        double precision :: B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
        call field(r,phi,z,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
        bmod=sqrt(B_r**2+B_p**2+B_z**2)*bmod_multiplier
!
        A_r=0.d0
        A_phi=psif
        A_z=-rtf*btf*log(r)
!
        B_phi=B_p*R
!
        if(ipert.gt.0) then
          A_r=A_r+ar*ampl
          A_z=A_z+az*ampl
          !A_phi = A_phi + eps * cos(m*theta +n phi)
        endif
!
      end subroutine vector_potential_rphiz
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine vector_potential_sthetaphi(s,theta,phi,ipert,bmod_multiplier,A_s,A_theta,A_phi,B_s,B_theta,B_phi,bmod, &
                                            & q,dR_ds,dZ_ds)
!
        use field_eq_mod,                 only : rtf,btf,psif
        use field_mod,                    only : ampl
        use getout_vector_potentials_mod, only : ar,az
        use magdata_in_symfluxcoor_mod,   only : psitor_max
        use gorilla_settings_mod,         only : boole_helical_pert, helical_pert_eps_Aphi, helical_pert_m_fourier, &
                                               & helical_pert_n_fourier
!
        implicit none
!
        integer, intent(in) :: ipert
        integer :: inp_label,m,n
        double precision, intent(in) :: bmod_multiplier
        double precision :: s,theta,phi,A_s,A_theta,A_phi,B_phi,bmod,bmod1
        double precision :: psi_pol,q,dq_ds,sqrtg,dbmod_dtheta,R,dR_ds,dR_dtheta
        double precision :: Z, dZ_ds, dZ_dtheta,B_s,B_theta,B_theta1
        double precision :: B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                            ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
        !inp_label - input switch: 1 for s and 2 for psi
        inp_label = 1
        call magdata_in_symfluxcoord_ext(inp_label,s,psi_pol,theta,q,dq_ds, &
                                         sqrtg,bmod1,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
                                         Z,dZ_ds,dZ_dtheta)
!
        call field(R,phi,Z,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ  &
                  ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
        bmod=sqrt(B_r**2+B_p**2+B_z**2)*bmod_multiplier
!
        !Calculate covarient components of the magnetic field B
        B_phi=B_p*R
        B_s = B_r*dR_ds + B_z*dZ_ds
        B_theta = B_r*dR_dtheta + B_z*dZ_dtheta
        B_theta1 = (dR_dtheta**2+dZ_dtheta**2)*B_p/(q*R)
!
        !Calculate vector potentials
        A_s=0.d0
        A_theta = s * psitor_max
        A_phi= psi_pol !Eventuell Vorzeichen
!
        !Analytical perturbation
        if(boole_helical_pert) then
            A_phi = A_phi + A_phi*helical_pert_eps_Aphi * cos(helical_pert_m_fourier*theta +helical_pert_n_fourier*phi)
        endif
!
      end subroutine vector_potential_sthetaphi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine vector_potential_sthetaphi_vmec(s,theta,phi,ipert,bmod_multiplier,A_s,A_theta,A_phi, &
                                                 & B_s,B_vartheta,B_varphi,bmod,sqg,dR_ds,dZ_ds)
!
        implicit none
!
        integer, intent(in) :: ipert
        double precision, intent(in) :: bmod_multiplier
        double precision, intent(in) :: s,theta,phi
        double precision, intent(out) :: bmod
        double precision :: A_s, A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,        &
                            alam,dl_ds,dl_dt,dl_dp
        double precision :: Bctrvr_vartheta,Bctrvr_varphi,Bcovar_vartheta,Bcovar_varphi,sqg
        double precision :: Bcovar_r,B_s,B_vartheta,B_varphi
!
        !VMEC field quantities
        double precision :: A_s_V,A_theta_V,A_phi_V,Bcovar_s_V,Bcovar_theta_V,Bcovar_phi_V
        double precision :: R,Z,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp
!       
        call splint_vmec_data(s,theta,phi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                        R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
        call vmec_field(s,theta,phi,A_theta,A_phi,dA_theta_ds,dA_phi_ds,aiota,     &
                        sqg,alam,dl_ds,dl_dt,dl_dp,Bctrvr_vartheta,Bctrvr_varphi,     &
                        Bcovar_r,Bcovar_vartheta,Bcovar_varphi)
!
        !Squareroot(g) is for some reason a negative quantity ???
        sqg = abs(sqg)

        !Physical component of B-Field
        bmod = sqrt(Bctrvr_vartheta*Bcovar_vartheta + Bctrvr_varphi*Bcovar_varphi)*bmod_multiplier
!
        !Covariant s-component of vector potential in symmetry flux coordinates
        A_s = 0.d0
!
        B_s = Bcovar_r
        B_vartheta = Bcovar_vartheta
        B_varphi = Bcovar_varphi
!
      end subroutine vector_potential_sthetaphi_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      function psi_rphiz(ind_tetr,x) result(psi_linear)
!
        implicit none
!
        integer :: ind_tetr
        double precision :: psi_linear
        double precision, dimension(3) :: x
!
        psi_linear = tetra_physics(ind_tetr)%Aphi1 + &
                   & tetra_physics(ind_tetr)%curlA(3)*(x(1)-tetra_physics(ind_tetr)%x1(1)) - &
                   & tetra_physics(ind_tetr)%curlA(1)*(x(3)-tetra_physics(ind_tetr)%x1(3))
!
      end function psi_rphiz
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      function dt_dtau(ind_tetr,x_enter,x_exit)
!
        implicit none
!
        integer, intent(in ) :: ind_tetr
        double precision :: dt_dtau,met_det_xenter,met_det_xexit
        double precision, dimension(3), intent(in) :: x_enter,x_exit
        double precision, dimension(3) :: x1
!
        ! This function calculates the constant factor dt_dtau for a particle orbit inside a tetrahedron,
        ! while considering the points, where the particle enters and exits the tetrahedron.
        ! The metric tensor sqrtg is adapted to the kind of coordinate system.
!
        !Coordinates at first vertex
        x1 = tetra_physics(ind_tetr)%x1
!
        !Metric determinants at entering and exit point
        met_det_xenter = metric_determinant(ind_tetr,x_enter)
        met_det_xexit = metric_determinant(ind_tetr,x_exit)
        dt_dtau = 0.5d0 * ( tetra_physics(ind_tetr)%bmod1*(met_det_xenter + met_det_xexit) &
                        & + sum(tetra_physics(ind_tetr)%gB * ((x_enter-x1)*met_det_xenter + (x_exit-x1)*met_det_xexit) ) )
!
      end function dt_dtau
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      function metric_determinant(ind_tetr,x)
!
        use magdata_in_symfluxcoor_mod,   only : psitor_max
        use tetra_grid_settings_mod, only: grid_kind
!
        implicit none
!
        integer :: ind_tetr
        double precision :: metric_determinant
        double precision, dimension(3),intent(in) :: x
        double precision, dimension(3) :: x1
!
        ! Calculate metric_determinant dependening on coordinate system and position
        select case(coord_system)
          case(1) !R,Phi,Z --> Cylindrical coordinate system
            metric_determinant = x(1)
          case(2) !s,theta,phi --> Symmetry flux coordinate system
            x1 = tetra_physics(ind_tetr)%x1
!
            select case(grid_kind)
              case(2) !EFIT
                metric_determinant = ( (tetra_physics(ind_tetr)%R1+sum(tetra_physics(ind_tetr)%gR*(x-x1)) )**2 * psitor_max ) / &
                                       & ( (tetra_physics(ind_tetr)%h3_1 + sum(tetra_physics(ind_tetr)%gh3*(x-x1))) * &
                                       &   (tetra_physics(ind_tetr)%bmod1 + sum(tetra_physics(ind_tetr)%gB*(x-x1))) )
                            
              case(3) !VMEC
                metric_determinant = tetra_physics(ind_tetr)%sqg1 + sum(tetra_physics(ind_tetr)%gsqg * (x-x1))
 !
            end select  
        end select
!
      end function metric_determinant
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine check_tetra_overlaps
!
        use tetra_grid_mod, only: ntetr, tetra_grid, set_neighbour_face, verts_sthetaphi
        use constants, only: pi
!
        implicit none
!
        integer :: i,j,neighbour_tetra,neighbour_face, counter_wrong_cases
!
        counter_wrong_cases = 0
!
      !$OMP PARALLEL &
      !$OMP& DO DEFAULT(NONE) &
      !$OMP& SHARED(ntetr,tetra_grid,tetra_physics) &
      !$OMP& FIRSTPRIVATE(counter_wrong_cases) &
      !$OMP& PRIVATE(i,j,neighbour_tetra,neighbour_face)
        do i = 1,ntetr
!
          do j = 1,4
            neighbour_tetra = tetra_grid(i)%neighbour_tetr(j)
!
            if(neighbour_tetra.eq.-1) cycle
            neighbour_face = tetra_grid(i)%neighbour_face(j)
!
            !If opposing normal vectors point in the same direction, tetrahedra overlap
            if ( sum(tetra_physics(i)%anorm(:,j) * tetra_physics(neighbour_tetra)%anorm(:,neighbour_face)).ge.0.d0) then
              !The respective face is not valid anymore
              call set_neighbour_face(i,j,-1) !External setter: tetra_grid(i)%neighbour_face(j) = -1
              !$omp critical
                counter_wrong_cases = counter_wrong_cases +1
              !$omp end critical
              print *, 'tetr i', i
            endif
          enddo
        enddo
        !$OMP END PARALLEL DO
!
        if(counter_wrong_cases.eq.0) then
          print *, 'Normal vectors are valid: &
          & All opposing normal vectors of neighbouring tetrahedra with shared face are pointing in a different direction. '
        else
          print *, 'Some normal vectors are invalid: &
          & Some opposing normal vectors of neighbouring tetrahedra with shared face are pointing in the "same" direction.'
          print *, 'Error counter: ', counter_wrong_cases
        endif
!        
      end subroutine check_tetra_overlaps
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  end module tetra_physics_mod

