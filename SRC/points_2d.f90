module points_2d
implicit none
private
double precision, parameter :: pi = 3.14159265358979d0
public :: create_points_2d, scaling_func, create_points_2d_vmec
!double precision, public, protected :: s_min =  1.d-1
double precision, dimension(:),allocatable,public,protected :: r_frac_mod !in each element of n_theta is written how many theta fragments for that r-cut exist
integer, dimension(:),allocatable,public,protected :: n_theta_vec
interface
    function scaling_func(x) result(x_scaled)
        double precision, dimension(:), intent(in) :: x
        double precision, dimension(size(x)) :: x_scaled
    end function
end interface

contains

subroutine create_points_2d(inp_label,n_theta,points,points_s_theta_phi,r_scaling_func,theta_scaling_func,repeat_center_point)
!
    use tetra_grid_settings_mod, only: s_min => sfc_s_min,theta_geom_flux, n_extra_rings
    use magdata_in_symfluxcoor_mod, only: psipol_max
    use magdata_in_symfluxcoordinates_mod, only: magdata_in_symfluxcoord_ext
!
    integer, intent(in) :: inp_label
    integer, dimension(:), intent(in) :: n_theta
    double precision, dimension(:,:), intent(out) :: points
    double precision, dimension(:,:), intent(out) :: points_s_theta_phi
    procedure(scaling_func), optional :: r_scaling_func, theta_scaling_func
    logical, optional, intent(in) :: repeat_center_point ! set true if we are insym flux coords
    double precision, dimension(size(n_theta)) :: r_frac
    integer :: nlabel, i, j, isurf, point_idx, n_theta_current
!
    double precision, dimension(:), allocatable :: theta_frac,theta_flux
!
    double precision :: s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta
!
    integer :: n_center_point

    integer :: n
    double precision :: s_second_ring
!
    n_center_point = 1
    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_center_point = n_theta(1)
        end if
    end if
!
    nlabel = size(n_theta) ! nlabel = n_s = grid_size(1)
    allocate(theta_frac(maxval(n_theta)))
    allocate(theta_flux(maxval(n_theta)))
!
!    ! calculates points of a fine grid of the core region by integrating along field lines and writes the output to some files
!    ! TODO: dont call if files are up to date
!    call preload_for_SYNCH()
!!
!    ! loads points that are calculated in preload_for_SYNCH into globals from magdata_in_symfluxcoor_mod and spline interpolats them
!    ! this is done so we can then call magdata_in_symfluxcoord_ext on this interpolated fine grid.
!    call load_magdata_in_symfluxcoord()

! start: make first few entries of r_frac non-equidistant in order to avoid very thin triangles close to the magnetic axis
    n = n_extra_rings

    if (n.gt.0) then
        s_second_ring = s_min + (1.d0-s_min)/(dble(nlabel-n+1))
        do i = 1,n
            r_frac(i) = exp(log(s_min) + dble(i)*(log(s_second_ring)-log(s_min))/dble(n+1))
        enddo
    endif
! end: make first few entries of r_frac non-equidistant in order to avoid very thin triangles close to the magnetic axis
!      when using this version, be sure to comment out the line sstarting with "r_frac =" instead of "r_frac(n:nlabel) ="
!
    r_frac(n+1:nlabel) = s_min + [(dble(i)*(1.d0-s_min), i=1, nlabel-n, 1)] / (dble(nlabel-n))

    ! r_frac = s_min + [(dble(i)*(1.d0-s_min), i=1, nlabel, 1)] / (dble(nlabel))
!
    !r_frac = [(i, i = 1, nlabel,1)] / dfloat(nlabel)
    if (present(r_scaling_func)) r_frac = r_scaling_func(r_frac)
!
    if(.not. allocated(r_frac_mod)) then
        allocate(r_frac_mod(nlabel))
        r_frac_mod = r_frac
    endif
    if(.not. allocated(n_theta_vec)) then
        allocate(n_theta_vec(nlabel))
        n_theta_vec = n_theta
    endif
! 
    point_idx = 1
!     
    do isurf=0, nlabel !for inner stuff
        if (isurf == 0) then
!            n_theta_current = n_center_point
!            s = s_min !0.0
!            if (s .gt. r_frac(1)) then
!                print *,"Error in points_2d.f90: inner-most value for s must be smaller than r_frac(1)!"
!                stop
!            endif
            n_theta_current = n_center_point
            select case (inp_label)
            case (1)  ! normalized toroidal flux s
               s = s_min !0.0
               if (s .gt. r_frac(1)) then
                  print *,"Error in points_2d.f90: inner-most value for s (s_min) must be smaller than r_frac(1)!"
                  error stop
               endif
            case (2)  ! poloidal flux psi
               psi = s_min * psipol_max
               if (abs(psi) .gt. abs(r_frac(1))) then
                  print *,"Error in points_2d.f90: inner-most value for psi (s_min * psipol_max) must be smaller than r_frac(1)!"
                  error stop
               endif
            end select
        else
!            n_theta_current = n_theta(isurf)
!            s = r_frac(isurf)
            n_theta_current = n_theta(isurf)
            select case (inp_label)
            case (1)
               s = r_frac(isurf)
            case (2)
               psi = r_frac(isurf)
            end select
        end if
!
        theta_frac(:n_theta_current) = [(i, i=0, n_theta_current-1, 1)] / dfloat(n_theta_current)
        if (present(theta_scaling_func)) theta_frac(:n_theta_current) = &
                         theta_scaling_func(theta_frac(:n_theta_current))        
!
        select case(theta_geom_flux)
          case(1) !Theta scaling in flux coordinates
!            theta_flux = theta_frac * 2.d0 * pi
            theta_flux(:n_theta_current) = theta_frac(:n_theta_current) * 2.d0 * pi
!
          case(2) !Theta scaling in geometrical theta
            !Transform geometric theta to symmetry flux theta (!De-normalize theta_frac --> real angle)
            !call theta_geom2theta_flux(s,theta_frac*2.d0*pi,theta_flux)
            call theta_geom2theta_flux(inp_label, s, psi, &
               theta_frac(:n_theta_current) * 2.d0 * pi, theta_flux(:n_theta_current))
        end select   
!              
        do j = 1, n_theta_current
            theta = theta_flux(j)
!            call magdata_in_symfluxcoord_ext(1, s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta, &
!                                             R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta)
            call magdata_in_symfluxcoord_ext(inp_label, s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta, &
                                              R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta)
            points(:, point_idx) = [R, 0.d0, Z]
            points_s_theta_phi(:, point_idx) = [s, theta, 0.d0]
!
            point_idx = point_idx + 1
        end do
    end do
!
    deallocate(theta_frac)
    deallocate(theta_flux)
!
end subroutine create_points_2d
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine create_points_2d_vmec(n_theta, points_sthetaphi, r_scaling_func, theta_scaling_func, repeat_center_point)
!
    use tetra_grid_settings_mod, only: s_min => sfc_s_min
!
    integer, dimension(:), intent(in) :: n_theta
    !double precision, dimension(:,:), intent(out) :: points_rphiz
    double precision, dimension(:,:), intent(out) :: points_sthetaphi
    procedure(scaling_func), optional :: r_scaling_func, theta_scaling_func
    logical, optional, intent(in) :: repeat_center_point ! set true if we are insym flux coords
    double precision, dimension(size(n_theta)) :: r_frac
    integer :: nlabel, i, j, isurf, point_idx,n_theta_current
    integer :: n_center_point
    
    double precision, dimension(:), allocatable :: theta_frac

    double precision :: s,theta,vartheta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                        R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
    n_center_point = 1
    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_center_point = n_theta(1)
        end if
    end if
!
    nlabel = size(n_theta)
    allocate(theta_frac(maxval(n_theta)))
!
    r_frac = s_min + [(dble(i)*(1.d0-s_min), i=1, size(r_frac), 1)] / (dble(nlabel))
    if (present(r_scaling_func)) r_frac = r_scaling_func(r_frac)
!
    if(.not. allocated(r_frac_mod)) then
        allocate(r_frac_mod(nlabel))
        r_frac_mod = r_frac
    endif
    if(.not. allocated(n_theta_vec)) then
        allocate(n_theta_vec(nlabel))
        n_theta_vec = n_theta
    endif
!
    point_idx = 1
    do isurf=0, nlabel
!
        varphi = 0.d0
        if (isurf == 0) then
            n_theta_current = n_center_point
            s = s_min
            if (s .gt. r_frac(1)) then
                print *,"Error in points_2d.f90: inner-most value for s must be smaller than r_frac(1)!"
                error stop
            endif
        else
            n_theta_current = n_theta(isurf)
            s = r_frac(isurf)
        end if
!
        theta_frac(:n_theta_current) = [(i, i=0, n_theta_current-1, 1)] / dfloat(n_theta_current)
        if (present(theta_scaling_func)) theta_frac(:n_theta_current) = &
                         theta_scaling_func(theta_frac(:n_theta_current))
!        
        do j = 1, n_theta_current
            theta = theta_frac(j) * 2.d0 * pi
            points_sthetaphi(:, point_idx) = [s, theta, varphi]
            point_idx = point_idx + 1
        end do
    end do
!
    deallocate(theta_frac)
!
  end subroutine create_points_2d_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  !subroutine theta_geom2theta_flux(s,theta_geom_vec,theta_flux_vec)
  subroutine theta_geom2theta_flux(inp_label, s, psi, theta_geom_vec,theta_flux_vec)
!
    use magdata_in_symfluxcoor_mod, only : raxis,zaxis
    use field_line_integration_mod, only: theta0
    use tetra_grid_settings_mod, only: theta0_at_xpoint
    use binsrc_mod, only: binsrc
    use plag_coeff_mod, only: plag_coeff
    use magdata_in_symfluxcoordinates_mod, only: magdata_in_symfluxcoord_ext
!
    implicit none
!
    double precision :: s, diff, diff1, diff2
    double precision, dimension(:), intent(in) :: theta_geom_vec
    double precision, dimension(:), intent(out) :: theta_flux_vec
!    integer :: ntheta, i, i_theta_geom
    integer :: ntheta, i, i_theta_geom, inp_label
    integer, parameter :: ntheta_interp = 500, nplag = 10
    double precision, dimension(nplag) :: coef
    double precision :: psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta,R0,Z0
    double precision, dimension(ntheta_interp+nplag) :: theta_geom_interp, theta_flux_interp
!
    ntheta = size(theta_geom_vec)
    theta_flux_interp = dfloat([(i, i=-nplag/2,ntheta_interp-1+nplag/2, 1)])/ntheta_interp* 2.d0*pi 
    
!     print *, "theta_interp:",theta_flux_interp
!     
    !Calculate magnetic axis
!    call magdata_in_symfluxcoord_ext(1,0.d0,psi,0.d0,q,dq_ds, &
!                                           & sqrtg,bmod,dbmod_dtheta,R0,dR_ds,dR_dtheta, &
!                                           & Z0,dZ_ds,dZ_dtheta)
    R0=raxis  !<=sergei10.11.19
    Z0=zaxis  !<=sergei10.11.19
!
    !Calculate sampling points for interpolation
    do i = 1,ntheta_interp+nplag !size theta_flux_interp
      theta = modulo(theta_flux_interp(i),2.d0*pi)
!       print *, 'i,s,theta',i,s,theta
!      call magdata_in_symfluxcoord_ext(1, s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta, &
!                                             R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta)
      call magdata_in_symfluxcoord_ext(inp_label, s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta, &
                                              R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta)
!      theta_geom_interp(i) = modulo(atan2((Z-Z0),(R-R0)),2.d0*pi)-2.d0*pi !to shift it in the same regime
      if (theta0_at_xpoint) then
          !to shift it to the same regime
          theta_geom_interp(i) = modulo(atan2(Z - Z0, R - R0) - &
               theta0, 2.d0 * pi) - 2.d0 * pi
      else
          !to shift it to the same regime
          theta_geom_interp(i) = modulo(atan2(Z - Z0, R - R0), 2.d0 * pi) - 2.d0 * pi
      endif
!       print *, 'i, theta, theta_geom_interp', i, theta,theta_geom_interp(i)
      if(i.gt.1) then
        if(theta_geom_interp(i).lt.theta_geom_interp(i-1)) then
          theta_geom_interp(i) = theta_geom_interp(i)+ 2.d0*pi*(ceiling((theta_geom_interp(i-1)-(theta_geom_interp(i)))/(2.d0*pi))) !in case it intersects with 2pi-periodicity
!         print *,'Error in theta_geom2theta_flux: Not monotonously ascending transformation function.'
        endif  
! 
        if(theta_geom_interp(i)-theta_geom_interp(i-1) .gt. pi) then
          print *, 'Error in theta_geom2theta_flux -> big jump in theta_geom_interp curve'
          !This error was likely caused by adding 2pi -> this means that the curve was not monotonously increasing

!         print *,'Normalized theta_geom_interp values =',theta_geom_interp/(2.d0*pi)!(size(theta_geom_interp)-nplag/2)
!         print *,'Normalized theta_flux_interp values =',theta_flux_interp/(2.d0*pi)!(size(theta_flux_interp)-nplag/2)
          error stop
        endif  
!          
      endif
    enddo
!
    !Interpolate with Lagrange polynomials
    do i = 1,ntheta
!      if (theta_geom_vec(i) .eq. 0.d0) then
!        theta_flux_vec(i) = 0.d0
!        cycle
!        endif
!         
      !Find indices of theta_geom_vec(i) among sampling points
      !call binsrc(theta_geom_interp(nplag/2+1:ntheta_interp+nplag/2),1,ntheta_interp,theta_geom_vec(i),i_theta_geom)
      theta = theta_geom_vec(i)                                       !<=sergei10.11.19
      if(theta.lt.theta_geom_interp(nplag/2+1)) then                  !<=sergei10.11.19
         theta=theta+2.d0*pi                                           !<=sergei10.11.19
      elseif(theta.gt.theta_geom_interp(ntheta_interp+nplag/2)) then  !<=sergei10.11.19
         theta=theta-2.d0*pi                                           !<=sergei10.11.19
      endif                                                           !<=sergei10.11.19
      call binsrc(theta_geom_interp(nplag/2+1:ntheta_interp+nplag/2),1,ntheta_interp,theta,i_theta_geom)  !<=sergei10.11.19
      i_theta_geom = i_theta_geom + nplag/2
!       print *, 'binsrc, i_theta_geom', i_theta_geom
!       print *, 'theta_geom_vec(i)',theta_geom_vec(i)
!       print *, 'theta_geom_interp(i_theta_geom-1:i_theta_geom)',theta_geom_interp(i_theta_geom-1:i_theta_geom)
!      call plag_coeff(nplag,0,theta_geom_vec(i),theta_geom_interp(i_theta_geom-nplag/2:i_theta_geom+nplag/2-1),coef)  !<=sergei10.11.19
       call plag_coeff(nplag,0,theta,theta_geom_interp(i_theta_geom-nplag/2:i_theta_geom+nplag/2-1),coef)  !<=sergei10.11.19
      theta_flux_vec(i) = modulo(sum(coef*theta_flux_interp(i_theta_geom-nplag/2:i_theta_geom+nplag/2-1)),2.d0*pi)   
!       print *, 'geom values'
!       print *, theta_geom_interp(i_theta_geom-nplag/2:i_theta_geom-nplag/2+1), &
!      & theta_geom_vec(i),theta_geom_interp(i_theta_geom:i_theta_geom+1)
!       print *, 'flux values'
!       print *, theta_flux_interp(i_theta_geom-nplag/2:i_theta_geom-nplag/2+1), &
!      & theta_flux_vec(i),theta_flux_interp(i_theta_geom:i_theta_geom+1)
    enddo
!     
    do i = 1,ntheta
      if (theta_flux_vec(i).eq.0.d0) cycle
      theta = theta_flux_vec(i)
      call magdata_in_symfluxcoord_ext(1, s,psi,theta,q,dq_ds,sqrtg,bmod,dbmod_dtheta, &
                                             R,dR_ds,dR_dtheta,Z,dZ_ds,dZ_dtheta)
!
      diff1 = abs(modulo(atan2((Z-Z0),(R-R0)),2.d0*pi) - theta_geom_vec(i))
      diff2 = min(abs(modulo(atan2((Z-Z0),(R-R0)),2.d0*pi)), theta_geom_vec(i))+ 2.d0*pi - & !see if the distance over periodic
      &max(abs(modulo(atan2((Z-Z0),(R-R0)),2.d0*pi)), theta_geom_vec(i))                     !boundary is shorter than directly
!              
      diff = min(diff1,diff2)  
! 
    enddo
!    
  end subroutine theta_geom2theta_flux
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module points_2d

    
