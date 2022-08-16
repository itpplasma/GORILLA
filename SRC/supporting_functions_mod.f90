!
  module supporting_functions_mod
!
    implicit none
!
    public
!
  contains
!
    subroutine sym_flux_in_cyl(sym_flux_filename,cyl_filename,n_additional_col)
!
      use tetra_physics_mod, only: coord_system
      use tetra_grid_settings_mod, only: grid_kind
!
      implicit none
!
      character(len=*), intent(in) :: sym_flux_filename,cyl_filename
      integer, intent(in) :: n_additional_col
      integer :: n_rows,io_error,i
      logical :: file_exist
      double precision :: x1,x2,x3,add1,add2,add3
!
      !Variables for magdata in symflux coordinates
      integer :: inp_label
      double precision :: psi_pol,q,dq_ds,sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta
      double precision :: Z, dZ_ds, dZ_dtheta,B_s,B_theta,B_theta1
      double precision :: theta_vmec,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                          alam,dR_dt,dR_dp,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
      !This subroutine is designed for maximum three additional columns
      if( (n_additional_col.lt.0).or.(n_additional_col.gt.3) ) then
        print *, 'Error: The subroutine sym_flux_in_cyl is desgined for maximum three additional columns.'
        print *, 'Please reduce the number of additional columns or adapt the subroutine.'
        stop
      endif
!
      !Check, if orbits were calculated in symmetry flux coordinate system
      if(coord_system.ne.2) then
        print *, "Error in sym_flux_in_cyl: Wrong coordinate system!"
        stop
      endif
!
      !Compute number of rows of sym_flux file
      n_rows = 0
      open(unit = 29, file=sym_flux_filename, status='old',action = 'read')
      do
        read(29,*,iostat = io_error)
        if(io_error.gt.0) exit
        n_rows = n_rows + 1
      enddo
      close(unit=29)
!
      !Open and read symmetry flux coordinates
      open(unit = 30, file=sym_flux_filename, status='old',action = 'read')
!
      !Read symmetry flux coordinates and write cylindrical coordinates
      inquire(file=cyl_filename,exist = file_exist)
      if(file_exist) then
        open(unit=31,file=cyl_filename,status='old',iostat=io_error)
        if(io_error.eq.0) close(unit=21,status='delete')
      else
        !print *,'Error: Could not open and delete existing file'
        print *,'No existing file. Created new one to store data.'
      endif
      open(unit=31,file=cyl_filename,status='new',action='write', &
           & iostat=io_error)
!
      do i = 1,n_rows-2
        select case(n_additional_col)
          case(0)
            read(30,*) x1,x2,x3
          case(1)
            read(30,*) x1,x2,x3,add1
          case(2)
            read(30,*) x1,x2,x3,add1,add2
          case(3)
            read(30,*) x1,x2,x3,add1,add2,add3
        end select
!
        select case(grid_kind)
          case(2) !EFIT field-aligned grid
            inp_label = 1
            call magdata_in_symfluxcoord_ext(inp_label,x1,psi_pol,x2,q,dq_ds, &
                                             sqrtg,bmod,dbmod_dtheta,R,dR_ds,dR_dtheta,       &
                                             Z,dZ_ds,dZ_dtheta)
          case(3) !VMEC field-aligned grid     
!
            !Find VMEC theta for sym-flux-theta
            x2 = theta_sym_flux2theta_vmec(x1,x2,x3)
!            
            call splint_vmec_data(x1,x2,x3,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                  R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
        end select !grid_kind
!        
        select case(n_additional_col)
          case(0)
            write(31,*) R,x3,Z
          case(1)
            write(31,*) R,x3,Z,add1
          case(2)
            write(31,*) R,x3,Z,add1,add2
          case(3)
            write(31,*) R,x3,Z,add1,add2,add3
        end select
!
      enddo
!
      close(unit=30)
      close(unit=31)
!
    end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function theta_sym_flux2theta_vmec(s,theta_sym_flux,varphi) result(theta_vmec)
!
      implicit none
!
      double precision :: s,varphi,theta_sym_flux
      double precision :: theta_vmec
      double precision :: dl_dt,alam,deltheta
      integer :: iter
      double precision, parameter :: epserr = 1.d-14
!      
!
      ! Begin Newton iteration to find VMEC theta
!
      theta_vmec = theta_sym_flux
!
      do iter=1,100
!
        call splint_lambda(s,theta_vmec,varphi,alam,dl_dt)
!
        deltheta = (theta_sym_flux-theta_vmec-alam)/(1.d0+dl_dt)
        theta_vmec = theta_vmec + deltheta
        if(abs(deltheta).lt.epserr) exit
      enddo
!      
    end function theta_sym_flux2theta_vmec
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    function theta_vmec2theta_sym_flux(s,theta_vmec,varphi) result(theta_sym_flux) 
!
        use constants, only: pi
!
        implicit none
!        
        double precision :: s,theta_vmec,theta_sym_flux,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                            R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
        call splint_vmec_data(s,theta_vmec,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!                                  
        theta_sym_flux = modulo(theta_vmec + alam,2.d0*pi)
!
    end function theta_vmec2theta_sym_flux
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_eigval(N,A,WR,VR)
!
      implicit none
!
      double precision, dimension(3) :: eigval
!
      INTEGER          N
      INTEGER          LDA, LDVL, LDVR
      !PARAMETER        ( LDA = N, LDVL = N, LDVR = N )
      INTEGER          LWMAX
      PARAMETER        ( LWMAX = 1000 )
!
!     .. Local Scalars ..
      INTEGER          INFO, LWORK
!
!     .. Local Arrays ..
      DOUBLE PRECISION A( N, N ), VL( N, N ), VR( N, N ), &
                        & WR( N ), WI( N ), WORK( LWMAX )
    
      external DGEEV
!
      LDA = N
      LDVL = N
      LDVR = N
!
      LWORK = -1
      CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,&
                 & VR, LDVR, WORK, LWORK, INFO )
      LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!
!     Solve eigenproblem.
!
      CALL DGEEV( 'Vectors', 'Vectors', N, A, LDA, WR, WI, VL, LDVL,&
                & VR, LDVR, WORK, LWORK, INFO )
    
!       print *, 'eigenvalues'
!       print '(E12.5)', WR
!      
    end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    elemental function logical2integer(boole_in)
!
        implicit none
!
        logical,intent(in)  :: boole_in
        integer             :: logical2integer
!
        if(boole_in) then
            logical2integer = 1
        else
            logical2integer = 0
        endif
!    
    end function logical2integer
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine alphas_runov2sym_flux(filename_runov,filename_symflux)
!
        use constants, only: pi
        use tetra_grid_settings_mod, only: n_field_periods
!
        implicit none
!
        character(len=*), intent(in) :: filename_runov,filename_symflux
        integer :: n_rows,io_error,i, i_tmp
        logical :: file_exist
        double precision, dimension(:,:), allocatable :: sthetaphilambda
        double precision :: dummy, x_tmp, phi_tmp, z_tmp, alpha_pitch,u_vmec
!
        !Compute number of rows of Runov's file
        n_rows = 0
        open(unit = 29, file=filename_runov, status='old',action = 'read')
!
        do
            read(29, '(A)', iostat=io_error)
            if (io_error /= 0) exit
            n_rows = n_rows + 1
        end do
!
        print*, "Runov's file with starting positions and starting pitch parameter contains ", n_rows, "starting values."
        print *, 'number of field periods', n_field_periods
!
        allocate(sthetaphilambda(4,n_rows))
!
        rewind(29)
!
        do i = 1, n_rows
            read(29,*) dummy, i_tmp, x_tmp, phi_tmp, z_tmp, alpha_pitch
!
            sthetaphilambda(1,i) = 5.0d-1
            sthetaphilambda(3,i) = modulo(phi_tmp,2.d0*pi/dble(n_field_periods))
            u_vmec = atan2(z_tmp, x_tmp)
            u_vmec = modulo(u_vmec, 2.d0*pi)
            sthetaphilambda(2,i) = theta_vmec2theta_sym_flux(sthetaphilambda(1,i),u_vmec,sthetaphilambda(3,i))
            sthetaphilambda(4,i) = alpha_pitch
        end do
!
        close(29)
!
        open(unit = 30, file=filename_symflux, status='unknown',action = 'write')
        do i = 1, n_rows
            write(30,*) sthetaphilambda(:,i)
        end do
        close(30)
!
    end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  end module
