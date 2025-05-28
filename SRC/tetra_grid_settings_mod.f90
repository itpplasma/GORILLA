!
  module tetra_grid_settings_mod
!
    implicit none
!
    private
!
    !Grid options
    logical,public,protected    :: boole_n_field_periods
    integer,public,protected    :: n_field_periods_manual
    integer,public,protected    :: n_field_periods
!    
    !Grid Size
    integer,public,protected  :: n1
    integer,public,protected  :: n2
    integer,public,protected  :: n3
!
    !In certain regions (e.g. close to the magnetic axis and / or close to rational surfaces), extra grid points might be added in 
    ! radial direction if grid_kind = 2 or 3
    integer,public,protected  :: n_extra_rings
!
    !Grid Size - for practicability - NO input quantity
    !Rectangular grid: [n1+n_extra_rings,n2,n3] = [nR,nphi,nZ], Field-aligned grid: [n1+n_extra_rings,n2,n3] = [nr,nphi,ntheta]
    integer, dimension(3),public,protected  :: grid_size

!
    !Grid kind - Selection of grid kind
    !(1 ... rectangular grid, 2 ... field-aligned grid EFIT, 3 ... field-aligned grid VMEC, 4 ... SOLEDGE3X-EIRENE grid)
    integer,public,protected :: grid_kind
!
    !MHD equilibrium filename
    character*64,public,protected :: g_file_filename
    character*64,public,protected :: convex_wall_filename
    character*100,public,protected :: netcdf_filename
!
    !MESH_SOLEDGE3X_EIRENE filename
    character*64,public,protected :: knots_SOLEDGE3X_EIRENE_filename
    character*64,public,protected :: triangles_SOLEDGE3X_EIRENE_filename
!
    !Symmetry Flux Coordinates Annulus
    double precision,public,protected :: sfc_s_min
!
    !Scaling of $\theta$-variable
    integer, public, protected :: theta_geom_flux
!
    !Option for $\theta$-variable origin
    logical, public, protected :: theta0_at_xpoint
!
    !Axisymmetric equilibrium type (reader variable)
    integer, public, protected :: iaxieq_in
!
    !Object file with mesh data
    logical,public,protected :: boole_write_mesh_obj
    character(50),public,protected :: filename_mesh_rphiz,filename_mesh_sthetaphi
!
    !Namelist for Tetrahedronal Grid input
    NAMELIST /TETRA_GRID_NML/ n1, n2, n3, grid_kind,boole_n_field_periods,n_field_periods_manual,sfc_s_min, &
                            & boole_write_mesh_obj,filename_mesh_rphiz,filename_mesh_sthetaphi,theta_geom_flux,theta0_at_xpoint, &
                            & g_file_filename,convex_wall_filename,netcdf_filename, &
                            & knots_SOLEDGE3X_EIRENE_filename, triangles_SOLEDGE3X_EIRENE_filename
!
    public :: load_tetra_grid_inp,set_grid_size,set_n_field_periods
!
    contains
!
        subroutine load_tetra_grid_inp()
!
            open(unit=9, file='tetra_grid.inp', status='unknown')
            read(9,nml=tetra_grid_nml)
            close(9)
!
            !Set axisymmetric equilibrium type
            select case(grid_kind)
                case(4) !Axisymmetric equilibrium of WEST
                    iaxieq_in = 1
                case default
                    iaxieq_in = 0
            end select
!
            print *,'Tetrahedronal Grid: Loaded input data from tetra_grid.inp'
!            
        end subroutine load_tetra_grid_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine set_n_field_periods(n_field_periods_in)
!
            implicit none
!
            integer,intent(in)    :: n_field_periods_in
!
            n_field_periods = n_field_periods_in
!
        end subroutine set_n_field_periods
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine set_grid_size(grid_size_in)
!
            implicit none
!
            integer, dimension(3),intent(in)  :: grid_size_in
            double precision :: s_ratio
!
            if (grid_kind.eq.2) then
                !s_ratio = ratio of innermost s value (= sfc_s_min) and s value at second innermost ring
                s_ratio = (sfc_s_min + (1.d0-sfc_s_min)/(dble(grid_size_in(1)))) / sfc_s_min
                n_extra_rings = int(abs(log(s_ratio)/log(10d0)))
            else
                n_extra_rings = 0
            endif

            grid_size = grid_size_in
            grid_size(1) = grid_size(1) + n_extra_rings
!
        end subroutine set_grid_size
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  end module tetra_grid_settings_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
