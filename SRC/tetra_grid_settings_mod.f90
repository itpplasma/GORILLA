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
    !(1 ... rectangular, 2 ... field-aligned EFIT, 3 ... field-aligned VMEC,
    !  4 ... SOLEDGE3X-EIRENE, 5 ... analytic circular tokamak,
    !  6 ... field-aligned Boozer chartmap)
    integer,public,protected :: grid_kind
!
    !For VMEC input, retain only harmonics with toroidal mode number n=0.
    logical,public,protected :: boole_axisymmetric_only = .false.
!
    !MHD equilibrium filename
    character*64,public,protected :: g_file_filename
    character*64,public,protected :: convex_wall_filename
    character*100,public,protected :: netcdf_filename
!
    !MESH_SOLEDGE3X_EIRENE filename
    character(256),public,protected :: knots_SOLEDGE3X_EIRENE_filename
    character(256),public,protected :: triangles_SOLEDGE3X_EIRENE_filename
!
    !Symmetry Flux Coordinates Annulus
    double precision,public,protected :: sfc_s_min
!
    !Outer s-bound of the radial mesh (default 1.0d0 = full plasma).  Set
    !sfc_s_min < sfc_s_max <= 1.0d0 to build only an annular shell in s.
    double precision,public,protected :: sfc_s_max = 1.0d0
    !Radial grid spacing option
    !1 ... equidistant in s (normalized toroidal flux) - default
    !2 ... equidistant in sqrt(s) (approximately equidistant in minor radius)
    integer, public, protected :: i_radial_spacing = 1
!
    !Scaling of $\theta$-variable
    integer, public, protected :: theta_geom_flux
!
    !Option for $\theta$-variable origin
    logical, public, protected :: theta0_at_xpoint
!
    !grid_kind = 2: cluster the poloidal grid points toward the X-point (needs theta0_at_xpoint).
    !The local point spacing (as a fraction of theta) is proportional to 1 - a*exp(-(t/w)^2), t = distance
    !to the X-point direction in units of a full turn, a = theta_xpoint_clustering in [0,1),
    !w = theta_xpoint_clustering_width in (0,0.25]. Spacing at the X-point: (1-a)/(1-a*w*sqrt(pi))
    !times the mean spacing; far from it: 1/(1-a*w*sqrt(pi)) times the mean spacing.
    double precision, public, protected :: theta_xpoint_clustering = 0.d0
    double precision, public, protected :: theta_xpoint_clustering_width = 0.1d0
!
    !Axisymmetric equilibrium type (reader variable)
    integer, public, protected :: iaxieq_in
!
    !grid_kind = 4 (external knots/triangles mesh): .false. -> WEST equilibrium format (iaxieq = 1, default),
    !.true. -> EFIT g-file (iaxieq = 0), e.g. an EFIT mesh extended beyond the separatrix to the limiter
    logical, public, protected :: boole_ext_mesh_efit = .false.
!
    !grid_kind = 2 extended to the wall (g-file limiter contour): the flux-aligned grid_kind 2 core plus
    !an edge mesh out to the wall, built once by PREPROCESSING/EFIT_EDGE_MESH/make_edge_mesh.sh into
    !edge_mesh_dir. The run checks that the stored settings (edge_mesh_dir/edge_mesh.nml) match the
    !grid_kind 2 settings of tetra_grid.inp, the g-file and edge_mesh_options, and then uses the mesh via
    !the external-mesh path (internally grid_kind = 4 with an EFIT g-file). Needs coord_system = 1.
    logical, public, protected :: boole_extend_to_wall = .false.
    character(256), public, protected :: edge_mesh_dir = 'edge_mesh'
    character(512), public, protected :: edge_mesh_options = ''
    !grid_kind requested in tetra_grid.inp (differs from grid_kind if boole_extend_to_wall)
    integer, public, protected :: grid_kind_input
!
    !Analytic field circular axisymetric geometry (reader variable)
    integer, public, protected :: ianalytic_circ = 0
    double precision, public, protected :: R0_analytic_circ = 0.d0
    double precision, public, protected :: a_analytic_circ  = 0.d0
    double precision, public, protected :: B0_analytic_circ = 0.d0
    double precision, public, protected :: q0_analytic_circ = 1.d0
    double precision, public, protected :: q1_analytic_circ = 0.d0
    integer, public, protected :: q_exp_analytic_circ = 2   ! q=q0+q1*(r/a)^p exponent (2 or 4)
    character(256), public, protected :: q_profile_file_analytic_circ = ''  ! tabulated q(r) (overrides analytic when set)
    ! Prescribed static helical perturbation (grid_kind = 5, coord_system = 1):
    ! A_phi += amp Re[psihat(rho) e^{i alpha}], Phi += amp Re[Phihat(rho) e^{i alpha}],
    ! alpha = m theta* + n phi (theta* straight-field-line angle). Table columns:
    ! rho [cm], Re psihat, Im psihat [G cm^2], Re Phihat, Im Phihat [statV]. '' = off.
    character(256), public, protected :: helical_pert_file_analytic_circ = ''
    double precision, public, protected :: helical_pert_amp_analytic_circ = 1.d0
    integer, public, protected :: helical_pert_m_analytic_circ = 0
    integer, public, protected :: helical_pert_n_analytic_circ = 0
!
    !Object file with mesh data
    logical,public,protected :: boole_write_mesh_obj
    character(50),public,protected :: filename_mesh_rphiz,filename_mesh_sthetaphi
!
    !Namelist for Tetrahedronal Grid input
    NAMELIST /TETRA_GRID_NML/ n1, n2, n3, grid_kind,boole_n_field_periods,n_field_periods_manual, &
                            & sfc_s_min,sfc_s_max, &
                            & i_radial_spacing, &
                            & boole_axisymmetric_only, boole_ext_mesh_efit, &
                            & boole_extend_to_wall, edge_mesh_dir, edge_mesh_options, &
                            & boole_write_mesh_obj,filename_mesh_rphiz,filename_mesh_sthetaphi,theta_geom_flux,theta0_at_xpoint, &
                            & theta_xpoint_clustering, theta_xpoint_clustering_width, &
                            & g_file_filename,convex_wall_filename,netcdf_filename, &
                            & knots_SOLEDGE3X_EIRENE_filename, triangles_SOLEDGE3X_EIRENE_filename, &
                            & R0_analytic_circ, a_analytic_circ, B0_analytic_circ, q0_analytic_circ, q1_analytic_circ, &
                            & q_exp_analytic_circ, q_profile_file_analytic_circ, &
                            & helical_pert_file_analytic_circ, &
                            & helical_pert_amp_analytic_circ, &
                            & helical_pert_m_analytic_circ, helical_pert_n_analytic_circ
!
    public :: load_tetra_grid_inp,set_grid_size,set_n_field_periods,set_n2,set_n3
    public :: write_edge_mesh_core_record
!
    contains
!
        subroutine load_tetra_grid_inp()
!
            integer :: iunit
!
            open(newunit=iunit, file='tetra_grid.inp', status='unknown')
            read(iunit,nml=tetra_grid_nml)
            close(iunit)
!
            if (sfc_s_max .le. sfc_s_min .or. sfc_s_max .gt. 1.0d0) then
                print *, 'Error: sfc_s_max must satisfy sfc_s_min < sfc_s_max <= 1.0d0'
                print *, '  sfc_s_min =', sfc_s_min, '  sfc_s_max =', sfc_s_max
                error stop
            endif
!
            if (theta_xpoint_clustering .ne. 0.d0) then
                if (theta_xpoint_clustering .lt. 0.d0 .or. theta_xpoint_clustering .ge. 1.d0) then
                    print *, 'Error: theta_xpoint_clustering must satisfy 0 <= theta_xpoint_clustering < 1'
                    error stop
                endif
                if (theta_xpoint_clustering_width .le. 0.d0 .or. theta_xpoint_clustering_width .gt. 0.25d0) then
                    print *, 'Error: theta_xpoint_clustering_width must satisfy 0 < width <= 0.25'
                    error stop
                endif
                if (grid_kind .ne. 2) then
                    print *, 'Error: theta_xpoint_clustering is only implemented for grid_kind = 2'
                    error stop
                endif
                if (.not. theta0_at_xpoint) then
                    print *, 'Error: theta_xpoint_clustering needs theta0_at_xpoint = .true. ', &
                             '(theta = 0 must point to the X-point)'
                    error stop
                endif
            endif
!
            grid_kind_input = grid_kind
!
            if (boole_extend_to_wall) then
                if (grid_kind .ne. 2) then
                    print *, 'Error: boole_extend_to_wall = .true. needs grid_kind = 2'
                    error stop
                endif
                call check_edge_mesh_record()
                !Use the extended mesh through the external-mesh reader with the EFIT g-file
                grid_kind = 4
                boole_ext_mesh_efit = .true.
                knots_SOLEDGE3X_EIRENE_filename = trim(edge_mesh_dir)//'/knots_aug_edge.dat'
                triangles_SOLEDGE3X_EIRENE_filename = trim(edge_mesh_dir)//'/triangles_aug_edge.dat'
                print *, 'grid_kind = 2 extended to the wall: using the mesh in ', trim(edge_mesh_dir)
            endif
!
            !Set axisymmetric equilibrium type
            select case(grid_kind)
                case(4) !External mesh: WEST equilibrium or EFIT g-file
                    if (boole_ext_mesh_efit) then
                        iaxieq_in = 0
                    else
                        iaxieq_in = 1
                    endif
                case(5) !Analytic field circular axisymmetric geometry tokamak
                    ianalytic_circ = 1
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
                s_ratio = (sfc_s_min + (sfc_s_max-sfc_s_min)/(dble(grid_size_in(1)))) / sfc_s_min
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
        subroutine set_n2(n2_in)
!
            implicit none
!
            integer,intent(in)    :: n2_in
!
            n2 = n2_in
!
        end subroutine set_n2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine set_n3(n3_in)
!
            implicit none
!
            integer,intent(in)    :: n3_in
!
            n3 = n3_in
!
        end subroutine set_n3
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  !
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine gfile_signature(header, nbytes)
            !First line and size of the g-file: identifies the equilibrium the edge mesh was built for
            character(256), intent(out) :: header
            integer(8), intent(out) :: nbytes
            integer :: iu, ios
            logical :: exists
!
            header = ''
            nbytes = -1
            inquire(file=trim(g_file_filename), exist=exists, size=nbytes)
            if (.not. exists) return
            open(newunit=iu, file=trim(g_file_filename), status='old', action='read', iostat=ios)
            if (ios .ne. 0) return
            read(iu, '(a)', iostat=ios) header
            close(iu)
            header = quote_free(header)
        end subroutine gfile_signature
!
        function quote_free(str) result(out)
            !replace quotes (namelist string delimiters) by blanks
            character(*), intent(in) :: str
            character(len(str)) :: out
            integer :: i
            out = str
            do i = 1, len(out)
                if (out(i:i) .eq. "'" .or. out(i:i) .eq. '"') out(i:i) = ' '
            enddo
        end function quote_free
!
        function collapse_blanks(str) result(out)
            !trim and collapse runs of blanks to one blank (for comparing option strings)
            character(*), intent(in) :: str
            character(len(str)) :: out
            integer :: i, k
            logical :: prev_blank
            out = ''
            k = 0
            prev_blank = .true.
            do i = 1, len_trim(str)
                if (str(i:i) .eq. ' ') then
                    if (prev_blank) cycle
                    prev_blank = .true.
                else
                    prev_blank = .false.
                endif
                k = k + 1
                out(k:k) = str(i:i)
            enddo
            out = trim(out)
        end function collapse_blanks
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine write_edge_mesh_core_record(filename)
            !Written by the core dump of make_edge_mesh.sh: grid_kind 2 settings the core was built with.
            !The edge generator adds generator_options and stores it as edge_mesh_dir/edge_mesh.nml.
            character(*), intent(in) :: filename
            character(256) :: header
            integer(8) :: nbytes
            integer :: iu
!
            call gfile_signature(header, nbytes)
            open(newunit=iu, file=filename, status='replace', action='write')
            write(iu, '(a)') '&EDGE_MESH_NML'
            write(iu, '(a,i0,a)') ' n1 = ', n1, ' ,'
            write(iu, '(a,i0,a)') ' n3 = ', n3, ' ,'
            write(iu, '(a,es25.17,a)') ' sfc_s_min = ', sfc_s_min, ' ,'
            write(iu, '(a,es25.17,a)') ' sfc_s_max = ', sfc_s_max, ' ,'
            write(iu, '(a,i0,a)') ' i_radial_spacing = ', i_radial_spacing, ' ,'
            write(iu, '(a,i0,a)') ' theta_geom_flux = ', theta_geom_flux, ' ,'
            write(iu, '(a,l1,a)') ' theta0_at_xpoint = ', theta0_at_xpoint, ' ,'
            write(iu, '(a,es25.17,a)') ' theta_xpoint_clustering = ', theta_xpoint_clustering, ' ,'
            write(iu, '(a,es25.17,a)') ' theta_xpoint_clustering_width = ', theta_xpoint_clustering_width, ' ,'
            write(iu, '(a,a,a)') " gfile_header = '", trim(header), "' ,"
            write(iu, '(a,i0,a)') ' gfile_size = ', nbytes, ' ,'
            write(iu, '(a)') '/'
            close(iu)
        end subroutine write_edge_mesh_core_record
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_edge_mesh_record()
            !Compare edge_mesh_dir/edge_mesh.nml with the settings of this run; stop with the command
            !that builds a matching mesh if the files are missing or differ.
            integer :: r_n1, r_n3, r_i_radial_spacing, r_theta_geom_flux
            double precision :: r_sfc_s_min, r_sfc_s_max, r_theta_xpoint_clustering, r_theta_xpoint_clustering_width
            logical :: r_theta0_at_xpoint
            character(256) :: r_gfile_header, header
            integer(8) :: r_gfile_size, nbytes
            character(512) :: r_generator_options
            namelist /EDGE_MESH_NML/ r_n1, r_n3, r_sfc_s_min, r_sfc_s_max, r_i_radial_spacing, r_theta_geom_flux, &
                                     r_theta0_at_xpoint, r_theta_xpoint_clustering, r_theta_xpoint_clustering_width, &
                                     r_gfile_header, r_gfile_size, r_generator_options
            character(1024) :: line
            character(512) :: fname
            logical :: exists, ok
            integer :: iu, ios, nmis
!
            fname = trim(edge_mesh_dir)//'/edge_mesh.nml'
            inquire(file=trim(fname), exist=exists)
            ok = exists
            nmis = 0
            if (exists) then
                !the record uses the names of tetra_grid.inp; read it through a copy with r_ prefixes
                call read_record(trim(fname), ios)
                if (ios .ne. 0) then
                    print *, 'Error: cannot read ', trim(fname)
                    ok = .false.
                else
                    call gfile_signature(header, nbytes)
                    call cmp_i('n1', r_n1, n1)
                    call cmp_i('n3', r_n3, n3)
                    call cmp_d('sfc_s_min', r_sfc_s_min, sfc_s_min)
                    call cmp_d('sfc_s_max', r_sfc_s_max, sfc_s_max)
                    call cmp_i('i_radial_spacing', r_i_radial_spacing, i_radial_spacing)
                    call cmp_i('theta_geom_flux', r_theta_geom_flux, theta_geom_flux)
                    if (r_theta0_at_xpoint .neqv. theta0_at_xpoint) call mismatch('theta0_at_xpoint')
                    call cmp_d('theta_xpoint_clustering', r_theta_xpoint_clustering, theta_xpoint_clustering)
                    if (theta_xpoint_clustering .gt. 0.d0) call cmp_d('theta_xpoint_clustering_width', &
                        r_theta_xpoint_clustering_width, theta_xpoint_clustering_width)
                    if (trim(adjustl(r_gfile_header)) .ne. trim(adjustl(header)) .or. r_gfile_size .ne. nbytes) &
                        call mismatch('g-file (header line or size of '//trim(g_file_filename)//')')
                    if (collapse_blanks(r_generator_options) .ne. collapse_blanks(edge_mesh_options)) then
                        call mismatch('edge_mesh_options')
                        print *, '    mesh: ', trim(collapse_blanks(r_generator_options))
                        print *, '    run : ', trim(collapse_blanks(edge_mesh_options))
                    endif
                    ok = (nmis .eq. 0)
                endif
            else
                print *, 'Error: boole_extend_to_wall = .true., but there is no edge mesh in ', trim(edge_mesh_dir)
            endif
!
            if (ok) return
!
            if (theta_geom_flux .ne. 2) print *, '  Note: the extended mesh runs in cylindrical coordinates; ', &
                'use theta_geom_flux = 2 (theta_geom_flux = 1 folds near the separatrix).'
            print *, '  Build a matching mesh once (about 10 s), from this run directory:'
            write(line, '(a)') '    GFILE=$PWD/'//trim(g_file_filename)
            if (sfc_s_max .ne. 1.d0 .or. i_radial_spacing .ne. 1 .or. theta_xpoint_clustering .gt. 0.d0 &
                .or. theta_geom_flux .ne. 2 .or. .not. theta0_at_xpoint) then
                write(line, '(a)') trim(line)//' CORE_EXTRA="'
                if (sfc_s_max .ne. 1.d0) write(line, '(a,g0,a)') trim(line)//' sfc_s_max = ', sfc_s_max, ' ,'
                if (i_radial_spacing .ne. 1) write(line, '(a,i0,a)') trim(line)//' i_radial_spacing = ', &
                    i_radial_spacing, ' ,'
                if (theta_geom_flux .ne. 2) write(line, '(a,i0,a)') trim(line)//' theta_geom_flux = ', &
                    theta_geom_flux, ' ,'
                if (.not. theta0_at_xpoint) write(line, '(a)') trim(line)//' theta0_at_xpoint = .false. ,'
                if (theta_xpoint_clustering .gt. 0.d0) write(line, '(a,g0,a,g0,a)') &
                    trim(line)//' theta_xpoint_clustering = ', theta_xpoint_clustering, &
                    ' , theta_xpoint_clustering_width = ', theta_xpoint_clustering_width, ' ,'
                write(line, '(a)') trim(line)//'"'
            endif
            print *, trim(line)//' \'
            write(line, '(a,a,i0,a,i0,a,g0,a,a)') '    <GORILLA>/PREPROCESSING/EFIT_EDGE_MESH/make_edge_mesh.sh ', &
                trim(edge_mesh_dir)//' ', n1, ' ', n3, ' ', sfc_s_min, ' -- ', trim(collapse_blanks(edge_mesh_options))
            print *, trim(line)
            print *, '  (see PREPROCESSING/EFIT_EDGE_MESH/README.md; DUMP_X and PYTHON may need to be set)'
            error stop
!
            contains
!
            subroutine read_record(fn, ios_out)
                !copy the record to a scratch unit with r_-prefixed names and read the namelist from it
                character(*), intent(in) :: fn
                integer, intent(out) :: ios_out
                integer :: iin, iscr, ios2, k
                character(1024) :: ln
                open(newunit=iin, file=fn, status='old', action='read', iostat=ios_out)
                if (ios_out .ne. 0) return
                open(newunit=iscr, status='scratch', action='readwrite')
                do
                    read(iin, '(a)', iostat=ios2) ln
                    if (ios2 .ne. 0) exit
                    ln = adjustl(ln)
                    k = index(ln, '=')
                    if (ln(1:1) .ne. '&' .and. ln(1:1) .ne. '/' .and. k .gt. 1) ln = 'r_'//trim(ln)
                    write(iscr, '(a)') trim(ln)
                enddo
                close(iin)
                rewind(iscr)
                r_generator_options = ''
                r_gfile_header = ''
                r_theta_xpoint_clustering = 0.d0
                r_theta_xpoint_clustering_width = 0.1d0
                read(iscr, nml=EDGE_MESH_NML, iostat=ios_out)
                close(iscr)
            end subroutine read_record
!
            subroutine mismatch(what)
                character(*), intent(in) :: what
                if (nmis .eq. 0) print *, 'Error: the edge mesh in ', trim(edge_mesh_dir), &
                    ' was built with different settings:'
                nmis = nmis + 1
                print *, '  differs: ', what
            end subroutine mismatch
!
            subroutine cmp_i(what, a, b)
                character(*), intent(in) :: what
                integer, intent(in) :: a, b
                if (a .ne. b) then
                    call mismatch(what)
                    print *, '    mesh:', a, '  run:', b
                endif
            end subroutine cmp_i
!
            subroutine cmp_d(what, a, b)
                character(*), intent(in) :: what
                double precision, intent(in) :: a, b
                if (abs(a - b) .gt. 1.d-12 * max(abs(a), abs(b), 1.d-300)) then
                    call mismatch(what)
                    print *, '    mesh:', a, '  run:', b
                endif
            end subroutine cmp_d
!
        end subroutine check_edge_mesh_record
!
end module tetra_grid_settings_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
