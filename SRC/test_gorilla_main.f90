!
  program test_gorilla_main
!
!-------------------------------------------------------------------------------------------!
!
! Program:  GORILLA - Guiding-center ORbit Integration with Local Linearization Approach
!
! Authors:  Michael Eder^(1)
!           Christopher G. Albert^(1,2)
!           Lukas M. P. Bauer^(1)
!           Sergei V. Kasilov^(1,3)
!           Winfried Kernbichler^(1)
!           Markus Meisterhofer^(1)
!
! Affilliations: (1) Fusion@OEAW, Institut für Theoretische Physik - Computational Physics,
!                    Technische Universität Graz, Petersgasse 16, 8010 Graz, Austria
!
!                (2) Max-Planck-Institut für Plasmaphysik
!                    Boltzmannstr. 2, 85748 Garching, Germany
!
!                (3) Institute of Plasma Physics, National Science Center, “Kharkov Institute of Physics and Technology”,
!                    Akademicheskaya str. 1, 61108 Kharkov, Ukraine
!
! Copyright/License: ...
!
! Theoretical background for this computer program:
!   M. Eder, et al., "Quasi-geometric integration of guiding-center orbits in piecewise linear toroidal fields"
!   Physics of Plasmas 27,  122508 (2020)
!   https://doi.org/10.1063/5.0022117
!
!-------------------------------------------------------------------------------------------!
!
    use tetra_grid_settings_mod, only:      load_tetra_grid_inp
    use gorilla_settings_mod, only:         load_gorilla_inp
    use orbit_timestep_gorilla_mod, only:   initialize_gorilla,orbit_timestep_gorilla
    use gorilla_plot_mod, only:             gorilla_plot
!
!use supporting_functions_mod, only: alphas_runov2sym_flux
!
    implicit none
!
    !Settings for demonstration of GORILLA
    integer                           :: i_option           ! Demonstration options (hard-coded)
!
    !Variables for single orbit time step
    integer                           :: ind_tetr           ! index of tetrahedron
    integer                           :: iface              ! index of tetrahedron face
    logical                           :: boole_initialized  ! switch for initialization of particle in tetrahedron
    double precision                  :: vpar               ! parallel velocity of particle (in cm/s)
    double precision                  :: vperp              ! perpendicular velocity of particle (in cm/s)
    double precision                  :: t_step             ! length of orbit time step (in s)
    double precision, dimension(3)    :: x                  ! particle positions x_1, x_2, x_3 depending on coordinate system
!
!-------------------------------------------------------------------------------------------!
!------------------------- Load input files and initialization -----------------------------!
!
    print *, ''
    print *, 'GORILLA - Guiding-center ORbit Integration with Local Linearization Approach'
    print *, ''
!
    !Load tetrahedronal grid input data
    call load_tetra_grid_inp()
!
    !Load GORILLA settings input data
    call load_gorilla_inp()
!
    !Initialization:
    ! - Load and spline electromagnetic field (MHD equilbira or analytical field)
    ! - Generation of 3D Tetrahedronal grid
    ! - Piecewise linearization of electromagnetic field in tetrahedronal grid
    ! - Particle initialization (species, energy)
    call initialize_gorilla()
!
!-------------------------------------------------------------------------------------------!
!---------------------------- Demonstration of GORILLA -------------------------------------!
!
    !Hard-coded options:
    !1 ... GORILLA plot: Demonstration of various GORILLA settings (see gorilla_plot.inp)
    !2 ... GORILLA orbit time step: For implementation of GORILLA in applications
    i_option = 1
!
    select case(i_option)
!
!-------------------------------------------------------------------------------------------!
!
        case(1) !GORILLA plot
!
            print *, ''
            print *, 'Demonstration of GORILLA: Plotting'
            print *, ''
!
            !Perform orbit integration, plot Poincaré cuts/sections and evolution of invariants of motion
            ! - All plotting settings (incl. starting conditions) are selected in gorilla_plot.inp
            call gorilla_plot()
!
!-------------------------------------------------------------------------------------------!
!
        case(2) !GORILLA orbit time step
!
            print *, ''
            print *, 'Demonstration of GORILLA: Single orbit time step'
            print *, ''
!
            ! Particle must be initialized in grid at first call of orbit_timestep_gorilla
            !(NOT needed for further time steps of the same particle - see below)
            boole_initialized = .false.
!
            ! Particle position
            x = [0.5d0, 0.1d0, 0.63d0]
!
            ! Parallel velocity of particle (in cm/s)
            vpar = 37525024.533239894d0
!
            ! Perpendicular velocity of particle (in cm/s)
            vperp = 38283182.426206760d0
!
            ! Time step (in s) for orbit integration with GORILLA
            t_step = 0.1d0
!
            !First orbit timestep with GORILLA
            !Particle will be initialized in grid:
            ! - ind_tetr and iface are NOT set
            ! - ind_tetr and iface will be set (boole_initialized = FALSE)
            call orbit_timestep_gorilla(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface)
!
            print *, ''
            print *, 'AFTER 1. orbit time step with GORILLA:'
            print *, ''
            print *, 'Particle position',x
            print *, 'Parallel velocity of particle (in cm/s)',vpar
            print *, 'Perpendicular velocity of particle (in cm/s)',vperp
            print *, 'Particle initialized in tetrahedral grid',boole_initialized
            print *, 'Index of tetrahedron', ind_tetr
            print *, 'Index of tetrahedron face',iface
!
            !Second orbit timestep with GORILLA
            !Particle is already initialized in grid:
            ! - ind_tetr and iface are set
            ! - ind_tetr and iface will NOT be set (boole_initialized = TRUE)
            call orbit_timestep_gorilla(x,vpar,vperp,t_step,boole_initialized,ind_tetr,iface)
!
            print *, ''
            print *, 'AFTER 2. orbit time step with GORILLA:'
            print *, ''
            print *, 'Particle position',x
            print *, 'Parallel velocity of particle (in cm/s)',vpar
            print *, 'Perpendicular velocity of particle (in cm/s)',vperp
            print *, 'Particle initialized in tetrahedral grid',boole_initialized
            print *, 'Index of tetrahedron', ind_tetr
            print *, 'Index of tetrahedron face',iface
!
    end select
!-------------------------------------------------------------------------------------------!
!
! call alphas_runov2sym_flux('start_t.confined','orbit_start_sthetaphilambda.dat')
!
  end program test_gorilla_main
!
