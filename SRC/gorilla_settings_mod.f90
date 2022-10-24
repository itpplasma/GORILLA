!
  module gorilla_settings_mod
!
    implicit none
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    !optional quantities filled during a tetrahedron pushing
    type optional_quantities_type
        sequence
        double precision :: t_hamiltonian     !real time of tetrahedron passing
        double precision :: gyrophase         !gyrophase of particle
        double precision :: vpar_int
        double precision :: vpar2_int
    end type optional_quantities_type
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    private
!
    !Electrostatic potential
    !Change in the electrostatic potential within the plasma volume in Gaussian units
    double precision,public,protected  :: eps_Phi
!
    !Coordinate system
    integer, public, protected :: coord_system
!    
    !particle species
    integer, public, protected :: ispecies
!
    !periodic coordinate particle re-location
    logical, public, protected :: boole_periodic_relocation
!
    !Gorilla pusher options
    integer, public, protected :: ipusher
!
    !Numerical pusher options
    logical, public, protected :: boole_pusher_ode45
    double precision, public, protected :: rel_err_ode45
!    
    logical, public, protected :: boole_dt_dtau
    logical, public, protected :: boole_newton_precalc
    integer, public, protected :: poly_order
    integer, public, protected :: i_precomp
    logical, public, protected :: boole_guess
!
    !Time tracing options
    integer, public, protected :: i_time_tracing_option
!
    !additional optional orbit quantities
    logical, public, protected :: boole_time_Hamiltonian
    logical, public, protected :: boole_gyrophase
    logical, public, protected :: boole_vpar_int
    logical, public, protected :: boole_vpar2_int
    logical, public, protected, dimension(4) :: boole_array_optional_quantities
!
    !Processing of particle handover to tetrahedron neighbour
    integer, public, protected :: handover_processing_kind
!
    !Manipulation of the axisymmetric electromagnetic field with noise
    logical, public, protected :: boole_axi_noise_vector_pot
    logical, public, protected :: boole_axi_noise_elec_pot
    logical, public, protected :: boole_non_axi_noise_vector_pot
    double precision, public, protected :: axi_noise_eps_A
    double precision, public, protected :: axi_noise_eps_Phi
    double precision, public, protected :: non_axi_noise_eps_A
!
    !Manipulation of the axisymmetric electromagnetic field (Tokamak) with helical harmonic perturbation
    logical, public  :: boole_helical_pert
    double precision, public, protected :: helical_pert_eps_Aphi
    integer, public, protected :: helical_pert_m_fourier
    integer, public, protected :: helical_pert_n_fourier
!
    !Dividing orbit integration into intermediate steps (adaptive) to reduce error made by finit polynomial
    logical, public, protected :: boole_adaptive_time_steps
    double precision, public, protected :: desired_delta_energy
    integer, public, protected :: max_n_intermediate_steps
!
    !Including additional terms in case of strong electric fields (cylindrical coordinates only)
    logical, public, protected :: boole_strong_electric_field
!
    logical, public, protected :: boole_save_electric
    character*64,public,protected :: filename_electric_field
    character*64,public,protected :: filename_electric_drift
!
    !Boolean for precalculation of rectangular grid to improve find_tetra (sensible for n_particles >> 1)
    logical, public, protected :: boole_grid_for_find_tetra
!
    !Namelist for Gorilla input
    NAMELIST /GORILLANML/ eps_Phi, coord_system, ispecies, ipusher, &
                        & boole_pusher_ode45, boole_dt_dtau, boole_newton_precalc, poly_order, i_precomp, boole_guess, &
                        & rel_err_ode45,boole_periodic_relocation,handover_processing_kind, boole_axi_noise_vector_pot, &
                        & boole_axi_noise_elec_pot, boole_non_axi_noise_vector_pot, axi_noise_eps_A, axi_noise_eps_Phi, &
                        & non_axi_noise_eps_A, boole_helical_pert, helical_pert_eps_Aphi, helical_pert_m_fourier, &
                        & helical_pert_n_fourier, boole_time_Hamiltonian, boole_gyrophase, boole_vpar_int, boole_vpar2_int, &
                        & boole_adaptive_time_steps, desired_delta_energy, max_n_intermediate_steps, boole_grid_for_find_tetra, &
                        & i_time_tracing_option, &
                        & boole_strong_electric_field, boole_save_electric, filename_electric_field, filename_electric_drift
!
    public :: load_gorilla_inp,set_eps_Phi, optional_quantities_type
!
    contains
!
        subroutine load_gorilla_inp()
!
            open(unit=10, file='gorilla.inp', status='unknown')
            read(10,nml=gorillanml)
            close(10)
!
            print *,'GORILLA: Loaded input data from gorilla.inp'
!
            call load_boole_array_optional_quantities
!
            !Dependencies of input parameters
            if( i_time_tracing_option.eq.2 ) then
                if (ipusher.ne.2) then
                    print *, 'ERROR: When Hamiltonian time tracing is activated, set ipusher to 2 in gorilla.inp.'
                    stop
                endif
            endif
!
            !Gyrophase requires Hamiltonian time evolution (At least, the current implementation requires that)
            if ( boole_gyrophase.and.(.not.boole_time_Hamiltonian)) then
                print *, 'ERROR: When gyro-phase computation is activated, set boole_time_Hamiltonian to TRUE in gorilla.inp.'
                stop
            endif
!
            !The additional terms for a strong electric field are only implemented in case of cylindrical coordinates (WEST)
            !Also there is currently not an implementation using the extended precomputation mode of either pusher mode
            if ( boole_strong_electric_field.and.((coord_system.ne.1).or.(i_precomp.ne.0).or.(boole_newton_precalc)) ) then
                print *, 'ERROR: When strong electric fields are included, set coord_system=1 and ipusher=2 in gorilla.inp.'
                print *, '       Also have i_precomp=0 and boole_newton_precalc=.false. in gorilla.inp!'
                stop
            endif
!
            !The explicite electric field and drift due to it can only be given in case of using the strong electric field mode
            if (boole_save_electric.and.(.not.boole_strong_electric_field)) then
                print *, 'ERROR: boole_strong_electric_field has to be set true to save electric field (boole_save_electric=true).'
                print *, '       Switch boole_save_electric to false, if not using the strong electric field mode!'
                stop
            endif
!
        end subroutine load_gorilla_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine load_boole_array_optional_quantities
            implicit none

            boole_array_optional_quantities(1) = boole_time_Hamiltonian
            boole_array_optional_quantities(2) = boole_gyrophase
            boole_array_optional_quantities(3) = boole_vpar_int
            boole_array_optional_quantities(4) = boole_vpar2_int

        end subroutine load_boole_array_optional_quantities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine set_eps_Phi(eps_Phi_in)
!
            implicit none
!
            double precision, intent(in) :: eps_Phi_in
!
            eps_Phi = eps_Phi_in
!
        end subroutine set_eps_Phi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  end module gorilla_settings_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
