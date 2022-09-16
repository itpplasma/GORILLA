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
    !additional optional orbit quantities
    logical, public, protected :: boole_time_Hamiltonian
    logical, public, protected :: boole_gyrophase
    logical, public, protected, dimension(2) :: boole_array_optional_quantities
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
    !Namelist for Gorilla input
    NAMELIST /GORILLANML/ eps_Phi, coord_system, ispecies, ipusher, &
                        & boole_pusher_ode45, boole_dt_dtau, boole_newton_precalc, poly_order, i_precomp, boole_guess, &
                        & rel_err_ode45,boole_periodic_relocation,handover_processing_kind, boole_axi_noise_vector_pot, &
                        & boole_axi_noise_elec_pot, boole_non_axi_noise_vector_pot, axi_noise_eps_A, axi_noise_eps_Phi, &
                        & non_axi_noise_eps_A, boole_helical_pert, helical_pert_eps_Aphi, helical_pert_m_fourier, &
                        & helical_pert_n_fourier, boole_time_Hamiltonian, boole_gyrophase, &
                        & boole_adaptive_time_steps, desired_delta_energy, max_n_intermediate_steps
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
        end subroutine load_gorilla_inp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine load_boole_array_optional_quantities
            implicit none

            boole_array_optional_quantities(1) = boole_time_Hamiltonian
            boole_array_optional_quantities(2) = boole_gyrophase

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
