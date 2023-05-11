!
module poly_without_precomp_mod
!
    implicit none
!
    double precision, dimension(4)      :: b
    double precision, dimension(4,4)    :: amat,amat2,amat3,amat4
    double precision, dimension(4)      :: amat_in_z,amat2_in_z,amat3_in_z,amat_in_b,amat2_in_b
    double precision, dimension(4)      :: amat4_in_z,amat3_in_b
    double precision                    :: perpinv3,perpinv4
!    
    !$OMP THREADPRIVATE(b,amat,amat2,amat3,amat4,amat_in_z,amat2_in_z,amat3_in_z,amat4_in_z, & 
    !$OMP& amat_in_b,amat2_in_b,amat3_in_b,perpinv3,perpinv4)
!    
end module poly_without_precomp_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module pusher_tetra_poly_mod

    implicit none
!
    interface poly_multiplication
    pure function scalar_multiplication_without_precomp(scalar_coef_1,scalar_coef_2)
        implicit none
        double precision, dimension(:), intent(in)           :: scalar_coef_1,scalar_coef_2
        double precision, dimension(:), allocatable          :: scalar_multiplication_without_precomp
    end function scalar_multiplication_without_precomp
!
    pure function vector_multiplication_without_precomp(scalar_coef,vector_coef)
        implicit none
        double precision, dimension(:), intent(in)           :: scalar_coef
        double precision, dimension(:,:), intent(in)         :: vector_coef
        double precision, dimension(:,:), allocatable        :: vector_multiplication_without_precomp
        end function vector_multiplication_without_precomp
!
    pure function tensor_multiplication_without_precomp(vector_coef_1,vector_coef_2)
        implicit none
        double precision, dimension(:,:), intent(in)         :: vector_coef_1,vector_coef_2
        double precision, dimension(:,:,:), allocatable      :: tensor_multiplication_without_precomp
        end function tensor_multiplication_without_precomp
    end interface poly_multiplication
!
    interface moment_integration
    pure function scalar_integral_without_precomp(poly_order,tau,scalar_coef)
        implicit none
        integer, intent(in)                                  :: poly_order
        double precision, intent(in)                         :: tau
        double precision, dimension(:), intent(in)           :: scalar_coef
        double precision                                     :: scalar_integral_without_precomp
    end function scalar_integral_without_precomp
!
    pure function vector_integral_without_precomp(poly_order,tau,vector_coef)
        implicit none
        integer, intent(in)                                   :: poly_order
        double precision, intent(in)                          :: tau
        double precision, dimension(:,:), intent(in)          :: vector_coef
        double precision, dimension(:), allocatable           :: vector_integral_without_precomp
    end function vector_integral_without_precomp
!
    pure function tensor_integral_without_precomp(poly_order,tau,tensor_coef)
        implicit none
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(:,:,:), intent(in)  :: tensor_coef
        double precision, dimension(:,:), allocatable   :: tensor_integral_without_precomp
    end function tensor_integral_without_precomp
    end interface moment_integration
!
    !change those for adaptive step sizes, probably allocatable
    double precision, dimension(:), allocatable, public, protected            :: tau_steps_list
    double precision, dimension(:,:), allocatable, public, protected          :: intermediate_z0_list
    integer, public, protected                                                :: number_of_integration_steps
!
    private
!    
    integer                             :: iface_init
    integer, public, protected          :: ind_tetr
    integer, public, protected          :: sign_rhs
    double precision                    :: perpinv2,vmod0
    double precision, public, protected :: perpinv,dt_dtau_const,bmod0
    double precision                    :: t_remain
    double precision, dimension(3)      :: x_init
    double precision, dimension(4), public, protected :: z_init
    double precision                    :: k1, k3
    double precision,parameter          :: eps_tau = 100.d0
    double precision, dimension(4,4)    :: unity_matrix4 = reshape( [ 1.d0, 0.d0, 0.d0, 0.d0, &
                                0.d0, 1.d0, 0.d0, 0.d0, & 
                                0.d0, 0.d0, 1.d0, 0.d0, &
                                0.d0, 0.d0, 0.d0, 1.d0 ], [4,4])
!
    !$OMP THREADPRIVATE(ind_tetr,iface_init,perpinv,perpinv2,dt_dtau_const,bmod0,t_remain,x_init,  &
    !$OMP& z_init,k1,k3,vmod0,tau_steps_list,intermediate_z0_list,number_of_integration_steps,sign_rhs)
!
    public :: pusher_tetra_poly,initialize_const_motion_poly, &
        & Quadratic_Solver2, Cubic_Solver, Quartic_Solver,analytic_integration_without_precomp, &
        & poly_multiplication_coef, moment_integration,set_integration_coef_manually
!   
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_const_motion_poly(perpinv_in,perpinv2_in)
!
            implicit none
!            
            double precision, intent(in)    :: perpinv_in,perpinv2_in
!
            perpinv = perpinv_in
            perpinv2 = perpinv2_in
!            
        end subroutine initialize_const_motion_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine manage_intermediate_steps_arrays(option)
!
            use gorilla_settings_mod, only: max_n_intermediate_steps, boole_adaptive_time_steps
            use gorilla_diag_mod, only: report_pusher_tetry_poly_adaptive
!
            implicit none
!            
            integer, intent(in) :: option
!
            integer             :: max_entries, k
!
            if (option .eq. 0) then
                !From one face to another there are max_n_intermediate_steps 
                !and that all possible two times as trajectory can be prolonged to a third face (turning on face) + buffer
                if (boole_adaptive_time_steps) then
                    max_entries = 3*(max_n_intermediate_steps)
                    allocate(tau_steps_list(max_entries),intermediate_z0_list(4,max_entries))
                else
                    !if no adaptive scheme -> just two steps in total
                    allocate(tau_steps_list(2),intermediate_z0_list(4,2))
                endif
!
if(report_pusher_tetry_poly_adaptive) then
open(124, file='./total_fluc_report.dat')
open(118, file='./step_fluc_report.dat')
endif

            elseif(option .eq. 1) then
                deallocate(tau_steps_list,intermediate_z0_list)
!
if(report_pusher_tetry_poly_adaptive) then
close(124)
close(118)
endif
!
            endif
!            
        end subroutine manage_intermediate_steps_arrays
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
            use tetra_physics_mod, only: tetra_physics, sign_sqg
            use gorilla_settings_mod, only: boole_strong_electric_field
            use supporting_functions_mod, only: bmod_func, phi_elec_func, v2_E_mod_func
!
            implicit none
!
            integer, intent(in)                        :: ind_tetr_inout,iface
            double precision, intent(in)               :: vpar
            double precision, dimension(3), intent(in) :: x
            double precision                           :: vperp2,vpar2,t_remain_in,phi_elec
!
            t_remain = t_remain_in
!    
            ind_tetr=ind_tetr_inout           !Save the index of the tetrahedron locally
!
            !Sign of the right hand side of ODE - ensures that tau is ALWAYS positive inside the algorithm
            sign_rhs = sign_sqg * int(sign(1.d0,t_remain))
!
            z_init(1:3)=x-tetra_physics(ind_tetr)%x1       !x is the entry point of the particle into the tetrahedron in (R,phi,z)-coordinates
!
            z_init(4)=vpar                         !Transform to z_init: 1st vertex of tetrahdron is new origin
!
            !Save initial orbit parameters
            x_init = x
!             vperp_init = vperp
            iface_init = iface
!
            !Tetrahedron constants
            dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const
!
            !Multiply with sign of rhs - ensures that tau is ALWAYS positive inside the algorithm
            dt_dtau_const = dt_dtau_const*dble(sign_rhs)
!    
            !Module of B at the entry point of the particle
            bmod0 = bmod_func(z_init(1:3),ind_tetr) !tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z_init(1:3))
!
            !Phi at the entry point of the particle
            phi_elec = phi_elec_func(z_init(1:3),ind_tetr)   !tetra_physics(ind_tetr)%Phi1+sum(tetra_physics(ind_tetr)%gPhi*z_init(1:3))
!
            !Auxiliary quantities
            vperp2 = -2.d0*perpinv*bmod0
            vpar2 = vpar**2
            !This is the total speed viewed in the MOVING FRAME of ExB drift (here it only acts as a coefficient for the EOM-set)
            !For physical estimates v_E is considered seperately anyway
            vmod0 = sqrt(vpar2+vperp2)
!
            k1 = vperp2+vpar2+2.d0*perpinv*tetra_physics(ind_tetr)%bmod1
            !Adding in +dx*grad(v2emod) to k1 helper-coefficiant for strong electric field case
            if (boole_strong_electric_field) k1 = k1 + (v2_E_mod_func(z_init(1:3),ind_tetr) - tetra_physics(ind_tetr)%v2Emod_1)
            k3 = tetra_physics(ind_tetr)%Phi1-phi_elec
!
        end subroutine initialize_pusher_tetra_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine pusher_tetra_poly(poly_order,ind_tetr_inout,iface,x,vpar,z_save,t_remain_in,t_pass, &
                        & boole_t_finished,iper_phi,optional_quantities)
!
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            use gorilla_diag_mod,  only: diag_pusher_tetry_poly
            use pusher_tetra_func_mod, only: pusher_handover2neighbour
            use gorilla_settings_mod, only: i_precomp, boole_guess, optional_quantities_type, boole_array_optional_quantities, &
                                    &  boole_adaptive_time_steps, i_time_tracing_option
!
            implicit none
!
            integer, intent(in)                                   :: poly_order
            double precision, intent(in)                          :: t_remain_in
!
            integer,intent(inout)                                 :: ind_tetr_inout,iface
            double precision, dimension(3), intent(inout)         :: x
            double precision, intent(inout)                       :: vpar
!
            double precision, dimension(3), intent(out)           :: z_save
            double precision, intent(out)                         :: t_pass
            logical, intent(out)                                  :: boole_t_finished
            integer,intent(out)                                   :: iper_phi
            type(optional_quantities_type), intent(out),optional  :: optional_quantities
!
            logical, dimension(4)                                 :: boole_faces
            integer                                               :: i,j,k
            double precision, dimension(4)                        :: z,operator_b_in_b,z_dummy
            integer                                               :: iface_new,i_scaling, max_entries
            double precision                                      :: tau,vperp2,tau_save,tau_est,tau_max, energy_init,energy_current
            logical                                               :: boole_analytical_approx,boole_face_correct
            logical                                               :: boole_trouble_shooting
            double precision, dimension(4,4)                      :: operator_b,operator_z_init
            double precision, allocatable                         :: t_hamiltonian
            double precision, dimension(:), allocatable           :: t_hamiltonian_list
            integer                                               :: i_step_root
            double precision                                      :: t_remain_new
!         
            call initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
            !In case of first call of orbit integration -> initialize the intermediate_steps_array
            if (.not.(allocated(tau_steps_list).OR.allocated(intermediate_z0_list))) call manage_intermediate_steps_arrays(0)
!
            !initialise the module variables number_of_integration_steps, tau_steps_list and intermediate_z0_list
            ! max_entries = 2*(max_n_intermediate_steps)
            ! allocate(tau_steps_list(max_entries),intermediate_z0_list(4,max_entries))
            number_of_integration_steps = 0
!
            if (any(boole_array_optional_quantities)) call initialise_optional_quantities(optional_quantities)
!
            !Initial computation values
            z = z_init
! 
if(diag_pusher_tetry_poly) then
    print *, 'z_init', z_init
    print *, 'iface init', iface_init
    print *, 'norm at start'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i)
        print *,i, 'normal velocity',normal_velocity_func(z,i)
    enddo
!
    call physical_estimate_tau(tau_est)
    print *, 'tau estimate',tau_est
!
endif
!
            !Initialize iper_phi (for the case that handover2neighbour does not set this value)
            iper_phi = 0
!
            !Initialize scaling for quartic solver
            i_scaling = 0
!
            !Initialize boole_t_finished
            boole_t_finished = .false.
!
            !Initialize boole trouble shooting
            boole_trouble_shooting = .true.
!
            !Set iface_new start value for polynomial approximation
            iface_new = iface_init !instead of iface_new = iface
!
            !boole_faces ... Boolean array for switching on/off root computation for individual face
            !Initialize boole_faces, such that roots for all 4 faces are computed
            boole_faces = .true.
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!FIRST ATTEMPT WITH SECOND ORDER GUESS!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            !Analytical calculation of orbit parameter to guess exit face and estimate of tau
            call analytic_approx(2,i_precomp,boole_faces, &
            & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
            !Define boundary for tau_max
            tau_max = tau*eps_tau
!
if(diag_pusher_tetry_poly) print *, 'boole_analytical_approx',boole_analytical_approx
!
            !use face prediction of second order for higher order computation (unneccessary if poly_order = 2)
            if(boole_guess .and. boole_analytical_approx .and. (poly_order.gt.2)) then               
                !Higher order polynomial root is computed only for previously predicted face in second order
                boole_faces = .false.                !Disable all 4 faces 
                boole_faces(iface_new) = .true.      !Enable guessed face
!
if(diag_pusher_tetry_poly) print *, 'tau',tau
if(diag_pusher_tetry_poly) print *, 'iface guess', iface_new
!
            endif   
!           
            !calculate exit time and exit face in higher order
            !if a successful guess was made in second order, boole_faces only allows the guessed face 
            if(poly_order.gt.2) then
                iface_new = iface_init !instead of iface_new = iface 
!                
                !Analytical calculation of orbit parameter to pass tetrahdron
                call analytic_approx(poly_order,i_precomp,boole_faces, &
                                    & i_scaling,z,iface_new,tau,boole_analytical_approx)   
            endif
!
            !Initialize face error recognition procedure
            boole_face_correct = .true.
!           
            !Detect face error: Particle is not lying in the correct exit face and/or velocity points inwards the tetrahedron    
!
            !Analytical result does not exist.
            if(.not.boole_analytical_approx) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: Analytic approximation'
            endif
!
            !Integrate trajectory analytically, if root exists
            if(boole_face_correct) then
                call analytic_integration(poly_order,i_precomp,z,tau)
!
                if (boole_adaptive_time_steps) then 
                    call overhead_adaptive_time_steps(poly_order, i_scaling, boole_guess, .true., &
                                                        &  iface_new, tau, z, boole_face_correct)
                endif !adaptive steps scheme
!
                call check_three_planes(z,iface_new,boole_face_correct)
                call check_face_convergence(z,iface_new,boole_face_correct)
                call check_exit_time(tau,tau_max,boole_face_correct,poly_order)
!                
                !If previous checks were fine, check for velocity and potentially prolong trajectory (only for 2nd order here)
                if (boole_face_correct) then
                    !Validation for vnorm
                    if(normal_velocity_func(z,iface_new).gt.0.d0) then
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: normal velocity'
                        if(poly_order.gt.2) then
                            boole_face_correct = .false.                       
                        else
!
                            call prolonged_trajectory(poly_order,i_scaling,z,tau,ind_tetr_inout,iface, iface_new, &
                            & boole_face_correct, boole_analytical_approx) 
                            if(.not.boole_analytical_approx) return
                        endif
                    endif !Normal velocity is positive at exit point   
                endif             
!
            endif !boole face correct
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!END OF FIRST ATTEMPT WITH SECOND ORDER GUESS!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(diag_pusher_tetry_poly) print *, 'boole_face_correct',boole_face_correct
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!SECOND ATTEMPT WITHOUT GUESSES PLUS RESCALING IN 2nd ORDER!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(diag_pusher_tetry_poly) print *, 'boole_face_correct',boole_face_correct

            !If error is detected, repeat orbit pushing but without low order polynomial guessing
            if(.not.boole_face_correct) then
!
                boole_faces = .true.    !Compute root for all 4 faces
                boole_face_correct = .true. !allow for consistency checks to set boole_face_correct .false. again
                iface_new = iface_init !instead of iface_new = iface
                z = z_init
                number_of_integration_steps = 0 !we only reset number_of_integration_steps, even though the other step quantities
                                        !contain wrong data (they will afterwards potentially be overwritten or else never read out)
!
                !Analytical calculation of orbit parameter to pass tetrahdron
                if(poly_order.eq.2) then
                    call analytic_approx(poly_order,i_precomp,boole_faces, &
                                        & 1,z,iface_new,tau,boole_analytical_approx)
                else        
                    call analytic_approx(poly_order,i_precomp,boole_faces, &
                                        & i_scaling,z,iface_new,tau,boole_analytical_approx)
                endif
!
                !Analytical result does not exist.
                if(.not.boole_analytical_approx) then
                    print *, 'Error: No analytical solution exists.'
                    ind_tetr_inout = -1
                    iface = -1
                    call manage_intermediate_steps_arrays(1)
                    return
                endif
!
                !Integrate trajectory analytically
                call analytic_integration(poly_order,i_precomp,z,tau)
!
                if (boole_adaptive_time_steps) then
                    !For this section the starting face is the initial face and recaling for second order
                    if (poly_order.eq.2) then  
                        call overhead_adaptive_time_steps(poly_order, 1, .false., .true., &
                                                        &  iface_new, tau, z, boole_face_correct)
                    else
                        call overhead_adaptive_time_steps(poly_order, i_scaling, .false., .true., &
                                                        &  iface_new, tau, z, boole_face_correct)
                    endif
                endif !adaptive steps scheme
!
                !call consistency checks for potential trouble shooting              
                call check_exit_time(tau,tau_max,boole_face_correct,poly_order) 
                call check_three_planes(z,iface_new,boole_face_correct)
                call check_face_convergence(z,iface_new,boole_face_correct)
!
                !If previous checks were fine, check for velocity and potentially prolong trajectory (here for all orders)
                if (boole_face_correct) then
                    !Validation for vnorm
                    if(normal_velocity_func(z,iface_new).gt.0.d0) then
!
if(diag_pusher_tetry_poly) print *, 'prolonged trajectory was called in second attempt'
!
                        call prolonged_trajectory(poly_order,i_scaling,z,tau,ind_tetr_inout,iface, iface_new, boole_face_correct, &
                        & boole_analytical_approx) 
                        if(.not.boole_analytical_approx) return
                    endif !Normal velocity is positive at exit point   
                endif
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!END OF SECOND ATTEMPT WITHOUT GUESSES PLUS RESCALE!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!THIRD ATTEMPT WITH RESCALE + ORDERREDUCTION (TROUBLESHOOTING)!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                !If any of the previous tests did not work (including those within prolonged trajectory), call trouble shooting
                if (.not.boole_face_correct) then
if(diag_pusher_tetry_poly) print*, 'Called Troubleshooting'
                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting)
!                       
                    if(.not.boole_trouble_shooting) then
                        print *, 'Error: Trouble shooting failed. Remove particle.'
                        ind_tetr_inout = -1
                        iface = -1
                        call manage_intermediate_steps_arrays(1)
                        return
                    endif
                endif
!
            endif !.not. boole_face_correct
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!END OF THRID ATTEMP (TROUBLESHOOTING)!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
if(diag_pusher_tetry_poly) then
    print *, 'After orbit pushing'
    print *, 'iface_new',iface_new
    print *, 'norm'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i) 
    enddo
endif
!
            !Final processing
            x=z(1:3)+tetra_physics(ind_tetr)%x1
            vpar=z(4)
!
            !Compute passing time dependent on time tracing option
            select case(i_time_tracing_option)
                !Time tracing in 0th order - dt_dtau = const.
                case(1)
                    t_pass = tau*dt_dtau_const
                !Hamiltonian time tracing with computation of polynomial
                case(2)
!
                    !Track Hamiltonian time for root finding operation
                    allocate(t_hamiltonian_list(number_of_integration_steps+1))
                    t_hamiltonian_list(1) = 0.d0
                    allocate(t_hamiltonian)
                    t_hamiltonian = 0.d0
!
                    !loop over number_of_integration_steps
                    do i = 1,number_of_integration_steps
                        !Calculate Hamiltonian time
                        call calc_t_hamiltonian(poly_order, intermediate_z0_list(:,i), tau_steps_list(i), t_hamiltonian)
!
                        !Track Hamiltonian time
                        t_hamiltonian_list(i+1) = t_hamiltonian
                    enddo
!
                    t_pass = t_hamiltonian
            end select
!
if(diag_pusher_tetry_poly) print *, 'tau total',tau
if(diag_pusher_tetry_poly) print *, 't_pass',t_pass
!if(diag_pusher_tetry_poly) then
!    print *, 't_remain',t_remain
!    if (t_remain .lt. 0) stop
!    if (t_pass .lt. 0) stop
!endif
!
            !Particle stops inside the tetrahedron - Absolute value is used, because negative time can be allowed
            if(abs(t_pass).ge.abs(t_remain)) then
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!FOURTH ATTEMPT IF PARTICLE DOES NOT LEAVE CELL IN REMAINING TIME!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
                !Case selection dependent on time tracing option
                select case(i_time_tracing_option)
                    !Time tracing in 0th order - dt_dtau = const.
                    case(1)
!
                        !Set z back to z_init
                        z = z_init
!
                        !If vnorm was corrected or intermediate steps taken, coefficients need to be computed again.
                        if(number_of_integration_steps .gt. 1) then
                            iface_new = iface_init
!
                            call set_integration_coef_manually(poly_order,z)
                        endif
!
                        !And reset the integration step counter
                        number_of_integration_steps = 0
!
                        !Compute orbit parameter tau from t_remain
                        tau = t_remain/dt_dtau_const
!
                    !Hamiltonian time tracing with computation of polynomial
                    case(2)
!
                        !Find integration step in which t_hamiltonian exceeds t_remain
                        !Absolute value is used, because negative time can be allowed
                        i_step_root =  findloc(abs(t_hamiltonian_list).gt.abs(t_remain),.true.,dim=1)
!
                        !Set orbit back to i-th integration step
                        z = intermediate_z0_list(:,i_step_root - 1)
                        call set_integration_coef_manually(poly_order,z)
                        number_of_integration_steps = i_step_root - 2
                        iface_new = iface_init ! (just in case)
!
                        !Find tau corresponding to root of t_Hamiltonian
                        t_remain_new = t_remain - t_hamiltonian_list(i_step_root-1)
                        call get_t_hamiltonian_root(poly_order,z,t_remain_new,tau)
!
                        !Fail save: consistency check due to order reduction in root solving operation in the case of 5th order
                        if(tau.gt.tau_steps_list(i_step_root - 1)) then
                            print *, 'Warning: Final t_hamiltonian step in 0th order for avoiding inconsistency.'
                            tau = t_remain_new/dt_dtau_const
                        endif
!
!print *, 'Warning: tau', tau, 'tau_steps_list',tau_steps_list(i_step_root - 1)
!print *, 'number_of_integration_steps', number_of_integration_steps
!print *, 'tau_steps_list',tau_steps_list
!print *, 'tau_steps_list * dt_dtau',tau_steps_list*dt_dtau_const
!print *, 't_hamiltonian_list',t_hamiltonian_list
!print *, 't_remain',t_remain
!print *, 't_remain_new',t_remain_new
!print *, 'i_step_root', i_step_root
!print *, 'intermediate_z0_list',intermediate_z0_list
!
                end select
                
if(diag_pusher_tetry_poly) print *, 'tau until t finished',tau
!
                !Integrate trajectory analytically from start until t_remain
                call analytic_integration(poly_order,i_precomp,z,tau)
!
                if (boole_adaptive_time_steps) then
                    !For this section the starting face is the initial face
                    call overhead_adaptive_time_steps(poly_order, i_scaling, .false., .false., &
                                                        &  iface_new, tau, z, boole_face_correct)
                endif !adaptive steps scheme
!
                ind_tetr_inout = ind_tetr
                iface = 0
!
                boole_face_correct = .true.
                !Particle must be inside the tetrahedron
                do i = 1,4
                    if(normal_distance_func(z(1:3),i).lt.0.d0) then
                        boole_face_correct = .false.
                    endif
                enddo
!
                if(boole_face_correct) then
                    boole_t_finished = .true.
                    
                    z_save = z(1:3)
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    vpar=z(4)
!
                    !Difference in between t_pass and t_hamiltonian:
!
                    ! t_pass: expresses the used time of integration including possible error due to order reduction (only relevant
                    ! if poly_order = 4) which is necessary for analytical root finding in the 5th order
                    !
                    ! t_hamiltonian: order consistent computation of Hamiltonian time elapsed
                    ! Due to above described order reduction, t_hamiltonian might sligthly differ from t_step
                    t_pass = t_remain
!                    
                else    !Time step finishes while particle is outside the tetrahedron (Wrong face was predicted, but logically true)
                    !Compute correct root by honestly computing higher order root for all four faces
if(diag_pusher_tetry_poly) then
    print *, 'Error: Particle is outside the tetrahedron when time is finished.'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i)
    enddo
endif
!
                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting)
!                       
                    if(.not.boole_trouble_shooting) then
                        print *, 'Error: Trouble shooting failed. Remove particle.'
                        ind_tetr_inout = -1
                        iface = -1
                        call manage_intermediate_steps_arrays(1)
                        return
                    endif
!
!
                    !Final processing
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    vpar=z(4)
!
                    !Case selection dependent on time tracing option
                    select case(i_time_tracing_option)
                        !Time tracing in 0th order - dt_dtau = const.
                        case(1)
                            t_pass = tau*dt_dtau_const
                        !Hamiltonian time tracing with computation of polynomial
                        case(2)
                            !Compute Hamiltonian time for the final step
                            t_hamiltonian = 0.d0
!
                            !Redo loop over number_of_integration_steps
                            do i = 1,number_of_integration_steps
                                call calc_t_hamiltonian(poly_order, intermediate_z0_list(:,i), tau_steps_list(i), t_hamiltonian)
                            enddo
!
                            t_pass = t_hamiltonian
                    end select
!
                    !Save relative coordinates after pushing
                    z_save = z(1:3)
!
                    iface = iface_new
                    call pusher_handover2neighbour(ind_tetr,ind_tetr_inout,iface,x,iper_phi)
                endif
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!END OF FOURTH ATTEMPT IF PARTICLE DOES NOT LEAVE CELL!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
            !Normal orbit that passes the whole tetrahedron
            else
!            
                !Save relative coordinates after pushing
                z_save = z(1:3)
!
                iface = iface_new
                call pusher_handover2neighbour(ind_tetr,ind_tetr_inout,iface,x,iper_phi)
!                
            endif
!
            if(any(boole_array_optional_quantities)) then
                !loop over number_of_integration_steps
                do i = 1,number_of_integration_steps
                    call calc_optional_quantities(poly_order, intermediate_z0_list(:,i), tau_steps_list(i), optional_quantities)
                enddo
            endif
!
            !If finished particle integration, deallocate the intermediate_steps_array (get also deallocated in previous instances if lost particle)
            if(boole_t_finished) call manage_intermediate_steps_arrays(1)
!
            !Deallocation
            if(i_time_tracing_option.eq.2) deallocate(t_hamiltonian_list,t_hamiltonian)
!
        end subroutine pusher_tetra_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_three_planes(z,iface_new,boole_face_correct)

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            double precision, dimension(4), intent(in)         :: z
            integer, intent(in)                                :: iface_new
            logical                                            :: boole_face_correct
            integer                                            :: j,k

            !Validation loop ('3-planes'-control)
            do j=1,3    !Just consider the planes without the "exit-plane"
            k=modulo(iface_new+j-1,4)+1
                if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error: three planes'
if(diag_pusher_tetry_poly) print *, 'face', k,'normal_distance',normal_distance_func(z(1:3),k)
                endif        
            enddo

        end subroutine check_three_planes
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_face_convergence(z,iface_new,boole_face_correct)

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            double precision, dimension(4), intent(in)         :: z
            integer, intent(in)                                :: iface_new
            logical                                            :: boole_face_correct

        if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
            boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error: distance'
        endif

        
        end subroutine check_face_convergence
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_velocity(z,iface_new,boole_face_correct,poly_order)

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            double precision, dimension(4), intent(in)         :: z
            integer, intent(in)                                :: iface_new
            logical                                            :: boole_face_correct
            integer, intent(in)                                :: poly_order

        if(normal_velocity_func(z,iface_new).gt.0.d0) then
if(diag_pusher_tetry_poly) print *, 'Error: normal velocity'
                boole_face_correct = .false.
        endif

        end subroutine check_velocity
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine check_exit_time(tau,tau_max,boole_face_correct,poly_order) !checks only for poly_order > 2

            use gorilla_diag_mod,only: diag_pusher_tetry_poly

            implicit none
            logical                                            :: boole_face_correct
            integer, intent(in)                                :: poly_order
            double precision                                   :: tau,tau_max

            !Higher order polynomial result for tau is out of safety boundary
            if(poly_order.gt.2) then
                if(tau.gt.tau_max) then
if(diag_pusher_tetry_poly)  print *, 'Error: Tau is out of safety boundary'
                    boole_face_correct = .false.
                endif
            endif
            
        end subroutine check_exit_time
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine prolonged_trajectory(poly_order,i_scaling,z,tau,ind_tetr_inout,iface, iface_new, boole_face_correct, &
            & boole_analytical_approx)
!
        use gorilla_settings_mod, only: i_precomp, boole_adaptive_time_steps
!
        implicit none
!
        integer, intent(in)                                 :: poly_order, i_scaling
        double precision, dimension(4), intent(inout)       :: z
        double precision, intent(inout)                     :: tau
        integer,intent(inout)                               :: ind_tetr_inout,iface, iface_new
        logical, intent(inout)                              :: boole_face_correct
        logical, intent(out)                                :: boole_analytical_approx
!
        logical, dimension(4)                               :: boole_faces
        integer                                             :: iface_new_save
        double precision                                    :: tau_save,tau_max

!                   
        !Particle orbit turns exactly at the face. Find next exit point of same tetrahedron.
        boole_faces = .true.
        tau_save = tau
        iface_new_save = iface_new !saving is done to retain the value for second analytic_approx

!                   
        if (poly_order.gt.2) then
            !Analytical calculation of orbit parameter for safety boundary
            call analytic_approx(2,i_precomp,boole_faces, &
                                & i_scaling,z,iface_new,tau,boole_analytical_approx)
            !Define boundary for tau_max
            tau_max = tau*eps_tau
        endif  
!                    
        iface_new = iface_new_save   !<-- Jonatan, Georg, 07.07.2022: above, iface_new is changed, this must be undone
        call analytic_approx(poly_order,i_precomp,boole_faces, &
                            & i_scaling,z,iface_new,tau,boole_analytical_approx)                                   
!                    
        !Analytical result does not exist.
        if(.not.boole_analytical_approx) then
            print *, 'Error in prolonged trajectory: No analytical solution exists.'
            ind_tetr_inout = -1
            iface = -1
            call manage_intermediate_steps_arrays(1)
            return
        endif
!                    
        !Integrate trajectory analytically
        call analytic_integration(poly_order,i_precomp,z,tau)
!
        if (boole_adaptive_time_steps) then
            call overhead_adaptive_time_steps(poly_order, i_scaling, .false., .true., &
                                                &  iface_new, tau, z, boole_face_correct)
        endif !adaptive steps scheme
!
        !Only executed if poly_order > 2, checks if exit time is within upper limit
        call check_exit_time(tau,tau_max,boole_face_correct,poly_order)
        !Validation loop ('3-planes'-control) 
        call check_three_planes(z,iface_new,boole_face_correct)        
        call check_velocity(z,iface_new,boole_face_correct,poly_order)
        call check_face_convergence(z,iface_new,boole_face_correct)

        !Integration in two parts (Therefore time must be added)
        tau = tau + tau_save !potentially alter this for energy adaption

    end subroutine prolonged_trajectory
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine overhead_adaptive_time_steps(poly_order, i_scaling, boole_guess_adaptive, boole_passing, &
                                            &  iface_inout_adaptive, tau, z, boole_face_correct)
!
        use gorilla_settings_mod, only: desired_delta_energy, max_n_intermediate_steps
        use gorilla_diag_mod, only: diag_pusher_tetry_poly, report_pusher_tetry_poly_adaptive
        use supporting_functions_mod, only: energy_tot_func
!
        implicit none
!
        integer, intent(in)                                   :: poly_order, i_scaling
        logical, intent(in)                                   :: boole_guess_adaptive
        logical, intent(in)                                   :: boole_passing
!
        integer, intent(inout)                                :: iface_inout_adaptive
        double precision, intent(inout)                       :: tau
        double precision, dimension(4), intent(inout)         :: z
        logical, intent(inout)                                :: boole_face_correct
!
        double precision                                      :: energy_start, energy_current, delta_energy_current, tau_save_report
        double precision, dimension(4)                        :: z0, z_save_report
        integer                                               :: number_of_integration_steps_save_report, &
                                                                & iface_inout_adaptive_save_report
        logical                                               :: boole_face_correct_save_report                        
!
        if ((desired_delta_energy .le. 0)) then
            print*, 'Error: The control setting desired_delta_energy is invalid! Check the limits in gorilla.inp!'
            stop
        elseif ((max_n_intermediate_steps .lt. 2)) then
            print*, 'Error: The control setting max_n_intermediate_steps is invalid! Check the limits in gorilla.inp!'
            stop
        endif
!
if(diag_pusher_tetry_poly) print*, 'Sanity check before checking energy conservation:'
if(diag_pusher_tetry_poly) print*, 'z', z
        !At least the check of proper inside (three planes + convergence in case of passing) have to be fullfilled
        if (boole_passing) then
            call check_three_planes(z,iface_inout_adaptive,boole_face_correct)
            call check_face_convergence(z,iface_inout_adaptive,boole_face_correct)
        elseif (.not.boole_passing) then
            call check_three_planes(z,0,boole_face_correct)
        else
            print*, 'Error: Non valid value for boole_passing of overhead_adaptive_time_step()!' 
            stop
        endif
        if (.not.boole_face_correct) return
!
        z0 = intermediate_z0_list(:,number_of_integration_steps)
        energy_start = energy_tot_func(z0,perpinv,ind_tetr)
        energy_current = energy_tot_func(z,perpinv,ind_tetr)
        delta_energy_current = abs(1-energy_current/energy_start)
!
if(diag_pusher_tetry_poly) print*, 'z0', z0
if(diag_pusher_tetry_poly) print*, 'energy_current', energy_current
if(diag_pusher_tetry_poly) print*, 'energy_start', energy_start
if(report_pusher_tetry_poly_adaptive) then
!Save the non-adaptive end results -> in report mode, the actual orbit is not changed by th adaptive scheme!
z_save_report = z
tau_save_report = tau
number_of_integration_steps_save_report = number_of_integration_steps
boole_face_correct_save_report = boole_face_correct
iface_inout_adaptive_save_report = iface_inout_adaptive
endif
!
        !If energy fluctuating too strong -> start recursive adaptive scheme
        if (delta_energy_current .gt. desired_delta_energy) then
if(diag_pusher_tetry_poly) print*, 'Adaptive Stepsize equidistant was called'
            call adaptive_time_steps_equidistant(poly_order, i_scaling, boole_guess_adaptive, boole_passing, &
                                & delta_energy_current, iface_inout_adaptive, tau, z, boole_face_correct)
        endif
!
if(report_pusher_tetry_poly_adaptive) then
!Restore non-adaptive results (report mode does not change orbit!)
z = z_save_report
tau = tau_save_report
number_of_integration_steps = number_of_integration_steps_save_report
tau_steps_list(number_of_integration_steps) = tau
boole_face_correct = boole_face_correct_save_report
iface_inout_adaptive = iface_inout_adaptive_save_report
endif
!
    end subroutine overhead_adaptive_time_steps
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine adaptive_time_steps_equidistant(poly_order, i_scaling, boole_guess_adaptive, boole_passing, &
        &  delta_energy_current, iface_out_adaptive, tau, z, boole_face_correct)
!
        use gorilla_settings_mod, only: i_precomp, desired_delta_energy, max_n_intermediate_steps
        use gorilla_diag_mod,only: diag_pusher_tetry_poly_adaptive, report_pusher_tetry_poly_adaptive
        use supporting_functions_mod, only: energy_tot_func
!
        implicit none
!
        integer, intent(in)                                   :: poly_order, i_scaling
        logical, intent(in)                                   :: boole_guess_adaptive
        logical, intent(in)                                   :: boole_passing
!
        double precision, intent(inout)                       :: tau, delta_energy_current
        double precision, dimension(4), intent(inout)         :: z
        logical, intent(inout)                                :: boole_face_correct
!
        integer, intent(out)                                  :: iface_out_adaptive
!
        logical, dimension(4)                                 :: boole_faces
        double precision, dimension(4)                        :: z_start_adaptive
        integer                                               :: i, k, eta, eta_extended, eta_limit,&
                                                                & eta_minimum, eta_buffer
        integer                                               :: iface_new_adaptive, number_of_integration_steps_start_adaptive
        logical                                               :: boole_analytical_approx, boole_exit_tetrahedron, & 
                                                                & boole_energy_check, boole_reached_minimum
        double precision                                      :: energy_start_adaptive, energy_current, &
                                                                & delta_energy_start, delta_energy_minimum, & 
                                                                & tau_prime, tau_collected, tau_exit, tau_minimum, tau_buffer
        double precision                                      :: delta_energy_step, energy_current_step, energy_previous_step
        !Diagnostics for adaptive step scheme (only useable for one particle calculation)
        double precision, dimension(:,:), allocatable         :: total_fluc_report, step_fluc_report
        integer                                               :: report_index
!
        !Save current z0 and stepped back global integration counter
        z_start_adaptive = intermediate_z0_list(:,number_of_integration_steps)
        number_of_integration_steps_start_adaptive = number_of_integration_steps - 1
        energy_start_adaptive = energy_tot_func(z_start_adaptive,perpinv,ind_tetr)
if(diag_pusher_tetry_poly_adaptive) print*, '------------------------'
!
        !Set up for Partition procedure
        eta = 1
        eta_extended = 1
        boole_reached_minimum = .false.
        delta_energy_start = delta_energy_current
        delta_energy_minimum = delta_energy_current
!
if (report_pusher_tetry_poly_adaptive) then
!Intialize report quantities for this instance of adaptive scheme
if((.not. allocated(step_fluc_report)) .OR. (.not. allocated(total_fluc_report))) then
allocate(total_fluc_report(max_n_intermediate_steps+100,2), step_fluc_report(max_n_intermediate_steps+100,2))
endif
report_index = 1
step_fluc_report(report_index,1) = eta
step_fluc_report(report_index,2) = delta_energy_current
total_fluc_report(report_index,1) = eta
total_fluc_report(report_index,2) = delta_energy_current
endif
!
        !Loop over possible equidistant splittings (eta changing depending on the current energy error, with uppper limit a priori set)
        PARTITION: do while (eta .lt. max_n_intermediate_steps)
!
            !If did not succeed (aka not left the loop) -> try to update number of steps
            !and try again, if not yet at maximal number of intermediate steps or passed minimum
            call adaptive_time_steps_update_eta(poly_order,delta_energy_current,eta)
!
            !Set up for redo of the minimum run
            if (boole_reached_minimum) then
                eta = eta_minimum
                tau = tau_minimum
            endif
!
            !Set up for every partition trial; boole_face_correct can be set true here as adaptive is before the consistency checks
            z = z_start_adaptive
            number_of_integration_steps = number_of_integration_steps_start_adaptive
            tau_prime = tau/eta
if(diag_pusher_tetry_poly_adaptive) print*, 'eta', eta
if(diag_pusher_tetry_poly_adaptive) print *, 'tau_prime', tau_prime
            tau_collected = 0
            boole_face_correct = .true.
            boole_exit_tetrahedron = .false.
            boole_energy_check = .false.
            if (boole_passing) then
                !As the intermediate steps change our orbit, we may need longer than the original tau -> more steps 
                eta_limit = ceiling(max_n_intermediate_steps*1.1d0)
            else
                !We do not expect a passing
                eta_limit = eta
            endif !boole_passing for eta_limit
!
            STEPWISE: do i = 1, (eta_limit -1) !used to be 1,eta-1
                !recalculate polynomial coefficients (tensors) as at every intermediate step the poly_coef change
                call set_integration_coef_manually(poly_order,z)
                call analytic_integration(poly_order,i_precomp,z,tau_prime)
!
                !GC must be still inside the tetrahedron
                !If fail -> return to last step and close orbit then
                CONTROL: do k = 1,4
                    if(normal_distance_func(z(1:3),k).lt.0.d0) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Adaptive steps: Position left tetrahedron'
if(diag_pusher_tetry_poly_adaptive) print *, k, 'norm', normal_distance_func(z(1:3),k)
if(diag_pusher_tetry_poly_adaptive) print *, 'steps taken', i, '/' ,eta
                        !Exception: if the first step is already outside (order-inconsistency if used guess) -> return failure
                        if (i .eq. 1) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps: Left tetrahedron in first step (order-inconsistency)!'
                            boole_face_correct = .false.
                            return
                        endif
                        z = intermediate_z0_list(:,number_of_integration_steps)
                        number_of_integration_steps = number_of_integration_steps - 1
                        boole_exit_tetrahedron = .true.
                        exit STEPWISE
                    endif
                end do CONTROL !Is the orbit still inside
!
                !If it was a valid step we add it to the total used time
                tau_collected = tau_collected + tau_prime
                eta_extended = i
!
if (report_pusher_tetry_poly_adaptive) then
!Collecting single step fluctuation data
if (i.eq.1) then
report_index = report_index + 1
step_fluc_report(report_index,2) = 0
energy_current_step = energy_start_adaptive
endif
energy_previous_step = energy_current_step
energy_current_step = energy_tot_func(z,perpinv,ind_tetr)
delta_energy_step = abs(1 - energy_current_step/energy_previous_step)
step_fluc_report(report_index,2) = step_fluc_report(report_index,2) + delta_energy_step
endif
!
                !The only reason to not use the guessing scheme is to avoid nipping out of orbits (ensure order-consistency)
                !For that we additionally need to check the remaining tau after each step when using the adaptive scheme
                if(.not.boole_guess_adaptive) then
                    iface_new_adaptive = 0
                    call adaptive_time_steps_exit_time(poly_order,i_scaling,z,.false., & 
                                                            & iface_new_adaptive,tau_exit,boole_analytical_approx)
                    !Analytical result does not exist can therefore not close cell orbit
                    ! -> "return" leaves orbit (falsely, therefore boole set to false) inside of tetrahedron
                    if(.not.boole_analytical_approx) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps: no analytical solution'
                        boole_face_correct = .false.
                        return
                    endif
                    if(tau_exit .le. tau_prime) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Adaptive step order-consistency: tau_exit',tau_exit,'smaller than tau_prime',tau_prime
                        tau_prime = tau_exit
                        exit STEPWISE !boole_exit_tetrahedron is still FALSE, as avoid leaving tetrahedron all together
                    endif
                endif !oberving remaining tau to avoid nipping out of orbits
!
            end do STEPWISE ! stepwise integration
!
            !If orbit exited, the now remaining time (tau_prime) from the last intermediate point to the exit face has to be recalculated
            !Alternatively stepwise integration might reveal a not passing orbit to be a passing one, which also has to be closed
            if (boole_exit_tetrahedron) then
                !As we already preformed at least one step before leaving the tetrahedron, the orbit is now INSIDE
                iface_new_adaptive = 0
                call adaptive_time_steps_exit_time(poly_order,i_scaling,z,boole_guess_adaptive, & 
                                                        & iface_new_adaptive,tau_exit,boole_analytical_approx)           
                !Analytical result does not exist can therefore not close cell orbit
                ! -> "return" leaves orbit (falsely, therefore boole set to false) inside of tetrahedron
                if(.not.boole_analytical_approx) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps: no analytical solution'
                    boole_face_correct = .false.
                    return
                endif
if(diag_pusher_tetry_poly_adaptive) print *, 'Adaptive step closure: tau_exit', tau_exit
                tau_prime = tau_exit !The last step will close the orbit to the exit
            else
                !If the orbit does not pass, the same timestep as for the other intermediate steps is used
                !Should it then end up outside, the higher level checks will detect it
                !We only have to manually update again the new position in the coefficiants, no new tau_prime needed
                call set_integration_coef_manually(poly_order,z)
            end if !boole_passing
!
            !Closing integration (either to the exit face or to the final position inside of tetrahedron)
            call analytic_integration(poly_order,i_precomp,z,tau_prime)
            tau_collected = tau_collected + tau_prime
!
            !Check if energy fluctuation was successfully decreased
            energy_current = energy_tot_func(z,perpinv,ind_tetr)
            delta_energy_current = abs(1-energy_current/energy_start_adaptive)
!
            !Update tau/eta actually needed to travers tetrahedron (for final result or for better approx in another iteration)
            !The buffer variables are used to 1:1 save the minimum run, as it was achieved now (even if the tau/eta was off, only the result matters not how we got to it!)
            eta_buffer = eta
            tau_buffer = tau
            eta = eta_extended + 1
            tau = tau_collected
!
if (report_pusher_tetry_poly_adaptive) then
!Finish up on data of this partition (average of step-fluctuation)
!And circumvent the exit mechanisms to collect more data
boole_energy_check = .true.
energy_previous_step = energy_current_step
energy_current_step = energy_tot_func(z,perpinv,ind_tetr)
delta_energy_step = abs(1 - energy_current_step/energy_previous_step)
step_fluc_report(report_index,2) = step_fluc_report(report_index,2) + delta_energy_step
step_fluc_report(report_index,1) = eta_extended + 1 
step_fluc_report(report_index,2) = step_fluc_report(report_index,2) / (eta_extended + 1)
total_fluc_report(report_index,2) = delta_energy_current
total_fluc_report(report_index,1) = eta_extended + 1
eta = eta_buffer
tau = tau_buffer
cycle PARTITION
endif
!
            !There is an effective minimum that can be achieved by decreasing the timesteps
            !After that the error actually increases/osscilates -> we only use the monoton decrease
if (diag_pusher_tetry_poly_adaptive) print*, 'delta_energy_current', delta_energy_current
if (diag_pusher_tetry_poly_adaptive) print*, 'delta_energy_minimum', delta_energy_minimum
            if (boole_reached_minimum) then
                exit PARTITION
            elseif (delta_energy_current.lt.delta_energy_minimum) then
                delta_energy_minimum = delta_energy_current
                eta_minimum = eta_buffer
                tau_minimum = tau_buffer
                if (delta_energy_current.le.desired_delta_energy) then
                    boole_energy_check = .true.
                    exit PARTITION
                endif
            elseif(delta_energy_current.gt.delta_energy_minimum) then
                !If we see that our new step choice has INCREASED the energy, we redo the minimum partition and stop
                boole_reached_minimum = .true.
                cycle PARTITION
            endif !energy check
!
        end do PARTITION
!
        iface_out_adaptive = iface_new_adaptive
!
        !If could not satisfy energy conservation with scheme, notify user
        if (.not.boole_energy_check) then
            print *, 'Error in adaptive steps equidistant: energy conservation could not be fullfilled!'
            print *, 'delta_energy_minimum', delta_energy_minimum
            print *, 'eta_minimum', eta_minimum
            print *, 'tau_minimum', tau_minimum
if(diag_pusher_tetry_poly_adaptive) stop
        endif 
!
if (report_pusher_tetry_poly_adaptive) then
!Write data of this instance of adaptive scheme into report files
report_index = report_index + 1
step_fluc_report(report_index,:) = -1
total_fluc_report(report_index,:) = -1
do k = 1, report_index
write(124,*) total_fluc_report(k,:)
end do
do k = 1, report_index
write(118,*) step_fluc_report(k,:)
end do
deallocate(total_fluc_report,step_fluc_report)
endif
!
    end subroutine adaptive_time_steps_equidistant
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine adaptive_time_steps_update_eta(poly_order,delta_energy_current,eta)
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly_adaptive, report_pusher_tetry_poly_adaptive
        use gorilla_settings_mod, only: desired_delta_energy, max_n_intermediate_steps
!
        implicit none
!
        integer, intent(in)                         :: poly_order
        double precision, intent(in)                :: delta_energy_current
!
        integer, intent(inout)                      :: eta
!
        double precision                            :: scale_factor, max_scale_factor
        double precision, parameter                 :: threshold = 1.d0, min_step_error = 1E-15, additive_increase = 1 !default 1,1E-15,10
!
if (report_pusher_tetry_poly_adaptive)  then 
    !To get more data points, do not optimal scaling
    if (delta_energy_current/eta .lt. min_step_error) then
        eta = ceiling(eta*2.d0**(1.d0/4.d0))
    else
        eta = ceiling(eta*2.d0**(1.d0/2.d0))
    endif
    eta = min(eta,max_n_intermediate_steps)
    return
endif
        scale_factor = (delta_energy_current/desired_delta_energy)**(1.0d0/poly_order)
        max_scale_factor = (delta_energy_current/(min_step_error*eta))**(1.0d0/(poly_order+1))
if (diag_pusher_tetry_poly_adaptive) print*, 'scale factor', scale_factor
if (diag_pusher_tetry_poly_adaptive) print*, 'max scale factor', max_scale_factor
        if ((scale_factor .gt. threshold) .AND. (max_scale_factor .gt. threshold)) then
            scale_factor = min(scale_factor,max_scale_factor)
            eta = min(ceiling(eta*scale_factor),max_n_intermediate_steps)
        elseif (scale_factor .gt. threshold) then
            eta = min(ceiling(eta*scale_factor),max_n_intermediate_steps)
        else
            eta = eta + additive_increase
        endif !Choice of next number of steps
!
    end subroutine adaptive_time_steps_update_eta
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    subroutine adaptive_time_steps_exit_time(poly_order,i_scaling,z,boole_guess_adaptive, & 
                                                & iface_new_adaptive,tau_exit,boole_analytical_approx)
!
        use gorilla_settings_mod, only: i_precomp
        use gorilla_diag_mod, only: diag_pusher_tetry_poly_adaptive
!
        implicit none
!
        integer, intent(in)                                 :: poly_order,i_scaling
        double precision, dimension(4),intent(in)           :: z
        logical, intent(in)                                 :: boole_guess_adaptive
!
        integer,intent(inout)                               :: iface_new_adaptive
!
        double precision, intent(out)                       :: tau_exit
        logical,intent(out)                                 :: boole_analytical_approx
!
        logical, dimension(4)                               :: boole_faces
!
        boole_faces = .true.
!
        !use face prediction of second order for higher order computation in the adaptive scheme if wished
        if (boole_guess_adaptive .and. (poly_order.gt.2)) then
            call analytic_approx(2,i_precomp,boole_faces, &
            & i_scaling,z,iface_new_adaptive,tau_exit,boole_analytical_approx)
            !Higher order polynomial root is computed only for previously predicted face in second order
            if(boole_analytical_approx) then               
                boole_faces = .false.                !Disable all 4 faces 
                boole_faces(iface_new_adaptive) = .true.      !Enable guessed face
            endif
            !Reset starting face (to INSIDE) for actual correct order calculation down below
            iface_new_adaptive = 0 
        endif  
        !calculate exit time and exit face in correct order
        !if a successful guess was made above in second order, boole_faces only allows the guessed face 
        call analytic_approx(poly_order,i_precomp,boole_faces, &
                            & i_scaling,z,iface_new_adaptive,tau_exit,boole_analytical_approx)   
!   
    end subroutine adaptive_time_steps_exit_time
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine analytic_approx(poly_order,i_precomp,boole_faces,i_scaling,z,iface_inout,dtau,boole_approx)
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
!
        implicit none
!
        integer, intent(in)                                 :: poly_order,i_precomp
        integer,intent(inout)                               :: iface_inout
        logical, dimension(4), intent(in)                   :: boole_faces
        double precision,intent(out)                        :: dtau
        double precision, dimension(4),intent(in)           :: z
        logical,intent(out)                                 :: boole_approx
        integer                                             :: i,iface,i_scaling
        double precision, dimension(4)                      :: dtau_vec
        logical, dimension(4)                               :: boole_valid_dtau
        double precision, dimension(4,poly_order+1)         :: coef_mat
!
        !Additional quantities for coefs = 0
        integer                                             :: poly_order_solver
        double precision                                    :: quart_a, quart_b, quart_c, quart_d, quart_e
        double precision                                    :: cub_a, cub_b, cub_c, cub_d
        double precision                                    :: quad_a, quad_b, quad_c
        double precision                                    :: lin_a, lin_b
!
        !Calculation of root of polynomial equation
        select case(i_precomp)
            case(0)
                call analytic_coeff_without_precomp(poly_order,boole_faces,z,coef_mat)
            case(1,2)
                call analytic_coeff_with_precomp(poly_order,i_precomp,boole_faces,z,coef_mat)
        end select
!     
        iface = iface_inout   
!
        dtau_vec = huge(0.d0) !Only results with positive times smaller than 'huge' will be evaluated
!
        !loop over faces
        faces_loop: do i = 1,4
            if(.not.boole_faces(i)) cycle
!
            !Initialization of poly_order_solver
            poly_order_solver = 0

            !Classification of coefficients
            select case(poly_order)
                case(1)
                    if( (i.eq.iface).or.(coef_mat(i,1).eq.0.d0) ) then
                        dtau_vec(i) = 0.d0
                        cycle !No further computations are needed below
                    else
                        poly_order_solver = 1
                        lin_a = coef_mat(i,2)
                        lin_b = coef_mat(i,1)
!
                        !Handling, if coefficients are exactly ZERO
                        if(lin_a.eq.0.d0) then
                            dtau_vec(i) = 0.d0
                            cycle
                        endif ! lin_a.eq.0.d0
!
                    endif
                case(2)
                    if( (i.eq.iface).or.(coef_mat(i,1).eq.0.d0) ) then ! coef_mat(i,1) = 0.d0 --> Lowest coefficient vanishes
                        poly_order_solver = 1
                        lin_a = coef_mat(i,3)/2.d0
                        lin_b = coef_mat(i,2)
!
                        !Handling, if coefficients are exactly ZERO
                        if(lin_a.eq.0.d0) then
                            dtau_vec(i) = 0.d0
                            cycle
                        endif ! lin_a.eq.0.d0
!
                    else !coef_mat(i,1) != 0.d0
                        poly_order_solver = 2
                        quad_a = coef_mat(i,3)
                        quad_b = coef_mat(i,2)
                        quad_c = coef_mat(i,1)
!
                        !Handling, if coefficients are exactly ZERO
                        if(quad_a.eq.0.d0) then
                            if(quad_b.ne.0.d0) then
                                poly_order_solver = 1
                                lin_a = quad_b
                                lin_b = quad_c
                            else
                                dtau_vec(i) = 0.d0
                                cycle
                            endif
                        endif ! quad_a.eq.0.d0
!
                    endif !coef_mat(i,1) = 0.d0
!
                case(3)        
                    if( (i.eq.iface).or.(coef_mat(i,1).eq.0.d0) ) then ! coef_mat(i,1) = 0.d0 --> Lowest coefficient vanishes
                        poly_order_solver = 2
                        quad_a = coef_mat(i,4)/3.d0
                        quad_b = coef_mat(i,3)/2.d0
                        quad_c = coef_mat(i,2)
!
                        !Handling, if coefficients are exactly ZERO
                        if(quad_a.eq.0.d0) then
                            if(quad_b.ne.0.d0) then
                                poly_order_solver = 1
                                lin_a = quad_b
                                lin_b = quad_c
                            else
                                dtau_vec(i) = 0.d0
                                cycle
                            endif
                        endif ! quad_a.eq.0.d0
!
                    else !coef_mat(i,1) != 0.d0
                        poly_order_solver = 3
                        cub_a = coef_mat(i,4)
                        cub_b = coef_mat(i,3)
                        cub_c = coef_mat(i,2)
                        cub_d = coef_mat(i,1)
!
                        !Handling, if coefficients are exactly ZERO
                        if(cub_a.eq.0.d0) then
                            if(cub_b.ne.0.d0) then
                                poly_order_solver = 2
                                quad_a = cub_b
                                quad_b = cub_c
                                quad_c = cub_d
                            elseif(cub_c.ne.0.d0) then
                                poly_order_solver = 1
                                lin_a = cub_c
                                lin_b = cub_d
                            else
                                dtau_vec(i) = 0.d0
                                cycle
                            endif
                        endif ! cub_a.eq.0.d0
!
                    endif !coef_mat(i,1) = 0.d0
!
                case(4) 
                    if( (i.eq.iface).or.(coef_mat(i,1).eq.0.d0) ) then ! coef_mat(i,1) = 0.d0 --> Lowest coefficient vanishes
                        poly_order_solver = 3
                        cub_a = coef_mat(i,5)/4.d0
                        cub_b = coef_mat(i,4)/3.d0
                        cub_c = coef_mat(i,3)/2.d0
                        cub_d = coef_mat(i,2)
!
                        !Handling, if coefficients are exactly ZERO
                        if(cub_a.eq.0.d0) then
                            if(cub_b.ne.0.d0) then
                                poly_order_solver = 2
                                quad_a = cub_b
                                quad_b = cub_c
                                quad_c = cub_d
                            elseif(cub_c.ne.0.d0) then
                                poly_order_solver = 1
                                lin_a = cub_c
                                lin_b = cub_d
                            else
                                dtau_vec(i) = 0.d0
                                cycle
                            endif
                        endif ! cub_a.eq.0.d0
!
                    else !coef_mat(i,1) != 0.d0
                        poly_order_solver = 4
                        quart_a = coef_mat(i,5)
                        quart_b = coef_mat(i,4)
                        quart_c = coef_mat(i,3)
                        quart_d = coef_mat(i,2)
                        quart_e = coef_mat(i,1)
!
                        !Handling, if coefficients are exactly ZERO
                        if(quart_a.eq.0.d0) then
                            if(quart_b.ne.0.d0) then
                                poly_order_solver = 3
                                cub_a = quart_b
                                cub_b = quart_c
                                cub_c = quart_d
                                cub_d = quart_e
                            elseif(quart_c.ne.0.d0) then
                                poly_order_solver = 2
                                quad_a = quart_c
                                quad_b = quart_d
                                quad_c = quart_e
                            elseif(quart_d.ne.0.d0) then
                                poly_order_solver = 1
                                lin_a = quart_d
                                lin_b = quart_e
                            else
                                dtau_vec(i) = 0.d0
                                cycle
                            endif
                        endif ! quart_a.eq.0.d0
!
                    endif !coef_mat(i,1) = 0.d0
            end select
!
            !Computation of roots with appropriate solver
            select case(poly_order_solver)
                case(1)
                    call Linear_Solver(lin_a,lin_b,dtau_vec(i))
                case(2)
                    call Quadratic_Solver(i_scaling,quad_a,quad_b,quad_c,dtau_vec(i))
                case(3)
                    call Cubic_Solver(cub_a,cub_b,cub_c,cub_d,dtau_vec(i))
                case(4)
                    call Quartic_Solver(i_scaling,quart_a,quart_b,quart_c,quart_d,quart_e,dtau_vec(i))
            end select
!
        enddo faces_loop
!
        boole_valid_dtau = (dtau_vec.lt.huge(0.d0)).and.(dtau_vec.gt.0.d0) !dtau_max
if(diag_pusher_tetry_poly) print *, 'boole_valid_dtau',boole_valid_dtau
        if( any( boole_valid_dtau) ) then
            boole_approx = .true.
!
            iface_inout = minloc(dtau_vec,1,boole_valid_dtau)
            dtau = dtau_vec(iface_inout)
        else
            boole_approx = .false.
        endif
        
if(diag_pusher_tetry_poly) print *, 'boole',boole_approx,'dtau',dtau,'iface_new',iface_inout
! 
    end subroutine analytic_approx
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_coeff_without_precomp(poly_order,boole_faces,z,coef_mat)
!
            use tetra_physics_mod, only: tetra_physics,cm_over_e
            use gorilla_settings_mod, only: boole_strong_electric_field
            use constants, only: clight
            use poly_without_precomp_mod
!
            implicit none
!
            integer, intent(in)                 :: poly_order
            logical, dimension(4), intent(in)   :: boole_faces
            logical, dimension(4)               :: boole_faces_not
            integer                             :: n
            double precision, dimension(4)      :: z
            double precision                    :: dist1
            double precision, dimension(4,poly_order+1)      :: coef_mat
!
            !ODE coefficients
            b(1:3)= ( tetra_physics(ind_tetr)%curlh*(k1) &
                    + perpinv*tetra_physics(ind_tetr)%gBxh1 )*cm_over_e &
                    - clight*(2.d0*(k3)*tetra_physics(ind_tetr)%curlh &
                    + tetra_physics(ind_tetr)%gPhixh1)
!
            b(4)=perpinv*tetra_physics(ind_tetr)%gBxcurlA-clight/cm_over_e*tetra_physics(ind_tetr)%gPhixcurlA
!
            amat = 0.d0
!
            amat(1:3,1:3)=perpinv*cm_over_e*tetra_physics(ind_tetr)%alpmat &
                        - clight* tetra_physics(ind_tetr)%betmat
            amat(4,4) = perpinv*cm_over_e*tetra_physics(ind_tetr)%spalpmat &
                        - clight* tetra_physics(ind_tetr)%spbetmat
            amat(1:3,4) = tetra_physics(ind_tetr)%curlA
!
            if (boole_strong_electric_field) then
                b(1:3) = b(1:3) - 0.5d0*cm_over_e*tetra_physics(ind_tetr)%gv2Emodxh1 !and also the modification of k1 done by initialization
                b(4) = b(4) + cm_over_e*perpinv*tetra_physics(ind_tetr)%gBxcurlvE - clight*tetra_physics(ind_tetr)%gPhixcurlvE &
                            & - 0.5d0*cm_over_e*tetra_physics(ind_tetr)%gv2EmodxcurlvE - 0.5d0*tetra_physics(ind_tetr)%gv2EmodxcurlA
                amat(1:3,1:3) = amat(1:3,1:3) - 0.5d0*cm_over_e*tetra_physics(ind_tetr)%gammat
                amat(4,4) = amat(4,4) - 0.5d0*cm_over_e*tetra_physics(ind_tetr)%spgammat
                amat(1:3,4) = amat(1:3,4) + cm_over_e*tetra_physics(ind_tetr)%curlvE
            endif !boole_stron_electric_fields (adding additional terms to the coefficiants)
!
            !Multiply amat and b with appropriate sign (which ensures that tau remains positive inside the algorithm)
            amat = amat * dble(sign_rhs)
            b = b * dble(sign_rhs)
!
            dist1= -tetra_physics(ind_tetr)%dist_ref
!
            boole_faces_not = .not.boole_faces
!
            if(poly_order.ge.0) then
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,1)=sum(tetra_physics(ind_tetr)%anorm(:,n)*z(1:3))
                enddo
                coef_mat(1,1)=coef_mat(1,1)-dist1
            endif
!
            if(poly_order.ge.1) then
                amat_in_z = matmul(amat,z)
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,2) = sum(tetra_physics(ind_tetr)%anorm(:,n)*amat_in_z(1:3)) + &
                                    & sum(tetra_physics(ind_tetr)%anorm(:,n)*b(1:3))
                enddo
            endif
!
            if(poly_order.ge.2) then
                amat2 = matmul(amat,amat)
                amat2_in_z = matmul(amat2,z)
                amat_in_b = matmul(amat,b)
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,3) = sum(tetra_physics(ind_tetr)%anorm(:,n)*amat2_in_z(1:3)) + &
                                    & sum(tetra_physics(ind_tetr)%anorm(:,n)*amat_in_b(1:3))
                enddo
            endif
!
            if(poly_order.ge.3) then
                amat3 = matmul(amat,amat2)
                amat3_in_z = matmul(amat3,z)
                amat2_in_b = matmul(amat2,b)
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,4) = sum(tetra_physics(ind_tetr)%anorm(:,n)*amat3_in_z(1:3)) + &
                                    & sum(tetra_physics(ind_tetr)%anorm(:,n)*amat2_in_b(1:3))
                enddo
            endif
!
            if(poly_order.ge.4) then
                amat4 = matmul(amat,amat3)
                amat4_in_z = matmul(amat4,z)
                amat3_in_b = matmul(amat3,b)
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,5) = sum(tetra_physics(ind_tetr)%anorm(:,n)*amat4_in_z(1:3)) + &
                                    & sum(tetra_physics(ind_tetr)%anorm(:,n)*amat3_in_b(1:3))
                enddo
            endif
!
        end subroutine analytic_coeff_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_coeff_with_precomp(poly_order,i_precomp,boole_faces,z,coef_mat)
!
            use tetra_physics_mod, only: tetra_physics,cm_over_e
            use constants, only: clight
            use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4
            use poly_without_precomp_mod, only: b,perpinv3,perpinv4
!
            implicit none
!
            integer, intent(in)                 :: poly_order,i_precomp
            logical, dimension(4), intent(in)   :: boole_faces
            logical, dimension(4)               :: boole_faces_not
            integer                             :: n
            double precision, dimension(4)      :: z
            double precision                    :: dist1
            double precision, dimension(4,poly_order+1)      :: coef_mat
!
            if(i_precomp.eq.1) then
                b(1:3)= ( tetra_physics(ind_tetr)%curlh*(k1) &
                        + perpinv*tetra_physics(ind_tetr)%gBxh1 )*cm_over_e &
                        - clight*(2.d0*(k3)*tetra_physics(ind_tetr)%curlh &
                        + tetra_physics(ind_tetr)%gPhixh1)
!
                b(4)=perpinv*tetra_physics(ind_tetr)%gBxcurlA-clight/cm_over_e*tetra_physics(ind_tetr)%gPhixcurlA
            endif
!
            dist1= -tetra_physics(ind_tetr)%dist_ref
!
            boole_faces_not = .not.boole_faces
!
            if(poly_order.ge.0) then
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    coef_mat(n,1)=sum(tetra_physics(ind_tetr)%anorm(:,n)*z(1:3))
                enddo
                coef_mat(1,1)=coef_mat(1,1)-dist1
            endif
!
            if(poly_order.ge.1) then
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    select case(i_precomp)
                        case(1)
                            coef_mat(n,2) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)) * z) + &
                                            & sum(tetra_physics(ind_tetr)%anorm(:,n) * b(1:3))
                        case(2)
                            coef_mat(n,2) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)) * z) + &
!
                                            & tetra_physics_poly4(ind_tetr)%anorm_in_b0(n) + &
                                            & k1 * tetra_physics_poly4(ind_tetr)%anorm_in_b1(n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_b2(n) + &
                                            & k3 * tetra_physics_poly4(ind_tetr)%anorm_in_b3(n)
                    end select
                enddo
            endif
!
            if(poly_order.ge.2) then
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    select case(i_precomp)
                        case(1)
                            coef_mat(n,3) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,n)) * z) + &
!
                                            & sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)) * b)
                        case(2)
                            coef_mat(n,3) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,n)) * z) + &
!
                                            & tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b0(n) + &
                                            & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b0(n) + &
                                            & k1 * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b1(n) + &
                                            & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b1(n)) + &
                                            & perpinv * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b2(n) + &
                                            & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b2(n)) + &
                                            & k3 * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b3(n) + &
                                            & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b3(n))
                    end select
                enddo
            endif
!
            if(poly_order.ge.3) then
!
                perpinv3 = perpinv2*perpinv
!
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    select case(i_precomp)
                        case(1)
                            coef_mat(n,4) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat3_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_2(:,n) + &
                                            & perpinv3 * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_3(:,n)) * z) + &
!
                                            & sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,n)) * b)
                    end select
                enddo
            endif
!
            if(poly_order.ge.4) then
!
                perpinv4 = perpinv2*perpinv2
!
                do n = 1,4
                    if(boole_faces_not(n)) cycle
                    select case(i_precomp)
                        case(1)
                            coef_mat(n,5) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat4_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat4_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat4_2(:,n) + &
                                            & perpinv3 * tetra_physics_poly4(ind_tetr)%anorm_in_amat4_3(:,n) + &
                                            & perpinv4 * tetra_physics_poly4(ind_tetr)%anorm_in_amat4_4(:,n)) * z) + &
!
                                            & sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat3_0(:,n) + &
                                            & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_1(:,n) + &
                                            & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_2(:,n) + &
                                            & perpinv3 * tetra_physics_poly4(ind_tetr)%anorm_in_amat3_3(:,n)) * b)
                    end select
                enddo
            endif
!
        end subroutine analytic_coeff_with_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_coeff_with_precomp_perpinv(poly_order,z,coef_mat)
!
            use tetra_physics_mod, only: tetra_physics
            use tetra_physics_poly_precomp_mod, only: tetra_physics_poly_perpinv,tetra_physics_poly4
!
            implicit none
!
            integer, intent(in)                 :: poly_order
            integer                             :: n
            double precision, dimension(4)      :: z
            double precision                    :: dist1
            double precision, dimension(4,poly_order+1)      :: coef_mat
!
                !Precompute quatities dependent on perpinv
!             if(.not.boole_precomp_analytic(ind_tetr)) call make_precomp_poly_perpinv(ind_tetr)
!
            dist1= -tetra_physics(ind_tetr)%dist_ref
!
            if(poly_order.ge.0) then
                coef_mat(:,1)=matmul(z(1:3),tetra_physics(ind_tetr)%anorm)
                coef_mat(1,1)=coef_mat(1,1)-dist1
            endif
!
            if(poly_order.ge.1) then
                do n = 1,4
                    coef_mat(n,2) = sum(tetra_physics_poly_perpinv(ind_tetr)%bcoef_pre_z(:,n) * z) + &
                                    & tetra_physics_poly_perpinv(ind_tetr)%bcoef_pre_k0(n) + &
                                    & k1 * tetra_physics_poly4(ind_tetr)%anorm_in_b1(n) + &
                                    & k3 * tetra_physics_poly4(ind_tetr)%anorm_in_b3(n)
                enddo
            endif
!
            if(poly_order.ge.2) then
                do n = 1,4
                    coef_mat(n,3) = sum(tetra_physics_poly_perpinv(ind_tetr)%acoef_pre_z(:,n)* z) + &
                                    & tetra_physics_poly_perpinv(ind_tetr)%acoef_pre_k0(n) + &
                                    & tetra_physics_poly_perpinv(ind_tetr)%acoef_pre_k1(n) * k1 + &
                                    & tetra_physics_poly_perpinv(ind_tetr)%acoef_pre_k3(n) * k3
                enddo
            endif
!
        end subroutine analytic_coeff_with_precomp_perpinv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Linear_Solver(a,b,dtau)
!
        !Find the root of the equation
        !f(tau) = a * tau + b
!
        implicit none
!
        double precision,intent(out)                :: dtau
        double precision,intent(in)                 :: a,b
!
        if(a.eq.0.d0) then
            dtau = huge(0.d0)
        else
            dtau = -b/a    
        endif        
!
    end subroutine Linear_Solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Quadratic_Solver(i_scaling,a,b,c,dtau)
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
!
        implicit none
!
        integer, intent(in) :: i_scaling
        double precision, intent(in)  :: a,b,c
        double precision,intent(out)                :: dtau
!
        select case(i_scaling)
            case(0)
                call Quadratic_Solver1(a,b,c,dtau)
            case DEFAULT
if(diag_pusher_tetry_poly) print *, 'New quadratic solver is called.'
                call Quadratic_Solver2(a,b,c,dtau)
        end select
                    
    end subroutine Quadratic_Solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Quadratic_Solver1(acoef,bcoef,ccoef,dtau)
!
        !Find the root of the equation
        !f(tau) = a/2 * tau**2 + b * tau + c
!        
        use tetra_physics_mod, only: tetra_physics
!         use orbit_timestep_gorilla_mod, only: counter_tetrahedron_passes
        use constants, only: eps
!
        implicit none
!
        double precision,intent(out)                :: dtau
        double precision                            :: discr, dummy
        double precision, intent(in)                :: acoef,bcoef,ccoef
!    
        !Initialize dtau        
        dtau = huge(0.d0)
!                   
        if(ccoef.gt.0.d0) then
            if(acoef.gt.0.d0) then
                if(bcoef.lt.0.d0) then
                    discr=bcoef**2-2.d0*acoef*ccoef
                    if(discr.gt.0.d0) then
                        dummy = (-bcoef +sqrt(discr))
                        if( abs(dummy).gt.eps ) then
                            dtau = 2.d0 * ccoef / dummy
                        else
                            dtau=(-sqrt(discr)-bcoef)/acoef !Numerically unstable
                        endif
                    elseif(discr.eq.0.d0) then
                        dtau=-bcoef/acoef
                    else !discr < 0
                        return
                    endif
                else ! b >= 0
                    return
                endif
            elseif(acoef.lt.0.d0) then
                discr=bcoef**2-2.d0*acoef*ccoef
                dummy = (-bcoef +sqrt(discr))
                if( abs(dummy).gt.eps ) then
                    dtau = 2.d0 * ccoef / dummy
                else
                    dtau=(-sqrt(discr)-bcoef)/acoef !Numerically unstable
                endif
            else !acoef = 0.d0
                if(bcoef.lt.0.d0) then
                    dtau = -ccoef/bcoef
                else
                    return
                endif
            endif
        elseif(ccoef.lt.0.d0) then
            if(acoef.lt.0.d0) then
                if(bcoef.gt.0.d0) then
                    discr=bcoef**2-2.d0*acoef*ccoef
                    if(discr.gt.0.d0) then
                        dtau=(sqrt(discr)-bcoef)/acoef
                    elseif(discr.eq.0.d0) then
                        dtau=-bcoef/acoef
                    else !discr < 0
                        return 
                    endif
                else ! b <= 0
                    return
                endif
            elseif(acoef.gt.0.d0) then
                discr=bcoef**2-2.d0*acoef*ccoef
                dtau=(sqrt(discr)-bcoef)/acoef
            else !acoef = 0.d0
                if(bcoef.gt.0.d0) then
                    dtau = -ccoef/bcoef
                else
                    return
                endif
            endif                
        else !ccoef = 0.d0
            if(((acoef.gt.0d0).and.(bcoef.lt.0d0)).or.((acoef.lt.0d0).and.(bcoef.gt.0d0))) then !If a and b have opposing signs
                dtau = -2.d0 * bcoef/acoef
            else
                return
            endif    
        endif
!
    end subroutine Quadratic_Solver1    
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       
    subroutine Quadratic_Solver2(a,b,c,min_tau_real)
!       WRAPPER FOR Polynomial234RootSolvers with only positive real roots
        !Find the root of the equation
        !f(tau) = a/2 * tau**2 + b * tau + c
!
        use Polynomial234RootSolvers, only: quadraticRoots
!        
        implicit none
!        
        double precision, intent(in)  :: a,b,c
        double precision              :: lambda,q1,q0
        integer                       :: nReal
        logical,dimension(2)          :: boole_positive_real
!         complex*16 , dimension(3)     :: tau
        double precision, intent(out) :: min_tau_real
        double precision              :: root(1:2,1:2)
!
        root = 0.d0
!         
        !Rescaling factor for coefficients
        lambda = b/c
!       
        q0 = 2.d0 * b**2 / (a * c)
        q1 = q0 
!         
        call quadraticRoots ( q1, q0, nReal, root)
        
        !Rescale root (Only if scaling parameter lambda is used)
        root = root/lambda
! 
        boole_positive_real = (abs(root(:,2)).eq.0.d0).and.(root(:,1).gt.0.d0)
        min_tau_real = minval(root(:,1),1,boole_positive_real)
!
    end subroutine Quadratic_Solver2
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!       
    subroutine Cubic_Solver(a,b,c,d,min_tau_real)
!       WRAPPER FOR Polynomial234RootSolvers with only positive real roots
!       f(tau) = a/6 * tau**3 + b/2 * tau**2 + c * tau + d
! 
        use Polynomial234RootSolvers, only: cubicRoots
!        
        implicit none
!        
        double precision, intent(in)  :: a,b,c,d
        double precision              :: lambda,c2,c1,c0
        integer                       :: nReal
        logical,dimension(3)          :: boole_positive_real
!         complex*16 , dimension(3)     :: tau
        double precision, intent(out) :: min_tau_real
        double precision              :: root(1:3,1:2)
!
        root = 0.d0
!         
        !Rescaling factor for coefficients
        lambda = b/(2.d0 * c)
!         
        c2 = 3.d0 * lambda * b / a
        c1 = 6.d0 * c * lambda**2 / a
        c0 = 6.d0 * d * lambda**3 / a
!         
        call cubicRoots (c2, c1, c0, nReal, root)
        
        !Rescale root (Only if scaling parameter lambda is used)
        root = root/lambda
! 
        boole_positive_real = (abs(root(:,2)).eq.0.d0).and.(root(:,1).gt.0.d0)
        min_tau_real = minval(root(:,1),1,boole_positive_real)
!
    end subroutine Cubic_Solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine Quartic_Solver(i_scaling,a,b,c,d,e,min_tau_real)
!       WRAPPER FOR QUARTIC SOLVER with only positive real roots
!       f(tau) = a/24 * tau**4 + b/6 * tau**3 + c/2 * tau**2 + d * tau + e
!
        use Polynomial234RootSolvers, only: quarticRoots
!
        implicit none
!        
        integer, intent(in)           :: i_scaling
        double precision, intent(in)  :: a,b,c,d,e
        logical,dimension(4)          :: boole_positive_real
!         complex*16 , dimension(4)     :: tau
        double precision, intent(out) :: min_tau_real
        integer                       :: nReal
        double precision              :: lambda,q3,q2,q1,q0
        double precision              :: root(1:4,1:2)
!
        !Scaling parameter lambda
        select case(i_scaling)
            case(0)
                lambda = sqrt(abs(b/(6.d0*d)))
            case(1)
                lambda = b/(3.d0*c)
            case(2)
                lambda = abs(b/(6.d0*e))**(1.d0/3.d0)
            case(3)
                lambda = c/(2.d0*d)
            case(4)
                lambda = sqrt(abs(c/(2.d0*e)))
            case(5)
                lambda = d/e
            case(6)
                lambda = abs(a/(24.d0*e))**(1.d0/4.d0)
        end select
!        
        q3 = 4.d0*b*lambda/a 
        q2 = 12.d0*c*lambda**2/a
        q1 = 24.d0*d*lambda**3/a
        q0 = 24.d0*e*lambda**4/a
!
        root = 0.d0
! 
        call  quarticRoots (q3, q2, q1, q0, nReal, root)
!         
        !Rescale root (Only if scaling parameter lambda is used)
        root = root/lambda
! 
        boole_positive_real = (abs(root(:,2)).eq.0.d0).and.(root(:,1).gt.0.d0)
        min_tau_real = minval(root(:,1),1,boole_positive_real)
!
    end subroutine Quartic_Solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration(poly_order,i_precomp,z,tau)
!
            implicit none
!
            integer, intent(in)                             :: poly_order,i_precomp
            double precision, intent(in)                    :: tau
            double precision, dimension(4),intent(inout)    :: z
!
            !Integrate trajectory analytically
            select case(i_precomp)
                case(0)
                    call analytic_integration_without_precomp(poly_order,z,tau)
                case(1,2)
                    call analytic_integration_with_precomp(poly_order,i_precomp,z,tau)
                case(3)
                    !call analytic_integration_with_precomp_perpinv(z,tau)
            end select
!
        end subroutine analytic_integration
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration_without_precomp(poly_order,z,tau)
!
            use poly_without_precomp_mod
!
            implicit none
!
            integer, intent(in)                          :: poly_order
            double precision, intent(in)                 :: tau
            double precision, dimension(4),intent(inout) :: z
            double precision                             :: tau2_half,tau3_sixth,tau4_twentyfourth
!
            !for each time integration performed, tau and z are saved for later computation of optional quantities
            !it is important to handle number_of_integration_steps appropriately when intrgration steps are discarded
            number_of_integration_steps = number_of_integration_steps + 1
            tau_steps_list(number_of_integration_steps) = tau
            intermediate_z0_list(:,number_of_integration_steps) = z
!
            if(poly_order.ge.1) then
                z = z + tau*(b+amat_in_z)
            endif
!
            if(poly_order.ge.2) then
                tau2_half = tau**2*0.5d0
                z = z + tau2_half*(amat_in_b + amat2_in_z)
            endif
!
            if(poly_order.ge.3) then
                tau3_sixth = tau**3/6.d0
                z = z + tau3_sixth*(amat2_in_b + amat3_in_z)
            endif

            if(poly_order.ge.4) then
                tau4_twentyfourth = tau**4/24.d0
                z = z + tau4_twentyfourth *(amat3_in_b + amat4_in_z)
            endif
!
        end subroutine analytic_integration_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine set_integration_coef_manually(poly_order,z0)
!
            use poly_without_precomp_mod, only: amat,amat2,amat3,amat4,b,&
                                                & amat_in_b,amat2_in_b,amat3_in_b,&
                                                & amat_in_z,amat2_in_z,amat3_in_z,amat4_in_z
!
            implicit none
!
            integer, intent(in)                          :: poly_order
            double precision, dimension(4),intent(in)    :: z0
!
            if(poly_order.ge.1) then
                amat_in_z = matmul(amat,z0)
            endif
            if(poly_order.ge.2) then
                amat2_in_z = matmul(amat2,z0)
                amat_in_b = matmul(amat,b)
            endif
            if(poly_order.ge.3) then
                amat3_in_z = matmul(amat3,z0)
                amat2_in_b = matmul(amat2,b)
            endif
            if(poly_order.ge.4) then
                amat4_in_z = matmul(amat4,z0)
                amat3_in_b = matmul(amat3,b)
            endif
        end subroutine set_integration_coef_manually
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialise_optional_quantities(optional_quantities)
!
            use gorilla_settings_mod, only: optional_quantities_type
!
            implicit none
!
            type(optional_quantities_type)    :: optional_quantities
!
            optional_quantities%t_hamiltonian = 0.d0
            optional_quantities%gyrophase = 0.d0
            optional_quantities%vpar_int = 0.d0
            optional_quantities%vpar2_int = 0.d0
!
        end subroutine initialise_optional_quantities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_optional_quantities(poly_order,z0,tau,optional_quantities)
!
        use gorilla_settings_mod, only: boole_time_Hamiltonian, boole_gyrophase, boole_vpar_int, boole_vpar2_int, &
                                        & optional_quantities_type
        use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(4), intent(in)      :: z0
        type(optional_quantities_type), intent(inout)   :: optional_quantities
!
        double precision                                :: delta_t_hamiltonian
        double precision                                :: tau2_half,tau3_sixth,tau4_twentyfourth
        double precision, dimension(3)                  :: x0
        double precision                                :: vpar0
        double precision, dimension(:,:), allocatable   :: x_coef,x_vpar_coef
        double precision, dimension(:), allocatable     :: vpar_coef,res_poly_coef
        double precision, dimension(4,poly_order+1)     :: dummy_coef_mat !exists because of analytic_coeff_without_precomp syntax
        logical, dimension(4)                           :: boole_faces_off = .false. !avoid root coefficients computation
!
        !recalculate polynomial coefficients (tensors) if more than one integration step was performed
        if (number_of_integration_steps .gt. 1) call set_integration_coef_manually(poly_order,z0)
!
            allocate(x_coef(3,poly_order+1))
            allocate(vpar_coef(poly_order+1))
!
            call z_series_coef(poly_order,z0,x_coef,vpar_coef)
!
        if(boole_time_hamiltonian.or.boole_vpar_int.or.boole_vpar2_int) then
!
            allocate(x_vpar_coef(3,poly_order+1))
!
            call poly_multiplication_coef(x_coef(1,:),vpar_coef(:),x_vpar_coef(1,:))
            call poly_multiplication_coef(x_coef(2,:),vpar_coef(:),x_vpar_coef(2,:))
            call poly_multiplication_coef(x_coef(3,:),vpar_coef(:),x_vpar_coef(3,:))
        endif
!
        !Optional computation of Hamiltonian time
        if(boole_time_hamiltonian) then
            delta_t_hamiltonian = hamiltonian_time(ind_tetr)%h1_in_curlA * tau + &
            & cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh * moment_integration(poly_order,tau,vpar_coef)+&
            & sum( hamiltonian_time(ind_tetr)%vec_mismatch_der * moment_integration(poly_order,tau,x_coef) ) + &
            & cm_over_e * sum(hamiltonian_time(ind_tetr)%vec_parcurr_der *  &
            & moment_integration(poly_order,tau,x_vpar_coef))
!
            !Multiply delta_t_hamiltonian with appropriate sign (We require that tau remains positive inside the algorithm)
            delta_t_hamiltonian = delta_t_hamiltonian * dble(sign_rhs)
!
!            call calc_t_hamiltonian(poly_order,z0,tau,delta_t_hamiltonian)
!
            optional_quantities%t_hamiltonian = optional_quantities%t_hamiltonian + delta_t_hamiltonian
!
            !Optional computation of gyrophase
            if(boole_gyrophase) then
                optional_quantities%gyrophase = optional_quantities%gyrophase - ( &
!
                    !Zeroth term
                    & 1.d0/cm_over_e * tetra_physics(ind_tetr)%bmod1 * delta_t_hamiltonian + &
!
                    !First term
                    & 1.d0/cm_over_e * sum( tetra_physics(ind_tetr)%gb * moment_integration(poly_order,tau,x_coef) ) *&
                    & hamiltonian_time(ind_tetr)%h1_in_curlA * dble(sign_rhs) + &
!
                    !Second term
                    & sum(moment_integration(poly_order,tau,x_vpar_coef) * tetra_physics(ind_tetr)%gb ) * &
                    & hamiltonian_time(ind_tetr)%h1_in_curlh * dble(sign_rhs)+ &
!
                    !Third term
                    & 1.d0/cm_over_e * sum( matmul( moment_integration(poly_order,tau,poly_multiplication(x_coef,x_coef)) , &
                    & hamiltonian_time(ind_tetr)%vec_mismatch_der ) * tetra_physics(ind_tetr)%gb) * dble(sign_rhs) + &
!
                    !Fourth term
                    & sum( matmul( moment_integration(poly_order,tau,poly_multiplication(x_coef,x_vpar_coef)), &
                    & hamiltonian_time(ind_tetr)%vec_parcurr_der ) *  tetra_physics(ind_tetr)%gb) * dble(sign_rhs) &
                & )
!

            endif ! gyrophase
!
        endif ! time_Hamiltonian
!
        if (boole_vpar_int) then 
            ! optional_quantities%vpar_int = optional_quantities%vpar_int + moment_integration(poly_order,tau,vpar_coef)* &
            !                                & dt_dtau_const
            optional_quantities%vpar_int = optional_quantities%vpar_int + &
!
                !First term
                & hamiltonian_time(ind_tetr)%h1_in_curlA * moment_integration(poly_order,tau,vpar_coef) + &
!
                !Second term
                & cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh * &
                & moment_integration(poly_order,tau,poly_multiplication(vpar_coef,vpar_coef))+ &
!
                !Third term
                & sum( hamiltonian_time(ind_tetr)%vec_mismatch_der * &
                & moment_integration(poly_order,tau,x_vpar_coef) ) + &
!
                !Fourth term
                & cm_over_e * sum(hamiltonian_time(ind_tetr)%vec_parcurr_der *  &
                & moment_integration(poly_order,tau,poly_multiplication(vpar_coef,x_vpar_coef)))
!
!print*, 'vpar:',vpar_coef, tau, moment_integration(poly_order,tau,vpar_coef)
!
        endif
!
        if (boole_vpar2_int) then
            ! optional_quantities%vpar2_int = optional_quantities%vpar2_int  + dt_dtau_const*&
            !                                 & moment_integration(poly_order,tau,poly_multiplication(vpar_coef,vpar_coef))
            optional_quantities%vpar2_int = optional_quantities%vpar2_int + &
!
                !First term
                & hamiltonian_time(ind_tetr)%h1_in_curlA * &
                & moment_integration(poly_order,tau,poly_multiplication(vpar_coef,vpar_coef)) + &
!
                !Second term
                & cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh * &
                & moment_integration(poly_order,tau,poly_multiplication(vpar_coef,poly_multiplication(vpar_coef,vpar_coef)))+ &
!
                !Third term
                & sum( hamiltonian_time(ind_tetr)%vec_mismatch_der * &
                & moment_integration(poly_order,tau,poly_multiplication(poly_multiplication(vpar_coef,vpar_coef),x_coef)) ) + &
!
                !Fourth term
                & cm_over_e * sum(hamiltonian_time(ind_tetr)%vec_parcurr_der *  &
                & moment_integration(poly_order,tau,poly_multiplication(poly_multiplication(vpar_coef,vpar_coef),x_vpar_coef)))
!

!print*, 'vpar2', poly_multiplication(vpar_coef,vpar_coef), tau, & 
!moment_integration(poly_order,tau,poly_multiplication(vpar_coef,vpar_coef))
!
        endif
!
        deallocate(vpar_coef)
        if(boole_time_hamiltonian.or.boole_vpar_int.or.boole_vpar2_int)  deallocate(x_coef,x_vpar_coef)
!
    end subroutine calc_optional_quantities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_t_hamiltonian(poly_order,z0,tau,t_hamiltonian)
!
        use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(4), intent(in)      :: z0
        double precision, intent(inout)                 :: t_hamiltonian
!
        double precision                                :: delta_t_hamiltonian
        double precision, dimension(:,:), allocatable   :: x_coef,x_vpar_coef
        double precision, dimension(:), allocatable     :: vpar_coef
!
        !recalculate polynomial coefficients (tensors) if more than one integration step was performed
        if (number_of_integration_steps .gt. 1) call set_integration_coef_manually(poly_order,z0)

        allocate(x_coef(3,poly_order+1))
        allocate(vpar_coef(poly_order+1))
        allocate(x_vpar_coef(3,poly_order+1))
!
        call z_series_coef(poly_order,z0,x_coef,vpar_coef)
!
        call poly_multiplication_coef(x_coef(1,:),vpar_coef(:),x_vpar_coef(1,:))
        call poly_multiplication_coef(x_coef(2,:),vpar_coef(:),x_vpar_coef(2,:))
        call poly_multiplication_coef(x_coef(3,:),vpar_coef(:),x_vpar_coef(3,:))
!
        delta_t_hamiltonian = hamiltonian_time(ind_tetr)%h1_in_curlA * tau + &
        & cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh * moment_integration(poly_order,tau,vpar_coef)+&
        & sum( hamiltonian_time(ind_tetr)%vec_mismatch_der * moment_integration(poly_order,tau,x_coef) ) + &
        & cm_over_e * sum(hamiltonian_time(ind_tetr)%vec_parcurr_der *  &
        & moment_integration(poly_order,tau,x_vpar_coef))
!
        !Multiply delta_t_hamiltonian with appropriate sign (We require that tau remains positive inside the algorithm)
        delta_t_hamiltonian = delta_t_hamiltonian * dble(sign_rhs)
!
        t_hamiltonian = t_hamiltonian + delta_t_hamiltonian
!
        deallocate(x_coef,vpar_coef,x_vpar_coef)
!
    end subroutine calc_t_hamiltonian
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_t_hamiltonian_by_order(poly_order,z0,tau,delta_t_hamiltonian)
!
        use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, dimension(4), intent(in)      :: z0
        double precision, intent(in)                    :: tau
        double precision, intent(out)                   :: delta_t_hamiltonian
        double precision, dimension(:,:), allocatable   :: x_coef,x_vpar_coef
        double precision, dimension(:), allocatable     :: vpar_coef
        double precision                                :: a_coef, b_coef, c_coef, d_coef, e_coef

!
        !recalculate polynomial coefficients (tensors) if more than one integration step was performed
        !if (number_of_integration_steps .gt. 1) call set_integration_coef_manually(poly_order,z0)
!
        allocate(x_coef(3,poly_order+1))
        allocate(vpar_coef(poly_order+1))
        allocate(x_vpar_coef(3,poly_order+1))
!
        call z_series_coef(poly_order,z0,x_coef,vpar_coef)
!
        call poly_multiplication_coef(x_coef(1,:),vpar_coef(:),x_vpar_coef(1,:))
        call poly_multiplication_coef(x_coef(2,:),vpar_coef(:),x_vpar_coef(2,:))
        call poly_multiplication_coef(x_coef(3,:),vpar_coef(:),x_vpar_coef(3,:))
!
        !Input for Root Solver
        !f(tau) = a/24 * tau**4 + b/6 * tau**3 + c/2 * tau**2 + d * tau + e
!
        a_coef = (1.d0/5.d0) * ( &
        & vpar_coef(5) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum ( x_coef(:,5) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum (x_vpar_coef(:,5) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        b_coef = (1.d0/4.d0) * ( &
        & vpar_coef(4) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum ( x_coef(:,4) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum (x_vpar_coef(:,4) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        c_coef = (1.d0/3.d0) * ( &
        & vpar_coef(3) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum ( x_coef(:,3) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum (x_vpar_coef(:,3) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        d_coef = 0.5d0 *(  &
        & vpar_coef(2) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum ( x_coef(:,2) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum (x_vpar_coef(:,2) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        e_coef = ( &
        & hamiltonian_time(ind_tetr)%h1_in_curlA + &
        & vpar_coef(1) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum( x_coef(:,1) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum( x_vpar_coef(:,1) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        delta_t_hamiltonian = &
!        & tau**5 * a_coef + &
        & tau**4 * b_coef + &
        & tau**3 * c_coef + &
        & tau**2 * d_coef + &
        & tau * e_coef
!
        deallocate(x_coef,vpar_coef,x_vpar_coef)
!
    end subroutine calc_t_hamiltonian_by_order

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine get_t_hamiltonian_root(poly_order,z0,t_remain,tau_t_hamiltonian_root)
!
        use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
!
        implicit none
!
        integer, intent(in)                                     :: poly_order
        double precision, dimension(4), intent(in)              :: z0
        double precision, intent(in)                            :: t_remain
        double precision, intent(out)                           :: tau_t_hamiltonian_root
        double precision                                        :: a_coef, b_coef, c_coef, d_coef, e_coef
        double precision, dimension(:,:), allocatable           :: x_coef,x_vpar_coef
        double precision, dimension(:), allocatable             :: vpar_coef
!
        !recalculate polynomial coefficients (tensors) if more than one integration step was performed
        !if (number_of_integration_steps .gt. 1) call set_integration_coef_manually(poly_order,z0)
!
        allocate(x_coef(3,poly_order+1))
        allocate(vpar_coef(poly_order+1))
        allocate(x_vpar_coef(3,poly_order+1))
!
        call z_series_coef(poly_order,z0,x_coef,vpar_coef)
!
        call poly_multiplication_coef(x_coef(1,:),vpar_coef(:),x_vpar_coef(1,:))
        call poly_multiplication_coef(x_coef(2,:),vpar_coef(:),x_vpar_coef(2,:))
        call poly_multiplication_coef(x_coef(3,:),vpar_coef(:),x_vpar_coef(3,:))
!
        !Coefficient for Root Solver
        !t_hamiltonian(tau) = a * tau**5 + b * tau**4 + c * tau**3 + d * tau**2 + e * tau - t_remain == 0
!
        !We neglect the 5th order term in the root solving operation in order to be analytical.
        !This is only relevant in the "final" tetrahedron of a time step.
!
        !t_hamiltonian(tau) = b * tau**4 + c * tau**3 + d * tau**2 + e * tau - t_remain == 0
!
        e_coef = ( &
        & hamiltonian_time(ind_tetr)%h1_in_curlA + &
        & vpar_coef(1) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
        & sum( x_coef(:,1) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
        & cm_over_e * sum( x_vpar_coef(:,1) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
        & )
!
        if(poly_order.ge.1) then
            d_coef = 0.5d0 *(  &
            & vpar_coef(2) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
            & sum ( x_coef(:,2) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
            & cm_over_e * sum (x_vpar_coef(:,2) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
            & )
        endif
!
        if(poly_order.ge.2) then
            c_coef = (1.d0/3.d0) * ( &
            & vpar_coef(3) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
            & sum ( x_coef(:,3) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
            & cm_over_e * sum (x_vpar_coef(:,3) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
            & )
        endif
!
        if(poly_order.ge.3) then
            b_coef = (1.d0/4.d0) * ( &
            & vpar_coef(4) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
            & sum ( x_coef(:,4) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
            & cm_over_e * sum (x_vpar_coef(:,4) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
            & )
        endif
!
!        a_coef = (1.d0/5.d0) * ( &
!        & vpar_coef(5) * cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh + &
!        & sum ( x_coef(:,5) * hamiltonian_time(ind_tetr)%vec_mismatch_der(:) ) + &
!        & cm_over_e * sum (x_vpar_coef(:,5) * hamiltonian_time(ind_tetr)%vec_parcurr_der(:) ) &
!        & )
!
        deallocate(x_coef,vpar_coef,x_vpar_coef)
!
        !Manipulate coefficients for root solvers
        !f(tau) = b/24 * tau**4 + c/6 * tau**3 + d/2 * tau**2 + e * tau - t_remain == 0
!
        !t_remain = t_remain
        !e_coef = e_coef
        if(poly_order.ge.1) d_coef = 2.d0 * d_coef
        if(poly_order.ge.2) c_coef = 6.d0 * c_coef
        if(poly_order.ge.3) b_coef = 24.d0 * b_coef
!
        !Multiply delta_t_hamiltonian with appropriate sign (We require that tau remains positive inside the algorithm)
        e_coef = e_coef * sign_rhs
        if(poly_order.ge.1) d_coef = d_coef * sign_rhs
        if(poly_order.ge.2) c_coef = c_coef * sign_rhs
        if(poly_order.ge.3) b_coef = b_coef * sign_rhs
!
        !Find root
        select case(poly_order)
            case(1)
                call Quadratic_Solver2(d_coef,e_coef,-t_remain,tau_t_hamiltonian_root)
            case(2)
                call Cubic_Solver(c_coef,d_coef,e_coef,-t_remain,tau_t_hamiltonian_root)
            case(3,4)
                call Quartic_Solver(0,b_coef,c_coef,d_coef,e_coef,-t_remain,tau_t_hamiltonian_root)
        end select
!
    end subroutine get_t_hamiltonian_root
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine z_series_coef(poly_order,z0,x_coef,vpar_coef)
!
        !This function shall only be called within the subroutine "calc_optional_quantities"
        !Values in modules are used that need to be precomputed/set in that subroutine.
!
        use poly_without_precomp_mod, only: b, amat_in_z,amat2_in_z,amat3_in_z,amat_in_b,amat2_in_b, &
                                            & amat4_in_z,amat3_in_b
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, dimension(4), intent(in)      :: z0
!
        double precision, dimension(:,:),intent(out)    :: x_coef
        double precision, dimension(:),intent(out)      :: vpar_coef
!
        x_coef(:,1) = z0(1:3)
        vpar_coef(1) = z0(4)
!
        if(poly_order.ge.1) then
            x_coef(:,2) = b(1:3) + amat_in_z(1:3)
            vpar_coef(2) = b(4) + amat_in_z(4)
        endif
!
        if(poly_order.ge.2) then
            x_coef(:,3) = 0.5d0 * (amat_in_b(1:3) + amat2_in_z(1:3))
            vpar_coef(3) = 0.5d0 * (amat_in_b(4) + amat2_in_z(4))
        endif
!
        if(poly_order.ge.3) then
            x_coef(:,4) = 1.d0 / 6.d0 * (amat2_in_b(1:3) + amat3_in_z(1:3))
            vpar_coef(4) = 1.d0 / 6.d0 * (amat2_in_b(4) + amat3_in_z(4))
        endif
!
        if(poly_order.ge.4) then
            x_coef(:,5) = 1.d0 / 24.d0 * (amat3_in_b(1:3) + amat4_in_z(1:3))
            vpar_coef(5) = 1.d0 / 24.d0 * (amat3_in_b(4) + amat4_in_z(4))
        endif
!
    end subroutine z_series_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    pure subroutine poly_multiplication_coef(poly_coef_1,poly_coef_2,res_poly_coef)
!
        implicit none
!
        double precision, dimension(:),intent(in) :: poly_coef_1,poly_coef_2
!
        double precision, dimension(:),intent(out) :: res_poly_coef
!
        integer :: max_order,i,j,k,sz_1,sz_2,cur_order
        integer,dimension(:),allocatable :: poly_order_vec_1,poly_order_vec_2
!
        sz_1 = size(poly_coef_1)
        sz_2 = size(poly_coef_2)
!
        !Maximum order of the two polynomials - No higher order than that will be computed
        max_order = maxval([sz_1,sz_2])
!
        allocate(poly_order_vec_1(sz_1),poly_order_vec_2(sz_2))
!
        poly_order_vec_1 = [(i, i = 0,sz_1-1,1)]
        poly_order_vec_2 = [(i, i = 0,sz_2-1,1)]
!
        !Initialize coefficient result vector
        res_poly_coef = 0.d0
!
        !Loop over coefficients of first polynomial
        do j = 1,sz_1
!
            !Loop over coefficients of first polynomial
            do k = 1,sz_2
!
                !Current order of term multiplication
                cur_order = (poly_order_vec_1(j) + poly_order_vec_2(k))
!
                !Check order of the multiplied term, and exit if max order is exceeded
                if( cur_order.gt.(max_order-1) ) exit
!
                res_poly_coef(cur_order+1) = res_poly_coef(cur_order+1) + &
                    & poly_coef_1(j) * poly_coef_2(k)

            enddo
!
        enddo
!
        deallocate(poly_order_vec_1,poly_order_vec_2)
        
    end subroutine poly_multiplication_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration_with_precomp(poly_order,i_precomp,z,tau)
!
            use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4
            use poly_without_precomp_mod, only: b,perpinv3,perpinv4
!
            implicit none
!
            integer, intent(in)                 :: poly_order,i_precomp
            double precision, intent(in) :: tau
            double precision, dimension(4),intent(inout) :: z
            double precision ::tau2_half,tau3_sixth,tau4_twentyfourth
            double precision,dimension(4,4) :: operator_b,operator_z_init
            double precision,dimension(4) :: operator_b_in_b
!
            if(poly_order.ge.2) then
                tau2_half = tau**2*0.5d0
            endif
!
            if(poly_order.ge.3) then
                tau3_sixth = tau**3/6.d0
            endif
!
            if(poly_order.ge.4) then
                tau4_twentyfourth = tau**4/24.d0
            endif
!
            select case(poly_order)
                case(2)
!
                    operator_z_init = tau*( perpinv*tetra_physics_poly4(ind_tetr)%amat1_1(:,:) + &
                                    & tetra_physics_poly4(ind_tetr)%amat1_0(:,:) ) + &
                                    & tau2_half* (perpinv2*tetra_physics_poly4(ind_tetr)%amat2_2(:,:) + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat2_1(:,:) + &
                                    & tetra_physics_poly4(ind_tetr)%amat2_0(:,:) )
!

                    select case(i_precomp)
                        case(1)
                            operator_b = tau*unity_matrix4 + tau2_half*(perpinv*tetra_physics_poly4(ind_tetr)%amat1_1(:,:) + &
                                        & tetra_physics_poly4(ind_tetr)%amat1_0(:,:))
!
                            !Integration of trajectory (analytic 2nd order)
                            z = z + matmul(operator_b,b) + matmul(operator_z_init,z)
                    case(2)
                        operator_b_in_b = tau* (tetra_physics_poly4(ind_tetr)%b0 + k1*tetra_physics_poly4(ind_tetr)%b1 &
                                        & + perpinv*tetra_physics_poly4(ind_tetr)%b2 + k3*tetra_physics_poly4(ind_tetr)%b3 )&
                                        & + tau2_half*( (tetra_physics_poly4(ind_tetr)%amat1_0_in_b0 &
                                        & + perpinv*tetra_physics_poly4(ind_tetr)%amat1_1_in_b0) &
                                        & + k1*(tetra_physics_poly4(ind_tetr)%amat1_0_in_b1 &
                                        & + perpinv*tetra_physics_poly4(ind_tetr)%amat1_1_in_b1) &
                                        & + perpinv*(tetra_physics_poly4(ind_tetr)%amat1_0_in_b2 &
                                        & + perpinv*tetra_physics_poly4(ind_tetr)%amat1_1_in_b2) &
                                        & + k3*(tetra_physics_poly4(ind_tetr)%amat1_0_in_b3 &
                                        & + perpinv*tetra_physics_poly4(ind_tetr)%amat1_1_in_b3) )
!
                        !Integration of trajectory (analytic 2nd order)
                        z = z + operator_b_in_b + matmul(operator_z_init,z)
                    end select
!
                case(3)

                    operator_z_init = tau*( perpinv*tetra_physics_poly4(ind_tetr)%amat1_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat1_0 ) + &
                                    & tau2_half* (perpinv2*tetra_physics_poly4(ind_tetr)%amat2_2 + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat2_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat2_0 ) + &
                                    & tau3_sixth* (perpinv3*tetra_physics_poly4(ind_tetr)%amat3_3 + &
                                    & perpinv2*tetra_physics_poly4(ind_tetr)%amat3_2 + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat3_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat3_0 )
!
                    select case(i_precomp)
                        case(1)
                            operator_b = tau*unity_matrix4 + tau2_half*(perpinv*tetra_physics_poly4(ind_tetr)%amat1_1 + &
                                        & tetra_physics_poly4(ind_tetr)%amat1_0) + &
                                        & tau3_sixth*(perpinv2*tetra_physics_poly4(ind_tetr)%amat2_2 + &
                                        & perpinv*tetra_physics_poly4(ind_tetr)%amat2_1 + &
                                        & tetra_physics_poly4(ind_tetr)%amat2_0)
!
                            !Integration of trajectory (analytic 3rd order)
                            z = z + matmul(operator_b,b) + matmul(operator_z_init,z)
                    end select
!
                case(4)

                    operator_z_init = tau*( perpinv*tetra_physics_poly4(ind_tetr)%amat1_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat1_0 ) + &
                                    & tau2_half* (perpinv2*tetra_physics_poly4(ind_tetr)%amat2_2 + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat2_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat2_0 ) + &
                                    & tau3_sixth* (perpinv3*tetra_physics_poly4(ind_tetr)%amat3_3 + &
                                    & perpinv2*tetra_physics_poly4(ind_tetr)%amat3_2 + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat3_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat3_0) + &
                                    & tau4_twentyfourth*(perpinv4*tetra_physics_poly4(ind_tetr)%amat4_4 + &
                                    & perpinv3*tetra_physics_poly4(ind_tetr)%amat4_3 + &
                                    & perpinv2*tetra_physics_poly4(ind_tetr)%amat4_2 + &
                                    & perpinv*tetra_physics_poly4(ind_tetr)%amat4_1 + &
                                    & tetra_physics_poly4(ind_tetr)%amat4_0 )
!
                    select case(i_precomp)
                        case(1)
                            operator_b = tau*unity_matrix4 + tau2_half*(perpinv*tetra_physics_poly4(ind_tetr)%amat1_1 + &
                                        & tetra_physics_poly4(ind_tetr)%amat1_0) + &
                                        & tau3_sixth*(perpinv2*tetra_physics_poly4(ind_tetr)%amat2_2 + &
                                        & perpinv*tetra_physics_poly4(ind_tetr)%amat2_1 + &
                                        & tetra_physics_poly4(ind_tetr)%amat2_0) + &
                                        & tau4_twentyfourth*(perpinv3*tetra_physics_poly4(ind_tetr)%amat3_3 + &
                                        & perpinv2*tetra_physics_poly4(ind_tetr)%amat3_2 + &
                                        & perpinv*tetra_physics_poly4(ind_tetr)%amat3_1 + &
                                        & tetra_physics_poly4(ind_tetr)%amat3_0)
!
                            !Integration of trajectory (analytic 4th order)
                            z = z + matmul(operator_b,b) + matmul(operator_z_init,z)
                    end select
            end select
!
        end subroutine analytic_integration_with_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration_with_precomp_perpinv(z,tau)
!
            use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4,tetra_physics_poly_perpinv
!
            implicit none
!
            double precision, intent(in) :: tau
            double precision, dimension(4),intent(inout) :: z
            double precision ::tau2_half
            double precision,dimension(4,4) :: unity_matrix4
            double precision,dimension(4,4) :: operator_b,operator_z_init
            double precision,dimension(4) :: operator_b_in_b
!
!             unity_matrix4 = reshape( [ 1.d0, 0.d0, 0.d0, 0.d0, &
!                                       0.d0, 1.d0, 0.d0, 0.d0, &
!                                       0.d0, 0.d0, 1.d0, 0.d0, &
!                                       0.d0, 0.d0, 0.d0, 1.d0 ], shape(unity_matrix4))
!
            tau2_half = tau**2*0.5d0
!
!             operator_b = tau*unity_matrix4 + tau2_half*perpinv*tetra_physics_poly4(ind_tetr)%amat1_1(:,:) + &
!                         & tau2_half*tetra_physics_poly4(ind_tetr)%amat1_0(:,:)
!
            operator_z_init = tau*tetra_physics_poly_perpinv(ind_tetr)%opz_pre_tau + &
                            & tau2_half*tetra_physics_poly_perpinv(ind_tetr)%opz_pre_tau2
!
            operator_b_in_b = tau* (tetra_physics_poly4(ind_tetr)%b0 + k1*tetra_physics_poly4(ind_tetr)%b1 &
                            & + perpinv*tetra_physics_poly4(ind_tetr)%b2 + k3*tetra_physics_poly4(ind_tetr)%b3 )&
                            & + tau2_half*(tetra_physics_poly_perpinv(ind_tetr)%opb_pre_k0 + &
                            & tetra_physics_poly_perpinv(ind_tetr)%opb_pre_k1*k1 + &
                            & tetra_physics_poly_perpinv(ind_tetr)%opb_pre_k3*k3 )
!
            !Integration of trajectory (analytic 2nd order)
            z = z + operator_b_in_b + matmul(operator_z_init,z)
!
        end subroutine analytic_integration_with_precomp_perpinv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    function normal_distance_func(z123,iface)
!
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!
        integer,intent(in)                          :: iface
        double precision, dimension(3),intent(in)   :: z123
        double precision                            :: normal_distance_func,dist1
!
        dist1= -tetra_physics(ind_tetr)%dist_ref
!
        normal_distance_func=sum(z123*tetra_physics(ind_tetr)%anorm(:,iface))
        if(iface.eq.1) normal_distance_func=normal_distance_func-dist1   !Correction for the one plane that is not lying in the first vertex
!        
    end function normal_distance_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function normal_velocity_func(z,iface)
!        
        use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4,tetra_physics_poly1
        use tetra_physics_mod, only: tetra_physics, cm_over_e
        use poly_without_precomp_mod, only: b
        use gorilla_settings_mod, only: i_precomp, boole_strong_electric_field
        use constants, only: clight
!
        implicit none
!        
        integer, intent(in)                         :: iface
        double precision, dimension(4), intent(in)  :: z
        double precision                            :: normal_velocity_func
!
        if(i_precomp.eq.0) then
            normal_velocity_func = (sum((-clight*tetra_physics_poly1(ind_tetr)%anorm_in_betmat(:,iface) + &
                        & perpinv*cm_over_e*tetra_physics_poly1(ind_tetr)%anorm_in_alpmat(:,iface)) * z(1:3)) + &
                        & tetra_physics_poly1(ind_tetr)%anorm_in_betvec(iface) * z(4)) * dble(sign_rhs) + &
                        & sum(tetra_physics(ind_tetr)%anorm(:,iface) * b(1:3))
            if(boole_strong_electric_field) then
                normal_velocity_func = normal_velocity_func + & 
                                    &   (sum(- 0.5d0*cm_over_e*tetra_physics_poly1(ind_tetr)%anorm_in_gammat(:,iface) * z(1:3)) &
                                    &   + cm_over_e*tetra_physics_poly1(ind_tetr)%anorm_in_gamvec(iface) * z(4)) * dble(sign_rhs) 
            endif
        else
            normal_velocity_func = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,iface) + &
                                    & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,iface)) * z) * dble(sign_rhs) + &
                                    & sum(tetra_physics(ind_tetr)%anorm(:,iface) * b(1:3))
        endif
!        
    end function normal_velocity_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine physical_estimate_tau(tau_est)
!
        use tetra_physics_mod, only: tetra_physics,cm_over_e
        use tetra_grid_settings_mod, only: grid_size
        use constants, only : clight
        use gorilla_settings_mod, only: boole_strong_electric_field
!
        implicit none
!
        double precision                           :: vperp2,tetra_dist_ref,vd_ExB,tau_est
        logical                                    :: boole_vd_ExB,boole_vperp
        double precision, parameter                :: eps_modulation = 0.1d0
!                      
        !Tetrahdron reference distance
        tetra_dist_ref = abs(tetra_physics(ind_tetr)%tetra_dist_ref)
!        
        vperp2 = -2.d0*perpinv*bmod0
!
        !ExB drift velocity
        !Only in case of strong electric field mode, ExB drift directly calculated during tetra_physics-setup
        if (boole_strong_electric_field) then
            vd_ExB = tetra_physics(ind_tetr)%v_E_mod_average
        else
            vd_ExB = abs(clight/bmod0*tetra_physics(ind_tetr)%Er_mod)
        endif
        boole_vd_ExB = .true.
        if(vd_ExB.eq.0.d0) boole_vd_ExB = .false.
!
        boole_vperp = .true.
        if(vperp2.eq.0.d0) boole_vperp = .false.
!
        !Reference time to pass the tetrahedron
        !Actually dt_ref (Save on additional variable dt_ref):
        if(boole_vd_ExB.and.boole_vperp) then
            tau_est = minval([abs(tetra_dist_ref/z_init(4)), &
            & sqrt(tetra_dist_ref*vmod0*tetra_physics(ind_tetr)%R1 / &
            & (vperp2*grid_size(2)*eps_modulation) ), &
            & tetra_dist_ref/vd_ExB],1)         
        elseif(boole_vd_ExB.and.(.not.boole_vperp)) then
            tau_est = minval([abs(tetra_dist_ref/z_init(4)), &
            & tetra_dist_ref/vd_ExB],1)        
        elseif((.not.boole_vd_ExB).and.boole_vperp) then
            tau_est = minval([abs(tetra_dist_ref/z_init(4)), &
            & sqrt(tetra_dist_ref*vmod0*tetra_physics(ind_tetr)%R1 / &
            & (vperp2*grid_size(2)*eps_modulation) )],1)          
        else
            tau_est = abs(tetra_dist_ref/z_init(4))            
        endif
!
        !Transform t to tau            
        tau_est = abs(tau_est/dt_dtau_const)
!
    end subroutine physical_estimate_tau
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
        use gorilla_settings_mod, only: boole_adaptive_time_steps
!
        implicit none
!
        integer, intent(in)                                 :: poly_order,i_precomp
        integer, intent(out)                                :: iface_new
        logical, intent(out)                                :: boole_trouble_shooting
        integer                                             :: j,k,i,poly_order_new
        integer                                             :: i_scaling
        double precision, dimension(4),intent(inout)        :: z
        double precision, intent(inout)                     :: tau
        double precision                                    :: tau_max,tau_max_est
        logical                                             :: boole_analytical_approx,boole_face_correct
        logical, dimension(4)                               :: boole_faces
!
        !Initialize boole_faces and i_scaling
        boole_faces = .true.
        i_scaling = 0
        boole_trouble_shooting = .true.
        iface_new = iface_init    !<-- Jonatan, Georg 07.07.2022: we think the face should be reseted because of analytic_approx belo
        z = z_init               !<-- Jonatan, Georg 07.07.2022. we think z should be reseted to allow for an appropriate tau_max
!
        !Analytical calculation of orbit parameter for safety boundary
        call analytic_approx(2,i_precomp,boole_faces, &
                            & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
        !Define boundary for tau_max
        tau_max = tau*eps_tau    
if(diag_pusher_tetry_poly) print *, 'quadratic tau',tau
!       
        !Initialization while loop
        boole_face_correct = .false.

        i = 0
        do while((.not.boole_face_correct).and.(poly_order.eq.4))
            i = i+1

            !Repeat root and orbit computation with different scaling for quartic solver
            i_scaling = i
            
            boole_faces = .true.    !Compute root for all 4 faces
            iface_new = iface_init
            z = z_init
            number_of_integration_steps = 0 !Resetting the number of integration steps because we start from scratch #1
!
            !Analytical calculation of orbit parameter to pass tetrahdron
            call analytic_approx(poly_order,i_precomp,boole_faces, &
                                & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
            !Analytical result does not exist.
            if(.not.boole_analytical_approx) then
                print *, 'Error: No analytical solution exists.'
                boole_trouble_shooting = .false.
                return                
            endif
!
            !Integrate trajectory analytically
            call analytic_integration(poly_order,i_precomp,z,tau)
!
            if (boole_adaptive_time_steps) then
                !For this section the starting face is the initial face
                !The i_scaling gets changed repeatedly
                call overhead_adaptive_time_steps(poly_order, i_scaling, .false., .true., &
                                                    &  iface_new, tau, z, boole_face_correct)
            endif !adaptive steps scheme
!
            !Control section
            boole_face_correct = .true.
!        
            call check_three_planes(z,iface_new,boole_face_correct)
            call check_face_convergence(z,iface_new,boole_face_correct)
            call check_velocity(z,iface_new,boole_face_correct,poly_order)    
!
            !Tau outside safety margin
            if(tau.gt.tau_max) then
                !Compare with physically estimated time
                call physical_estimate_tau(tau_max_est)
                tau_max_est = tau_max_est*eps_tau
!                
                if(tau.gt.tau_max_est) then
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: Safety margin'
                endif
            endif
!            
            if(i.eq.6) exit
        enddo
!
if(diag_pusher_tetry_poly) print *, 'scaling', i
!
        !3rd order polynomial pushing: (1) For Quadratic and Cubic solver, (2) if quartic solver fails
        if(.not.boole_face_correct) then
if(diag_pusher_tetry_poly) print *, 'Trouble shooting routine failed: Use 3rd order polynomial'
            boole_faces = .true.    !Compute root for all 4 faces
            iface_new = iface_init
            z = z_init
            number_of_integration_steps = 0 !Resetting the number of integration steps because we start from scratch #2
!
            select case(poly_order)
                case(2)
                    poly_order_new = 2
                    i_scaling = 1
                case(3)
                    poly_order_new = 3
                    i_scaling = 1
                case(4)
                    !Decrease poly_order
                    poly_order_new = 3
                    i_scaling = 0
            end select
!
            !Analytical calculation of orbit parameter to pass tetrahdron
            call analytic_approx(poly_order_new,i_precomp,boole_faces, &
                                & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
            !Analytical result does not exist.
            if(.not.boole_analytical_approx) then
                print *, 'Error: No analytical solution exists.'
                boole_trouble_shooting = .false.
                return                
            endif
!
            !Integrate trajectory analytically
            call analytic_integration(poly_order,i_precomp,z,tau)
!           
            !Control section
            boole_face_correct = .true.
!
            if (boole_adaptive_time_steps) then
                !For this section the starting face is the initial face
                !i_scaling and the poly_order were changed 
                call overhead_adaptive_time_steps(poly_order_new, i_scaling, .false., .true., &
                                                    &  iface_new, tau, z, boole_face_correct)
            endif !adaptive steps scheme    
!    
            !consistency checks
            call check_three_planes(z,iface_new,boole_face_correct)
            call check_face_convergence(z,iface_new,boole_face_correct)
            call check_velocity(z,iface_new,boole_face_correct,poly_order)            
!        
            !Tau outside safety margin
            if(tau.gt.tau_max) then
                !Compare with physically estimated time
                call physical_estimate_tau(tau_max_est)
                tau_max_est = tau_max_est*eps_tau
!                
                if(tau.gt.tau_max_est) then
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: Safety margin'
                endif
            endif
!            
            if(.not.boole_face_correct) then
                print *, 'Trouble shooting routine failed: Also 3rd order polynomial pushing failed'
                boole_trouble_shooting = .false.
                return                
            endif
!            
        endif
    
    end subroutine trouble_shooting_polynomial_solver
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module pusher_tetra_poly_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
module par_adiab_inv_poly_mod
!
    implicit none
!
    private
!
    integer,public,protected    :: counter_banana_mappings = 0
    double precision            :: par_adiab_inv
    integer                     :: nskip
!
    public :: par_adiab_inv_tetra_poly
!
    !$OMP THREADPRIVATE(counter_banana_mappings,par_adiab_inv)
!
    contains
!    
    subroutine par_adiab_inv_tetra_poly(poly_order,t_pass,vpar_in,vpar_end,file_id_vpar_0,file_id_J_par,n_skip_vpar_0, &
                                        & boole_J_par,boole_poincare_vpar_0,boole_e_tot,file_id_e_tot)
!
        !Attention: This subroutine CAN ONLY BE CALLED directly after pusher, when modules still contain
        !           the quantities from last pushing
!        
        use poly_without_precomp_mod, only: amat,b
        use supporting_functions_mod, only: energy_tot_func
        use pusher_tetra_poly_mod, only: dt_dtau_const,z_init,ind_tetr, perpinv, &
                                                number_of_integration_steps, intermediate_z0_list, tau_steps_list, &
                                                & set_integration_coef_manually
        use tetra_physics_mod, only: particle_mass,tetra_physics
        use gorilla_diag_mod, only: diag_pusher_tetry_poly_adaptive
!  
        implicit none
!
        integer, intent(in)                 :: poly_order,file_id_vpar_0,file_id_J_par,file_id_e_tot,n_skip_vpar_0
        double precision, intent(in)        :: t_pass,vpar_in,vpar_end
        logical, intent(in)                 :: boole_J_par,boole_poincare_vpar_0,boole_e_tot
        double precision                    :: a44,b4,tau_part1
        double precision, dimension(4)      :: z
        double precision, dimension(3)      :: x
        integer                             :: i, turning_index
!
        !If particle was lost and the recording quantities deallocated -> no further calculations (here for safety)
        if (.not.allocated(intermediate_z0_list)) return
!
        nskip = n_skip_vpar_0
!
!        !Convert t_pass to tau
!        tau = t_pass/dt_dtau_const
!
        !Select matrix and vector elements that are use for v_parallel integration
        a44 = amat(4,4)
        b4 = b(4)
!
if(diag_pusher_tetry_poly_adaptive) then
if(number_of_integration_steps.gt.1) then
print*, '----------------------------------'
do i = 1, number_of_integration_steps
print*, 'z0 number', i
print*, intermediate_z0_list(:,i)
print*, tau_steps_list(i)
enddo
print*, 'J_par before', par_adiab_inv
endif
endif
        !Trace banana tips
        if((vpar_end.gt.0.d0).and.(vpar_in.lt.0.d0)) then
            !Find exact root of vpar(tau)  --> tau_part1 is the "time" inside the tetrahedron until $v=\parallel$ = 0
            !turning_index = findloc(intermediate_z0_list(4,1:number_of_integration_steps) .gt. 0.0d0, .TRUE.) - 1
            turning_index = findloc(intermediate_z0_list(4,1:number_of_integration_steps) .gt. 0.0d0, .TRUE., dim = 1) - 1
            if (turning_index.eq.0) then
                !as then vpar > 0 at tetrahedron entry, which conflicts vpar_in < 0 of above
                print*, 'Error in par_adiab_inv_poly_mod: orbit already bounced before entering!'
                stop
            endif
            !If findloc does not find any match, aka orbit turned in LAST step -> returns 0 -> turning_index = -1
            if (turning_index.eq.-1) turning_index = number_of_integration_steps
!
            !Find now the root on the relevant orbit section determined above, to get the order-consistent result
            tau_part1 = tau_vpar_root(poly_order,a44,b4,intermediate_z0_list(4,turning_index))
!    
            if(tau_part1.gt.tau_steps_list(turning_index)) print *, 'Error in par_adiab_inv_poly_mod: bounce happens after section?'
!
            !Integrate par_adiab_inv exactly until root of vpar(tau)
            do i = 1, turning_index - 1
                par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau_steps_list(i), & 
                                                                & intermediate_z0_list(4,i))*dt_dtau_const
            enddo !integrating UP to turning point (before turn)
            !And finish from the turning index to the actual bounce
            par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau_part1, & 
                                                            & intermediate_z0_list(4,turning_index))*dt_dtau_const
!
            if(counter_banana_mappings.gt.1) then
                if(counter_banana_mappings/nskip*nskip.eq.counter_banana_mappings) then
!                    
                    !Poloidal projection of Poincaré sections at $v=\parallel$ = 0
                    !Need again to integrate from turning point to actual bounce
                    z = intermediate_z0_list(:,turning_index)
                    call set_integration_coef_manually(poly_order,z)
                    call analytic_integration_external(poly_order,z,tau_part1)
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
!
                    ! write coordinates for poincare cuts
                    !$omp critical
                        if(boole_poincare_vpar_0) then
                            write(file_id_vpar_0,*) x
                        endif
                    !$omp end critical
!
                    !$omp critical
                        if(boole_J_par) then
                            write(file_id_J_par,*) counter_banana_mappings,par_adiab_inv
                        endif
                    !$omp end critical
!
                    !$omp critical
                        if(boole_e_tot) then
                            write(file_id_e_tot,*) counter_banana_mappings,energy_tot_func(z,perpinv,ind_tetr)
                        endif
                    !$omp end critical
                    !print *, 'banana bounces:', counter_banana_mappings
                endif
            endif
            counter_banana_mappings = counter_banana_mappings + 1
!
            !Start to integrate par_adiab_inv for new bounce period
            par_adiab_inv = 0.d0
            !The new period starts at the bounce, so we finish the rest of the turning section
            par_adiab_inv = par_adiab_inv + &
                            & par_adiab_tau(poly_order,a44,b4,tau_steps_list(turning_index)-tau_part1,0.d0)*dt_dtau_const
            !And we then the remaining steps of total tetrahedron orbit
            do i = turning_index + 1, number_of_integration_steps
                par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau_steps_list(i), & 
                                                                & intermediate_z0_list(4,i))*dt_dtau_const
            enddo !integrating FROM end of turning section to the end of the total tetrahedron orbit
        else    
            !Compute parallel adiabatic invariant as a function of time
            do i = 1, number_of_integration_steps
                par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau_steps_list(i), & 
                                                                & intermediate_z0_list(4,i))*dt_dtau_const
            enddo !stepwise integration over whole tetrahedron orbit    
        endif
!
if(diag_pusher_tetry_poly_adaptive) then
if(number_of_integration_steps .gt. 1) then
print*, 'J_par after', par_adiab_inv
print*, 'turning_index', turning_index
print*, 'Bounce?', (vpar_end.gt.0.d0).and.(vpar_in.lt.0.d0)
endif  
endif 
!  
    end subroutine par_adiab_inv_tetra_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
    function par_adiab_tau(poly_order,a44,b4,tau,vpar_in)
!
        use pusher_tetra_poly_mod, only: sign_rhs
!
        implicit none
!
        integer                  :: poly_order
        double precision         :: tau,vpar_in,par_adiab_tau,a44,b4
!
        select case(poly_order)
            case(2)
                par_adiab_tau = tau*vpar_in**2 + 0.5d0*tau**2 * (2.d0*b4*vpar_in + 2.d0*a44*vpar_in**2) &
                            & + 1.d0/3.d0*tau**3*(b4**2 + 3.d0*a44*b4*vpar_in + 2.d0*a44**2* vpar_in**2)
!                            
            case(3)
                par_adiab_tau = tau*vpar_in**2 + 0.5d0*tau**2 * (2.d0*b4*vpar_in + 2.d0*a44*vpar_in**2) &
                            & + 1.d0/3.d0*tau**3*(b4**2 + 3.d0*a44*b4*vpar_in + 2.d0*a44**2*vpar_in**2) &
                            & + 1.d0/4.d0*tau**4*(a44*b4**2 + 7.d0/3.d0*a44**2*b4*vpar_in+(4.d0*a44**3*vpar_in**2)/3.d0)
!                
            case(4)
                par_adiab_tau = tau*vpar_in**2 + 0.5d0*tau**2 * (2.d0*b4*vpar_in + 2.d0*a44*vpar_in**2) &
                            & + 1.d0/3.d0*tau**3*(b4**2+3.d0*a44*b4*vpar_in+2.d0*a44*2.d0*vpar_in**2) &
                            & + 1.d0/4.d0*tau**4*(a44*b4**2 + 7.d0/3.d0*a44**2*b4*vpar_in+(4.d0*a44**3*vpar_in**2)/3.d0) &
                            & + 1.d0/60.d0*tau**5*(7.d0*a44**2*b4**2+15.d0*a44**3*b4*vpar_in + 8.d0*a44**4*vpar_in**2)
! 
        end select
!
    end function par_adiab_tau
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function tau_vpar_root(poly_order,a44,b4,vpar_in)
!
        use pusher_tetra_poly_mod, only: Quadratic_Solver2, Cubic_Solver, Quartic_Solver
!
        implicit none
!
        integer                                     :: poly_order,i_scaling
        double precision                            :: tau,vpar_in,tau_vpar_root,a44,b4
        double precision, dimension(poly_order+1)   :: coef_vec
!

        !Compute coefficients for vpar(t)==0 root finding
        !f(tau)  = a + b*tau + c/2 * tau^2 + d/6 * tau^3 + e/24 * tau^4
!
        i_scaling = 0
!
        !a
        coef_vec(1) = vpar_in
!
        if(poly_order.ge.1) then
            !b
            coef_vec(2) = (b4 + a44*vpar_in)
        end if
!
        if(poly_order.ge.2) then
            !c
            coef_vec(3) = (a44*b4 + a44**2*vpar_in)
        end if
!
        if(poly_order.ge.3) then
            !d
            coef_vec(4) = (a44**2*b4 + a44**3*vpar_in)
        end if
!
        if(poly_order.ge.4) then
            !e
            coef_vec(5) = (a44**3*b4 + a44**4*vpar_in)
        end if
!
        !Find root
        select case(poly_order)
            case(1)
                print *, 'Error: vpar root finding in linear order not available'
            case(2)
                call Quadratic_Solver2(coef_vec(3),coef_vec(2),coef_vec(1),tau_vpar_root)               
            case(3)        
                call Cubic_Solver(coef_vec(4),coef_vec(3),coef_vec(2),coef_vec(1),tau_vpar_root)
            case(4) 
                call Quartic_Solver(i_scaling,coef_vec(5),coef_vec(4),coef_vec(3), &
                                    & coef_vec(2),coef_vec(1),tau_vpar_root)
        end select  
!
    end function tau_vpar_root
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine analytic_integration_external(poly_order,z,tau)
!
        use poly_without_precomp_mod
!
        implicit none
!
        integer, intent(in)                          :: poly_order
        double precision, intent(in)                 :: tau
        double precision, dimension(4),intent(inout) :: z
        double precision                             :: tau2_half,tau3_sixth,tau4_twentyfourth
!
        !In contrast to the poly_pusher_mod counterpart, this integrator DOES NOT change the orbit recording quantities (for safety)
        if(poly_order.ge.1) then
            z = z + tau*(b+amat_in_z)
        endif
!
        if(poly_order.ge.2) then
            tau2_half = tau**2*0.5d0
            z = z + tau2_half*(amat_in_b + amat2_in_z)
        endif
!
        if(poly_order.ge.3) then
            tau3_sixth = tau**3/6.d0
            z = z + tau3_sixth*(amat2_in_b + amat3_in_z)
        endif

        if(poly_order.ge.4) then
            tau4_twentyfourth = tau**4/24.d0
            z = z + tau4_twentyfourth *(amat3_in_b + amat4_in_z)
        endif
!
    end subroutine analytic_integration_external
!
end module par_adiab_inv_poly_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
pure function scalar_multiplication_without_precomp(scalar_coef_1,scalar_coef_2)
use pusher_tetra_poly_mod, only: poly_multiplication_coef
implicit none
double precision, dimension(:), intent(in)           :: scalar_coef_1,scalar_coef_2
double precision, dimension(:), allocatable          :: scalar_multiplication_without_precomp
integer                                              :: k
!
k = maxval((/size(scalar_coef_1),size(scalar_coef_2)/))
allocate(scalar_multiplication_without_precomp(k))
call poly_multiplication_coef(scalar_coef_1,scalar_coef_2,scalar_multiplication_without_precomp)
!
end function scalar_multiplication_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
pure function vector_multiplication_without_precomp(scalar_coef,vector_coef)
use pusher_tetra_poly_mod, only: poly_multiplication_coef
implicit none
double precision, dimension(:), intent(in)           :: scalar_coef
double precision, dimension(:,:), intent(in)         :: vector_coef
double precision, dimension(:), allocatable          :: scalar_coef_res
integer                                              :: i,j,k
double precision, dimension(:,:), allocatable        :: vector_multiplication_without_precomp
!
j = size(vector_coef(:,1))
k = maxval((/size(vector_coef(1,:)),size(scalar_coef)/))
allocate(scalar_coef_res(k))
allocate(vector_multiplication_without_precomp(j,k))
!
do i = 1,j
        call poly_multiplication_coef(scalar_coef,vector_coef(i,:),scalar_coef_res)
        vector_multiplication_without_precomp(i,:) = scalar_coef_res
enddo
deallocate(scalar_coef_res)
end function vector_multiplication_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
pure function tensor_multiplication_without_precomp(vector_coef_1,vector_coef_2)
use pusher_tetra_poly_mod, only: poly_multiplication_coef
implicit none
double precision, dimension(:,:), intent(in)         :: vector_coef_1,vector_coef_2
integer                                              :: i,j,k,l,m
double precision, dimension(:), allocatable          :: scalar_coef_res
double precision, dimension(:,:,:), allocatable      :: tensor_multiplication_without_precomp
!
k = size(vector_coef_1(:,1))
l = size(vector_coef_2(:,1))
m = maxval((/size(vector_coef_1(1,:)),size(vector_coef_2(1,:))/))
allocate(scalar_coef_res(m))
allocate(tensor_multiplication_without_precomp(k,l,m))
!
        do i = 1,k
            do j = 1,l
                call poly_multiplication_coef(vector_coef_1(i,:),vector_coef_2(j,:),scalar_coef_res)
                tensor_multiplication_without_precomp(i,j,:) = scalar_coef_res
            enddo
        enddo
deallocate(scalar_coef_res)
!
end function tensor_multiplication_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
pure function scalar_integral_without_precomp(poly_order,tau,scalar_coef)
!
        !This function shall only be called within the subroutine "calc_optional_quantities"
        !Values in modules are used that need to be precomputed/set in that subroutine.
!
        implicit none
!
        integer, intent(in)                                  :: poly_order
        double precision, intent(in)                         :: tau
        double precision, dimension(:), intent(in)           :: scalar_coef

        double precision                                     :: scalar_integral_without_precomp
!
        if(poly_order.ge.1) then
            scalar_integral_without_precomp = scalar_coef(1)*tau + tau**2 * 0.5d0 * scalar_coef(2)
        endif
!
        if(poly_order.ge.2) then
            scalar_integral_without_precomp = scalar_integral_without_precomp + tau**3/3.d0 * scalar_coef(3)
        endif
!
        if(poly_order.ge.3) then
            scalar_integral_without_precomp = scalar_integral_without_precomp + tau**4/4.d0 * scalar_coef(4)
        endif
!
        if(poly_order.ge.4) then
            scalar_integral_without_precomp = scalar_integral_without_precomp + tau**5/5.d0 * scalar_coef(5)
        endif
!
    end function scalar_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
    pure function vector_integral_without_precomp(poly_order,tau,vector_coef)
!
        !This function shall only be called within the subroutine "calc_optional_quantities"
        !Values in modules are used that need to be precomputed/set in that subroutine.
    !
    use pusher_tetra_poly_mod, only: moment_integration
!
        implicit none

        
!
        integer, intent(in)                                   :: poly_order
        double precision, intent(in)                          :: tau
        double precision, dimension(:,:), intent(in)          :: vector_coef
!
        double precision, dimension(:), allocatable           :: vector_integral_without_precomp
!
        integer                                               :: i,m
!
        m = size(vector_coef(:,1))
        allocate(vector_integral_without_precomp(m))
        do i = 1,m
        vector_integral_without_precomp(i) = moment_integration(poly_order,tau,vector_coef(i,:))
        enddo
!
    end function vector_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    pure function tensor_integral_without_precomp(poly_order,tau,tensor_coef)
!
        !This function shall only be called within the subroutine "calc_optional_quantities"
        !Values in modules are used that need to be precomputed/set in that subroutine.
!
    use pusher_tetra_poly_mod, only: moment_integration
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(:,:,:), intent(in)  :: tensor_coef
!
        double precision, dimension(:,:), allocatable   :: tensor_integral_without_precomp
        integer :: j,k,m,n
!
        m = size(tensor_coef(:,1,1))
        n = size(tensor_coef(1,:,1))
        allocate(tensor_integral_without_precomp(m,n))
!
        do j = 1,m
            do k = 1,n
                tensor_integral_without_precomp(j,k) = moment_integration(poly_order,tau,tensor_coef(j,k,:))
            enddo
        enddo
!
    end function tensor_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!