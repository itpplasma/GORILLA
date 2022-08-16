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
    private
!    
    integer                             :: iface_init
    integer, public, protected          :: ind_tetr
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
    !change those for adaptive step sizes, probably allocatable
    double precision, dimension(:), allocatable           :: tau_steps_list
    double precision, dimension(:,:), allocatable         :: intermediate_z0_list
    integer                                               :: number_of_integration_steps
!
    !Diagnostics for adaptive step scheme (only useable for one particle calculation)
    double precision, dimension(:,:), allocatable         :: total_fluctuation_report, single_step_fluctuation_report, &
                                                           & closure_fluctuation_report
    integer                                               :: report_total_entry_index, report_single_step_entry_index, & 
                                                            & report_entry_index
    integer, parameter                                    :: max_n_report_entries = 1, minimum_partition_number = 9000
    logical                                               :: boole_collect_data
!
    !$OMP THREADPRIVATE(ind_tetr,iface_init,perpinv,perpinv2,dt_dtau_const,bmod0,t_remain,x_init,  &
    !$OMP& z_init,k1,k3,vmod0,tau_steps_list,intermediate_z0_list,number_of_integration_steps)
!
    public :: pusher_tetra_poly,initialize_const_motion_poly, initialize_intermediate_steps_arrays, &
        & Quadratic_Solver2, Cubic_Solver, Quartic_Solver,analytic_integration_without_precomp,energy_tot_func
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
        subroutine initialize_intermediate_steps_arrays(status)
!
            use gorilla_settings_mod, only: max_n_intermediate_steps, boole_adaptive_time_steps
            use gorilla_diag_mod, only: report_pusher_tetry_poly_adaptive
!
            implicit none
!            
            integer, intent(in) :: status
!
            integer             :: max_entries, k
!
            if (status .eq. 0) then
                !From one face to another there are max_n_intermediate_steps 
                !and that all possible two times as trajectory can be prolonged to a third face (turning on face) + buffer
                if (boole_adaptive_time_steps) then
                    max_entries = 3*(max_n_intermediate_steps)
                    allocate(tau_steps_list(max_entries),intermediate_z0_list(4,max_entries))
                else
                    !if no adaptive scheme -> just two steps in total
                    allocate(tau_steps_list(2),intermediate_z0_list(4,2))
                endif
if (report_pusher_tetry_poly_adaptive) then
allocate(total_fluctuation_report(max_n_report_entries*(max_n_intermediate_steps+1),2), & 
& single_step_fluctuation_report(max_n_report_entries*(max_n_intermediate_steps+1),3), &
& closure_fluctuation_report(max_n_report_entries*(max_n_intermediate_steps+1),3))
report_single_step_entry_index = 0
report_total_entry_index = 0
report_entry_index = 0
endif
            elseif (status .eq. 1) then
!
if (report_pusher_tetry_poly_adaptive) then
open(123, file='./total_fluctuation_report.dat')
do k = 1, report_total_entry_index
write(123,*) total_fluctuation_report(k,:)
end do
close(123)
open(123, file='./single_step_fluctuation_report.dat')
do k = 1, report_single_step_entry_index
write(123,*) single_step_fluctuation_report(k,:)
end do
close(123)
open(123, file='./closure_fluctuation_report.dat')
do k = 1, report_single_step_entry_index
write(123,*) closure_fluctuation_report(k,:)
end do
close(123)
deallocate(total_fluctuation_report,single_step_fluctuation_report,closure_fluctuation_report)
endif
!
                deallocate(tau_steps_list,intermediate_z0_list)
            else
                print*, 'Error: invalid input for initialization of intermediate_steps_arrays!'
                stop
            endif
!            
        end subroutine initialize_intermediate_steps_arrays
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
            use tetra_physics_mod, only: tetra_physics
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
            !Module of B at the entry point of the particle
            bmod0 = bmod_func(z_init(1:3)) !tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z_init(1:3))
!
            !Phi at the entry point of the particle
            phi_elec = phi_elec_func(z_init(1:3))   !tetra_physics(ind_tetr)%Phi1+sum(tetra_physics(ind_tetr)%gPhi*z_init(1:3))
!
            !Auxiliary quantities
            vperp2 = -2.d0*perpinv*bmod0
            vpar2 = vpar**2
            vmod0 = sqrt(vpar2+vperp2)
!
            k1 = vperp2+vpar2+2.d0*perpinv*tetra_physics(ind_tetr)%bmod1
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
            use gorilla_diag_mod,only: diag_pusher_tetry_poly
            use pusher_tetra_func_mod, only: pusher_handover2neighbour
            use gorilla_settings_mod, only: i_precomp, boole_guess, optional_quantities_type, boole_array_optional_quantities, &
                                    &  boole_adaptive_time_steps
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
!         
            call initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
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
            !!!!!!!!!!!!!!!!SECOND ATTEMPT WITHOUT GUESSES PLUS RESCALE!!!!!!!!!!!!!!!!!!
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
            t_pass = tau*dt_dtau_const
!
if(diag_pusher_tetry_poly) print *, 'tau total',tau
if(diag_pusher_tetry_poly) print *, 't_pass',t_pass
if(diag_pusher_tetry_poly) then
    print *, 't_remain',t_remain
    if (t_remain .lt. 0) stop
    if (t_pass .lt. 0) stop
endif

!
            !Particle stops inside the tetrahedron
            if(t_pass.ge.t_remain) then
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!FORTH ATTEMPT IF PARTICLE DOES NOT LEAVE CELL IN REMAINING TIME!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
                        return
                    endif
!
!
                    !Final processing
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    vpar=z(4)
                    t_pass = tau*dt_dtau_const
!
                    !Save relative coordinates after pushing
                    z_save = z(1:3)
!
                    iface = iface_new
                    call pusher_handover2neighbour(ind_tetr,ind_tetr_inout,iface,x,iper_phi)
                endif
!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!END OF FORTH ATTEMPT IF PARTICLE DOES NOT LEAVE CELL!!!!!!!!!!!!!
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
            ! deallocate(tau_steps_list,intermediate_z0_list)
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
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: three planes'
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
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: distance'
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
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: normal velocity'
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
            print *, 'Error: No analytical solution exists.'
            ind_tetr_inout = -1
            iface = -1
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
        use gorilla_settings_mod, only: delta_energy, max_n_intermediate_steps, min_tau_intermediate_steps
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
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
        double precision                                      :: energy_start, energy_current, delta_energy_current, tau_min
        double precision, dimension(4)                        :: z0
        integer                                               :: dummy
!
        if ((delta_energy .le. 0)) then
            print*, 'Error: The control setting delta_energy is invalid! Check the limits in gorilla.inp!'
            stop
        elseif ((max_n_intermediate_steps .lt. 1)) then
            print*, 'Error: The control setting max_n_intermediate_steps is invalid! Check the limits in gorilla.inp!'
            stop
        elseif ((min_tau_intermediate_steps .gt. 0.5)) then
            print*, 'Error: The control setting min_tau_intermediate_steps is invalid! Check the limits in gorilla.inp!'
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
        dummy = 1 !Needed to initialize the number of intermediate steps
        tau_min = tau*min_tau_intermediate_steps !Needed for limitting the splitting of steps
        z0 = intermediate_z0_list(:,number_of_integration_steps)
if(diag_pusher_tetry_poly) print*, 'z0', z0
        energy_start = energy_tot_func(z0)
if(diag_pusher_tetry_poly) print*, 'energy_start', energy_start
        energy_current = energy_tot_func(z)
if(diag_pusher_tetry_poly) print*, 'energy_end', energy_current
        delta_energy_current = abs(1-energy_current/energy_start)
        !If energy fluctuating too strong -> start recursive adaptive scheme
        if (delta_energy_current .gt. delta_energy) then
            if(.false.) then
if(diag_pusher_tetry_poly) print*, 'Adaptive Stepsize recursive was called'
                call adaptive_time_steps_recursive(poly_order, i_scaling, boole_guess_adaptive, boole_passing, tau_min, &
                                &  iface_inout_adaptive, dummy, tau, z, boole_face_correct)
            else
if(diag_pusher_tetry_poly) print*, 'Adaptive Stepsize equidistant was called'
                call adaptive_time_steps_equidistant(poly_order, i_scaling, boole_guess_adaptive, boole_passing, &
                                & delta_energy_current, iface_inout_adaptive, tau, z, boole_face_correct)
            endif
        endif
!
    end subroutine overhead_adaptive_time_steps
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    recursive subroutine adaptive_time_steps_recursive(poly_order, i_scaling, boole_guess_adaptive, boole_passing, tau_min, &
                                    &  iface_inout_adaptive, n_intermediate_steps, tau, z, boole_face_correct)
!
        use gorilla_settings_mod, only: i_precomp, delta_energy, max_n_intermediate_steps
        use gorilla_diag_mod,only: diag_pusher_tetry_poly
!
        implicit none 
!
        integer, intent(in)                                   :: poly_order, i_scaling
        logical, intent(in)                                   :: boole_guess_adaptive
        logical, intent(in)                                   :: boole_passing
        double precision, intent(in)                          :: tau_min
!
        integer, intent(inout)                                :: iface_inout_adaptive, n_intermediate_steps
        double precision, intent(inout)                       :: tau
        double precision, dimension(4), intent(inout)         :: z
        logical, intent(inout)                                :: boole_face_correct
!
        logical, dimension(4)                                 :: boole_faces
        integer                                               :: i
        integer                                               :: iface_new_adaptive
        logical                                               :: boole_analytical_approx
        double precision                                      :: energy_init, energy_current, tau_half, tau_prime
!
        !Start with previous face, set z back to current z0
        iface_new_adaptive = iface_inout_adaptive
if(diag_pusher_tetry_poly) print *, 'integrated z', normal_distance_func(z(1:3),2)
        z = intermediate_z0_list(:,number_of_integration_steps)
if(diag_pusher_tetry_poly) print *, 'starting z', normal_distance_func(z(1:3),2)
!
        !Then redo the integration with half timestep
        number_of_integration_steps = number_of_integration_steps - 1
        tau_half = tau/2
        call analytic_integration(poly_order,i_precomp,z,tau_half)
!
        !Needed for recursion (can be set to true, as the whole adaptive scheme takes place before the other consistency checks)
        boole_face_correct = .true.
!
if(diag_pusher_tetry_poly) print *, 'Adaptive step begin: tau_half', tau_half
if(diag_pusher_tetry_poly) print *, 'Adaptive step begin: minimal tau allowed', tau_min
!        
        !Particle must be still inside the tetrahedron
        !If fail -> return to upper level without closure or another call -> climibing out the recusion scheme
        do i = 1,4
            if(normal_distance_func(z(1:3),i).lt.0.d0) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error in adaptive steps: position outside tetrahedron'
if(diag_pusher_tetry_poly) print *, i, 'norm', normal_distance_func(z(1:3),i)
                return
            endif
        enddo
        !Particle is now not at any of the exit faces
        if (boole_face_correct) iface_new_adaptive = 0
!
        !If this intermediate step still not fullfills energy conservation (or left tetrahedron if at top level guess was made)
        !The routine calls itself (recursivly), again redoing the step with additional halfed step if not at MINIMUM STEP SIZE
        energy_init = energy_tot_func(z_init)
        energy_current = energy_tot_func(z)
        if ((abs(1-energy_current/energy_init) .gt. delta_energy) .AND. (tau_half/2 .ge. tau_min)) then
if(diag_pusher_tetry_poly) print *, 'DOWN'
            call adaptive_time_steps_recursive(poly_order, i_scaling, boole_guess_adaptive, boole_passing, tau_min, &
                                    &  iface_new_adaptive, n_intermediate_steps, tau_half, z, boole_face_correct)
!
        !If everything fine (or at least inside if reach stepsize limit), we try to finish the orbit
        else
            !If orbit passes, the now remaining time (tau_prime) from the intermediate point to the exit face has to be recalcula
            if (boole_passing) then
!
                boole_faces = .true.
!
                !use face prediction of second order for higher order computation in the adaptive scheme if wished
                if (boole_guess_adaptive .and. (poly_order.gt.2)) then
                    call analytic_approx(2,i_precomp,boole_faces, &
                    & i_scaling,z,iface_new_adaptive,tau_prime,boole_analytical_approx)
                    !Higher order polynomial root is computed only for previously predicted face in second order
                    if(boole_analytical_approx) then               
                        boole_faces = .false.                !Disable all 4 faces 
                        boole_faces(iface_new_adaptive) = .true.      !Enable guessed face
                    endif
                    !Reset starting face for actual correct order calculation down below
                    iface_new_adaptive = iface_inout_adaptive 
                endif  
!           
                !calculate exit time and exit face in correct order
                !if a successful guess was made above in second order, boole_faces only allows the guessed face 
                call analytic_approx(poly_order,i_precomp,boole_faces, &
                                    & i_scaling,z,iface_new_adaptive,tau_prime,boole_analytical_approx)   
!
if(diag_pusher_tetry_poly) print *, 'Adaptive step closure: tau_prime', tau_prime
if(diag_pusher_tetry_poly) print *, 'Adaptive step closure: minimal tau allowed', tau_min
!           
                !Analytical result does not exist can therefore not close cell orbit
                ! -> "return" leaves orbit (falsely therefore boole set to false) inside of tetrahedron at this recursion level
                !Then the upper recursion levels get resolved with the orbit still there and bool_face_correct set to false
                if(.not.boole_analytical_approx) then
if(diag_pusher_tetry_poly) print *, 'Error in adaptive steps: no analytical solution'
                    boole_face_correct = .false.
                    return
                endif
!
            !If the orbit does not pass, the remaining time after the intermediate one is just used
            else
!
                tau_prime = tau - tau_half
!
            end if !boole_passing
!
            !Integrate trajectory analytically, if root exists or not passing anyway
            call analytic_integration(poly_order,i_precomp,z,tau_prime)
!
            !If this closing step not fullfills energy conservation and we have not hit the MAXIMUM NUMBER OF INTERMEDIATE STEPS
            !The routine calls itself (recursivly), again redoing the step now with the remaining time (tau_prime) halfed
            !One can not detect if the orbit ventured outside, as we rerooted it to finish (if energy fits this is the last step)
            !Exception is in the case of the none passing orbit, but this is checked independantly in the top level routine
            energy_init = energy_tot_func(z_init)
            energy_current = energy_tot_func(z)
            if ((abs(1-energy_current/energy_init).gt.delta_energy).AND.(n_intermediate_steps.lt.max_n_intermediate_steps) &
                    & .AND. (tau_prime/2 .ge. tau_min)) then
                !This creates an ADDITIONAL intermediate step (in contrast to the above which just changes the size of the step)
if(diag_pusher_tetry_poly) print *, 'DOWN'
                n_intermediate_steps = n_intermediate_steps + 1
                call adaptive_time_steps_recursive(poly_order, i_scaling, boole_guess_adaptive, boole_passing, tau_min, &
                                        &  iface_new_adaptive, n_intermediate_steps, tau_prime, z, boole_face_correct)
            else
                !Gives the same duration of total step if no rerooting had to be done (aka not passing orbit)
                tau = tau_half + tau_prime
            endif ! energy fluctuation > delta_energy in closing step
!
        endif ! energy fluctuation > delta_energy or outside of tetrahedron in normal step
!
        iface_inout_adaptive = iface_new_adaptive
if(diag_pusher_tetry_poly) print *, 'UP'
! if(diag_pusher_tetry_poly) print *, 'ending z', normal_distance_func(z(1:3),2)
! if(diag_pusher_tetry_poly) print *, 'leaving face', iface_inout_adaptive
!
    end subroutine adaptive_time_steps_recursive
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine adaptive_time_steps_equidistant(poly_order, i_scaling, boole_guess_adaptive, boole_passing, &
        &  delta_energy_current, iface_out_adaptive, tau, z, boole_face_correct)
!
        use gorilla_settings_mod, only: i_precomp, delta_energy, max_n_intermediate_steps
        use gorilla_diag_mod,only: diag_pusher_tetry_poly, diag_pusher_tetry_poly_adaptive, report_pusher_tetry_poly_adaptive
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
        integer                                               :: i, k, n_intermediate_steps, limit, &
                                                                & n_intermediate_steps_minimum
        integer                                               :: iface_new_adaptive, number_of_integration_steps_start_adaptive
        logical                                               :: boole_analytical_approx, boole_exit_tetrahedron, & 
                                                               & boole_energy_check, boole_reached_minimum
        double precision                                      :: energy_start_adaptive, energy_current, tau_prime, tau_collected, &
                                                               & scale_factor, max_scale_factor, delta_energy_start, &
                                                               & delta_energy_minimum
        double precision, parameter                           :: freshhold = 1, min_step_error = 1E-15, additive_increase = 10!10,2,1
        double precision, dimension(4,poly_order+1)           :: dummy_coef_mat !only here for analytic_coeff_without_precomp syntax
        logical, dimension(4)                                 :: boole_faces_off = .false. !avoid root coefficients computation

!
        !Set z back to current z0, step back in global integration counter and save
if(diag_pusher_tetry_poly_adaptive) print *, 'integrated z', normal_distance_func(z(1:3),2)
        z = intermediate_z0_list(:,number_of_integration_steps)
if(diag_pusher_tetry_poly_adaptive) print *, 'previous z', normal_distance_func(z(1:3),2)
        number_of_integration_steps = number_of_integration_steps - 1
        number_of_integration_steps_start_adaptive = number_of_integration_steps
        z_start_adaptive = z
if(report_pusher_tetry_poly_adaptive) boole_collect_data = .false.
!
        !Set up for Partition procedure
        n_intermediate_steps = 1
        boole_reached_minimum = .false.
        delta_energy_start = delta_energy_current
        delta_energy_minimum = delta_energy_current
        !A bit more than the maximal "desired" number of steps (need more/less 
        !steps as stepwise integration changes path and total dwell time, escpecially in relevant cases)
        limit = max_n_intermediate_steps*1.5d0
!
        !Loop over possible equidistant splittings with increasing number of steps
        PARTITION: do while (n_intermediate_steps .lt. max_n_intermediate_steps)
!
            !If did not succeed (aka not leave the loop) -> try to increase number of steps
            !and try again, if not yet at maximal number of intermediate steps or passed minimum
            scale_factor = (delta_energy_current/delta_energy)**(1.0d0/poly_order)
            max_scale_factor = (delta_energy_current/(min_step_error*n_intermediate_steps))**(1.0d0/(poly_order+1))
if (diag_pusher_tetry_poly_adaptive) print*, 'scale factor', scale_factor
if (diag_pusher_tetry_poly_adaptive) print*, 'max scale factor', max_scale_factor
           if (max_scale_factor .gt. freshhold) then
                scale_factor = min(scale_factor,max_scale_factor)
                n_intermediate_steps = min(ceiling(n_intermediate_steps*scale_factor),max_n_intermediate_steps)
            elseif (scale_factor .gt. freshhold) then
                n_intermediate_steps = min(ceiling(n_intermediate_steps*scale_factor),max_n_intermediate_steps)
            else
                n_intermediate_steps = n_intermediate_steps + additive_increase
            endif !Choice of next number of steps
!
            !Set up for redo of the minimum run
            if (boole_reached_minimum) n_intermediate_steps = n_intermediate_steps_minimum
!
if (report_pusher_tetry_poly_adaptive) then
if(.not.boole_collect_data.AND.(n_intermediate_steps.gt.minimum_partition_number)) then
report_entry_index = report_entry_index + 1
boole_collect_data = .true.
n_intermediate_steps = 1
delta_energy_current = delta_energy_start
report_single_step_entry_index = report_single_step_entry_index + 1
report_total_entry_index = report_total_entry_index + 1
single_step_fluctuation_report(report_single_step_entry_index,1) = n_intermediate_steps
single_step_fluctuation_report(report_single_step_entry_index,2) = delta_energy_current
single_step_fluctuation_report(report_single_step_entry_index,3) = tau
total_fluctuation_report(report_total_entry_index,1) = n_intermediate_steps
total_fluctuation_report(report_total_entry_index,2) = delta_energy_current
closure_fluctuation_report(report_single_step_entry_index,1) = n_intermediate_steps
closure_fluctuation_report(report_single_step_entry_index,2) = delta_energy_current
closure_fluctuation_report(report_single_step_entry_index,3) = 1
cycle PARTITION
endif
endif
!
            !Set up for every partition trial; boole_face_correct can be set true here as adaptive is before the consistency checks
            z = z_start_adaptive
            tau_prime = tau/n_intermediate_steps
if(diag_pusher_tetry_poly_adaptive) print *, 'tau_prime', tau_prime
!if(diag_pusher_tetry_poly_adaptive) print *, 'n_intermediate_steps', n_intermediate_steps
            number_of_integration_steps = number_of_integration_steps_start_adaptive
            tau_collected = 0
            boole_face_correct = .true.
            boole_exit_tetrahedron = .false.
            boole_energy_check = .false.
            if (.not.boole_passing) limit = n_intermediate_steps
!
            STEPWISE: do i = 1, (limit -1) !used to be n_intermediate_steps
!
                !recalculate polynomial coefficients (tensors) as at every intermediate step the poly_coeff change
                call analytic_coeff_without_precomp(poly_order,boole_faces_off,z,dummy_coef_mat)
                call analytic_integration(poly_order,i_precomp,z,tau_prime)
!
!if(diag_pusher_tetry_poly_adaptive) print *, 'loop z exit', normal_distance_func(z(1:3),1)
!if(diag_pusher_tetry_poly_adaptive) print *, 'loop z compare', normal_distance_func(z(1:3),3)
                !Particle must be still inside the tetrahedron
                !If fail -> return to last step and close orbit then
                CONTROL: do k = 1,4
!
                    if(normal_distance_func(z(1:3),k).lt.0.d0) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps: position outside tetrahedron'
if(diag_pusher_tetry_poly_adaptive) print *, k, 'norm', normal_distance_func(z(1:3),k)
if(diag_pusher_tetry_poly_adaptive) print *, 'steps taken', i, '/' ,n_intermediate_steps
                        z = intermediate_z0_list(:,number_of_integration_steps)
                        !Exception: if the first step is already outside -> need smaller steps
                        if (i .eq. 1) cycle PARTITION
                        boole_exit_tetrahedron = .true.
                        exit STEPWISE
                    endif
!
                end do CONTROL !Is the orbit still inside
!
                !If it was a valid step we add it to the total used time
                tau_collected = tau_collected + tau_prime
!
if (report_pusher_tetry_poly_adaptive) then
if (boole_collect_data ) then
if (i.eq.1) then
report_single_step_entry_index = report_single_step_entry_index + 1
single_step_fluctuation_report(report_single_step_entry_index,1) = 0
single_step_fluctuation_report(report_single_step_entry_index,2) = 0
single_step_fluctuation_report(report_single_step_entry_index,3) = tau_prime
endif
single_step_fluctuation_report(report_single_step_entry_index,2) = &
& single_step_fluctuation_report(report_single_step_entry_index,2) + abs(1-energy_tot_func(z)/ & 
& energy_tot_func(intermediate_z0_list(:,number_of_integration_steps)))
single_step_fluctuation_report(report_single_step_entry_index,1) = & 
& single_step_fluctuation_report(report_single_step_entry_index,1) + 1
endif
endif
!
            end do STEPWISE ! stepwise integration
!
            !If orbit passes, the now remaining time (tau_prime) from the intermediate point to the exit face has to be recalculatet
            !Alternatively stepwise integration might reveal a not passing orbit to be a passing one, which also has to be closed
            if (boole_passing .OR. boole_exit_tetrahedron) then
!
                !If not the first step brought us outside (else this part of the code is not reached), the orbit is now INSIDE
                iface_new_adaptive = 0
                boole_faces = .true.
!
                !use face prediction of second order for higher order computation in the adaptive scheme if wished
                if (boole_guess_adaptive .and. (poly_order.gt.2)) then
                    call analytic_approx(2,i_precomp,boole_faces, &
                    & i_scaling,z,iface_new_adaptive,tau_prime,boole_analytical_approx)
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
                                    & i_scaling,z,iface_new_adaptive,tau_prime,boole_analytical_approx)   
!
if(diag_pusher_tetry_poly_adaptive) print *, 'Adaptive step closure: tau_prime', tau_prime
!           
                !Analytical result does not exist can therefore not close cell orbit
                ! -> "return" leaves orbit (falsely, therefore boole set to false) inside of tetrahedron
                if(.not.boole_analytical_approx) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps: no analytical solution'
                    boole_face_correct = .false.
                    return
                endif
!
            !If the orbit does not pass, the same timestep as for the other intermediate steps is used
            !Should it then end up outside, the higher level checks will detect it
            !We however have to manually update again the new position in the coefficiants
            else
                call analytic_coeff_without_precomp(poly_order,boole_faces_off,z,dummy_coef_mat)
            end if !boole_passing
!
            !Closing integration (either to the exit face or to the final position inside of tetrahedron)
            call analytic_integration(poly_order,i_precomp,z,tau_prime)
            !Adding up the last step to the total
            tau_collected = tau_collected + tau_prime
!
            !Check if energy fluctuation was successfully decreased
            energy_start_adaptive = energy_tot_func(z_start_adaptive)
            energy_current = energy_tot_func(z)
            delta_energy_current = abs(1-energy_current/energy_start_adaptive)
!
if (diag_pusher_tetry_poly_adaptive) print*, 'current delta_energy', delta_energy_current
if (diag_pusher_tetry_poly_adaptive) print*, 'delta_energy_minimum', delta_energy_minimum
if (diag_pusher_tetry_poly_adaptive) print*, 'n_intermediate_steps', n_intermediate_steps
if (report_pusher_tetry_poly_adaptive) delta_energy_minimum = delta_energy_current*2 !Circumvent the minimum approach for more data
if (report_pusher_tetry_poly_adaptive) then
if (boole_collect_data) then
single_step_fluctuation_report(report_single_step_entry_index,2) = &
& single_step_fluctuation_report(report_single_step_entry_index,2) + abs(1-energy_tot_func(z)/ & 
& energy_tot_func(intermediate_z0_list(:,number_of_integration_steps)))
single_step_fluctuation_report(report_single_step_entry_index,1) = & 
& single_step_fluctuation_report(report_single_step_entry_index,1) + 1
!Forming the average of the single step error done in this partition  
single_step_fluctuation_report(report_single_step_entry_index,2) = & 
& single_step_fluctuation_report(report_single_step_entry_index,2) / & 
& single_step_fluctuation_report(report_single_step_entry_index,1)
report_total_entry_index = report_total_entry_index + 1
total_fluctuation_report(report_total_entry_index,2) = delta_energy_current
total_fluctuation_report(report_total_entry_index,1) = single_step_fluctuation_report(report_single_step_entry_index,1)
closure_fluctuation_report(report_single_step_entry_index,1) = single_step_fluctuation_report(report_single_step_entry_index,1)
closure_fluctuation_report(report_single_step_entry_index,2) = abs(1-energy_tot_func(z)/ & 
& energy_tot_func(intermediate_z0_list(:,number_of_integration_steps)))
closure_fluctuation_report(report_single_step_entry_index,3) = & 
& tau_prime/single_step_fluctuation_report(report_single_step_entry_index,3)
endif
endif
!
            !There is an effective minimum that can be achieved by decreasing the timesteps
            !After that the error actually increases/osscilates -> we only use the monoton decrease
            if (boole_reached_minimum) then
                exit PARTITION
            elseif (delta_energy_current.lt.delta_energy_minimum) then
                delta_energy_minimum = delta_energy_current
                n_intermediate_steps_minimum = n_intermediate_steps
                if (delta_energy_current.le.delta_energy) then
                    boole_energy_check = .true.
                    exit PARTITION
                endif
            else
                !If we see that our new step choice has INCREASED the energy, we redo the last partition and stop
                boole_reached_minimum = .true.
                cycle PARTITION
            endif !energy check
!
        end do PARTITION
!
if (report_pusher_tetry_poly_adaptive) then
if (boole_collect_data) then
if (report_entry_index.ge.max_n_report_entries) then
if (diag_pusher_tetry_poly_adaptive) print*, n_intermediate_steps
if (diag_pusher_tetry_poly_adaptive) print*, boole_exit_tetrahedron
call initialize_intermediate_steps_arrays(1)
print *, 'Collected enough data for adaptive scheme report!'
stop
endif
report_single_step_entry_index = report_single_step_entry_index + 1
report_total_entry_index = report_total_entry_index + 1
single_step_fluctuation_report(report_single_step_entry_index,:) = -1
total_fluctuation_report(report_total_entry_index,:) = -1
closure_fluctuation_report(report_single_step_entry_index,:) = -1
endif
endif
!
        !If none of the partitions could fully satisfy the energy condition, failure is given back to the higher level routine
        if (.not.boole_energy_check) then
if(diag_pusher_tetry_poly_adaptive) print *, 'Error in adaptive steps equidistant: energy conservation could not be fullfilled!'
            !boole_face_correct = .false.
            !print*, delta_energy_current
if(diag_pusher_tetry_poly_adaptive .AND. (delta_energy_minimum.gt.1.0d-10)) stop
        endif !adaptive sucess check
!
        iface_out_adaptive = iface_new_adaptive
        tau = tau_collected
! if (diag_pusher_tetry_poly_adaptive) print*, iface_out_adaptive
! if (diag_pusher_tetry_poly_adaptive) print*, tau
!
    end subroutine adaptive_time_steps_equidistant
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
        integer, intent(in) :: i_scaling
        double precision, intent(in)  :: a,b,c
        double precision,intent(out)                :: dtau
!
        select case(i_scaling)
            case(0)
                call Quadratic_Solver1(a,b,c,dtau)
            case DEFAULT
                print *, 'New quadratic solver is called.'
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
            use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
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
            optional_quantities%t_hamiltonian = 0
            optional_quantities%gyrophase = 0
!
        end subroutine initialise_optional_quantities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine calc_optional_quantities(poly_order,z0,tau,optional_quantities)
!
        use gorilla_settings_mod, only: boole_time_Hamiltonian, boole_gyrophase, optional_quantities_type
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
        !Save coordinates before analytical integration
        x0 = z0(1:3)
        vpar0 = z0(4)

        !Optional computation of Hamiltonian time
        if(boole_time_Hamiltonian) then
!
            allocate(x_coef(3,poly_order+1))
            allocate(vpar_coef(poly_order+1))
            allocate(x_vpar_coef(3,poly_order+1))
!
            call x_series_coef(x_coef,poly_order,x0)
            call vpar_series_coef(vpar_coef,poly_order,vpar0)
!
            call poly_multiplication_coef(x_coef(1,:),vpar_coef(:),x_vpar_coef(1,:))
            call poly_multiplication_coef(x_coef(2,:),vpar_coef(:),x_vpar_coef(2,:))
            call poly_multiplication_coef(x_coef(3,:),vpar_coef(:),x_vpar_coef(3,:))
!
            delta_t_hamiltonian = hamiltonian_time(ind_tetr)%h1_in_curlA * tau + &
            & cm_over_e * hamiltonian_time(ind_tetr)%h1_in_curlh * vpar_integral_without_precomp(poly_order,tau,vpar_coef)+&
            & sum( hamiltonian_time(ind_tetr)%vec_mismatch_der * x_integral_without_precomp(poly_order,tau,x_coef) ) + &
            & cm_over_e * sum(hamiltonian_time(ind_tetr)%vec_parcurr_der *  &
            & x_vpar_integral_without_precomp(poly_order,tau,x_vpar_coef))
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
                    & 1.d0/cm_over_e * sum( tetra_physics(ind_tetr)%gb * x_integral_without_precomp(poly_order,tau,x_coef) ) * &
                    & hamiltonian_time(ind_tetr)%h1_in_curlA + &
!
                    !Second term
                    & sum(x_vpar_integral_without_precomp(poly_order,tau,x_vpar_coef) * tetra_physics(ind_tetr)%gb ) * &
                    & hamiltonian_time(ind_tetr)%h1_in_curlh + &
!
                    !Third term
                    & 1.d0/cm_over_e * sum( matmul( xx_integral_without_precomp(poly_order,tau,x_coef) , &
                    & hamiltonian_time(ind_tetr)%vec_mismatch_der ) * tetra_physics(ind_tetr)%gb) + &
!
                    !Fourth term
                    & sum( matmul( xxvpar_integral_without_precomp( poly_order,tau,x_coef,x_vpar_coef), &
                    & hamiltonian_time(ind_tetr)%vec_parcurr_der ) *  tetra_physics(ind_tetr)%gb) &
                & )

            endif ! gyrophase
!            
            deallocate(x_coef,vpar_coef,x_vpar_coef)
!
        endif ! time_Hamiltonian
!
    end subroutine calc_optional_quantities
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine x_series_coef(x_coef,poly_order,x0)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        use poly_without_precomp_mod, only: b, amat_in_z,amat2_in_z,amat3_in_z,amat_in_b,amat2_in_b, &
                                            & amat4_in_z,amat3_in_b
!
        implicit none
!
        double precision, dimension(:,:),intent(out)    :: x_coef
        integer, intent(in)                             :: poly_order
        double precision, dimension(3), intent(in)      :: x0
!
        x_coef(:,1) = x0
!
        if(poly_order.ge.1) then
            x_coef(:,2) = b(1:3) + amat_in_z(1:3)
        endif
!
        if(poly_order.ge.2) then
            x_coef(:,3) = 0.5d0 * (amat_in_b(1:3) + amat2_in_z(1:3))
        endif
!
        if(poly_order.ge.3) then
            x_coef(:,4) = 1.d0 / 6.d0 * (amat2_in_b(1:3) + amat3_in_z(1:3))
        endif
!
        if(poly_order.ge.4) then
            x_coef(:,5) = 1.d0 / 24.d0 * (amat3_in_b(1:3) + amat4_in_z(1:3))
        endif
!
    end subroutine x_series_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vpar_series_coef(vpar_coef,poly_order,vpar0)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        use poly_without_precomp_mod, only: amat,b
!
        implicit none
!
        double precision, dimension(:),intent(out)      :: vpar_coef
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: vpar0
        double precision                                :: a44,b4
!
        a44 = amat(4,4)
        b4 = b(4)
!
        vpar_coef(1) = vpar0
!
        if(poly_order.ge.1) then
            vpar_coef(2) = vpar0*a44 + b4
        endif
!
        if(poly_order.ge.2) then
            vpar_coef(3) = 0.5d0 * (vpar0 * a44**2 + a44 * b4)
        endif
!
        if(poly_order.ge.3) then
            vpar_coef(4) = 1.d0 / 6.d0 * (vpar0 * a44**3 + a44**2 * b4)
        endif
!
        if(poly_order.ge.4) then
            vpar_coef(5) = 1.d0 / 24.d0 * (vpar0 * a44**4 + a44**3 * b4)
        endif

    end subroutine vpar_series_coef
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function x_integral_without_precomp(poly_order,tau,x_coef)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(3)                  :: x_integral_without_precomp
        double precision, dimension(:,:), intent(in)    :: x_coef
!
        if(poly_order.ge.1) then
            x_integral_without_precomp = x_coef(:,1) * tau + x_coef(:,2) * 0.5d0 * tau**2
        endif
!
        if(poly_order.ge.2) then
            x_integral_without_precomp = x_integral_without_precomp + tau**3/3.d0 * x_coef(:,3)
        endif
!
        if(poly_order.ge.3) then
            x_integral_without_precomp = x_integral_without_precomp + tau**4/4.d0 * x_coef(:,4)
        endif
!
        if(poly_order.ge.4) then
            x_integral_without_precomp = x_integral_without_precomp + tau**5/5.d0 * x_coef(:,5)
        endif

    end function x_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function vpar_integral_without_precomp(poly_order,tau,vpar_coef)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        implicit none
!
        integer, intent(in)                         :: poly_order
        double precision, intent(in)                :: tau
        double precision, dimension(:),intent(in)   :: vpar_coef
        double precision                            :: vpar_integral_without_precomp
!
        if(poly_order.ge.1) then
            vpar_integral_without_precomp = vpar_coef(1)*tau + tau**2 * 0.5d0 * vpar_coef(2)
        endif
!
        if(poly_order.ge.2) then
            vpar_integral_without_precomp = vpar_integral_without_precomp + tau**3/3.d0 * vpar_coef(3)
        endif
!
        if(poly_order.ge.3) then
            vpar_integral_without_precomp = vpar_integral_without_precomp + tau**4/4.d0 * vpar_coef(4)
        endif
!
        if(poly_order.ge.4) then
            vpar_integral_without_precomp = vpar_integral_without_precomp + tau**5/5.d0 * vpar_coef(5)
        endif
!
    end function vpar_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine poly_multiplication_coef(poly_coef_1,poly_coef_2,res_poly_coef)
!
        implicit none
!
        double precision, dimension(:),intent(in) :: poly_coef_1,poly_coef_2
        double precision, dimension(:),intent(out) :: res_poly_coef
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
    function x_vpar_integral_without_precomp(poly_order,tau,x_vpar_coef)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(3)                  :: x_vpar_integral_without_precomp
        double precision, dimension(:,:), intent(in)    :: x_vpar_coef
!
        if(poly_order.ge.1) then
            x_vpar_integral_without_precomp = x_vpar_coef(:,1) * tau + x_vpar_coef(:,2) * 0.5d0 * tau**2
        endif
!
        if(poly_order.ge.2) then
            x_vpar_integral_without_precomp = x_vpar_integral_without_precomp + tau**3/3.d0 * x_vpar_coef(:,3)
        endif
!
        if(poly_order.ge.3) then
            x_vpar_integral_without_precomp = x_vpar_integral_without_precomp + tau**4/4.d0 * x_vpar_coef(:,4)
        endif
!
        if(poly_order.ge.4) then
            x_vpar_integral_without_precomp = x_vpar_integral_without_precomp + tau**5/5.d0 * x_vpar_coef(:,5)
        endif
!
    end function x_vpar_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function xx_integral_without_precomp(poly_order,tau,x_coef)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(3,3)                :: xx_integral_without_precomp
        double precision, dimension(:,:), intent(in)    :: x_coef
        double precision, dimension(:), allocatable     :: xx_coef_component
        integer :: j,k
!
        allocate(xx_coef_component(poly_order+1))
!
        do j = 1,3
            do k = 1,3
!
                !Compute convolution of coefficient for one component
                call poly_multiplication_coef(x_coef(j,:),x_coef(k,:),xx_coef_component)
!
                xx_integral_without_precomp(j,k) = tau * xx_coef_component(1)
!
                if(poly_order.ge.1) then
                    xx_integral_without_precomp(j,k) = xx_integral_without_precomp(j,k) +  &
                        & tau**2 * 0.5d0 * xx_coef_component(2)
                endif
!
                if(poly_order.ge.2) then
                    xx_integral_without_precomp(j,k) = xx_integral_without_precomp(j,k) + &
                        & tau**3 / 3.d0 * xx_coef_component(3)
                endif
!
                if(poly_order.ge.3) then
                    xx_integral_without_precomp(j,k) = xx_integral_without_precomp(j,k) + &
                        & tau**4 / 4.d0 * xx_coef_component(4)
                endif
!
                if(poly_order.ge.4) then
                    xx_integral_without_precomp(j,k) = xx_integral_without_precomp(j,k) + &
                        & tau**5 / 5.d0 * xx_coef_component(5)
                endif
!
            enddo
        enddo
!
        deallocate(xx_coef_component)
!
    end function xx_integral_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function xxvpar_integral_without_precomp(poly_order,tau,x_coef,x_vpar_coef)
!
        !This function shall only be called within the subroutine "analytic_integration_without_precomp"
        !Values in modules are used that need to be precomputed in that subroutine.
!
        implicit none
!
        integer, intent(in)                             :: poly_order
        double precision, intent(in)                    :: tau
        double precision, dimension(3,3)                :: xxvpar_integral_without_precomp
        double precision, dimension(:,:), intent(in)    :: x_coef,x_vpar_coef
        double precision, dimension(:), allocatable     :: xxvpar_coef_component
        integer :: j,k
!
        allocate(xxvpar_coef_component(poly_order+1))
!
        do j = 1,3
            do k = 1,3
!
                !Compute convolution of coefficient for one component
                call poly_multiplication_coef(x_vpar_coef(j,:),x_coef(k,:),xxvpar_coef_component)
!
                xxvpar_integral_without_precomp(j,k) = tau * xxvpar_coef_component(1)
!
                if(poly_order.ge.1) then
                    xxvpar_integral_without_precomp(j,k) = xxvpar_integral_without_precomp(j,k) +  &
                        & tau**2 * 0.5d0 * xxvpar_coef_component(2)
                endif
!
                if(poly_order.ge.2) then
                    xxvpar_integral_without_precomp(j,k) = xxvpar_integral_without_precomp(j,k) + &
                        & tau**3 / 3.d0 * xxvpar_coef_component(3)
                endif
!
                if(poly_order.ge.3) then
                    xxvpar_integral_without_precomp(j,k) = xxvpar_integral_without_precomp(j,k) + &
                        & tau**4 / 4.d0 * xxvpar_coef_component(4)
                endif
!
                if(poly_order.ge.4) then
                    xxvpar_integral_without_precomp(j,k) = xxvpar_integral_without_precomp(j,k) + &
                        & tau**5 / 5.d0 * xxvpar_coef_component(5)
                endif
!
            enddo
        enddo
!
        deallocate(xxvpar_coef_component)
!
    end function xxvpar_integral_without_precomp
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
        function bmod_func(z123)
!
            use tetra_physics_mod, only: tetra_physics
!
            implicit none
!
            double precision :: bmod_func
            double precision, dimension(3),intent(in) :: z123
!
            bmod_func = tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z123)
!        
        end function bmod_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        
        function phi_elec_func(z123)
!
            use tetra_physics_mod, only: tetra_physics
!
            implicit none
!
            double precision :: phi_elec_func
            double precision, dimension(3),intent(in) :: z123
!
            phi_elec_func = tetra_physics(ind_tetr)%Phi1 + sum(tetra_physics(ind_tetr)%gPhi*z123)
!        
        end function phi_elec_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
        function energy_tot_func(z)
!
            use tetra_physics_mod, only: tetra_physics,particle_mass,particle_charge
!
            implicit none
!
            double precision                :: energy_tot_func,vperp
            double precision, dimension(4)  :: z
!
            !Compute vperp
            vperp=sqrt(2.d0*abs(perpinv)*( tetra_physics(ind_tetr)%bmod1+sum(tetra_physics(ind_tetr)%gb*z(1:3)) ))
!
            energy_tot_func = particle_mass/2.d0*(vperp**2 + z(4)**2) + particle_charge*phi_elec_func(z(1:3))
!       
        end function energy_tot_func
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
        use tetra_physics_mod, only: tetra_physics
        use poly_without_precomp_mod, only: b
        use gorilla_settings_mod, only: i_precomp
!
        implicit none
!        
        integer, intent(in)                         :: iface
        double precision, dimension(4), intent(in)  :: z
        double precision                            :: normal_velocity_func
!
        if(i_precomp.eq.0) then
            normal_velocity_func = sum((tetra_physics_poly1(ind_tetr)%anorm_in_amat1_0(:,iface) + &
                        & perpinv * tetra_physics_poly1(ind_tetr)%anorm_in_amat1_1(:,iface)) * z) + &
                        & sum(tetra_physics(ind_tetr)%anorm(:,iface) * b(1:3))
        else
            normal_velocity_func = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,iface) + &
                                    & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,iface)) * z) + &
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
        vd_ExB = abs(clight/bmod0*tetra_physics(ind_tetr)%Er_mod)
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
        iface_new = face_init    !<-- Jonatan, Georg 07.07.2022: we think the face should be reseted because of analytic_approx belo
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
end module
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
        use pusher_tetra_poly_mod, only: dt_dtau_const,z_init,ind_tetr,analytic_integration_without_precomp, &
                                                & energy_tot_func
        use tetra_physics_mod, only: particle_mass,tetra_physics
!  
        implicit none
!
        integer, intent(in)                 :: poly_order,file_id_vpar_0,file_id_J_par,file_id_e_tot,n_skip_vpar_0
        double precision, intent(in)        :: t_pass,vpar_in,vpar_end
        logical, intent(in)                 :: boole_J_par,boole_poincare_vpar_0,boole_e_tot
        double precision                    :: tau,a44,b4,tau_part1
        double precision, dimension(4)      :: z
        double precision, dimension(3)      :: x
!
        nskip = n_skip_vpar_0
!
        !Convert t_pass to tau
        tau = t_pass/dt_dtau_const
!
        !Select matrix and vector elements that are use for v_parallel integration
        a44 = amat(4,4)
        b4 = b(4)
!
        !Trace banana tips
        if((vpar_end.gt.0.d0).and.(vpar_in.lt.0.d0)) then
            !Find exact root of vpar(tau)  --> tau_part1 is the "time" inside the tetrahedron until $v=\parallel$ = 0
            tau_part1 = tau_vpar_root(poly_order,a44,b4,vpar_in)
            
            if(tau_part1.gt.tau) print *, 'Error'
!
            !Integrate par_adiab_inv exactly until root of vpar(tau)
            par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau_part1,vpar_in)*dt_dtau_const  
        
            if(counter_banana_mappings.gt.1) then
                if(counter_banana_mappings/nskip*nskip.eq.counter_banana_mappings) then

!                    
                    !Poloidal projection of Poincaré sections at $v=\parallel$ = 0
                    z = z_init
                    call analytic_integration_without_precomp(poly_order,z,tau_part1)
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
                            write(file_id_e_tot,*) counter_banana_mappings,energy_tot_func(z)
                        endif
                    !$omp end critical
                    !print *, 'banana bounces:', counter_banana_mappings
                endif
            endif
            counter_banana_mappings = counter_banana_mappings + 1
!
            !Start to integrate par_adiab_inv for new bounce period
            par_adiab_inv = 0.d0
            par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau-tau_part1,0.d0)*dt_dtau_const 
        else    
            !Compute parallel adiabatic invariant as a function of time
            par_adiab_inv = par_adiab_inv + par_adiab_tau(poly_order,a44,b4,tau,vpar_in)*dt_dtau_const    
        endif     
!  
    end subroutine par_adiab_inv_tetra_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
    function par_adiab_tau(poly_order,a44,b4,tau,vpar_in)
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
        !f(tau)  = a + b*tau + c/2 * tau^2 + d/6 * tau^3 + e/24 * tau^5
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
end module
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module vperp2_integral_mod
!
!WARNING: This module and its functions and routines were never tested and only written conceptually.
!
    implicit none
!
    private
!
    public :: vperp2_integral_tetra
!
    contains
!    
    function vperp2_integral_tetra(poly_order,t_pass)
        !WARNING: This subroutine CAN ONLY BE CALLED directly after pusher, when modules still contain
        !           the quantities from last pushing
!        
        use pusher_tetra_poly_mod, only: ind_tetr,perpinv,dt_dtau_const,bmod0
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!        
        integer                         :: poly_order
        double precision                :: vperp2_integral_tetra,tau,t_pass
!
        !Convert t_pass to tau
        tau = t_pass/dt_dtau_const
!
        !Integral of vper^2(tau) over dwell time-orbit parameter
        vperp2_integral_tetra = -2.d0*perpinv*( bmod0*tau &
                            & + sum(tetra_physics(ind_tetr)%gB * x_integral_tetra(poly_order,tau)) )
!
        !Integral of vper^2(time)
        vperp2_integral_tetra = vperp2_integral_tetra*dt_dtau_const
! 
    end function vperp2_integral_tetra
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function x_integral_tetra(poly_order,tau)
!
        use poly_without_precomp_mod
        use pusher_tetra_poly_mod, only: z_init
!
        implicit none
!
        integer, parameter                  :: n_dim=3
        integer                             :: poly_order
        double precision, dimension(n_dim)  :: x_integral_tetra
        double precision                    :: tau,tau2_half,tau3_sixth,tau4_twelvth,tau5_hundredtwentyth
!
            if(poly_order.ge.1) then
                tau2_half = tau**2*0.5d0
                x_integral_tetra(1:n_dim) = z_init(1:n_dim)*tau + tau2_half*(b(1:n_dim)+amat_in_z(1:n_dim))
            endif
!
            if(poly_order.ge.2) then
                tau3_sixth = tau**3/6.d0
                x_integral_tetra(1:n_dim) = x_integral_tetra(1:n_dim) + tau3_sixth*(amat_in_b(1:n_dim) &
                                          & + amat2_in_z(1:n_dim))
            endif
!
            if(poly_order.ge.3) then
                tau4_twelvth = tau**4/12.d0
                x_integral_tetra(1:n_dim) = x_integral_tetra(1:n_dim) + tau4_twelvth*(amat2_in_b(1:n_dim) &
                                          & + amat3_in_z(1:n_dim))
            endif
!
            if(poly_order.ge.4) then
                tau5_hundredtwentyth = tau**5/120.d0
                x_integral_tetra(1:n_dim) = x_integral_tetra(1:n_dim) + tau5_hundredtwentyth *(amat3_in_b(1:n_dim) &
                                          & + amat4_in_z(1:n_dim))  
            endif
!
    end function x_integral_tetra
!
end module    
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
