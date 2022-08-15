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
    !$OMP THREADPRIVATE(ind_tetr,iface_init,perpinv,perpinv2,dt_dtau_const,bmod0,t_remain,x_init,  &
    !$OMP& z_init,k1,k3,vmod0)
!
    public :: pusher_tetra_poly,initialize_const_motion_poly, &
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
                                               & boole_t_finished,iper_phi,t_hamiltonian,gyrophase)
!
            use tetra_physics_mod, only: tetra_physics,particle_charge,particle_mass
            use gorilla_diag_mod,only: diag_pusher_tetry_poly
            use pusher_tetra_func_mod, only: pusher_handover2neighbour
            use gorilla_settings_mod, only: i_precomp, boole_guess, boole_time_Hamiltonian, boole_gyrophase
!
            implicit none
!
            integer                                             :: poly_order
            integer,intent(inout)                               :: ind_tetr_inout,iface
            integer,intent(out)                                 :: iper_phi
            double precision, dimension(3), intent(inout)       :: x
            double precision, intent(inout)                     :: vpar
            double precision, dimension(3), intent(out)         :: z_save
            double precision, intent(in)                        :: t_remain_in
            double precision, intent(out)                       :: t_pass
            logical, intent(out)                                :: boole_t_finished
            double precision, intent(out),optional              :: t_hamiltonian,gyrophase
            logical, dimension(4)                               :: boole_faces
            integer                                             :: i,j,k
            double precision, dimension(4)                      :: z,operator_b_in_b,z_dummy
            integer                                             :: iface_new,i_scaling
            double precision                                    :: tau,vperp2,tau_save,tau_est,tau_max
            logical                                             :: boole_analytical_approx,boole_face_correct,boole_vnorm_correction
            logical                                             :: boole_trouble_shooting
            double precision,dimension(4,4)                     :: operator_b,operator_z_init
!         
            call initialize_pusher_tetra_poly(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
            !Initial computation values
            z = z_init
!
            !Hamiltonian time (Only Delta-value within one pushing is used for output)
            if(boole_time_Hamiltonian) then
                if(present(t_hamiltonian)) then
                    t_hamiltonian = 0.d0
                else
                    print *, 'Error: Hamiltonian-time dynamics computation is on, but NOT used. ', &
                             & 'Please, set boole_time_Hamiltonian to .false.'
                    stop
                endif
            endif
!
            !Gyrophase (Only Delta-value within one pushing is used for output)
            if(boole_gyrophase) then
                if(present(gyrophase)) then
                    gyrophase = 0.d0
                else
                    print *, 'Error: Gyrophyse computation is on, but NOT used. ', &
                             & 'Please, set boole_gyrophase to .false.'
                    stop
                endif
            endif
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
            !Initialize vnorm correction (Particle turns exactly at the face)
            boole_vnorm_correction = .false.
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
            iface_new = iface
!
            !boole_faces ... Boolean array for switching on/off root computation for individual face
            !Initialize boole_faces, such that roots for all 4 faces are computed
            boole_faces = .true.
!
            !Low order polynomial face prediction
            if(boole_guess) then
!
                !Analytical calculation of orbit parameter to guess exit face
                call analytic_approx(2,i_precomp,boole_faces, &
                                    & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
                !Define boundary for tau_max
                tau_max = tau*eps_tau
                                    
if(diag_pusher_tetry_poly) print *, 'boole_analytical_approx',boole_analytical_approx
!               
                !Higher order polynomial root is computed only for quadratically predicted face
                if(poly_order.gt.2) then
!    
                    if(boole_analytical_approx) then
                        boole_faces = .false.                !Disable all 4 faces 
                        boole_faces(iface_new) = .true.      !Enable guessed face
                        iface_new = iface                    !Set iface_new to entering iface
!                    
                        !Analytical calculation of orbit parameter to pass tetrahdron
                        call analytic_approx(poly_order,i_precomp,boole_faces, &
                                            & i_scaling,z,iface_new,tau,boole_analytical_approx)
! 
if(diag_pusher_tetry_poly) print *, 'tau',tau
if(diag_pusher_tetry_poly) print *, 'iface guess', iface_new
!
                    else !No analytical approximation exists
                        iface_new = iface   
! 
                        !Analytical calculation of orbit parameter to pass tetrahdron
                        call analytic_approx(poly_order,i_precomp,boole_faces, &
                                            & i_scaling,z,iface_new,tau,boole_analytical_approx)                  
                    endif
!                
                endif
!    
            !No face guessing, instead roots are solved for all 4 faces
            else
                !Analytical calculation of orbit parameter for safety boundary
                call analytic_approx(2,i_precomp,boole_faces, &
                                    & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
                !Define boundary for tau_max
                tau_max = tau*eps_tau  
                
                if(poly_order.gt.2) then
                    iface_new = iface  
!                
                    !Analytical calculation of orbit parameter to pass tetrahdron
                    call analytic_approx(poly_order,i_precomp,boole_faces, &
                                        & i_scaling,z,iface_new,tau,boole_analytical_approx)   
                endif
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
                if(boole_time_Hamiltonian) then
                    if(boole_gyrophase) then
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                    else
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                    endif !gyrophase
                else
                    call analytic_integration(poly_order,i_precomp,z,tau)
                endif
!
                !Validation loop ('3-planes'-control)
                three_planes_loop: do j=1,3                             !Just consider the planes without the "exit-plane"
                    k=modulo(iface_new+j-1,4)+1
                    if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                        boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: three planes'
                    endif        
                enddo three_planes_loop
!
                !Accuracy on face
                if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: distance'
                endif
!                
                !Validation for vnorm
                if(normal_velocity_func(z,iface_new).gt.0.d0) then
if(diag_pusher_tetry_poly) print *, 'Error in predicted integrator: normal velocity'
                    if(poly_order.gt.2) then
                        boole_face_correct = .false.
!                        
                    else !Error handling (only for quadratic solver): Particle turns 'exactly' on the face
!                   
                        boole_vnorm_correction = .true.
!
                        !Particle orbit turns exactly at the face. Find next exit point of same tetrahedron.
                        boole_faces = .true.
                        tau_save = tau
!                    
                        !Analytical calculation of orbit parameter for safety boundary
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
                        if(boole_time_Hamiltonian) then
                            if(boole_gyrophase) then
                                call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                            else
                                call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                            endif !gyrophase
                        else
                            call analytic_integration(poly_order,i_precomp,z,tau)
                        endif
!
                        !Integration in two parts (Therefore time must be added)
                        tau = tau + tau_save
!
                        !Validation loop ('3-planes'-control)
                        do j=1,3                             !Just consider the planes without the "exit-plane"
                            k=modulo(iface_new+j-1,4)+1
                            if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                                boole_face_correct = .false.
                            endif        
                        enddo
!                    
                        !Accuracy on face
                        if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
                            boole_face_correct = .false.
                        endif
!                        
                        !Validation for vnorm
                        if(normal_velocity_func(z,iface_new).gt.0.d0) then
                            boole_face_correct = .false.
                        endif
!       
                    endif !vnorm is positive
                endif
!                
                !Higher order polynomial result for tau is out of safety boundary
                if(poly_order.gt.2) then
                    if(tau.gt.tau_max) then
if(diag_pusher_tetry_poly)                        print *, 'Error: Tau is out of safety boundary (1)'
                        boole_face_correct = .false.
                    endif
                endif               
!
            endif !boole face correct
!
if(diag_pusher_tetry_poly) print *, 'boole_face_correct',boole_face_correct

            !If error is detected, repeat orbit pushing but without low order polynomial guessing
            if(.not.boole_face_correct) then
!
!                if(poly_order.eq.2) then
!                    print *, 'Error: Wrong orbit pushing in quadratic order'
!                    stop
!                endif
!
                boole_faces = .true.    !Compute root for all 4 faces
                iface_new = iface 
                z = z_init
!
                !Hamiltonian time
                if(boole_time_Hamiltonian) t_hamiltonian = 0.d0
!
                !Gyrophase
                if(boole_gyrophase) gyrophase = 0.d0
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
                if(boole_time_Hamiltonian) then
                    if(boole_gyrophase) then
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                    else
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                    endif ! gyrophase
                else
                    call analytic_integration(poly_order,i_precomp,z,tau)
                endif
!
                !Higher order polynomial result for tau is out of safety boundary
                if(poly_order.gt.2) then
                    if(tau.gt.tau_max) then
!                       print *, 'Error: Tau is out of safety boundary (2)'
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (1)'
                        if(boole_time_Hamiltonian) then
                            if(boole_gyrophase) then
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian,gyrophase)
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian)
                            endif !gyrophase
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                        endif
!                       
                        if(.not.boole_trouble_shooting) then
                            print *, 'Error: Trouble shooting failed. Remove particle.'
                            ind_tetr_inout = -1
                            iface = -1
                            return
                        endif
!
                    endif
                endif
!
                !Validation loop ('3-planes'-control)
                do j=1,3                             !Just consider the planes without the "exit-plane"
                    k=modulo(iface_new+j-1,4)+1

                    if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (2)'
                        if(boole_time_Hamiltonian) then
                            if(boole_gyrophase) then
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian,gyrophase)
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian)
                            endif !gyrophase
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                        endif
!                       
                        if(.not.boole_trouble_shooting) then
                            print *, 'Error: Trouble shooting failed. Remove particle.'
                            ind_tetr_inout = -1
                            iface = -1
                            return
                        endif
!
                    endif        
                enddo
!
                !Accuracy on face
                if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
!                   print *, 'Error in accuracy'
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (3)'
                    if(boole_time_Hamiltonian) then
                        if(boole_gyrophase) then
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting,t_hamiltonian,gyrophase)
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting,t_hamiltonian)
                        endif !gyrophase
                    else
                        call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                    endif
!                       
                    if(.not.boole_trouble_shooting) then
                        print *, 'Error: Trouble shooting failed. Remove particle.'
                        ind_tetr_inout = -1
                        iface = -1
                        return
                    endif
!
                endif
!
                !Validation for vnorm
                if(normal_velocity_func(z,iface_new).gt.0.d0) then
!                   
                    boole_vnorm_correction = .true.
!
                    !Particle orbit turns exactly at the face. Find next exit point of same tetrahedron.
                    boole_faces = .true.
                    tau_save = tau
!                    
                    !Analytical calculation of orbit parameter for safety boundary
                    call analytic_approx(2,i_precomp,boole_faces, &
                                        & i_scaling,z,iface_new,tau,boole_analytical_approx)
!
                    !Define boundary for tau_max
                    tau_max = tau*eps_tau     
!                    
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
                    if(boole_time_Hamiltonian) then
                        if(boole_gyrophase) then
                            call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                        else
                            call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                        endif !gyrophase
                    else
                        call analytic_integration(poly_order,i_precomp,z,tau)
                    endif
!
                    !Higher order polynomial result for tau is out of safety boundary
                    if(poly_order.gt.2) then
                        if(tau.gt.tau_max) then
!                             print *, 'Error: Tau is out of safety boundary (3)'
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (4)'
                            if(boole_time_Hamiltonian) then
                                if(boole_gyrophase) then
                                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                        & boole_trouble_shooting,t_hamiltonian,gyrophase)
                                else
                                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                        & boole_trouble_shooting,t_hamiltonian)
                                endif !gyrophase
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                            endif
!                       
                            if(.not.boole_trouble_shooting) then
                                print *, 'Error: Trouble shooting failed. Remove particle.'
                                ind_tetr_inout = -1
                                iface = -1
                                return
                            endif
!
                        endif
                    endif  
!
                    !Integration in two parts (Therefore time must be added)
                    tau = tau + tau_save
!
                    !Validation loop ('3-planes'-control)
                    do j=1,3                             !Just consider the planes without the "exit-plane"
                        k=modulo(iface_new+j-1,4)+1
                        if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (5)'
                            if(boole_time_Hamiltonian) then
                                if(boole_gyrophase) then
                                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                        & boole_trouble_shooting,t_hamiltonian,gyrophase)
                                else
                                    call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                        & boole_trouble_shooting,t_hamiltonian)
                                endif !gyrophase
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                            endif
!                       
                            if(.not.boole_trouble_shooting) then
                                print *, 'Error: Trouble shooting failed. Remove particle.'
                                ind_tetr_inout = -1
                                iface = -1
                                return
                            endif
!
                        endif        
                    enddo
!                    
                    !Validation for vnorm
                    if(normal_velocity_func(z,iface_new).gt.0.d0) then
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (6)'
                        if(boole_time_Hamiltonian) then
                            if(boole_gyrophase) then
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian,gyrophase)
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian)
                            endif !gyrophase
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                        endif
!                       
                        if(.not.boole_trouble_shooting) then
                            print *, 'Error: Trouble shooting failed. Remove particle.'
                            ind_tetr_inout = -1
                            iface = -1
                            return
                        endif
!
                    endif
!
                    !Accuracy on face
                    if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
if(diag_pusher_tetry_poly) print *, 'Trouble shooting (7)'
                        if(boole_time_Hamiltonian) then
                            if(boole_gyrophase) then
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian,gyrophase)
                            else
                                call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                    & boole_trouble_shooting,t_hamiltonian)
                            endif !gyrophase
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                        endif
!                       
                        if(.not.boole_trouble_shooting) then
                            print *, 'Error: Trouble shooting failed. Remove particle.'
                            ind_tetr_inout = -1
                            iface = -1
                            return
                        endif
!
                    endif
                endif !Normal velocity is positive at exit point

            endif !.not. boole_face_correct
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

            !Particle stops inside the tetrahedron
            if(t_pass.ge.t_remain) then
!
                !Set z back to z_init
                z = z_init
!
                !Hamiltonian time
                if(boole_time_Hamiltonian) t_hamiltonian = 0.d0
!
                !Gyrophase
                if(boole_gyrophase) gyrophase = 0.d0
!
                !If vnorm was corrected, coefficients need to be computed again.
                if(boole_vnorm_correction) then
                    iface_new = iface_init
                    boole_faces = .true.
!
                    call analytic_approx(poly_order,i_precomp,boole_faces, &
                                         & i_scaling,z,iface_new,tau,boole_analytical_approx)
                endif
!
                !Compute orbit parameter tau from t_remain
                tau = t_remain/dt_dtau_const
if(diag_pusher_tetry_poly) print *, 'tau until t finished',tau
!
                !Integrate trajectory analytically from start until t_remain
                if(boole_time_Hamiltonian) then
                    if(boole_gyrophase) then
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                    else
                        call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                    endif !gyrophase
                else
                    call analytic_integration(poly_order,i_precomp,z,tau)
                endif
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

                if(boole_face_correct) then
                    boole_t_finished = .true.
                    
                    z_save = z(1:3)
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    vpar=z(4)
                    t_pass = t_remain
                    
                else    !Time step finishes while particle is outside the tetrahedron (Wrong face was predicted, but logically true)
                    !Compute correct root by honestly computing higher order root for all four faces
if(diag_pusher_tetry_poly) then
    print *, 'Error: Particle is outside the tetrahedron when time is finished.'
    do i = 1,4
        print *,i, 'norm', normal_distance_func(z(1:3),i) 
    enddo
endif
                    if(boole_time_Hamiltonian) then
                        if(boole_gyrophase) then
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting,t_hamiltonian,gyrophase)
                        else
                            call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new, &
                                                                & boole_trouble_shooting,t_hamiltonian)
                        endif !gyrophase
                    else
                        call trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting)
                    endif
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
        end subroutine pusher_tetra_poly
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
        subroutine analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
!
            implicit none
!
            integer, intent(in)                             :: poly_order,i_precomp
            double precision, intent(in)                    :: tau
            double precision, dimension(4),intent(inout)    :: z
            double precision, intent(inout),optional        :: t_hamiltonian,gyrophase
!
            !Integrate trajectory analytically
            select case(i_precomp)
                case(0)
                    call analytic_integration_without_precomp(poly_order,z,tau,t_hamiltonian,gyrophase)
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
        subroutine analytic_integration_without_precomp(poly_order,z,tau,t_hamiltonian,gyrophase)
!
            use poly_without_precomp_mod
            use gorilla_settings_mod, only: boole_time_Hamiltonian, boole_gyrophase
            use tetra_physics_mod, only: hamiltonian_time,cm_over_e,tetra_physics
!
            implicit none
!
            integer, intent(in)                          :: poly_order
            double precision, intent(in)                 :: tau
            double precision, dimension(4),intent(inout) :: z
            double precision, intent(inout),optional     :: t_hamiltonian,gyrophase
            double precision                              :: delta_t_hamiltonian
            double precision                             :: tau2_half,tau3_sixth,tau4_twentyfourth
!
            double precision, dimension(3)                :: x0
            double precision                              :: vpar0
            double precision, dimension(:,:), allocatable :: x_coef,x_vpar_coef
            double precision, dimension(:), allocatable   :: vpar_coef,res_poly_coef
!
            !Save coordinates before analytical integration
            x0 = z(1:3)
            vpar0 = z(4)
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
                t_hamiltonian = t_hamiltonian + delta_t_hamiltonian
!
                !Optional computation of gyrophase
                if(boole_gyrophase) then
                    gyrophase = gyrophase - ( &
!
!                        !Version with dt_dtau_const
!                        & 1.d0/cm_over_e * hamiltonian_time(ind_tetr)%dt_dtau_const_save * &
!                        & ( tetra_physics(ind_tetr)%bmod1 * tau + &
!                        & sum( tetra_physics(ind_tetr)%gb * x_integral_without_precomp(poly_order,tau,x_coef) ) ) &
!
!
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
        end subroutine analytic_integration_without_precomp
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine analytic_integration_without_precomp_stepwise(poly_order,z,tau)
!
            use poly_without_precomp_mod, only: amat,amat2,amat3,amat4,b
!
            implicit none
!
            integer, intent(in)                          :: poly_order
            double precision, intent(in)                 :: tau
            double precision, dimension(4),intent(inout) :: z
            double precision                             :: tau2_half,tau3_sixth,tau4_twentyfourth
            double precision, dimension(4)      :: amat_in_z,amat2_in_z,amat3_in_z,amat_in_b,amat2_in_b
            double precision, dimension(4)      :: amat4_in_z,amat3_in_b, z_save
!
            z_save = z
!
            if(poly_order.ge.1) then
                amat_in_z = matmul(amat,z_save)
                z = z + tau*(b+amat_in_z)
            endif
!
            if(poly_order.ge.2) then
                amat2_in_z = matmul(amat2,z_save)
                amat_in_b = matmul(amat,b)
                tau2_half = tau**2*0.5d0
                z = z + tau2_half*(amat_in_b + amat2_in_z)
            endif
!
            if(poly_order.ge.3) then
                amat3_in_z = matmul(amat3,z_save)
                amat2_in_b = matmul(amat2,b)
                tau3_sixth = tau**3/6.d0
                z = z + tau3_sixth*(amat2_in_b + amat3_in_z)
            endif

            if(poly_order.ge.4) then
                amat4_in_z = matmul(amat4,z_save)
                amat3_in_b = matmul(amat3,b)
                tau4_twentyfourth = tau**4/24.d0
                z = z + tau4_twentyfourth *(amat3_in_b + amat4_in_z)
            endif
!
        end subroutine analytic_integration_without_precomp_stepwise
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
    subroutine trouble_shooting_polynomial_solver(poly_order,i_precomp,z,tau,iface_new,boole_trouble_shooting,t_hamiltonian, &
        & gyrophase)
!
        use gorilla_diag_mod, only: diag_pusher_tetry_poly
        use gorilla_settings_mod, only: boole_time_Hamiltonian,boole_gyrophase
!
        implicit none
!
        integer, intent(in)                                 :: poly_order,i_precomp
        integer, intent(out)                                :: iface_new
        logical, intent(out)                                :: boole_trouble_shooting
        double precision, intent(out),optional              :: t_hamiltonian,gyrophase
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
!
!        !Trouble shooting is only developed for fourth order
!        if(poly_order.eq.2) then
!            print *, 'Error: Trouble shooting is only developed for third and fourth order'
!            stop
!        endif
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
!
            !Hamiltonian time
            if(boole_time_Hamiltonian) t_hamiltonian = 0.d0
!
            !Gyrophase
            if(boole_gyrophase) gyrophase = 0.d0
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
            if(boole_time_Hamiltonian) then
                if(boole_gyrophase) then
                    call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                else
                    call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                endif !gyrophase
            else
                call analytic_integration(poly_order,i_precomp,z,tau)
            endif
!
            !Control section
            boole_face_correct = .true.
!        
            !Validation loop ('3-planes'-control)
            three_planes_loop: do j=1,3                             !Just consider the planes without the "exit-plane"
                k=modulo(iface_new+j-1,4)+1
                if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: three planes'
                endif        
            enddo three_planes_loop
!
            !Accuracy on face
            if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: distance'
            endif
!
            !Validation for vnorm
            if(normal_velocity_func(z,iface_new).gt.0.d0) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: normal velocity'
            endif
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
!
            !Hamiltonian time
            if(boole_time_Hamiltonian) t_hamiltonian = 0.d0
!
            !Gyrophase
            if(boole_gyrophase) gyrophase = 0.d0
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
            if(boole_time_Hamiltonian) then
                if(boole_gyrophase) then
                    call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian,gyrophase)
                else
                    call analytic_integration(poly_order,i_precomp,z,tau,t_hamiltonian)
                endif !gyrophase
            else
                call analytic_integration(poly_order,i_precomp,z,tau)
            endif
            
            !Control section
            boole_face_correct = .true.
!        
            !Validation loop ('3-planes'-control)
            do j=1,3                             !Just consider the planes without the "exit-plane"
                k=modulo(iface_new+j-1,4)+1
                if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                    boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: three planes'
                endif        
            enddo
!
            !Accuracy on face
            if(abs(normal_distance_func(z(1:3),iface_new)).gt.1.d-11) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: distance'
            endif
!            
            !Validation for vnorm
            if(normal_velocity_func(z,iface_new).gt.0.d0) then
                boole_face_correct = .false.
if(diag_pusher_tetry_poly) print *, 'Trouble shooting: normal velocity'
            endif           
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
