!
module pusher_tetra_rk_mod
!
    private
!
    integer                          :: iface_init
    integer,public,protected         :: ind_tetr
    double precision                 :: B0,perpinv,perpinv2,vmod_init
    double precision,public,protected :: spamat,dt_dtau_const
    double precision                 :: dist1,dist_min,dist_max,t_remain
    double precision, dimension(3)   :: gradB,x_init
    double precision, dimension(3),public,protected   :: Bvec
    double precision, dimension(4),public,protected   :: b,z_init
    double precision, dimension(3,3),public,protected :: amat
    double precision, dimension(3,4) :: anorm
    double precision,parameter                 :: eps_distmin = 1.d-10
    double precision,parameter                 :: eps_distmax = 10.d0
    double precision                 :: dtau_ref,dtau_max,dtau_quad
    integer, parameter               :: kiter=48  
    double precision,parameter      :: eps_modulation = 0.1d0
    double precision,parameter      :: eps_dtau_max = 10.d0
    double precision,parameter      :: eps_dtau_quad = 1.5d0
    double precision                 :: k1, k3
!   
    !$OMP THREADPRIVATE(ind_tetr,spamat,B0,perpinv,perpinv2,iface_init,dt_dtau_const,k1,k3, &
    !$OMP& dist1,dist_min,dist_max,Bvec,gradB,x_init,b,z_init,amat,anorm,t_remain,vmod_init,dtau_ref,dtau_max,dtau_quad)
!    
    public :: pusher_tetra_rk,find_tetra,initialize_const_motion_rk,energy_tot_func
!    
    contains
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine initialize_const_motion_rk(perpinv_in,perpinv2_in)
!
        implicit none
!            
        double precision, intent(in)    :: perpinv_in,perpinv2_in
!
        perpinv = perpinv_in
        perpinv2 = perpinv2_in
!            
    end subroutine initialize_const_motion_rk
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine initialize_pusher_tetra_rk_mod(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
        use tetra_physics_mod, only: tetra_physics, cm_over_e,dt_dtau,coord_system
        use tetra_grid_settings_mod, only: grid_size
        use constants, only : clight,eps,pi
        use gorilla_settings_mod, only: boole_dt_dtau
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!
        implicit none
!
        integer, intent(in)                        :: ind_tetr_inout,iface
        double precision, intent(in)               :: vpar
        double precision, dimension(3), intent(in) :: x
        double precision                           :: bmod,vperp2,vpar2,phi_elec,t_remain_in,tetra_dist_ref,vd_ExB
        logical                                    :: boole_vd_ExB,boole_vperp
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
        iface_init = iface
!
        !Tetrahedron constants
        B0 = tetra_physics(ind_tetr)%bmod1
        gradB = tetra_physics(ind_tetr)%gb
        anorm = tetra_physics(ind_tetr)%anorm
        dt_dtau_const = tetra_physics(ind_tetr)%dt_dtau_const
!    
        !Module of B at the entry point of the particle
        bmod=bmod_func(z_init(1:3))
!
        !Phi at the entry point of the particle
        phi_elec=tetra_physics(ind_tetr)%Phi1+sum(tetra_physics(ind_tetr)%gPhi*z_init(1:3))
!
        !Auxiliary quantities
        vperp2 = -2.d0*perpinv*bmod
        vpar2 = vpar**2
!
        !velocity module
        vmod_init = sqrt(vpar2+vperp2)
!
        !ODE coefficients
        b(1:3)=( tetra_physics(ind_tetr)%curlh*(vperp2+vpar2+2.d0*perpinv*B0) &
            + perpinv*tetra_physics(ind_tetr)%gBxh1 )*cm_over_e &
            - clight*(2.d0*(tetra_physics(ind_tetr)%Phi1-phi_elec)*tetra_physics(ind_tetr)%curlh &
            + tetra_physics(ind_tetr)%gPhixh1)
!
        b(4)=perpinv*tetra_physics(ind_tetr)%gBxcurlA-clight/cm_over_e*tetra_physics(ind_tetr)%gPhixcurlA
!
!         k1 = vperp2+vpar2+2.d0*perpinv*B0
!         k3 = tetra_physics(ind_tetr)%Phi1-phi_elec
!
        Bvec=tetra_physics(ind_tetr)%curlA    !B-Vector
        amat=perpinv*cm_over_e*tetra_physics(ind_tetr)%alpmat &
            - clight* tetra_physics(ind_tetr)%betmat    ! a-Matrix (elements 1:3)
!                
        spamat=perpinv*cm_over_e*tetra_physics(ind_tetr)%spalpmat &
            - clight* tetra_physics(ind_tetr)%spbetmat    !tr(a-Matrix)
!
        !Reference distances
        dist1=-tetra_physics(ind_tetr)%dist_ref     ! Normal distance of first vertex to first plane
        dist_min=eps_distmin*abs(dist1)                     ! relative value for convergence criterion, because of the different volumes of tetrahedra
        dist_max=eps_distmax*abs(dist1)             ! relative maximum distance (to avoid huge integration distances for RK4 accuracy)
!
        !Tetrahdron reference distance
        tetra_dist_ref = abs(tetra_physics(ind_tetr)%tetra_dist_ref)
!        
        !ExB drift velocity
        vd_ExB = abs(clight/bmod*tetra_physics(ind_tetr)%Er_mod)
        boole_vd_ExB = .true.
        if(vd_ExB.eq.0.d0) boole_vd_ExB = .false.
!
        boole_vperp = .true.
        if(vperp2.eq.0.d0) boole_vperp = .false.
!
        !Reference time to pass the tetrahedron
        !Actually dt_ref (Save on additional variable dt_ref):
        if(boole_vd_ExB.and.boole_vperp) then
            dtau_ref = minval([abs(tetra_dist_ref/z_init(4)), &
            & sqrt(tetra_dist_ref*vmod_init*tetra_physics(ind_tetr)%R1 / &
            & (vperp2*grid_size(2)*eps_modulation) ), &
            & tetra_dist_ref/vd_ExB],1)         
        elseif(boole_vd_ExB.and.(.not.boole_vperp)) then
            dtau_ref = minval([abs(tetra_dist_ref/z_init(4)), &
            & tetra_dist_ref/vd_ExB],1)        
        elseif((.not.boole_vd_ExB).and.boole_vperp) then
            dtau_ref = minval([abs(tetra_dist_ref/z_init(4)), &
            & sqrt(tetra_dist_ref*vmod_init*tetra_physics(ind_tetr)%R1 / &
!
            & (vperp2*grid_size(2)*eps_modulation) )],1)          
        else
            dtau_ref = abs(tetra_dist_ref/z_init(4))            
        endif
!
        !Transform t to tau            
        if(boole_dt_dtau) then
            dtau_ref = abs(dtau_ref/dt_dtau_const)
        else
            dtau_ref = abs(dtau_ref/dt_dtau(ind_tetr,x_init,x_init))
        endif   
!        
        !Maximum integration orbit parameter tau that is allowed
        dtau_max = eps_dtau_max*dtau_ref
        dtau_quad = eps_dtau_quad*dtau_ref !Reference length, if quadratic analytical approximation should start from initial values or not
!
    end subroutine initialize_pusher_tetra_rk_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine pusher_tetra_rk(ind_tetr_inout,iface,x,vpar,z_final,t_remain_in,t_pass,boole_t_finished,iper_phi)
!
        use tetra_grid_mod, only: tetra_grid, verts_sthetaphi
        use tetra_physics_mod, only: tetra_physics, cm_over_e,coord_system,dt_dtau,isinside
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
        use constants, only : clight,amp,echarge,eps,pi,ev2erg
        use supporting_functions_mod, only: calc_eigval
        use odeint_mod, only: odeint_allroutines
        use gorilla_settings_mod, only: boole_pusher_ode45
!
        implicit none
!
        logical,intent(out) :: boole_t_finished
        integer,intent(inout) :: ind_tetr_inout,iface
        integer,intent(out) :: iper_phi
        integer :: iper_theta
        double precision,intent(inout) :: vpar
        double precision, intent(in) :: t_remain_in
        double precision, intent(out) :: t_pass
        double precision,dimension(3),intent(inout) :: x
!         double precision,dimension(3),intent(out),optional :: x_final
!
!
        logical :: finish,llod_outside_logical,boole_particle_converged,boole_particle_turned
        logical :: boole_particle_turned_tangential,boole_acoef0,boole_newton1_tangential,boole_newton2_tangential
        logical, dimension(4) :: allowed_faces,boole_vec_acoef0
        integer :: iface_new,iface_new1,iface_bisec,i_outside_plane,i,j,k,l,m,n,o
        integer :: llod,extr_tang_err,eps_err,llod_inside_init,llod_iface_outside,llod_iface_outside_counter,fp_iface_outside
        integer, dimension(:), allocatable :: p
        double precision :: tau,tau2,bmod,bmod_beg,vperp2,v2,dtau,dtau_bisec,phi_elec,t_pass_inout,p_phi,tau_save
        double precision :: discr,dtau_new,dist,dist_new,vnorm,tau_bisec
        double precision :: dummy,dtau_save
        double precision, dimension(3)   :: x_beg, z_final
        double precision, dimension(4)   :: z,acoef,bcoef,ccoef,zbeg,zbeg2,z_bisec,z_llod_init,z_rk4, z_save
!
        double precision, dimension(4)   :: dzdtau,dzdtau_save
!
        !Testing variables
        double precision, dimension(3)   :: x_save
        double precision, dimension(3) :: delta_x,x_analytic
        !double precision :: dr_save
        double precision :: tau_ode45,eigval_lambda,trace_full_a, vmod,vdrift,vdrift2
        double precision, dimension(4) :: z_analytic,z_term5,z_term,normal_distances
        logical                      :: boole_quad_approx,boole_newton_converged,boole_converged
        logical                      :: boole_last_line_defense_converged
        logical                      :: boole_accuracy_ode45,boole_final_processing_converged
!  
    ! ----------------------------------------------------------------------------------!

        !Initialize starting values and 'working' constants in module for tetrahedron
        call initialize_pusher_tetra_rk_mod(ind_tetr_inout,x,iface,vpar,t_remain_in)
!
!
if(diag_pusher_tetra_rk) call print_starting_condition()
!
        !Initial computation values
        z = z_init
        tau = 0.d0
        iface_new = iface
!
        !Initialize boole_t_finished
        boole_t_finished = .false.
!
        !Initialize accuracy for integration step
        boole_accuracy_ode45 = boole_pusher_ode45
!
        !Set all allowed faces to 'possibly' true
        allowed_faces = .true.
!
!        ! Important testing feature: Is the particle outside the tetrahedron regarding more than one plane?
!        k = 0
!        normal_distances = normal_distances_func(z(1:3))
!if(diag_pusher_tetra_rk) print *,'normal_distances at beginning',normal_distances
!        call rk4_step(z,0.d0,dzdtau)
!        do i = 1,4
!            if(normal_distances(i).lt.0.d0) k = k+1
!        enddo
!
!        if(k.gt.1) then
!            print *, 'Error: The particle is outside the tetrahedron regarding more than one plane.'
!            print *, 'normal_distances',normal_distances
!            call print_starting_condition()
!            stop
!        endif


!
!if(iface_init.ne.0) then
!    if(abs(normal_distance_func(z(1:3),iface_init)).gt.(abs(dist_min)*1.d1)) then
!        print *, 'Error - Wrong starting condition: particle is initially not converged at iface'
!        call print_starting_condition()
!        print *, 'normal_distance(iface_init)', abs(normal_distance_func(z(1:3),iface_init))
!        print *, 'dist_min',dist_min
!        stop
!    endif
!endif

!

        !First quadratic analytical approximation
        call quad_analytic_approx(z,allowed_faces,iface_new,dtau,boole_quad_approx)
!
        !Analytical approximation exists
        if(boole_quad_approx) then
if(diag_pusher_tetra_rk) print *, 'Quadratic approximation exists'
!           
            !First RK4 step with analytical approximation
            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
            tau = tau + dtau
!
        !No analytical approximation exists
        else  
if(diag_pusher_tetra_rk) print *, 'no analytical approximation exists'
!
            !If no analytical approximation exists, use dtau_ref as dtau
            dtau = dtau_ref
!
            !First RK4 step with guessed dtau
            call rk4_step(z,dtau,dzdtau)
            tau = tau + dtau
!
            !Guess iface_new: minimum normal distance
            iface_new = minloc(abs(normal_distances_func(z(1:3))) ,1)
        endif
!
        !Extreme distance handling - If particle is already extremely far outside, use last line defense
        if(any(normal_distances_func(z(1:3)).gt.dist_max)) then
            call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                    & boole_last_line_defense_converged)
            if(.not.boole_last_line_defense_converged) then
                print *, 'Initial last line defense did not converge.'
                ind_tetr_inout = -1
                iface = -1
                call print_starting_condition()
                return
            endif
        endif
!
!
!----------------- Convergence and validation loop -----------------------!
!
        !Loop over all 4 faces to scan for the right converged exit point
        conv_val_loop: do i = 1,5 !number of faces + 1, because a particle is only converged, if it passes the whole loop
!
if(diag_pusher_tetra_rk) print *, 'conv_val_loop i',i
!
if(diag_pusher_tetra_rk) then
    print *, 'normal distances before Newton', normal_distances_func(z(1:3))
endif

            boole_converged = .true.        
!
            !Linear and quadratic Newton's method to converge at selected iface_new
            call newton_face_convergence(z,tau,iface_new,dzdtau,boole_accuracy_ode45,boole_newton_converged)
!
if(diag_pusher_tetra_rk) print *, 'Boolean Newton converged', boole_newton_converged
!
            !Propose new quadratic analytical approximation in case Newton's method didn't converge
            if(.not.boole_newton_converged) then
if(diag_pusher_tetra_rk) print *, 'Newton NOT converged'
                allowed_faces(iface_new) = .false.
                if(all(.not.allowed_faces)) then                    !All planes are forbidden -> try LLOD                       
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)
                    boole_converged = .false.  
                    cycle conv_val_loop
                endif
                
if(diag_pusher_tetra_rk) print *, 'normal distances after Newton:', normal_distances_func(z(1:3))


                !Extreme distance handling: Set back to initial values, if tau is too far outside
                if(tau.gt.dtau_quad) then
if(diag_pusher_tetra_rk) print *, 'Extreme distance handling using tau.'
                    z = z_init
                    tau = 0.d0
                endif
!
                call quad_analytic_approx(z,allowed_faces,iface_new,dtau,boole_quad_approx)
!                
                if(.not.boole_quad_approx) then
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged) 
                    boole_converged = .false.
                    cycle conv_val_loop
                endif 
!                 
                !Extreme distance prevention
                if(any(normal_distances_func(z(1:3)).gt.dist_max)) then
                    dtau = tau + dtau
                    tau = 0.d0
                    z = z_init
                endif
! 
                call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                tau = tau + dtau
!              
                boole_converged = .false.
                cycle conv_val_loop    
            endif
!            
            !Validation loop ('3-planes'-control)
            three_planes_loop: do j=1,3                             !Just consider the planes without the "exit-plane"
if(diag_pusher_tetra_rk) print *, 'three_planes loop i',j
if(diag_pusher_tetra_rk) print *, 'normal distances', normal_distances_func(z(1:3))
if(diag_pusher_tetra_rk) print *, 'iface_new',iface_new
                k=modulo(iface_new+j-1,4)+1
                if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside the tetrahedron
                    allowed_faces(iface_new) = .false.              !The respective face cannot be an exit plane
                    if(all(.not.allowed_faces)) then                !All planes are forbidden -> try LLOD                       
                        call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                                & boole_last_line_defense_converged)
                        boole_converged = .false.
                        cycle conv_val_loop
                    endif
!          
                    !Check, if the proposed j is forbidden or not
                    if(allowed_faces(k)) then                       !j is not forbidden
                        iface_new=k                                 !This means the orbit has to cross the k-plane before --> new iface
                        boole_converged = .false.
                        cycle conv_val_loop
                    else    
                        call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                            & boole_last_line_defense_converged)                    
                        boole_converged = .false.
                        cycle conv_val_loop
                    endif    
                endif        
            enddo three_planes_loop
!
            !Validation loop for vnorm
            if(normal_velocity_func(iface_new,dzdtau,z).gt.0.d0) then
if(diag_pusher_tetra_rk) print *, 'normal velocity = ', normal_velocity_func(iface_new,dzdtau,z)
                allowed_faces(iface_new) = .false.                  !The respective face cannot be an exit plane
                if(all(.not.allowed_faces)) then                    !All planes are forbidden -> try LLOD                       
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)
                    boole_converged = .false.
                    cycle conv_val_loop
                endif
!
                !Extreme distance handling: Set back to initial values, if tau is too far outside
                if(tau.gt.dtau_quad) then
if(diag_pusher_tetra_rk) print *, 'Extreme distance handling using tau.'
                    z = z_init
                    tau = 0.d0
                endif
!                
                call quad_analytic_approx(z,allowed_faces,iface_new,dtau,boole_quad_approx)
!                
                if(.not.boole_quad_approx) then
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)  
                    boole_converged = .false.
                    cycle conv_val_loop
                endif 
!                
                !Extreme distance prevention
                if(any(normal_distances_func(z(1:3)).gt.dist_max)) then
                    dtau = tau + dtau
                    tau = 0.d0
                    z = z_init
                endif
!                
                call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                tau = tau + dtau                
!
                boole_converged = .false.
                cycle conv_val_loop 
!                
            endif
!
            if(tau.le.0.d0) then
if(diag_pusher_tetra_rk) print *, 'tau < 0.d0 case'
                allowed_faces(iface_new) = .false.                  !The respective face cannot be an exit plane
                z = z_init
                tau = 0.d0
!                
                if(all(.not.allowed_faces)) then                    !All planes are forbidden -> try LLOD                       
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)
                    boole_converged = .false.
                    cycle conv_val_loop
                endif
!
                call quad_analytic_approx(z,allowed_faces,iface_new,dtau,boole_quad_approx)
!                
                if(.not.boole_quad_approx) then
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)
                    boole_converged = .false.
                    cycle conv_val_loop
                endif 
!                
                call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                tau = tau + dtau                
!
                boole_converged = .false.
                cycle conv_val_loop 
!
            endif
!
            exit conv_val_loop
!                            
        enddo conv_val_loop
!
        if(.not.boole_converged) then
            print *, 'Particle convergence failed'
            ind_tetr_inout = -1
            iface = -1
            call print_starting_condition()
            return
        endif
!
!-------------------------------------------------------------------!







!--------------------------------------------------------------------!
if(diag_pusher_tetra_rk) then
    print *, 'coordinates before final processing but after convergence procedure'
    print *, z
    print *, 'iface_new', iface_new
    print *, 'normal distances = ',normal_distances_func(z(1:3))
    print *, 'tau',tau
endif

!
        iface = iface_new
!
!         boole_accuracy_ode45 = .false.
        call final_processing(z,tau,iface,ind_tetr_inout,iper_phi,x,vpar,t_pass,boole_t_finished, &
                               & boole_accuracy_ode45,boole_final_processing_converged)
!
        if(.not.boole_final_processing_converged) then
            print *, 'Error in final processing: Could not converge even with ODE45'
            ind_tetr_inout = -1
            iface = -1
            call print_starting_condition()
            return
        endif

        
        !x_final without periodic boundary handling is needed for the calculation of FluxTV (time dt_dtau)
        z_final = z(1:3)
!         vperp = vperp_func(z_final)
        
!         if(present(x_final)) then
!             x_final = z(1:3)+tetra_physics(ind_tetr)%x1
!         endif
!

! if(diag_pusher_tetra_rk) then
        if((t_pass.ge.t_remain).and.(.not.boole_t_finished)) then
            print *, 'boole_t_finished was not set'
            ind_tetr_inout = -1
            iface = -1
            call print_starting_condition()
            return
        endif
        


        if(t_pass.le.0.d0) then
            print *, 't_pass is <= Zero'
            ind_tetr_inout = -1
            iface = -1
            call print_starting_condition()
            return
        endif
! endif        
!
        if(diag_pusher_tetra_rk)    print *, 'Particle is successfully converged.'
!
    end subroutine pusher_tetra_rk
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine analytic_coeff(norder,z,coef_mat)
!        
        use tetra_physics_mod, only: tetra_physics
        use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4
        use gorilla_diag_mod, only:diag_pusher_tetra_rk
!
        implicit none
!
        integer, intent(in)                 :: norder
        integer                             :: n
        double precision, dimension(4)      :: z
        double precision                    :: k1,k3
        double precision, dimension(4,norder+1)      :: coef_mat
!
        if(norder.ge.0) then
            coef_mat(:,1)=matmul(z(1:3),anorm)
            coef_mat(1,1)=coef_mat(1,1)-dist1
        endif
!
        if(norder.ge.1) then
            do n = 1,4
                coef_mat(n,2) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) + &
                                & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)) * z) + &
!                                    
!                                 & tetra_physics_poly4(ind_tetr)%anorm_in_b0(n) + &
!                                 & k1 * tetra_physics_poly4(ind_tetr)%anorm_in_b1(n) + &
!                                 & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_b2(n) + &
!                                 & k3 * tetra_physics_poly4(ind_tetr)%anorm_in_b3(n)
                                & sum(anorm(:,n) * b(1:3))
            enddo
        endif
!
        if(norder.ge.2) then
            perpinv2 = perpinv**2
            do n = 1,4
                coef_mat(n,3) = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,n) + &
                                & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,n) + &
                                & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,n)) * z) + &
!                    
!                                 & tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b0(n) + &
!                                 & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b0(n) + &
!                                 & k1 * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b1(n) + &
!                                 & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b1(n)) + &
!                                 & perpinv * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b2(n) + &
!                                 & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b2(n)) + &
!                                 & k3 * (tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b3(n) + &
!                                 & perpinv*tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b3(n))
!
                                & sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) + &
                                & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)) * b)
            enddo
        endif
!            
    end subroutine analytic_coeff
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine quad_analytic_approx(z,allowed_faces,iface_inout,dtau,boole_quad_approx)
!
        use tetra_physics_mod, only: tetra_physics
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!         use pusher_tetra_poly_mod, only: analytic_coeff
        use gorilla_settings_mod, only: boole_newton_precalc
        use constants, only: eps
!
        implicit none
!
        integer,intent(inout)                       :: iface_inout
        double precision,intent(out)                :: dtau
        double precision, dimension(4),intent(in)   :: z
        logical, dimension(4),intent(in)            :: allowed_faces
        logical,intent(out)                         :: boole_quad_approx
        integer                                     :: i,iface
        double precision                            :: discr,dummy
        double precision, dimension(4)              :: acoef,bcoef,ccoef,dtau_vec
        logical, dimension(4)                       :: boole_valid_dtau
!
        double precision, dimension(4,3)            :: coef_mat
        !Calculation of the quadratic equation constants

        if(boole_newton_precalc) then
            call analytic_coeff(2,z,coef_mat)
            ccoef = coef_mat(:,1)
            bcoef = coef_mat(:,2)
            acoef = coef_mat(:,3)
        else
            acoef= tetra_physics(ind_tetr)%acoef_pre
            bcoef=z(4)*acoef+matmul(b(1:3),anorm)
            acoef=acoef*(b(4)+spamat*z(4))
            ccoef = normal_distances_func(z(1:3))
        endif
!     
if(diag_pusher_tetra_rk) then
    print *, 'Quadratic analytical approximation'
    print *, 'acoef', acoef
    print *, 'bcoef', bcoef
    print *, 'ccoef', ccoef
    print *, 'normal_distances',normal_distances_func(z(1:3))
    print *, 'dist_min',dist_min
endif

        iface = iface_inout
if(diag_pusher_tetra_rk) print *, 'iface',iface
!        
        dtau_vec = dtau_max !Only results with positive times smaller than 'huge' will be evaluated
!
        !loop over faces
        faces_loop: do i = 1,4 
            if(.not.allowed_faces(i)) then  !Exclude forbidden faces
                cycle
            endif    
            if(iface.eq.i) then             !Particle starts from face i (ccoef(i) = 0.d0)
if(diag_pusher_tetra_rk) print *, 'i',i,'ccoef = 0'
                if(acoef(i).gt.0.d0) then   
                    if(bcoef(i).lt.0.d0) then
                        dtau_vec(i)=-2.d0*bcoef(iface)/acoef(iface)
                    else !acoef and bcoef have the same sign
                        cycle
                    endif
                elseif(acoef(i).lt.0.d0) then   
                    if(bcoef(i).gt.0.d0) then
                        dtau_vec(i)=-2.d0*bcoef(iface)/acoef(iface)
                    else !acoef and bcoef have the same sign
                        cycle
                    endif
                else ! acoef = 0
                    cycle !Orbit is already converged at iface. Linear solution with same iface would not yield any result.
                endif
            else !ccoef != 0
if(diag_pusher_tetra_rk) print *, 'abs(ccoef(i))',abs(ccoef(i))
                if(abs(ccoef(i)).gt.dist_min) then
                if(ccoef(i).gt.0.d0) then
if(diag_pusher_tetra_rk) print *, 'i',i,'ccoef > 0'
                    if(acoef(i).gt.0.d0) then
                        if(bcoef(i).lt.0.d0) then
                            discr=bcoef(i)**2-2.d0*acoef(i)*ccoef(i)
                            if(discr.gt.0.d0) then
                                dummy = (-bcoef(i) +sqrt(discr))
                                if( abs(dummy).gt.eps ) then
                                    dtau_vec(i)= 2.d0 * ccoef(i) / dummy
                                else
                                    dtau_vec(i)=(-sqrt(discr)-bcoef(i))/acoef(i) !Numerically unstable
                                endif
                            elseif(discr.eq.0.d0) then
                                dtau_vec(i)=-bcoef(i)/acoef(i)
                            else !discr < 0
                                cycle
                            endif
                        else ! b >= 0
                            cycle
                        endif
                    elseif(acoef(i).lt.0.d0) then
                        discr=bcoef(i)**2-2.d0*acoef(i)*ccoef(i)
                        dummy = (-bcoef(i) +sqrt(discr))
                        if( abs(dummy).gt.eps ) then
                            dtau_vec(i)= 2.d0 * ccoef(i) / dummy
                        else
                            dtau_vec(i)=(-sqrt(discr)-bcoef(i))/acoef(i) !Numerically unstable
                        endif
                    else !acoef = 0.d0
                        if(bcoef(i).lt.0.d0) then
                            dtau_vec(i) = -ccoef(i)/bcoef(i)
                        else
                            cycle
                        endif
                    endif
                elseif(ccoef(i).lt.0.d0) then
if(diag_pusher_tetra_rk) print *, 'i',i,'ccoef < 0'
                    if(acoef(i).lt.0.d0) then
                        if(bcoef(i).gt.0.d0) then
                            discr=bcoef(i)**2-2.d0*acoef(i)*ccoef(i)
                            if(discr.gt.0.d0) then
                                dtau_vec(i)=(sqrt(discr)-bcoef(i))/acoef(i)
                            elseif(discr.eq.0.d0) then
                                dtau_vec(i)=-bcoef(i)/acoef(i)
                            else !discr < 0
                                cycle
                            endif
                        else ! b <= 0
                            cycle
                        endif
                    elseif(acoef(i).gt.0.d0) then
                        discr=bcoef(i)**2-2.d0*acoef(i)*ccoef(i)
                        dtau_vec(i)=(sqrt(discr)-bcoef(i))/acoef(i)
                    else !acoef = 0.d0
                        if(bcoef(i).gt.0.d0) then
                            dtau_vec(i) = -ccoef(i)/bcoef(i)
                        else
                            cycle
                        endif
                    endif
                else
                    print *, 'Error: Should not happen'
                    stop
                endif
                else !ccoef = 0.d0
                    if(((acoef(i).gt.0d0).and.(bcoef(i).lt.0d0)).or.((acoef(i).lt.0d0).and.(bcoef(i).gt.0d0))) then !If a and b have opposing signs
                        dtau_vec(i) = -2.d0 * bcoef(i)/acoef(i)
                    else
                        cycle
                    endif    
                endif
            endif
        enddo faces_loop
!
        boole_valid_dtau = (dtau_vec.lt.dtau_max).and.(dtau_vec.gt.0.d0) !dtau_max
        if( any( boole_valid_dtau) ) then
            boole_quad_approx = .true.
!
            iface_inout = minloc(dtau_vec,1,boole_valid_dtau)
            dtau = dtau_vec(iface_inout)
        else
            boole_quad_approx = .false.
        endif
!
if(diag_pusher_tetra_rk) then
    print *, 'dtau_vec', dtau_vec
    print *, 'dtau_max', dtau_max
    print *, 'selected dtau', dtau    
endif    

!
    end subroutine quad_analytic_approx
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
    subroutine rhs_pusher_tetra_rk4(z,dzdtau)
!
        implicit none
!
        ! In/out variables
        ! z(1:3) ... coordinates of the particle
        ! z(4) ... parallel velocity of the particle (vpar)
!        
        ! Output variables
        ! dzdtau ... left hand side of ODE set
!  
        ! Constant module variables
        ! b ... const. vector of differential equation
        ! amat ... const. Matrix of differential equation
        ! spamat ... trace of amat
        ! Bvec ... vector of magnetic field B
!  
        double precision, dimension(4),intent(inout) :: z
        double precision, dimension(4),intent(out) :: dzdtau
!
        dzdtau(1:3)=b(1:3)+matmul(amat,z(1:3))+Bvec*z(4)
        dzdtau(4)=b(4)+spamat*z(4)
!
    end subroutine rhs_pusher_tetra_rk4
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine rk4_step(y,h,dzdtau)
!
        implicit none
!
        ! In/out variables
        ! y(1:3) ... coordinates of the particle
        ! y(4) ... parallel velocity of the particle (vpar)
        ! h ... dtau (time for integration)
!  
        ! Constant module variables
        ! b ... const. vector of differential equation
        ! amat ... const. Matrix of differential equation
        ! spamat ... trace of amat
        ! Bvec ... vector of magnetic field B
!  
        ! Output variables
        ! dzdtau ... velocity vector of the particle at y(t_end)
!
        integer, parameter :: n=4  
!
        double precision, dimension(n), intent(inout) :: y
        double precision, dimension(n)                :: y_save
        double precision, intent(in) :: h
        double precision             :: h_save
        double precision, dimension(4),intent(out) :: dzdtau
!
        !RK variables:
        double precision :: hh,h6
        double precision, dimension(n) :: dydx,yt,dyt,dym
!
        !Save input variables
        y_save = y
        h_save = h
!
        hh=h*0.5d0
        h6=h/6.d0
!
        call rhs_pusher_tetra_rk4(y,dydx)
!
        yt=y+hh*dydx
!
        call rhs_pusher_tetra_rk4(yt,dyt)
!
        yt=y+hh*dyt
!
        call rhs_pusher_tetra_rk4(yt,dym)
!
        yt=y+h*dym
        dym=dyt+dym
!
        call rhs_pusher_tetra_rk4(yt,dyt)
!
        y=y+h6*(dydx+dyt+2.d0*dym)
!
        dzdtau=dyt
!  
    end subroutine rk4_step
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine rhs_pusher_tetra_rk45(dummy,z,dzdtau)
!
      implicit none
!
      double precision, dimension(4) :: z,dzdtau
      double precision :: dummy
!
      dzdtau(1:3) = b(1:3)+matmul(amat,z(1:3))+Bvec*z(4)
      dzdtau(4) = b(4)+spamat*z(4)
!
    end subroutine rhs_pusher_tetra_rk45
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine newton_face_convergence(z,tau,iface,dzdtau,boole_accuracy_ode45_in,boole_newton_converged,boole_start_quadratic)
!
        !This function is a wrapper for Newton's method to assure fast ODE45 accuracy
!
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!
        implicit none
!
        integer,intent(in)                              :: iface
        double precision, intent(inout)                 :: tau
        double precision, dimension(4),intent(inout)    :: dzdtau
        double precision, dimension(4),intent(inout)    :: z   
        double precision, dimension(4)                  :: z_save,dzdtau_save
        logical, intent(in)                             :: boole_accuracy_ode45_in
        logical,intent(out)                             :: boole_newton_converged
        logical,intent(in),optional                     :: boole_start_quadratic
        double precision                                :: tau_save,dtau
        logical                                         :: boole_accuracy_ode45
!
        z_save = z
        tau_save = tau
        dzdtau_save = dzdtau
!
        !Call Newton with RK4 method
        boole_accuracy_ode45 = .false.
      
        call newton_face_convergence_wrapped(z,tau,iface,dzdtau,boole_accuracy_ode45,boole_newton_converged, &
                                               & boole_start_quadratic)                                             
!
        if(.not.boole_newton_converged) then
            z = z_save
            tau = tau_save
            dzdtau = dzdtau_save
            return
        endif    
!
        if(boole_accuracy_ode45_in) then
!
            z = z_save
            dtau = tau - tau_save
!
if(diag_pusher_tetra_rk) then
    print *, 'original tau', tau_save
    print *, 'tau after RK4', tau
    print *, 'Newton Face convergence: dtau of RK4', dtau
    print *, 'normal distance before ODE45 step', normal_distances_func(z(1:3))
endif    
            

            !Repeat RK4-Newton with ODE45
            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45_in)
            !(tau = tau_save + dtau)
if(diag_pusher_tetra_rk) print *, 'normal distance after ODE45 step', normal_distances_func(z(1:3))
!
!             !Extreme distance prevention
!             if(any(abs(normal_distances_func(z(1:3))).gt.dist_max)) then
! if(diag_pusher_tetra_rk) print *, 'Extreme distance prevention in wrapper'
!                 z = z_save
!                 tau = tau_save
!                 dzdtau = dzdtau_save
!             endif
!
            !Converge with ODE45-Newton
            call newton_face_convergence_wrapped(z,tau,iface,dzdtau,boole_accuracy_ode45_in,boole_newton_converged, &
                                               & .false.)
        endif
!    
    end subroutine newton_face_convergence
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine newton_face_convergence_wrapped(z,tau,iface,dzdtau,boole_accuracy_ode45,boole_newton_converged, &
                                               & boole_start_quadratic_in)
!
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!
        implicit none
!
        integer,intent(in)                              :: iface
        double precision, intent(inout)                 :: tau
        double precision, dimension(4),intent(inout)    :: dzdtau
        double precision, dimension(4)                  :: dzdtau_save,dzdtau_start
        double precision, dimension(4),intent(inout)    :: z   
        double precision, dimension(4)                  :: z_save,z_start
        logical, intent(in)                             :: boole_accuracy_ode45
        logical,intent(out)                             :: boole_newton_converged
        logical,intent(in),optional                     :: boole_start_quadratic_in
        integer                                         :: k,j,i
        double precision                                :: dtau,dtau_save,tau_save,tau_start,dist,dist_new,discr,dummy
        double precision                                :: normal_velocity,normal_acceleration
        logical                                         :: boole_start_quadratic
!
if(diag_pusher_tetra_rk) print *, 'Newton - Accuracy:', boole_accuracy_ode45
!     
select case(iface)
    case(1:4)
!
        !Save initial values, before this routine tries to converge particle on face
        z_start = z
        tau_start = tau
        dzdtau_start = dzdtau
!
        !Default value for convergence
        boole_newton_converged = .false.
!        
        !Default value for boole_start_quadratic              
        !For an optional special convergence procedure that needs to START with Newton2 in order to avoid finding the wrong intersections in final processing.
        !The second intersection will have a (1) too large time which is not allowed and (2) will fly inside.
        boole_start_quadratic = .false.
        if(present(boole_start_quadratic_in)) boole_start_quadratic = boole_start_quadratic_in
!
        ! Calculate parameters for Newton's method
        dist = normal_distance_func(z(1:3),iface)        
!        
        k = 0
!
        convergence_loop: do while(abs(dist).gt.dist_min)       
!
            k = k+1
! if(diag_pusher_tetra_rk) print *, 'k',k
!
            !Save variable 
            z_save=z
            dzdtau_save=dzdtau
!        
            !Linear Newton's guess
            normal_velocity = normal_velocity_func(iface,dzdtau,z)
            if(normal_velocity.ne.0.d0) then
                dtau = -dist/normal_velocity !minus sign, because flying outwards has negative velocity sign
            else
if(diag_pusher_tetra_rk) print *, 'Error in Newtons method: normal velocity equals zero.'
                return
            endif
!
            !Save dtau for quadratic Newton (in case of extreme distance prevention)
            dtau_save = dtau
            tau_save = tau
!
            !Extreme distance prevention
            if(any(normal_distances_func(z(1:3)).gt.dist_max)) then
if(diag_pusher_tetra_rk) print *, 'Newton extreme distance prevention'
                dtau = tau + dtau
                tau = 0.d0
                z = z_init
            endif
!
            !Integration step of prevention in case of ODE45
!            if(boole_accuracy_ode45) then
                !Extreme dtau prevention OR integration step
                if(abs(dtau).gt.dtau_max) then
if(diag_pusher_tetra_rk) print *, 'Extreme dtau prevention'
                    boole_start_quadratic = .true.
                else
if(diag_pusher_tetra_rk) print *, 'linear Newton step'

! if(diag_pusher_tetra_rk) then
!     print *, ''
!     print *, 'boole_accuracy', boole_accuracy_ode45
!     print *, 'dtau quadtratic', dtau
!     print *, 'dtau_max',dtau_max
!     print *, 'normal distances', normal_distances_func(z(1:3))
!     print *, 'dtau through normal_velocity', min_dist_ref / maxval([(abs(normal_velocity_func(i,dzdtau,z)),i=1,4)],1)
!     print *, 'total tau after step', tau + dtau
! endif  
                    !RK4 step with linear Newton's guess  
                    call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
    !                
                    !Distance to plane after RK4-step with Newton's guess
                    dist_new = normal_distance_func(z(1:3),iface)
                endif
                !Integration step, because dtau is within range (dtau_max)
!            else
!                !RK4 step with linear Newton's guess
!                call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
!!
!                !Distance to plane after RK4-step with Newton's guess
!                dist_new = normal_distance_func(z(1:3),iface)
!            endif
! 
            !Quadratic Newton procedure, if Newton-step yielded a worse result (larger normal distance)
            if((abs(dist_new).ge.abs(dist)).or.boole_start_quadratic) then
!
if(diag_pusher_tetra_rk) print *, 'quadratic Newton'
                boole_start_quadratic = .false.    
!
                !Use z value, before Newton yielded a worse result
                z = z_save
                dzdtau = dzdtau_save !just to make sure that the value for the output is correct (if discriminant is negative)
                tau = tau_save
!
                !Normal acceleration
                normal_acceleration = 0.5d0*normal_acceleration_func(iface,dzdtau,z)
!
                !Discriminant of quadratic equation
                discr = normal_velocity**2-4.d0*normal_acceleration*dist
!
                !Discriminant positive: intersectionany(normal_distances_func(z(1:3)).gt.dist_max)
                if(discr.gt.0.d0) then 
                    !If the distance to the plane is quadratic in time, discriminant has to be eq. o. gt. than zero
!
                    if(normal_acceleration.lt.0.d0) then
                        dtau = (-normal_velocity - sqrt(discr))/(2.d0*normal_acceleration)
                    elseif(normal_acceleration.gt.0.d0) then
                        dtau = (-normal_velocity + sqrt(discr))/(2.d0*normal_acceleration)
                    else ! normal acceleration equals zero
!                        print *, 'This case should not happen.'
!                        stop
                        dtau = -dist/normal_velocity
                    endif
!
                    !Extreme distance prevention
                    if(any(normal_distances_func(z(1:3)).gt.dist_max)) then
if(diag_pusher_tetra_rk) print *, 'extreme distance prevention'
                        dtau = tau + dtau
                        tau = 0.d0
                        z = z_init
                    endif
!
                    if(abs(dtau).gt.dtau_max) then
if(diag_pusher_tetra_rk) print *, 'Extreme dtau in quadratic approximation'
                        z = z_start
                        tau = tau_start
                        dzdtau = dzdtau_start
                        return
                    endif
!
                    !RK4 step with quadratic Newton's guess
                    call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                    tau = tau + dtau
                    dist = normal_distance_func(z(1:3),iface)
!                
                !Discriminant negative or zero: No intersection
                else     
!                    
                    !z and tau values from 'before' the step, where Newton yielded worse result, are returned.
if(diag_pusher_tetra_rk) print *, 'Newton: negative discriminant'
                    return !Exit the subroutine
                endif
!                
            !Newton iteration yielded a better result (smaller normal distance)
            else
if(diag_pusher_tetra_rk) print *, 'Linear Newton'
                tau = tau + dtau
                dist = dist_new                                       
            endif
!
            !Number of Newton iterations exceed a certain number
            if(k.gt.kiter) then
if(diag_pusher_tetra_rk)   print *, 'Error: Newton convergence procedure fails: Maximum number of Newton iterations exceed.'
!                 call print_starting_condition()
                return
            endif
!
        enddo convergence_loop
!       
        !Negative tau convergence --> Return values before Newton converged on wrong plane
        if(tau.le.0.d0) then
if(diag_pusher_tetra_rk) print *, 'Newton: negative tau'
            z = z_start
            tau = tau_start
            dzdtau = dzdtau_start
            return
        endif
!
        !Particle is converged on face
        boole_newton_converged = .true.
!
case default
!
    print *, 'Error in Newton face convergence. iface ne 1,2,3,4'
    call print_starting_condition()
    stop
!
end select
!
    end subroutine newton_face_convergence_wrapped
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!

    subroutine newton_face_convergence1(z,tau,iface,dzdtau,boole_accuracy_ode45_in,boole_newton_converged,boole_start_quadratic_in)
!
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!
        implicit none
!
        integer,intent(in)                              :: iface
        double precision, intent(inout)                 :: tau
        double precision, dimension(4),intent(inout)    :: dzdtau
        double precision, dimension(4)                  :: dzdtau_save,dzdtau_start
        double precision, dimension(4),intent(inout)    :: z   
        double precision, dimension(4)                  :: z_save,z_start
        logical, intent(in)                             :: boole_accuracy_ode45_in
        logical,intent(out)                             :: boole_newton_converged
        logical,intent(in),optional                     :: boole_start_quadratic_in
        integer                                         :: k,j
        double precision                                :: dtau,dtau_save,tau_save,tau_start,dist,dist_new,discr,dummy
        double precision                                :: normal_velocity
        logical                                         :: boole_newton1_tangential,boole_newton2_tangential,boole_start_quadratic
        logical                                         :: boole_accuracy_ode45
if(diag_pusher_tetra_rk) then
print *, 'dist_min', dist_min
    z_save = z_init
    j = 1000
    tau_save = 0.d0
    dtau_save = 1.d-20/dble(j)
write(19,*) tau_save, normal_distance_func(z_save(1:3),iface)    
    do k = 1,j
        call integration_step(z_save,dtau_save,dzdtau_save,.true.)
        tau_save = tau_save + dtau_save
        write(19,*) tau_save, normal_distance_func(z_save(1:3),iface)
    enddo    
     stop
endif        !
select case(iface)
    case(1:4)
!
boole_accuracy_ode45 = .false.
!
        !Save initial values, before this routine tries to converge particle on face
        z_start = z
        tau_start = tau
        dzdtau_start = dzdtau
!
        !Default value for convergence
        boole_newton_converged = .false.
!
        !Default values for handling of extreme distance cases (tangential)
        boole_newton1_tangential = .false.
        boole_newton2_tangential = .false.
!
        !Default value for boole_start_quadratic              
        !For an optional special convergence procedure that needs to START with Newton2 in order to avoid finding the wrong intersections in final processing.
        !The second intersection will have a (1) too large time which is not allowed and (2) will fly inside.
        boole_start_quadratic = .false.
        if(present(boole_start_quadratic_in)) boole_start_quadratic = boole_start_quadratic_in
!
        ! Calculate parameters for Newton's method
        dist = normal_distance_func(z(1:3),iface)        
!        
        k = 0
        convergence_loop: do while(abs(dist).gt.dist_min)       
!
if(diag_pusher_tetra_rk) print *, 'k',k
            !Save variable 
            z_save=z
            dzdtau_save=dzdtau
!        
if(diag_pusher_tetra_rk) then
    print *, 'vnorm'
    print *, normal_velocity_func(iface,dzdtau,z)
endif    
!
if(diag_pusher_tetra_rk) then
    if(normal_velocity_func(iface,dzdtau,z).eq.0.d0) stop
endif

            !Linear Newton's guess
            normal_velocity = normal_velocity_func(iface,dzdtau,z)
            if(normal_velocity.ne.0.d0) then
                dtau = -dist/normal_velocity !minus sign, because flying outwards has negative velocity sign
            else
if(diag_pusher_tetra_rk) print *, 'Error in Newtons method: normal velocity equals zero.'
                return
            endif
!
            !Save dtau for quadratic Newton (in case of extreme distance prevention)
            dtau_save = dtau
            tau_save = tau
!
            !Extreme distance prevention
            boole_newton1_tangential = any(normal_distances_func(z(1:3)).gt.dist_max)
            if(boole_newton1_tangential) then
if(diag_pusher_tetra_rk) print *, 'Newton extreme distance prevention'
                dtau = tau + dtau
                tau = 0.d0
                z = z_init
            endif
!
            !RK4 step with linear Newton's guess  
            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
!
if(diag_pusher_tetra_rk) write(17,*) tau+dtau,normal_distance_func(z(1:3),iface)
            !Distance to plane after RK4-step with Newton's guess
            dist_new = normal_distance_func(z(1:3),iface)
!
            !Quadratic Newton procedure, if Newton-step yielded a worse result (larger normal distance)
            !if(abs(dist_new).ge.abs(dist)) then
            
            
            if((abs(dist_new).ge.abs(dist)).or.(boole_newton1_tangential.and.(abs(dist_new).gt.dist_min)) &
                & .or.boole_start_quadratic) then
!
if(diag_pusher_tetra_rk) print *, 'quadtratic Newton'
                if(present(boole_start_quadratic_in)) boole_start_quadratic = .false.
!            
                !Use z value, before Newton yielded a worse result
                z = z_save
                dzdtau = dzdtau_save !just to make sure that the value for the output is correct (if discriminant is negative)
                tau = tau_save
!
                !Discriminant of quadratic equation
                discr=dist*(dist-4.d0*dist_new)
!
                !Discriminant positive: intersection
                !if((discr.gt.0.d0)) then
                if((discr.gt.0.d0).and.(.not.boole_newton2_tangential)) then 
                  !If the distance to the plane is quadratic in time, discriminant has to be eq. o. gt. than zero
                    dummy=0.5d0*dtau_save/dist_new                 ! This means that the orbit has an intersection with the plane
                    !Distinguish if the acoef of the quadratic equation is positive or negative --> sign of dist.
                    if(dist.lt.0.d0) then
                        !Distinguish cases depending on the sign of dtau. (By Simplifying the formula the sign of the solution changes.)
                        if(dummy.gt.0.d0) then
                            !dtau > 0, dist_new > 0
                            dtau=dummy*(dist-sqrt(discr))
                        else
                            !dtau < 0, dist_new > 0
                            dtau=dummy*(dist+sqrt(discr))
                        endif
                    else
                        !Here the sign is the otherway arround (depending on dummy)
                        if(dummy.gt.0.d0) then
                            !dtau < 0, dist_new < 0
                            dtau=dummy*(dist+sqrt(discr))
                        else
                            !dtau > 0, dist_new < 0
                            dtau=dummy*(dist-sqrt(discr))
                        endif
                    endif
!                    
                    !!!! ATTENTION: Super sloppy workaround
                    boole_newton2_tangential = boole_newton1_tangential
!
                    !Extreme distance prevention
                    boole_newton1_tangential = any(normal_distances_func(z(1:3)).gt.dist_max)
                    if(boole_newton1_tangential) then
if(diag_pusher_tetra_rk) print *, 'extreme distance prevention'
                        dtau = tau + dtau
                        tau = 0.d0
                        z = z_init
                    endif
!
                    !RK4 step with quadratic Newton's guess
                    call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                    tau = tau + dtau
                    dist = normal_distance_func(z(1:3),iface)
!                
! write(18,*)     tau, dist
                !Discriminant negative or zero: No intersection
                else     
!                    
                    !z and tau values from 'before' the step, where Newton yielded worse result, are returned.
if(diag_pusher_tetra_rk) print *, 'Newton: negative discriminant'
                    return !Exit the subroutine
                endif
!                
            !Newton iteration yielded a better result (smaller normal distance)
            else
if(diag_pusher_tetra_rk) print *, 'Linear Newton'
                tau = tau + dtau
                dist = dist_new           
! write(18,*)     tau, dist                
            endif
!
            !Number of Newton iterations exceed a certain number
            if(k.gt.kiter) then
                print *, 'Error: Newton convergence procedure fails: Maximum number of Newton iterations exceed.'
                call print_starting_condition()
                return
            endif
!
        enddo convergence_loop
!       
        !Negative tau convergence --> Return values before Newton converged on wrong plane
        if(tau.le.0.d0) then
if(diag_pusher_tetra_rk) print *, 'Newton: negative tau'
            z = z_start
            tau = tau_start
            dzdtau = dzdtau_start
            return
        endif
!
        !Particle is converged on face
        boole_newton_converged = .true.
!        

case default

    print *, 'Error in Newton face convergence. iface ne 1,2,3,4'
    call print_starting_condition()
    stop

end select
    end subroutine newton_face_convergence1
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine bisection_face_convergence(z,tau_inout,dtau_in,iface,dzdtau,boole_accuracy_ode45,boole_bisection_converged)
!
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
        use supporting_functions_mod, only: logical2integer
!
        implicit none
!
        double precision, dimension(4), intent(inout)       :: z
        double precision, intent(inout)                     :: tau_inout
        double precision, intent(in)                        :: dtau_in
        integer, intent(out)                                :: iface
        double precision, dimension(4), intent(out)         :: dzdtau
        logical, intent(out)                                :: boole_bisection_converged  
        logical, intent(in)                                 :: boole_accuracy_ode45
        double precision                                    :: dtau,tau,dtau_save, tau_save, tau_ref
        integer                                             :: i,j,k,l,n_faces_converged
        double precision, dimension(4)                      :: normal_distances, z_save
!
if(diag_pusher_tetra_rk) print *, 'Bisection convergence procedure: '


! if(diag_pusher_tetra_rk) then
! 
!      j = 1000
!
!
!      dtau_save = tau_inout/dble(j)
!      tau_save = 0.d0
!      z_save = z_init
!      write(18,*) tau_save, normal_distances_func(z_save(1:3))
!      do i = 1,j
!          call integration_step(z_save,dtau_save,dzdtau,.true.)
!          tau_save = tau_save + dtau_save
!          write(18,*) tau_save, normal_distances_func(z_save(1:3))
!      enddo
!
!      call bisection_search_start(tau_inout,1000,boole_accuracy_ode45,z_save,dtau,tau_save)
!      write(19,*) tau_save, normal_distances_func(z_save(1:3))
!
!  stop
! 
! 
! ! Check, if particle crossed a plane in previous dtau step
!     z_save = z
!     normal_distances = sign([1.d0,1.d0,1.d0,1.d0],normal_distances_func(z_save(1:3)))
!     dtau = -dtau_in
!     
!     call integration_step(z_save, dtau,dzdtau,.true.)
!     if(all((sign([1.d0,1.d0,1.d0,1.d0],normal_distances_func(z_save(1:3))) * normal_distances).eq.1.d0)) then
!         print *, 'Particle did not cross the plane in previous step.'
!         stop
!     endif
! endif


!
        tau = tau_inout     !Cumulative integration time (tau) until call of bisection procedure
        tau_ref = tau       !Tau at the beginning of bisection
        dtau = dtau_in      !Last integration step size before bisection procedure
        z_save = z          !Save z for adaptive convergence region loop
!

!
if(diag_pusher_tetra_rk) print *, 'dtau', dtau
if(diag_pusher_tetra_rk) print *, 'tau', tau
if(diag_pusher_tetra_rk) print *, 'dist_min',dist_min



if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(93,*) -1,tau, normal_distances_func(z(1:3))


        boole_bisection_converged = .false.
!
        !Outer loop for possibly changing the convergence region (Only in ODE45 accuracy mode)
        convergence_region_loop: do l = 1,2
!        
            k = 0
            do while (.not.boole_bisection_converged)
!        
if(diag_pusher_tetra_rk) then
    print *, 'k = ',k
    print *, 'dtau',dtau
    print *, 'normal_distances', normal_distances_func(z(1:3))
    call rk4_step(z,0.d0,dzdtau)
    print *, 'normal velocites', [(normal_velocity_func(i,dzdtau,z),i=1,4)]
    print *, 'tau',tau
endif    


                normal_distances = normal_distances_func(z(1:3))
!
                if (minval(normal_distances(:),1) .lt. -dist_min) then
if(diag_pusher_tetra_rk) print *, 'Bisection: Outside the convergence region and outside the tetrahedron'
if(diag_pusher_tetra_rk) print *, 'tau-tau_ref', tau-tau_ref
                    dtau = -abs(dtau/2.d0)
!                
                    !Extreme distance prevention
                    if(any(normal_distances.gt.dist_max)) then
if(diag_pusher_tetra_rk) print *, 'extreme distance'
                        dtau_save = dtau
                        dtau = tau + dtau
                        tau = 0.d0
                        z = z_init
!                    
                        !Runge-Kutta step
                        call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                        tau = tau + dtau
                        dtau = dtau_save
if(diag_pusher_tetra_rk) write(93,*) k,tau, normal_distances_func(z(1:3))
                    else                 
                        !Runge-Kutta step    
                        call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                        tau = tau + dtau
if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(93,*) k, tau, normal_distances_func(z(1:3) )
                    endif
                    
                    
                elseif(minval(normal_distances(:),1) .gt. dist_min) then
                    dtau = +abs(dtau/2.d0)
!
                    !Extreme distance prevention is not needed in this case. (Particle is inside the tetrahedron.)
!                
                    !Runge-Kutta step
                    call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                    tau = tau + dtau
if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(93,*) k,tau, normal_distances_func(z(1:3)     )
                endif
!            
                normal_distances = normal_distances_func(z(1:3))
!            
                if (abs(minval(normal_distances(:),1)) .lt. dist_min) then

if(diag_pusher_tetra_rk) print *, 'Bisection: Inside convergence region'
!  
                    if(normal_velocity_func(Minloc(normal_distances(:),1),dzdtau,z).gt.0.d0) then
                        dtau = +abs(dtau/2.d0)
                        !Runge-Kutta step
                        call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                        tau = tau + dtau
if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(93,*) k,tau, normal_distances_func(z(1:3) )
                    else
                        j = 0
                        do i = 1,4
                            if(normal_distances(i).lt.0.d0) j = j+1
                        enddo
                        if(j.le.1) then
                            iface = Minloc(normal_distances(:),1)
                            boole_bisection_converged = .true.
                        else
                            dtau = -abs(dtau/2.d0)
                            !Runge-Kutta step 
                            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                            tau = tau + dtau
if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(93,*) k,tau, normal_distances_func(z(1:3))
                        endif
                    endif
                endif
!            
                k = k+1
                if(k.gt.(kiter)) then
!
                    if(l.eq.1) then
                        !if(boole_accuracy_ode45) then
!                    
                            !Check for higher order convergence problem:
                            !The particle is inside the convergence region regarding more than one face
!                            
                            !Number of particles outside the tetrahdron, but inside the convergence region
                            n_faces_converged = sum(logical2integer( &
                                            & (normal_distances.lt.0.d0).and.(abs(normal_distances).lt.dist_min)),1)
!                                            
                            if(n_faces_converged.gt.1) then
                                dist_min = 2.d0 * dist_min
                                tau = tau_inout
                                dtau = dtau_in     
                                z = z_save
if(diag_pusher_tetra_rk) then
                            print *, 'Higher order convergence problem.'
                            print *, 'Number of particles outside the tetrahdron, &
                                      & but inside the convergence region:',n_faces_converged
endif                                
                                exit
                            else
!                                print *, 'Bisection failed: Particle did not converge.'
!                                print *, 'Number of particles outside the tetrahdron, &
!                                        & but inside the convergence region:',n_faces_converged
!                                call print_starting_condition()
!                                stop
!
if(diag_pusher_tetra_rk)  print *, 'Bisection: Find new starting position'
                                !Search appropriate start position for bisection procedure
                                call bisection_search_start(tau_inout,1000,boole_accuracy_ode45,z,dtau,tau)
!                                
                                exit
                            endif
                        !else
                        !    exit convergence_region_loop
                        !endif
                    else
                        exit convergence_region_loop
                    endif
!
                endif
!            
            enddo !while (.not.boole_bisection_converged)
!
            if(boole_bisection_converged) exit convergence_region_loop
!
        enddo convergence_region_loop
!
        tau_inout = tau
!
    end subroutine bisection_face_convergence
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine bisection_search_start(tau_in,n_steps,boole_accuracy_ode45,z_start,dtau,tau_out)
!
    !This subroutine should be used, when the standard bisection procedure fails,
    !because the wrong starting point was chosen.
    !In such a case, the orbit is intgrated stepwise for a total of 'n_steps'.
    !An appropriate starting position for a consecutive bisection procedure is then chosen.
    !
    !   Input:
    !
    !   tau_in  ... Integration time, for which the previous bisection method DID NOT find a correct exit point
    !   n_steps ... Number of steps for which orbit integration should be split
    !   boole_accuracy_ode45 ... Accuracy of orbit integration
    !
    !   Output:
    !
    !   z_start ... Proposal for new starting position of bisection procedure
    !               (Starting point is chosen, when orbit left the tetrahedron.
    !               Particle is still inside in the previous integration step)
    !   dtau    ... Integration time between time steps
    !   tau_out ... Cumulative integration time from entry point until z_start
!
        use supporting_functions_mod, only: logical2integer
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
!
        integer, intent(in)                         :: n_steps
        double precision, intent(in)                :: tau_in
        logical, intent(in)                         :: boole_accuracy_ode45
        double precision, dimension(4), intent(out) :: z_start
        double precision, intent(out)               :: dtau,tau_out
        integer                                     :: i,start_index
        logical, dimension(n_steps)                 :: boole_particle_inside
        double precision, dimension(n_steps)        :: tau_vec
        double precision, dimension(4)              :: z_save,dzdtau
        double precision, dimension(4,n_steps)      :: z_mat, normal_distances_mat
!
        !Save starting values in matrix
        z_save = z_init
        normal_distances_mat(:,1) = normal_distances_func(z_save(1:3))
        z_mat(:,1) = z_save
        tau_vec(1) = 0.d0
!
        !Define step length
        dtau = tau_in/dble(n_steps)
!
        !Integrate orbit and save normal distances and position
        do i = 2,n_steps
            call integration_step(z_save,dtau,dzdtau,boole_accuracy_ode45)
            z_mat(:,i) = z_save
            normal_distances_mat(:,i) = normal_distances_func(z_save(1:3))
            tau_vec(i) = tau_vec(i-1) + dtau
            
            write(20,*) tau_vec(i), normal_distances_mat(:,i)
        enddo
!
        !Boolean, true, if particle is inside the tetrahedron
        boole_particle_inside = all(normal_distances_mat.gt.0.d0,DIM=1)
!
        !Find starting index for bisection:
        !Bisection should start, when orbit just left the tetrahedron.
        !Particle is still inside in the previous integration step (+1)
        start_index = MaxLoc( [(i ,i = 1,n_steps)] , DIM = 1,MASK =boole_particle_inside) + 1
!
        z_start = z_mat(:,start_index)
        tau_out = tau_vec(start_index)
!
    end subroutine bisection_search_start
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine last_line_defense(z,tau,iface,dzdtau,boole_accuracy_ode45,boole_last_line_defense_converged)
!    
        use supporting_functions_mod, only: logical2integer
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
        use tetra_physics_mod, only: tetra_physics
!
        implicit none
!
        double precision, dimension(4), intent(out)             :: z
        double precision, intent(out)                           :: tau
        integer, intent(out)                                    :: iface
        double precision, dimension(4), intent(out)             :: dzdtau
        logical,intent(in)                                      :: boole_accuracy_ode45
        double precision                                        :: dtau,vnorm,tau_save,dtau_save
        double precision, dimension(4)                          :: z_save,normal_distances
        logical, dimension(4)                                   :: allowed_faces, l_normal_distances
        integer, dimension(4)                                   :: i_normal_distances
        logical                                                 :: boole_quad_approx,boole_turned_tangential
        logical                                                 :: boole_distance_bisection
        logical                                                 :: boole_converged,boole_bisection_converged
        logical                                                 :: boole_newton_converged,boole_dtau_decreased
        logical,intent(out)                                     :: boole_last_line_defense_converged
        integer                                                 :: iface_new,iface_init_outside,n_faces_outside
        integer                                                 :: i,j,k,l
!
        !Initialize booleans
        boole_turned_tangential = .false.       !Particle is outside the initial face and has tangentially crossed another face
        boole_converged = .false.               !Particle has converged on face and flies outside
        boole_last_line_defense_converged = .true.
        boole_distance_bisection = .false.

        !Set all allowed faces to true 
        allowed_faces = .true.
!
        !Set to initial values
        tau = 0.d0
        z = z_init
        iface_new = iface_init
!
if(diag_pusher_tetra_rk) then
    print *, ''
    print *, 'Last line of defense'
    print *, 'iface_init', iface_init
    print *, 'normal_distances', normal_distances_func(z(1:3))
    !call rk4_step(z,0.d0,dzdtau)
    !print *, 'vnorm(iface_init)', normal_velocity_func(iface_init,dzdtau,z)
    print *, 'boole_accuracy_ode45', boole_accuracy_ode45
endif

!
        !Store, if initial normal distance to iface is negative (--> Is the particle inside or outside regarding the four planes?)
        iface_init_outside = 0     !Default value, if particle is inside regarding all faces
        if(iface_init.ne.0) then
            if(normal_distance_func(z(1:3),iface_init).lt.0.d0) iface_init_outside = iface_init
        endif
!
if(diag_pusher_tetra_rk) print *, 'Initial quadratic analytical approximation'
        !First quadratic analytical approximation
        call quad_analytic_approx(z,allowed_faces,iface_new,dtau,boole_quad_approx)
!
        !Analytical approximation exists
        if(boole_quad_approx) then
!
            !First RK4 step  
            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
            tau = tau + dtau
!            
        !No analytical approximation exists
        else  
if(diag_pusher_tetra_rk) print *, 'no analytical approximation exists'
!
            !If no analytical approximation exists, use dtau_ref as dtau
if(diag_pusher_tetra_rk) print *, 'dtau_ref', dtau_ref
            dtau = dtau_ref
!
            !First RK4 step 
            call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
            tau = tau + dtau
!
            !Guess iface_new: minimum normal distance
            iface_new = minloc(normal_distances_func(z(1:3)) ,1)
        endif
!
if(diag_pusher_tetra_rk) then
    print *, 'After analytical approximation - normal_distances'
    print *, normal_distances_func(z(1:3))
    print *, 'dtau',dtau
endif    
!
!
        !Extreme distance prevention Bisection: Push particle with bisection inside dist_max sphere
        k = 0
        do       
            k = k+1
            normal_distances = normal_distances_func(z(1:3))
            if(any(normal_distances.gt.dist_max)) then
                boole_distance_bisection = .true.
if(diag_pusher_tetra_rk)   print *, 'Extreme distance prevention Bisection is used'
                dtau = tau - 0.5d0*abs(dtau)
                tau = 0.d0
                z = z_init
!                
                call integration_step(z,dtau,dzdtau,.false.) !Set accuracy to FALSE, always!
                tau = tau + dtau
            else
                if(boole_distance_bisection) then
                    !Recompute the orbit position inside dist_max with higher accuracy
                    if(boole_accuracy_ode45) then
                        z = z_init
                        dtau = tau
                        call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                        boole_distance_bisection = .false.
                    else
                        exit
                    endif
                else
                    exit 
                endif
            endif

            if(k.gt.kiter) then
if(diag_pusher_tetra_rk)       print *,"'Error in last line of defense: Extreme distance prevention bisection. &
                                & No intersection within ",kiter, " RK-steps."
!                         call print_starting_condition() 
                boole_last_line_defense_converged = .false.
if(diag_pusher_tetra_rk)       print *, 'boole_accuracy_ode45', boole_accuracy_ode45
                return                             
            endif          
        enddo
!
        
!
        !Classification after analytical quadratic step
        if(iface_init_outside.ne.0) then
            if(normal_distance_func(z(1:3),iface_init_outside).lt.0.d0) then
                if(normal_velocity_func(iface_init_outside,dzdtau,z).lt.0.d0) then
                    !Check, if particle has crossed another plane
                    do i=1,3
                        j=modulo(iface_init_outside+i-1,4)+1
                        if(normal_distance_func(z(1:3),j).lt.0.d0) then
                            boole_turned_tangential = .true. !Particle has turned and also crossed another plane
                        endif
                    enddo
                    iface_new = iface_init_outside
                    if(abs(normal_distance_func(z(1:3),iface_new)).lt.dist_min) then
                        boole_converged = .true.
                    endif
                endif
            endif  
        endif
!
        !Treat all cases, that are not converged yet
        if(.not.boole_converged) then
!
            !Treat tangential cases (that have crossed another plane) with bisection convergence procedure
            if(boole_turned_tangential) then
!
if(diag_pusher_tetra_rk) print *, 'Turned tangential '
!
                call bisection_face_convergence(z,tau,dtau,iface_new,dzdtau,boole_accuracy_ode45,boole_bisection_converged)
!
                if(.not.boole_bisection_converged) then
if(diag_pusher_tetra_rk)   print *, 'Error in last line of defense: Bisection procedure fails. (1)'
!                     call print_starting_condition() 
                    boole_last_line_defense_converged = .false.
if(diag_pusher_tetra_rk)   print *, 'boole_accuracy_ode45', boole_accuracy_ode45
                    return
                endif
!
            !Treat all other cases (Non tangential)
            else
!
if(diag_pusher_tetra_rk)    print *, 'Inside special bisection:'
                !Special bisection procedure until particle is outside with respect to ONLY ONE plane and flies outside. 
                boole_dtau_decreased = .false.
!   

if(diag_pusher_tetra_rk) print *, 'dist_min',dist_min

if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(94,*) tau, normal_distances_func(z(1:3) )
                k = 0
                do 
                    k = k+1
!
                    normal_distances = normal_distances_func(z(1:3))
                    l_normal_distances = normal_distances.lt.0d0
                    i_normal_distances = logical2integer(l_normal_distances) !conversion from LOGICAL to INTEGER
                    n_faces_outside = sum(i_normal_distances,1)
!                    
                    if(n_faces_outside.eq.0) then
                        if(boole_dtau_decreased) then
                            dtau = 0.5d0*abs(dtau)
                        else
                            dtau = 2.d0*abs(dtau)  !increase time
                        endif  
                    elseif(n_faces_outside.eq.1) then
                        iface_new = MaxLoc(i_normal_distances,1)
                        if(normal_velocity_func(iface_new,dzdtau,z).ge.0.d0) then
                            if(iface_init_outside.ne.iface_new) then !Particle is outside another plane, than it was intially. Therefore, the particle must have crossed a plane.
                                dtau = - 0.5d0*abs(dtau) !decrease time
                                boole_dtau_decreased = .true.
                            else
                                if(boole_dtau_decreased) then
                                    dtau = 0.5d0*abs(dtau)
                                else
                                    dtau = 2.d0*abs(dtau)  !increase time
                                endif    
                            endif
                        else
                            !Exit loop: Particle is outside only ONE plane and flies outside
                            exit
                        endif
                    else !n_faces_outside >= 2
!                       
                        l = 0 !Counter for particles that are outside, but within convergence region
                        !Check, if particles, that are outside are within convergence region
                        do i = 1,4
                            if(.not.l_normal_distances(i)) cycle    !Skip particles, that are inside the tetrahedron
                            if(abs(normal_distances(i)).lt.dist_min) l = l +1
                        enddo
                        
                        !Distinguish treatment, whether all outside-particles are inside the convergence region or not
                        if(l.eq.n_faces_outside) then
                            j = 0
                            do i = 1,4
                                if(.not.l_normal_distances(i)) cycle    !Skip particles, that are inside the tetrahedron
                                if(abs(normal_distances(i)).ge.dist_min) cycle ! Skip particles, that are not inside the convergence region
                                if(normal_velocity_func(i,dzdtau,z).gt.0.d0) j = j +1 !Count particles, that fly inside
                            enddo
!                            
                            if(j.gt.0) then !At least one particle that is outside, but within convergence region, flies inside
                                dtau = 2.d0*abs(dtau)  !increase time
                            else !All particles that are outside, but within convergence region, fly outside
                                !Formally, this must be the solution. (But due to second order convergence problem requirements are relaxed.)
                                !dtau = -0.5d0*abs(dtau)  !decrease time
                                !boole_dtau_decreased = .true.
                                !
                                iface_new = minloc(abs(normal_distances),1,MASK = l_normal_distances) !Take the minimum of those particles, that are outside but within convergence region
if(diag_pusher_tetra_rk) print *, 'iface_new', iface_new, 'l_normal_distances',l_normal_distances, &
                            & 'normal_distances', normal_distances
                                !Exit loop: More than one particle are outside, but within convergence region and all fly outside
                                exit
                            endif
                        else
                            dtau = -0.5d0*abs(dtau)  !decrease time
                            boole_dtau_decreased = .true.
                        endif
                    endif
!                   

!
                    !Extreme distance prevention
                    if(any(normal_distances.gt.dist_max)) then
                        dtau = tau + dtau
                        tau = 0.d0
                        z = z_init
                    endif
!                    
                    !Runge-Kutta step   
                    call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
                    tau = tau + dtau                    
!
if(diag_pusher_tetra_rk) then
    print *, 'k = ',k
    print *, 'dtau',dtau
    print *, 'normal_distances'
    print *, normal_distances_func(z(1:3))
endif 

if(diag_pusher_tetra_rk.and.boole_accuracy_ode45) write(94,*) tau, normal_distances_func(z(1:3) )

                    if(k.gt.kiter) then
if(diag_pusher_tetra_rk)  print *,"'Error in last line of defense: Special bisection. No intersection within ",kiter, " RK-steps."
!                         call print_starting_condition() 
                        boole_last_line_defense_converged = .false.
if(diag_pusher_tetra_rk)  print *, 'boole_accuracy_ode45', boole_accuracy_ode45
                        return                             
                    endif
                enddo                
!


                !Save position and time after special bisection
                z_save = z
                dtau_save = dtau
                tau_save = tau

if(diag_pusher_tetra_rk) then
    !Negative tau prevention
    if(tau.le.0.d0) then
        print *, 'Error in last line of defense: tau is negative after special bisection'
        print *, 'tau',tau
        call print_starting_condition()
        stop
    endif
endif
!
if(diag_pusher_tetra_rk) print *, 'Before Newton - normal_distances', normal_distances_func(z(1:3)), 'tau',tau

                !Try to converge with linear or quadratic Newton
                call newton_face_convergence(z,tau,iface_new,dzdtau,boole_accuracy_ode45,boole_newton_converged)
!
if(diag_pusher_tetra_rk) print *, 'After Newton - normal_distances', normal_distances_func(z(1:3)), 'tau',tau
!
                !Three-planes control
                do i=1,3
                    j=modulo(iface_new+i-1,4)+1
                    if(normal_distance_func(z(1:3),j).lt.0.d0) then
                        boole_newton_converged = .false. !Particle has converged with Newton, but is outside another plane
if(diag_pusher_tetra_rk) print *, 'After Newton - Three planes'
                    endif
                enddo
!
if(diag_pusher_tetra_rk) print *, 'vnorm(iface)', normal_velocity_func(iface_new,dzdtau,z)
!

                !If linear or quadratic Newton fails, use bisection
                if((.not.boole_newton_converged).or.(normal_velocity_func(iface_new,dzdtau,z).ge.0.d0)) then
!
                    !Use position and time before Newton convergence procedure
                    z = z_save
                    tau = tau_save
                    dtau = dtau_save
!                    
                    call bisection_face_convergence(z,tau,dtau,iface_new,dzdtau,boole_accuracy_ode45,boole_bisection_converged)
!
                    if(.not.boole_bisection_converged) then
if(diag_pusher_tetra_rk)   print *, 'Error in last line of defense: Bisection procedure fails. (2)'
!                          call print_starting_condition() 
if(diag_pusher_tetra_rk)   print *, 'boole_accuracy_ode45', boole_accuracy_ode45
                        boole_last_line_defense_converged = .false.
                        return
                    endif                    
                endif 
!            
            endif !boole_turned_tangential
!        
        endif !not boole converged
!
        iface = iface_new
!        
if(diag_pusher_tetra_rk) then
!
        print *, 'LLOD finished: iface', iface
!
        ! Important testing feature: Is the particle outside the tetrahedron regarding more than one plane?
        k = 0
        normal_distances = normal_distances_func(z(1:3))
        call rk4_step(z,0.d0,dzdtau)
        do i = 1,4
            if(normal_distances(i).lt.0.d0) k = k+1
        enddo

        if(k.gt.1) then
            print *, 'Error: The particle is outside the tetrahedron regarding more than one plane.'
            print *, 'normal_distances',normal_distances
!             call print_starting_condition()
!             stop
        endif
!
        if(any(normal_distances.lt.-dist_min)) then
            print *, 'Error: Particle is outside the tetrahedron and outside the convergence region'
!             call print_starting_condition()
!             stop
        endif
!
        if(normal_velocity_func(iface,dzdtau,z).ge.0.d0) then
            print *, 'vnorm(iface)',iface
        endif
!        
        if(tau.le.0.d0) then
            print *, 'Tau is negative:', tau
        endif
endif      
!
    end subroutine last_line_defense
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine final_processing(z,tau,iface_inout,ind_tetr_out,iper_phi,x,vpar,t_pass,boole_t_finished, &
                               & boole_accuracy_ode45, boole_final_processing_converged)
!
        use tetra_physics_mod, only: tetra_physics,dt_dtau
        use gorilla_diag_mod, only: diag_pusher_tetra_rk
        use gorilla_settings_mod, only: boole_dt_dtau
        use pusher_tetra_func_mod, only: pusher_handover2neighbour
!
        implicit none
!
        double precision, dimension(4), intent(inout)           :: z
        double precision, intent(inout)                         :: tau
        integer, intent(inout)                                  :: iface_inout
        integer, intent(out)                                    :: ind_tetr_out,iper_phi
        double precision, dimension(3), intent(out)             :: x
        double precision, intent(out)                           :: vpar,t_pass
        logical, intent(out)                                    :: boole_t_finished, boole_final_processing_converged
        logical, intent(in)                                     :: boole_accuracy_ode45
        double precision, dimension(4)                          :: dzdtau,dzdtau_save
        double precision, dimension(4)                          :: z_save,normal_distances,normal_distances_save
        double precision, dimension(4)                          :: delta_normal_distances
        double precision                                        :: eps_time = 1.d-8
        double precision                                        :: dtau,tau_save,dtau_save
        integer                                                 :: i,j,k,iface_new,iface_outside
        integer                                                 :: i_outside_plane,iper_theta
        logical                                                 :: boole_newton_converged,boole_bisection_converged
        logical                                                 :: boole_last_line_defense_converged 
!        
        boole_final_processing_converged = .true.
!
        iface_new = iface_inout
!
        !Calculate final values
        x=z(1:3)+tetra_physics(ind_tetr)%x1
!
        !Calculate time to to pass the tetrahedron
        if(boole_dt_dtau) then
            t_pass = tau*dt_dtau_const
        else
            t_pass = tau*dt_dtau(ind_tetr,x_init,x)
        endif
!
!
if(diag_pusher_tetra_rk) then
    if(boole_accuracy_ode45) then
        print *, 'Start Final processing with ODE45'
    else
        print *, 'Start Final processing with RK4'
    endif
endif

if(diag_pusher_tetra_rk) print *, 't_pass',t_pass
if(diag_pusher_tetra_rk) print *, 't_remain',t_remain
! if(diag_pusher_tetra_rk) print *, 'x_init',x_init
! if(diag_pusher_tetra_rk) print *, 'x',x
! if(diag_pusher_tetra_rk) print *, 'tau',tau
!
        !Time handling - Treat case, if orbit stops within the cell
        if(t_remain.lt.t_pass) then
!
            !Start from the original coordinates
            z = z_init
            tau=0.d0
!
            !New guess for dtau
            if(boole_dt_dtau) then
                dtau = t_remain/dt_dtau_const
            else
                dtau = t_remain/dt_dtau(ind_tetr,x_init,x)! !x(1) is used, because it is the "best" guess
            endif
!

            !Iterative RK-procedure until t (not tau)! equals t_remain
            k = 0
            do
                k = k+1                
!
                !Save z for extreme distance prevention
                z_save = z
                normal_distances_save = normal_distances_func(z(1:3))
!
                call integration_step(z,dtau,dzdtau,boole_accuracy_ode45)
!
                !Extreme distance prevention (otherwise negative dt_dtau might occur)
                if(any(normal_distances_save.gt.dist_max)) then
!
                    z = z_save
!                    
                    !Use last line of defense solver 
                    call last_line_defense(z,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                        & boole_last_line_defense_converged)                  
!
                    !Validity check of LLOD regarding other planes
                    three_planes_loop: do j=1,3                             !Just consider the planes without the "exit-plane"
                        k=modulo(iface_new+j-1,4)+1
                        if(normal_distance_func(z(1:3),k).lt.0.d0) then     !If distance is negative, exitpoint of the considered plane is outside 
if(diag_pusher_tetra_rk)           print *, 'Error final processing: LLOD yields wrong result.'
                            boole_final_processing_converged = .false.
                            return
                        endif
                    enddo three_planes_loop
!
                    !Validity check of LLOD regarding normal velocity
                    if(normal_velocity_func(iface_new,dzdtau,z).gt.0.d0) then 
if(diag_pusher_tetra_rk)       print *, 'Error final processing: LLOD yields wrong result.'
                        boole_final_processing_converged = .false.
                        return
                    endif
!                    
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    if(boole_dt_dtau) then
                        t_pass = tau*dt_dtau_const
                    else
                        t_pass = tau*dt_dtau(ind_tetr,x_init,x)
                    endif
!                   
                    !Desired convergence point MUST be reached within t_remain
                    if(t_pass.le.t_remain) then
!                        
                        boole_t_finished = .false.
                        vpar=z(4)
!
                        iface_inout = iface_new
                        call pusher_handover2neighbour(ind_tetr,ind_tetr_out,iface_inout,x,iper_phi)
!                        
                        exit
                    else
if(diag_pusher_tetra_rk)       print *, 'Error final processing: LLOD does not yield desired convergence point.'
                        boole_final_processing_converged = .false.
                        return
                    endif
!                    
                endif
!
                tau = tau + dtau
                x=z(1:3)+tetra_physics(ind_tetr)%x1
                normal_distances = normal_distances_func(z(1:3))
!                
                if(boole_dt_dtau) then
                    t_pass = tau*dt_dtau_const
                    exit !If dt_dtau is a constant, loop should only be run once
                else
                    t_pass = tau*dt_dtau(ind_tetr,x_init,x)
                endif    
!                
                
                delta_normal_distances = abs(normal_distances - normal_distances_save)
!                
                if(Maxval(delta_normal_distances,1).lt.dist_min ) then
                    exit !particle is converged
                endif
!
                if(k.gt.kiter) then
if(diag_pusher_tetra_rk)   print *,"Error - Final processing: Can't calculate position that stops within cell for kiter steps."
!                     call print_starting_condition()              
!                     stop
                    boole_final_processing_converged = .false.
                    return
                endif
!
                !Calculate new dtau
                if(boole_dt_dtau) then
                    dtau = t_remain/dt_dtau_const-tau
                else
                    dtau = t_remain/dt_dtau(ind_tetr,x_init,x)-tau
                endif
            enddo
!
            !Check, if orbit time is finished, but particle is inside another tetrahedron.
            !This might happen close to the banana tip.
            i_outside_plane = 0
            iface_outside = 0
            do i = 1,4
                if(normal_distances(i).lt.0.d0) then 
                    i_outside_plane = i_outside_plane + 1
                    iface_outside = i
                endif  
            enddo
!     
            !Check, if orbit is converged at iface
            if(any(abs(normal_distances).lt.dist_min)) then
!
                if(i_outside_plane.eq.1) then !If particle is outside of face iface_outside of the tetrahedron
                    iface_new = iface_outside
                else
                    do i = 1,4
                        if(abs(normal_distances(i)).lt.dist_min) iface_new = i !find converged face iface_new
                    enddo
                endif
!          
                if(normal_velocity_func(iface_new,dzdtau,z).lt.0.d0) then
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    if(boole_dt_dtau) then
                        t_pass = tau*dt_dtau_const
                    else
                        t_pass = tau*dt_dtau(ind_tetr,x_init,x) 
                    endif
                    boole_t_finished = .true.
                    vpar=z(4)
!
                    iface_inout = iface_new
                    call pusher_handover2neighbour(ind_tetr,ind_tetr_out,iface_inout,x,iper_phi)
                else
                    !Particle flies inside, but is already converged at t_remain = 0 flight time inside the convergence region
                    !Particle is handed to the same tetrahedron once more for a new t.
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
                    if(boole_dt_dtau) then
                        t_pass = tau*dt_dtau_const
                    else
                        t_pass = tau*dt_dtau(ind_tetr,x_init,x) 
                    endif
                    boole_t_finished = .true.
                    vpar=z(4)
!                        
                    ind_tetr_out = ind_tetr
                    iface_inout = iface_new      
                endif
!        
            elseif(i_outside_plane.ne.0) then
!
                if(i_outside_plane.eq.1) then
!
                    !Try to converge particle with Newton1 and Newton2 at iface_new = m
                    iface_new = iface_outside
!                    
if(diag_pusher_tetra_rk) print *, 'iface outside', iface_new
!
                    !Save values before Newton's method
                    tau_save = tau
                    z_save = z
!
                    !Linear and quadratic Newton's method to converge at selected iface_new - Try immediately quadratic solution, due to the fact that this specific case mainly occurs with tangential trajectories
                    call newton_face_convergence(z,tau,iface_new,dzdtau,boole_accuracy_ode45,boole_newton_converged,.true.)
!
                    !Propose another Newton's method iteration, where also linear Newton is allowed
                    !This might be a case close to the banana tip
                    if(.not.boole_newton_converged) then
!                    
                        z = z_save
                        tau = tau_save
!
                        call rk4_step(z,0.d0,dzdtau)
                        call newton_face_convergence(z,tau,iface_new,dzdtau,boole_accuracy_ode45,boole_newton_converged)
!                                
                        if(.not.boole_newton_converged) then
                            z = z_save
                            tau = tau_save
if(diag_pusher_tetra_rk) print *, 'bisection (0)'
                            call bisection_face_convergence(z,tau,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                    & boole_bisection_converged)
!                            
                            if(.not.boole_bisection_converged) then
if(diag_pusher_tetra_rk)           print *,"Error in final processing. - Bisection convergence. (0)"
                                boole_final_processing_converged = .false.
                                return
                            endif 
                        endif
                    endif                  
!                    
                    if(normal_velocity_func(iface_new,dzdtau,z).gt.0.d0) then
!
                        z = z_save
                        tau = tau_save
if(diag_pusher_tetra_rk) print *, 'bisection (1)'
                        call bisection_face_convergence(z,tau,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                                & boole_bisection_converged)
!                            
                        if(.not.boole_bisection_converged) then
if(diag_pusher_tetra_rk)           print *,"Error in final processing. - Bisection convergence. (1)"
                            boole_final_processing_converged = .false.
                            return
                        endif                       
                    endif
!                            
                else !particle is outside the tetrahedron with respect to more than one plane
!                     print *, 'Error - Final processing: Particle is outside with respect to more than one plane.'
! 
if(diag_pusher_tetra_rk) print *, 'bisection (2)'
                    call bisection_face_convergence(z,tau,tau,iface_new,dzdtau,boole_accuracy_ode45, &
                            & boole_bisection_converged)
!                            
                    if(.not.boole_bisection_converged) then
if(diag_pusher_tetra_rk)               print *,"Error in final processing. - Bisection convergence. (2)"
!                                 call print_starting_condition
!                                 stop
                        boole_final_processing_converged = .false.
                        return
                    endif
!
                endif
!
                x=z(1:3)+tetra_physics(ind_tetr)%x1
                if(boole_dt_dtau) then
                    t_pass = tau*dt_dtau_const
                else
                    t_pass = tau*dt_dtau(ind_tetr,x_init,x) 
                endif
                boole_t_finished = .false.
                vpar=z(4)
!
                iface_inout = iface_new
                call pusher_handover2neighbour(ind_tetr,ind_tetr_out,iface_inout,x,iper_phi)
!
            !Orbit time is finished inside the tetrahedron
            else                  
!
                vpar=z(4)
                boole_t_finished = .true.
!                    
                ind_tetr_out = ind_tetr
                iface_inout = 0
                iper_phi = 0
                iper_theta = 0
            endif
! 
        !Continue final processing for completely passing cases
        else
!        
            boole_t_finished = .false.    
!
            vpar=z(4)
!
            iface_inout = iface_new
            call pusher_handover2neighbour(ind_tetr,ind_tetr_out,iface_inout,x,iper_phi)
!
        endif
!   
    end subroutine final_processing
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function bmod_func(z123)
!
        implicit none
!
        double precision :: bmod_func
        double precision, dimension(3),intent(in) :: z123
!
        bmod_func = B0+sum(gradB*z123)
!        
    end function bmod_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    function vperp_func(z123)
!    
        implicit none
!
        double precision :: vperp_func
        double precision, dimension(3),intent(in) :: z123

            if(perpinv.ne.0.d0) then
                vperp_func=sqrt(2.d0*abs(perpinv)*bmod_func(z123))
            else
                vperp_func = 0.d0
            endif
!
    end function vperp_func
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
        energy_tot_func = particle_mass/2.d0*(vperp_func(z(1:3))**2 + z(4)**2) + particle_charge*phi_elec_func(z(1:3))
!       
    end function energy_tot_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    function normal_distances_func(z123)
!
        implicit none
!
        double precision, dimension(3),intent(in)   :: z123
        double precision, dimension(4)              :: normal_distances_func
!        
        normal_distances_func=matmul(z123,anorm)
        normal_distances_func(1)=normal_distances_func(1)-dist1   !Correction for the one plane that is not lying in the first vertex
!        
    end function normal_distances_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!    
    function normal_distance_func(z123,iface)
!
        implicit none
!
        integer,intent(in)                          :: iface
        double precision, dimension(3),intent(in)   :: z123
        double precision                            :: normal_distance_func
!        
        normal_distance_func=sum(z123*anorm(:,iface))
        if(iface.eq.1) normal_distance_func=normal_distance_func-dist1   !Correction for the one plane that is not lying in the first vertex
!        
    end function normal_distance_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function normal_velocity_func(iface,dzdtau,z)
!
        use gorilla_settings_mod, only: boole_newton_precalc
!
        implicit none
!
        integer,intent(in)                              :: iface
        double precision, dimension(4),intent(in)       :: dzdtau,z
        double precision                                :: normal_velocity_func
!        
        if(boole_newton_precalc) then
            normal_velocity_func = normal_velocity_analytic(iface,z)
        else
            normal_velocity_func=sum(dzdtau(1:3)*anorm(:,iface))
        endif
!
    end function normal_velocity_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function normal_acceleration_func(iface,dzdtau,z)
!
        use gorilla_settings_mod, only: boole_newton_precalc
!
        implicit none
!
        integer,intent(in)                              :: iface
        double precision, dimension(4),intent(in)       :: dzdtau,z
        double precision                                :: normal_acceleration_func
!       
        if(boole_newton_precalc) then
            normal_acceleration_func = normal_acceleration_analytic(iface,z)
        else
            normal_acceleration_func = sum(matmul(anorm(:,iface),amat)*dzdtau(1:3))+sum(anorm(:,iface)*bvec)*dzdtau(4)
        endif
!
    end function normal_acceleration_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function normal_velocity_analytic(iface,z)
!        
        use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4
!
        implicit none
!        
        integer, intent(in)                         :: iface
        double precision, dimension(4), intent(in)  :: z
        double precision                            :: normal_velocity_analytic
!            
        normal_velocity_analytic = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,iface) + &
                                & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,iface)) * z) + &
                                & sum(anorm(:,iface) * b(1:3))
!        
    end function normal_velocity_analytic
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    function normal_acceleration_analytic(iface,z)
!        
        use tetra_physics_poly_precomp_mod, only: tetra_physics_poly4
!
        implicit none
!        
        integer, intent(in)                         :: iface
        double precision, dimension(4), intent(in)  :: z
        double precision                            :: normal_acceleration_analytic
!
        normal_acceleration_analytic = sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,iface) + &
                                        & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,iface) + &
                                        & perpinv2 * tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,iface)) * z) + &
!                    
                                        & sum((tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,iface) + &
                                        & perpinv * tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,iface)) * b)
!                                            
    end function normal_acceleration_analytic    
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine print_starting_condition()
!
        implicit none
!
        print *, 'Starting conditions'
        print *, 't_remain', t_remain
        print *, 'ind_tetr', ind_tetr
        print *, 'x_init', x_init
        print *, 'iface_init', iface_init
        print *, 'perpinv', perpinv
        print *, 'vperp_init', sqrt(2.d0*abs(perpinv)*bmod_func(z_init(1:3)))
        print *, 'vpar_init', z_init(4)
!        
    end subroutine print_starting_condition
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine find_tetra(x,vpar,vperp,ind_tetr_out,iface)
!
        use tetra_grid_mod, only : Rmin,Rmax,Zmin,Zmax,ntetr,tetra_grid, verts_sthetaphi
        use tetra_grid_settings_mod, only: grid_kind, grid_size, n_field_periods
        use tetra_physics_mod, only: tetra_physics,cm_over_e,isinside,coord_system
        use pusher_tetra_func_mod, only: pusher_handover2neighbour
        use constants, only: pi,clight,eps
        use supporting_functions_mod, only: logical2integer
!
        implicit none
!
        double precision, dimension(3), intent(inout) :: x
        double precision, intent(in) :: vpar,vperp
!
        integer, intent(out) :: ind_tetr_out,iface
!
        integer :: ir,iphi,iz,ind_search_tetra, indtetr_start,indtetr_end, ind_normdist, ind_tetr_save
        integer :: ind_plane_tetra_start, ntetr_in_plane, numerical_corr
        integer :: nr, nphi, nz
        integer :: iper_phi
        integer :: n_plane_conv, l, counter_vnorm_pos, iface_new, iface_new_save, i_tetra_try
        double precision ::hr,hphi,hz,vnorm,vperp2
        double precision, dimension(4):: cur_dist_value, z
        double precision, dimension(3) :: x_save
        double precision, dimension(4)   :: dzdtau
        logical, dimension(4) :: boole_plane_conv,boole_plane_conv_temp
        integer, dimension(:), allocatable :: ind_tetr_tried
!
        ! Initialize numerical correction
        numerical_corr = 0
!
        !Calculation of search domain depending on grid_kind
        select case(grid_kind)
            case(1) !rectangular grid
!
                nr = grid_size(1)
                nphi = grid_size(2)
                nz = grid_size(3)
!
                hr=(Rmax-Rmin)/nr !calculating the discretization step sizes in these directions
                hphi=(2.d0 * pi)/nphi
                hz=(Zmax-Zmin)/nz
!
                ir=int((x(1)-Rmin)/hr) +1 !get the reference coordinate in the (indexwise: 1-based) grid, so for instance the particle is in box (ir,iphi,iz)=(4,2,7)
                iphi=int(x(2)/hphi) +1
                iz=int((x(3)-Zmin)/hz) +1
!
                !Diagnostic to check, if out of domain
                if(ir.lt.1.or.ir.gt.nr.or.iphi.lt.1.or.iphi.gt.nphi.or.iz.lt.1.or.iz.gt.nz) then
                    ind_tetr_out=-1
                    iface=-1
                    print *, 'Error in start_tetra_rect: Particle is outside of the domain!'
                    stop
                endif
!
                indtetr_start = int((dble(iz)-1.d0)*6.d0 + 6.d0*dble(nz)*(dble(ir)-1.d0) & 
                & + 6.d0*(dble(iphi)-1.d0)*dble(nr)*dble(nz) +1.d0)
                    !This calculates the starting tetrahedron index as a function of (ir,iphi,iz)
                    !+1, since formula is 0-based (like istarttetr in make_grid_rect)
                indtetr_end = indtetr_start + 5
!
            case(2,4) !EFIT field-aligned grid or SOLEDGE3X_EIRENE
!
            select case(coord_system)
                case(1)
                    nphi = grid_size(2)
                    ind_plane_tetra_start = int(x(2)*nphi/(2.d0*pi/n_field_periods))
                    if(abs(x(2)*nphi/(2.d0*pi/n_field_periods) - dble(ind_plane_tetra_start)).gt.(1.d0-eps)) numerical_corr = 1
                case(2)
                    nphi = grid_size(2)
                    ind_plane_tetra_start = int(x(3)*nphi/(2.d0*pi/n_field_periods))
                    if(abs(x(3)*nphi/(2.d0*pi/n_field_periods) - dble(ind_plane_tetra_start)).gt.(1.d0-eps)) numerical_corr = 1
                end select
                ntetr_in_plane = ntetr/nphi !number of tetrahedra in a slice which is delimited by two phi=const. planes (with phi2-phi1 = 1*delta_phi)
                ! index of the slice, 0 if in on phi=0 plane
                indtetr_start = ind_plane_tetra_start*ntetr_in_plane +1
                indtetr_end = (ind_plane_tetra_start+1+numerical_corr)*ntetr_in_plane !num_corr ... 0 for normal cases, 1 for numerical computer precision rounding errors
!
            case(3) !VMEC field-aligned gird
                nphi = grid_size(2)
                ind_plane_tetra_start = int(x(3)*nphi/(2.d0*pi/n_field_periods))
                if(abs(x(3)*nphi/(2.d0*pi/n_field_periods) - dble(ind_plane_tetra_start)).gt.(1.d0-eps)) numerical_corr = 1
    !
                ntetr_in_plane = ntetr/nphi !number of tetrahedra in a slice which is delimited by two phi=const. planes (with phi2-phi1 = 1*delta_phi)
                ! index of the slice, 0 if in on phi=0 plane
                indtetr_start = ind_plane_tetra_start*ntetr_in_plane +1
                indtetr_end = (ind_plane_tetra_start+1+numerical_corr)*ntetr_in_plane !num_corr ... 0 for normal cases, 1 for numerical computer precision rounding errors
!
        end select
!
        indtetr_end = Minval([indtetr_end,ntetr])
        do ind_search_tetra = indtetr_start, indtetr_end
!
            if (isinside(ind_search_tetra,x)) then !inside tetrahedron
!
            ind_tetr_out = ind_search_tetra
            iface = 0
!
            do ind_normdist = 1,4 !calculate distances
                if (ind_normdist .ne. 1) then
                cur_dist_value(ind_normdist) = &
                    & sum(tetra_physics(ind_search_tetra)%anorm(:,ind_normdist)*(x-tetra_physics(ind_search_tetra)%x1))
                else ! ind_normdist .eq. 1
                cur_dist_value(ind_normdist) = sum(tetra_physics(ind_search_tetra)%anorm(:,ind_normdist)*&
                &(x-tetra_physics(ind_search_tetra)%x1))+tetra_physics(ind_search_tetra)%dist_ref
                endif
            enddo
!
            !Check, if particle is in vicinity of a plane (converged on plane)
            boole_plane_conv = abs(cur_dist_value) .le. (eps*abs(tetra_physics(ind_search_tetra)%dist_ref))
            n_plane_conv = sum(logical2integer(boole_plane_conv),1)
!
!print *, 'n_plane_conv',n_plane_conv
!
            if ( n_plane_conv.gt.0 ) then !if it is inside and next to a plane, set iface to the index of where it is 0
!
                !Initialize starting values and 'working' constants in module for tetrahedron
                !Squared perpendicular velocity
                vperp2 = vperp**2
!
                !Parallel velocity and position
                z(1:3) = x-tetra_physics(ind_tetr_out)%x1
                z(4) = vpar
!
                !Temporary iface
                iface_new = minloc(abs(cur_dist_value),1)
!
                ! Allocate and initialize vector with tried indices
                allocate(ind_tetr_tried(2*n_plane_conv))
                ind_tetr_tried = 0
!
                !Compute perpendicular invariant of particle
                perpinv=-0.5d0*vperp2/(tetra_physics(ind_tetr_out)%bmod1+sum(tetra_physics(ind_tetr_out)%gb*z(1:3)))
                call initialize_const_motion_rk(perpinv,perpinv**2)
!
                !Loop over all possible tetrahedra options
                try_loop: do i_tetra_try = 1,(2*n_plane_conv)
!
!print *, 'i_tetra_try',i_tetra_try
!print *, 'ind_tetr',ind_tetr_out
!print *, 'x',x
!
                    !Write 'tried' tetrahedron in 'memory'-vector
                    ind_tetr_tried(i_tetra_try) = ind_tetr_out
!
                    z(1:3) = x-tetra_physics(ind_tetr_out)%x1
!
                    call initialize_pusher_tetra_rk_mod(ind_tetr_out,x,iface_new,vpar,-1.d0)
!
                    cur_dist_value = normal_distances_func(z(1:3))
!print *, 'norm',cur_dist_value
!
                    !Check, if particle is inside the tetrahedron
                    if (any(cur_dist_value.lt.(-dist_min) )) then
                        print *, 'Error in find_tetra: Particle is not inside the tetrahedron while searching.'
                        stop
                    endif
!
                    boole_plane_conv_temp = abs(cur_dist_value) .le. (eps*abs(tetra_physics(ind_tetr_out)%dist_ref))
!
                    !call rk4 at z with dtau = 0 to get the velocity
                    call rk4_step(z,0.d0,dzdtau)
!
                    !Now, ALL normal velocites, where particle is converged  on plane, must be positive
                    counter_vnorm_pos = 0
                    do l = 1,4
                        if(.not.boole_plane_conv_temp(l)) cycle
                        vnorm = normal_velocity_func(l,dzdtau,z)
!print *, 'l',l,'norm(l)', normal_distance_func(z(1:3),l),'vnorm(l)',vnorm, 'boole',boole_plane_conv_temp(l)
                        if (vnorm .gt. 0.d0) then
                            counter_vnorm_pos = counter_vnorm_pos + 1
                        endif
                    enddo
!
!print *, 'counter_vnorm_pos',counter_vnorm_pos

                    if (counter_vnorm_pos.eq.n_plane_conv) then
                        iface = iface_new
!print *, 'gotcha'
                        exit try_loop!this means the tetrahedron was found
                    else
                        ind_tetr_save = ind_tetr_out
                        x_save = x
                        iface_new_save = iface_new
!
                        do l = 1,4
                            if(.not.boole_plane_conv_temp(l)) cycle
                            vnorm = normal_velocity_func(l,dzdtau,z)
                            if (vnorm.gt.0.d0) cycle
                            iface_new = l
                            call pusher_handover2neighbour(ind_tetr_save,ind_tetr_out,iface_new,x,iper_phi)
!
                            if( any(ind_tetr_out.eq.ind_tetr_tried)) then
!print *, 'Denied ind_tetr',ind_tetr_out
!print *, 'Denied iface_new', l
                                x = x_save
                                iface_new = iface_new_save
                            else
!print *, 'Succesful iface_new',l
                                exit !New try
                            endif
                        enddo
                        
                    endif

                enddo try_loop !l = i_tetra_try,n_exactly_conv
!
                deallocate(ind_tetr_tried)
!
            endif ! n_plane_conv.gt.0
!
            exit !this means the tetrahedron was found
!
            else !not inside
                ind_tetr_out = -1
                iface = -1
            endif !isinside
        enddo
!
        if (ind_tetr_out .eq. -1) then
            print *, ' in start_tetra: Starting tetrahedron was not found.'
            return
        endif
!
    end subroutine find_tetra
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine integration_step(z,dtau,dzdtau,boole_accuracy)
!
        use odeint_mod, only: odeint_allroutines
        use gorilla_settings_mod, only: rel_err_ode45
!
        implicit none
!
        double precision, dimension(4), intent(inout)       :: z
        double precision, intent(inout)                     :: dtau
        double precision, dimension(4) ,intent(out)         :: dzdtau
        logical, intent(in)                                 :: boole_accuracy
        integer                                             :: i,n_rk_steps,accuracy_option
        double precision                                    :: dtau_new
!
        accuracy_option = 1 ! 1 ... adaptive ODE45, 2 ... multiple RK4 steps with prescribed step-number
        n_rk_steps = 500
!
            if(boole_accuracy) then
                select case(accuracy_option)
                    case(1) !Adaptive ODE45
                        call odeint_allroutines(z,4,0.d0,dtau,rel_err_ode45,rhs_pusher_tetra_rk45)
                        call rk4_step(z,0.d0,dzdtau)
                    case(2) !RK4 with multiple steps
                        dtau_new = dtau/dble(n_rk_steps)
                        do i = 1,n_rk_steps
                            call rk4_step(z,dtau_new,dzdtau)
                        enddo
                end select
            else !Single RK4 step
                call rk4_step(z,dtau,dzdtau)
            endif    
!
    end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module pusher_tetra_rk_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module par_adiab_inv_rk_mod
!
    implicit none
!
    private
!
    integer,public,protected    :: counter_banana_mappings = 0
    double precision            :: par_adiab_inv
    integer                     :: nskip
!
    public :: par_adiab_inv_tetra_rk
!
    !$OMP THREADPRIVATE(counter_banana_mappings,par_adiab_inv)
!
    contains
!    
    subroutine par_adiab_inv_tetra_rk(t_pass,vpar_in,vpar_end,file_id_vpar_0,file_id_J_par,n_skip_vpar_0, &
                                     & boole_J_par,boole_poincare_vpar_0,boole_e_tot,file_id_e_tot)
!
        !Attention: This subroutine CAN ONLY BE CALLED directly after pusher, when modules still contain
        !           the quantities from last pushing
!        
        use pusher_tetra_rk_mod, only: dt_dtau_const,z_init,ind_tetr,energy_tot_func
        use tetra_physics_mod, only: particle_mass,tetra_physics
!  
        implicit none
!
        integer, intent(in)                 :: file_id_vpar_0,file_id_J_par,n_skip_vpar_0,file_id_e_tot
        double precision, intent(in)        :: t_pass,vpar_in,vpar_end
        logical, intent(in)                 :: boole_J_par,boole_poincare_vpar_0,boole_e_tot
        double precision                    :: tau,tau_part1,par_adiab_tau
        double precision, dimension(4)      :: z
        double precision, dimension(3)      :: x
!
        nskip = n_skip_vpar_0
!
        !Convert t_pass to tau
        tau = t_pass/dt_dtau_const
!
        z = z_init
!
        !Trace banana tips
        if((vpar_end.gt.0.d0).and.(vpar_in.lt.0.d0)) then
!
            !Compute integral of $v_\parallel^2$ and position z until $vpar = 0$
            call calc_par_adiab_until_root(tau,z,par_adiab_tau,tau_part1)
            par_adiab_inv = par_adiab_inv + par_adiab_tau*dt_dtau_const          
!
            if(counter_banana_mappings.gt.1) then
                if(counter_banana_mappings/nskip*nskip.eq.counter_banana_mappings) then
                    !$omp critical
                        if(boole_J_par) then
                            write(file_id_J_par,*) counter_banana_mappings,par_adiab_inv
                        endif
                    !$omp end critical
                    !Poloidal projection of Poincaré sections at $v=\parallel$ = 0
                    x=z(1:3)+tetra_physics(ind_tetr)%x1
!
                    !$omp critical
                        if(boole_poincare_vpar_0) then
                            write(file_id_vpar_0,*) x    ! write coordinates for poincare cuts
                        endif
                    !$omp end critical
!
                    !$omp critical
                        if(boole_e_tot) then
                            write(file_id_e_tot,*) counter_banana_mappings,energy_tot_func(z)
                        endif
                    !$omp end critical
                endif
            endif
            counter_banana_mappings = counter_banana_mappings + 1
!
            !Start to integrate par_adiab_inv for new bounce period
            par_adiab_inv = 0.d0
            call calc_par_adiab_tau(tau-tau_part1,z,par_adiab_tau)
            par_adiab_inv = par_adiab_inv + par_adiab_tau*dt_dtau_const  
        else    
            !Compute parallel adiabatic invariant as a function of time
            call calc_par_adiab_tau(tau,z,par_adiab_tau)
            par_adiab_inv = par_adiab_inv + par_adiab_tau*dt_dtau_const    
        endif     
!  
    end subroutine par_adiab_inv_tetra_rk
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!   
    subroutine calc_par_adiab_tau(dtau,z_inout,par_adiab_tau)
!
        use odeint_mod, only: odeint_allroutines  
        use gorilla_settings_mod, only: rel_err_ode45
!
        implicit none
!
        double precision                :: dtau,par_adiab_tau
        double precision, dimension(4)  :: z_inout
        double precision, dimension(5)  :: z
!
        !z(5) = par_adiab_inv
        z(1:4) = z_inout
        z(5) = 0
!
        call odeint_allroutines(z,5,0.d0,dtau,rel_err_ode45,rhs_par_adiab_ode45)
!
        par_adiab_tau = z(5)
        z_inout = z(1:4)
!        
    end subroutine calc_par_adiab_tau
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine calc_par_adiab_until_root(tau_in,z_inout,par_adiab_tau,tau_out)
!
        use odeint_mod, only: odeint_allroutines  
        use gorilla_settings_mod, only: rel_err_ode45
!
        implicit none
!
        double precision                :: dtau,par_adiab_tau,tau_in,tau_out
        double precision, dimension(4)  :: z_inout
        double precision, dimension(5)  :: z
        double precision                :: vpar_save,signum_dtau,vpar_min = 1.d1
        integer                         :: i,n_iter
!
        !z(5) = par_adiab_inv
        z(1:4) = z_inout
        z(5) = 0
!
        dtau = tau_in
        tau_out = 0.d0
!
        n_iter = 100
        i = 0
        do while(abs(z(4)).gt.vpar_min)
            i = i+1
!         
            vpar_save = z(4)
!
            call odeint_allroutines(z,5,0.d0,dtau,rel_err_ode45,rhs_par_adiab_ode45)
            tau_out = tau_out+dtau
!            
            !vpar_save > 0
            if(vpar_save.gt.0.d0) then
                !vpar > 0
                if(z(4).gt.0.d0) then
                    if(dtau.gt.0.d0) then
                        dtau = abs(dtau/2)
                    else
                        dtau = -abs(dtau/2)
                    endif
                !vpar < 0
                else
                    if(dtau.gt.0.d0) then
                        dtau = -abs(dtau/2)
                    else
                        dtau = +abs(dtau/2)
                    endif
                endif
            !vpar_save < 0
            else 
                !vpar > 0
                if(z(4).gt.0.d0) then
                    if(dtau.gt.0.d0) then
                        dtau = -abs(dtau/2)
                    else
                        dtau = +abs(dtau/2)
                    endif
                !vpar < 0
                else
                    if(dtau.gt.0.d0) then
                        dtau = +abs(dtau/2)
                    else
                        dtau = -abs(dtau/2)
                    endif
                endif
            endif
!            
            if(i.gt.n_iter) then
                print *, 'Root of vpar could not be found within ', n_iter, 'steps'
                stop
            endif
        enddo
!
        par_adiab_tau = z(5)
        z_inout = z(1:4)
!        
    end subroutine calc_par_adiab_until_root
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
    subroutine rhs_par_adiab_ode45(dummy,z,dzdtau)
!
        use pusher_tetra_rk_mod, only:b,amat,Bvec,spamat
!
        implicit none
!
        double precision, dimension(5) :: z,dzdtau
        double precision :: dummy
!
        dzdtau(1:3) = b(1:3)+matmul(amat,z(1:3))+Bvec*z(4)
        dzdtau(4) = b(4)+spamat*z(4)
        !Differential equation for integral over $v_\parallel^2$
        dzdtau(5) = z(4)**2
!
    end subroutine rhs_par_adiab_ode45
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! 
end module par_adiab_inv_rk_mod
