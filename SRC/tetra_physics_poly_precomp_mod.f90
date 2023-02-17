!
module tetra_physics_poly_precomp_mod
!
    implicit none
!
    private
!
    type tetrahedron_physics_precomp_poly1
        sequence
!
        double precision, dimension(4)              :: anorm_in_betvec, anorm_in_gamvec
!
        double precision, dimension(3,4)            :: anorm_in_alpmat,anorm_in_betmat,anorm_in_gammat
!
    end type
!
    type(tetrahedron_physics_precomp_poly1), dimension(:),   allocatable, public, protected :: tetra_physics_poly1
!
!--------------------------------------------------------------------------------------------------------------------------------!
!
    type tetrahedron_physics_precomp_poly4
        sequence
!
        double precision, dimension(4,4)            :: amat1_0,amat1_1
        double precision, dimension(4,4)            :: amat2_0,amat2_1,amat2_2
        double precision, dimension(4,4)            :: amat3_0,amat3_1,amat3_2,amat3_3
        double precision, dimension(4,4)            :: amat4_0,amat4_1,amat4_2,amat4_3,amat4_4
!
        double precision, dimension(4,4)            :: anorm_in_amat1_0,anorm_in_amat1_1
        double precision, dimension(4,4)            :: anorm_in_amat2_0,anorm_in_amat2_1,anorm_in_amat2_2
        double precision, dimension(4,4)            :: anorm_in_amat3_0,anorm_in_amat3_1,anorm_in_amat3_2,anorm_in_amat3_3
        double precision, dimension(4,4)            :: anorm_in_amat4_0,anorm_in_amat4_1,anorm_in_amat4_2,anorm_in_amat4_3
        double precision, dimension(4,4)            :: anorm_in_amat4_4
!
        double precision, dimension(4)              :: b0,b1,b2,b3
        double precision, dimension(4)              :: amat1_0_in_b0,amat1_0_in_b1,amat1_0_in_b2,amat1_0_in_b3
        double precision, dimension(4)              :: amat1_1_in_b0,amat1_1_in_b1,amat1_1_in_b2,amat1_1_in_b3
!
        double precision, dimension(4)              :: anorm_in_b0,anorm_in_b1,anorm_in_b2,anorm_in_b3
        double precision, dimension(4)              :: anorm_in_amat1_0_in_b0,anorm_in_amat1_0_in_b1
        double precision, dimension(4)              :: anorm_in_amat1_0_in_b2,anorm_in_amat1_0_in_b3
        double precision, dimension(4)              :: anorm_in_amat1_1_in_b0,anorm_in_amat1_1_in_b1
        double precision, dimension(4)              :: anorm_in_amat1_1_in_b2,anorm_in_amat1_1_in_b3

    end type
!   
    type(tetrahedron_physics_precomp_poly4), dimension(:),   allocatable, public, protected :: tetra_physics_poly4
!
!--------------------------------------------------------------------------------------------------------------------------------!
!
    type tetrahedron_physics_precomp_poly_perpinv
        sequence
!
        double precision,dimension(4)   :: bcoef_pre_k0,acoef_pre_k0,acoef_pre_k1,acoef_pre_k3
        double precision,dimension(4)   :: opb_pre_k0,opb_pre_k1,opb_pre_k3
        double precision,dimension(4,4) :: bcoef_pre_z,acoef_pre_z,opz_pre_tau,opz_pre_tau2
!
    end type
!
    type(tetrahedron_physics_precomp_poly_perpinv), dimension(:),allocatable, public, protected :: tetra_physics_poly_perpinv
    logical, dimension(:),allocatable :: boole_precomp_poly_perpinv
!
!--------------------------------------------------------------------------------------------------------------------------------!
!
    !$OMP THREADPRIVATE(tetra_physics_poly_perpinv,boole_precomp_poly_perpinv)
!
    public :: make_precomp_poly,make_precomp_poly_perpinv, &
            & initialize_boole_precomp_poly_perpinv,alloc_precomp_poly_perpinv
!    
    contains
!
        subroutine make_precomp_poly
!
            use gorilla_settings_mod, only: i_precomp,ipusher,boole_newton_precalc
!
            implicit none
!
            select case(ipusher)
                case(1) !numerical RK pusher
!
                    if(boole_newton_precalc) then
                        call make_precomp_poly4()
                    endif
!
                case(2) !polynomial pusher
!
                    if(i_precomp.eq.0) then
                        call make_precomp_poly1()
                    else
                        call make_precomp_poly4()
                    endif
!
            end select !ipusher
!
        end subroutine make_precomp_poly
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine make_precomp_poly1()
!
            use tetra_grid_mod, only: ntetr
            use tetra_physics_mod, only: tetra_physics
            use gorilla_settings_mod, only: boole_strong_electric_fields
!
            implicit none
!
            integer                                     :: ind_tetr,n
            double precision, dimension(3)              :: n_vec
!
            allocate(tetra_physics_poly1(ntetr))
!
            !$OMP PARALLEL &
            !$OMP& DO DEFAULT(NONE) &
            !$OMP& SHARED(ntetr,tetra_physics,tetra_physics_poly1,boole_strong_electric_fields) &
            !$OMP& PRIVATE(ind_tetr,n,n_vec)

            !Loop over all tetrahedra
            do ind_tetr = 1,ntetr
!
                !Loop over all normal vectors
                do n = 1,4
!
                    n_vec(1) = tetra_physics(ind_tetr)%anorm(1,n)
                    n_vec(2) = tetra_physics(ind_tetr)%anorm(2,n)
                    n_vec(3) = tetra_physics(ind_tetr)%anorm(3,n)
!
                    !anorm in amat:
!
                    tetra_physics_poly1(ind_tetr)%anorm_in_alpmat(:,n) = &
                                & matmul(n_vec,tetra_physics(ind_tetr)%alpmat)
!
                    tetra_physics_poly1(ind_tetr)%anorm_in_betmat(:,n) = &
                                & matmul(n_vec,tetra_physics(ind_tetr)%betmat)
!
                    tetra_physics_poly1(ind_tetr)%anorm_in_betvec(n) = &
                                & sum(n_vec*tetra_physics(ind_tetr)%curlA)

                    if (boole_strong_electric_fields) then
                        tetra_physics_poly1(ind_tetr)%anorm_in_gammat(:,n) = &
                                & matmul(n_vec,tetra_physics(ind_tetr)%gammat)
                        tetra_physics_poly1(ind_tetr)%anorm_in_gamvec(n) = &
                                & sum(n_vec*tetra_physics(ind_tetr)%curlvE)
                    endif
!
                enddo !(n=1,4)
!
            enddo !(ind_tetr = 1,ntetr)
            !$OMP END PARALLEL DO
!
            print *, 'Precomputation of coefficients for polynomial solution finished.'
!
        end subroutine make_precomp_poly1
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine make_precomp_poly4()
!
            use tetra_grid_mod, only: ntetr
            use tetra_physics_mod, only: tetra_physics, cm_over_e
            use constants, only: pi,clight
!
            implicit none
!
            integer         :: ind_tetr,n
            double precision, dimension(3,3)            :: alpmat
            double precision, dimension(3,3)            :: betmat
            double precision, dimension(4,4)            :: alp,bet,alp_alp,bet_bet,alp_bet,bet_alp
            double precision, dimension(4,4)            :: alp_alp_alp,alp_alp_bet,alp_bet_alp,alp_bet_bet
            double precision, dimension(4,4)            :: bet_alp_alp,bet_alp_bet,bet_bet_alp,bet_bet_bet 
            double precision, dimension(4,4)            :: alp_alp_alp_alp,alp_alp_alp_bet,alp_alp_bet_alp,alp_alp_bet_bet
            double precision, dimension(4,4)            :: alp_bet_alp_alp,alp_bet_alp_bet,alp_bet_bet_alp,alp_bet_bet_bet 
            double precision, dimension(4,4)            :: bet_alp_alp_alp,bet_alp_alp_bet,bet_alp_bet_alp,bet_alp_bet_bet
            double precision, dimension(4,4)            :: bet_bet_alp_alp,bet_bet_alp_bet,bet_bet_bet_alp,bet_bet_bet_bet 
            
            double precision                            :: spalpmat
            double precision                            :: spbetmat
            double precision, dimension(3)              :: bvec
            double precision, dimension(4)              :: n_vec            
!
            allocate(tetra_physics_poly4(ntetr))
!
            !$OMP PARALLEL &
            !$OMP& DO DEFAULT(NONE) &
            !$OMP& SHARED(ntetr,tetra_physics,tetra_physics_poly4,cm_over_e) &
            !$OMP& PRIVATE(ind_tetr,alpmat,spalpmat,betmat,spbetmat,bvec,n,n_vec, &
            !$OMP& alp,bet,alp_alp,bet_bet,alp_bet,bet_alp,alp_alp_alp,alp_alp_bet,alp_bet_alp, &
            !$OMP& alp_bet_bet,bet_alp_alp,bet_alp_bet,bet_bet_alp,bet_bet_bet, &
            !$OMP& alp_alp_alp_alp,alp_alp_alp_bet,alp_alp_bet_alp,alp_alp_bet_bet, &
            !$OMP& alp_bet_alp_alp,alp_bet_alp_bet,alp_bet_bet_alp,alp_bet_bet_bet, & 
            !$OMP& bet_alp_alp_alp,bet_alp_alp_bet,bet_alp_bet_alp,bet_alp_bet_bet, &
            !$OMP& bet_bet_alp_alp,bet_bet_alp_bet,bet_bet_bet_alp,bet_bet_bet_bet)

            !Loop over all tetrahedra
            do ind_tetr = 1,ntetr
!
                alpmat = cm_over_e*tetra_physics(ind_tetr)%alpmat
                spalpmat = cm_over_e*tetra_physics(ind_tetr)%spalpmat
                betmat = clight*tetra_physics(ind_tetr)%betmat
                spbetmat = clight*tetra_physics(ind_tetr)%spbetmat
                bvec = tetra_physics(ind_tetr)%curla
!
                !Define 4x4 matrices alpha, beta and powers
!
                !alpha:
                alp(1,:) = [alpmat(1,1), alpmat(1,2),alpmat(1,3),0.d0]
                alp(2,:) = [alpmat(2,1), alpmat(2,2),alpmat(2,3),0.d0]
                alp(3,:) = [alpmat(3,1), alpmat(3,2),alpmat(3,3),0.d0]
                alp(4,:) = [0.d0, 0.d0, 0.d0, spalpmat]
!
                !beta:
                bet(1,:) = [-betmat(1,1), -betmat(1,2), -betmat(1,3), bvec(1)]
                bet(2,:) = [-betmat(2,1), -betmat(2,2), -betmat(2,3), bvec(2)]
                bet(3,:) = [-betmat(3,1), -betmat(3,2), -betmat(3,3), bvec(3)]
                bet(4,:) = [0.d0, 0.d0, 0.d0, -spbetmat]
!
                !alpha in alpha
                alp_alp = matmul(alp,alp)
!
                !alpha in beta
                alp_bet = matmul(alp,bet)
!
                !beta in beta
                bet_bet = matmul(bet,bet)
!
                !beta in alpha
                bet_alp = matmul(bet,alp)
!
                !alpha in alpha in alpha
                alp_alp_alp = matmul(alp,alp_alp)
!
                !alpha in alpha in beta
                alp_alp_bet = matmul(alp,alp_bet)
!
                !alpha in beta in alpha
                alp_bet_alp = matmul(alp,bet_alp)
!
                !alpha in beta in beta
                alp_bet_bet = matmul(alp,bet_bet)        
!
                !beta in alpha in alpha
                bet_alp_alp = matmul(bet,alp_alp)
!
                !beta in alpha in beta
                bet_alp_bet = matmul(bet,alp_bet)
!
                !beta in beta in alpha
                bet_bet_alp = matmul(bet,bet_alp)
!
                !beta in beta in beta
                bet_bet_bet = matmul(bet,bet_bet)
!
                !alpha in alpha in alpha in alpha
                alp_alp_alp_alp = matmul(alp,alp_alp_alp)
!
                !alpha in alpha in alpha in beta
                alp_alp_alp_bet = matmul(alp,alp_alp_bet)
!
                !alpha in alpha in beta in alpha
                alp_alp_bet_alp = matmul(alp,alp_bet_alp)
!
                !alpha in alpha in beta in beta
                alp_alp_bet_bet = matmul(alp,alp_bet_bet)
!
                !alpha in beta in alpha in alpha
                alp_bet_alp_alp = matmul(alp,bet_alp_alp)
!
                !alpha in beta in alpha in beta
                alp_bet_alp_bet = matmul(alp,bet_alp_bet)
!
                !alpha in beta in beta in alpha
                alp_bet_bet_alp = matmul(alp,bet_bet_alp)
!
                !alpha in beta in beta in beta
                alp_bet_bet_bet = matmul(alp,bet_bet_bet)
!                
                !beta in alpha in alpha in alpha
                bet_alp_alp_alp = matmul(bet,alp_alp_alp)
!
                !beta in alpha in alpha in beta
                bet_alp_alp_bet = matmul(bet,alp_alp_bet)
!
                !beta in alpha in beta in alpha
                bet_alp_bet_alp = matmul(bet,alp_bet_alp)
!
                !beta in alpha in beta in beta
                bet_alp_bet_bet = matmul(bet,alp_bet_bet)
!
                !beta in beta in alpha in alpha
                bet_bet_alp_alp = matmul(bet,bet_alp_alp)
!
                !beta in beta in alpha in beta
                bet_bet_alp_bet = matmul(bet,bet_alp_bet)
!
                !beta in beta in beta in alpha
                bet_bet_bet_alp = matmul(bet,bet_bet_alp)
!
                !beta in beta in beta in beta
                bet_bet_bet_bet = matmul(bet,bet_bet_bet)
!
                !Compute a^k matrices factorized for initial condition c1^n
!
                !a^1 (c1^0)
                tetra_physics_poly4(ind_tetr)%amat1_0 = bet
!                
                !a^1 (c1^1)
                tetra_physics_poly4(ind_tetr)%amat1_1 = alp
!
                !a^2 (c1^0)
                tetra_physics_poly4(ind_tetr)%amat2_0 = bet_bet
!
                !a^2 (c1^1)
                tetra_physics_poly4(ind_tetr)%amat2_1 = bet_alp + alp_bet
!
                !a^2 (c1^2)
                tetra_physics_poly4(ind_tetr)%amat2_2 = alp_alp
!
                !a^3 (c1^0)
                tetra_physics_poly4(ind_tetr)%amat3_0 = bet_bet_bet
!
                !a^3 (c1^1)
                tetra_physics_poly4(ind_tetr)%amat3_1 = bet_bet_alp + bet_alp_bet + alp_bet_bet
!
                !a^3 (c1^2)
                tetra_physics_poly4(ind_tetr)%amat3_2 = bet_alp_alp + alp_bet_alp + alp_alp_bet
!                
                !a^3 (c1^3)
                tetra_physics_poly4(ind_tetr)%amat3_3 = alp_alp_alp
!                
                !a^4 (c1^0)
                tetra_physics_poly4(ind_tetr)%amat4_0 = bet_bet_bet_bet
!
                !a^4 (c1^1)
                tetra_physics_poly4(ind_tetr)%amat4_1 = bet_bet_bet_alp + bet_bet_alp_bet + bet_alp_bet_bet + alp_bet_bet_bet
!
                !a^4 (c1^2)
                tetra_physics_poly4(ind_tetr)%amat4_2 = bet_bet_alp_alp + bet_alp_bet_alp + bet_alp_alp_bet + &
                                                         & alp_bet_bet_alp + alp_bet_alp_bet + alp_alp_bet_bet
!                
                !a^4 (c1^3)
                tetra_physics_poly4(ind_tetr)%amat4_3 = alp_bet_alp_alp + alp_alp_bet_alp + alp_alp_alp_bet + bet_alp_alp_alp
!                
                !a^4 (c1^4)
                tetra_physics_poly4(ind_tetr)%amat4_4 = alp_alp_alp_alp
!
                !Compute n in a^k
!
                !Loop over all normal vectors
                do n = 1,4
!
                    n_vec(1) = tetra_physics(ind_tetr)%anorm(1,n)
                    n_vec(2) = tetra_physics(ind_tetr)%anorm(2,n)
                    n_vec(3) = tetra_physics(ind_tetr)%anorm(3,n)
                    n_vec(4) = 0.d0
!                
                    !anorm in amat:
!
                    !c1^0
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat1_0)
!
                    !c1^1 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat1_1)
!
                    !anorm in amat^2:
!
                    !c1^0
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat2_0(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat2_0)
!
                    !c1^1 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat2_1(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat2_1)
!
                    !c1^2 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat2_2(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat2_2)
!
                    !anorm in amat^3:
!
                    !c1^0
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat3_0(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat3_0)
!
                    !c1^1 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat3_1(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat3_1)
!
                    !c1^2 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat3_2(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat3_2)
!
                    !c1^3 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat3_3(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat3_3)
!
                    !anorm in amat^4:
!
                    !c1^0
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat4_0(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat4_0)
!
                    !c1^1 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat4_1(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat4_1)
!
                    !c1^2 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat4_2(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat4_2)
!
                    !c1^3 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat4_3(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat4_3)
!
                    !c1^4 
                    tetra_physics_poly4(ind_tetr)%anorm_in_amat4_4(:,n) = matmul(n_vec,tetra_physics_poly4(ind_tetr)%amat4_4)
                enddo !(n=1,4)
!                
                !Precomputation of factorized b-vector
                tetra_physics_poly4(ind_tetr)%b0(1:3) = -clight*tetra_physics(ind_tetr)%gPhixh1
                tetra_physics_poly4(ind_tetr)%b0(4) = -clight/cm_over_e*tetra_physics(ind_tetr)%gPhixcurlA
!
                tetra_physics_poly4(ind_tetr)%b1(1:3) = cm_over_e*tetra_physics(ind_tetr)%curlh
                tetra_physics_poly4(ind_tetr)%b1(4) = 0.d0
!
                tetra_physics_poly4(ind_tetr)%b2(1:3) = cm_over_e*tetra_physics(ind_tetr)%gBxh1
                tetra_physics_poly4(ind_tetr)%b2(4) = tetra_physics(ind_tetr)%gBxcurlA
!
                tetra_physics_poly4(ind_tetr)%b3(1:3) = - 2.d0*clight*tetra_physics(ind_tetr)%curlh
                tetra_physics_poly4(ind_tetr)%b3(4) = 0.d0
!
                !Precomputation of amat1 in b-vectors
                tetra_physics_poly4(ind_tetr)%amat1_0_in_b0 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_0,tetra_physics_poly4(ind_tetr)%b0)
                tetra_physics_poly4(ind_tetr)%amat1_0_in_b1 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_0,tetra_physics_poly4(ind_tetr)%b1)
                tetra_physics_poly4(ind_tetr)%amat1_0_in_b2 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_0,tetra_physics_poly4(ind_tetr)%b2)
                tetra_physics_poly4(ind_tetr)%amat1_0_in_b3 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_0,tetra_physics_poly4(ind_tetr)%b3)
                tetra_physics_poly4(ind_tetr)%amat1_1_in_b0 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_1,tetra_physics_poly4(ind_tetr)%b0)
                tetra_physics_poly4(ind_tetr)%amat1_1_in_b1 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_1,tetra_physics_poly4(ind_tetr)%b1)
                tetra_physics_poly4(ind_tetr)%amat1_1_in_b2 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_1,tetra_physics_poly4(ind_tetr)%b2)
                tetra_physics_poly4(ind_tetr)%amat1_1_in_b3 = &
                    & matmul(tetra_physics_poly4(ind_tetr)%amat1_1,tetra_physics_poly4(ind_tetr)%b3)
!             
                !Precomputation of anorm in b-vectors
                tetra_physics_poly4(ind_tetr)%anorm_in_b0 = [(sum(tetra_physics(ind_tetr)%anorm(:,n) &
                    & * tetra_physics_poly4(ind_tetr)%b0(1:3)),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_b1 = [(sum(tetra_physics(ind_tetr)%anorm(:,n) &
                    & * tetra_physics_poly4(ind_tetr)%b1(1:3)),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_b2 = [(sum(tetra_physics(ind_tetr)%anorm(:,n) &
                    & * tetra_physics_poly4(ind_tetr)%b2(1:3)),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_b3 = [(sum(tetra_physics(ind_tetr)%anorm(:,n) &
                    & * tetra_physics_poly4(ind_tetr)%b3(1:3)),n=1,4)]
!
                !Precomputation of anorm_in_amat1_0 in b-vectors
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b0 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n)*tetra_physics_poly4(ind_tetr)%b0),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b1 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n)*tetra_physics_poly4(ind_tetr)%b1),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b2 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n)*tetra_physics_poly4(ind_tetr)%b2),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0_in_b3 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_0(:,n)*tetra_physics_poly4(ind_tetr)%b3),n=1,4)]
!
                !Precomputation of anorm_in_amat1_1 in b-vectors
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b0 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)*tetra_physics_poly4(ind_tetr)%b0),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b1 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)*tetra_physics_poly4(ind_tetr)%b1),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b2 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)*tetra_physics_poly4(ind_tetr)%b2),n=1,4)]
                tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1_in_b3 = [( &
                    & sum(tetra_physics_poly4(ind_tetr)%anorm_in_amat1_1(:,n)*tetra_physics_poly4(ind_tetr)%b3),n=1,4)]
!                
            enddo !(ind_tetr = 1,ntetr)
            !$OMP END PARALLEL DO
!        
            print *, 'Precomputation of coefficients for polynomial solution finished.'
!
        end subroutine make_precomp_poly4
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine alloc_precomp_poly_perpinv(i_alloc,ntetr)
!
            implicit none
!
            integer,intent(in) :: i_alloc,ntetr
!
            if(i_alloc.eq.1) then
                allocate(tetra_physics_poly_perpinv(ntetr),boole_precomp_poly_perpinv(ntetr))
            else
                deallocate(tetra_physics_poly_perpinv,boole_precomp_poly_perpinv)
            endif
!
        end subroutine alloc_precomp_poly_perpinv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine initialize_boole_precomp_poly_perpinv()
!
            implicit none
!
            boole_precomp_poly_perpinv = .false.
!
        end subroutine initialize_boole_precomp_poly_perpinv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
        subroutine make_precomp_poly_perpinv(perpinv,perpinv2)
!
            use tetra_grid_mod, only: ntetr
!
            implicit none
!
            double precision, intent(in) :: perpinv,perpinv2
!
            integer                 :: i
!
            do i = 1,ntetr
                tetra_physics_poly_perpinv(i)%bcoef_pre_k0 = &
                    & tetra_physics_poly4(i)%anorm_in_b0(:) + perpinv*tetra_physics_poly4(i)%anorm_in_b2(:)
                tetra_physics_poly_perpinv(i)%acoef_pre_k0 = &
                    & tetra_physics_poly4(i)%anorm_in_amat1_0_in_b0(:) + &
                    & perpinv*tetra_physics_poly4(i)%anorm_in_amat1_1_in_b0(:) + &
                    & perpinv * (tetra_physics_poly4(i)%anorm_in_amat1_0_in_b2(:) + &
                    & perpinv*tetra_physics_poly4(i)%anorm_in_amat1_1_in_b2(:))
                tetra_physics_poly_perpinv(i)%acoef_pre_k1 = &
                    & tetra_physics_poly4(i)%anorm_in_amat1_0_in_b1(:) + &
                    & perpinv*tetra_physics_poly4(i)%anorm_in_amat1_1_in_b1(:)
                tetra_physics_poly_perpinv(i)%acoef_pre_k3 = &
                    & tetra_physics_poly4(i)%anorm_in_amat1_0_in_b3(:) + &
                    & perpinv*tetra_physics_poly4(i)%anorm_in_amat1_1_in_b3(:)
                tetra_physics_poly_perpinv(i)%bcoef_pre_z = &
                    & tetra_physics_poly4(i)%anorm_in_amat1_0(:,:) + &
                    & perpinv * tetra_physics_poly4(i)%anorm_in_amat1_1(:,:)
                tetra_physics_poly_perpinv(i)%acoef_pre_z = &
                    & tetra_physics_poly4(i)%anorm_in_amat2_0(:,:) + &
                    & perpinv * tetra_physics_poly4(i)%anorm_in_amat2_1(:,:) + &
                    & perpinv2 * tetra_physics_poly4(i)%anorm_in_amat2_2(:,:)
!
                tetra_physics_poly_perpinv(i)%opz_pre_tau =  &
                    & perpinv*tetra_physics_poly4(i)%amat1_1(:,:) + &
                    & tetra_physics_poly4(i)%amat1_0(:,:)
!
                tetra_physics_poly_perpinv(i)%opz_pre_tau2 = &
                    & perpinv2*tetra_physics_poly4(i)%amat2_2(:,:) + &
                    & perpinv*tetra_physics_poly4(i)%amat2_1(:,:) + &
                    & tetra_physics_poly4(i)%amat2_0(:,:)
!
                tetra_physics_poly_perpinv(i)%opb_pre_k0 = &
                    & tetra_physics_poly4(i)%amat1_0_in_b0 + &
                    & perpinv*tetra_physics_poly4(i)%amat1_1_in_b0 + &
                    & perpinv*(tetra_physics_poly4(i)%amat1_0_in_b2 + &
                    & perpinv*tetra_physics_poly4(i)%amat1_1_in_b2)
!
                tetra_physics_poly_perpinv(i)%opb_pre_k1 = &
                    & tetra_physics_poly4(i)%amat1_0_in_b1 + &
                    & perpinv*tetra_physics_poly4(i)%amat1_1_in_b1
!
                tetra_physics_poly_perpinv(i)%opb_pre_k3 = &
                    & tetra_physics_poly4(i)%amat1_0_in_b3 + &
                    & perpinv*tetra_physics_poly4(i)%amat1_1_in_b3
!
            enddo
!
            boole_precomp_poly_perpinv(i) = .true.
!
        end subroutine make_precomp_poly_perpinv
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module tetra_physics_poly_precomp_mod

