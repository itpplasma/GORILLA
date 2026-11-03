module orbit_timestep_gorilla_batch_mod
    use iso_c_binding, only: c_ptr, c_null_ptr, c_associated, c_int, c_int64_t, c_double
    use iso_fortran_env, only: real64, int64
    use ieee_arithmetic, only: ieee_is_finite
    use omp_lib, only: omp_get_wtime
    use gorilla_status_mod
    use orbit_timestep_gorilla_mod, only: orbit_timestep_gorilla, check_coordinate_domain, &
        gorilla_mesh_generation, gorilla_mesh_ready
    use tetra_grid_mod, only: tetra_grid, ntetr
    use tetra_physics_mod, only: tetra_physics, cm_over_e, sign_sqg
    use tetra_grid_settings_mod, only: n_field_periods
    use gorilla_settings_mod, only: coord_system, ipusher, poly_order, i_precomp, i_time_tracing_option, &
        handover_processing_kind, boole_adaptive_time_steps, boole_strong_electric_field, boole_array_optional_quantities
    implicit none
    private
    public :: gorilla_batch_t, gorilla_batch_profile_t, initialize_gorilla_batch, finalize_gorilla_batch
    public :: orbit_timestep_gorilla_batch
    integer, parameter, public :: GORILLA_BATCH_OK=0, GORILLA_BATCH_BAD_ARG=1, GORILLA_BATCH_UNSUPPORTED=2, &
        GORILLA_BATCH_NO_MESH=3, GORILLA_BATCH_NO_CUDA=4, GORILLA_BATCH_LAUNCH=5, GORILLA_BATCH_STALE_MESH=6

    type gorilla_batch_t
        private
        type(c_ptr) :: handle = c_null_ptr
        integer :: backend=0, generation=-1, cells=0, coordinates=0, periods=0, metric_sign=0
        real(real64) :: species=0, setup_seconds=0
    end type
    ! Timings are wall seconds. Setup is reported only on the first batch call.
    ! Cell crossings count handovers, including physical exits, rather than partial cells.
    type gorilla_batch_profile_t
        real(real64) :: setup_seconds=0, upload_seconds=0, kernel_seconds=0, download_seconds=0
        real(real64) :: replay_seconds=0, total_seconds=0
        integer(int64) :: cell_crossings=0, replay_cells=0, numerical_failures=0, physical_losses=0, completed=0
    end type
#ifdef GORILLA_HAVE_CUDA
    interface
        integer(c_int) function cuda_create(n, mesh, topo, coord, nper, sqg, species, light, eps, handle) &
                bind(C,name='gorilla_cuda_create')
            import c_int, c_int64_t, c_double, c_ptr
            integer(c_int64_t), value :: n
            real(c_double), intent(in) :: mesh(*)
            integer(c_int64_t), intent(in) :: topo(*)
            integer(c_int), value :: coord,nper,sqg
            real(c_double), value :: species,light,eps
            type(c_ptr), intent(out) :: handle
        end function
        integer(c_int) function cuda_advance(handle,n,rs,ist,crossings,budget,ms) bind(C,name='gorilla_cuda_advance')
            import c_ptr,c_int,c_int64_t,c_double
            type(c_ptr), value :: handle
            integer(c_int64_t), value :: n,budget
            real(c_double), intent(inout) :: rs(*)
            integer(c_int64_t), intent(inout) :: ist(*),crossings(*)
            real(c_double), intent(out) :: ms(*)
        end function
        integer(c_int) function cuda_destroy(handle) bind(C,name='gorilla_cuda_destroy')
            import c_int,c_ptr
            type(c_ptr), intent(inout) :: handle
        end function
    end interface
#endif
contains
    subroutine initialize_gorilla_batch(batch,backend,ierr)
        use constants, only: clight, eps
        type(gorilla_batch_t), intent(inout) :: batch
        character(*), intent(in) :: backend
        integer, intent(out) :: ierr
#ifdef GORILLA_HAVE_CUDA
        real(c_double), allocatable :: mesh(:,:)
        integer(c_int64_t), allocatable :: topo(:,:)
        integer :: i
#endif
        real(real64) :: start
        start=omp_get_wtime()
        call finalize_gorilla_batch(batch)
        ierr=GORILLA_BATCH_NO_MESH
        if (.not.gorilla_mesh_ready) return
        if (.not.allocated(tetra_physics).or..not.allocated(tetra_grid)) return
        if (size(tetra_physics)/=ntetr.or.ntetr<1) return
        select case(backend)
        case('cpu')
            ! Only the polynomial pusher exposes numerical failure separately from loss.
            ierr=GORILLA_BATCH_UNSUPPORTED
            if (ipusher/=2) return
            batch%backend=1
        case('cuda')
            ierr=GORILLA_BATCH_UNSUPPORTED
            if (.not.cuda_settings_supported()) return
#ifndef GORILLA_HAVE_CUDA
            ierr=GORILLA_BATCH_NO_CUDA
            return
#else
            allocate(mesh(ntetr,59),topo(ntetr,16))
            do i=1,ntetr
                associate(p=>tetra_physics(i),g=>tetra_grid(i))
                    mesh(i,1:3)=p%x1
                    mesh(i,4:15)=reshape(p%anorm,[12])
                    mesh(i,16)=p%dist_ref
                    mesh(i,17)=p%dt_dtau_const
                    mesh(i,18)=p%bmod1
                    mesh(i,19:21)=p%gB
                    mesh(i,22)=p%Phi1
                    mesh(i,23:25)=p%gPhi
                    mesh(i,26:28)=p%curlh
                    mesh(i,29:31)=p%gBxh1
                    mesh(i,32:34)=p%gPhixh1
                    mesh(i,35)=p%gBxcurlA
                    mesh(i,36)=p%gPhixcurlA
                    mesh(i,37:45)=reshape(p%alpmat,[9])
                    mesh(i,46:54)=reshape(p%betmat,[9])
                    mesh(i,55)=p%spalpmat
                    mesh(i,56)=p%spbetmat
                    mesh(i,57:59)=p%curlA
                    topo(i,1:4)=g%neighbour_tetr
                    topo(i,5:8)=g%neighbour_face
                    topo(i,9:12)=g%neighbour_perbou_phi
                    topo(i,13:16)=g%neighbour_perbou_theta
                end associate
            enddo
            ierr=cuda_create(int(ntetr,c_int64_t),mesh,topo,int(coord_system,c_int), &
                int(n_field_periods,c_int),int(sign_sqg,c_int),cm_over_e,clight,eps,batch%handle)
            if (ierr/=0) return
            batch%backend=2
#endif
        case default
            ierr=GORILLA_BATCH_BAD_ARG
            return
        end select
        batch%generation=gorilla_mesh_generation
        batch%cells=ntetr
        batch%species=cm_over_e
        batch%coordinates=coord_system
        batch%periods=n_field_periods
        batch%metric_sign=sign_sqg
        batch%setup_seconds=omp_get_wtime()-start
        ierr=GORILLA_BATCH_OK
    end subroutine

    subroutine finalize_gorilla_batch(batch)
        type(gorilla_batch_t), intent(inout) :: batch
#ifdef GORILLA_HAVE_CUDA
        integer(c_int) :: ierr
        if(c_associated(batch%handle)) ierr=cuda_destroy(batch%handle)
#endif
        batch%handle=c_null_ptr
        batch%backend=0
        batch%generation=-1
        batch%cells=0
        batch%setup_seconds=0
    end subroutine

    logical function cuda_settings_supported() result(ok)
        ok=ipusher==2.and.poly_order==2.and.i_precomp==0.and.i_time_tracing_option==1.and. &
            handover_processing_kind==1.and..not.boole_adaptive_time_steps.and..not.boole_strong_electric_field.and. &
            .not.any(boole_array_optional_quantities).and.(coord_system==1.or.coord_system==2)
    end function

    ! Reinitialize handles after rebuilding the mesh or changing the species/settings.
    subroutine orbit_timestep_gorilla_batch(batch,x,vpar,vperp,t_step,initialized,ind_tetr,iface,status,profile,ierr)
        type(gorilla_batch_t), intent(inout) :: batch
        real(real64), intent(inout) :: x(:,:),vpar(:),vperp(:)
        real(real64), intent(in) :: t_step
        logical, intent(inout) :: initialized(:)
        integer, intent(inout) :: ind_tetr(:),iface(:)
        integer, intent(out) :: status(:)
        type(gorilla_batch_profile_t), intent(out), optional :: profile
        integer, intent(out), optional :: ierr
        integer :: error_code
        type(gorilla_batch_profile_t) :: p
        integer :: n,i
        integer(int64) :: crossings,marker_crossings
        real(real64) :: start
        start=omp_get_wtime()
        p%setup_seconds=batch%setup_seconds
        batch%setup_seconds=0
        status=GORILLA_NUMERICAL_FAILURE
        error_code=GORILLA_BATCH_BAD_ARG
        n=size(vpar)
        if(size(x,1)/=3.or.size(x,2)/=n.or.size(vperp)/=n.or.size(initialized)/=n.or. &
            size(ind_tetr)/=n.or.size(iface)/=n.or.size(status)/=n) goto 900
        if(batch%backend==0) goto 900
        error_code=GORILLA_BATCH_STALE_MESH
        if(batch%generation/=gorilla_mesh_generation.or..not.gorilla_mesh_ready) goto 900
        if(.not.allocated(tetra_physics).or..not.allocated(tetra_grid)) goto 900
        if(batch%cells/=ntetr.or.batch%species/=cm_over_e.or.ipusher/=2.or. &
            batch%coordinates/=coord_system.or.batch%periods/=n_field_periods.or.batch%metric_sign/=sign_sqg) goto 900
        error_code=GORILLA_BATCH_BAD_ARG
        if(.not.ieee_is_finite(t_step)) goto 900
        error_code=GORILLA_BATCH_OK
        if(batch%backend==1) then
            crossings=0
            !$omp parallel do private(i,marker_crossings) reduction(+:crossings) schedule(static)
            do i=1,n
                if(.not.all(ieee_is_finite(x(:,i))).or..not.ieee_is_finite(vpar(i)).or. &
                    .not.ieee_is_finite(vperp(i))) cycle
                call orbit_timestep_gorilla(x(:,i),vpar(i),vperp(i),t_step,initialized(i),ind_tetr(i),iface(i), &
                    status=status(i),crossing_count=marker_crossings)
                crossings=crossings+marker_crossings
            enddo
            !$omp end parallel do
            p%cell_crossings=crossings
        else
#ifdef GORILLA_HAVE_CUDA
            if (.not.cuda_settings_supported()) then
                error_code=GORILLA_BATCH_UNSUPPORTED
                goto 900
            endif
            call advance_cuda(batch,x,vpar,vperp,t_step,initialized,ind_tetr,iface,status,p,error_code)
#endif
        endif
900     p%completed=count(status==GORILLA_COMPLETE,kind=int64)
        p%physical_losses=count(status==GORILLA_PHYSICAL_LOSS,kind=int64)
        p%numerical_failures=count(status==GORILLA_NUMERICAL_FAILURE,kind=int64)
        p%total_seconds=omp_get_wtime()-start
        if(present(profile)) profile=p
        if(present(ierr)) ierr=error_code
    end subroutine
#ifdef GORILLA_HAVE_CUDA
    subroutine advance_cuda(batch,x,vpar,vperp,t_step,initialized,ind_tetr,iface,status,p,error_code)
        use supporting_functions_mod, only: bmod_func,vperp_func
        use find_tetra_mod, only: find_tetra
        use pusher_tetra_poly_mod, only: initialize_const_motion_poly,pusher_tetra_poly
        type(gorilla_batch_t), intent(in) :: batch
        real(real64), intent(inout) :: x(:,:),vpar(:),vperp(:)
        real(real64), intent(in) :: t_step
        logical, intent(inout) :: initialized(:)
        integer, intent(inout) :: ind_tetr(:),iface(:)
        integer, intent(out) :: status(:)
        type(gorilla_batch_profile_t), intent(inout) :: p
        integer, intent(out) :: error_code
        real(c_double), allocatable :: rs(:,:)
        integer(c_int64_t), allocatable :: ist(:,:),crossings(:)
        real(c_double) :: ms(5)
        real(real64) :: z(3),entry_x(3),entry_v,t_pass,bmag,start
        integer :: i,n,cell,face,entry_cell,cell_status,iper,ierr,launch
        logical :: finished
        error_code=GORILLA_BATCH_OK
        n=size(vpar)
        allocate(rs(n,9),ist(n,4),crossings(n))
        rs=0; ist=0; crossings=0
        rs(:,1:3)=transpose(x); rs(:,4)=vpar; rs(:,6)=t_step
        ist(:,1)=ind_tetr; ist(:,2)=iface
        ist(:,4)=GORILLA_NUMERICAL_FAILURE
        do i=1,n
            if(.not.all(ieee_is_finite(x(:,i))).or..not.ieee_is_finite(vpar(i)).or. &
                .not.ieee_is_finite(vperp(i))) cycle
            if(.not.initialized(i)) then
                call check_coordinate_domain(x(:,i))
                call find_tetra(x(:,i),vpar(i),vperp(i),ind_tetr(i),iface(i),int(sign(1.d0,t_step)))
                initialized(i)=ind_tetr(i)>0
            endif
            rs(i,1:3)=x(:,i); rs(i,4)=vpar(i); rs(i,6)=t_step
            ist(i,1)=ind_tetr(i); ist(i,2)=iface(i)
            if(ind_tetr(i)<1) then
                ist(i,4)=GORILLA_PHYSICAL_LOSS
                cycle
            endif
            if(ind_tetr(i)>ntetr.or.iface(i)<0.or.iface(i)>4) cycle
            cell=ind_tetr(i)
            z=x(:,i)-tetra_physics(cell)%x1
            rs(i,7:9)=z; ist(i,3)=cell
            bmag=bmod_func(z,cell)
            if(.not.ieee_is_finite(bmag).or.bmag<=0) cycle
            rs(i,5)=-0.5d0*vperp(i)**2/bmag
            ist(i,4)=GORILLA_ACTIVE
            if(t_step==0) ist(i,4)=GORILLA_COMPLETE
        enddo
        ! Each uncertain cell is replayed once; failed native solves are terminal.
        launch=0
        do while(any(ist(:,4)==GORILLA_ACTIVE))
            launch=launch+1
            if(launch>1000000) then
                error_code=GORILLA_BATCH_LAUNCH
                exit
            endif
            ierr=cuda_advance(batch%handle,int(n,c_int64_t),rs,ist,crossings,256_c_int64_t,ms)
            if(ierr/=0) then
                error_code=ierr
                exit
            endif
            p%upload_seconds=p%upload_seconds+ms(2)*1.d-3
            p%kernel_seconds=p%kernel_seconds+ms(3)*1.d-3
            p%download_seconds=p%download_seconds+ms(4)*1.d-3
            start=omp_get_wtime()
            do i=1,n
                if(ist(i,4)/=GORILLA_REPLAY_CELL) cycle
                p%replay_cells=p%replay_cells+1
                cell=int(ist(i,1)); face=int(ist(i,2)); entry_cell=cell
                if(cell<1.or.cell>ntetr) then
                    ist(i,4)=GORILLA_NUMERICAL_FAILURE
                    cycle
                endif
                entry_x=rs(i,1:3); entry_v=rs(i,4)
                z=entry_x-tetra_physics(cell)%x1
                call initialize_const_motion_poly(rs(i,5),rs(i,5)**2)
                call pusher_tetra_poly(poly_order,cell,face,rs(i,1:3),rs(i,4),z,rs(i,6), &
                    t_pass,finished,iper,status=cell_status)
                if(cell_status==GORILLA_NUMERICAL_FAILURE) then
                    rs(i,1:3)=entry_x; rs(i,4)=entry_v
                    rs(i,7:9)=entry_x-tetra_physics(entry_cell)%x1
                    ist(i,3)=entry_cell
                    ist(i,4)=GORILLA_NUMERICAL_FAILURE
                    cycle
                endif
                rs(i,6)=rs(i,6)-t_pass
                rs(i,7:9)=z
                ist(i,1)=cell; ist(i,2)=face; ist(i,3)=entry_cell; ist(i,4)=cell_status
                if(cell_status/=GORILLA_COMPLETE) crossings(i)=crossings(i)+1
            enddo
            p%replay_seconds=p%replay_seconds+omp_get_wtime()-start
        enddo
        where(ist(:,4)==GORILLA_ACTIVE.or.ist(:,4)==GORILLA_REPLAY_CELL) ist(:,4)=GORILLA_NUMERICAL_FAILURE
        p%cell_crossings=sum(crossings)
        do i=1,n
            x(:,i)=rs(i,1:3); vpar(i)=rs(i,4)
            ind_tetr(i)=int(ist(i,1)); iface(i)=int(ist(i,2)); status(i)=int(ist(i,4))
            cell=int(ist(i,3))
            if(cell>0.and.cell<=ntetr) vperp(i)=vperp_func(rs(i,7:9),rs(i,5),cell)
        enddo
    end subroutine
#endif
end module
