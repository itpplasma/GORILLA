program test_orbit_batch_driver
    use iso_fortran_env, only: real64
    use gorilla_status_mod
    use orbit_timestep_gorilla_mod, only: initialize_gorilla,orbit_timestep_gorilla
    use orbit_timestep_gorilla_batch_mod
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use tetra_grid_mod, only: tetra_grid,verts_rphiz,ntetr
    use tetra_physics_mod, only: tetra_physics
    use pusher_tetra_poly_mod, only: initialize_const_motion_poly,pusher_tetra_poly
    implicit none
    type(gorilla_batch_t) :: batch
    type(gorilla_batch_profile_t) :: profile
    real(real64) :: x(3,3),xs(3,3),vp(3),vps(3),vt(3),vts(3),dt,z(3),passed,remaining
    integer :: cells(3),faces(3),cs(3),fs(3),status(3),ss(3),i,j,u,v,ios,ierr,iper
    logical :: initialized(3),is(3),finished
    character(1024) :: root,line

    call get_command_argument(1,root)
    ! Reuse the public example's complete namelist, with analytic geometry below.
    open(newunit=u,file=trim(root)//'/EXAMPLES/example_1/gorilla.inp',status='old')
    open(newunit=v,file='gorilla.inp',status='replace')
    do
        read(u,'(A)',iostat=ios) line
        if(ios/=0) exit
        if(index(line,'coord_system =')>0) line='coord_system = 1 ,'
        if(index(line,'poly_order =')>0) line='poly_order = 2 ,'
        write(v,'(A)') trim(line)
    enddo
    close(u); close(v)
    open(newunit=u,file='tetra_grid.inp',status='replace')
    write(u,'(A)') '&TETRA_GRID_NML', 'grid_kind=5,n1=4,n2=12,n3=4,', &
        'boole_n_field_periods=.true.,boole_write_mesh_obj=.false.,', &
        'R0_analytic_circ=1000.,a_analytic_circ=100.,B0_analytic_circ=10000.,', &
        'q0_analytic_circ=2.,q1_analytic_circ=0.,/'
    close(u)
    open(newunit=u,file='field_divB0.inp',status='replace')
    write(u,'(A)') '0','1','0.0','1','0.0','0', &
        '"unused"','"unused"','"unused"','"."','1','1','0'
    close(u)
    call load_gorilla_inp()
    call load_tetra_grid_inp()
    call initialize_gorilla()
    call initialize_gorilla_batch(batch,'cpu',ierr)
    if(ierr/=GORILLA_BATCH_OK) error stop 'CPU initialization'

    do j=-1,1
        call seed()
        xs=x; vps=vp; vts=vt; cs=cells; fs=faces; is=initialized
        dt=j*1.d-6
        do i=1,3
            call orbit_timestep_gorilla(xs(:,i),vps(i),vts(i),dt,is(i),cs(i),fs(i),status=ss(i))
        enddo
        call orbit_timestep_gorilla_batch(batch,x,vp,vt,dt,initialized,cells,faces,status,profile,ierr)
        if(ierr/=0.or.any(status/=ss)) error stop 'scalar/batch statuses'
        if(any(status/=GORILLA_COMPLETE)) error stop 'valid marker completion'
        if(maxval(abs(x-xs))>1.d-11.or.maxval(abs(vp-vps))>1.d-7.or. &
            maxval(abs(vt-vts))>1.d-7) error stop 'scalar/batch trajectories'
        if(any(cells/=cs).or.any(faces/=fs)) error stop 'scalar/batch topology'
        if(profile%completed/=3.or.profile%numerical_failures/=0) error stop 'profile outcomes'
        if(j>-1.and.profile%setup_seconds/=0) error stop 'setup must be charged once'
        if(j==0.and.profile%cell_crossings/=0) error stop 'zero step has no handoffs'
    enddo

    call seed()
    cells(1)=-1
    cells(2)=ntetr+1
    call orbit_timestep_gorilla_batch(batch,x,vp,vt,0.d0,initialized,cells,faces,status,ierr=ierr)
    if(status(1)/=GORILLA_PHYSICAL_LOSS.or.status(2)/=GORILLA_NUMERICAL_FAILURE.or. &
        status(3)/=GORILLA_COMPLETE) error stop 'initial cell classification'

    ! A stationary trajectory has no exit root: it must not consume time or become a physical loss.
    call seed()
    vp=0; vt=0
    xs=x; cs=cells; fs=faces
    call initialize_const_motion_poly(0.d0,0.d0)
    call pusher_tetra_poly(2,cells(1),faces(1),x(:,1),vp(1),z,1.d-6,passed,finished,iper,status=ss(1))
    if(ss(1)/=GORILLA_NUMERICAL_FAILURE.or.passed/=0.or.finished.or.iper/=0) &
        error stop 'safe failed cell outputs'
    if(any(z/=xs(:,1)-tetra_physics(cs(1))%x1)) error stop 'safe failed cell coordinates'
    x=xs; cells=cs; faces=fs
    call orbit_timestep_gorilla(x(:,1),vp(1),vt(1),1.d-6,initialized(1),cells(1),faces(1),remaining,ss(1))
    if(ss(1)/=GORILLA_NUMERICAL_FAILURE.or.remaining/=1.d-6.or.cells(1)/=cs(1).or. &
        faces(1)/=fs(1).or.any(x(:,1)/=xs(:,1))) error stop 'failed scalar state preservation'
    call orbit_timestep_gorilla(x(:,1),vp(1),vt(1),1.d-6,initialized(1),cells(1),faces(1))
    if(cells(1)/=-1) error stop 'legacy numerical failure indicator'

    call seed()
    call seed_boundary_exit()
    if(ss(1)/=GORILLA_COMPLETE.or.cells(1)<1) error stop 'boundary seed starts inside the mesh'
    xs=x; vps=vp; vts=vt; cs=cells; fs=faces; is=initialized
    call orbit_timestep_gorilla(xs(:,1),vps(1),vts(1),1.d-5,is(1),cs(1),fs(1),remaining,ss(1))
    call orbit_timestep_gorilla_batch(batch,x,vp,vt,1.d-5,initialized,cells,faces,status,profile,ierr)
    if(status(1)/=GORILLA_PHYSICAL_LOSS.or.cells(1)>=0.or.profile%physical_losses<1.or. &
        profile%cell_crossings<1) error stop 'integrated physical boundary exit'
    if(ss(1)/=status(1).or.remaining<=0.or.remaining>=1.d-5.or. &
        maxval(abs(x(:,1)-xs(:,1)))>1.d-11.or.abs(vp(1)-vps(1))>1.d-7.or. &
        abs(vt(1)-vts(1))>1.d-7) error stop 'physical exit agrees with scalar'

    call seed()
    call initialize_gorilla(i_option=2)
    xs=x
    call orbit_timestep_gorilla_batch(batch,x,vp,vt,0.d0,initialized,cells,faces,status,ierr=ierr)
    if(ierr/=GORILLA_BATCH_STALE_MESH.or.any(status/=GORILLA_NUMERICAL_FAILURE).or.any(x/=xs)) &
        error stop 'stale mesh rejection'
    call initialize_gorilla_batch(batch,'cpu',ierr)
    if(ierr/=0) error stop 'reinitialize rebuilt mesh'

    call initialize_gorilla(i_option=1)
    call initialize_gorilla_batch(batch,'cpu',ierr)
    if(ierr/=GORILLA_BATCH_NO_MESH) error stop 'grid without matching physics'
    call initialize_gorilla(i_option=2)
    call initialize_gorilla_batch(batch,'cuda',ierr)
    if(ierr/=GORILLA_BATCH_NO_CUDA.and.ierr/=GORILLA_BATCH_OK) error stop 'explicit CUDA admission'
    call finalize_gorilla_batch(batch)
    call finalize_gorilla_batch(batch)
    call seed()
    xs=x
    call orbit_timestep_gorilla_batch(batch,x,vp,vt,0.d0,initialized,cells,faces,status,ierr=ierr)
    if(ierr==0.or.any(status/=GORILLA_NUMERICAL_FAILURE).or.any(x/=xs)) error stop 'finalized handle'
    print *, 'orbit batch oracle: PASS'
contains
    subroutine seed_boundary_exit()
        real(real64) :: vertices(3,4),centroid(3),face_centre(3),candidate(3),trial(3)
        real(real64) :: trial_vpar,trial_vperp,trial_remaining
        integer :: cell,face,direction,trial_cell,trial_face,trial_status
        logical :: trial_initialized

        ! Boundary coordinates must follow the mesh, including circular variants.
        do cell=1,ntetr
            vertices=verts_rphiz(:,tetra_grid(cell)%ind_knot)
            centroid=sum(vertices,dim=2)/4.d0
            do face=1,4
                if(tetra_grid(cell)%neighbour_tetr(face)>=0) cycle
                face_centre=(sum(vertices,dim=2)-vertices(:,face))/3.d0
                candidate=face_centre+1.d-4*(centroid-face_centre)
                do direction=-1,1,2
                    trial=candidate; trial_vpar=direction*1.d8; trial_vperp=1.d7
                    trial_cell=0; trial_face=0; trial_initialized=.false.
                    call orbit_timestep_gorilla(trial,trial_vpar,trial_vperp,0.d0, &
                        trial_initialized,trial_cell,trial_face,status=trial_status)
                    if(trial_status/=GORILLA_COMPLETE.or.trial_cell<1) cycle
                    call orbit_timestep_gorilla(trial,trial_vpar,trial_vperp,1.d-5, &
                        trial_initialized,trial_cell,trial_face,trial_remaining,trial_status)
                    if(trial_status/=GORILLA_PHYSICAL_LOSS) cycle
                    if(trial_remaining<=0.or.trial_remaining>=1.d-5) cycle
                    x(:,1)=candidate; vp(1)=direction*1.d8; vt(1)=1.d7
                    initialized(1)=.false.
                    call orbit_timestep_gorilla(x(:,1),vp(1),vt(1),0.d0,initialized(1), &
                        cells(1),faces(1),status=ss(1))
                    return
                enddo
            enddo
        enddo
        error stop 'mesh has no integrated boundary-exit seed'
    end subroutine

    subroutine seed()
        integer :: k,cell
        do k=1,3
            cell=ntetr/2+6*k
            cells(k)=cell
            x(:,k)=sum(verts_rphiz(:,tetra_grid(cell)%ind_knot),dim=2)/4.d0
        enddo
        vp=[1.d7,-1.d7,5.d6]; vt=1.d7
        faces=0; initialized=.true.
    end subroutine
end program
