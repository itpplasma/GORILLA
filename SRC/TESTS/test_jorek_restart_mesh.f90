program test_jorek_restart_mesh
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_mesh_mod, only: build_jorek_mesh_arrays
    use jorek_restart, only: jorek_restart_t, load_jorek_restart

    implicit none

    type(jorek_restart_t) :: data
    real(dp), allocatable :: points(:, :)
    integer, allocatable :: verts(:, :), neighbours(:, :), faces(:, :), perbou(:, :)
    integer, allocatable :: owner(:)
    real(dp), allocatable :: owner_st(:, :)
    character(len=1024) :: filename
    integer :: ierr, n_plane, n_tetras

    if (command_argument_count() /= 1) error stop 'restart path argument is required'
    call get_command_argument(1, filename)
    call load_jorek_restart(trim(filename), data, ierr)
    if (ierr /= 0) error stop 'golden JOREK restart load failed'

    call build_jorek_mesh_arrays(data, 3, 100.0_dp, points, verts, &
        neighbours, faces, perbou, ierr)
    if (ierr /= 0) then
        write (*, '(A, I0)') 'golden JOREK mesh construction error: ', ierr
        error stop 'golden JOREK mesh construction failed'
    end if
    n_plane = 6167
    n_tetras = 109287
    if (any(shape(points) /= [3, 3*n_plane]) &
            .or. any(shape(verts) /= [4, n_tetras])) &
        error stop 'golden JOREK mesh has wrong dimensions'
    if (.not. all(ieee_is_finite(points))) &
        error stop 'golden JOREK mesh contains non-finite coordinates'
    if (any(verts < 1) .or. any(verts > size(points, 2))) &
        error stop 'golden JOREK mesh contains an invalid vertex index'
    if (maxval(abs(points([1, 3], n_plane + 1:2*n_plane) &
            - points([1, 3], 1:n_plane))) > 1.0e-12_dp) &
        error stop 'golden JOREK mesh extrusion changed the poloidal plane'
    if (count(perbou /= 0) == 0 .or. any(abs(perbou) > 1)) &
        error stop 'golden JOREK mesh periodic flags are invalid'

    call build_jorek_mesh_arrays(data, 3, 100.0_dp, points, verts, &
        neighbours, faces, perbou, ierr, owner, owner_st, &
        poloidal_subdivisions=2)
    if (ierr /= 0) then
        write (*, '(A, I0)') 'refined golden JOREK mesh error: ', ierr
        error stop 'refined golden JOREK mesh construction failed'
    end if
    if (any(shape(points) /= [3, 73788]) &
            .or. any(shape(verts) /= [4, 439218])) &
        error stop 'refined golden JOREK mesh has wrong dimensions'
    if (count(owner == 0) /= 3 .or. any(owner < 0) &
            .or. any(owner > data%n_elements)) &
        error stop 'refined golden JOREK mesh owners are invalid'
    if (any(owner_st < 0.0_dp) .or. any(owner_st > 1.0_dp)) &
        error stop 'refined golden JOREK owner coordinates are invalid'

    print '(A, I0, A, I0)', 'PASS: golden JOREK mesh with ', &
        size(points, 2), ' vertices and ', size(verts, 2)
end program test_jorek_restart_mesh
