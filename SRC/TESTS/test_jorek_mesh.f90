program test_jorek_mesh
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use jorek_mesh_mod, only: build_jorek_mesh_arrays
    use jorek_restart, only: jorek_restart_t

    implicit none

    type(jorek_restart_t) :: data
    real(dp), allocatable :: points(:, :)
    integer, allocatable :: verts(:, :), neighbours(:, :), faces(:, :), perbou(:, :)
    integer, allocatable :: owner(:)
    real(dp), allocatable :: owner_st(:, :)
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: r_node(4) = [10.0_dp, 11.0_dp, 11.0_dp, 10.0_dp]
    real(dp), parameter :: z_node(4) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    integer :: face, ierr, neighbour, node, tetra

    data%jorek_model = 303
    data%n_order = 3
    data%n_degrees = 4
    data%n_coord_tor = 1
    data%n_dim = 2
    data%n_tor = 1
    data%n_period = 1
    data%n_var = 7
    data%n_nodes = 4
    data%n_elements = 1
    data%n_vertex_max = 4
    allocate(data%x(4, 1, 4, 2), data%values(4, 1, 4, 7))
    allocate(data%vertex(1, 4), data%size(1, 4, 4))
    data%x = 0.0_dp
    data%values = 0.0_dp
    data%vertex(1, :) = [1, 2, 3, 4]
    data%size = 1.0_dp
    data%size(:, 2:3, 2) = -1.0_dp
    data%size(:, 3:4, 3) = -1.0_dp
    data%size(:, 2, 4) = -1.0_dp
    data%size(:, 4, 4) = -1.0_dp
    do node = 1, 4
        data%x(node, 1, 1, 1) = r_node(node)
        data%x(node, 1, 1, 2) = z_node(node)
        data%x(node, 1, 2, 1) = 1.0_dp/3.0_dp
        data%x(node, 1, 3, 2) = 1.0_dp/3.0_dp
        data%values(node, 1, 1, 1) = r_node(node)*z_node(node)
        data%values(node, 1, 2, 1) = z_node(node)/3.0_dp
        data%values(node, 1, 3, 1) = r_node(node)/3.0_dp
        data%values(node, 1, 4, 1) = 1.0_dp/9.0_dp
    end do

    call build_jorek_mesh_arrays(data, 3, 100.0_dp, points, verts, &
        neighbours, faces, perbou, ierr, owner, owner_st)
    if (ierr /= 0) then
        write (*, '(A, I0)') 'synthetic JOREK mesh construction error: ', ierr
        error stop 'synthetic JOREK mesh construction failed'
    end if
    if (any(shape(points) /= [3, 12]) .or. any(shape(verts) /= [4, 18])) &
        error stop 'synthetic JOREK mesh has wrong dimensions'
    if (maxval(abs(points(1, 1:4) - 100.0_dp*r_node)) > 1.0e-12_dp &
            .or. maxval(abs(points(3, 1:4) - 100.0_dp*z_node)) > 1.0e-12_dp) &
        error stop 'JOREK mesh length conversion failed'
    if (maxval(abs(points(2, 1:4))) > 1.0e-14_dp &
            .or. maxval(abs(points(2, 5:8) - 2.0_dp*pi/3.0_dp)) > 1.0e-14_dp &
            .or. maxval(abs(points(2, 9:12) - 4.0_dp*pi/3.0_dp)) > 1.0e-14_dp) &
        error stop 'JOREK mesh toroidal extrusion failed'
    if (any(verts < 1) .or. any(verts > 12) .or. count(perbou /= 0) == 0) &
        error stop 'JOREK tetrahedral connectivity is invalid'
    if (any(owner /= 1)) error stop 'JOREK mesh vertex owner is wrong'
    if (maxval(abs(owner_st(:, 1:4) - reshape([0.0_dp, 0.0_dp, &
            1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp], [2, 4]))) &
            > 1.0e-14_dp) error stop 'JOREK mesh owner coordinates are wrong'
    if (maxval(abs(owner_st(:, 5:8) - owner_st(:, 1:4))) > 1.0e-14_dp &
            .or. maxval(abs(owner_st(:, 9:12) - owner_st(:, 1:4))) &
            > 1.0e-14_dp) error stop 'JOREK mesh owner coordinates changed by extrusion'

    do tetra = 1, size(neighbours, 2)
        do face = 1, 4
            neighbour = neighbours(face, tetra)
            if (neighbour == -1) cycle
            if (neighbours(faces(face, tetra), neighbour) /= tetra) &
                error stop 'JOREK mesh neighbour relation is not reciprocal'
            if (perbou(faces(face, tetra), neighbour) /= -perbou(face, tetra)) &
                error stop 'JOREK mesh periodic orientation is inconsistent'
        end do
    end do

    call build_jorek_mesh_arrays(data, 3, 100.0_dp, points, verts, &
        neighbours, faces, perbou, ierr, owner, owner_st, &
        poloidal_subdivisions=2)
    if (ierr /= 0) error stop 'refined JOREK mesh construction failed'
    if (any(shape(points) /= [3, 27]) .or. any(shape(verts) /= [4, 72])) &
        error stop 'refined JOREK mesh has wrong dimensions'
    if (any(owner /= 1)) error stop 'refined JOREK mesh owner is wrong'
    node = minloc(sum((owner_st(:, 1:9) - 0.5_dp)**2, dim=1), dim=1)
    if (maxval(abs(points([1, 3], node) - [1050.0_dp, 50.0_dp])) &
            > 1.0e-12_dp) error stop 'refined JOREK center is wrong'

    call build_jorek_mesh_arrays(data, 3, 100.0_dp, points, verts, &
        neighbours, faces, perbou, ierr, owner, owner_st, &
        poloidal_subdivisions=4)
    if (ierr /= 9) error stop 'unsafe JOREK volume refinement was accepted'

    call build_jorek_mesh_arrays(data, 3, 0.0_dp, points, verts, &
        neighbours, faces, perbou, ierr, owner, owner_st)
    if (ierr /= 5) error stop 'invalid JOREK mesh scale was accepted'

    print '(A)', 'PASS: JOREK-derived tetrahedral mesh'
end program test_jorek_mesh
