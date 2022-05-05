module circular_mesh_SOLEDGE3X_EIRENE
implicit none
private
integer, dimension(:,:), allocatable :: triangles_SOLEDGE3X_EIRENE
integer :: n_triangles
double precision, parameter :: pi = 3.14159265358979d0
public :: calc_mesh_SOLEDGE3X_EIRENE, create_points_SOLEDGE3X_EIRENE

contains

subroutine create_points_SOLEDGE3X_EIRENE(n_slices, points_rphiz, verts_per_slice)
!
    integer, intent(in) :: n_slices
    integer, intent(out) :: verts_per_slice
    double precision, dimension(:, :), allocatable, intent(out) :: points_rphiz !(r, phi, z)
!
    integer :: file_id_knots
    integer, dimension(2) :: shape_knots
    character(50) :: filename_knots
    double precision, dimension(:,:),allocatable :: knots_SOLEDGE3X_EIRENE
!
    integer :: n_verts,phi_position, i
!           
    !Load SOLEDGE3X_EIRENE mesh data
!
    !Define file ids and filenames for SOLEDGE3X_EIRENE mesh data
    file_id_knots = 501
    filename_knots = './MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/knots.dat'
!
    !Load knots
    open(unit=file_id_knots, file=filename_knots, status='unknown')
    read(file_id_knots,*) shape_knots
    allocate(knots_SOLEDGE3X_EIRENE(shape_knots(1),shape_knots(2)))
    do i = 1,shape_knots(1)
        read(file_id_knots,*) knots_SOLEDGE3X_EIRENE(i,:)
    enddo
    close(file_id_knots)
!
    !Number of vertices per 2d-plane
    verts_per_slice = shape_knots(1)
!
    !Total number of vertices  
    n_verts = verts_per_slice*n_slices
!
    if (allocated(points_rphiz)) deallocate(points_rphiz)
!
    allocate(points_rphiz(3, n_verts))
!    
    !Fill first slice of point matrix with 2d-SOLEDGE3X_EIRENE-knots
    points_rphiz(1, :verts_per_slice) = knots_SOLEDGE3X_EIRENE(:,1)
    points_rphiz(2, :verts_per_slice) = 0.d0
    points_rphiz(3, :verts_per_slice) = knots_SOLEDGE3X_EIRENE(:,2)
!      
    phi_position = 2 !In coordinate system (R,Phi,Z) --> Phi is at second position
    call extrude_points(verts_per_slice, n_slices,phi_position, points_rphiz)
!
    deallocate(knots_SOLEDGE3X_EIRENE)
!
end subroutine create_points_SOLEDGE3X_EIRENE
!
subroutine calc_mesh_SOLEDGE3X_EIRENE(n_slices, points_rphiz, verts_per_slice, n_tetras, & 
                     verts, neighbours, neighbour_faces, perbou_phi)
!
    use circular_mesh, only: wrap_idx_inplace
!
    integer, intent(in) :: n_slices,verts_per_slice
    double precision, dimension(:,:), intent(in) :: points_rphiz 
    integer, intent(out) :: n_tetras
    integer :: i,j,f
    integer, dimension(:, :), allocatable, intent(out) :: verts, neighbours, neighbour_faces, perbou_phi
!
    double precision, DIMENSION(verts_per_slice) :: A_phi_vec
!
    integer :: file_id_triangles
    integer, dimension(2) :: shape_triangles
    character(50) :: filename_triangles
    integer, dimension(:,:), allocatable :: triangle_type
!
    integer, dimension(4, 6) :: tetra_conf, mask_theta, mask_phi, mask_r
    integer, dimension(4, 3) :: slice_offset, ring_offset, cur_triangle_offset, no_offset
    integer :: n_verts, mask_idx,prism_orientation, slice,&
                & cur_triangle, cur_stand_alone_vertex, tetra_idx, prism_i, prism_j, tetras_per_slice, &
                vertex_I, vertex_II, vertex_III
    integer, dimension(:), allocatable :: count_connected
    logical :: match
!         
    !Compute A_phi (psif) for every vertex in 2d-plane
    do i = 1,verts_per_slice
        call vector_potential_rz(points_rphiz(1,i),points_rphiz(3,i),A_phi_vec(i))
    enddo
!
    file_id_triangles = 502
    filename_triangles = './MHD_EQUILIBRIA/MESH_SOLEDGE3X_EIRENE/triangles.dat'
!
    !Load triangles
    open(unit=file_id_triangles, file=filename_triangles, status='unknown')
    read(file_id_triangles,*) shape_triangles
    allocate(triangles_SOLEDGE3X_EIRENE(shape_triangles(1),shape_triangles(2)))
    do i = 1,shape_triangles(1)
        read(file_id_triangles,*) triangles_SOLEDGE3X_EIRENE(i,:)
    enddo
    close(file_id_triangles)
!
    !Number of triagles in 2d-plane
    n_triangles = shape_triangles(1)
!  
    !Compute triangle types (top/bottom with respect to psi) and "stand-alone" vertex
    allocate(triangle_type(n_triangles,2))
    call calc_triangle_type(A_phi_vec,triangle_type)
! 
    if ( n_slices < 3) then
        print *, "Invalid parameters for function calc_mesh_SOLEDGE3X_EIRENE"
        stop
    end if
!    
    n_verts = verts_per_slice * n_slices
!    
    n_tetras = n_triangles*3*n_slices
!
    if (.not. size(verts, dim=2) == n_tetras .and. size(neighbours, dim=2) == n_tetras &
            .and. size(neighbour_faces, dim=2) == n_tetras) then
        print *, "Invalid size of output arrays for function calc_mesh_SOLEDGE3X_EIRENE"
        stop
    end if
!
    tetras_per_slice = n_tetras / n_slices

    if (allocated(verts)) deallocate(verts)
    if (allocated(neighbours)) deallocate(neighbours)
    if (allocated(neighbour_faces)) deallocate(neighbour_faces)
    if (allocated(perbou_phi)) deallocate(perbou_phi)
    allocate(verts(4, n_tetras), neighbours(4, n_tetras), neighbour_faces(4, n_tetras), &
             perbou_phi(4, n_tetras))
!
! 
!     7------6
!    /|     /|
!   3------2 |
!   | 5----|-4           r     
!   |/     |/            ^ ,phi
!   1------0 *           |/    
!               theta <--*
!
    ! top facing prism
    ! top edge from 2 - 7
    tetra_conf(:, :3) = reshape((/ &
        0, 2, 3, 7,  & !1 
        0, 2, 6, 7,  & !1
        0, 4, 6, 7  & !1
        /), (/4, 3/))
!    
    ! calculate the bottom facing prism tetraeder configurations from the top facing one one
    ! bottom facing prism
    tetra_conf(:, 4:) = tetra_conf(:, :3)
    where (iand(tetra_conf(:, 4:) , 2) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 1)
    end where
    where (iand(tetra_conf(:, 4:) , 1) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 2)
    end where
!
    do i = 1, 6
        print *, tetra_conf(:, i)
        if (mod(i, 3) == 0) print *, ''
    end do
!
    ! extract masks from tetra_conf
    mask_theta = iand(tetra_conf, 1)
    mask_r = ishft(iand(tetra_conf, 2), -1)
    mask_phi = ishft(iand(tetra_conf, 4), -2)
!
    ! initialize neighbour_faces
    neighbour_faces = - 1
!
    ! first slice
    tetra_idx = 1

    do cur_triangle = 1, n_triangles
!
        prism_orientation = triangle_type(cur_triangle,1)
        cur_stand_alone_vertex = triangle_type(cur_triangle,2)
!
        vertex_I = triangles_SOLEDGE3X_EIRENE(cur_triangle,cur_stand_alone_vertex)
        vertex_II = triangles_SOLEDGE3X_EIRENE(cur_triangle,mod(cur_stand_alone_vertex ,3) + 1)
        vertex_III = triangles_SOLEDGE3X_EIRENE(cur_triangle,mod(cur_stand_alone_vertex + 1,3) + 1)
!      
        ! --- get the correct offset masks for up- or down facing prism ---
        mask_idx = 1 + prism_orientation * 3
!
        slice_offset = mask_phi(:, mask_idx:mask_idx + 2) * verts_per_slice
        ring_offset = abs(mask_r(:, mask_idx:mask_idx + 2)-prism_orientation) &
                    & *abs(mask_theta(:, mask_idx:mask_idx + 2) - 1 + prism_orientation) &
                    & *vertex_II
        cur_triangle_offset = abs(mask_r(:, mask_idx:mask_idx + 2) - prism_orientation) &
                    & *abs(mask_theta(:, mask_idx:mask_idx + 2)-prism_orientation) &
                    & *vertex_III
        no_offset = abs(mask_r(:, mask_idx:mask_idx + 2) - 1 + prism_orientation) &
                    & *abs(mask_theta(:, mask_idx:mask_idx + 2) - 1) &
                    & *vertex_I
!
! --- with offset masks calculate verts ---
        tetra_idx = (cur_triangle - 1) * 3 + 1
        verts(:, tetra_idx:tetra_idx + 2) = no_offset + slice_offset + ring_offset + cur_triangle_offset
!
    end do !cur_triangles
!
    if (allocated(count_connected)) deallocate(count_connected)
    allocate(count_connected(n_triangles))
    count_connected = 0
!
    !$OMP PARALLEL DEFAULT(NONE) &
    !$OMP& SHARED(n_triangles,tetras_per_slice,neighbours,neighbour_faces,verts,count_connected) &
    !$OMP& PRIVATE(prism_i,tetra_idx,prism_j,match)
    !$OMP DO
    do prism_i = 1, n_triangles
!
        ! connect with prisms in neighbouring slices
        tetra_idx = (prism_i - 1) * 3 + 1
        neighbours(4, tetra_idx) = tetra_idx + 2 - tetras_per_slice
        neighbour_faces(4, tetra_idx) = 1
        neighbours(1, tetra_idx + 2) = tetra_idx + tetras_per_slice
        neighbour_faces(1, tetra_idx + 2) = 4
!
        do prism_j = prism_i, n_triangles
!
            ! --- connect neighbouring tetras  ---
            ! connect neighbouring tetras in this prism
            call connect_prisms_SOLEDGE3X_EIRENE(prism_i, prism_j, verts, neighbours, neighbour_faces,match)
!
            !$omp critical
            if (match .and. (prism_i /= prism_j)) then
                count_connected(prism_i) = count_connected(prism_i) + 1
                count_connected(prism_j) = count_connected(prism_j) + 1
            end if
            !$omp end critical
!
            !If already three prisms were found that connect to the current prism_i
            !-> stop further search
            if (count_connected(prism_i) > 2) exit
!
        end do !prism_j
!
    end do !prism_i
    !$OMP END DO
    !$OMP END PARALLEL
!
    ! checks for missmatched prims and tries to solve it by changing classification LOCALLY
    call repair(count_connected,triangle_type,mask_r,mask_phi,mask_theta,verts_per_slice,verts,neighbours,neighbour_faces)
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
if(.true.) then
!if (allocated(count_connected)) deallocate(count_connected)
!allocate(count_connected(n_triangles))
count_connected = 0
!
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP& SHARED(n_triangles,tetras_per_slice,neighbours,neighbour_faces,verts,count_connected) &
!$OMP& PRIVATE(prism_i,tetra_idx,prism_j,match)
!$OMP DO
do prism_i = 1, n_triangles
!
do prism_j = prism_i, n_triangles
!
! --- connect neighbouring tetras  ---
! connect neighbouring tetras in this prism
call connect_prisms_SOLEDGE3X_EIRENE(prism_i, prism_j, verts, neighbours, neighbour_faces,match)
!
!$omp critical
if (match .and. (prism_i /= prism_j)) then
count_connected(prism_i) = count_connected(prism_i) + 1
count_connected(prism_j) = count_connected(prism_j) + 1
end if
!$omp end critical
!
!If already three prisms were found that connect to the current prism_i
!-> stop further search
if (count_connected(prism_i) > 2) exit
!
end do !prism_j
!
end do !prism_i
!$OMP END DO
!$OMP END PARALLEL
!
end if
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
    ! for all the other slices we can calculate the verts by incrementing the
    ! vert indices of the first(idx = 0) slice by slice * verts_per_slice
    ! and the neighbours by incrementing the neighbours of the first slice by slice * tetras_per_slice
    ! neighbour_faces and perbou_theta can just be copied

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice
!
        verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = verts(:, 1:tetras_per_slice) + slice * verts_per_slice
!
        ! wrap indices around in the last slice
        if (slice == n_slices - 1) then
            call wrap_idx_inplace(verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1), n_verts)
        end if
    end do
!

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        neighbours(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = &
                neighbours(:, 1:tetras_per_slice) + slice * tetras_per_slice

        ! wrap indices around in the last slice
        if (slice == n_slices - 1) then
            call wrap_idx_inplace(neighbours(:, tetra_idx:tetra_idx + tetras_per_slice - 1), n_tetras)
        end if
    end do

    ! we dont wrap the neighbours index of the first slice earlier otherwise we would need to wrap the 
    ! index of every incrementally calculated slice
    call wrap_idx_inplace(neighbours(:, 1:tetras_per_slice), n_tetras)

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        neighbour_faces(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = neighbour_faces(:, 1:tetras_per_slice)
    end do
    
    ! remove outer faces after modulo operation
    where (neighbour_faces == -1)
        neighbours = -1
    end where

    ! periodic boundary
    perbou_phi = 0
    perbou_phi(4, :tetras_per_slice:3) = -1
    perbou_phi(1, n_tetras - tetras_per_slice + 3::3) = 1

    do i = 1, n_tetras
        do f = 1, 4
            if (.not. (neighbours(f, i) == -1 .and. neighbour_faces(f, i) == -1)) then
                if (.not. neighbours(neighbour_faces(f, i), neighbours(f, i)) == i) then
                    print *, "!!! CONSISTENCY CHECK FOR NEIGHBOURS FAILED, MESH IS BROKEN !!!"
                    stop
                end if
                if (.not. perbou_phi(neighbour_faces(f, i), neighbours(f, i)) == - perbou_phi(f, i)) then
                    print *, "!!! CONSISTENCY CHECK FOR PERBOU FAILED, MESH IS BROKEN !!!"
                    stop
                end if
            end if
        end do
    end do
    print *, "mesh consistency check ok"
!
open(123, file='./MESH_CHECK/border_triangles.dat')
do cur_triangle = 1, n_triangles
    if ((count_connected(cur_triangle) < 3)) then
        write(123,*) triangles_SOLEDGE3X_EIRENE(cur_triangle,:)
    end if
end do
close(123)
!
open(123, file='./MESH_CHECK/border_triangles_types.dat')
do cur_triangle = 1, n_triangles
    if ((count_connected(cur_triangle) < 3)) then
        write(123,*) triangle_type(cur_triangle,1)
    end if
end do
close(123)
!
open(123, file='./MESH_CHECK/border_triangles_error.dat')
do cur_triangle = 1, n_triangles
    if ((count_connected(cur_triangle) < 3)) then
        if(is_not_border(cur_triangle)) then
            write(123,*) triangles_SOLEDGE3X_EIRENE(cur_triangle,:)
        end if
    end if
end do
close(123)
!
open(123, file='./MESH_CHECK/triangles_types.dat')
do cur_triangle = 1, n_triangles
    write(123,*) triangle_type(cur_triangle,1)
end do
close(123)
!
    deallocate(triangles_SOLEDGE3X_EIRENE, triangle_type,count_connected)
!
stop
!
end subroutine calc_mesh_SOLEDGE3X_EIRENE
!
subroutine extrude_points(verts_per_slice, n_slices,phi_position,  points)
!
    use tetra_grid_settings_mod, only: n_field_periods
!
    integer, intent(in) :: verts_per_slice, n_slices,phi_position
    double precision, dimension(:, :), intent(inout) :: points !(r, phi, z)
!
    integer :: slice, vert_idx
    double precision :: phi
    ! copy and rotate slice
    do slice = 2, n_slices
        vert_idx = (slice- 1) * verts_per_slice + 1
        phi = (2.d0 * pi / n_field_periods * (slice - 1)) / n_slices
        points(:, vert_idx:vert_idx + verts_per_slice - 1) = points(:, 1:verts_per_slice)
        points(phi_position, vert_idx:vert_idx + verts_per_slice - 1) = phi
    end do
end subroutine
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine vector_potential_rz(r,z,A_phi)
!
    use field_eq_mod, only : rtf,btf,psif
!
    implicit none
!
    double precision,intent(in) :: r,z
    double precision, intent(out) :: A_phi
    double precision :: phi
    double precision :: B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ    &
                        ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ
!
    phi = 0.d0
!
    call field(r,phi,z,B_r,B_p,B_z,dBrdR,dBrdp,dBrdZ  &
                ,dBpdR,dBpdp,dBpdZ,dBzdR,dBzdp,dBzdZ)
!
    A_phi=psif
!
end subroutine vector_potential_rz
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine calc_triangle_type(A_phi_vec,triangle_type)
!
    double precision, dimension(:), intent(in) :: A_phi_vec
    integer, dimension(n_triangles,2), intent(out) :: triangle_type
    double precision, dimension(3) :: diff_Aphi_vec
    integer :: i, i_min_diff, vertex_1, vertex_2, vertex_3
!
    do i = 1, n_triangles
!
        vertex_1 = triangles_SOLEDGE3X_EIRENE(i,1)
        vertex_2 = triangles_SOLEDGE3X_EIRENE(i,2)
        vertex_3 = triangles_SOLEDGE3X_EIRENE(i,3)
!
        diff_Aphi_vec(1) = abs(A_phi_vec(vertex_1)-A_phi_vec(vertex_2))
        diff_Aphi_vec(2) = abs(A_phi_vec(vertex_2)-A_phi_vec(vertex_3))
        diff_Aphi_vec(3) = abs(A_phi_vec(vertex_3)-A_phi_vec(vertex_1))
!
        i_min_diff = minloc(diff_Aphi_vec,1)
!
if (i == 372) then
print*, diff_Aphi_vec
end if
!
        select case(i_min_diff)
            case(1)
                if( (A_phi_vec(vertex_1)).gt.(A_phi_vec(vertex_3))) then
                    triangle_type(i,1) = 0
                else
                    triangle_type(i,1) = 1
                endif
                triangle_type(i,2) = 3
            case(2)
                if( (A_phi_vec(vertex_2)).gt.(A_phi_vec(vertex_1))) then
                    triangle_type(i,1) = 0
                else
                    triangle_type(i,1) = 1
                endif
                triangle_type(i,2) = 1
            case(3)
                if( (A_phi_vec(vertex_3)).gt.(A_phi_vec(vertex_2))) then
                    triangle_type(i,1) = 0
                else
                    triangle_type(i,1) = 1
                endif
                triangle_type(i,2) = 2
        end select
    end do
!
end subroutine calc_triangle_type
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine connect_prisms_SOLEDGE3X_EIRENE(prism_1_idx, prism_2_idx, verts, neighbours, neighbour_faces,match)
    integer, intent(in) :: prism_1_idx, prism_2_idx
    integer, dimension(:, :), intent(in) :: verts
    integer, dimension(:, :), intent(inout) :: neighbours, neighbour_faces
    logical, intent(out) :: match
!
    integer :: tetra_1_base, tetra_2_base, tetra_1_idx, tetra_2_idx, tetra_1_off, tetra_2_off, &
               tetra_1_face, tetra_2_face, i
!   
    logical, dimension(4) :: same_vert_1, same_vert_2
!
    match = .false.
!
    tetra_1_base = (prism_1_idx - 1) * 3 + 1
    tetra_2_base = (prism_2_idx - 1) * 3 + 1
!
    do tetra_1_off = 0, 2
        tetra_1_idx = tetra_1_base + tetra_1_off
        do tetra_2_off = 0, 2
            tetra_2_idx = tetra_2_base + tetra_2_off
!
            if (tetra_1_idx == tetra_2_idx) cycle
!
            do i = 1, 4
                same_vert_1(i) = any(verts(:, tetra_2_idx) == verts(i, tetra_1_idx))
                same_vert_2(i) = any(verts(i, tetra_2_idx) == verts(:, tetra_1_idx))
            end do
!
            if (count(same_vert_1) == 3) then
                tetra_1_face = minloc(abs(transfer(same_vert_1 , 1, size=4)), dim=1)
                tetra_2_face = minloc(abs(transfer(same_vert_2 , 1, size=4)), dim=1)
!
! if ((prism_1_idx /= prism_2_idx) .and. &
!     & (neighbour_faces(tetra_1_face, tetra_1_idx) /= -1 .or. (neighbour_faces(tetra_2_face, tetra_2_idx) /= -1))) then
!     print*, 'Try to connect faces that are already link to others -> broken'
!     print*, 'when connecting ', prism_1_idx, ' with ', prism_2_idx
!     stop
! end if
!
                neighbours(tetra_1_face, tetra_1_idx) = tetra_2_idx
                neighbours(tetra_2_face, tetra_2_idx) = tetra_1_idx
!
                neighbour_faces(tetra_1_face, tetra_1_idx) = tetra_2_face
                neighbour_faces(tetra_2_face, tetra_2_idx) = tetra_1_face
!
                match = .true.
!
            end if
        end do
    end do
!
end subroutine connect_prisms_SOLEDGE3X_EIRENE
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    pure function is_not_border(triangle_idx)
        integer, intent(in) :: triangle_idx
        logical :: is_not_border
!
        integer :: n_neighbours
        integer :: i
!
        is_not_border = .false.
        n_neighbours = 0
!
        do i = 1, n_triangles
            if (i == triangle_idx) cycle
!
            if(are_neighbours(triangle_idx,i)) n_neighbours = n_neighbours + 1
!
            if (n_neighbours == 3) then
                is_not_border = .true.
                exit
            end if
!
        end do !n_triangle
!
    end function is_not_border
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    pure function are_neighbours(triangle_1_idx,triangle_2_idx)
        integer, intent(in) :: triangle_1_idx, triangle_2_idx
        logical :: are_neighbours
!
        integer, dimension(3) :: triangle_1, triangle_2
        integer :: n_common_vertices
        integer :: j,k
!
        are_neighbours = .false.
        n_common_vertices = 0
        triangle_1 = triangles_SOLEDGE3X_EIRENE(triangle_1_idx,:)
        triangle_2 = triangles_SOLEDGE3X_EIRENE(triangle_2_idx,:)
!
        do j = 1, 3  
!           
            do k = 1,3
                if (triangle_1(j) == triangle_2(k)) then
                    n_common_vertices = n_common_vertices + 1
                    exit
                end if
            end do !triangle_2 vertices
!
            if (n_common_vertices == 2) then
                are_neighbours = .true.
                exit
            end if
!           
        end do !triangle_1 vertices
!
    end function are_neighbours
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    pure function find_neighbour_triangles(triangle_idx)
        integer, intent(in) :: triangle_idx
        integer, dimension(3) :: find_neighbour_triangles
!
        integer :: n_neighbours, i
!
        n_neighbours = 0
        find_neighbour_triangles = -1
        do i = 1, n_triangles
            if (i == triangle_idx) cycle
!
            if(are_neighbours(triangle_idx,i)) then
                n_neighbours = n_neighbours + 1
                find_neighbour_triangles(n_neighbours) = i
            end if !are_neighbours
!
            if (n_neighbours == 3) then
                exit
            end if
!
        end do !n_triangle
!
    end function find_neighbour_triangles
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine change_verts(cur_triangle,triangle_type,mask_r,mask_phi,mask_theta,verts_per_slice,verts)
        integer, intent(in) :: cur_triangle, verts_per_slice
        integer, dimension(:,:), intent(in) :: triangle_type
        integer, dimension(:,:), intent(in) :: mask_r, mask_phi, mask_theta
        integer, dimension(:,:), intent(inout) :: verts
!
        integer :: prism_orientation, cur_stand_alone_vertex, vertex_I, vertex_II, vertex_III, mask_idx, tetra_idx
        integer, dimension(4, 3) :: slice_offset, ring_offset, cur_triangle_offset, no_offset
!
            prism_orientation = triangle_type(cur_triangle,1)
            cur_stand_alone_vertex = triangle_type(cur_triangle,2)
!
            vertex_I = triangles_SOLEDGE3X_EIRENE(cur_triangle,cur_stand_alone_vertex)
            vertex_II = triangles_SOLEDGE3X_EIRENE(cur_triangle,mod(cur_stand_alone_vertex ,3) + 1)
            vertex_III = triangles_SOLEDGE3X_EIRENE(cur_triangle,mod(cur_stand_alone_vertex + 1,3) + 1)
!      
            ! --- get the correct offset masks for up- or down facing prism ---
            mask_idx = 1 + prism_orientation * 3
!
            slice_offset = mask_phi(:, mask_idx:mask_idx + 2) * verts_per_slice
            ring_offset = abs(mask_r(:, mask_idx:mask_idx + 2)-prism_orientation) &
                        & *abs(mask_theta(:, mask_idx:mask_idx + 2) - 1 + prism_orientation) &
                        & *vertex_II
            cur_triangle_offset = abs(mask_r(:, mask_idx:mask_idx + 2) - prism_orientation) &
                        & *abs(mask_theta(:, mask_idx:mask_idx + 2)-prism_orientation) &
                        & *vertex_III
            no_offset = abs(mask_r(:, mask_idx:mask_idx + 2) - 1 + prism_orientation) &
                        & *abs(mask_theta(:, mask_idx:mask_idx + 2) - 1) &
                        & *vertex_I
            ! --- with offset masks calculate verts ---
            tetra_idx = (cur_triangle - 1) * 3 + 1
            verts(:, tetra_idx:tetra_idx + 2) = no_offset + slice_offset + ring_offset + cur_triangle_offset
!
    end subroutine change_verts
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine repair(count_connected,triangle_type,mask_r,mask_phi,mask_theta,verts_per_slice,verts,neighbours,neighbour_faces)
!
        integer, intent(in) :: verts_per_slice
        integer, dimension(:), intent(in) :: count_connected
        integer, dimension(:,:), intent(in) :: mask_r, mask_phi, mask_theta
        integer, dimension(:,:), intent(inout) :: triangle_type, verts, neighbours, neighbour_faces
!
        integer, dimension(:), allocatable :: error_triangles
        integer, dimension(:,:), allocatable :: neighbour_triangles, auxillary_array
        integer :: n_error, n_repair, k, q, l, m, n, cur_triangle, tetra_idx, n_match
        integer, dimension(3) :: cur_neighbours
        integer, dimension(4,12) :: old_neighbours, old_neighbour_faces
        logical :: match
!
        print*, 'Start repair routine!'
!
        !Not pretty but faster workaround, to not have to determine n_error twice (allocation size of error_triangles -> then fill)
        !n_error can only be maximal the number of not fully connected triangles (so including the proper borders as well)
        n_error = count(count_connected<3)
        allocate(auxillary_array(n_error,4))
!
        n_error = 0
        do k = 1, n_triangles
            if((count_connected(k)<3)) then
                if(is_not_border(k)) then
                    n_error = n_error + 1
                    auxillary_array(n_error,1) = k
                    auxillary_array(n_error,2:4) = find_neighbour_triangles(k)
                end if
            end if
        end do !n_error
!       
        !the actual arrays storing the problem triangles and its neighbours -> deallocate the auxillary array after filling them
        if (allocated(error_triangles)) deallocate(error_triangles)
        if (allocated(neighbour_triangles)) deallocate(neighbour_triangles)
        allocate(error_triangles(n_error),neighbour_triangles(n_error,3))
        error_triangles = auxillary_array(:n_error,1)
        neighbour_triangles = auxillary_array(:n_error,2:4)
        deallocate(auxillary_array)
!
        !We go over all problem triangles seperatly and try to fix them LOCALLY
        n_repair = 0
        PARENT: do k = 1, n_error
            cur_triangle = error_triangles(k)
            cur_neighbours = neighbour_triangles(k,:)
!           
            !save original state of triangle and neighbours just in case the repair DOES fail
            tetra_idx = (cur_triangle-1)*3 + 1
            old_neighbours(:,1:3) = neighbours(:,tetra_idx:tetra_idx+2)
            old_neighbour_faces(:,1:3) = neighbour_faces(:,tetra_idx:tetra_idx+2)
            do q = 1,3
                tetra_idx = (cur_neighbours(q)-1)*3 + 1
                old_neighbours(:,q*3+1:q*3+3) = neighbours(:,tetra_idx:tetra_idx+2)
                old_neighbour_faces(:,q*3+1:q*3+3) = neighbour_faces(:,tetra_idx:tetra_idx+2)
            end do
!
            !Here all possible flip (topface-bottomface) and rotation (change of loose vertex) combinations are tried:
            !By the nature of the prism-cuts there is ALWAYS at least one (rotated) cutting of the prism, that connects with
            !three surrounding (arbitrally cut) surfaces properly (see the figure 2.10 in Bauer et al.)
            FLIP: do l = 0,1
                triangle_type(cur_triangle,1) = modulo(triangle_type(cur_triangle,1)+1,2)
!
                ROTATE: do m = 1, 3
                triangle_type(error_triangles(k),2) = modulo(triangle_type(error_triangles(k),2), 3) + 1
!
                call change_verts(cur_triangle,triangle_type,mask_r,mask_phi,mask_theta,verts_per_slice,verts)
                    CHILD: do n = 1, 3
                        call connect_prisms_SOLEDGE3X_EIRENE(cur_triangle,cur_neighbours(n),verts,neighbours,neighbour_faces,match)
                        if (.not. match) cycle ROTATE
                    end do CHILD
!
                n_repair = n_repair + 1
                cycle PARENT !if all match checks were past -> now it sits right -> go for next problem triangle
!                
                end do ROTATE !permuts through the lose vertex of the triangle and back if no sucess in fitting it
            end do FLIP !flips between top and bottom face and back if not sucess in fitting it
!
            !If repair failed on that triangle (e.g. never hit cycle PARENT) -> reset original state
            tetra_idx = (cur_triangle-1)*3 + 1
            neighbours(:,tetra_idx:tetra_idx+2) = old_neighbours(:,1:3)
            neighbour_faces(:,tetra_idx:tetra_idx+2) = old_neighbour_faces(:,1:3)
            do q = 1,3
                tetra_idx = (cur_neighbours(q)-1)*3 + 1
                neighbours(:,tetra_idx:tetra_idx+2) = old_neighbours(:,q*3+1:q*3+3)
                neighbour_faces(:,tetra_idx:tetra_idx+2) = old_neighbour_faces(:,q*3+1:q*3+3)
            end do
        end do PARENT !n_error
!
        write(*,'(A,I3,A,I3,A)') 'Repaired ', n_repair, ' / ', n_error, ' prisms!'
!
    end subroutine repair
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
end module circular_mesh_SOLEDGE3X_EIRENE
