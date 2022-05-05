module circular_mesh
implicit none
private
double precision, parameter :: pi = 3.14159265358979d0
public :: calc_mesh, create_points, calc_points_circular, calc_n_tetras, calc_n_verts, &
& wrap_idx_inplace
integer :: n_field_periods

contains

subroutine create_points(verts_per_ring, n_slices, points_rphiz, points_sthetaphi, efit_vmec,n_field_periods_in,n_verts, &
                         verts_theta_vmec,r_scaling_func, theta_scaling_func,repeat_center_point)
    use points_2d, only: create_points_2d, create_points_2d_vmec, scaling_func

    integer, intent(in) :: n_slices,efit_vmec,n_field_periods_in
    integer, dimension(:), intent(in) :: verts_per_ring ! without venter vert; e.g. (/6, 8, 10/)
    double precision, dimension(:, :), allocatable, intent(out) :: points_rphiz !(r, phi, z)
    double precision, dimension(:, :), allocatable, intent(out) :: points_sthetaphi !(s, theta, phi)
    double precision, dimension(:), allocatable, intent(out),optional :: verts_theta_vmec !(theta in VMEC coordinates)
    integer, intent(out) :: n_verts ! without venter vert; e.g. (/6, 8, 10/)
    procedure(scaling_func), optional :: r_scaling_func, theta_scaling_func
    logical, optional, intent(in) :: repeat_center_point ! set true if we are insym flux coords

    integer :: verts_per_slice,phi_position, i
    double precision :: s,theta,vartheta,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                        R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
    !Write field periodicity in module
    n_field_periods = n_field_periods_in
!                        
    n_verts = calc_n_verts(verts_per_ring, n_slices, repeat_center_point)
    verts_per_slice = n_verts / n_slices

    if (allocated(points_rphiz)) deallocate(points_rphiz)
    if (allocated(points_sthetaphi)) deallocate(points_sthetaphi)
!
    if (present(verts_theta_vmec)) then
        if (allocated(verts_theta_vmec)) deallocate(verts_theta_vmec)
    endif

    allocate(points_rphiz(3, n_verts), points_sthetaphi(3, n_verts))

    !Distinguish inbetween axisymmetric (EFIT) and non-axisymmetric (VMEC) data
!    
    select case(efit_vmec)
      case(1) !EFIT (axisymmetric)
!      
!        call create_points_2d(verts_per_ring, points_rphiz(:, :verts_per_slice), &
        call create_points_2d(1,verts_per_ring, points_rphiz(:, :verts_per_slice), &
                               & points_sthetaphi(:, :verts_per_slice), &
                               & r_scaling_func, theta_scaling_func, repeat_center_point)
        phi_position = 2 !In coordinate system (R,Phi,Z) --> Phi is at second position
        call extrude_points(verts_per_slice, n_slices,phi_position, points_rphiz)
!
        phi_position = 3 !In coordinate system (s,theta,Phi) --> Phi is at third position
        call extrude_points(verts_per_slice, n_slices,phi_position, points_sthetaphi)
!      
      case(2) !VMEC (non-axisymmetric)
!        
        allocate(verts_theta_vmec(n_verts))
!
        call create_points_2d_vmec(verts_per_ring,&
                                &  points_sthetaphi(:, :verts_per_slice), &
                               & r_scaling_func, theta_scaling_func,repeat_center_point)
!
        !For the VMEC field theta is not vartheta !! vartheta = theta + alam
        !But VMEC-theta must be kept in order to compute respective R and Z values
        !VMEC-theta is stored in separate vector verts_theta_vmec
        !Theta value in points_sthetaphi IS SFC theta
!        
        phi_position = 3 !In coordinate system (s,theta,Phi) --> Phi is at third position
        call extrude_points(verts_per_slice, n_slices,phi_position, points_sthetaphi)
!        
        do i = 1,n_verts
!    
            s = points_sthetaphi(1,i)
            theta = points_sthetaphi(2,i)
            varphi = points_sthetaphi(3,i)        
!         
            !Find corresponding VMEC-theta for SFC-theta
            verts_theta_vmec(i) = theta_sym_flux2theta_vmec(s,theta,varphi)
!
            !Read VMEC Spline for Cylindrical variables
            call splint_vmec_data(s,verts_theta_vmec(i),varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                  R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!
            points_rphiz(1,i) = R
            points_rphiz(2,i) = varphi
            points_rphiz(3,i) = Z
!          
        enddo
!
    end select    
        
!    open(123, file='points_rphiz.dat')
!    do i=1, verts_per_slice
!        write(123,101) points_rphiz(:, i)
!    end do
!    close(123)
!
!    101 format(1000(e21.14,x))

end subroutine create_points

subroutine calc_points_circular(verts_per_ring, n_slices, R0, Rmax, Z0,  points, repeat_center_point)

    double precision, intent(in) :: R0, Z0, Rmax
    integer, intent(in) :: n_slices
    integer, dimension(:), intent(in) :: verts_per_ring ! without venter vert; e.g. (/6, 8, 10/)
    double precision, dimension(:, :), intent(inout) :: points !(r, phi, z)
    logical, optional, intent(in) :: repeat_center_point ! set true if we are insym flux coords

    integer :: vert, ring, n_rings, n_verts, verts_per_slice, vert_idx, phi_position, n_verts_current
    double precision :: r, theta

    integer :: n_center_point

    n_center_point = 1
    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_center_point = verts_per_ring(1)
        end if
    end if

    n_rings = size(verts_per_ring)
    n_verts = calc_n_verts(verts_per_ring, n_slices, repeat_center_point)
    verts_per_slice = n_verts / n_slices

    vert_idx = 1
    do ring = 1, n_rings
        if (ring == 0) then
            n_verts_current = n_center_point
        else
            n_verts_current = verts_per_ring(ring)
        end if
        r = (Rmax * ring) / n_rings
        do vert = 1, n_verts_current
            theta = (2.d0 * pi * (vert - 1)) / n_verts_current
            points(:, vert_idx) = r*(/cos(theta), 0.d0, sin(theta)/)
            vert_idx = vert_idx + 1
        end do
    end do

    points(1, :verts_per_slice) = points(1, :verts_per_slice) + R0
    points(3, :verts_per_slice) = points(3, :verts_per_slice) + Z0

    phi_position = 2 !In coordinate system (R,Phi,Z) --> Phi is at second position
    call extrude_points(verts_per_slice, n_slices,phi_position, points)

end subroutine calc_points_circular

subroutine extrude_points(verts_per_slice, n_slices,phi_position,  points)
    integer, intent(in) :: verts_per_slice, n_slices,phi_position
    double precision, dimension(:, :), intent(inout) :: points !(r, phi, z)

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

function theta_vmec2theta_sym_flux(s,theta_vmec,varphi) result(theta_sym_flux) 
!
        implicit none
!        
        double precision :: s,theta_vmec,theta_sym_flux,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                            R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp
!
        call splint_vmec_data(s,theta_vmec,varphi,A_phi,A_theta,dA_phi_ds,dA_theta_ds,aiota,       &
                                R,Z,alam,dR_ds,dR_dt,dR_dp,dZ_ds,dZ_dt,dZ_dp,dl_ds,dl_dt,dl_dp)
!                                  
        theta_sym_flux = theta_vmec + alam
!
end function theta_vmec2theta_sym_flux
!
function theta_sym_flux2theta_vmec(s,theta_sym_flux,varphi) result(theta_vmec)
!
  implicit none
!
  double precision :: s,varphi,theta_sym_flux
  double precision :: theta_vmec
  double precision :: dl_dt,alam,deltheta
  integer :: iter
  double precision, parameter :: epserr = 1.d-14
!
!
  ! Begin Newton iteration to find VMEC theta
!
  theta_vmec = theta_sym_flux
!
  do iter=1,100
!
    call splint_lambda(s,theta_vmec,varphi,alam,dl_dt)
!
    deltheta = (theta_sym_flux-theta_vmec-alam)/(1.d0+dl_dt)
    theta_vmec = theta_vmec + deltheta
    if(abs(deltheta).lt.epserr) exit
  enddo
!
end function theta_sym_flux2theta_vmec
!
subroutine calc_mesh(verts_per_ring, n_slices, points, n_tetras, & 
                     verts, neighbours, neighbour_faces, perbou_phi, perbou_theta, repeat_center_point)
    integer, intent(in) :: n_slices
    integer, dimension(:), intent(in) :: verts_per_ring! without venter vert; e.g. (/6, 8, 10/)
    double precision, dimension(:, :), intent(in) :: points !points in first (r, phi, z) 
    integer, intent(out) :: n_tetras
    integer, dimension(:, :), allocatable, intent(out) :: verts, neighbours, neighbour_faces, perbou_phi, perbou_theta
    logical, optional, intent(in) :: repeat_center_point ! set true if we are insym flux coords

    integer, dimension(4, 6) :: tetra_conf, mask_theta, mask_phi, mask_r
    integer, dimension(4, 3) :: slice_offset, ring_offset, segment_offset
    integer, dimension(size(verts_per_ring)) :: prisms_per_ring
    integer, allocatable, dimension(:) :: top_facing_prisms

    integer :: n_rings, n_verts, verts_per_slice, tetras_per_slice, &
               mask_idx, n_verts_lower, n_verts_upper, prism_orientation, &
               slice, ring, segment, base_idx, prism_idx, tetra_idx, lower_off, upper_off, i, f

    double precision, dimension(2) :: u, v, p, q
    integer :: n_center_point

    n_center_point = 1
    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_center_point = verts_per_ring(1)
        end if
    end if

    n_rings = size(verts_per_ring)
    
    if (any(verts_per_ring < 3) .or. n_rings < 1 .or. n_slices < 3) then
        print *, "Invalid parameters for function calc_mesh"
        stop
    end if
    
    verts_per_slice = sum(verts_per_ring) + n_center_point ! verts per slice
    n_verts = verts_per_slice * n_slices
    
    n_tetras = calc_n_tetras(verts_per_ring, n_slices, repeat_center_point)
    if (.not. size(verts, dim=2) == n_tetras .and. size(neighbours, dim=2) == n_tetras &
            .and. size(neighbour_faces, dim=2) == n_tetras) then
        print *, "Invalid size of output arrays for function calc_mesh"
        stop
    end if

    tetras_per_slice = n_tetras / n_slices
    
    prisms_per_ring = verts_per_ring
    prisms_per_ring(2:) = prisms_per_ring(2:) + verts_per_ring(:n_rings - 1)
    if (n_center_point /= 1) then
        prisms_per_ring(1) = prisms_per_ring(1) * 2
    end if

    if (allocated(verts)) deallocate(verts)
    if (allocated(neighbours)) deallocate(neighbours)
    if (allocated(neighbour_faces)) deallocate(neighbour_faces)
    if (allocated(perbou_phi)) deallocate(perbou_phi)
    if (allocated(perbou_theta)) deallocate(perbou_theta)
    allocate(verts(4, n_tetras), neighbours(4, n_tetras), neighbour_faces(4, n_tetras), &
             perbou_phi(4, n_tetras), perbou_theta(4, n_tetras))
    allocate(top_facing_prisms(verts_per_slice - 1))

    ! 
    !     7------6
    !    /|     /|
    !   3------2 |
    !   | 5----|-4           r     
    !   |/     |/            ^ ,phi
    !   1------0 *           |/    
    !               theta <--*

    ! top facing prism
    ! top edge from 2 - 7
    tetra_conf(:, :3) = reshape((/ &
        0, 2, 3, 7,  & !1 
        0, 2, 6, 7,  & !1
        0, 4, 6, 7  & !1
        /), (/4, 3/))
    
    ! calculate the bottom facing prism tetraeder configurations from the top facing one one
    ! bottom facing prism
    tetra_conf(:, 4:) = tetra_conf(:, :3)
    where (iand(tetra_conf(:, 4:) , 2) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 1)
    end where
    where (iand(tetra_conf(:, 4:) , 1) > 0)
        tetra_conf(:, 4:) = ieor(tetra_conf(:, 4:) , 2)
    end where

    do i = 1, 6
        print *, tetra_conf(:, i)
        if (mod(i, 3) == 0) print *, ''
    end do

    ! extract masks from tetra_conf
    mask_theta = iand(tetra_conf, 1)
    mask_r = ishft(iand(tetra_conf, 2), -1)
    mask_phi = ishft(iand(tetra_conf, 4), -2)

    ! initialize neighbour_faces
    neighbour_faces = - 1

    ! initialize periodic boundary for theta
    perbou_theta = 0

    ! first slice
    prism_idx = 1
    tetra_idx = 1
    base_idx = 1
    do ring = 1, n_rings
        upper_off = 0
        lower_off = 0
        n_verts_upper = verts_per_ring(ring)
        
        if (ring == 1) then
            n_verts_lower = n_center_point
        else
            n_verts_lower = verts_per_ring(ring - 1)
        end if

        do segment = 1, prisms_per_ring(ring)

            u = points((/1, 3/), base_idx + modulo(lower_off, n_verts_lower))
            v = points((/1, 3/), base_idx + n_verts_lower + modulo(upper_off, n_verts_upper))
            p = points((/1, 3/), base_idx + n_verts_lower + modulo(upper_off + 1, n_verts_upper))
            q = points((/1, 3/), base_idx + modulo(lower_off + 1, n_verts_lower))
            
            ! if ((q(1) - v(1))*(p(2) - v(2)) - (q(2) - v(2))*(p(1) - v(1)) > 0.d0) then
            !     print *, "ERROR: points in outer ring are too far apart or meshed in wrong coordinate system"
            !     stop
            ! end if
            
            ! --- calculate prism orientation ---
            if (verts_per_ring(ring) == n_verts_lower) then
                prism_orientation = mod((segment - 1), 2)
            else if (lower_off > n_verts_lower) then
                prism_orientation = 0
            else if (upper_off > n_verts_upper) then
                prism_orientation = 1
            else
                if ((ring == 1 .and. n_center_point == 1).or. delaunay_condition(u, v, p, q, .true.)) then
                    ! triangle [0, 2, 3] satisfies delaunay condition => prism is up facing
                    prism_orientation  = 0
                else
                    ! prism is down facing
                    prism_orientation = 1
                end if
            end if

            
            ! --- get the correct offset masks for up- or down facing prism ---
            mask_idx = 1 + prism_orientation * 3

            slice_offset = mask_phi(:, mask_idx:mask_idx + 2) * verts_per_slice
            ring_offset = mask_r(:, mask_idx:mask_idx + 2)
            segment_offset = mask_theta(:, mask_idx:mask_idx + 2)

            where (ring_offset /= 0) segment_offset = modulo(segment_offset + upper_off, n_verts_upper)
            where (ring_offset == 0) segment_offset = modulo(segment_offset + lower_off, n_verts_lower)
            ring_offset = ring_offset * n_verts_lower
            
            ! --- with offset masks calculate verts ---
            verts(:, tetra_idx:tetra_idx + 2) = base_idx + slice_offset + ring_offset + segment_offset

            ! --- connect neighbouring tetras  ---
            ! connect neighbouring tetras in this prism
            call connect_prisms(prism_idx, prism_idx, verts, neighbours, neighbour_faces)

            ! connect with prisms in neighbouring slices
            neighbours(4, tetra_idx) = tetra_idx + 2 - tetras_per_slice
            neighbour_faces(4, tetra_idx) = 1
            neighbours(1, tetra_idx + 2) = tetra_idx + tetras_per_slice
            neighbour_faces(1, tetra_idx + 2) = 4

            ! connect with previous prism in the same ring
            if (.not. segment == 1) then
                call connect_prisms(prism_idx, prism_idx - 1, verts, neighbours, neighbour_faces)
            end if
            
            if (prism_orientation == 0) then ! top facing prism
                ! add idx of this prism to array of top facing prisms
                top_facing_prisms(base_idx + (n_verts_lower - 1) + upper_off) = prism_idx
                
                ! periodic boundary for theta
                if (segment == 1) then
                    perbou_theta(4, tetra_idx + 1) = -1
                    perbou_theta(4, tetra_idx + 2) = -1
                else if (segment == prisms_per_ring(ring)) then
                    perbou_theta(2, tetra_idx) = 1
                    perbou_theta(3, tetra_idx + 2) = 1
                end if
                
                upper_off = upper_off + 1
            else ! bottom facing prism
                ! connect with prism on previous ring
                if (.not. ring == 1) then
                    call connect_prisms(prism_idx, top_facing_prisms(base_idx - 1 + lower_off), &
                                        verts, neighbours, neighbour_faces)
                end if
                
                ! periodic boundary for theta
                if (segment == 1) then
                    perbou_theta(2, tetra_idx) = -1
                    perbou_theta(3, tetra_idx + 2) = -1
                else if (segment == prisms_per_ring(ring)) then
                    perbou_theta(1, tetra_idx) = 1
                    perbou_theta(1, tetra_idx + 1) = 1
                end if
                
                lower_off = lower_off + 1 
            end if
            
            prism_idx = prism_idx + 1
            tetra_idx = (prism_idx - 1) * 3 + 1

        end do

        ! connect first and last prism in ring
        call connect_prisms(prism_idx - 1, prism_idx - prisms_per_ring(ring), verts, neighbours, neighbour_faces)

        base_idx = base_idx + n_verts_lower

    end do

    deallocate(top_facing_prisms)

    ! for all the other slices we can calculate the verts by incrementing the
    ! vert indices of the first(idx = 0) slice by slice * verts_per_slice
    ! and the neighbours by incrementing the neighbours of the first slice by slice * tetras_per_slice
    ! neighbour_faces and perbou_theta can just be copied

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = verts(:, 1:tetras_per_slice) + slice * verts_per_slice

        ! wrap indices around in the last slice
        if (slice == n_slices - 1) then
            call wrap_idx_inplace(verts(:, tetra_idx:tetra_idx + tetras_per_slice - 1), n_verts)
        end if
    end do

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

    do slice = 1, n_slices - 1 ! slice is the 0-based index of the slice
        tetra_idx = 1 + slice * tetras_per_slice

        perbou_theta(:, tetra_idx:tetra_idx + tetras_per_slice - 1) = perbou_theta(:, 1:tetras_per_slice)
    end do

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
                if (.not. perbou_theta(neighbour_faces(f, i), neighbours(f, i)) == - perbou_theta(f, i)) then
                    print *, "!!! CONSISTENCY CHECK FOR PERBOU THETA FAILED, MESH IS BROKEN !!!"
                    stop
                end if
            end if
        end do
    end do
    print *, "mesh consistency check ok"

end subroutine calc_mesh

pure function calc_n_verts(verts_per_ring, n_slices, repeat_center_point) result(n_verts)
    ! calculate the number of vertices in a torus given verts_per_ring and n_slices
    integer, dimension(:), intent(in) :: verts_per_ring
    integer, intent(in) :: n_slices
    logical, optional, intent(in) ::repeat_center_point

    integer :: n_verts

    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_verts = (sum(verts_per_ring) + verts_per_ring(1)) * n_slices
            return
        end if
    end if

    n_verts = (sum(verts_per_ring) + 1) * n_slices

end function calc_n_verts

pure function calc_n_tetras(verts_per_ring, n_slices, repeat_center_point) result(n_tetras)
    ! calculate the number of tetraeders in a torus given verts_per_ring, n_rings and n_slices
    integer, dimension(:), intent(in) :: verts_per_ring
    integer, intent(in) :: n_slices
    logical, optional, intent(in) ::repeat_center_point

    integer :: n_tetras

    if (present(repeat_center_point)) then
        if (repeat_center_point .eqv. .true.) then
            n_tetras = (verts_per_ring(1) + &
                       sum(verts_per_ring) + sum(verts_per_ring(:size(verts_per_ring) - 1))) * 3 * n_slices
            return
        end if
    end if
        
    n_tetras = (sum(verts_per_ring) + sum(verts_per_ring(:size(verts_per_ring) - 1))) * 3 * n_slices

end function calc_n_tetras

pure function delaunay_condition(u, v, p, q, fixed_order) result(valid)
    ! check if the triangle [u, v, p] satisfies the delaunay condition with respect to point q
    ! returns: 
    !   valid: .true. the triangle [u, v, p] satisfies the delaunay condition 
    !          .false. otherwise

    double precision, dimension(2), intent(in) :: u, v, p, q
    logical, intent(in), optional :: fixed_order
    logical :: valid

    double precision :: a, b, c, d, e, f, g, h, i
    double precision :: delta, gamma

    a = u(1) - q(1)
    b = u(2) - q(2)
    c = (u(1) - q(1)) ** 2.d0 + (u(2) - q(2)) ** 2.d0
    d = v(1) - q(1)
    e = v(2) - q(2)
    f = (v(1) - q(1)) ** 2.d0 + (v(2) - q(2)) ** 2.d0
    g = p(1) - q(1)
    h = p(2) - q(2)
    i = (p(1) - q(1)) ** 2.d0 + (p(2) - q(2)) ** 2.d0

    delta = a*e*i + b*f*g + c*d*h - c*e*g - b*d*i - a*f*h

    if (.not.present(fixed_order) .or. fixed_order .eqv. .false.) then
        ! needed for result to be invariant to order of points
        gamma = ((u(1) - p(1)) * (v(2) - p(2))) - ((v(1) - p(1)) * (u(2) - p(2)))
    else 
        ! if the order of the points is counterclockwise we dont need to calculate gamma
        gamma = 1.d0
    end if

    valid = (delta * gamma < 0.0d0)

end function delaunay_condition

subroutine connect_prisms(prism_1_idx, prism_2_idx, verts, neighbours, neighbour_faces)
    integer, intent(in) :: prism_1_idx, prism_2_idx
    integer, dimension(:, :), intent(in) :: verts
    integer, dimension(:, :), intent(inout) :: neighbours, neighbour_faces

    integer :: tetra_1_base, tetra_2_base, tetra_1_idx, tetra_2_idx, tetra_1_off, tetra_2_off, &
               tetra_1_face, tetra_2_face, i
    
    logical, dimension(4) :: same_vert_1, same_vert_2

    tetra_1_base = (prism_1_idx - 1) * 3 + 1
    tetra_2_base = (prism_2_idx - 1) * 3 + 1

    do tetra_1_off = 0, 2
        tetra_1_idx = tetra_1_base + tetra_1_off
        do tetra_2_off = 0, 2
            tetra_2_idx = tetra_2_base + tetra_2_off

            if (tetra_1_idx == tetra_2_idx) cycle

            do i = 1, 4
                same_vert_1(i) = any(verts(:, tetra_2_idx) == verts(i, tetra_1_idx))
                same_vert_2(i) = any(verts(i, tetra_2_idx) == verts(:, tetra_1_idx))
            end do

            if (count(same_vert_1) == 3) then
                tetra_1_face = minloc(abs(transfer(same_vert_1 , 1, size=4)), dim=1)
                tetra_2_face = minloc(abs(transfer(same_vert_2 , 1, size=4)), dim=1)

                neighbours(tetra_1_face, tetra_1_idx) = tetra_2_idx
                neighbours(tetra_2_face, tetra_2_idx) = tetra_1_idx

                neighbour_faces(tetra_1_face, tetra_1_idx) = tetra_2_face
                neighbour_faces(tetra_2_face, tetra_2_idx) = tetra_1_face
            end if
        end do
    end do

end subroutine connect_prisms

elemental function wrap_idx(index, period) result(wrapped_index)
    ! wrap around the 1-based periodic index which may be larger or smaller than the period
    integer, intent(in) :: index
    integer, intent(in) :: period
    integer :: wrapped_index
    wrapped_index = modulo(index - 1, period) + 1
end function wrap_idx

elemental subroutine wrap_idx_inplace(index, period)
    integer, intent(inout) :: index
    integer, intent(in) :: period
    index = wrap_idx(index, period)
end subroutine wrap_idx_inplace

end module circular_mesh

! program test
!     use circular_mesh, only: calc_mesh, calc_points_circular, calc_n_tetras, calc_n_verts
!     implicit none
!       
!     integer, allocatable, dimension(:, :) :: v, n, nf, pb
!     double precision, allocatable, dimension(:, :) :: points
!     integer :: n_tetras, n_points, i, f
!     
!     integer :: n_slices = 3
!     !integer, dimension(3) :: vpr = (/6, 9, 12/)
!     integer, dimension(3) :: vpr
!     vpr = (/(i, i=6, 6 + size(vpr) * 2 - 1, 2)/)
! 
!     n_tetras = calc_n_tetras(vpr, n_slices)
!     n_points = calc_n_verts(vpr, n_slices)
! 
!     allocate(v(4, n_tetras))
!     allocate(n(4, n_tetras))
!     allocate(nf(4, n_tetras))
!     allocate(pb(4, n_tetras))
!     allocate(points(3, n_points))
! 
!     call calc_points_circular(vpr, n_slices, 171.d0, 96.d0, 0.d0, points)
!     call calc_mesh(vpr, n_slices, points(:, :n_points / n_slices), v, n, nf, pb)
!     !do i = 1, size(points, 2)
!     !    print *, i, points(:, i)
!     !end do
!     ! simple consistency check
! 
!     do i = 1, size(v, 2)
!         do f = 1, 4
!             if (.not. (n(f, i) == -1 .and. nf(f, i) == -1)) then
!                 if (.not. n(nf(f, i), n(f, i)) == i) then
!                     print *, "!!! CONCISTENCY CHECK FAILED, MESH IS BROKEN !!!"
!                     stop
!                 end if
!             end if
!         end do
!         !print *, v(:, i)
!         !if (mod(i, 3) == 0) print *,
!     end do
!     print *, "consistency check ok"
! end program test
! 
