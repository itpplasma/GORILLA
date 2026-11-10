module vertex_data_mod
!
! Poloidal triangle mesh with field data at its vertices (grid_kind = 7).
! The netCDF file is written e.g. by PREPROCESSING/JOREK/jorek2gorilla.py.
! All quantities are in CGS units and right-handed cylindrical coordinates (R,phi,Z):
!   triangles(three,n_triangles)  one-based vertex indices
!   R, Z                          vertex positions [cm]
!   A_R, A_phi, A_Z               covariant vector potential [G cm, G cm^2, G cm]
!   B_R, B_phi, B_Z, Bmod         physical magnetic field [G]
!   Phi                           electrostatic potential [statV]
!   n, T_i, T_e                   density [cm^-3] and temperatures [eV]
!   mesh_label (optional)         label constant along the poloidal mesh lines, orients the prisms
!                                 (default: A_phi)
!
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
!
    implicit none
    private
!
    integer, public, protected :: n_vertices_2d = 0
    integer, public, protected :: n_triangles_2d = 0
    integer, dimension(:,:), allocatable, public, protected :: triangles_2d
    real(dp), dimension(:), allocatable, public, protected :: R_2d, Z_2d
    real(dp), dimension(:), allocatable, public, protected :: A_R_2d, A_phi_2d, A_Z_2d
    real(dp), dimension(:), allocatable, public, protected :: B_R_2d, B_phi_2d, B_Z_2d, Bmod_2d
    real(dp), dimension(:), allocatable, public, protected :: Phi_2d, n_2d, T_i_2d, T_e_2d
    real(dp), dimension(:), allocatable, public, protected :: mesh_label_2d
!
    public :: load_vertex_data, vertex_data_rphiz, vertex_data_profiles
!
contains
!
    subroutine load_vertex_data(filename)
!
        character(len=*), intent(in) :: filename
        integer :: ncid, varid, dimid
!
        call clear_vertex_data()
        call check_nc(nf90_open(trim(filename), nf90_nowrite, ncid), 'open '//trim(filename))
!
        call check_nc(nf90_inq_dimid(ncid, 'n_vertices', dimid), 'find dimension n_vertices')
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=n_vertices_2d), 'read dimension n_vertices')
        call check_nc(nf90_inq_dimid(ncid, 'n_triangles', dimid), 'find dimension n_triangles')
        call check_nc(nf90_inquire_dimension(ncid, dimid, len=n_triangles_2d), 'read dimension n_triangles')
!
        allocate(triangles_2d(3, n_triangles_2d))
        call check_nc(nf90_inq_varid(ncid, 'triangles', varid), 'find variable triangles')
        call check_nc(nf90_get_var(ncid, varid, triangles_2d), 'read variable triangles')
!
        call read_vertex_variable(ncid, 'R', R_2d)
        call read_vertex_variable(ncid, 'Z', Z_2d)
        call read_vertex_variable(ncid, 'A_R', A_R_2d)
        call read_vertex_variable(ncid, 'A_phi', A_phi_2d)
        call read_vertex_variable(ncid, 'A_Z', A_Z_2d)
        call read_vertex_variable(ncid, 'B_R', B_R_2d)
        call read_vertex_variable(ncid, 'B_phi', B_phi_2d)
        call read_vertex_variable(ncid, 'B_Z', B_Z_2d)
        call read_vertex_variable(ncid, 'Bmod', Bmod_2d)
        call read_vertex_variable(ncid, 'Phi', Phi_2d)
        call read_vertex_variable(ncid, 'n', n_2d)
        call read_vertex_variable(ncid, 'T_i', T_i_2d)
        call read_vertex_variable(ncid, 'T_e', T_e_2d)
        if (nf90_inq_varid(ncid, 'mesh_label', varid) == nf90_noerr) then
            call read_vertex_variable(ncid, 'mesh_label', mesh_label_2d)
        else
            allocate(mesh_label_2d, source=A_phi_2d)
        end if
!
        call check_nc(nf90_close(ncid), 'close '//trim(filename))
!
        if (minval(triangles_2d) < 1 .or. maxval(triangles_2d) > n_vertices_2d) then
            error stop 'vertex_data_mod: triangle vertex index out of range'
        end if
        if (minval(R_2d) <= 0.0_dp) error stop 'vertex_data_mod: vertex with R <= 0'
!
        print *, 'Vertex data: loaded ', n_vertices_2d, ' vertices and ', n_triangles_2d, &
                 ' triangles from ', trim(filename)
!
    end subroutine load_vertex_data
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vertex_data_rphiz(iv_2d, A_r, A_phi, A_z, B_r, B_phi, B_z, bmod, phi_elec)
!
        ! Same output convention as vector_potential_rphiz: B_phi is the covariant component R*B_phi.
!
        integer, intent(in) :: iv_2d
        real(dp), intent(out) :: A_r, A_phi, A_z, B_r, B_phi, B_z, bmod, phi_elec
!
        A_r = A_R_2d(iv_2d)
        A_phi = A_phi_2d(iv_2d)
        A_z = A_Z_2d(iv_2d)
        B_r = B_R_2d(iv_2d)
        B_phi = B_phi_2d(iv_2d)*R_2d(iv_2d)
        B_z = B_Z_2d(iv_2d)
        bmod = Bmod_2d(iv_2d)
        phi_elec = Phi_2d(iv_2d)
!
    end subroutine vertex_data_rphiz
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vertex_data_profiles(ind_tetr, r, z, dens, temp_i, temp_e)
!
        ! Density [cm^-3] and temperatures [eV] at (r,z) inside tetrahedron ind_tetr, linearly
        ! interpolated over the triangle of the prism the tetrahedron belongs to (exact for the
        ! axisymmetric vertex data). Relies on the tetrahedron ordering of calc_mesh_from_triangles:
        ! three tetrahedra per prism, prisms in the order of the triangles, slice after slice.
!
        integer, intent(in) :: ind_tetr
        real(dp), intent(in) :: r, z
        real(dp), intent(out) :: dens, temp_i, temp_e
        integer :: i_triangle, k1, k2, k3
        real(dp) :: det, l1, l2, l3
!
        if (n_triangles_2d == 0) error stop 'vertex_data_profiles: no vertex data loaded (grid_kind = 7 required)'
!
        i_triangle = mod((ind_tetr-1)/3, n_triangles_2d) + 1
        k1 = triangles_2d(1, i_triangle)
        k2 = triangles_2d(2, i_triangle)
        k3 = triangles_2d(3, i_triangle)
!
        det = (R_2d(k2)-R_2d(k1))*(Z_2d(k3)-Z_2d(k1)) - (R_2d(k3)-R_2d(k1))*(Z_2d(k2)-Z_2d(k1))
        l2 = ((r-R_2d(k1))*(Z_2d(k3)-Z_2d(k1)) - (R_2d(k3)-R_2d(k1))*(z-Z_2d(k1)))/det
        l3 = ((R_2d(k2)-R_2d(k1))*(z-Z_2d(k1)) - (r-R_2d(k1))*(Z_2d(k2)-Z_2d(k1)))/det
        l1 = 1.0_dp - l2 - l3
!
        dens = l1*n_2d(k1) + l2*n_2d(k2) + l3*n_2d(k3)
        temp_i = l1*T_i_2d(k1) + l2*T_i_2d(k2) + l3*T_i_2d(k3)
        temp_e = l1*T_e_2d(k1) + l2*T_e_2d(k2) + l3*T_e_2d(k3)
!
    end subroutine vertex_data_profiles
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine read_vertex_variable(ncid, name, values)
!
        integer, intent(in) :: ncid
        character(len=*), intent(in) :: name
        real(dp), dimension(:), allocatable, intent(out) :: values
        integer :: varid
!
        allocate(values(n_vertices_2d))
        call check_nc(nf90_inq_varid(ncid, name, varid), 'find variable '//name)
        call check_nc(nf90_get_var(ncid, varid, values), 'read variable '//name)
!
    end subroutine read_vertex_variable
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine check_nc(status, operation)
!
        integer, intent(in) :: status
        character(len=*), intent(in) :: operation
!
        if (status /= nf90_noerr) then
            write(*, '(A)') trim(operation)//': '//trim(nf90_strerror(status))
            error stop 'vertex_data_mod: NetCDF error'
        end if
!
    end subroutine check_nc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine clear_vertex_data()
!
        if (allocated(triangles_2d)) deallocate(triangles_2d)
        if (allocated(R_2d)) deallocate(R_2d, Z_2d, A_R_2d, A_phi_2d, A_Z_2d, B_R_2d, B_phi_2d, &
                                        B_Z_2d, Bmod_2d, Phi_2d, n_2d, T_i_2d, T_e_2d, mesh_label_2d)
        n_vertices_2d = 0
        n_triangles_2d = 0
!
    end subroutine clear_vertex_data
!
end module vertex_data_mod
