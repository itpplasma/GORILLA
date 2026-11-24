program test_cyl_vertex_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use constants, only: pi
    use gorilla_settings_mod, only: load_gorilla_inp
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use tetra_grid_mod, only: nvert, verts_rphiz, ntetr
    use tetra_physics_mod, only: tetra_physics, cm_over_e
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla
    implicit none

    real(dp), parameter :: toroidal_field = -20000.0_dp, phi_gauge = 3.0e5_dp
    real(dp), parameter :: radial_potential = 25000.0_dp, step = 2.0e-6_dp
    real(dp) :: x(3), a(3), b(3), velocity, perpendicular, field_norm
    real(dp) :: expected_phi, expected_z, initial_r
    integer :: unit, iv, vertex, cell, face
    logical :: located
    character(len=24) :: mode
    character(len=1024) :: defaults

    call get_command_argument(1, mode)
    call get_command_argument(2, defaults)
    call write_inputs(defaults)
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla(i_option=1)

    open(newunit=unit, file='vertex-field.dat', status='replace')
    write(unit, '(a)') 'GORILLA_CYL_VERTEX_AB_V1'
    write(unit, '(i0)') nvert
    do iv=1,nvert
        vertex = iv
        if (mode == 'order' .and. iv == 1) vertex = 2
        x = verts_rphiz(:,vertex)
        call exact_field(x,a,b)
        if (mode == 'zero-field' .and. iv == 1) b = 0.0_dp
        write(unit, '(9(es25.17e3,1x))') x,a,b
    enddo
    if (mode == 'extra-row') write(unit, '(9(es25.17e3,1x))') x,a,b
    close(unit)
    call initialize_gorilla(i_option=2)

    if (mode /= 'field' .and. mode /= 'orbit') stop
    do cell=1,ntetr
        x = tetra_physics(cell)%x1
        call exact_field(x,a,b)
        field_norm = norm2(b)
        if (abs(tetra_physics(cell)%bmod1/field_norm - 1.0_dp) > 1.0e-14_dp) error stop 'B magnitude'
        if (abs(tetra_physics(cell)%Aphi1/phi_gauge - 1.0_dp) > 1.0e-14_dp) error stop 'A phi covariant'
        if (abs(tetra_physics(cell)%h2_1/x(1) - b(2)/field_norm) > 1.0e-14_dp) error stop 'h phi covariant'
        if (abs(tetra_physics(cell)%h3_1 - b(3)/field_norm) > 1.0e-14_dp) error stop 'h z'
        if (abs(tetra_physics(cell)%curlA(2)/toroidal_field - 1.0_dp) > 1.0e-13_dp) error stop 'A curl'
    enddo

    if (mode == 'orbit') then
        x = [182.0_dp, 0.03_dp, 12.0_dp]
        initial_r = x(1)
        velocity = 4.0e6_dp
        perpendicular = 1.0e6_dp
        located = .false.
        cell = -1
        face = 0
        expected_phi = modulo(x(2) - step*velocity/x(1), 2.0_dp*pi)
        expected_z = x(3) + step*cm_over_e*velocity**2/(toroidal_field*x(1))
        call orbit_timestep_gorilla(x,velocity,perpendicular,step,located,cell,face)
        if (.not. located .or. cell < 0) error stop 'toroidal orbit exited mesh'
        if (abs(x(1)-initial_r) > 1.0e-8_dp) error stop 'toroidal orbit radius'
        if (abs(x(2)-expected_phi) > 1.0e-9_dp) error stop 'toroidal transit across seam'
        if (abs(x(3)-expected_z) > 1.0e-8_dp) error stop 'toroidal curvature drift'
    endif
    print *, 'PASS'

contains

    subroutine exact_field(x, a, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: a(3), b(3)

        a = [0.0_dp, phi_gauge/x(1), -toroidal_field*x(1)]
        b = [0.0_dp, toroidal_field, 0.0_dp]
        if (mode == 'field') then
            a(1) = radial_potential*cos(3.0_dp*x(2))
            b(3) = 3.0_dp*radial_potential*sin(3.0_dp*x(2))/x(1)
        endif
    end subroutine exact_field

    subroutine write_inputs(template)
        character(len=*), intent(in) :: template
        character(len=1024) :: line
        integer :: input, output, io

        open(newunit=input, file=template, status='old')
        open(newunit=output, file='gorilla.inp', status='replace')
        do
            read(input,'(a)',iostat=io) line
            if (io /= 0) error stop 'invalid GORILLA default namelist'
            line = adjustl(line)
            if (line(1:1) == '/') exit
            write(output,'(a)') trim(line)
        enddo
        write(output,'(a)') ' coord_system=1, ispecies=2,'
        if (mode == 'orbit') then
            write(output,'(a)') ' ipusher=1, i_time_tracing_option=1, boole_dt_dtau=.false.,'
        else
            write(output,'(a)') ' ipusher=2, poly_order=4, i_time_tracing_option=2, boole_time_Hamiltonian=.true.,'
        endif
        write(output,'(a)') " vertex_field_file='vertex-field.dat', /"
        close(input)
        close(output)

        open(newunit=output, file='tetra_grid.inp', status='replace')
        write(output,'(a)') '&TETRA_GRID_NML'
        write(output,'(a)') ' grid_kind=5, n1=8, n2=12, n3=16, sfc_s_min=0.02,'
        write(output,'(a)') ' boole_n_field_periods=.true.,'
        write(output,'(a)') ' R0_analytic_circ=160.0, a_analytic_circ=60.0, B0_analytic_circ=19000.0,'
        write(output,'(a)') ' q0_analytic_circ=1.0, q1_analytic_circ=2.0, /'
        close(output)
    end subroutine write_inputs
end program test_cyl_vertex_field
