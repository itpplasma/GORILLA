module test_guiding_centre_step_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  use guiding_centre_step_mod
  use funit
  implicit none

  ! Background:
  ! GORILLA's charged guiding-centre pusher (orbit_timestep_gorilla) is
  ! exposed through the species-agnostic characteristic-step facade in
  ! guiding_centre_step_mod.  These tests
  !   * verify that invalid state dimensions, non-finite fields and
  !     unsupported characteristic models fail explicitly (mesh-free), and
  !   * reproduce the legacy charged guiding-centre endpoint and event
  !     sequence through the facade on a small analytic circular tokamak
  !     (grid_kind = 5), comparing the facade's result with the legacy
  !     orbit_timestep_gorilla stepping.

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Mesh-free rejection tests: validation happens before any mesh use.
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  @test
  subroutine test_reject_wrong_state_dimension_low()
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(4)
    integer :: ierr
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    state = 0.0_dp
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_INVALID, ierr)
  end subroutine test_reject_wrong_state_dimension_low

  @test
  subroutine test_reject_wrong_state_dimension_high()
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(6)
    integer :: ierr
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    state = 0.0_dp
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_INVALID, ierr)
  end subroutine test_reject_wrong_state_dimension_high

  @test
  subroutine test_reject_non_finite_position()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    state = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    state(1) = ieee_value(1.0_dp, ieee_quiet_nan)   ! non-finite position
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_INVALID, ierr)
  end subroutine test_reject_non_finite_position

  @test
  subroutine test_reject_non_finite_time()
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    state = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call trace_characteristic_step(state, ieee_value(1.0_dp, ieee_quiet_nan), eq, res, ierr)
    @assertEqual(CHAR_ERR_INVALID, ierr)
  end subroutine test_reject_non_finite_time

  @test
  subroutine test_reject_unsupported_model_neutral()
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr
    eq%model = CHAR_MODEL_NEUTRAL
    state = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_UNSUPPORTED, ierr)
  end subroutine test_reject_unsupported_model_neutral

  @test
  subroutine test_reject_unsupported_model_fluid()
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr
    eq%model = CHAR_MODEL_FLUID
    state = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_UNSUPPORTED, ierr)
  end subroutine test_reject_unsupported_model_fluid

  @test
  subroutine test_reject_unsupported_model_unknown()
    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr
    eq%model = 99
    state = [0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call trace_characteristic_step(state, 1.0_dp, eq, res, ierr)
    @assertEqual(CHAR_ERR_UNSUPPORTED, ierr)
  end subroutine test_reject_unsupported_model_unknown

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Legacy endpoint / event-sequence reproduction on a small analytic
  ! circular tokamak (grid_kind = 5, rectangular grid, no netcdf).
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  @test
  subroutine test_reproduces_legacy_guiding_centre_endpoint()
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use gorilla_settings_mod, only: load_gorilla_inp
    use orbit_timestep_gorilla_mod, only: initialize_gorilla, orbit_timestep_gorilla

    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp), parameter :: eps = 1.0d-10

    ! Initial charged guiding-centre phase space (R, phi, Z, vpar, vperp).
    real(dp) :: x0(3), vpar0, vperp0
    real(dp) :: t_step

    ! Legacy stepping bookkeeping.
    real(dp) :: x_leg(3), vpar_leg, vperp_leg, t_remain_leg
    logical :: boole_leg
    integer :: ind_tetr_leg, iface_leg

    ! Facade stepping bookkeeping.
    real(dp) :: state(CHAR_STATE_DIM)
    logical :: boole_fac
    integer :: ind_tetr_fac, iface_fac
    integer :: ierr_local
    integer :: k

    ! --- setUp: analytic circular tokamak mesh (no netcdf required) ---
    call write_analytic_input_files()
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()

    x0 = [172.0d0, 0.1d0, 0.0d0]
    vpar0 = 1.0d8
    vperp0 = 2.0d7
    t_step = 1.0d-8

    ! --- legacy charged guiding-centre sequence (2 steps) ---
    x_leg = x0
    vpar_leg = vpar0
    vperp_leg = vperp0
    boole_leg = .false.
    ind_tetr_leg = -1
    iface_leg = -1
    do k = 1, 2
      call orbit_timestep_gorilla(x_leg, vpar_leg, vperp_leg, t_step, &
                                  boole_leg, ind_tetr_leg, iface_leg, t_remain_leg)
    end do

    ! --- facade sequence from the same start (2 steps) ---
    state = [x0(1), x0(2), x0(3), vpar0, vperp0]
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    boole_fac = .false.
    ind_tetr_fac = -1
    iface_fac = -1
    do k = 1, 2
      call trace_characteristic_step(state, t_step, eq, res, ierr_local, &
                                     boole_initialized=boole_fac, ind_tetr=ind_tetr_fac, iface=iface_fac)
      @assertEqual(CHAR_OK, ierr_local)
    end do

    ! --- the facade must reproduce the legacy endpoint ---
    @assertEqual(x_leg(1), state(1), tolerance = eps)
    @assertEqual(x_leg(2), state(2), tolerance = eps)
    @assertEqual(x_leg(3), state(3), tolerance = eps)
    @assertEqual(vpar_leg, state(4), tolerance = eps)
    @assertEqual(vperp_leg, state(5), tolerance = eps)

    ! --- event sequence: both steps finished inside the same cell ---
    @assertTrue(res%finished)
    @assertEqual(CHAR_EVENT_NONE, res%event%kind)
    @assertEqual(ind_tetr_leg, res%event%tetrahedron)
    @assertEqual(iface_leg, res%event%face)
    @assertEqual(t_step, res%time, tolerance = 1.0d-12)
    @assertEqual(0.0_dp, res%event%remaining, tolerance = 1.0d-10)
    @assertEqual(t_remain_leg, res%event%remaining, tolerance = eps)
  end subroutine test_reproduces_legacy_guiding_centre_endpoint

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! A wall (material-boundary) event must report the exit tetrahedron,
  ! the crossed boundary face and its outward normal -- the contract the
  ! facade promises for CHAR_EVENT_WALL.  A guiding centre launched near
  ! the low-field side outer wall (R0+a = 250) drifts outward and leaves
  ! the domain there.
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  @test
  subroutine test_wall_event_reports_boundary_face_and_normal()
    use tetra_grid_settings_mod, only: load_tetra_grid_inp
    use gorilla_settings_mod, only: load_gorilla_inp
    use orbit_timestep_gorilla_mod, only: initialize_gorilla

    type(characteristic_equations_t) :: eq
    type(characteristic_result_t) :: res
    real(dp), parameter :: eps = 1.0d-6
    real(dp) :: state(CHAR_STATE_DIM)
    integer :: ierr_local
    logical :: boole
    integer :: ind_tetr, iface

    ! setUp: analytic circular tokamak mesh (no netcdf required).
    call write_analytic_input_files()
    call load_tetra_grid_inp()
    call load_gorilla_inp()
    call initialize_gorilla()

    ! Start near the outer wall (Rmax = R0 + a = 250) on the low-field side;
    ! the drift carries the guiding centre across the outer boundary.
    state = [220.0_dp, 3.0_dp, -70.0_dp, 1.0d8, 2.0d7]
    eq%model = CHAR_MODEL_GUIDING_CENTRE
    boole = .false.
    ind_tetr = -1
    iface = -1

    call trace_characteristic_step(state, 1.0d-3, eq, res, ierr_local, &
                                   boole_initialized=boole, ind_tetr=ind_tetr, iface=iface)

    @assertEqual(CHAR_OK, ierr_local)
    @assertEqual(CHAR_EVENT_WALL, res%event%kind)
    @assertFalse(res%finished)
    ! The boundary (exit) tetrahedron/face must be preserved (not -1).
    @assertTrue(res%event%tetrahedron > 0)
    @assertTrue(res%event%face >= 1 .and. res%event%face <= 4)
    ! The particle ends on the outer wall and the outward normal points +R.
    @assertEqual(250.0_dp, res%position(1), tolerance = eps)
    @assertEqual(1.0_dp, res%event%normal(1), tolerance = eps)
    @assertEqual(0.0_dp, res%event%normal(2), tolerance = eps)
    @assertEqual(0.0_dp, res%event%normal(3), tolerance = eps)
    @assertEqual(1.0_dp, norm2(res%event%normal), tolerance = eps)
  end subroutine test_wall_event_reports_boundary_face_and_normal

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Write the grid_kind=5 analytic circular tokamak input files that
  ! initialize_gorilla needs in the working directory.
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine write_analytic_input_files()
    integer :: iunit

    ! --- gorilla.inp: cylindrical coordinates, electron, RK pusher ---
    open(newunit=iunit, file='gorilla.inp', status='replace', action='write')
    write(iunit, '(a)') '&GORILLANML'
    write(iunit, '(a)') ' eps_Phi =  0.d0 ,'
    write(iunit, '(a)') ' coord_system = 1 ,'
    write(iunit, '(a)') ' ispecies = 1 ,'
    write(iunit, '(a)') ' boole_periodic_relocation = .true. ,'
    write(iunit, '(a)') ' ipusher = 1 ,'
    write(iunit, '(a)') ' boole_pusher_ode45 = .false. ,'
    write(iunit, '(a)') ' rel_err_ode45 = 1.E-8 ,'
    write(iunit, '(a)') ' boole_dt_dtau = .true. ,'
    write(iunit, '(a)') ' boole_newton_precalc = .false. ,'
    write(iunit, '(a)') ' poly_order = 2 ,'
    write(iunit, '(a)') ' i_precomp = 0 ,'
    write(iunit, '(a)') ' boole_guess = .false. ,'
    write(iunit, '(a)') ' handover_processing_kind = 1 ,'
    write(iunit, '(a)') ' boole_axi_noise_vector_pot = .false. ,'
    write(iunit, '(a)') ' boole_axi_noise_elec_pot = .false. ,'
    write(iunit, '(a)') ' boole_non_axi_noise_vector_pot = .false. ,'
    write(iunit, '(a)') ' axi_noise_eps_A = 0.d0 ,'
    write(iunit, '(a)') ' axi_noise_eps_Phi = 0.d0 ,'
    write(iunit, '(a)') ' non_axi_noise_eps_A = 0.d0 ,'
    write(iunit, '(a)') ' boole_helical_pert = .false. ,'
    write(iunit, '(a)') ' helical_pert_eps_Aphi = 0.d0 ,'
    write(iunit, '(a)') ' helical_pert_m_fourier = 0 ,'
    write(iunit, '(a)') ' helical_pert_n_fourier = 0 ,'
    write(iunit, '(a)') ' boole_time_Hamiltonian = .false. ,'
    write(iunit, '(a)') ' boole_gyrophase = .false. ,'
    write(iunit, '(a)') ' boole_vpar_int = .true. ,'
    write(iunit, '(a)') ' boole_vpar2_int = .true. ,'
    write(iunit, '(a)') ' boole_adaptive_time_steps = .false. ,'
    write(iunit, '(a)') ' desired_delta_energy = 0.d0 ,'
    write(iunit, '(a)') ' max_n_intermediate_steps = 0 ,'
    write(iunit, '(a)') ' boole_grid_for_find_tetra = .false. ,'
    write(iunit, '(a)') ' a_factor = 1 ,'
    write(iunit, '(a)') ' b_factor = 1 ,'
    write(iunit, '(a)') ' c_factor = 1 ,'
    write(iunit, '(a)') ' i_time_tracing_option = 0 ,'
    write(iunit, '(a)') ' boole_strong_electric_field = .false. ,'
    write(iunit, '(a)') ' boole_save_electric = .false. ,'
    write(iunit, '(a)') " filename_electric_field = 'E_field.dat' ,"
    write(iunit, '(a)') " filename_electric_drift = 'v_E_drift.dat' ,"
    write(iunit, '(a)') ' boole_pert_from_mephit = .false. ,'
    write(iunit, '(a)') '/'
    close(iunit)

    ! --- tetra_grid.inp: small analytic circular tokamak grid ---
    open(newunit=iunit, file='tetra_grid.inp', status='replace', action='write')
    write(iunit, '(a)') '&TETRA_GRID_NML'
    write(iunit, '(a)') ' n1 = 8 ,'
    write(iunit, '(a)') ' n2 = 6 ,'
    write(iunit, '(a)') ' n3 = 8 ,'
    write(iunit, '(a)') ' grid_kind = 5 ,'
    write(iunit, '(a)') ' boole_n_field_periods = .true. ,'
    write(iunit, '(a)') ' n_field_periods_manual = 1 ,'
    write(iunit, '(a)') ' sfc_s_min = 1.d-1 ,'
    write(iunit, '(a)') ' theta_geom_flux = 1 ,'
    write(iunit, '(a)') ' theta0_at_xpoint = .true. ,'
    write(iunit, '(a)') ' boole_write_mesh_obj = .false. ,'
    write(iunit, '(a)') ' R0_analytic_circ = 170.0d0 ,'
    write(iunit, '(a)') ' a_analytic_circ = 80.0d0 ,'
    write(iunit, '(a)') ' B0_analytic_circ = 10000.0d0 ,'
    write(iunit, '(a)') ' q0_analytic_circ = 1.0d0 ,'
    write(iunit, '(a)') ' q1_analytic_circ = 0.0d0 ,'
    write(iunit, '(a)') '/'
    close(iunit)
  end subroutine write_analytic_input_files

end module test_guiding_centre_step_mod
