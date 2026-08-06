!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! guiding_centre_step_mod
!
! Species-agnostic characteristic-step event contract and facade.
!
! This module exposes the existing charged guiding-centre integrator
! (orbit_timestep_gorilla) through a common characteristic-step facade.
! The facade
!
!   * takes a 5-dimensional phase-space state  [x1,x2,x3,vpar,vperp],
!   * advances it by a requested physical time,
!   * reports the elapsed time, the remaining time, the containing cell,
!     the crossed face, a material-surface (domain-boundary) event and the
!     final phase-space state in a plain data result -- without
!     application-global output state.
!
! Locally linear dynamics (for which GORILLA's exact geometric stepping
! is valid) continue to use the GORILLA integrator internally.  The
! equations of motion are supplied through the characteristic_equations_t
! data boundary; only the charged guiding-centre model is implemented on
! this branch, so requesting any other characteristic model fails
! explicitly with CHAR_ERR_UNSUPPORTED.
!
! The event/result contract (CHAR_EVENT_*, CHAR_ERR_*, characteristic_event_t,
! characteristic_result_t) is the shared, species-agnostic vocabulary that
! neutral and fluid-characteristic implementations are meant to reuse.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
module guiding_centre_step_mod
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  private

  ! --- event kinds ---
  integer, parameter, public :: CHAR_EVENT_NONE = 0
  integer, parameter, public :: CHAR_EVENT_FACE = 1
  integer, parameter, public :: CHAR_EVENT_WALL = 2
  integer, parameter, public :: CHAR_EVENT_COLLISION = 3

  ! --- characteristic model identifiers ---
  integer, parameter, public :: CHAR_MODEL_GUIDING_CENTRE = 1
  integer, parameter, public :: CHAR_MODEL_NEUTRAL = 2
  integer, parameter, public :: CHAR_MODEL_FLUID = 3

  ! --- error codes ---
  integer, parameter, public :: CHAR_OK = 0
  integer, parameter, public :: CHAR_ERR_INVALID = 1
  integer, parameter, public :: CHAR_ERR_GRAZING = 2
  integer, parameter, public :: CHAR_ERR_NOT_INSIDE = 3
  integer, parameter, public :: CHAR_ERR_UNSUPPORTED = 4

  ! Phase-space state dimension handled by the guiding-centre facade.
  integer, parameter, public :: CHAR_STATE_DIM = 5

  ! Data boundary supplying the local equations of motion.
  !
  ! For CHAR_MODEL_GUIDING_CENTRE the local EOM are those of the
  ! (already initialized) GORILLA magnetic field, linearized per
  ! tetrahedron, and the facade uses GORILLA's exact geometric stepping
  ! via orbit_timestep_gorilla.  The field is selected through the model
  ! identifier; further per-model data fields can be added later without
  ! changing the call boundary.
  type, public :: characteristic_equations_t
    integer :: model = CHAR_MODEL_GUIDING_CENTRE
  end type characteristic_equations_t

  ! A single characteristic event.
  type, public :: characteristic_event_t
    integer :: kind = CHAR_EVENT_NONE
    integer :: tetrahedron = -1    ! containing cell at the event
    integer :: face = -1           ! crossed / exit face index
    real(dp) :: time = 0.0_dp      ! elapsed time up to the event
    real(dp) :: position(3) = 0.0_dp
    real(dp) :: normal(3) = 0.0_dp ! outward boundary normal (if any)
    real(dp) :: remaining = 0.0_dp ! time still to be integrated
  end type characteristic_event_t

  ! Result of a characteristic step.
  type, public :: characteristic_result_t
    real(dp) :: position(3) = 0.0_dp  ! final phase-space position
    real(dp) :: time = 0.0_dp         ! elapsed time of the step
    type(characteristic_event_t) :: event
    logical :: finished = .false.     ! requested time fully integrated
  end type characteristic_result_t

  public :: trace_characteristic_step

contains

  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  ! trace_characteristic_step
  !
  ! Advance a characteristic by the requested time t_step.
  !
  ! Arguments:
  !   state    [inout] 5-vector [x1,x2,x3,vpar,vperp]; valid dimension is
  !                    CHAR_STATE_DIM.  Returns the final phase-space
  !                    state.
  !   t_step   [in]    requested integration time (>= 0).
  !   equations[in]    local equations of motion data boundary; selects
  !                    the characteristic model.
  !   result   [out]   event/result report (see characteristic_result_t).
  !   ierr     [out]   CHAR_OK on success, or
  !                    CHAR_ERR_INVALID   (invalid state dimension or
  !                                       non-finite state/time),
  !                    CHAR_ERR_UNSUPPORTED (unsupported characteristic
  !                                       model).
  !   boole_initialized [inout,opt] forward compat with the legacy
  !                    orbit_timestep_gorilla caller interface.
  !   ind_tetr [inout,opt] containing cell index (persistent across calls).
  !   iface    [inout,opt] exit/crossed face index (persistent across calls).
  !
  ! The underlying GORILLA mesh and species parameters must have been set
  ! up by a prior call to initialize_gorilla, exactly as for the legacy
  ! charged guiding-centre interface.
  !
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  subroutine trace_characteristic_step(state, t_step, equations, result, ierr, &
                                       boole_initialized, ind_tetr, iface)
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use orbit_timestep_gorilla_mod, only: orbit_timestep_gorilla
    real(dp), intent(inout) :: state(:)
    real(dp), intent(in) :: t_step
    type(characteristic_equations_t), intent(in) :: equations
    type(characteristic_result_t), intent(out) :: result
    integer, intent(out) :: ierr
    logical, intent(inout), optional :: boole_initialized
    integer, intent(inout), optional :: ind_tetr, iface

    real(dp) :: x(3), vpar, vperp, t_remain
    logical :: boole_init_local
    integer :: ind_tetr_local, iface_local

    ! Default: no event, nothing advanced.
    result = characteristic_result_t()
    ierr = CHAR_OK

    ! --- fail explicitly on invalid state dimensions ---
    if (size(state) /= CHAR_STATE_DIM) then
      ierr = CHAR_ERR_INVALID
      return
    end if

    ! --- fail explicitly on non-finite fields ---
    if (.not. (ieee_is_finite(t_step) .and. &
               ieee_is_finite(state(1)) .and. ieee_is_finite(state(2)) .and. &
               ieee_is_finite(state(3)) .and. ieee_is_finite(state(4)) .and. &
               ieee_is_finite(state(5)))) then
      ierr = CHAR_ERR_INVALID
      return
    end if

    ! --- fail explicitly on unsupported characteristic models ---
    if (equations%model /= CHAR_MODEL_GUIDING_CENTRE) then
      ierr = CHAR_ERR_UNSUPPORTED
      return
    end if

    ! Decode the 5-dimensional phase-space state.
    x = state(1:3)
    vpar = state(4)
    vperp = state(5)

    ! Persistent legacy interface state, defaulting to a fresh particle.
    if (present(boole_initialized)) then
      boole_init_local = boole_initialized
    else
      boole_init_local = .false.
    end if
    if (present(ind_tetr)) then
      ind_tetr_local = ind_tetr
    else
      ind_tetr_local = -1
    end if
    if (present(iface)) then
      iface_local = iface
    else
      iface_local = -1
    end if

    ! Advance the charged guiding centre with GORILLA's exact geometric
    ! stepping (locally linear dynamics).
    call orbit_timestep_gorilla(x, vpar, vperp, t_step, boole_init_local, &
                                ind_tetr_local, iface_local, t_remain)

    ! Encode the final phase-space state back into the caller's buffer.
    state(1:3) = x
    state(4) = vpar
    state(5) = vperp

    ! Persist the legacy bookkeeping, if requested.
    if (present(boole_initialized)) boole_initialized = boole_init_local
    if (present(ind_tetr)) ind_tetr = ind_tetr_local
    if (present(iface)) iface = iface_local

    ! --- assemble the species-agnostic result ---
    result%position = x
    result%time = t_step - t_remain
    result%event%position = x
    result%event%time = result%time
    result%event%remaining = t_remain
    result%event%tetrahedron = ind_tetr_local
    result%event%face = iface_local

    if (ind_tetr_local == -1) then
      ! Particle reached / lost at the domain (material) boundary.
      result%event%kind = CHAR_EVENT_WALL
      result%finished = .false.
    else
      ! Requested time integrated within the containing cell.
      result%event%kind = CHAR_EVENT_NONE
      result%finished = .true.
    end if

  end subroutine trace_characteristic_step

end module guiding_centre_step_mod
