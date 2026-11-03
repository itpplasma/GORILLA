!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
module gorilla_status_mod
!
    ! Terminal and non-terminal state of one marker (particle) after a call to
    ! pusher_tetra_poly, orbit_timestep_gorilla or orbit_timestep_gorilla_batch.
    !
    ! The values are part of the C ABI frozen in SRC/gorilla_cuda.h and must not
    ! be renumbered without updating that header.
    !
    implicit none
!
    private
!
    ! Particle is still inside the mesh and may be advanced further. The cell in
    ! which it currently resides and the entry face are returned in ind_tetr/iface.
    integer, parameter, public :: GORILLA_ACTIVE           = 0
!
    ! The requested orbital time step has been integrated completely, the particle
    ! stopped inside the mesh (ind_tetr is the cell of the final position, iface=0).
    integer, parameter, public :: GORILLA_COMPLETE         = 1
!
    ! The particle left the computation domain through a physical boundary
    ! (ind_tetr < 0). This is a legitimate physical outcome, not a failure of the
    ! integrator. The last valid cell and the last relative position are retained
    ! so that vperp can still be reconstructed.
    integer, parameter, public :: GORILLA_PHYSICAL_LOSS    = 2
!
    ! The device could not decide the cell traversal and hands the untouched entry
    ! state back to the host for exact CPU replay of that single cell. Only used by
    ! the batch/device backends, never by the scalar CPU pusher.
    integer, parameter, public :: GORILLA_REPLAY_CELL      = 3
!
    ! The integrator failed (no analytical root, failed consistency check, failed
    ! troubleshooting, non-finite value, host replay failure). Position, velocity and
    ! remaining time are those of the last accepted state; no elapsed time of the
    ! failed cell may be subtracted from the remaining time.
    integer, parameter, public :: GORILLA_NUMERICAL_FAILURE = 4
!
    integer, parameter, public :: GORILLA_STATUS_MAX = GORILLA_NUMERICAL_FAILURE
!
    public :: gorilla_status_name
!
    contains
!
        pure function gorilla_status_name(status) result(name)
!
            implicit none
!
            integer, intent(in) :: status
            character(len=20)   :: name
!
            select case(status)
                case(GORILLA_ACTIVE)
                    name = 'ACTIVE'
                case(GORILLA_COMPLETE)
                    name = 'COMPLETE'
                case(GORILLA_PHYSICAL_LOSS)
                    name = 'PHYSICAL_LOSS'
                case(GORILLA_REPLAY_CELL)
                    name = 'REPLAY_CELL'
                case(GORILLA_NUMERICAL_FAILURE)
                    name = 'NUMERICAL_FAILURE'
                case default
                    name = 'UNKNOWN'
            end select
!
        end function gorilla_status_name
!
end module gorilla_status_mod
