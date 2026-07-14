program test_jorek_orbit_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    real(dp), parameter :: final_time = 1.0e-5_dp
    real(dp), parameter :: time_tolerance = 1.0e-14_dp
    real(dp) :: time, value, first_value, last_time, previous_time
    real(dp) :: minimum, maximum, relative_energy_range
    real(dp) :: r, phi, z, first_pphi, last_pphi, last_pphi_time
    real(dp) :: minimum_pphi_response, relative_pphi_response
    character(len=64) :: argument
    integer :: iunit, status, count, energy_count

    minimum_pphi_response = 1.0e-6_dp
    if (command_argument_count() >= 1) then
        call get_command_argument(1, argument)
        read (argument, *) minimum_pphi_response
    end if

    open (newunit=iunit, file='e_tot.dat', status='old', action='read')
    count = 0
    previous_time = -huge(previous_time)
    do
        read (iunit, *, iostat=status) time, value
        if (status /= 0) exit
        if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(value)) &
            error stop 'non-finite orbit energy record'
        if (time <= previous_time) error stop 'orbit energy times are not increasing'
        count = count + 1
        if (count == 1) then
            first_value = value
            minimum = value
            maximum = value
        end if
        last_time = time
        previous_time = time
        minimum = min(minimum, value)
        maximum = max(maximum, value)
    end do
    close (iunit)
    if (count < 2) error stop 'orbit energy history is incomplete'
    energy_count = count
    if (abs(last_time - final_time) > time_tolerance) &
        error stop 'orbit did not reach the requested final time'
    relative_energy_range = (maximum - minimum)/abs(first_value)
    if (relative_energy_range > 1.0e-12_dp) &
        error stop 'orbit energy conservation regression'

    open (newunit=iunit, file='full_orbit_plot_rphiz.dat', &
        status='old', action='read')
    count = 0
    do
        read (iunit, *, iostat=status) r, phi, z
        if (status /= 0) exit
        if (.not. all(ieee_is_finite([r, phi, z]))) &
            error stop 'non-finite orbit position'
        if (r < 156.0_dp .or. r > 186.0_dp &
                .or. z < -15.0_dp .or. z > 15.0_dp) &
            error stop 'orbit left the core regression mesh'
        if (phi < 0.0_dp .or. phi >= 2.0_dp*acos(-1.0_dp)) &
            error stop 'orbit toroidal angle is outside its periodic interval'
        count = count + 1
    end do
    close (iunit)
    if (count /= energy_count) error stop 'orbit position history is incomplete'

    open (newunit=iunit, file='p_phi.dat', status='old', action='read')
    count = 0
    previous_time = -huge(previous_time)
    do
        read (iunit, *, iostat=status) time, value
        if (status /= 0) exit
        if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(value)) &
            error stop 'non-finite canonical momentum record'
        if (time <= previous_time) &
            error stop 'canonical-momentum times are not increasing'
        count = count + 1
        if (count == 1) first_pphi = value
        last_pphi = value
        last_pphi_time = time
        previous_time = time
    end do
    close (iunit)
    if (count /= energy_count) error stop 'canonical-momentum history is incomplete'
    if (abs(last_pphi_time - final_time) > time_tolerance) &
        error stop 'canonical-momentum history did not reach the final time'
    if (abs(first_pphi) <= tiny(first_pphi)) &
        error stop 'zero initial canonical momentum'
    relative_pphi_response = abs(last_pphi/first_pphi - 1.0_dp)
    if (relative_pphi_response < minimum_pphi_response) &
        error stop 'non-axisymmetric snapshot response was not resolved'

    print '(A, I0)', 'records=', energy_count
    print '(A, ES12.4)', 'final time=', last_time
    print '(A, ES12.4)', 'relative energy range=', relative_energy_range
    print '(A, ES12.4)', 'relative canonical-momentum response=', &
        relative_pphi_response
    print '(A)', 'PASS: JOREK snapshot orbit completion and invariants'
end program test_jorek_orbit_output
