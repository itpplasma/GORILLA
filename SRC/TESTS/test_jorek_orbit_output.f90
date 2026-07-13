program test_jorek_orbit_output
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite

    implicit none

    real(dp) :: time, value, first_value, last_time, minimum, maximum
    real(dp) :: r, phi, z, first_pphi, last_pphi
    integer :: iunit, status, count

    open (newunit=iunit, file='e_tot.dat', status='old', action='read')
    count = 0
    do
        read (iunit, *, iostat=status) time, value
        if (status /= 0) exit
        if (.not. ieee_is_finite(time) .or. .not. ieee_is_finite(value)) &
            error stop 'non-finite orbit energy record'
        count = count + 1
        if (count == 1) then
            first_value = value
            minimum = value
            maximum = value
        end if
        last_time = time
        minimum = min(minimum, value)
        maximum = max(maximum, value)
    end do
    close (iunit)
    if (count /= 18) error stop 'unexpected number of completed orbit steps'
    if (abs(last_time - 1.0e-5_dp) > 1.0e-14_dp) &
        error stop 'orbit did not reach the requested final time'
    if ((maximum - minimum)/abs(first_value) > 1.0e-12_dp) &
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
        count = count + 1
    end do
    close (iunit)
    if (count /= 18) error stop 'orbit position record is incomplete'

    open (newunit=iunit, file='p_phi.dat', status='old', action='read')
    count = 0
    do
        read (iunit, *, iostat=status) time, value
        if (status /= 0) exit
        count = count + 1
        if (count == 1) first_pphi = value
        last_pphi = value
    end do
    close (iunit)
    if (count /= 18) error stop 'canonical-momentum record is incomplete'
    if (.not. ieee_is_finite(first_pphi) .or. &
            .not. ieee_is_finite(last_pphi)) &
        error stop 'non-finite canonical momentum'
    if (abs(first_pphi) <= tiny(first_pphi)) &
        error stop 'zero initial canonical momentum'
    if (abs(last_pphi/first_pphi - 1.0_dp) < 1.0e-6_dp) &
        error stop 'non-axisymmetric snapshot response was not resolved'

    print '(A)', 'PASS: JOREK snapshot orbit completion and invariants'
end program test_jorek_orbit_output
