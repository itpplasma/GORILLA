program test_radial_electric_projection_driver
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use tetra_physics_mod, only: radial_electric_projection_flux
    implicit none
    real(dp) :: actual

    actual = radial_electric_projection_flux(3.0_dp, 2.0_dp, 0.0_dp, &
                                              0.5_dp, 0.0_dp)
    if (abs(actual - 6.0_dp) > 1.0e-14_dp) error stop "regular projection"

    actual = radial_electric_projection_flux(0.0_dp, 0.0_dp, 0.0_dp, &
                                              0.0_dp, 0.0_dp)
    if (actual /= 0.0_dp) error stop "zero field at coordinate singularity"

    actual = radial_electric_projection_flux(3.0_dp, 0.0_dp, 0.0_dp, &
                                              2.0_dp, 1.0_dp)
    if (actual /= 0.0_dp) error stop "undefined radial direction at axis"

    print *, "radial electric projection oracle: PASS"
end program test_radial_electric_projection_driver
