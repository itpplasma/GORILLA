module test_tetra_physics_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use funit
    use tetra_physics_mod, only: radial_electric_projection_flux
    implicit none

contains

    @test
    subroutine test_regular_flux_projection()
        real(dp) :: actual

        actual = radial_electric_projection_flux(3.0_dp, 2.0_dp, 0.0_dp, &
                                                  0.5_dp, 0.0_dp)
        @assertEqual(6.0_dp, actual, tolerance=1.0e-14_dp)
    end subroutine test_regular_flux_projection

    @test
    subroutine test_zero_field_at_coordinate_singularity()
        real(dp) :: actual

        actual = radial_electric_projection_flux(0.0_dp, 0.0_dp, 0.0_dp, &
                                                  0.0_dp, 0.0_dp)
        @assertEqual(0.0_dp, actual, tolerance=0.0_dp)
    end subroutine test_zero_field_at_coordinate_singularity

    @test
    subroutine test_axis_has_no_radial_direction()
        real(dp) :: actual

        actual = radial_electric_projection_flux(3.0_dp, 0.0_dp, 0.0_dp, &
                                                  2.0_dp, 1.0_dp)
        @assertEqual(0.0_dp, actual, tolerance=0.0_dp)
    end subroutine test_axis_has_no_radial_direction

end module test_tetra_physics_mod
