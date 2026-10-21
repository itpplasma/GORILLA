program test_marker_species
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gorilla_settings_mod, only: load_gorilla_inp
    use orbit_timestep_gorilla_mod, only: apply_marker_species
    use tetra_physics_mod, only: particle_mass, particle_charge, cm_over_e
    use constants, only: clight, echarge
    implicit none
    character(len=80) :: arg
    real(dp) :: expected_mass, expected_charge

    call get_command_argument(1, arg)
    read(arg, *) expected_mass
    call get_command_argument(2, arg)
    read(arg, *) expected_charge
    particle_mass = 1.0_dp
    particle_charge = 2.0_dp
    cm_over_e = 3.0_dp
    call load_gorilla_inp()
    call apply_marker_species()
    if (expected_mass <= 0.0_dp) then
        if (abs(particle_mass - 1.0_dp) > 0.0_dp) error stop 'Default mass changed'
        if (abs(particle_charge - 2.0_dp) > 0.0_dp) error stop 'Default charge changed'
        if (abs(cm_over_e - 3.0_dp) > 0.0_dp) error stop 'Default cm_over_e changed'
    else
        if (abs(particle_mass / expected_mass - 1.0_dp) > 1e-15_dp) &
            error stop 'Marker mass differs from native NEO-2 species'
        if (abs(particle_charge / (expected_charge * echarge) - 1.0_dp) > 1e-15_dp) &
            error stop 'Marker charge differs from native NEO-2 species'
        if (abs(cm_over_e * particle_charge / (clight * particle_mass) - 1.0_dp) > 1e-15_dp) &
            error stop 'Inconsistent pusher mass/charge ratio'
    end if
    print *, 'PASS: marker species and pusher ratio'
end program test_marker_species
