module test_vmec_asymmetry
  use funit
  use spline_vmec_data_mod, only: add_vmec_fourier_mode,vmec_mode_is_selected
  implicit none

contains

  @test
  subroutine test_asymmetric_fourier_mode_matches_direct_oracle()
    double precision, parameter :: theta=0.37d0,varphi=0.23d0
    double precision :: phase,cosphase,sinphase,R,Z,alam
    double precision :: expected_R,expected_Z,expected_alam

    phase=2.d0*theta-3.d0*varphi
    cosphase=cos(phase)
    sinphase=sin(phase)
    R=1.25d0
    Z=-0.5d0
    alam=0.75d0

    call add_vmec_fourier_mode(cosphase,sinphase,2.d0,-3.d0, &
                               5.d0,7.d0,-11.d0,13.d0,R,Z,alam)

    expected_R=1.25d0+2.d0*cos(phase)-3.d0*sin(phase)
    expected_Z=-0.5d0+5.d0*cos(phase)+7.d0*sin(phase)
    expected_alam=0.75d0-11.d0*cos(phase)+13.d0*sin(phase)
    @assertEqual(expected_R,R,tolerance=1.d-14)
    @assertEqual(expected_Z,Z,tolerance=1.d-14)
    @assertEqual(expected_alam,alam,tolerance=1.d-14)
  end subroutine test_asymmetric_fourier_mode_matches_direct_oracle

  @test
  subroutine test_symmetric_mode_matches_legacy_series()
    double precision, parameter :: phase=-0.83d0
    double precision :: R,Z,alam

    R=0.d0
    Z=0.d0
    alam=0.d0
    call add_vmec_fourier_mode(cos(phase),sin(phase),2.d0,0.d0, &
                               0.d0,7.d0,0.d0,13.d0,R,Z,alam)

    @assertEqual(2.d0*cos(phase),R,tolerance=1.d-14)
    @assertEqual(7.d0*sin(phase),Z,tolerance=1.d-14)
    @assertEqual(13.d0*sin(phase),alam,tolerance=1.d-14)
  end subroutine test_symmetric_mode_matches_legacy_series

  @test
  subroutine test_axisymmetric_filter_retains_only_n_zero()
    @assertTrue(vmec_mode_is_selected(0,.true.))
    @assertFalse(vmec_mode_is_selected(2,.true.))
    @assertFalse(vmec_mode_is_selected(-2,.true.))
    @assertTrue(vmec_mode_is_selected(2,.false.))
    @assertTrue(vmec_mode_is_selected(-2,.false.))
  end subroutine test_axisymmetric_filter_retains_only_n_zero

end module test_vmec_asymmetry
