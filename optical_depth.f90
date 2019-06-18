module optical_depth
  use iso_c_binding, only: c_double
  use const
  use photoabsorption
  implicit none
  !public :: 
  private
  integer, parameter :: n_BLR = 26
  ! BLR template from Finke 2016, Appendix
  ! wavelength in cm, converted from \AA
  real(dp), dimension(n_BLR), parameter :: &
    lambda_BLR = [937.80_dp, 949.74_dp, 977.02_dp, 990.69_dp, 1025.72_dp, &
      1033.83_dp, 1066.66_dp, 1215.67_dp, 1304.35_dp, 1306.82_dp, 1396.76_dp, &
      1402.06_dp, 1549.06_dp, 1718.55_dp, 1721.89_dp, 1908.73_dp, 2423.83_dp, &
      2798.75_dp, 3188.67_dp, 4102.89_dp, 4341.68_dp, 4687.02_dp, 4862.68_dp, &
      5539.43_dp, 5877.29_dp, 6564.61_dp]*1e-8_dp, &
    ! luminosity in units of L(Hβ)
    lumi_BLR = [0.24_dp, 0.24_dp, 0.60_dp, 0.60_dp, 1.1_dp, 1.1_dp, 0.094_dp, &
      12.0_dp, 0.23_dp, 0.23_dp, 1.0_dp, 1.0_dp, 2.9_dp, 0.030_dp, 0.030_dp, &
      1.8_dp, 0.051_dp, 1.7_dp, 0.051_dp, 0.12_dp, 0.30_dp, 0.016_dp, 1.0_dp, &
      0.039_dp, 0.092_dp, 3.6_dp], &
    ! radius in units of R(Hβ)
    rad_BLR = [2.7_dp, 2.8_dp, 0.83_dp, 0.85_dp, 1.2_dp, 1.2_dp, 4.5_dp, &
      0.27_dp, 4.0_dp, 4.0_dp, 0.83_dp, 0.83_dp, 0.83_dp, 3.8_dp, 3.8_dp, &
      0.46_dp, 5.8_dp, 0.45_dp, 4.3_dp, 3.4_dp, 3.2_dp, 0.63_dp, 1.0_dp, &
      4.8_dp, 0.39_dp, 1.3_dp]
contains
  subroutine optical_depth_blr(tau_gg_blr, eps_c, z_jet_c, redshift_c, &
      L_Hbeta_c, R_Hbeta_c) bind(C)
    implicit none
    real(c_double), intent(out) :: tau_gg_blr
    real(c_double), intent(in) :: eps_c, z_jet_c, redshift_c, L_Hbeta_c, &
      R_Hbeta_c
    real(dp) :: eps, z_jet, redshift, L_Hbeta, R_Hbeta
    real(dp), parameter :: M_8 = 1.0_dp, R_g = 1.5e13_dp*M_8, &
      L_disk = 1.26e46_dp*M_8, l_edd = 1.0_dp

    real(dp), dimension(n_BLR) :: eps_BLR, xi_BLR, R_BLR

    eps = real(eps_c, dp)
    z_jet = real(z_jet_c, dp)
    redshift = real(redshift_c, dp)
    L_Hbeta = real(L_Hbeta_c, dp)
    R_Hbeta = real(R_Hbeta_c, dp)

    eps_BLR = (h*c/lambda_BLR)/(me*c**2)
    xi_BLR = lumi_BLR*L_Hbeta/L_disk
    R_BLR = rad_BLR*R_Hbeta

    tau_gg_blr = real(tau_blr(eps, redshift, l_edd, z_jet, R_g, xi_BLR, R_BLR, &
      eps_BLR, n_BLR), c_double)
  end subroutine optical_depth_blr
end module optical_depth
