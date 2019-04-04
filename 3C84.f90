program test_absorption
  ! Reproduce figure 14 from Finke 2016 to test absorption routines
  use photoabsorption
  use const
  implicit none
  ! Parameters from Table 3
  real(dp), parameter :: M_8 = 2.9e0_dp, R_g = 1.5e13_dp*M_8, L_disk = 8.97e43_dp, &
    eta = 1/12.0_dp, Theta = k*1e3_dp/(me*c**2), xi_dt = 0.1_dp, &
    R_dt = 6.5e17_dp, redshift = 0.0_dp
  ! Derived parameter
  real(dp), parameter :: l_edd = L_disk/(1.26e46_dp*M_8)
  ! From Appendix
  real(dp), parameter :: L_H_beta = 9.28e40_dp, R_H_beta = 2.41e16_dp

  integer :: i, j
  integer, parameter :: steps = 10, n_BLR = 26
  real(dp) :: eps, energy, z_jet
  real(dp), dimension(n_BLR) :: eps_BLR, xi_BLR, lambda_BLR, lumi_BLR, &
    rad_BLR, R_BLR

  ! BLR template from Finke 2016, Appendix
  ! wavelength in cm, converted from \AA
  lambda_BLR = [937.80_dp, 949.74_dp, 977.02_dp, 990.69_dp, 1025.72_dp, &
    1033.83_dp, 1066.66_dp, 1215.67_dp, 1304.35_dp, 1306.82_dp, 1396.76_dp, &
    1402.06_dp, 1549.06_dp, 1718.55_dp, 1721.89_dp, 1908.73_dp, 2423.83_dp, &
    2798.75_dp, 3188.67_dp, 4102.89_dp, 4341.68_dp, 4687.02_dp, 4862.68_dp, &
    5539.43_dp, 5877.29_dp, 6564.61_dp]*1e-8_dp
  eps_BLR = (h*c/lambda_BLR)/(me*c**2)
  ! luminosity in units of L(Hβ)
  lumi_BLR = [0.24_dp, 0.24_dp, 0.60_dp, 0.60_dp, 1.1_dp, 1.1_dp, 0.094_dp, &
    12.0_dp, 0.23_dp, 0.23_dp, 1.0_dp, 1.0_dp, 2.9_dp, 0.030_dp, 0.030_dp, &
    1.8_dp, 0.051_dp, 1.7_dp, 0.051_dp, 0.12_dp, 0.30_dp, 0.016_dp, 1.0_dp, &
    0.039_dp, 0.092_dp, 3.6_dp]
  ! radius in units of R(Hβ)
  rad_BLR = [2.7_dp, 2.8_dp, 0.83_dp, 0.85_dp, 1.2_dp, 1.2_dp, 4.5_dp, &
    0.27_dp, 4.0_dp, 4.0_dp, 0.83_dp, 0.83_dp, 0.83_dp, 3.8_dp, 3.8_dp, &
    0.46_dp, 5.8_dp, 0.45_dp, 4.3_dp, 3.4_dp, 3.2_dp, 0.63_dp, 1.0_dp, 4.8_dp, &
    0.39_dp, 1.3_dp]

  xi_BLR = lumi_BLR*L_H_beta/L_disk
  R_BLR = rad_BLR*R_H_beta


  do i = -2*steps, 2*steps
    z_jet = R_BLR(8)*10.0_dp**(i/real(steps, dp)) ! R_BLR(8) = R(Ly α)
    do j = 0*steps, 5*steps
      energy = 1e-3_dp*TeV*10.0_dp**(j/real(steps, dp))
      eps = energy/(me*c**2)
      print *, z_jet/R_BLR(8), energy*1e3_dp/TeV, &
        tau_disk(eps, redshift, l_edd, M_8, R_g, eta, z_jet), &
        tau_blr(eps, redshift, l_edd, z_jet, R_g, xi_BLR, R_BLR, eps_BLR, &
          n_BLR), &
        tau_dust(eps, redshift, xi_dt, l_edd, Theta, R_dt, R_g, z_jet)
    end do
    print *,
  end do
end program
