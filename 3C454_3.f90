program fsrq
! SED of the FSRQ 3C 454.3 as modelled in arXiv:1811.07654
  use ssc
  use external_compton
  use photoabsorption
  use luminosity_distance
  implicit none
  real(dp), parameter :: redshift = 0.859_dp,  g_max = 5e3_dp, M_8 = 1.2e1_dp
  real(dp), parameter :: z_jet = 1e18_dp, R_blob = 4.13e16_dp, &
    L_disk = 2e46_dp, f_BLR = 3e-2_dp, R_Torus = 1e19_dp, &!, R_BLR = 1.6e18_dp
    f_Torus = 0.1_dp, delta = 39.3_dp, Gamm = delta/2
  real(dp), parameter :: B = 0.56_dp, g_pk = 180.0_dp, Ne_pk = 2.4e50_dp, &
    bb = 1.0_dp, g_min = 1.0_dp 

  real(dp) :: l_edd = L_disk/(1.26e46_dp*M_8), eta = 1/12.0_dp, &
    Theta = k*1e3_dp/(me*c**2)

  integer :: j
  integer, parameter :: steps = 10, n_BLR = 26
  real(dp) :: nu, eps, dL, tau, tau_internal, tau_external, R_g, absorption
  real(dp), dimension(n_BLR) :: eps_BLR, xi_BLR, lambda_BLR, lumi_BLR, &
    rad_BLR, R_BLR

  dL = lum_dist(redshift)
  R_g = 1.8e14

  ! BLR template from Finke 2016, Appendix
  ! wavelength in cm, converted from \AA
  lambda_BLR = [937.80_dp, 949.74_dp, 977.02_dp, 990.69_dp, 1025.72_dp, &
    1033.83_dp, 1066.66_dp, 1215.67_dp, 1304.35_dp, 1306.82_dp, 1396.76_dp, &
    1402.06_dp, 1549.06_dp, 1718.55_dp, 1721.89_dp, 1908.73_dp, 2423.83_dp, &
    2798.75_dp, 3188.67_dp, 4102.89_dp, 4341.68_dp, 4687.02_dp, 4862.68_dp, &
    5539.43_dp, 5877.29_dp, 6564.61_dp]*1e-10_dp
  eps_BLR = (h*c/lambda_BLR)/(me*c**2)
  ! luminosity in units of L(Hβ)
  lumi_BLR = [0.24_dp, 0.24_dp, 0.60_dp, 0.60_dp, 1.1_dp, 1.1_dp, 0.094_dp, &
    12.0_dp, 0.23_dp, 0.23_dp, 1.0_dp, 1.0_dp, 2.9_dp, 0.030_dp, 0.030_dp, &
    1.8_dp, 0.051_dp, 1.7_dp, 0.051_dp, 0.12_dp, 0.30_dp, 0.016_dp, 1.0_dp, &
    0.039_dp, 0.092_dp, 3.6_dp]
  xi_BLR = lumi_BLR*4.2e43_dp/L_disk
  ! radius in units of R(Hβ)
  rad_BLR = [2.7_dp, 2.8_dp, 0.83_dp, 0.85_dp, 1.2_dp, 1.2_dp, 4.5_dp, &
    0.27_dp, 4.0_dp, 4.0_dp, 0.83_dp, 0.83_dp, 0.83_dp, 3.8_dp, 3.8_dp, &
    0.46_dp, 5.8_dp, 0.45_dp, 4.3_dp, 3.4_dp, 3.2_dp, 0.63_dp, 1.0_dp, 4.8_dp, &
    0.39_dp, 1.3_dp]
  R_BLR = rad_BLR*4.3e17_dp
  !do j = 9*steps, 25*steps
  do j = 9*steps, 35*steps
    nu = 10**(real(j, dp)/steps)
    eps = h*nu/(me*c**2)
    tau_internal = tau_ssc(eps, delta, redshift, dL, Ne, B, R_blob, g_min, &
      g_max)
    tau_external = tau_disk(eps, redshift, l_edd, M_8, eta, z_jet) &
      + tau_blr(eps, redshift, l_edd, z_jet, R_g, xi_BLR, [R_BLR], eps_BLR, &
        n_BLR) &
      + tau_dust(eps, redshift, f_Torus, l_edd, Theta, R_Torus, R_g, z_jet)
    tau = tau_internal + tau_external + epsrel
    absorption = (1 - exp(-tau))/tau
    print *, nu, &
      f_syn(eps, delta, redshift, dL, Ne, B, g_min, g_max)*absorption, &
      f_disk(eps, redshift, dL, M_8, l_edd, eta)*absorption, &
      f_ssc(eps, delta, redshift, dL, Ne, B, R_blob, g_min, g_max)*absorption, &
      ec_disk(eps, redshift, dL, delta, Gamm, Ne, M_8, l_edd, eta, z_jet, &
        g_min)*absorption, &
      ec_blr(eps, redshift, dL, delta, Gamm, Ne, L_disk, z_jet, [f_BLR], &
        [R_BLR], [eps_BLR], n_BLR, g_min)*absorption, &
      f_dust(eps, redshift, dL, Theta, L_disk, f_Torus)*absorption, &
      ec_dust(eps, redshift, dL, delta, Gamm, Ne, L_disk, z_jet, f_Torus, &
        R_Torus, Theta, g_min)*absorption
    !print *, nu, &
    !  tau_disk(eps, redshift, l_edd, M_8, eta, z_jet), &
    !  tau_blr(eps, redshift, l_edd, z_jet, R_g, [f_BLR], [R_BLR], [eps_BLR], &
    !    n_BLR), &
    !  tau_dust(eps, redshift, f_Torus, l_edd, Theta, R_Torus, R_g, z_jet)
  end do
contains
  function Ne(g)
    implicit none
    real(dp) :: Ne, g

    if (g >= g_min .and. g <= g_max) then
      Ne = Ne_pk*(g/g_pk)**(-(2 + bb*log10(g/g_pk)))
    else
      Ne = 0.0_dp
    end if
  end function Ne
end program
