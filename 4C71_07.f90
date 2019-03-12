program fsrq
! SED of the FSRQ 4C +71.07 as modelled in arXiv:1811.07654
  use const
  use ssc
  use external_compton
  use photoabsorption
  use luminosity_distance
  implicit none
  real(dp), parameter :: redshift = 2.172_dp, M_8 = 50.0_dp
  real(dp), parameter :: z_jet = 1e18_dp, R_blob = 4e16_dp, &
    L_disk = 2.2e47_dp, R_BLR = 1.6e18_dp, f_BLR = 3e-2_dp, R_Torus = 1e19_dp, &
    f_Torus = 50e-2_dp, Gamm = 20.0_dp, delta = 27.0_dp, B = 0.5_dp!0.9_dp
  real(dp), parameter :: a1 = 2.1_dp, a2 = 3.5_dp, &!5.0_dp, &
    g_min = 20.0_dp, gb = 750.0_dp, g_max = 5e3_dp, &!5e3_dp, &
    Ke = 19*(4*pi*R_blob**3)/3

  real(dp) :: l_edd = L_disk/(1.26e46_dp*M_8), eta = 0.3_dp, &
    Theta = k*100/(me*c**2)

  integer :: j
  integer, parameter :: steps = 10, n_BLR = 1
  real(dp) :: nu, eps, dL, tau, tau_internal, tau_external, R_g, absorption
  real(dp), dimension(n_BLR) :: eps_BLR

  dL = lum_dist(redshift)
  R_g = 1.5e13_dp*M_8
  eps_BLR = [2e-5_dp]

  do j = 9*steps, 25*steps
  !do j = 9*steps, 35*steps
    nu = 10**(real(j, dp)/steps)
    eps = h*nu/(me*c**2)
    tau_internal = tau_ssc(eps, delta, redshift, dL, Ne, B, R_blob, g_min, &
      g_max)
    tau_external = tau_disk(eps, redshift, l_edd, M_8, eta, z_jet) &
      + tau_blr(eps, redshift, l_edd, z_jet, R_g, [f_BLR], [R_BLR], [eps_BLR], &
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
      Ne = Ke/gb/((g/gb)**a1 + (g/gb)**a2)
    else
      Ne = 0.0_dp
    end if
  end function Ne
end program
