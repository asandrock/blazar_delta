program W_Com
! Reproduce SED of W Comae (cf. Acciari et al.)
  use ssc
  use external_compton
  use photoabsorption
  use luminosity_distance
  implicit none
  real(dp), parameter :: delta_D = 20.0_dp, Gamm = delta_D, L_e = 5.7e43_dp, &
    R = 1e16_dp, g1 = 8e3_dp, g2 = 3e5_dp, q = 2.55_dp, B = 0.35_dp, &
    z = 0.102_dp, nu_ext = 1.5e14_dp, u_ext = 2.4e-4_dp, g_c = 39.5_dp, &
    eps_ext = h*nu_ext/(me*c**2)
  real(dp), parameter :: g_min = g_c, g_break = g1, g_max = g2, &
    N0 = L_e/(pi*R**2*Gamm**2*sqrt(1 - 1/Gamm**2)*me*c**3 &
      *(g_break**2*log(g_break/g_min) &
      + (g_break**2 - g_break**(q + 1)*g_max**(1 - q))/(q - 1))) &
      *4*pi*R**3/3

  integer :: j, steps = 10
  real(dp) :: nu, eps, d_L, tau

  d_L = lum_dist(z)
  do j = 9*steps, 28*steps
    nu = 10.0_dp**(real(j, dp)/steps)
    eps = h*nu/(me*c**2)
    tau = tau_ssc(eps, delta_D, z, d_L, Ne, B, R, g_min, g_max)
    print *, nu, &
      f_syn(eps, delta_D, z, d_L, Ne, B, g_min, g_max), &
      f_ssc(eps, delta_D, z, d_L, Ne, B, R, g_min, g_max)*absorption(tau), &
      ec_iso_mono(eps, z, d_L, delta_D, Ne, u_ext, eps_ext, g_min) &
        *absorption(tau)
  end do
contains
  function Ne(g)
    implicit none
    real(dp) :: Ne, g

    if (g < g_min .or. g > g_max) then
      Ne = 0.0_dp
    else if (g < g_break) then
      Ne = N0*(g/g_break)**(-2)
    else
      Ne = N0*(g/g_break)**(- q - 1)
    end if
  end function Ne

  function absorption(tau_gg)
    implicit none
    real(dp) :: absorption, tau_gg

    if (tau_gg > 1e-3_dp) then
      absorption = (1 - exp(-tau_gg))/tau_gg
    else
      absorption = 1.0_dp
    end if
  end function absorption
end program W_Com
