program dermer
! Reproduce SED of 3C279 from Dermer et al. (2014), ApJ 782, 82.
  use ssc
  use external_compton
  use photoabsorption
  use luminosity_distance
  implicit none
!  real(dp), parameter :: z = 0.5362_dp, Theta = k*560/(me*c**2), &
!    xi_dt = 0.2_dp, Gamm = 25.7_dp, delta_D = Gamm, B = 2.6_dp, &
!    g_pk = 31.6_dp, L_disk = 0.55e46_dp, v_Lv_syn_pk = 0.023e48_dp, &
!    t_var = 1e4_dp, R = c*delta_D*t_var/(1 + z)
!  real(dp), parameter :: eps_IR = 2e-7_dp, u_IR = 0.77e-3_dp, &
!    eps_BLR = 2e-5_dp, u_BLR = 2.16e-3_dp
  real(dp), parameter :: Ke = 1e49_dp, g1 = 1e2_dp, gc = 1e3_dp, z = 0.116_dp, &
    delta_D = 100.0_dp, B = 10e-3_dp, t_var = 300.0_dp, p = 1.8_dp, &
    R = c*delta_D*t_var/(1 + z)

  integer :: j, steps = 10
  real(dp) :: nu, eps, d_L!, g_min = 1.0_dp, g_max = 1e8_dp
  real(dp) :: g_min = g1, g_max = 100*gc

  d_L = lum_dist(z)
  do j = 10*steps, 24*steps
    nu = 10.0_dp**(real(j, dp)/steps)
    eps = h*nu/(me*c**2)
    print *, nu, &
      f_syn(eps, delta_D, z, d_L, Ne, B, g_min, g_max), &
      f_ssc(eps, delta_D, z, d_L, Ne, B, R, g_min, g_max)!, &
      !ec_iso_mono(eps, z, d_L, delta_D, Ne, u_IR, eps_IR, g_min), & 
      !ec_iso_mono(eps, z, d_L, delta_D, Ne, u_BLR, eps_BLR, g_min)
  end do
contains
  function Ne(g)
    implicit none
    real(dp) :: Ne, g!, x
    !real(dp), parameter :: bb = 0.8_dp, Ne_pk = 12*pi*v_Lv_syn_pk/(c*sigma_T &
    !  *B**2*delta_D**4*g_pk**2)
      !Ne_pk = v_Lv_syn_pk/3/(4/3.0_dp*c*sigma_t*B**2/(8*pi)*g_pk**2*delta_D**4)

    !x = g/g_pk
    !Ne = Ne_pk*x**(-2 - bb*log(x))
    if (g > g1) then
      Ne = Ke*g**(-p)*exp(-g/gc)
    else
      Ne = 0.0_dp
    end if
  end function Ne
end program
