program test_external_compton
  ! Reproduce fig. 8
  use external_compton
  implicit none
  ! Parameters of baseline FSRQ from Table 1
  real(dp), parameter :: Gamm = 40.0_dp, delta_D = 40.0_dp, B = 0.56_dp, &
    g1 = 20.0_dp, g2 = 5e7_dp, gb = 1e4_dp, p1 = 2.0_dp, p2 = 3.5_dp, &
    M_8 = 12.0_dp, R_g = 1.8e14_dp, L_disk = 2e46_dp, eta = 1/12.0_dp
  real(dp), parameter :: d_L = 1e28_dp, z = 0.0_dp, N0 = 1e50_dp, &
    Lumi_edd = 1.26e46_dp*M_8, l_edd = L_disk/Lumi_edd
  real(dp), parameter :: T_dt_3 = 1.0_dp, Theta = T_dt_3*k*1e3_dp/(me*c**2), xi_dt = 0.1_dp, &
    R_dt = 8e5*sqrt(l_edd/M_8)*T_dt_3**(-2.6_dp)*R_g
  real(dp), dimension(4), parameter :: rr = [1e18_dp, 1e19_dp, 1e20_dp, 1e21_dp]
  real(dp) :: nu, eps_s
  integer :: j, steps = 10

  do j = 16*steps, 29*steps
    nu = 10**(j/real(steps, dp))
    eps_s = h*nu/(me*c**2)
    print *, nu, &
      ec_dust(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr(1), xi_dt, R_dt, &
        Theta, g1), &
      ec_dust(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr(2), xi_dt, R_dt, &
        Theta, g1), &
      ec_dust(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr(3), xi_dt, R_dt, &
        Theta, g1), &
      ec_dust(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr(4), xi_dt, R_dt, &
        Theta, g1)
  end do
contains
  function Ne(g)
    implicit none
    real(dp) :: Ne, g

    if (g < g1 .or. g > g2) then
      Ne = 0.0_dp
    else if (g <= gb) then
      Ne = N0*(g/gb)**(-p1)
    else
      Ne = N0*(g/gb)**(-p2)
    end if
  end function Ne
end program test_external_compton
