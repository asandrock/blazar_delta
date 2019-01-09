program test_external_compton
  ! Reproduce fig. 8
  use external_compton
  implicit none
  ! Parameters of baseline FSRQ from Table 1
  real(dp), parameter :: Gamm = 40.0_dp, delta_D = 40.0_dp, B = 0.56_dp, &
    g1 = 20.0_dp, g2 = 5e7_dp, gb = 1e4_dp, p1 = 2.0_dp, p2 = 3.5_dp, &
    M_8 = 12.0_dp, R_g = 1.8e14_dp, L_disk = 2e46_dp, eta = 1/12.0_dp
  real(dp), parameter :: rr = 1e17_dp, d_L = 1e28_dp, z = 0.0_dp, &
    u = L_disk/(4*pi*rr**2*c), Lumi_edd = 1.26e46_dp*M_8, &
    l_edd = L_disk/Lumi_edd, N0 = 1e50_dp

  real(dp) :: nu, eps_s, e0
  integer :: j, steps = 10

  e0 = eps0(8.17_dp*R_g)
  do j = 18*steps, 29*steps
    nu = 10**(j/real(steps, dp))
    eps_s = h*nu/(me*c**2)
    print *, nu, &
      ec_disk(eps_s, z, d_L, delta_D, Gamm, Ne, M_8, l_edd, eta, rr), &
      ec_point_mono(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, e0, rr), &
      ec_iso_mono(eps_s, z, d_L, delta_D, Ne, u, e0)
  end do
contains
  function Ne(g)
    implicit none
    real(dp) :: Ne, g

    !if (g < g1 .or. g > g2) then
    !  Ne = 0.0_dp
    !else
 if (g <= gb) then
      Ne = N0*(g/gb)**(-p1)
    else
      Ne = N0*(g/gb)**(-p2)
    end if
  end function Ne

  function eps0(R)
    ! cf. eq. (64)
    implicit none
    real(dp) :: eps0, R

    eps0 = 2.7e-4_dp*(l_edd/(M_8*eta))**0.25_dp*(R/R_g)**(-0.75_dp)
  end function eps0
end program test_external_compton
