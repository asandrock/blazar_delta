module photoabsorption
  ! Photoabsorption by external radiation fields
  use const
  implicit none
contains
  function tau_disk(eps_1, z, l_edd, M_8, eta, rr)
    implicit none
    real(dp) :: tau_disk, eps_1, l_edd, M_8, eta, rr, z
    real(dp) :: l_err, l_save, prefactor, R_g, R_in, R_out
    integer :: l_neval, l_ier

    R_g = 1.5e13_dp*M_8
    R_in = 6*R_g ! Schwarzschild black hole
    R_out = 200*R_g ! simple assumption from Finke (2016)
    prefactor = 1e7_dp*l_edd**0.75_dp*M_8**0.25_dp/eta**0.75_dp
    call qagi(l_int, rr/R_g, 1, epsabs, epsrel, tau_disk, l_err, l_neval, l_ier)
  contains
    function l_int(l)
      implicit none
      real(dp) :: l_int, l
      real(dp) :: R_err
      integer :: R_neval, R_ier

      l_save = l
      call qng(R_int, R_in/R_g, R_out/R_g, epsabs, epsrel, l_int, R_err, &
        R_neval, R_ier)
    end function l_int

    function R_int(Rt)
      implicit none
      real(dp) :: R_int, Rt, mu, st

      mu = 1/sqrt(1 + Rt**2/l_save**2)
      st = eps0(Rt)*eps_1*(1 + z)*(1 - mu)/2
      R_int = phi(Rt)/(Rt**1.25_dp*(1 + Rt**2/l_save**2)**1.5_dp) &
        *sigma_gg(st)*(1 - mu)
    end function R_int

    function eps0(Rt)
      ! cf. eq. (64)
      implicit none
      real(dp) :: eps0, Rt

      eps0 = 2.7e-4_dp*(l_edd/(M_8*eta))**0.25_dp*Rt**(-0.75_dp)
    end function eps0

    function phi(Rt)
      implicit none
      real(dp) :: phi, Rt

      phi = sqrt(1 - R_in/(R_g*Rt))
    end function phi
  end function tau_disk

  function tau_blr(eps_1, z, l_edd, rr, R_g, xi_li, R_li, eps_li, n_li)
    implicit none
    real(dp) :: tau_blr, eps_1, z, l_edd, R_g, rr
    integer :: n_li
    real(dp), dimension(n_li) :: xi_li, R_li, eps_li
    real(dp) :: l_save, l_err
    integer :: l_neval, l_ier

    call qagi(l_int, rr/R_g, 1, epsabs, epsrel, tau_blr, l_err, l_neval, l_ier)
  contains
    function l_int(l)
      implicit none
      real(dp) :: l_int, l
      real(dp) :: mu_err
      integer :: mu_neval, mu_ier

      l_save = l
      call qng(mu_int, -1.0_dp, 1.0_dp, epsabs, epsrel, l_int, mu_err, &
        mu_neval, mu_ier)
    end function l_int

    function mu_int(mu_re)
      implicit none
      real(dp) :: mu_int, mu_re
      real(dp), dimension(n_li) :: integ, x2, mu_star, s

      x2 = (R_li**2 + l_save**2 - 2*l_save*R_li*mu_re)
      mu_star = 1 - (R_li/(R_g*sqrt(x2)))**2*(1 - mu_re**2)
      s = eps_li*eps_1*(1 + z)*(1 - mu_star)/2
      integ = 900*xi_li*l_edd/eps_li/x2*sigma(s)*(1 - mu_star)
      mu_int = sum(integ)
    end function mu_int

    function sigma(s)
      implicit none
      real(dp), dimension(n_li) :: sigma, s
      integer :: j

      do j = 1, n_li
        sigma(j) = sigma_gg(s(j))
      end do
    end function sigma
  end function tau_blr

  function tau_dust(eps_1, z, xi, l_edd, Theta, R_dt, R_g, rr)
    implicit none
    real(dp) :: tau_dust, eps_1, z, xi, l_edd, Theta, R_dt, R_g, rr
    real(dp) :: err
    integer :: neval, ier

    call qagi(integrand, rr/R_g, 1, epsabs, epsrel, tau_dust, err, neval, ier)
  contains
    function integrand(l)
      implicit none
      real(dp) :: integrand, l
      real(dp) :: x, s

      x = sqrt((R_dt/R_g)**2 + l**2)
      s = 2.7_dp*Theta*eps_1*(1 + z)*(1 - l/x)/2
      integrand = 900*xi*l_edd/(2.7_dp*Theta)*(1 - l/x)/x**2*sigma_gg(s)
    end function integrand
  end function tau_dust

  function sigma_gg(s)
    ! sigma_gg/sigma_T
    implicit none
    real(dp) :: sigma_gg, s, beta

    if (s > 4.0_dp) then
      beta = sqrt(1 - 1/s)
      sigma_gg = 3/(8*s)*((3 - beta**4)*log((1 + beta)**2*s) &
        - 2*beta*(2 - beta**2))
    else
      sigma_gg = 0.0_dp
    end if
  end function sigma_gg
end module photoabsorption
