module photoabsorption
  ! Photoabsorption by external radiation fields
  use const
  use multidim_integrate, only: dcuhre
  implicit none
contains
  function tau_disk(eps_1, z, l_edd, M_8, eta, rr)
    implicit none
    real(dp) :: tau_disk, eps_1, l_edd, M_8, eta, rr, z
    real(dp) :: l_err, l_save, prefactor, R_g, R_in, R_out
    integer :: l_neval, l_ier, nsub
    real(dp), dimension(2) :: aa, bb
    real(dp), dimension(1) :: res, err

    R_g = 1.5e13_dp*M_8
    R_in = 6*R_g ! Schwarzschild black hole
    R_out = 200*R_g ! simple assumption from Finke (2016)
    prefactor = 1e7_dp*l_edd**0.75_dp*M_8**0.25_dp/eta**0.75_dp
    !call qagi(l_int, rr/R_g, 1, epsabs, epsrel, tau_disk, l_err, l_neval, l_ier)
    aa = [0.0_dp, log(R_in/R_g)]
    bb = [R_g/rr, log(R_out/R_g)]
    call dcuhre(2, 1, aa, bb, int(1e3), int(1e5), funsub, epsabs, epsrel, 0, 0, &
      res, err, l_neval, l_ier, nsub)
    tau_disk = res(1)
    !if (l_ier /= 0) tau_disk = -1.0_dp
  contains
    subroutine funsub(ndim, z, nfun, f)
      implicit none
      integer, intent(in) :: ndim, nfun
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: f(:)

      l_save = 1/z(1)
      f = R_int(z(2))*l_save**2
    end subroutine funsub

    function l_int(l)
      implicit none
      real(dp) :: l_int, l
      real(dp) :: R_err
      integer :: R_neval, R_ier

      l_save = l
      call qng(R_int, log(R_in/R_g), log(R_out/R_g), epsabs, epsrel, l_int, &
        R_err, R_neval, R_ier)
    end function l_int

    function R_int(log_Rt)
      ! corrected according to Dermer et al. (2009)
      implicit none
      real(dp), intent(in) :: log_Rt
      real(dp) :: R_int, Rt, mu, st

      Rt = exp(log_Rt)

      mu = 1/sqrt(1 + Rt**2/l_save**2)
      st = eps0(Rt)*eps_1*(1 + z)*(1 - mu)/2
      R_int = prefactor*phi(Rt)**0.25_dp/(Rt**1.25_dp &
        *(1 + Rt**2/l_save**2)**1.5_dp) &
        *sigma_gg(st)*(1 - mu)/l_save**2 &
        *Rt
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

  function tau_blr(eps_1, z, l_edd, rr, R_g, xis_li, Rs_li, epss_li, n_li)
    implicit none
    real(dp) :: tau_blr, eps_1, z, l_edd, R_g, rr
    integer :: n_li
    real(dp), dimension(n_li) :: xis_li, Rs_li, epss_li
    real(dp) :: xi_li, R_li, eps_li
    real(dp) :: l_save, l_err
    integer :: l_neval, l_ier, nsub, j
    real(dp), dimension(2) :: aa, bb
    real(dp), dimension(1) :: res, err

    !call qagi(l_int, rr/R_g, 1, epsabs, epsrel, tau_blr, l_err, l_neval, l_ier)
    !call qng(l_int, 0.0_dp, R_g/rr, epsabs, epsrel, tau_blr, l_err, l_neval, l_ier)
    aa = [0.0_dp, -1.0_dp]
    bb = [R_g/rr, 1.0_dp]
    tau_blr = 0.0_dp
    do j = 1, n_li
      R_li = Rs_li(j)
      xi_li = xis_li(j)
      eps_li = epss_li(j)
      call dcuhre(2, 1, aa, bb, int(1e3), int(1e5), funsub, epsabs, epsrel, 0, 0, &
        res, err, l_neval, l_ier, nsub)
      tau_blr = tau_blr + res(1)
    end do
  contains
    subroutine funsub(ndim, z, nfun, f)
      implicit none
      integer, intent(in) :: ndim, nfun
      real(dp), intent(in) :: z(:)
      real(dp), intent(out) :: f(:)

      l_save = 1/z(1)
      f = mu_int(z(2))
    end subroutine funsub

    function l_int(l1)
      implicit none
      real(dp) :: l_int, l1
      real(dp) :: mu_err
      integer :: mu_neval, mu_ier

      l_save = 1/l1
      !call qags(mu_int, -1.0_dp, 1.0_dp, epsabs, epsrel, l_int, mu_err, &
      !  mu_neval, mu_ier)
      call qng(mu_int, -1.0_dp, 1.0_dp, epsabs, epsrel, l_int, mu_err, &
        mu_neval, mu_ier)
    end function l_int

    function mu_int(mu_re)
      implicit none
      real(dp) :: mu_int
      real(dp) :: mu_re
      real(dp) :: integ, x2, mu_star, s, Rg2_x2, aux

      x2 = ((R_li/R_g)**2 + l_save**2 - 2*l_save*R_li/R_g*mu_re)
      Rg2_x2 = R_li**2 + (R_g*l_save)**2 - 2*(l_save*R_g)*R_li*mu_re
      if (abs(mu_re) < 1) then
        aux = max(0.0_dp, 1 - R_li**2/Rg2_x2*(1 - mu_re**2))
      else
        aux = 0.0_dp
      end if
      mu_star = sqrt(aux)
      s = eps_li*eps_1*(1 + z)*(1 - mu_star)/2
      integ = 900*xi_li*l_edd/eps_li/x2*sigma_gg(s)*(1 - mu_star)
      mu_int = integ*l_save**2
    end function mu_int
  end function tau_blr

  function tau_dust(eps_1, z, xi, l_edd, Theta, R_dt, R_g, rr)
    implicit none
    real(dp) :: tau_dust, eps_1, z, xi, l_edd, Theta, R_dt, R_g, rr
    real(dp) :: err
    integer :: neval, ier

    call qagi(integrand, rr/R_g, 1, epsabs, epsrel, tau_dust, err, neval, ier)
    !if (ier /= 0) tau_dust = -1.0_dp
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

  elemental function sigma_gg(s)
    ! sigma_gg/sigma_T
    implicit none
    real(dp) :: sigma_gg, beta
    real(dp), intent(in) :: s

    if (s > 1.0_dp) then
      beta = sqrt(1 - 1/s)
      sigma_gg = 3/(8*s)*((3 - beta**4)*log((1 + beta)**2*s) &
        - 2*beta*(2 - beta**2))
    else
      sigma_gg = 0.0_dp
    end if
  end function sigma_gg
end module photoabsorption
