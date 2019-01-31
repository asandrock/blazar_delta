module ssc
  use const
  use gamma_avg
  use dilog
  !use quadpack
  implicit none
contains
  function f_syn(eps, delta_D, z, d_L, Ne, B, g1, g2)
    implicit none
    real(dp) :: f_syn, eps, delta_D, z, d_L, B, g1, g2
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: prefactor, eps_p, abserr
    integer :: neval, ier

    eps_p = (1 + z)*eps/delta_D
    prefactor = sqrt(3.0_dp)*delta_D**4*eps_p*e**3*B/(4*pi*h*d_L**2)
    call qng(integrand, log(g1), log(g2), epsabs, epsrel, f_syn, abserr, &
      neval, ier)
  contains
    function integrand(log_g)
      implicit none
      real(dp) :: integrand, log_g, gamm, x

      gamm = exp(log_g)
      x = 4*pi*eps_p*me**2*c**3/(3*e*B*h*gamm**2)
      integrand = prefactor*Ne(gamm)*G_syn(x)*gamm
    end function integrand

    function G_syn(x)
      ! Exact synchrotron function, Crusius & Schlickeiser (1986):
      !         x   pi              infty 
      ! G(x) = --- int dth sin(th)   int    dt K_{5/3}(t)
      !         2   0              x/sin(th)
      ! Parametrization from Aharonian, Kelner & Prosekin, PRD (2010) 82:043002.
      implicit none
      real(dp) :: G_syn, x, y
 
      y = x**(1/3.0_dp)
      G_syn = 1.808_dp*y/sqrt(1 + 3.4_dp*y**2)*(1 + 2.21_dp*y**2 &
        + 0.347_dp*y**4)/(1 + 1.353_dp*y**2 + 0.217_dp*y**4)*exp(-x)
    end function G_syn
  end function f_syn

  function f_ssc(eps_s, delta_D, z, d_L, Ne, B, R, g1, g2)
    implicit none
    real(dp) :: f_ssc, eps_s, delta_D, z, d_L, B, R, g1, g2
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: u_prefac, ssc_prefac, eps_sp, eps_B, abserr
    integer :: neval, ier

    eps_sp = (1 + z)*eps_s/delta_D
    eps_B = B/B_cr
    u_prefac = sigma_T*B**2/(8*pi)/(2*pi*R**2*eps_B)
    ssc_prefac = delta_D**4/(4*pi*d_L**2) &
      *c*sigma_T*eps_sp/(32*pi)*4*pi

    call qng(integrand, log(g1**2*eps_B), log(g2**2*eps_B), epsabs, epsrel, &
      f_ssc, abserr, neval, ier)
  contains
    function integrand(log_eps)
      implicit none
      real(dp) :: integrand, log_eps, eps, g_min, g_avg

      eps = exp(log_eps)
      g_avg = gamma_isotropic(eps_sp, eps)
      g_min = max(eps_sp/2*(1 + sqrt(1 + 1/(eps*eps_sp))), g1)
      if (g_min < g_avg) then
        !integrand = ssc_prefac*f_syn_delta(eps)/eps**4*Ne(g_avg)/g_avg &
        !  *M3(4*g_avg*eps) &
        !  *eps
        integrand = ssc_prefac*u(eps)/eps**3*Ne(g_avg)/g_avg*M3(4*g_avg*eps)*eps
      else
        !integrand = ssc_prefac*f_syn_delta(eps)/eps**4*Ne(g_min)/g_min &
        !  *M3(2*g_min*eps)*(g_avg/g_min)**1.5_dp &
        !  *eps
        integrand = ssc_prefac*u(eps)/eps**3*Ne(g_min)/g_min*M3(4*g_min*eps) &
          *(g_avg/g_min)**1.5_dp*eps
      end if
    end function integrand

    function u(eps)
      implicit none
      real(dp) :: u, eps, g_s

      g_s = sqrt(eps/eps_B)
      u = u_prefac*g_s*Ne(g_s)
    end function u
  end function f_ssc

  function tau_ssc(eps, delta_D, z, d_L, Ne, B, R, g1, g2)
    implicit none
    real(dp) :: tau_ssc, eps, delta_D, z, d_L, B, R, g1, g2
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: eps_p, eps_B, prefactor, syn_prefac, abserr
    integer :: neval, ier

    eps_p = eps*(1 + z)/delta_D
    eps_B = B/B_cr
    syn_prefac = delta_D**4/(6*pi*d_L**2)*c*sigma_T*B**2/(8*pi)
    prefactor = 9*d_L**2*sigma_T/(8*me*c**3*R*delta_D**5*eps_p**2)
    call qng(integrand, log(eps_B*g1**2), log(eps_B*g2**2), epsabs, epsrel, &
      tau_ssc, abserr, neval, ier)
  contains
    function integrand(log_e)
      implicit none
      real(dp) :: integrand
      real(dp), intent(in) :: log_e

      real(dp) :: e, s0

      e = exp(log_e)
      s0 = e*eps_p
      integrand = prefactor*f_syn_delta(e)/e**4*phi_bar(s0)*e
    end function integrand

    function phi_bar(s0)
      implicit none
      real(dp) :: phi_bar
      real(dp), intent(in) :: s0
 
      real(dp) :: beta0, w0, L
 
      if (s0 > 1) then
        beta0 = sqrt(1 - 1/s0)
        w0 = (1 + beta0)**2*s0
        L = -pi**2/12 - li2(-w0)
        phi_bar = (1 + beta0**2)*s0*log(w0) - beta0**2*log(w0) &
            - 4*beta0*s0 + 2*beta0 + 4*log(w0)*log(1 + w0) - 4*L
      else
        phi_bar = 0.0_dp
      end if
    end function phi_bar

    function f_syn_delta(eps)
      implicit none
      real(dp) :: f_syn_delta, eps, g_s

      g_s = sqrt(eps/eps_B)
      f_syn_delta = syn_prefac*g_s**3*Ne(g_s)
    end function f_syn_delta
  end function tau_ssc
end module ssc
