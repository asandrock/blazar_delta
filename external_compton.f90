module external_compton
  use compton
  use const
  use gamma_avg
  implicit none
contains
  function ec_iso(eps_s, z, d_L, delta_D, Ne, u, g1)
    ! vF_v for continuous isotropic external radiation field
    ! cf. eq. (54, 55)
    implicit none
    real(dp) :: ec_iso, eps_s, z, d_L, delta_D, g1
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne

      function u(eps)
        use const
        real(dp) :: u, eps
      end function u
    end interface

    real(dp) :: abserr, prefactor, eps_sp
    integer :: neval, ier

    eps_sp = eps_s*(1 + z)
    prefactor = c*sigma_T*delta_D**3*eps_sp/(32*pi*d_L**2)
    call qagi(integrand, 0.0_dp, 1, epsabs, epsrel, ec_iso, abserr, neval, ier)
  contains
    function integrand(eps)
      implicit none
      real(dp) :: integrand, eps, g_avg, g_min

      g_avg = gamma_isotropic(eps_s*(1 + z), eps)
      g_min = max(eps_sp/2*(1 + sqrt(1 + 1/(eps*eps_sp))), g1*delta_D)
      if (g_avg > g_min) then
        integrand = u(eps)/eps**3*Ne(g_avg/delta_D)/g_avg*M3(4*g_avg*eps)
      else
        integrand = u(eps)/eps**3*Ne(g_min/delta_D)/g_min*M3(4*g_min*eps) &
          *(g_avg/g_min)**1.5_dp
      end if
    end function integrand
  end function ec_iso

  function ec_iso_mono(eps_s, z, d_L, delta_D, Ne, u0, eps0, g1)
    ! same as ec_iso for monochromatic radiation field
    !   u(eps) = u0*delta(eps - eps0)
    ! cf. eq. (57, 58)
    implicit none
    real(dp) :: ec_iso_mono, eps_s, z, d_L, delta_D, u0, eps0, g_min, g1
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: g_avg, eps_sp

    eps_sp = eps_s*(1 + z)
    g_avg = gamma_isotropic(eps_sp, eps0)
    g_min = max(eps_sp/2*(1 + sqrt(1 + 1/(eps0*eps_sp))), g1*delta_D)
    if (g_avg > g_min) then
      ec_iso_mono = c*sigma_T*delta_D**3*eps_sp/(32*pi*d_L**2) &
        *u0/eps0**3*Ne(g_avg/delta_D)/g_avg*M3(4*g_avg*eps0)
    else
      ec_iso_mono = c*sigma_T*delta_D**3*eps_sp/(32*pi*d_L**2) &
        *u0/eps0**3*Ne(g_min/delta_D)/g_min*M3(4*g_min*eps0) &
        *(g_avg/g_min)**1.5_dp
    end if
  end function ec_iso_mono

  function ec_point_mono(eps_s, z, d_L, delta_D, Gamm, Ne, L0, eps0, rr, g1) 
    ! v F_v for external Compton scattering on accretion disk photons, assuming
    ! a monochromatic point source radially behind the jet
    implicit none
    real(dp) :: ec_point_mono, eps_s, z, d_L, delta_D, Gamm, L0, eps0, rr, g1
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: eps_sp, g_min, mu_s, g_avg

    eps_sp = eps_s*(1 + z)
    mu_s = (Gamm*delta_D - 1)/(Gamm*sqrt(1 - 1/Gamm**2)*delta_D)
    g_min = max(eps_sp/2*(1 + sqrt(1 + 2/(eps0*eps_sp*(1 - mu_s)))), g1*delta_D)
    g_avg = gamma_general(eps0, eps_sp, mu_s)
    if (g_avg > g_min) then
      ec_point_mono = eps_sp*L0*sigma_T*delta_D**3 &
        /(16*pi**2*rr**2*d_L**2*eps0**2)*S3(g_avg*eps0*(1 - mu_s)) &
        *Ne(g_avg/delta_D)
    else
      ec_point_mono = eps_sp*L0*sigma_T*delta_D**3 &
        /(16*pi**2*rr**2*d_L**2*eps0**2)*S3(g_min*eps0*(1 - mu_s)) &
        *Ne(g_min/delta_D)*(g_avg/g_min)**2
    end if
  end function ec_point_mono

  function ec_disk(eps_s, z, d_L, delta_D, Gamm, Ne, M_8, l_edd, eta, rr, g1)
    ! v F_v for external Compton scattering on accretion disk photons, assuming
    ! a Shakura-Sunyaev accretion disk
    ! cf. eq. (63-73)
    implicit none
    real(dp) :: ec_disk, eps_s, z, d_L, delta_D, Gamm, M_8, l_edd, eta, rr, g1
    ! ec_disk : vF_v
    ! eps_s : observed dimensionless photon energy
    ! z : redshift
    ! d_L : luminosity distance
    ! delta_D : Doppler factor
    ! Gamm : bulk Lorentz factor
    ! M_8 : black hole mass in 10^8 solar masses
    ! l_edd : Eddington ratio L_disk/L_Edd
    ! eta : accretion efficiency
    ! rr : black hole-blob distance
    ! g1 : minimum electron Lorentz factor
    interface
      function Ne(g)
        ! electron distribution
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: R_g, Lumi_Edd, R_in, R_out, mu_max, mu_min, eps_sp, mu_s
    real(dp) :: integral, abserr, ph_save, prefactor
    integer :: neval, ier

    eps_sp = (1 + z)*eps_s
    mu_s = (Gamm*delta_D - 1)/(Gamm*sqrt(1 - 1/Gamm**2)*delta_D)

    R_g = 1.5e13_dp*M_8
    Lumi_Edd = 1.26e46_dp*M_8
    R_in = 6*R_g ! Schwarzschild black hole
    R_out = 200*R_g ! simple assumption from Finke (2016)

    prefactor = 3*eps_s*(1 + z)*delta_D**3*l_edd*Lumi_Edd*R_g*sigma_T &
      /(64*pi**3*d_L**2*eta*rr**3)

    mu_min = 1/sqrt(1 + (R_out/rr)**2)
    mu_max = 1/sqrt(1 + (R_in/rr)**2)
    call qng(ph_integrand, 0.0_dp, 2*pi, epsabs, epsrel, integral, abserr, &
      neval, ier)
    ec_disk = integral
  contains
    function ph_integrand(ph)
      implicit none
      real(dp) :: ph_integrand, ph
      real(dp) :: ph_err
      integer :: ph_neval, ph_ier

      ph_save = ph
      call qng(mu_integrand, mu_min, mu_max, epsabs, epsrel, ph_integrand, &
        ph_err, ph_neval, ph_ier)
    end function ph_integrand

    function mu_integrand(mu)
      implicit none
      real(dp) :: mu_integrand, mu
      real(dp) :: cos_psi, g_avg, R, eps, g_min

      cos_psi = mu*mu_s + sqrt(1 - mu**2)*sqrt(1 - mu_s**2)*cos(ph_save)
      R = rr/mu*sqrt(1 - mu**2)
      eps = eps0(R)
      g_avg = gamma_general(eps, eps_sp, cos_psi)
      g_min = max(eps_sp/2*(1 + sqrt(1 + 2/(eps*eps_sp*(1 - cos_psi)))), &
        g1*delta_D)
      if (g_avg > g_min) then
        mu_integrand = prefactor*phi(R)/((1/mu**2 - 1)**1.5_dp*eps**2) &
          *S3(g_avg*eps*(1 - cos_psi))*Ne(g_avg/delta_D)
      else
        mu_integrand = prefactor*phi(R)/((1/mu**2 - 1)**1.5_dp*eps0(R)**2) &
          *S3(g_min*eps0(R)*(1 - cos_psi))*Ne(g_min/delta_D) &
          *(g_avg/g_min)**2
      end if
    end function mu_integrand

    function phi(R)
      ! cf. eq. (66)
      implicit none
      real(dp) :: phi, R

      phi = sqrt(1 - R_in/R)
    end function phi

    function eps0(R)
      ! cf. eq. (64)
      implicit none
      real(dp) :: eps0, R

      eps0 = 2.7e-4_dp*(l_edd/(M_8*eta))**0.25_dp*(R/R_g)**(-0.75_dp)
    end function eps0
  end function ec_disk

  function f_disk(eps_obs, z, d_L, M_8, l_edd, eta)
    ! emission of the accretion disk
    implicit none
    real(dp) :: f_disk, eps_obs, z, d_L, M_8, l_edd, eta
    real(dp) ::  eps, eps_max, Lumi_Edd, R_in, R_g

    eps = eps_obs*(1 + z)
    R_g = 1.5e13_dp*M_8
    Lumi_Edd = 1.26e46_dp*M_8
    R_in = 6*R_g ! Schwarzschild black hole
    eps_max = 2.7e-4*(l_edd/(M_8*eta))**0.25_dp*(R_in/R_g)**(-0.75_dp)
    f_disk = 1.12_dp*l_edd*Lumi_Edd/(4*pi*d_L**2)*(eps/eps_max)**(4/3.0_dp) &
      *exp(-eps/eps_max)
  end function f_disk

  function ec_blr(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr, xi_li, R_li, &
    eps_li, n_li, g1)
    implicit none
    real(dp) :: ec_blr, eps_s, z, d_L, delta_D, Gamm, L_disk, rr, g1
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    integer, intent(in) :: n_li
    real(dp), dimension(n_li) :: xi_li, R_li, eps_li

    real(dp) :: prefactor, cos_ph, ph_int, ph_err, mu_s, eps_sp
    integer :: ph_neval, ph_ier

    eps_sp = (1 + z)*eps_s
    mu_s = (Gamm*delta_D - 1)/(Gamm*sqrt(1 - 1/Gamm**2)*delta_D)
    prefactor = L_disk*sigma_T*delta_D**3/(80*pi**3*d_L**2)
    call qng(ph_integrand, 0.0_dp, 2*pi, epsabs, epsrel, ph_int, ph_err, &
      ph_neval, ph_ier)
    ec_blr = ph_int
  contains
    function ph_integrand(phi)
      implicit none
      real(dp) :: ph_integrand, phi, mu_err
      integer :: mu_neval, mu_ier

      cos_ph = cos(phi)
      call qng(mu_integrand, -1.0_dp, 1.0_dp, epsabs, epsrel, ph_integrand, &
        mu_err, mu_neval, mu_ier)
    end function ph_integrand

    function mu_integrand(mu_re)
      implicit none
      real(dp) :: mu_integrand, mu_re
      real(dp), dimension(n_li) :: g_avg, g_min, cos_psi, mu_star, x2, integ
      integer :: j

      x2 = R_li**2 + rr**2 - 2*rr*r_li*mu_re
      mu_star = 1 - R_li**2/x2*(1 - mu_re**2)
      cos_psi = mu_star*mu_s + sqrt(1 - mu_star**2)*sqrt(1 - mu_s**2) &
        *cos_ph

      do j = 1, n_li
        g_avg(j) = gamma_general(eps_li(j), eps_sp, cos_psi(j))
      end do
      g_min = max(eps_sp/2*(1 + sqrt(1 + 2/(eps_li*eps_sp*(1 - cos_psi)))), &
        g1*delta_D)

      where (g_avg > g_min)
        integ = prefactor*xi_li*(eps_sp/eps_li**2) &
          /x2*S3(g_avg*eps_li*(1 - cos_psi)) &
          *N_e(g_avg/delta_D)
      else where
        integ = prefactor*xi_li*(eps_sp/eps_li**2) &
          /x2*S3(g_min*eps_li*(1 - cos_psi)) &
          *N_e(g_min/delta_D)*(g_avg/g_min)**2
      end where
      mu_integrand = sum(integ)
    end function mu_integrand

    function N_e(g)
      implicit none
      real(dp), dimension(n_li) :: N_e, g
      integer :: k

      do k = 1, n_li
        N_e(k) = Ne(g(k))
      end do
    end function
  end function ec_blr

  function ec_dust(eps_s, z, d_L, delta_D, Gamm, Ne, L_disk, rr, xi_dt, R_dt, &
    Theta, g1)
    implicit none
    real(dp) :: ec_dust, eps_s, z, d_L, delta_D, Gamm, L_disk, rr, xi_dt, &
      R_dt, Theta, g1
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface

    real(dp) :: prefactor, mu_s, eps_sp, x2, abserr, eps0
    integer :: neval, ier

    eps_sp = (1 + z)*eps_s
    mu_s = (Gamm*delta_D - 1)/(Gamm*sqrt(1 - 1/Gamm**2)*delta_D)
    x2 = R_dt**2 + rr**2
    eps0 = 2.7_dp*Theta
    prefactor = eps_sp*xi_dt*L_disk*sigma_T*delta_D**3/(80*pi**3 &
      *eps0**2*x2*d_L**2)
    call qng(phi_int, 0.0_dp, 2*pi, epsabs, epsrel, ec_dust, abserr, neval, ier)
  contains
    function phi_int(phi)
      implicit none
      real(dp) :: phi_int, phi
      real(dp) :: g_avg, g_min, cos_psi

      cos_psi = rr*mu_s/sqrt(x2) + sqrt(1 - rr**2/x2)*sqrt(1 - mu_s**2)*cos(phi)
      g_avg = gamma_general(eps0, eps_sp, cos_psi)
      g_min = max(eps_sp/2*(1 + sqrt(1 + 2/(eps0*eps_sp*(1 - cos_psi)))), &
        g1*delta_D)
      if (g_min < g_avg) then
        phi_int = prefactor*S3(eps0*g_avg*(1 - cos_psi)) &
          *Ne(g_avg/delta_D)
      else
        phi_int = prefactor*S3(eps0*g_min*(1 - cos_psi)) &
          *Ne(g_min/delta_D)*(g_avg/g_min)**2
      end if
    end function phi_int
  end function ec_dust

  function f_dust(eps, z, d_L, Theta, L_disk, xi_dt)
    implicit none
    real(dp) :: f_dust, eps, z, d_L, Theta, L_disk, xi_dt, eps_bh

    eps_bh = eps*(1 + z)
    f_dust = L_disk*xi_dt/(4*pi*d_L**2)*15/pi**4*(eps/Theta)**4 &
      *exp(-eps/Theta)/(1 - exp(-eps/Theta))
  end function f_dust
end module external_compton
