module external_compton
  use compton
  use const
  use gamma_avg
  implicit none
contains
  function ec_iso(eps_s, z, d_L, delta_D, Ne, u)
    ! vF_v for continuous isotropic external radiation field
    ! cf. eq. (54, 55)
    implicit none
    real(dp) :: ec_iso, eps_s, z, d_L, delta_D
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
      g_min = eps_sp/2*(1 + sqrt(1 + 1/(eps*eps_sp)))
      if (g_avg > g_min) then
        integrand = u(eps)/eps**3*Ne(g_avg/delta_D)/g_avg*M3(4*g_avg*eps)
      else 
        integrand = u(eps)/eps**3*Ne(g_min/delta_D)/g_min*M3(4*g_min*eps) &
          *(g_avg/g_min)**1.5_dp
      end if
    end function integrand
  end function ec_iso

  function ec_iso_mono(eps_s, z, d_L, delta_D, Ne, u0, eps0)
    ! same as ec_iso for monochromatic radiation field
    !   u(eps) = u0*delta(eps - eps0)
    ! cf. eq. (57, 58)
    implicit none
    real(dp) :: ec_iso_mono, eps_s, z, d_L, delta_D, u0, eps0, g_min
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: g_avg, eps_sp

    eps_sp = eps_s*(1 + z)
    g_avg = gamma_isotropic(eps_sp, eps0)
    g_min = eps_sp/2*(1 + sqrt(1 + 1/(eps0*eps_sp)))
    if (g_avg > g_min) then
      ec_iso_mono = c*sigma_T*delta_D**3*eps_sp/(32*pi*d_L**2) &
        *u0/eps0**3*Ne(g_avg/delta_D)/g_avg*M3(4*g_avg*eps0)
    else
      ec_iso_mono = c*sigma_T*delta_D**3*eps_sp/(32*pi*d_L**2) &
        *u0/eps0**3*Ne(g_min/delta_D)/g_min*M3(4*g_min*eps0) &
        *(g_avg/g_min)**1.5_dp
    end if
  end function ec_iso_mono

  function ec_point_mono(eps_s, z, d_L, delta_D, Gamm, Ne, L0, eps0, rr) 
    ! v F_v for external Compton scattering on accretion disk photons, assuming
    ! a monochromatic point source radially behind the jet
    implicit none
    real(dp) :: ec_point_mono, eps_s, z, d_L, delta_D, Gamm, L0, eps0, rr
    interface
      function Ne(g)
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: eps_sp, g_min, mu_s, g_avg

    eps_sp = eps_s*(1 + z)
    mu_s = (Gamm*delta_D - 1)/(Gamm*sqrt(1 - 1/Gamm**2)*delta_D)
    g_min = eps_sp/2*(1 + sqrt(1 + 2/(eps0*eps_sp*(1 - mu_s))))
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

  function ec_disk(eps_s, z, d_L, delta_D, Gamm, Ne, M_8, l_edd, eta, rr)
    ! v F_v for external Compton scattering on accretion disk photons, assuming
    ! a Shakura-Sunyaev accretion disk
    ! cf. eq. (63-73)
    implicit none
    real(dp) :: ec_disk, eps_s, z, d_L, delta_D, Gamm, M_8, l_edd, eta, rr
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
    interface
      function Ne(g)
        ! electron distribution
        use const
        real(dp) :: Ne, g
      end function Ne
    end interface
    real(dp) :: R_g, Lumi_Edd, R_in, R_out, mu_max, mu_min, eps_sp, mu_s
    real(dp) :: abserr, ph_save, prefactor
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
    call qng(ph_integrand, 0.0_dp, 2*pi, epsabs, epsrel, ec_disk, abserr, &
      neval, ier)
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
      g_min = eps_sp/2*(1 + sqrt(1 + 2/(eps*eps_sp*(1 - cos_psi))))
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
end module external_compton
