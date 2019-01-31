module const
  use iso_fortran_env, only: dp=>real64
  implicit none
  real(dp), parameter :: pi = 4*atan(1.0_dp), &
    c = 2.99792458e10_dp, & ! cm/s, speed of light
    h = 6.62607004e-27_dp, & ! erg s Planck's constant
    hbar = 1.0545718e-27_dp, & ! erg s, reduced Planck's constant
    sigma_T = 6.652458558e-25, & ! cm^2, Thomson cross section
    e =  4.80320425e-10_dp, & ! statCoulomb, elementary charge
    TeV = 1.602176462_dp, & ! erg, conversion factor of TeV to erg
    me = 9.10938356e-28_dp, & ! g, electron mass
    mp = 1.672621777e-24_dp, & ! g, proton mass
    m_pion = (139.57018e-6_dp + 134.9766e-6_dp)/2*TeV/c**2, & ! g, pion mass
    B_cr = me**2*c**3/(e*hbar), & ! = 4.414e13_dp G, critical magnetic field
    me_c2 = 0.510998928e6_dp, & ! eV, conversion factor from eV to dimensionless
    ! energy
    re = 2.8179403227e-13_dp, & ! cm, classical electron radius
    lambda_C = h/(me*c), & ! electron Compton wavelength
    alpha = 1/137.035999074_dp, & ! fine structure constant
    mu = me/mp, & ! electron-proton mass ratio
    zeta3 = 1.20206_dp, & ! Riemann zeta function at 3
    zeta4 = 1.082323_dp, & ! Riemann zeta function at 4
    microbarn = 1e-30_dp ! cm^2
  real(dp), parameter :: k = 1.3806488e-16_dp, & ! erg/K, Boltzmann constant
    a_SB = 8*pi**5*k**4/(15*h**3*c**3) ! Stefan-Boltzmann radiation energy density constant
  real(dp), parameter :: Omega_m = 0.3089_dp, Omega_L = 0.6911_dp, &
    Omega_k = 1 - Omega_m - Omega_L, Omega_r = 0.0_dp, & ! Cosmology parameters
    Mpc = 3.08567758149137e24_dp, & ! cm, Megaparsec in cm
    H0 = 67.74e5_dp/Mpc ! cm/s/Mpc, Hubble constant
  real(dp), parameter :: epsabs = 1e-60_dp, epsrel = 1e-5_dp
end module const
