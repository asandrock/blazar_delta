module bh_pair
  ! Functions related to Bethe-Heitler pair production, cf. Chodorowski et al.,
  ! ApJ 400, 181
  use const
  !use quadpack
  implicit none
contains
  pure function S0(k)
    ! total cross section, divided by Thomson cross section to be analogous to
    ! Compton scattering, cf. eqq. (A1, A2).
    implicit none
    real(dp) :: S0, eta, ln2k, _2k
    real(dp), intent(in) :: k

    if (k < 2.0_dp) then
      S0 = 0.0_dp
    else if (k < 4.0_dp) then
      eta = (k - 2)/(k + 2)
      S0 = alpha/4*((k - 2)/k)**3*(1 + eta/12 + (23*eta**2)/40 &
        + (37*eta**3)/120 + (61*eta**4)/192)
    else
      ln2k = log(2*k)
      _2k = 2/k
      S0 = 3*alpha/(8*pi)*((28*ln2k)/9 - 218/27.0_dp + _2k**2*(6*ln2k &
        - 7/2.0_dp + 2/3.0_dp*ln2k**3 - ln2k**2 - pi**2/3*ln2k + 2*zeta3 &
        + pi**/6) - _2k**4*((3*ln2k)/16 + 1/8.0_dp) - _2k**6*(29/(9.0_dp*256) &
        *ln2k - 77/(27.0_dp*512)))
    end if
  end function S0

  pure function S1(k)
    ! K, cf. eq. (3.7-3.10)
    implicit none
    real(dp) :: S1, lnk_1, lnk1
    real(dp), intent(in) :: k
    real(dp), parameter :: a0 = 1.0_dp, a1 = 0.3958_dp, a2 = 0.1000_dp, &
      a3 = 0.00781_dp
    real(dp), parameter :: b0 = -8.778_dp, b1 = 5.512_dp, b2 = -1.614_dp, &
      b3 = 2/3.0_dp

    if (x < 2.0_dp) then
      S1 = 0.0_dp
    else if (x < 1e3_dp) then
      lnk_1 = log(k - 1)
      S1 = 4*me/mp/k*(a0 + a1*lnk_1 + a2*lnk_1**2 + a3*lnk_1**3)
    else
      lnk = log(k)
      S1 = 4*me/mp/k*(b0 + b1*lnk + b2*lnk**2 + b3*lnk**3) &
        /(28/9.0_dp*log(2*k) - 218/27.0_dp)
    end if
  end function S1

  function M0(x)
    implicit none
    real(dp) :: M0, x

    M0 = alpha*re**2*8/(3*sigma_T)*psi(x/2)
  end function M0

  function M1(x)
    implicit none
    real(dp) :: M1, x

    M1 = phi(x/2)/psi(x/2)
  end function M1

  function psi(kappa)
    implicit none
    real(dp) :: psi, kappa

    if (kappa < 2.0_dp) then
      psi = 0.0_dp
    else if (kappa < 4.0_dp) then
      psi = (2*pi)/3*(347/(40*kappa) - 130903/1440.0_dp - (19151 &
        *log(2.0_dp))/48 - (15163*kappa)/480 + (2593*kappa**2)/1920 &
        + 3904/(9*(2 + kappa)**3) - 9688/(15*(2 + kappa)**2) + 10676/(15 &
        *(2 + kappa)) + 1007/48.0_dp*log(kappa) + 189*log(2 + kappa))
    else
      psi = -67/(1728*kappa**4) + 29*log(2*kappa)/(144*kappa**4) &
        + 7/(4*kappa**2) + 3*log(2*kappa)/(2*kappa**2) - 2.71245 &
        + 2*(pi**2/3 - 7 + 4*zeta3)*log(2*kappa) + 2*(6 - pi**2/3) &
        *log(2*kappa)**2 - 4/3.0_dp*log(2*kappa)**3 + 2/3.0_dp &
        *log(2*kappa)**4 - (130*kappa**2)/27 + (14*kappa)**2/9*log(2*kappa)
    end if
  end function psi

  function phi(kappa)
    implicit none
    real(dp) :: phi, kappa, k_2
    real(dp), parameter :: c1 = 0.8048_dp, c2 = 0.1459_dp, c3 = 1.137e-3_dp, &
      c4 = -3.879e-6_dp
    real(dp), parameter :: d0 = -86.07_dp, d1 = 50.96_dp, d2 = -14.45_dp, &
      d3 = 8/3.0_dp, f1 = 2.910_dp, f2 = 78.35_dp, f3 = 1837.0_dp

    if (kappa < 2.0_dp) then
      phi = 0.0_dp
    else if (kappa < 25.0_dp) then
      k_2 = kappa - 2
      phi = pi/12*k_2**4/(1 + c1*k_2 + c2*k_2**2 + c3*k_2**3 + c4*k_2**4)
    else
      lnk = log(kappa)
      phi = kappa*(d0 + d1*lnk + d2*lnk**2 + d3*lnk**3)/(1 - f1/kappa &
        - f2/kappa**2 - f3/kappa**3)
    end if
  end function phi
end module bh_pair
