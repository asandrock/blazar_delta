module dilog
  ! Real part of dilogarithm of a real argument, following the convention
  ! $$
  !   Li_2(x) = -\Re \int_0^x \frac{\ln(1 - t)}{t} dt,
  ! $$
  ! which coincides with the convention of Mathematica and GNU Maxima. The
  ! Spence function used e.g. in Abramowitz & Stegun and SciPy is connected
  ! to this one as
  ! $$
  !   Spence(x) = Li_2(1 - x).
  ! $$
  ! The implementation is based on the implementation in the GSL.
  implicit none
  public :: li2
  private
    integer, parameter :: dp = kind(1.0d0)
    real(dp), parameter :: pi = 4*atan(1.0_dp)
contains
  elemental function li2(x)
    implicit none
    real(dp) :: li2
    real(dp), intent(in) :: x

    if (x >= 0.0_dp) then
      li2 = dilog_xge0(x)
    else
      li2 = -dilog_xge0(-x) + dilog_xge0(x**2)/2
    end if
  end function li2

  elemental function dilog_xge0(x)
    ! Dilogarithm for real x, where x >= 0.
    implicit none
    real(dp) :: dilog_xge0
    real(dp), intent(in) :: x
    real(dp) :: eps, lne
    real(dp), dimension(0:8) :: c
    integer :: k

    if (x > 2.0_dp) then
      dilog_xge0 = pi**2/3 - dilog_series_2(1/x) - log(x)**2/2
    else if (x > 1.01_dp) then
      dilog_xge0 = pi**2/6 + dilog_series_2(1 - 1/x) - log(x)*(log(1 - 1/x) &
        + log(x)/2)
    else if (x > 1.0_dp) then
      ! series around x = 1
      eps = x - 1
      lne = log(eps)
      c(0) = pi**2/6
      do k = 1, 8
        c(k) = (-1)**(1 + k)*(1 - k*lne)/k**2
      end do
      dilog_xge0 = c(0) + eps*(c(1) + eps*(c(2) + eps*(c(3) + eps*(c(4) &
        + eps*(c(5) + eps*(c(6) + eps*(c(7) + eps*c(8))))))))
    else if (x == 1.0_dp) then
      dilog_xge0 = pi**2/6
    else if (x > 0.5_dp) then
      dilog_xge0 = pi**2/6 - dilog_series_2(1 - x) - log(x)*log(1 - x)
    else if (x > 0.25_dp) then
      dilog_xge0 = dilog_series_2(x)
    else if (x > 0.0_dp) then
      dilog_xge0 = dilog_series_1(x)
    else ! x = 0
      dilog_xge0 = 0.0_dp
    end if
  end function dilog_xge0

  elemental function dilog_series_1(x)
    implicit none
    real(dp) :: dilog_series_1
    real(dp), intent(in) :: x
    ! Series for real dilogarithm, converging rapidly for |x| < 1/2.
    ! $$
    !   \sum_{k = 1}^\infty \frac{x^k}{k^2}
    ! $$
    integer, parameter :: kmax = 1000
    integer :: k
    real(dp) :: series

    series = 0.0_dp
    do k = 1, kmax
      series = series + x**k/k**2
    end do 
    dilog_series_1 = series
  end function dilog_series_1

  elemental function dilog_series_2(x)
    implicit none
    real(dp) :: dilog_series_2
    real(dp), intent(in) :: x
    ! Associated series representation of the dilogarithm for values
    ! -1 < x < 1
    ! $$
    !   Li_2(x) = 1 + \frac{1 - x}{x} \ln (1 - x)
    !           + \sum_{k = 1}^\infty \frac{x^k}{k^2 (1 + k)}
    ! $$
    integer, parameter :: kmax = 100
    integer :: k
    real(dp) :: series

    series = 0.0_dp
    do k = 1, kmax
      series = series + x**k/(k**2*(1 + k))
    end do
    dilog_series_2 = 1 + (1 - x)/x*log(1 - x) + series
  end function dilog_series_2
end module dilog
