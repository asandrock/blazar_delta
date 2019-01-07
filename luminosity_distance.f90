module luminosity_distance
  use iso_fortran_env, only: dp=>real64
  use const
  implicit none
contains
  function lum_dist(redshift)
    implicit none
    real(dp) :: lum_dist, redshift
    real(dp) :: result, abserr
    integer :: neval, ifail

    if (redshift <= 0.0_dp) then
      print *, "luminosity_distance: redshift <= 0 is unphysical"
      STOP
    else
      call qags(integrand, 0.0_dp, redshift, epsabs, epsrel, &
          result, abserr, neval, ifail)
      lum_dist = result
    end if
  contains
    function integrand(zeta)
      implicit none
      real(dp) :: integrand
      real(dp), intent(in) :: zeta

      integrand = c/H0*(1 + redshift)/sqrt(Omega_r*(1 + zeta)**4 &
          + Omega_m*(1 + zeta)**3 + Omega_k*(1 + zeta)**2 + Omega_L)
    end function integrand
  end function lum_dist
end module luminosity_distance
