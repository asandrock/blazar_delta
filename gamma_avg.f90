module gamma_avg
  ! Functions to determine <gamma> for general and isotropic seed radiation
  ! fields
  use compton
  !use zeroin
  implicit none
  interface
    function zeroin(ax,bx,f,tol)
      double precision :: ax, bx, tol, zeroin
      interface
        function f(x)
          double precision :: f, x
        end function
      end interface
    end function zeroin
  end interface
contains
  function gamma_general(eps, eps_s, cos_psi)
    ! numerical solution of eq. (30)
    implicit none
    real(dp) :: gamma_general, eps, eps_s, cos_psi
    real(dp) :: g_min, g_max, g_th, g_kn

    ! asymptotic solutions
    g_th = sqrt(4*eps_s/(5*eps*(1 - cos_psi)))
    g_kn = (4*eps_s)/5
    ! search region
    g_min = max(1.0_dp, min(g_th, g_kn)/10)
    g_max = 10*max(g_th, g_kn)

    gamma_general = zeroin(g_min, g_max, eq30, epsrel)
  contains
    function eq30(gg)
      implicit none
      real(dp) :: eq30, gg

      eq30 = eps_s/gg - 1.25_dp*S1(gg*eps*(1 - cos_psi))
    end function eq30    
  end function gamma_general

  function gamma_isotropic(eps_s, eps)
    ! numerical solution of eq. (48)
    implicit none
    real(dp) :: gamma_isotropic, eps_s, eps
    real(dp) :: g_min, g_max, g_th, g_kn

    ! asymptotic solutions
    g_th = sqrt(eps_s/(2*eps))
    g_kn = eps_s/0.9419_dp
    ! search region
    g_min = max(1.0_dp, min(g_th, g_kn)/10)
    g_max = 10*max(g_th, g_kn)

    gamma_isotropic = zeroin(g_min, g_max, eq48, epsrel)
  contains
    function eq48(gg)
      implicit none
      real(dp) :: eq48, gg

      eq48 = eps_s/gg - 1.5_dp*M1(4*gg*eps)
    end function eq48
  end function gamma_isotropic
end module gamma_avg
