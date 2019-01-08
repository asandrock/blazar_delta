module compton
  ! Functions related to Compton scattering
  use const
  use dilog
  !use quadpack
  implicit none
contains
  function S0(x)
    ! cf. eq. (11-13)
    implicit none
    real(dp) :: S0, x

    if (x < 1e-2_dp) then
      S0 = 1 - 2*x + (26*x**2)/5 - (133*x**3)/10 + (1444*x**4)/35
    else if (x > 1e3_dp) then
      S0 = 3/(8*x)*(log(x) + 0.5_dp*(1 + log(4.0_dp)))
    else
      S0 = 3/(8*x**2)*(4 + 2*x**2*(1 + x)/(1 + 2*x)**2 + (x**2 - 2*x - 2)/x &
        *log(1 + 2*x))
    end if
  end function S0

  function M0(x)
    ! cf. eq. (15-17)
    implicit none
    real(dp) :: M0, x

    if (x < 1e-2_dp) then
      M0 = x**2/3 - (2*x**3)/9 + (13*x**4)/60 - (133*x**5)/600
    else if (x > 1e2_dp) then
      M0 = x/2*(log(x) - 0.5_dp) + 3.114_dp*log(x) - log(x/2)**2 - 6.559_dp &
        + 2/x*(2*log(x) + 3.249_dp)
    else
      M0 = (2*(x + 1)*x*li2(-x) - (x/2 + 1)*x**2/2 + (x/2 + 4)*(x + 1)**2 &
        *log(x + 1))/(x*(x + 1)) - 4
    end if
  end function M0

  function S1(x)
    ! cf. eq. (18-22)
    implicit none
    real(dp) :: S1, x

    if (x < 1e-2_dp) then
      S1 = x
    else if (x > 1e2_dp) then
      S1 = 1 - 4/(3*log(2*x) + 0.5_dp)
    else
      S1 = (S0(x) - 3/(8*x**3)*(x**2/3*(((1 + 2*x)**3 - 1)/(1 + 2*x)**3) &
        + 2*x*(x**2 - x - 1)/(1 + 2*x) + log(1 + 2*x)))/S0(x)
    end if
  end function S1

  function M1(x)
    ! cf. eq. (24)
    implicit none
    real(dp) :: M1, x
    real(dp) :: err, aux1, aux2, y_max
    integer :: neval, ier

    y_max = x/(1 + x)
    if (x < 1e-3_dp) then
      M1 = x/3
    !else if (x > 1e3_dp) then
    !  M1 = 0.6279_dp
    else
      call qags(int_num, 0.0_dp, y_max, epsabs, epsrel, aux1, err, neval, ier)
      if (ier /= 0) aux1 = 0
      call qags(int_den, 0.0_dp, y_max, epsabs, epsrel, aux2, err, neval, ier)
      if (ier /= 0) aux2 = 0
      M1 = aux1/aux2
    end if
  contains
    function int_num(y)
      implicit none
      real(dp) :: int_num, y

      int_num = y*JC(x, y)
    end function int_num

    function int_den(y)
      implicit none
      real(dp) :: int_den, y

      int_den = JC(x, y)
    end function int_den

    function JC(x, y)
      ! cf. eq. (25, 26)
      implicit none
      real(dp) :: JC, x, y, w, xw
 
      w = y/(x*(1 - y))
      xw = y/(1 - y)
      JC = 2*w*log(w) + (1 + 2*w)*(1 - w) + xw**2/(2*(1 + xw))*(1 - w)
    end function JC
  end function M1

  function S2(x)
    ! cf. eq. (32)
    implicit none
    real(dp) :: S2, x

    if (x < 1e-3_dp) then
      S2 = 0.5_dp
    else
      S2 = 1/(1 + x/S1(x)*dS1_dx(x))
    end if
  contains
    function dS1_dx(x)
      ! calculated with GNU Maxima
      implicit none
      real(dp) :: dS1_dx, x

    dS1_dx =(2.0E+0*(192*x**7*log(2*x+1)**2+384*x**6*log(2*x+1)**2+144*x**5*l&
&og(2*x+1)**2-240*x**4*log(2*x+1)**2-300*x**3*log(2*x+1)**2-144*x*&
&*2*log(2*x+1)**2-33*x*log(2*x+1)**2-3*log(2*x+1)**2+112*x**8*log(&
&2*x+1)+896*x**7*log(2*x+1)+2268*x**6*log(2*x+1)+2852*x**5*log(2*x&
&+1)+1996*x**4*log(2*x+1)+780*x**3*log(2*x+1)+156*x**2*log(2*x+1)+&
&12*x*log(2*x+1)+128*x**9-664*x**8-2464*x**7-3296*x**6-2348*x**5-9&
&24*x**4-180*x**3-12*x**2))/(3.0E+0*(2*x+1)**2*(4*x**4*log(2*x+1)-&
&4*x**3*log(2*x+1)-15*x**2*log(2*x+1)-10*x*log(2*x+1)-2*log(2*x+1)&
&+2*x**4+18*x**3+16*x**2+4*x)**2) 
    end function dS1_dx
  end function S2

  function S3(x)
    implicit none
    real(dp) :: S3, x

    S3 = x*S0(x)*S2(x)
  end function S3

  function M2(x)
    ! cf. eq. (45)
    implicit none
    real(dp) :: M2, x

    M2 = 1/(1 + x/M1(x)*dM1_dx(x))
  end function M2

  function dM1_dx(x)
    implicit none
    real(dp) :: dM1_dx, x
    real(dp) :: aux_den, aux_J, aux_yJ, err, d_yJ, d_J, y_max
    integer :: ier, neval

    y_max = x/(1 + x)
    call qags(int_den, 0.0_dp, y_max, epsabs, epsrel, aux_den, err, neval, ier)
    call qags(int_J, 0.0_dp, y_max, epsabs, epsrel, aux_J, err, neval, ier)
    call qags(int_yJ, 0.0_dp, y_max, epsabs, epsrel, aux_yJ, err, neval, ier)
    d_yJ = x/(1 + x)**3*JC(x, y_max) + aux_yJ
    d_J = 1/(1 + x)**2*JC(x, y_max) + aux_J
    dM1_dx = d_yJ/aux_den - M1(x)/aux_den*d_J
  contains
    function int_J(y)
      implicit none
      real(dp) :: int_J, y

      int_J = dJ_dx(x, y)
    end function int_J

    function int_yJ(y)
      implicit none
      real(dp) :: int_yJ, y

      int_yJ = y*dJ_dx(x, y)
    end function int_yJ

    function dJ_dx(x, y)
      ! diff(J_C(x, y), x)
      implicit none
      real(dp) :: dJ_dx, x, y

      dJ_dx = -y/(x**2*(1 - y))*(2*log(y/(x*(1 - y))) - (x*y**2 + 6*x*y + 8*y &
        - 6*x)/(2*x*(1 - y)))
    end function dJ_dx

    function int_den(y)
      implicit none
      real(dp) :: int_den, y

      int_den = JC(x, y)
    end function int_den

    function JC(x, y)
      ! cf. eq. (25, 26)
      implicit none
      real(dp) :: JC, x, y, w, xw
 
      w = y/(x*(1 - y))
      xw = y/(1 - y)
      JC = 2*w*log(w) + (1 + 2*w)*(1 - w) + xw**2/(2*(1 + xw))*(1 - w)
    end function JC
  end function dM1_dx

  function M3(x)
    implicit none
    real(dp) :: M3, x

    M3 = M0(x)*M2(x)
  end function M3
end module compton
