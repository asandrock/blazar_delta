program test_compton
  ! Reproduce figg. 2-7
  use compton
  implicit none
  real(dp) :: x
  integer :: j

  open(unit=23, file='fig-2.dat', status='replace')
  open(unit=24, file='fig-3.dat', status='replace')
  open(unit=25, file='fig-4.dat', status='replace')
  open(unit=26, file='fig-5.dat', status='replace')
  open(unit=27, file='fig-6.dat', status='replace')
  open(unit=28, file='fig-7.dat', status='replace')
  do j = -500, 500
    x = 10**(j/100.0_dp)
    ! Fig. 2
    write (23,*) x, S0(x), 1 - 2*x + (26*x**2)/5 - (133*x**3)/10 &
      + (1444*x**4)/35, 3/(8*x)*(log(x) + 0.5_dp*(1 + log(4.0_dp)))
    ! Fig. 3
    write (24,*) x, M0(x)/x**2, (x**2/3 - (2*x**3)/9 + (13*x**4)/60 &
      - (133*x**5)/600)/x**2, (x/2*(log(x) - 0.5_dp) + 3.114_dp*log(x) &
      - log(x/2)**2 - 6.559_dp + 2/x*(2*log(x) + 3.249_dp))/x**2
    ! Fig. 4
    write (25,*) x, S1(x), x, 1 - 4/(3*(log(2*x) + 0.5_dp))
    ! Fig. 5
    write (26,*) x, M1(x), x/3, 0.6279_dp
    ! Fig. 6
    write (27,*) x, S2(x)
    ! Fig. 7
    write (28,*) x, M2(x)
    print *, x, S3(x), M3(x)
  end do
  close(unit=23, status='keep')
  close(unit=24, status='keep')
  close(unit=25, status='keep')
  close(unit=26, status='keep')
  close(unit=27, status='keep')
  close(unit=28, status='keep')
end program test_compton
