! =============================================================================
! test_erf.F90
!
! Unit tests for the error-function family in module `special`
! (src/special.F90): `casim_erfc` (complementary error function) and the
! `erfinv` generic interface (currently bound only to `erfinv1`, a truncated
! series approximation).
! =============================================================================
program test_erf

  use variable_precision, only: wp
  use special,            only: casim_erfc, erfinv
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol

  ! --- casim_erfc(x) = erfc(x) -----------------------------------------
  ! Reference values from the standard complementary error function.
  tol = 1.0e-6_wp
  call check_close('casim_erfc(0.01)', casim_erfc(0.01_wp), 0.9887166_wp, tol)
  call check_close('casim_erfc(0.5)',  casim_erfc(0.5_wp),   0.4795001_wp,    tol)
  call check_close('casim_erfc(1.0)',  casim_erfc(1.0_wp),   0.1572992_wp,    tol)
  call check_close('casim_erfc(2.0)',  casim_erfc(2.0_wp),   0.0046777_wp,    1.0e-5_wp)
  call check_close('casim_erfc(-1.0)', casim_erfc(-1.0_wp),  1.8427008_wp,    tol)

  ! Near-zero regression: erfg(x,c) previously had two bugs here (now fixed
  ! in src/special.F90): a 0/0 -> NaN at x==0 exactly, and a branch shared
  ! between erf/erfc that always returned 0.0 for |x| < 1e-10 (wrong for
  ! erfc, which should tend to 1.0).
  call check_close('casim_erfc(0.0)',    casim_erfc(0.0_wp),    1.0_wp, tol)
  call check_close('casim_erfc(1e-12)',  casim_erfc(1.0e-12_wp), 1.0_wp, tol)

  ! --- erfinv(x) = inverse error function --------------------------------
  ! erfinv1 is a truncated (11th order) series, documented in-source as
  ! "needs more work to get good accuracy" - it diverges from the true
  ! inverse erf function as |x| grows. Expected values below are the exact
  ! result of that same truncated series (independently evaluated), so this
  ! is a regression check on the series implementation, not a check against
  ! the true inverse error function.
  tol = 1.0e-6_wp
  call check_close('erfinv(0.0)',        erfinv(0.0_wp),        0.0_wp,               tol)
  call check_close('erfinv(0.5)',        erfinv(0.5_wp),        0.4769296236996508_wp, tol)
  call check_close('erfinv(-0.5)',       erfinv(-0.5_wp),      -0.4769296236996508_wp, tol)
  call check_close('erfinv(0.8427008)',  erfinv(0.8427008_wp),  0.9877954710770509_wp, tol)

  call test_summary()

end program test_erf
