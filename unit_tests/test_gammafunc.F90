! =============================================================================
! test_gammafunc.F90
!
! Pilot unit test proving out the unit_tests/ build+run infrastructure.
! Exercises `Gammafunc` (module `special`, src/special.F90) - a pure
! function computing the Gamma function via a Stirling-series expansion
! with a recursion relation for small arguments - against known exact
! values.
! =============================================================================
program test_gammafunc

  use variable_precision, only: wp
  use special,            only: Gammafunc
  use mphys_constants,    only: pi
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol

  ! Gammafunc uses a recursion relation to raise small arguments above 3,
  ! then applies Stirling's asymptotic series - this gives a small RELATIVE
  ! error (observed ~1e-4) rather than a fixed absolute error, so scale the
  ! tolerance with the expected magnitude.
  tol = 2.0e-4_wp

  ! Gamma(1) = 1
  call check_close('Gammafunc(1.0)', Gammafunc(1.0_wp), 1.0_wp, tol*1.0_wp)

  ! Gamma(2) = 1! = 1
  call check_close('Gammafunc(2.0)', Gammafunc(2.0_wp), 1.0_wp, tol*1.0_wp)

  ! Gamma(3) = 2! = 2
  call check_close('Gammafunc(3.0)', Gammafunc(3.0_wp), 2.0_wp, tol*2.0_wp)

  ! Gamma(4) = 3! = 6
  call check_close('Gammafunc(4.0)', Gammafunc(4.0_wp), 6.0_wp, tol*6.0_wp)

  ! Gamma(5) = 4! = 24
  call check_close('Gammafunc(5.0)', Gammafunc(5.0_wp), 24.0_wp, tol*24.0_wp)

  ! Gamma(0.5) = sqrt(pi)
  call check_close('Gammafunc(0.5)', Gammafunc(0.5_wp), sqrt(pi), tol*sqrt(pi))

  call test_summary()

end program test_gammafunc
