! =============================================================================
! test_lognormal.F90
!
! Unit tests for `MNtoRm(M, N, density, sigma)` in module `lognormal_funcs`
! (src/lognormal_funcs.F90): inverts lognormal-mode aerosol mass & number to
! a mean radius, with a zero-return guard below thresholds%aeromass_small
! (1e-25 kg/kg) / thresholds%aeronumber_small (1e-6 /kg).
! =============================================================================
program test_lognormal

  use variable_precision, only: wp
  use lognormal_funcs,    only: MNtoRm
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol

  ! M=1e-9 kg/kg, N=1e8 /kg, density=1777 kg/m3 (typical ammonium sulfate),
  ! sigma=1.5 (geometric std dev)
  tol = 1.0e-3_wp * 8.622695084822716e-8_wp
  call check_close('MNtoRm(1e-9,1e8,1777,1.5)', &
       MNtoRm(1.0e-9_wp, 1.0e8_wp, 1777.0_wp, 1.5_wp), 8.622695084822716e-8_wp, tol)

  ! Below the mass threshold (aeromass_small=1e-25) -> returns 0
  tol = 1.0e-30_wp
  call check_close('MNtoRm below mass threshold', &
       MNtoRm(1.0e-30_wp, 1.0e8_wp, 1777.0_wp, 1.5_wp), 0.0_wp, tol)

  ! Below the number threshold (aeronumber_small=1e-6) -> returns 0
  call check_close('MNtoRm below number threshold', &
       MNtoRm(1.0e-9_wp, 1.0e-10_wp, 1777.0_wp, 1.5_wp), 0.0_wp, tol)

  call test_summary()

end program test_lognormal
