! =============================================================================
! test_which_mode.F90
!
! Unit test for subroutine `which_mode` (module `which_mode_to_use`,
! src/which_mode_to_use.F90), which partitions an evaporated-aerosol
! mass/number increment between two lognormal modes.
!
! With the repo's default `mphys_switches::l_aeroproc_midway = .FALSE.`,
! `imethod` is a compile-time parameter fixed to `isimple_method`, so the
! only logic actually exercised at runtime is:
!
!   rm = (rftpi*dm/dn/density)^(1/3)      [rftpi = 1/(4/3*pi)]
!   if (rm >= r_thresh_fixed=0.5e-6) then all of dm,dn -> mode 2
!   else                                   all of dm,dn -> mode 1
!
! Reference dm values below are chosen (in Python, replicating the same
! formula) so rm lands cleanly either side of the 0.5 micron threshold.
! r1_in, r2_in and sigma are unused on this code path but must be supplied.
! =============================================================================
program test_which_mode

  use variable_precision, only: wp
  use which_mode_to_use,  only: which_mode
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: dm, dn, density, sigma, r1_in, r2_in
  real(wp) :: dm1, dm2, dn1, dn2
  real(wp) :: tol

  density = 1777.0_wp
  dn = 1.0e8_wp
  r1_in = 0.1e-6_wp
  r2_in = 2.0e-6_wp
  sigma = 1.5_wp

  ! Case A: rm = 0.3 micron < 0.5 micron threshold -> all mass/number to mode 1
  dm = 2.009739652354462e-08_wp
  call which_mode(dm, dn, r1_in, r2_in, density, sigma, dm1, dm2, dn1, dn2)
  tol = 1.0e-6_wp * abs(dm)
  call check_close('which_mode small-rm dm1', dm1, dm, tol)
  call check_close('which_mode small-rm dm2', dm2, 0.0_wp, 1.0e-20_wp)
  call check_close('which_mode small-rm dn1', dn1, dn, 1.0e-3_wp*dn)
  call check_close('which_mode small-rm dn2', dn2, 0.0_wp, 1.0e-6_wp)

  ! Case B: rm = 0.8 micron >= 0.5 micron threshold -> all mass/number to mode 2
  dm = 3.811061859279573e-07_wp
  call which_mode(dm, dn, r1_in, r2_in, density, sigma, dm1, dm2, dn1, dn2)
  tol = 1.0e-6_wp * abs(dm)
  call check_close('which_mode large-rm dm2', dm2, dm, tol)
  call check_close('which_mode large-rm dm1', dm1, 0.0_wp, 1.0e-20_wp)
  call check_close('which_mode large-rm dn2', dn2, dn, 1.0e-3_wp*dn)
  call check_close('which_mode large-rm dn1', dn1, 0.0_wp, 1.0e-6_wp)

  ! Case C: dm*dn <= 0 -> all outputs left at zero
  dm = -1.0e-8_wp
  call which_mode(dm, dn, r1_in, r2_in, density, sigma, dm1, dm2, dn1, dn2)
  call check_close('which_mode dm*dn<=0 dm1', dm1, 0.0_wp, 1.0e-20_wp)
  call check_close('which_mode dm*dn<=0 dm2', dm2, 0.0_wp, 1.0e-20_wp)
  call check_close('which_mode dm*dn<=0 dn1', dn1, 0.0_wp, 1.0e-20_wp)
  call check_close('which_mode dm*dn<=0 dn2', dn2, 0.0_wp, 1.0e-20_wp)

  call test_summary()

end program test_which_mode
