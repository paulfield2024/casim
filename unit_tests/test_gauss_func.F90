! =============================================================================
! test_gauss_func.F90
!
! Unit test for `gauss_casim_func(a, b)` in module `gauss_casim_micro`
! (src/gauss_4A_func.F90) - the brute-force double numerical integral used
! in the snow-aggregation collision kernel. Reference values were computed
! independently in Python replicating the same rectangle-rule integration
! (tmax=18, dt=0.08), so this is a regression check on the Fortran
! implementation matching that exact algorithm.
! =============================================================================
program test_gauss_func

  use variable_precision, only: wp
  use gauss_casim_micro,  only: gauss_casim_func
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol

  tol = 1.0e-3_wp
  call check_close('gauss_casim_func(0.0,1.0)', &
       gauss_casim_func(0.0_wp, 1.0_wp), 11.996579679342574_wp, tol)

  call check_close('gauss_casim_func(1.0,2.0)', &
       gauss_casim_func(1.0_wp, 2.0_wp), 314.9546636807022_wp, tol*100.0_wp)

  call check_close('gauss_casim_func(0.5,0.5)', &
       gauss_casim_func(0.5_wp, 0.5_wp), 6.6421780915186845_wp, tol)

  call test_summary()

end program test_gauss_func
