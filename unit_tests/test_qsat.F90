! =============================================================================
! test_qsat.F90
!
! Unit tests for module `qsat_funs` (src/qsat_casim_func.F90): Tetens-formula
! saturation mixing ratio over water/ice, and the Clausius-Clapeyron-derived
! d(qsat)/dT helper.
! =============================================================================
program test_qsat

  use variable_precision, only: wp
  use qsat_funs,          only: Qsaturation, Qisaturation, dqwsatdt
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol, qs

  tol = 1.0e-4_wp

  ! Qsaturation(T, p): T in Kelvin, p in mb.
  call check_close('Qsaturation(273.15, 1013.25)', &
       Qsaturation(273.15_wp, 1013.25_wp), 3.7730566e-3_wp, tol*1.0e-2_wp)

  call check_close('Qsaturation(293.15, 1013.25)', &
       Qsaturation(293.15_wp, 1013.25_wp), 1.4696304e-2_wp, tol*1.0e-1_wp)

  ! Below the valid domain (T <= qsa3=35.86), the function returns the
  ! documented sentinel value of 999.0
  call check_close('Qsaturation(30.0, 1013.25) sentinel', &
       Qsaturation(30.0_wp, 1013.25_wp), 999.0_wp, tol)

  ! Qisaturation(T, p): saturation mixing ratio over ice.
  ! At T=273.15K, the ice and water saturation curves coincide (both use
  ! T=0C as their reference).
  call check_close('Qisaturation(273.15, 1013.25)', &
       Qisaturation(273.15_wp, 1013.25_wp), 3.7730566e-3_wp, tol*1.0e-2_wp)

  ! At colder temperatures, qsat over ice is lower than qsat over water
  ! (this is the basis of the Wegener-Bergeron-Findeisen process).
  call check_close('Qisaturation(253.15, 1013.25) < Qsaturation', &
       Qisaturation(253.15_wp, 1013.25_wp), 6.3175045e-4_wp, tol*1.0e-3_wp)

  ! dqwsatdt(qsat, T): derivative should be positive (qsat increases with T)
  qs = Qsaturation(293.15_wp, 1013.25_wp)
  call check_close('dqwsatdt(qsat(20C), 293.15)', &
       dqwsatdt(qs, 293.15_wp), 9.3123676e-4_wp, 1.0e-6_wp)

  call test_summary()

end program test_qsat
