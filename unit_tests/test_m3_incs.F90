! =============================================================================
! test_m3_incs.F90
!
! Unit tests for the algebraic 3rd-moment increment subroutines in module
! `m3_incs` (src/m3_incs.F90): `m3_inc_type1`, `m3_inc_type3`, `m3_inc_type4`.
! (`m3_inc_type2` is not tested directly here - it dispatches to
! `m3_inc_type3` as a fallback and otherwise needs realistic pre-existing
! moment triples m1,m2,m3 to exercise its main branch meaningfully.)
! =============================================================================
program test_m3_incs

  use variable_precision, only: wp
  use m3_incs,            only: m3_inc_type1, m3_inc_type3, m3_inc_type4
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol, dm2, dm3, dm3_x

  ! m3_inc_type1(m1, m2, dm1, dm2): dm2 = dm1 * m1 / m2
  tol = 1.0e-10_wp
  call m3_inc_type1(1.0e-3_wp, 1.0e6_wp, 1.0e-5_wp, dm2)
  call check_close('m3_inc_type1', dm2, 1.0e-14_wp, tol)

  ! m3_inc_type4(dm3_y, c_x, c_y, p3, dm3_x): dm3_x = -(c_y/c_x)^(p3/3) * dm3_y
  ! c_x=1000, c_y=900, p3=3 -> -(900/1000)^1 * dm3_y = -0.9 * dm3_y
  call m3_inc_type4(2.0e-10_wp, 1000.0_wp, 900.0_wp, 3.0_wp, dm3_x)
  call check_close('m3_inc_type4', dm3_x, -1.8e-10_wp, tol)

  ! m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu): with p1=0,p2=3,p3=6,mu=0
  ! (see test_lookup.F90): k1=(p3-p2)/(p2-p1)=1, k2=(p1-p3)/(p2-p1)=-2
  ! dm3 = Gfunc(0,0,3,6) * dm1^-1 * dm2^-(-2) = 20 * dm1^-1 * dm2^2
  ! dm1=2.0, dm2=3.0 -> dm3 = 20 * 0.5 * 9.0 = 90.0
  tol = 2.0e-3_wp * 90.0_wp
  call m3_inc_type3(0.0_wp, 3.0_wp, 6.0_wp, 2.0_wp, 3.0_wp, dm3, 0.0_wp)
  call check_close('m3_inc_type3', dm3, 90.0_wp, tol)

  call test_summary()

end program test_m3_incs
