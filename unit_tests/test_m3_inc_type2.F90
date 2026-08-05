! =============================================================================
! test_m3_inc_type2.F90
!
! Unit test for subroutine `m3_inc_type2` (module `m3_incs`, src/m3_incs.F90):
! tendency of a 3rd moment given pre-existing moments m1,m2,m3 and increments
! dm1,dm2, assuming the shape parameter mu does not vary.
!
! Two branches are exercised:
!   - m1 > 0: main algebraic branch. With p1=3, p2=0, p3=6 the code takes
!     the hardwired shortcut fac=((m1+dm1)/m1)^2*(m2/(m2+dm2))-1, dm3=m3*fac
!     (reference value computed in Python replicating this exact formula).
!   - m1 <= 0: falls back to m3_inc_type3(p1,p2,p3,dm1,dm2,dm3,mu_init).
!     Reusing the p1=0,p2=3,p3=6,mu=0 case already verified in
!     test_m3_incs.F90 (Gfunc(0,0,3,6)=20 -> dm3=90.0) confirms the dispatch.
! =============================================================================
program test_m3_inc_type2

  use variable_precision, only: wp
  use m3_incs,            only: m3_inc_type2
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: dm3, tol

  ! m1>0 branch, p3=6 hardwired shortcut:
  ! m1=1, m2=1, m3=2, dm1=0.1, dm2=-0.05
  ! fac = ((1.1)/1)^2 * (1/(1-0.05)) - 1 = 0.273684210526316
  ! dm3 = m3*fac = 0.547368421052632
  tol = 1.0e-9_wp
  call m3_inc_type2(1.0_wp, 1.0_wp, 2.0_wp, 3.0_wp, 0.0_wp, 6.0_wp, &
       0.1_wp, -0.05_wp, dm3)
  call check_close('m3_inc_type2 (m1>0, p3=6 shortcut)', dm3, &
       0.547368421052632_wp, tol)

  ! m1<=0 branch: falls back to m3_inc_type3(p1,p2,p3,dm1,dm2,dm3,mu_init)
  ! p1=0,p2=3,p3=6,mu=0,dm1=2,dm2=3 -> dm3 = Gfunc(0,0,3,6)*dm1^-1*dm2^2 = 90
  tol = 2.0e-3_wp * 90.0_wp
  call m3_inc_type2(0.0_wp, 1.0_wp, 1.0_wp, 0.0_wp, 3.0_wp, 6.0_wp, &
       2.0_wp, 3.0_wp, dm3, mu_init=0.0_wp)
  call check_close('m3_inc_type2 (m1<=0 fallback to type3)', dm3, 90.0_wp, tol)

  call test_summary()

end program test_m3_inc_type2
