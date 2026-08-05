! =============================================================================
! test_lookup.F90
!
! Unit tests for the moment/gamma-combination pure functions in module
! `lookup` (src/lookup.F90): `Gfunc` and `moment`. These take plain real(wp)
! arguments (no derived types), so are usable as-is in isolation.
! (`Hfunc` also exists in this module but is not in its public list, so it
! is not directly testable from outside the module.)
!
! Reference values are chosen so p1=0, p2=3, p3=6 (mass/number/6th-moment
! powers for spherical particles) giving k1=1, k2=-2, which combined with
! mu=0 makes the expected Gamma-function arguments land on small integers
! (Gamma(1)=1, Gamma(4)=6, Gamma(7)=720) for easy hand-verification.
!
! Also tests the subroutines `get_n0` (plain real(wp) args, no derived
! type needed) and `get_mu` (needs a (mu, Gfunc-value) lookup table, built
! here locally via the public `set_mu_lookup` helper rather than depending
! on module-level state).
! =============================================================================
program test_lookup

  use variable_precision, only: wp
  use lookup,             only: Gfunc, moment, get_n0, get_mu, set_mu_lookup, nmu
  use test_utils,         only: check_close, test_summary

  implicit none

  real(wp) :: tol
  real(wp) :: n0, mu
  real(wp), allocatable :: mu_i_arr(:), mu_g_arr(:)

  ! Gfunc(mu, p1, p2, p3) with mu=0, p1=0, p2=3, p3=6:
  !   k1=(p3-p2)/(p2-p1)=1, k2=(p1-p3)/(p2-p1)=-2
  !   Gfunc = Gamma(1)^1 * Gamma(4)^-2 * Gamma(7) = 1 * (1/36) * 720 = 20
  tol = 2.0e-3_wp * 20.0_wp
  call check_close('Gfunc(0,0,3,6)', Gfunc(0.0_wp, 0.0_wp, 3.0_wp, 6.0_wp), 20.0_wp, tol)

  ! moment(n0, lam, mu, p) = n0 * Gamma(1+mu+p) * lam^(-p) / Gamma(1+mu)
  ! n0=1e6, lam=1000, mu=0, p=3 -> 1e6 * Gamma(4) * 1000^-3 / Gamma(1)
  !                              = 1e6 * 6 * 1e-9 / 1 = 6e-3
  tol = 2.0e-3_wp * 6.0e-3_wp
  call check_close('moment(1e6,1000,0,3)', &
       moment(1.0e6_wp, 1000.0_wp, 0.0_wp, 3.0_wp), 6.0e-3_wp, tol)

  ! get_n0(m, p, mu, lam, n0) = m/(Gamma(1+mu+p)*lam^-p/Gamma(1+mu))
  ! mu=0, p=3, lam=1, m=6 -> 6/(Gamma(4)*1^-3/Gamma(1)) = 6/(6*1/1) = 1.0
  tol = 2.0e-3_wp * 1.0_wp
  call get_n0(6.0_wp, 3.0_wp, 0.0_wp, 1.0_wp, n0)
  call check_close('get_n0(m=6,p=3,mu=0,lam=1)', n0, 1.0_wp, tol)

  ! get_mu(m1,m2,m3,p1,p2,p3,mu,mu_g_o,mu_i_o): searches a (mu -> Gfunc)
  ! lookup table for where Hfunc(m1,m2,m3,p1,p2,p3) falls, and interpolates
  ! back to mu. Build the table locally (min_mu=0, max_mu=35, nmu=501
  ! points) via the public set_mu_lookup helper for p1=0,p2=3,p3=6, then
  ! pick m1=1, m2=1, m3=Gfunc(2.5,0,3,6) so Hfunc(m1,m2,m3,..)=m3 exactly
  ! equals the target Gfunc value, i.e. the table should recover mu~=2.5.
  allocate(mu_i_arr(nmu), mu_g_arr(nmu))
  call set_mu_lookup(0.0_wp, 3.0_wp, 6.0_wp, mu_i_arr, mu_g_arr)
  call get_mu(1.0_wp, 1.0_wp, 4.783549783549783_wp, 0.0_wp, 3.0_wp, 6.0_wp, &
       mu, mu_g_o=mu_g_arr, mu_i_o=mu_i_arr)
  ! Table spacing is (35-0)/500 = 0.07, so allow a bit more than that for
  ! linear-interpolation error.
  call check_close('get_mu(->2.5)', mu, 2.5_wp, 0.1_wp)

  call test_summary()

end program test_lookup
