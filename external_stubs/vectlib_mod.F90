! Minimal stand-in for the UM/vector-math-library "vectlib_mod" module.
! CASIM's casim_moments_mod.F90 (only exercised when l_kfsm = .true., which
! is not the default and is not used by the standalone unit tests or 1D
! column model) needs simple vectorised elementwise exp/reciprocal/power.
! These plain elementwise implementations are numerically equivalent to the
! optimised UM versions.
module vectlib_mod
  use variable_precision, only: wp
  implicit none

contains

  subroutine exp_v(n, x, y)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(out) :: y(n)
    y(1:n) = exp(x(1:n))
  end subroutine exp_v

  subroutine oneover_v(n, x, y)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(out) :: y(n)
    y(1:n) = 1.0_wp / x(1:n)
  end subroutine oneover_v

  subroutine powr_v(n, x, p, y)
    integer, intent(in) :: n
    real(wp), intent(in) :: x(n)
    real(wp), intent(in) :: p
    real(wp), intent(out) :: y(n)
    y(1:n) = x(1:n) ** p
  end subroutine powr_v

end module vectlib_mod
