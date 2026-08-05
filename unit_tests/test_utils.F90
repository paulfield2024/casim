! =============================================================================
! test_utils.F90
!
! Minimal shared assertion helpers for CASIM unit tests. Each unit test
! program `use`s this module, calls `check_close`/`check_true` for each
! individual assertion, prints a PASS/FAIL line per assertion (which also
! becomes part of the KGO text baseline), and calls `test_summary` at the
! end to get a non-zero process exit status if anything failed.
! =============================================================================
module test_utils

  use variable_precision, only: wp

  implicit none
  private

  integer :: n_checks = 0
  integer :: n_failed = 0

  public :: check_close, check_true, test_summary

contains

  ! Assert that `actual` is within `tol` of `expected` (absolute tolerance).
  subroutine check_close(name, actual, expected, tol)
    character(len=*), intent(in) :: name
    real(wp), intent(in) :: actual, expected, tol

    n_checks = n_checks + 1
    if (abs(actual - expected) <= tol) then
      write(*, '(A,A,A,ES16.7,A,ES16.7)') 'PASS ', trim(name), &
           ': actual=', actual, ' expected=', expected
    else
      n_failed = n_failed + 1
      write(*, '(A,A,A,ES16.7,A,ES16.7,A,ES16.7)') 'FAIL ', trim(name), &
           ': actual=', actual, ' expected=', expected, ' tol=', tol
    end if
  end subroutine check_close

  ! Assert that a logical condition holds.
  subroutine check_true(name, condition)
    character(len=*), intent(in) :: name
    logical, intent(in) :: condition

    n_checks = n_checks + 1
    if (condition) then
      write(*, '(A,A)') 'PASS ', trim(name)
    else
      n_failed = n_failed + 1
      write(*, '(A,A)') 'FAIL ', trim(name)
    end if
  end subroutine check_true

  ! Print a final summary line and stop with exit code 1 if any check failed.
  subroutine test_summary()
    write(*, '(A,I0,A,I0,A)') 'SUMMARY: ', n_checks - n_failed, '/', n_checks, ' checks passed'
    if (n_failed > 0) then
      stop 1
    end if
  end subroutine test_summary

end module test_utils
