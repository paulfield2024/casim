! Minimal stand-in for the UM/ECMWF "yomhook" DrHook profiling module.
! `lhook` is kept permanently false so `dr_hook` is never actually invoked
! by CASIM, and dr_hook itself is a no-op provided purely to satisfy the
! link step.
module yomhook
  use parkind1, only: jprb, jpim
  implicit none
  logical, parameter :: lhook = .false.

contains

  subroutine dr_hook(name, code, handle)
    character(len=*), intent(in) :: name
    integer(kind=jpim), intent(in) :: code
    real(kind=jprb), intent(inout) :: handle
    ! No-op: profiling is disabled (lhook = .false.)
  end subroutine dr_hook

end module yomhook
