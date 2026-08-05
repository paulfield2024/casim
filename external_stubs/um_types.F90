! Minimal stand-in for the UM "um_types" module.
! CASIM's sedimentation.F90 only needs the `real_lsprec` kind, used to
! declare variables passed to the (also stubbed) lsp_sedim_eulexp routine.
module um_types
  implicit none
  integer, parameter :: real_lsprec = selected_real_kind(13, 300)
end module um_types
