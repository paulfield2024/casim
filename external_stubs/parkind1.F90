! Minimal stand-in for the UM/ECMWF "parkind1" kind-definition module.
! CASIM only uses jprb (working real precision) and jpim (working integer
! precision) for its DrHook instrumentation calls, so only those two kinds
! are provided here.
module parkind1
  implicit none
  integer, parameter :: jprb = selected_real_kind(13, 300)
  integer, parameter :: jpim = selected_int_kind(9)
end module parkind1
