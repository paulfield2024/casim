! Minimal stand-in for the UM "mphys_radar_mod" module (radar-reflectivity
! constants). CASIM's casim_reflec_mod.F90 defines identical local constants
! of the same name/value for its own use in setup_reflec_constants; this
! stub simply re-exports those same standard values so that casim_reflec's
! `use mphys_radar_mod` resolves for standalone builds outside the UM.
module mphys_radar_mod
  use precision, only: wp
  implicit none
  real(wp), parameter :: kliq  = 0.93_wp
  real(wp), parameter :: kice  = 0.174_wp
  real(wp), parameter :: mm6m3 = 1.0e18_wp
  real(wp), parameter :: ref_lim = -35.0_wp
  real(wp), parameter :: ref_lim_lin = 3.1623e-4_wp
  real(wp), parameter :: mr_lim = 1.0e-8_wp
end module mphys_radar_mod
