! Minimal stand-in for the UM diagnostic (Smith-type) cloud fraction scheme
! module used by condensation.F90. The real scheme lives in the UM and
! cannot be redistributed with CASIM. It is only called when
! `l_cfrac_casim_diag_scheme .and. casim_parent == parent_um`; CASIM's
! standalone unit tests and 1D column model run with
! l_cfrac_casim_diag_scheme = .false. (the CASIM default), so this stub is
! never actually executed - it exists purely so the code links.
module cloud_frac_scheme
  use variable_precision, only: wp
  implicit none

contains

  subroutine cloud_frac_casim_mphys(k, p, T, T_liq, rhcrit, qsat, qt,       &
       cloud_mass, rain_mass, cloud_mass_new)

    integer, intent(in) :: k
    real(wp), intent(in) :: p, T, T_liq, rhcrit, qsat, qt
    real(wp), intent(in) :: cloud_mass, rain_mass
    real(wp), intent(out) :: cloud_mass_new

    ! Fallback all-or-nothing behaviour (not used in practice - see above).
    if (qt > qsat) then
      cloud_mass_new = qt - qsat
    else
      cloud_mass_new = 0.0_wp
    end if

  end subroutine cloud_frac_casim_mphys

end module cloud_frac_scheme
