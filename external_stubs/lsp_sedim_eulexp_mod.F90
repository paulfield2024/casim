! Minimal stand-in for the UM large-scale-precipitation Eulerian/exponential
! sedimentation routine (lsp_sedim_eulexp_mod). The real UM implementation
! cannot be redistributed with CASIM. This routine is only ever called by
! CASIM's sedimentation.F90 when `casim_parent == parent_um`; for standalone
! unit tests and the 1D column model we run with casim_parent == parent_monc,
! so this stub is never actually executed - it exists purely so the code
! links. It implements a simple upstream (donor-cell) explicit fallout so
! that, in case it were ever exercised, results remain physically sane.
module lsp_sedim_eulexp_mod
  use um_types, only: real_lsprec
  implicit none

contains

  subroutine lsp_sedim_eulexp(points, m0, dhi, dhir, rhoin, rhor,           &
       flux_fromabove, fallspeed_thislayer, mixratio_thislayer,            &
       fallspeed_fromabove, total_flux_out)

    integer, intent(in) :: points
    real(real_lsprec), intent(in) :: m0
    real(real_lsprec), intent(in) :: dhi(points), dhir(points)
    real(real_lsprec), intent(in) :: rhoin(points), rhor(points)
    real(real_lsprec), intent(in) :: flux_fromabove(points)
    real(real_lsprec), intent(inout) :: fallspeed_thislayer(points)
    real(real_lsprec), intent(in) :: mixratio_thislayer(points)
    real(real_lsprec), intent(in) :: fallspeed_fromabove(points)
    real(real_lsprec), intent(out) :: total_flux_out(points)

    integer :: i

    do i = 1, points
      total_flux_out(i) = flux_fromabove(i) +                              &
           mixratio_thislayer(i) * rhoin(i) * fallspeed_thislayer(i) * dhir(i)
    end do

  end subroutine lsp_sedim_eulexp

end module lsp_sedim_eulexp_mod
