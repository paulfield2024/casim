! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Flags which columns/levels have any microphysics work to do, to skip empty ones (preconditioner).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Decides which grid columns/levels actually need the full
!   microphysics calculation (i.e. contain hydrometeors or are
!   sub-saturated), to skip unnecessary work elsewhere.
!
! Paper reference:
!   CASIM software/performance infrastructure; not discussed in
!   Field et al. (2023).
!
MODULE Preconditioning
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: qws
  USE mphys_switches, ONLY: i_qv, i_ql, i_qr, i_qi, i_qs, i_qg , cloud_params, &
       rain_params, ice_params, snow_params, graupel_params, l_cfrac_casim_diag_scheme
  USE thresholds, ONLY: thresh_tidy
  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='PRECONDITIONING'

  LOGICAL, ALLOCATABLE :: precondition(:,:)

!$OMP THREADPRIVATE(precondition)

  PUBLIC precondition, preconditioner
CONTAINS

  SUBROUTINE preconditioner(ixy_inner, qfields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='PRECONDITIONER'

    INTEGER, INTENT(IN) :: ixy_inner

    REAL(wp), INTENT(IN) :: qfields(:,:)

    INTEGER :: k
    LOGICAL :: l_temp

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k=ubound(precondition,1)-1, 1, -1
      ! Do we have any existing hydrometeor mass?
      l_temp=.FALSE.
      IF (cloud_params%l_1m) l_temp=l_temp .OR. qfields(k, i_ql) > thresh_tidy(i_ql)
      IF (rain_params%l_1m) l_temp=l_temp .OR. qfields(k, i_qr) > thresh_tidy(i_qr)
      IF (ice_params%l_1m) l_temp=l_temp .OR. qfields(k, i_qi) > thresh_tidy(i_qi)
      IF (snow_params%l_1m) l_temp=l_temp .OR. qfields(k, i_qs) > thresh_tidy(i_qs)
      IF (graupel_params%l_1m) l_temp=l_temp .OR. qfields(k, i_qg) > thresh_tidy(i_qg)
      ! Do we have supersaturation
      l_temp=l_temp .OR. qws(k,ixy_inner) < qfields(k, i_qv)
      ! Do we meet heterogeneous freezing condtion
      ! Need to add this for ice phase...
      ! l_temp = l_temp .or. Si > 0.25
      ! Do we have something above which might fall down
      l_temp=l_temp .OR. precondition(k+1,ixy_inner)
      l_temp = l_temp .OR. l_cfrac_casim_diag_scheme !DPG - To prevent early quit from
        !CASIM if we are sub-saturated since want to allow cloud scheme to operate even
        !if we are subsaturated.
      ! qsat doesn't work at very low pressures,
      ! so if qsaturation is 0.0 then don't do microphysics
      IF (qws(k,ixy_inner) <= 1.0e-6) l_temp=.FALSE.
      IF (qfields(k,i_qv) <= 3.0e-6) l_temp=.FALSE.

      ! OK, that's all...
      precondition(k,ixy_inner)=l_temp
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE preconditioner
END MODULE Preconditioning
