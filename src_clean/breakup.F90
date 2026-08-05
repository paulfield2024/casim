! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Mechanical breakup of snow aggregates (ice_breakup).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE breakup
  USE variable_precision, ONLY: wp, iwp
  USE process_routines, ONLY: process_rate, process_name, i_sbrk
  USE mphys_parameters, ONLY: hydro_params, DSbrk, tau_sbrk
  USE thresholds, ONLY: thresh_small
  USE distributions, ONLY: dist_lambda, dist_mu

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='BREAKUP'

  PRIVATE

  PUBLIC ice_breakup
CONTAINS
  !< Subroutine to determine the breakup of large particles
  !< This code is specified for just snow, but could be used for
  !< other species (e.g. rain)
  !< For triple moment species there is a corresponding change in the
  !< 3rd moment assuming shape parameter is not changed
  !< NB: Aerosol mass is not modified by this process
  !
  !< OPTIMISATION POSSIBILITIES: strip out shape parameters
  SUBROUTINE ice_breakup(nz, l_Tcold, params, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:)
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! Local variables
    TYPE(process_name) :: iproc ! processes selected depending on which species we're modifying
    REAL(wp) :: dnumber
    REAL(wp) :: num, mass
    REAL(wp) :: lam, mu
    REAL(wp) :: Dm ! Mass-weighted mean diameter
    
    INTEGER :: k
    
    CHARACTER(len=*), PARAMETER :: RoutineName='ICE_BREAKUP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (l_Tcold(k)) THEN

          SELECT CASE (params%id)
          CASE (4_iwp) !snow
             iproc=i_sbrk
          END SELECT

          mass=qfields(k, params%i_1m)

          IF (mass > thresh_small(params%i_1m) .AND. params%l_2m) THEN ! if no existing ice, we don't bother
             num=qfields(k, params%i_2m)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)
             Dm=(1.0 + params%d_x + mu)/lam
             
             IF (Dm > DSbrk) THEN ! Mean size exceeds threshold
                dnumber=(Dm/DSbrk - 1.0)**params%d_x * num / tau_sbrk
                procs(params%i_2m, iproc%id)%column_data(k)=dnumber
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ice_breakup
END MODULE breakup
