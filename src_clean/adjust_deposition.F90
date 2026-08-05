! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Rescales ice/snow/graupel vapour deposition rates to stay physically consistent with the available vapour.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE adjust_deposition
  ! As in Harrington et al. 1995 move some of the depositional
  ! growth on ice into the snow category

  USE variable_precision, ONLY: wp
  USE mphys_parameters, ONLY: DImax, snow_params, ice_params
  USE distributions, ONLY: dist_lambda
  USE process_routines, ONLY: process_rate, i_idep, i_sdep, i_saut

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ADJUST_DEPOSITION'

  PRIVATE

  PUBLIC adjust_dep
CONTAINS

  SUBROUTINE adjust_dep(nz, l_Tcold, procs)
    ! only grow ice which is not autoconverted to
    ! snow, c.f. Harrington et al (1995)
    ! This assumes that mu_ice==0, so the fraction becomes
    ! P(mu+2, lambda*DImax) (see Abramowitz & Stegun 6.5.13)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! local variables

    !type(process_rate), pointer :: ice_dep, snow_dep, ice_aut
    REAL(wp) :: lam, frac, dmass

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='ADJUST_DEP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (l_Tcold(k)) THEN
          IF (procs(ice_params%i_1m, i_saut%id)%column_data(k) > 0) THEN
             lam=dist_lambda(k,ice_params%id)
             frac=1.0-exp(-lam*DImax)*(1.0+lam*DImax)
             dmass=frac*procs(ice_params%i_1m, i_idep%id)%column_data(k)

             procs(ice_params%i_1m, i_idep%id)%column_data(k)= &
                  procs(ice_params%i_1m, i_idep%id)%column_data(k)-dmass
             procs(snow_params%i_1m,i_sdep%id)%column_data(k)= &
                  procs(snow_params%i_1m,i_sdep%id)%column_data(k)+dmass
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE adjust_dep
END MODULE adjust_deposition
