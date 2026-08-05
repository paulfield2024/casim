! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Ventilation factor calculation for depositional/evaporative growth (ventilation_3M, ventilation_1M_2M).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Ventilation factor for diffusional growth/evaporation of
!   falling hydrometeors, integrated over the particle size
!   distribution.
!
! Paper reference:
!   Field et al. (2023) Appendix A.8.1, Eq. (A23)-(A24): F is
!   the per-particle ventilation coefficient; the integrated
!   ventilation factor chi_x combines capacitance, F and the PSD.
!
MODULE ventfac
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho
  USE special, ONLY: pi, GammaFunc
  USE mphys_parameters, ONLY: hydro_params, vent_1, vent_2
  USE mphys_constants, ONLY: visair, rho0, Dv

!prf
  USE mphys_switches, ONLY: l_kfsm
!prf


  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='VENTFAC'

  PUBLIC ventilation_3M, ventilation_1M_2M
CONTAINS

  ! Integrated ventilation factor for a triple-moment PSD,
  ! Field et al. (2023) Eq. (A23)-(A24).
  SUBROUTINE ventilation_3M(ixy_inner, k, V, n0, lam, mu, params)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='VENTILATION_3M'

    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: k
    REAL(wp), INTENT(IN) :: n0, lam, mu
    REAL(wp), INTENT(OUT) :: V  ! bulk ventilation factor
    TYPE(hydro_params), INTENT(IN) :: params

    REAL(wp) :: T1, T2
    REAL(wp) :: Sc ! Schmidt number
    REAL(wp) :: a_x, b_x, f_x

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x

    Sc=visair/Dv

    T1=vent_1*(mu+1.0)/lam
    T2=vent_2*Sc**(1.0/3.0)*(a_x*(rho0/rho(k,ixy_inner))**(0.5)*rho(k,ixy_inner)/visair)**(0.5)

!changing to rho on top give V units of 1/m2 which is correct to combine with melting, wet growth that have units (m2/s *V) to give 1/s for mass change rate
    V=2.0*pi*n0*rho(k,ixy_inner)*(T1+T2*GammaFunc(0.5*b_x+mu+2.5)/GammaFunc(1.0+mu) &
         *(1.0 + 0.5*f_x/lam)**(-(0.5*b_x + mu + 2.5))*lam**(-0.5*b_x-1.5))
!prf for single moment use capacitance of 0.5*sphere -i.e. ~plate or disk or aggregate
    IF (l_kfsm) V=0.5*V
!prf

    ! for computational efficiency/accuracy changed from...
    ![' ']    !*(lam + 0.5*f_x)**(-(0.5*b_x + mu + 2.5))*lam**(1.+mu)) &

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ventilation_3M

  ! Integrated ventilation factor for single-/double-moment
  ! PSDs, Field et al. (2023) Eq. (A23)-(A24).
  SUBROUTINE ventilation_1M_2M(ixy_inner, k, V, n0, lam, mu, params)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='VENTILATION_1M_2M'

    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: k
    REAL(wp), INTENT(IN) :: n0, lam, mu
    REAL(wp), INTENT(OUT) :: V  ! bulk ventilation factor
    TYPE(hydro_params), INTENT(IN) :: params

    REAL(wp) :: T1, T2
    REAL(wp) :: Sc ! Schmidt number
    REAL(wp) :: a_x, b_x, f_x

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    a_x=params%a_x
    b_x=params%b_x
    f_x=params%f_x

    Sc=visair/Dv

    T1=vent_1*(mu+1.0)/lam
    T2=vent_2*Sc**(1.0/3.0)*SQRT(a_x*SQRT(rho0*rho(k,ixy_inner))/visair)

    V=2.0*pi*n0*rho(k,ixy_inner)*(T1+T2*params%gam_0p5bx_mu_2p5/params%gam_1_mu &
         *(1.0 + 0.5*f_x/lam)**(-(0.5*b_x + mu + 2.5))*lam**(-0.5*b_x-1.5))
!prf for single moment use capacitance of 0.5*sphere -i.e. ~plate or disk or aggregate
    IF (l_kfsm) V=0.5*V
!prf

    ! for computational efficiency/accuracy changed from...
    ![' ']    !*(lam + 0.5*f_x)**(-(0.5*b_x + mu + 2.5))*lam**(1.+mu)) &

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ventilation_1M_2M
END MODULE ventfac
