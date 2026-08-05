! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Numerical double-integral (gauss_casim_func) used in the snow self-aggregation collision kernel.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE gauss_casim_micro

  !< Calculate the integral needed in the aggregation calculations
  !<
  !< OPTIMISATION POTENTIAL - LOOKUP

  USE variable_precision, ONLY: wp
  USE mphys_switches, ONLY: max_mu
  IMPLICIT NONE
  PRIVATE

  INTEGER, PARAMETER :: maxq = 5
  REAL(wp) :: gaussfunc_save(maxq) ! max 5 values - I.e. this is only used with 2m schemes
  ! should be extended
  INTEGER, PARAMETER, PRIVATE :: nbins_a = 50
  REAL(wp) :: gaussfunc_save_2D(nbins_a,maxq)
  LOGICAL :: l_save_2D(nbins_a,maxq) = .FALSE.

!$OMP THREADPRIVATE(gaussfunc_save, gaussfunc_save_2D, l_save_2D)

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='GAUSS_CASIM_MICRO'

  PUBLIC gauss_casim_func, gaussfunclookup, gaussfunclookup_2d
CONTAINS

  SUBROUTINE gaussfunclookup(iq, gauss_value, a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: iq !< parameter index relating to variable we're considering
    REAL(wp), INTENT(OUT) :: gauss_value !< returned value
    REAL(wp), INTENT(IN), OPTIONAL :: a, b  !< Value of a and b to use. Only required if initializing

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='GAUSSFUNCLOOKUP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF ( present(a) .AND. present(b)) THEN
      gauss_value=gauss_casim_func(a,b)
      gaussfunc_save(iq)=gauss_casim_func(a,b)
    ELSE
      gauss_value=gaussfunc_save(iq)
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE gaussfunclookup

  SUBROUTINE gaussfunclookup_2d(iq, gauss_value, a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: iq !< parameter index relating to variable we're considering
    REAL(wp), INTENT(OUT) :: gauss_value !< returned value
    REAL(wp), INTENT(IN) :: a, b  !< Value of a and b to use.  (a is mu, b is b_x)

    ! Local variables
    INTEGER :: ibin ! mu(a) bin in which we sit

    CHARACTER(len=*), PARAMETER :: RoutineName='GAUSSFUNCLOOKUP_2D'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ibin=int((a/max_mu)*(nbins_a-1))+1
    IF (.NOT. l_save_2D(ibin,iq)) THEN
      gauss_value=gauss_casim_func(a,b)
      gaussfunc_save_2d(ibin,iq)=gauss_casim_func(a,b)
      l_save_2D(ibin,iq)=.TRUE.
    ELSE
      gauss_value=gaussfunc_save_2d(ibin,iq)
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE gaussfunclookup_2d

  FUNCTION gauss_casim_func(a, b)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(wp), INTENT(IN) :: a, b !< function arguments

    ! Local variables
    REAL(wp), PARAMETER ::   tmax = 18.0    !< Limit of integration
    REAL(wp), PARAMETER ::   dt   = 0.08   !< step size
    REAL(wp) :: gauss_sum, t1, t2
    REAL(wp) :: Gauss_casim_Func

    CHARACTER(len=*), PARAMETER :: RoutineName='GAUSS_CASIM_FUNC'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    gauss_sum=0.0
    t1=0.5*dt
    Gaus_t1: DO WHILE (t1 <= tmax)
      t2=0.5*dt
      Gaus_t2: DO WHILE (t2 <= tmax)
        gauss_sum=gauss_sum+(t1+t2)**2*abs((t1**b)-(t2**b))*(t1**a)*(t2**a)*exp(-(t1+t2))
        t2=t2+dt
      END DO Gaus_t2
      t1=t1+dt
    END DO Gaus_t1
    Gauss_casim_Func=gauss_sum*dt*dt

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION Gauss_casim_Func
END MODULE gauss_casim_micro
