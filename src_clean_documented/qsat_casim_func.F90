! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Saturation mixing ratio (Tetens formula) over liquid/ice and its temperature derivative.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Saturation mixing ratio with respect to liquid water
!   (Qsaturation) and ice (Qisaturation), used by evaporation,
!   deposition/sublimation and melting rate calculations.
!
! Paper reference:
!   Field et al. (2023) Appendix A.8.1-A.8.4 use qsat/qisat
!   (Eq. A21-A22, A25-A26, A28-A29) computed by this module; the
!   saturation-vapour-pressure formula itself is not spelled out
!   in the paper.
!
MODULE qsat_funs
  USE variable_precision, ONLY: wp
  IMPLICIT NONE
  PRIVATE

  PUBLIC Qsaturation, Qisaturation, dqwsatdt

CONTAINS
  ! Function to return the saturation mr over water
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  ! Saturation mixing ratio w.r.t. liquid water, used as
  ! qwsat in Field et al. (2023) Eq. (A25)-(A26).
  FUNCTION Qsaturation (T, p)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: T, p
    REAL(wp) :: Qsaturation
    ! Temperature in Kelvin
    ! Pressure in mb
    REAL(wp), PARAMETER ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation
    !

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    IF (T > qsa3 .AND. p * exp (qsa2 * (t - tk0c) / (T - qsa3)) > qsa4) THEN
      Qsaturation=qsa1/(p*exp(qsa2*(t-tk0c)/(T-qsa3))-qsa4)
    ELSE
      qsaturation=999.0
    END IF

  END FUNCTION Qsaturation

  ! Function to return the saturation mr over ice
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-21.8745584*(T-273.15)/(T-7.66))-6.109)
  ! Saturation mixing ratio w.r.t. ice, used as qisat in
  ! Field et al. (2023) Eq. (A28)-(A29), (A38).
  FUNCTION Qisaturation(T, p)

    IMPLICIT NONE

    REAL(wp), INTENT(IN) ::  T, p
    REAL(wp) :: qisaturation
    ! Temperature in Kelvin
    ! Pressure in mb
    REAL(wp), PARAMETER :: tk0c = 273.15, qis1 = 3.8, qis2 = -21.8745584 , qis3 = 7.66  , qis4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qisat equation
    ! Constant in qisat equation
    ! Constant in qisat equation

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    IF (T > qis3 .AND. p * exp (qis2 * (t - tk0c) / (T - qis3)) > qis4) THEN
       qisaturation = qis1/(p*exp(qis2*(T - tk0c)/(T - qis3)) - qis4)
    ELSE
       qisaturation = 999.0
    END IF

  END FUNCTION Qisaturation

  ! Function to return the rate of change with temperature
  ! of saturation mixing ratio over liquid water.
  !
  ! Based on tetans formular
  ! QS=3.8/(P*EXP(-17.2693882*(T-273.15)/(T-35.86))-6.109)
  FUNCTION dqwsatdt (qsat, T)

    IMPLICIT NONE

    REAL(wp) , INTENT(IN) ::qsat, T
    REAL(wp) ::dqwsatdt
    ! Saturatio mixing ratio
    ! Temperature in Kelvin

    REAL(wp), PARAMETER ::tk0c = 273.15, qsa1 = 3.8, qsa2 = - 17.2693882, qsa3 = 35.86, qsa4 = 6.109
    ! Temperature of freezing in Kelvin
    ! Top in equation to calculate qsat
    ! Constant in qsat equation
    ! Constant in qsat equation
    ! Constant in qsat equation

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------

    dqwsatdt=-qsa2*(TK0C-qsa3)*(1.0+qsa4*qsat/qsa1)*qsat*(T-qsa3)**(- 2.0)

  END FUNCTION dqwsatdt
END MODULE qsat_funs
