! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Error function/complementary error function used by the Shipway activation scheme.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE shipway_erf

  USE variable_precision, ONLY: wp
  USE shipway_constants, ONLY: pi

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SHIPWAY_ERF'

CONTAINS

  FUNCTION erfg(x,c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE
    
    REAL(wp), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: c ! 0 gives erf(x)
                             ! 1 gives erfc(x)
    REAL(wp) :: erfg, f, z
    INTEGER :: j, cc

    REAL(wp) :: t, a1, a2, a3, a4, a5, p
    INTEGER :: sign_x

    CHARACTER(len=*), PARAMETER :: RoutineName='ERFG'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    z=x
    cc=c

    IF (abs(z) < 1.5)THEN 
       j=3 + int(9*abs(z))
       f=1
       DO WHILE(j /= 0)
          f=1.+f*z**2*(.5-j)/(j*(.5+j))
          j=j-1
       END DO
       f=cc+f*z*(2.-4.*cc)/sqrt(pi)
    ELSE
      ! Something wrong with this for x<~-2
      ! so use A&S instead
      ! cc=cc*sign_x(1.,z)
      ! j=3 + int(32/abs(z))
      ! f=0
      ! do while(j /= 0)
      !    f=1./(f*j+sqrt(2*z*z))
      !    j=j-1
      ! end do
      ! f=f*(cc*cc+cc-1.)*sqrt(2./pi)*exp(-z*z) + (1-cc)
      !---------------------
      ! Following from A&S 
      ! save the sign_x of x
      sign_x = 1
      IF (x <= 0) sign_x = -1
      z = abs(x)

      !# constants
      a1 =  0.254829592
      a2 = -0.284496736
      a3 =  1.421413741
      a4 = -1.453152027
      a5 =  1.061405429
      p  =  0.3275911

      ! A&S formula 7.1.26
      t = 1.0/(1.0 + p*x)
      f = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x)
      IF (c==1)f=1-f
    END IF
    
    erfg=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION erfg

END MODULE shipway_erf
