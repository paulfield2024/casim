! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Special mathematical functions: Gamma function (lookup table), error function/complementary error function, and their inverses.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Special mathematical functions (Gamma function, error
!   function/complementary error function and its inverse) used
!   throughout the PSD and ventilation/deposition calculations.
!
! Paper reference:
!   Field et al. (2023) does not describe these generic special-
!   function implementations directly; they underpin the
!   Gamma(.) terms used throughout Appendix A (e.g. Eq. A5,
!   A8-A9, A12-A18, A36-A38, A44).
!
MODULE special
  USE variable_precision, ONLY: wp
  USE mphys_constants, ONLY: pi
  !  Use solvers, only: brent

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SPECIAL'

  REAL(wp), PARAMETER :: euler=0.57721566

  ! pi is set to the same value as that used in the UM 

  INTERFACE erfinv
     MODULE PROCEDURE erfinv1
  END INTERFACE erfinv

  INTERFACE gammafunc
     !     module procedure gammafunc1
     MODULE PROCEDURE gammalookup
     !     module procedure intrinsic_gamma
  END INTERFACE gammafunc

  REAL(wp), ALLOCATABLE :: gammalookup_arg(:)
  REAL(wp), ALLOCATABLE :: gammalookup_val(:)
  REAL(wp) :: gammalookup_xmin, gammalookup_xmax, gammalookup_dx
  LOGICAL :: l_gammalookup_set=.FALSE.
  
  PUBLIC pi, Gammafunc, casim_erfc, erfinv
CONTAINS
  ! NB The following should provide sufficient range
  ! and density of points for linear interpolation
  ! to provide appropriate accuracy for any values required.
  SUBROUTINE set_gammalookup(xmin, xmax, dx)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SET_GAMMALOOKUP'

    REAL(wp), INTENT(IN) :: xmin !< Minimum value of argument
    REAL(wp), INTENT(IN) :: xmax !< Maximum value of argument
    REAL(wp), INTENT(IN) :: dx   !< spacing of argument calculations

    ! Local variables
    REAL(wp) :: arg
    INTEGER :: nargs
    INTEGER :: i

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    gammalookup_xmin=xmin
    gammalookup_xmax=xmax
    gammalookup_dx=dx

    nargs=ceiling((xmax - xmin)/dx + 1)
    ALLOCATE(gammalookup_arg(nargs))
    ALLOCATE(gammalookup_val(nargs))

    arg=xmin
    DO i=1, nargs
      gammalookup_arg(i)=arg
      gammalookup_val(i)=gammaFunc1(arg)
      arg=arg+dx
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE set_gammalookup

  FUNCTION gammalookup(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GAMMALOOKUP'

    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: gammalookup

    REAL(wp) :: xmin=1e-12, xmax=100.0, dx=.0001
    INTEGER :: i_minus

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (.NOT. l_gammalookup_set) THEN
      CALL set_gammalookup(xmin, xmax, dx)
      l_gammalookup_set=.TRUE.
    END IF
    ! Locate x in table
    i_minus=int((x - gammalookup_xmin)/gammalookup_dx)+1
    gammalookup=0.5*(gammalookup_val(i_minus)+gammalookup_val(i_minus+1))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION gammalookup

  !================!
  ! Gamma function !
  !================!
  FUNCTION gammafunc1(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GAMMAFUNC1'

    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: gammafunc1

    REAL(wp) :: f,g,z

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    f=huge(x)
    g=1
    z=x
    IF ((z+int(abs(z))) /= 0) THEN
      DO WHILE(z < 3)
        ! Lets use a recursion relation for Gamma functions
        ! to get a large argument and use Stirlings formula
        g=g*z
        z=z+1
      END DO

      ! This is just stirlings formula...
      f=(1.0-2.0*(1-2.0/(3.0*z*z))/(7.0*z*z))/(30.0*z*z)
      f=(1.0-f)/(12.0*z)+z*(log(z)-1)
      f=(exp(f)/g)*sqrt(2.0*pi/z)
    END IF
    gammafunc1=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION gammafunc1

  FUNCTION erfg(x,c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ERFG'

    REAL(wp), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: c ! 0 gives erf(x)
    ! 1 gives erfc(x)
    REAL(wp) :: erfg, f, z
    INTEGER :: j, cc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    z=x
    cc=c

    IF (abs(z) < 1e-10) THEN
      f=0.0
    ELSE IF (abs(z) < 1.5) THEN
      j=3+int(9*abs(z))
      f=1
      DO WHILE(j /= 0)
        f=1.0+f*z**2*(.5-j)/(j*(.5+j))
        j=j-1
      END DO
      f=cc+f*z*(2.0-4.0*cc)/sqrt(pi)
    ELSE
      cc=cc*int(abs(z)/z)
      j=3+int(32/abs(z))
      f=0.0
      DO WHILE(j /= 0)
        f=1.0/(f*j + sqrt(2.0*z*z))
        j=j-1
      END DO
      f=f*(cc*cc+cc-1.0)*sqrt(2.0/pi)*exp(-z*z)+(1.0-cc)
    END IF

    ! quick fix, but should do this properly...
    f=f*((1-c)*abs(z)/z +c)
    erfg=f

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION erfg

  FUNCTION casim_erfc(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='CASIM_ERFC'

    REAL(wp), INTENT(IN) :: x
    INTEGER, PARAMETER :: c=1
    REAL(wp) :: casim_erfc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    casim_erfc=erfg(x,c)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION casim_erfc

  ! Inverse of error function
  !
  ! This needs more work to get good accuracy
  FUNCTION erfinv1(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ERFINV1'

    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: erfinv1

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    erfinv1=.5*sqrt(pi)*(x+pi/12.0*x*x*x+7.0/480.0*pi*pi*x**5 &
         +127.0/40320*pi**3*x**7+4369.0/5806080*pi**4*x**9 &
         +34807.0/182476800.0*pi**5*x**11)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION erfinv1

  ! Inverse of error function
  !
  ! Alternative version solves equation
  FUNCTION erfinv2(x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ERFINV2'

    REAL(wp), INTENT(IN) :: x
    REAL(wp) :: erfinv2

    REAL(wp) :: work, work_old, diff, erfx, erfx_old

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    diff=9999.0
    work_old=.2
    work=1.0
    DO WHILE(abs(diff) > 1e-3)
      erfx=erf(work)-x
      erfx_old=erf(work_old)-x
      diff=-erfx*(work_old-work)/(erfx_old-erfx)
      work_old=work
      work=work + diff
    END DO

    erfinv2=work

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION erfinv2

  ! Inverse of error function
  !
  ! Alternative version solves equation
  FUNCTION erfinv3(x, tol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ERFINV3'

    REAL(wp), INTENT(IN) :: x
    REAL(wp), OPTIONAL, INTENT(IN) :: tol
    REAL(wp) :: erfinv3

    REAL(wp) :: work, diff, erfx, derfx, tolval

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    tolval=1e-3
    IF (present(tol)) tolval=tol

    IF (abs(x) > .95) THEN
      ! don't converge well for abs(x)->1
      ! should treat this properly
      ! (i.e. more sophisticated solver)
      ! but don't really care too much
      ! about these values for now
      work=erfinv1(x)
    ELSE

      diff=9999.0
      work=erfinv1(x)
      DO WHILE(abs(diff) > tolval)
        erfx=erf(work)-x
        derfx=2.0*exp(-work*work)/sqrt(pi)
        diff=-erfx/derfx
        work=work+diff
      END DO
    END IF
    erfinv3=work

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION erfinv3

END MODULE special
