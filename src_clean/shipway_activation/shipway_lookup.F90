! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Lookup-table parameters/grid used by the Shipway activation scheme.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE shipway_lookup

  USE variable_precision, ONLY: wp

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SHIPWAY_LOOKUP'

  REAL(wp), PARAMETER :: tol=1.e-9

  ! parameters for lookup tables
  REAL(wp), PARAMETER :: xmin=1.e-4  ! May need to adjust if parameters go out of range
  REAL(wp), PARAMETER :: xmax=20000.

  REAL(wp), PARAMETER :: ymin=.1
  REAL(wp), PARAMETER :: ymax=3.

  INTEGER, PARAMETER :: nx=500
  INTEGER, PARAMETER :: ny=10

  INTEGER, PARAMETER :: nterms=10000 ! We can afford to be generous

  LOGICAL :: l_generate_tables=.TRUE.

  REAL(wp) :: x_param, y_param
  REAL(wp) :: xvalues(nx)
  REAL(wp) :: yvalues(ny)
  REAL(wp) :: J_table(nx,ny)
  
  PUBLIC lookup_I, xmax, xmin, ymax, ymin, generate_tables

CONTAINS

  SUBROUTINE generate_tables(filename)
    ! NB The lookup table is stored as I(x,y)/x^2, then the x^2 
    ! term is restored after the interpolation is done in lookup_I  

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(*), OPTIONAL, INTENT(IN) :: filename
    REAL(wp) :: dlx,dly
    LOGICAL :: fexist
    LOGICAL :: l_calculate=.TRUE., l_writeout=.FALSE.
    INTEGER :: i,j 

    CHARACTER(len=*), PARAMETER :: RoutineName='GENERATE_TABLES'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    IF (present(filename))THEN
      INQUIRE(FILE=filename, EXIST=fexist)
      IF (fexist)THEN
        OPEN(200, file=trim(filename), STATUS='old')
        READ(200, *) xvalues
        READ(200, *) yvalues
        READ(200, *) J_table
        l_calculate=.FALSE.
      ELSE
        l_writeout=.TRUE.
      END IF
    END IF

    IF (l_calculate)THEN

      dlx=(log(xmax)-log(xmin))/(nx-1)
      dly=(log(ymax)-log(ymin))/(ny-1)
      DO i=1,nx
        xvalues(i) = xmin*exp((i-1)*dlx)
      END DO

      DO j=1,ny
        yvalues(j) = ymin*exp((j-1)*dly)
      END DO

      DO j=1,ny
        DO i=1,nx
          J_table(i,j) = J1(xvalues(i),yvalues(j),nterms)/(xvalues(i)*xvalues(i))
        END DO
      END DO

      ! Quick hack to ensure J1 is monotonic increasing.
      DO j=1,ny
        DO i=2,nx
          IF (J_table(i,j) <= J_table(i-1,j)) J_table(i,j) = J_table(i-1,j) + tiny(xmin)
        END DO
      END DO

      IF (l_writeout)THEN
        OPEN(200, file=trim(filename), STATUS='new')
        WRITE(200, *) xvalues
        WRITE(200, *) yvalues
        WRITE(200, *) J_table
      END IF
    END IF

    l_generate_tables=.FALSE.

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE generate_tables

  REAL(wp) FUNCTION J1_integrand(t)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! See equation 21 of Shipway (2015)    
    REAL(wp), INTENT(IN) :: t

    REAL(wp) :: x, y, expnt

    ! local variables

    CHARACTER(len=*), PARAMETER :: RoutineName='J1_INTEGRAND'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    x=x_param
    y=y_param

    ! Limit exponent to prevent underflow - if it's small enough
    ! we don't care about the loss of accuracy
    ! This allows for the factor multiplying the exponential term to be
    ! a minimum of e**(-20.) = 2.e-9
    ! Not that we still expect some floating underflow in the
    ! part of the parameter space where the integral < tiny(1.0)
    ! This doesn't really matter since the left had side will presumabely be > tiny(1.0)!
!    expnt=-0.5*log(t)*log(t)/(y*y)
    expnt=max(-0.5*log(t)*log(t)/(y*y), (minexponent(t))*log(2.0) + 20.0)
    J1_integrand = x*sqrt(x*x-t*t)/t*exp(expnt) &
       /sqrt(((x**3-t**3)/x**3)**0.6)/sqrt(0.5)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)


  END FUNCTION J1_integrand

  REAL(wp) FUNCTION J1(x, y, nterms)
    
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: x, y
    INTEGER, INTENT(IN) :: nterms
    CHARACTER(len=*), PARAMETER :: RoutineName='J1'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    x_param=x
    y_param=y
    
    J1 = simpson(J1_integrand, tol, x - tol, nterms)

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION J1

  SUBROUTINE lookup_I(x,y,J1)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! NB The lookup table is stored as I(x,y)/x^2, then the x^2 
    ! term is restored after the interpolation is done in lookup_I  

    REAL(wp), INTENT(IN) :: x, y
    REAL(wp), INTENT(OUT) :: J1
    INTEGER :: m(1)
    INTEGER :: ix, iy
    REAL(wp) :: dx1, dy1, dx2, dy2

    CHARACTER(len=*), PARAMETER :: RoutineName='LOOKUP_I'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    IF (l_generate_tables)THEN
      CALL generate_tables()
    END IF

    m=minloc((xvalues(:) - x)*(xvalues(:) - x))
    ix=m(1)
    m=minloc((yvalues(:) - y)*(yvalues(:) - y))
    iy=m(1)

    IF (x < xvalues(ix)) ix=ix-1
    IF (y < yvalues(iy)) iy=iy-1
    
    dx1=x-xvalues(ix)
    dx2=xvalues(ix+1)-x
    dy1=y-yvalues(iy)
    dy2=yvalues(iy+1)-y

    J1=x*x*(1./((dx1+dx2)*(dy1+dy2))) &
       *(J_table(ix,iy)*dx2*dy2 &
       + J_table(ix+1,iy)*dx1*dy2 &
       + J_table(ix,iy+1)*dx2*dy1 &
       + J_table(ix+1,iy+1)*dx1*dy1)

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE lookup_I


  FUNCTION simpson(f,a,b,N)
    
    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp) :: f
    REAL(wp), INTENT(IN) :: a,b
    INTEGER, INTENT(IN) :: N
    
    REAL(wp) :: s, simpson, fk1,fk2,fk3
    REAL(wp) :: h
    INTEGER  :: i
    REAL(wp) :: tol ! to mitigate rounding issues at single precision

    CHARACTER(len=*), PARAMETER :: RoutineName='SIMPSON'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    s=0
    
    tol=(b-a)*1e-7
    
    h = (b-a-tol)/(2*N)
    fk3=f(a)
    DO i=0,N-1
       fk1=fk3
       fk2=f(a+(2*i+1)*h)
       fk3=f(a+2*(i+1)*h)
       s=s+fk1+4*fk2+fk3
    END DO
    
    simpson=h*s/3.

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION simpson

END MODULE shipway_lookup
