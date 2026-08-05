! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Algebraic third-moment increments on phase change/process rates (m3_inc_type1-4).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Increments to a hydrometeor's third prognostic moment (e.g.
!   reflectivity-related moment for triple-moment species) when
!   mass and/or number are incremented by a process rate.
!
! Paper reference:
!   Field et al. (2023) Sec. 1 notes rain/snow/graupel can carry
!   a third prognostic moment; no numbered equation is given for
!   these specific moment-increment helper routines.
!
MODULE m3_incs
  USE variable_precision, ONLY: wp
  USE lookup, ONLY: Gfunc
  USE mphys_die, ONLY: throw_mphys_error, warn, std_msg
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='M3_INCS'

  PUBLIC m3_inc_type4, m3_inc_type3, m3_inc_type2, m3_inc_type1
CONTAINS

  ! for changes of phase from category y to category x
  ! dM3_x = -(c_y/c_x)**(p3/3) * dM3_y
  ! c_x, c_y are densities of x and y respectively
  ! p3 is the third moment (must be the same for both x and y)
  ! This assumes both categories are spherical.
  SUBROUTINE m3_inc_type4(dm3_y, c_x, c_y, p3, dm3_x)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='M3_INC_TYPE4'

    REAL(wp), INTENT(IN) :: dm3_y
    REAL(wp), INTENT(IN) :: c_x,c_y
    REAL(wp), INTENT(IN) :: p3
    REAL(wp), INTENT(OUT) :: dm3_x

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm3_x = -(c_y/c_x)**(p3/3.0) * dm3_y

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE m3_inc_type4

  ! tendency of m3 as a function of dm1, dm2
  ! assuming initial value for mu
  SUBROUTINE m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='M3_INC_TYPE3'

    REAL(wp), INTENT(IN) :: dm1, dm2
    REAL(wp), INTENT(IN) :: p1, p2, p3
    REAL(wp), INTENT(IN) :: mu
    REAL(wp), INTENT(OUT) :: dm3
    
    REAL(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3 = p2-p1
    k1 = (p3-p2)/k3
    k2 = (p1-p3)/k3

    dm3 = (Gfunc(mu, p1, p2, p3)*dm1**(-k1)*dm2**(-k2))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE m3_inc_type3

  ! tendency of m3 as a function of m1,m2,dm1,dm2
  ! assuming shape parameter does not vary
  SUBROUTINE m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3, mu_init)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='M3_INC_TYPE2'
   
    REAL(wp), INTENT(IN) :: m1, m2, m3, dm1, dm2
    REAL(wp), INTENT(IN) :: p1, p2, p3
    REAL(wp), INTENT(OUT) :: dm3
    REAL(wp), INTENT(IN), OPTIONAL :: mu_init ! initial mu for type3
    
    REAL(wp) :: k1, k2, k3
    REAL(wp) :: fac

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (m1 > 0.0) THEN
      k1 = p3-p2
      k2 = p1-p3
      k3 = p2-p1

      IF (abs(p3-6.0)< spacing(p3)) THEN ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to...
        fac=((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2))-1.0
      ELSE IF (abs(p3-4.5)< spacing(p3)) THEN ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to.
        fac=sqrt(((m1+dm1)/m1)*((m1+dm1)/m1)*((m1+dm1)/m1)*(m2/(m2+dm2)))-1.0
      ELSE IF (abs(p3-1.5)< spacing(p3)) THEN ! p3=6
        ! since we've hardwired p1=3 and p2=0, we can reduce the calculation to...
        fac=sqrt(((m1+dm1)/m1)*((m2+dm2)/m2))-1.0
      ELSE
        !      fac=((m1/(m1+dm1))**(k1/k3)*(m2/(m2+dm2))**(k2/k3)-1.)
        fac=exp((k1*log(m1/(m1+dm1)) + k2*log(m2/(m2+dm2)))/k3)-1.0
      END IF
      IF (fac < -1.0) THEN

        WRITE(std_msg, *) 'm3inc_2 ERROR:', fac, m1,m2,m3,dm1,dm2,k1,k2,k3, &
                            sqrt(fac)

        CALL throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)

      END IF
      dm3 = m3*fac
    ELSE ! If there is no pre-existing mass, then use type 3
      CALL m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_init)
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE m3_inc_type2

  ! tendency of m2 as a function of m1,dm1
  ! assuming shape parameter and slope do not vary
  SUBROUTINE m3_inc_type1(m1, m2, dm1, dm2)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='M3_INC_TYPE1'

    REAL(wp), INTENT(IN) :: m1, m2, dm1
    REAL(wp), INTENT(OUT) :: dm2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm2 = dm1*m1/m2

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE m3_inc_type1
END MODULE m3_incs
