! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Single- and two-species hydrometeor collection/sweep-out rate kernels used by the accretion/aggregation routines.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE sweepout_rate
  USE variable_precision, ONLY: wp
  USE special, ONLY: pi, Gammafunc
  USE mphys_constants, ONLY: rho0
  USE mphys_parameters, ONLY: hydro_params
  
  IMPLICIT NONE
  PRIVATE
  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SWEEPOUT_RATE'

  PUBLIC sweepout, binary_collection, sweepout_1M2M, binary_collection_1M2M
CONTAINS

  ! Calculate the sweepout rate given the distribution
  ! and fallspeed parameters
  FUNCTION sweepout(n0, lam, mu, params, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SWEEPOUT'

    REAL(wp), INTENT(IN) :: mu, n0, lam
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted sweepout
    LOGICAL, INTENT(IN), OPTIONAL :: mass_weight
    REAL(wp) :: sweepout

    ! local variables
    REAL(wp) :: G3_b_mu, G1_mu !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: arg3 !< argument for gamma function
    REAL(wp) :: coef !< coefficient = c_x if mass-weighting

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg3=3.0+params%b_x+mu
    coef=1.0
    IF (present(mass_weight)) THEN
      IF (mass_weight) THEN
        arg3=arg3+params%d_x
        coef=params%c_x
      END IF
    END IF

    G3_b_mu=GammaFunc(arg3)
    G1_mu=GammaFunc(1.0+mu)
    sweepout=coef*(pi*n0*params%a_x/4.0)*G3_b_mu/G1_mu* (1.0 + params%f_x/lam)**(-arg3)&
         *lam**(1 + mu - arg3)* (rho0/rho)**(params%g_x)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION sweepout

 ! Calculate the sweepout rate given the distribution
  ! and fallspeed parameters
  FUNCTION sweepout_1M2M(n0, lam, params, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SWEEPOUT_1M2M'

    REAL(wp), INTENT(IN) :: n0, lam
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted sweepout
    LOGICAL, INTENT(IN), OPTIONAL :: mass_weight
    REAL(wp) :: sweepout_1M2M

    ! local variables
    REAL(wp) :: G3_b_mu, G1_mu !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: arg3 !< argument for gamma function
    REAL(wp) :: coef !< coefficient = c_x if mass-weighting

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg3=3.0+params%b_x+params%fix_mu
    G3_b_mu=params%gam_3_bx_mu
    coef=1.0

    IF (present(mass_weight)) THEN
      IF (mass_weight) THEN
         arg3=arg3+params%d_x
         G3_b_mu=params%gam_3_bx_mu_dx
         coef=params%c_x
      END IF
    END IF

    G1_mu=params%gam_1_mu

    sweepout_1M2M=coef*(pi*n0*params%a_x/4.0)*G3_b_mu/G1_mu* (1.0 + params%f_x/lam)**(-arg3)&
         *lam**(1 + params%fix_mu - arg3)* (rho0/rho)**(params%g_x)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION sweepout_1M2M

  ! Calculate the number of species Y collected by species X through
  ! binary collisions as both species sediment
  FUNCTION binary_collection(n0_Xin, lam_X, mu_X, n0_Yin, lam_Y, mu_Y,   &
       params_X, params_Y, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='BINARY_COLLECTION'

    REAL(wp), INTENT(IN) :: mu_X, n0_Xin, lam_X
    REAL(wp), INTENT(IN) :: mu_Y, n0_Yin, lam_Y
    REAL(wp) :: n0_X, n0_Y
    TYPE(hydro_params), INTENT(IN) :: params_X
    TYPE(hydro_params), INTENT(IN) :: params_Y
    REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted value (i.e. total mass accreted)
    LOGICAL, INTENT(IN), OPTIONAL :: mass_weight
    REAL(wp) :: binary_collection

    ! local variables
    REAL(wp) :: G1_X, G2_X, G3_X  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: G1_Y, G2_Y, G3_Y  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: arg1_X !< argument for gamma function
    REAL(wp) :: arg2_X !< argument for gamma function
    REAL(wp) :: arg3_X !< argument for gamma function
    REAL(wp) :: arg1_Y !< argument for gamma function
    REAL(wp) :: arg2_Y !< argument for gamma function
    REAL(wp) :: arg3_Y !< argument for gamma function
    REAL(wp) :: l_X1 !< lam_X^-arg1_X
    REAL(wp) :: l_X2 !< lam_X^-arg2_X
    REAL(wp) :: l_X3 !< lam_X^-arg3_X
    REAL(wp) :: l_Y1 !< lam_Y^-arg1_Y
    REAL(wp) :: l_Y2 !< lam_Y^-arg2_Y
    REAL(wp) :: l_Y3 !< lam_Y^-arg3_Y

    REAL(wp) :: coef !< coefficient = c_x if mass-weighting
    REAL(wp) :: V_X  !< mass-weighted fall velocity for X
    REAL(wp) :: V_Y  !< mass-weighted fall velocity for Y
    REAL(wp) :: delV !< bulk fall-speed differential

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    arg1_X=1.0+mu_X
    arg2_X=2.0+mu_X
    arg3_X=3.0+mu_X
    arg1_Y=1.0+mu_Y
    arg2_Y=2.0+mu_Y
    arg3_Y=3.0+mu_Y
    coef=1.0
    IF (present(mass_weight)) THEN
      IF (mass_weight) THEN
        arg1_Y=arg1_Y+params_Y%d_x
        arg2_Y=arg2_Y+params_Y%d_x
        arg3_Y=arg3_Y+params_Y%d_x
        coef=params_Y%c_x
      END IF
    END IF

    G1_X=GammaFunc(arg1_X)
    G2_X=GammaFunc(arg2_X)
    G3_X=GammaFunc(arg3_X)
    G1_Y=GammaFunc(arg1_Y)
    G2_Y=GammaFunc(arg2_Y)
    G3_Y=GammaFunc(arg3_Y)

    l_X1=lam_X**(-arg1_X)
    l_X2=lam_X**(-arg2_X)
    l_X3=lam_X**(-arg3_X)
    l_Y1=lam_Y**(-arg1_Y)
    l_Y2=lam_Y**(-arg2_Y)
    l_Y3=lam_Y**(-arg3_Y)

    n0_X=n0_Xin *lam_X**(mu_X+1.0)/GammaFunc(1.0+mu_X)
    n0_Y=n0_Yin *lam_Y**(mu_Y+1.0)/GammaFunc(1.0+mu_Y)
    
    V_X=params_X%a_x * lam_X**(-params_X%b_x)*(rho0/rho)**(params_X%g_x)     &
         *GammaFunc(1.0+mu_X+params_X%d_x+params_X%b_x)/GammaFunc(1.0+mu_X+params_X%d_x)

    V_Y=params_Y%a_x * lam_Y**(-params_Y%b_x)*(rho0/rho)**(params_Y%g_x)     &
         *GammaFunc(1.0+mu_Y+params_Y%d_x+params_Y%b_x)/GammaFunc(1.0+mu_Y+params_Y%d_x)

    delV=max(max(V_X,V_Y)/4.0, abs(V_X-V_Y))

    binary_collection=coef*0.25*pi*n0_X*n0_Y*delV*(l_X3*l_Y1*G1_Y*G3_X+2.0*l_X2*l_Y2*G2_Y*G2_X+l_X1*l_Y3*G3_Y*G1_X)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION binary_collection

  FUNCTION binary_collection_1M2M(n0_Xin, lam_X, n0_Yin, lam_Y,   &
       params_X, params_Y, rho, mass_weight)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='BINARY_COLLECTION_1M2M'

    REAL(wp), INTENT(IN) :: n0_Xin, lam_X
    REAL(wp), INTENT(IN) :: n0_Yin, lam_Y
    REAL(wp) :: n0_X, n0_Y, mu_X, mu_Y
    
    TYPE(hydro_params), INTENT(IN) :: params_X
    TYPE(hydro_params), INTENT(IN) :: params_Y
    REAL(wp), INTENT(IN) :: rho ! air density
    ! if present and true use mass-weighted value (i.e. total mass accreted)
    LOGICAL, INTENT(IN), OPTIONAL :: mass_weight
    REAL(wp) :: binary_collection_1M2M

    ! local variables
    REAL(wp) :: G1_X, G2_X, G3_X  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: G1_Y, G2_Y, G3_Y  !< gamma functions (these may be taken directly from
    !< params once code has been made more efficent)
    REAL(wp) :: arg1_X !< argument for gamma function
    REAL(wp) :: arg2_X !< argument for gamma function
    REAL(wp) :: arg3_X !< argument for gamma function
    REAL(wp) :: arg1_Y !< argument for gamma function
    REAL(wp) :: arg2_Y !< argument for gamma function
    REAL(wp) :: arg3_Y !< argument for gamma function
    REAL(wp) :: l_X1 !< lam_X^-arg1_X
    REAL(wp) :: l_X2 !< lam_X^-arg2_X
    REAL(wp) :: l_X3 !< lam_X^-arg3_X
    REAL(wp) :: l_Y1 !< lam_Y^-arg1_Y
    REAL(wp) :: l_Y2 !< lam_Y^-arg2_Y
    REAL(wp) :: l_Y3 !< lam_Y^-arg3_Y

    REAL(wp) :: coef !< coefficient = c_x if mass-weighting
    REAL(wp) :: V_X  !< mass-weighted fall velocity for X
    REAL(wp) :: V_Y  !< mass-weighted fall velocity for Y
    REAL(wp) :: delV !< bulk fall-speed differential

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    mu_X=params_X%fix_mu
    mu_Y=params_Y%fix_mu

    arg1_X=1.0+params_X%fix_mu
    arg2_X=2.0+params_X%fix_mu
    arg3_X=3.0+params_X%fix_mu
    arg1_Y=1.0+params_Y%fix_mu
    arg2_Y=2.0+params_Y%fix_mu
    arg3_Y=3.0+params_Y%fix_mu
    
    G1_X=params_X%gam_1_mu
    G2_X=params_X%gam_2_mu
    G3_X=params_X%gam_3_mu
    G1_Y=params_Y%gam_1_mu
    G2_Y=params_Y%gam_2_mu
    G3_Y=params_Y%gam_3_mu
    coef=1.0
    IF (present(mass_weight)) THEN
       IF (mass_weight) THEN
          arg1_Y=arg1_Y+params_Y%d_x
          arg2_Y=arg2_Y+params_Y%d_x
          arg3_Y=arg3_Y+params_Y%d_x
          G1_Y=params_Y%gam_1_mu_dx
          G2_Y=params_Y%gam_2_mu_dx
          G3_Y=params_Y%gam_3_mu_dx
          coef=params_Y%c_x
       END IF
    END IF
        
    l_X1=lam_X**(-arg1_X)
    l_X2=lam_X**(-arg2_X)
    l_X3=lam_X**(-arg3_X)
    l_Y1=lam_Y**(-arg1_Y)
    l_Y2=lam_Y**(-arg2_Y)
    l_Y3=lam_Y**(-arg3_Y)

    n0_X=n0_Xin *lam_X**(mu_X+1.0)/GammaFunc(1.0+mu_X)
    n0_Y=n0_Yin *lam_Y**(mu_Y+1.0)/GammaFunc(1.0+mu_Y)

    V_X=params_X%a_x * lam_X**(-params_X%b_x)*(rho0/rho)**(params_X%g_x)     &
         *params_X%gam_1_mu_dx_bx/params_X%gam_1_mu_dx
         !*GammaFunc(1.0+mu_X+params_X%d_x+params_X%b_x)/GammaFunc(1.0+mu_X+params_X%d_x)

    V_Y=params_Y%a_x * lam_Y**(-params_Y%b_x)*(rho0/rho)**(params_Y%g_x)     &
         *params_Y%gam_1_mu_dx_bx/params_Y%gam_1_mu_dx
         !*GammaFunc(1.0+mu_Y+params_Y%d_x+params_Y%b_x)/GammaFunc(1.0+mu_Y+params_Y%d_x)

    delV=max(max(V_X,V_Y)/4.0, abs(V_X-V_Y))

    binary_collection_1M2M=coef*0.25*pi*n0_X*n0_Y*delV*(l_X3*l_Y1*G1_Y*G3_X+2.0*l_X2*l_Y2*G2_Y*G2_X+l_X1*l_Y3*G3_Y*G1_X)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION binary_collection_1M2M

END MODULE sweepout_rate
