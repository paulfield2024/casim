! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Gamma-function lookup helpers and size-distribution-parameter inversion (Gfunc, moment, get_n0, get_mu, get_slope_generic family).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE lookup
  USE mphys_die, ONLY: throw_mphys_error, bad_values, std_msg
  USE variable_precision, ONLY: wp
  USE mphys_switches, ONLY: max_mu, l_kfsm
! use mphys_switches, only: l_passive3m
  USE mphys_parameters, ONLY: hydro_params, a_i, b_i, &
                              a_s, b_s
  USE special, ONLY: GammaFunc
  USE passive_fields, ONLY: rho
!#if DEF_MODEL==UM
  USE casim_moments_mod, ONLY: casim_moments
!#endif
  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='LOOKUP'

  INTEGER, PARAMETER :: nmu=501
  REAL(wp) :: min_mu = 0.0  ! ( max_mu = 35 )
  REAL(wp), ALLOCATABLE :: mu_g(:), mu_i(:)
  REAL(wp), ALLOCATABLE :: mu_g_sed(:), mu_i_sed(:)

  INTERFACE get_lam_n0
     MODULE PROCEDURE get_lam_n0_3M, get_lam_n0_2M, get_lam_n0_1M, get_lam_n0_1M_KF
  END INTERFACE get_lam_n0

  PUBLIC Gfunc, get_slope_generic, get_slope_generic_kf, moment, get_n0,        &
         get_mu, get_lam_n0, set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
CONTAINS

  FUNCTION Gfunc(mu, p1, p2, p3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GFUNC'

    REAL(wp), INTENT(IN) :: mu, p1, p2, p3

    REAL(wp) :: Gfunc
!   real(wp) :: GfuncL
    REAL(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    !    GfuncL=exp(k1*log(GammaFunc(1.+mu+p1)) &
    !         +k2*log(GammaFunc(1.+mu+p2)) &
    !         +log(GammaFunc(1.+mu+p3)) &
    !         )

    Gfunc=GammaFunc(1.0+mu+p1)**k1*GammaFunc(1.0+mu+p2)**k2*GammaFunc(1.0+mu+p3)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION Gfunc

  FUNCTION Hfunc(m1,m2,m3, p1, p2, p3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='HFUNC'

    REAL(wp), INTENT(IN) :: m1,m2,m3
    REAL(wp), INTENT(IN) :: p1,p2,p3
    REAL(wp) :: Hfunc
    REAL(wp) :: k1, k2, k3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    k3=p2-p1
    k1=(p3-p2)/k3
    k2=(p1-p3)/k3

    Hfunc=exp(k1*log(m1)+k2*log(m2)+log(m3))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION Hfunc

  SUBROUTINE set_mu_lookup(p1, p2, p3, ind, val)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SET_MU_LOOKUP'

    REAL(wp), INTENT(IN) :: p1, p2, p3
    REAL(wp), INTENT(INOUT) :: ind(:), val(:)

    INTEGER :: i, lb, ub

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    lb=lbound(ind,1)
    ub=ubound(ind,1)

    DO i=lb, ub
      ind(i)=min_mu+(max_mu-min_mu)/(nmu-1.0)*(i-1)
      val(i)=Gfunc(ind(i), p1, p2, p3)
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE set_mu_lookup

  SUBROUTINE get_slope_generic(params, n0, lam, mu, mass, num, m3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_SLOPE_GENERIC'

    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(OUT) :: n0, lam, mu
    REAL(wp), INTENT(IN) :: mass
    REAL(wp), INTENT(IN), OPTIONAL :: num, m3
                                    !! m3 NEEDED for 3rd moment code

    REAL(wp) :: m1, m2, p1, p2, p3

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    m1=mass/params%c_x
    ! p1=params%p1
    ! p2=params%p2
    ! p3=params%p3

    ! if (params%l_3m) then
    !   if (l_passive3m) then
    !     m2=num
    !     mu=params%fix_mu
    !      call get_lam_n0(m1, m2, params, lam, n0, params%id)
    !   else
    !     m2=num
    !     call get_mu(m1, m2, m3, p1, p2, p3, mu)
    !     call get_lam_n0(m1, m2, m3, params, mu, lam, n0)
    !   end if
    ! elseif (params%l_2m) then
    IF (params%l_2m) THEN
      m2=num
      mu=params%fix_mu
      !call get_lam_n0(m1, m2, p1, p2, mu, lam, n0)
      CALL get_lam_n0(m1, m2, params, lam, n0, params%id)
    ELSE
      mu=params%fix_mu
      n0=params%fix_n0
      !call get_lam_n0(m1, p1, mu, lam, n0)
      CALL get_lam_n0(m1, params, lam, n0)
    END IF
    IF (lam <= 0) THEN
      WRITE(std_msg, *) 'ERROR in lookup', params%id, params%i_2m, m1, num, lam
      CALL throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_slope_generic

  SUBROUTINE get_slope_generic_kf(ixy_inner, k, params, n0, lam, mu, lams, mass, Tk, num, &
                               m3)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_SLOPE_GENERIC_KF'

    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: k
    TYPE(hydro_params), INTENT(INOUT) :: params
    REAL(wp), INTENT(OUT) :: n0, lam, mu
    REAL(wp), INTENT(OUT) :: lams(2)
    REAL(wp), INTENT(IN) :: mass, Tk
    REAL(wp), INTENT(IN), OPTIONAL :: num, m3

    REAL(wp) :: m1, ms, p1, p2, p3
    REAL(wp) :: n_p, cficei(1), m_s(1), Ta(1), rhoa(1), qcf(1)
    REAL(wp) :: j1, j2
    REAL(wp) :: na, nb

    INTEGER  :: points

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    m1=mass/params%c_x
    p1=params%p1
    p2=params%p2
    p3=params%p3

    IF (l_kfsm) THEN

      j1=p1+params%b_x
      j2=p1+b_i
      ! modify p2 so it is equal to j1
      params%p2 = j1

      points = 1
      n_p = 0.0
      m_s(1)=0.0
      cficei(1) = 1.0
      rhoa(1) = rho(k,ixy_inner)
      qcf(1) = mass
      Ta(1) = Tk

    END IF

      mu=params%fix_mu
      n0=params%fix_n0
!#if DEF_MODEL==UM
!AH-KID changes - is this needed?
!prf make n0 fixed param consistent with Field 2017
      IF (params%id==3 .AND. l_kfsm) THEN
        qcf(1)=max(qcf(1),1e-8) !stop lam going to zero
        CALL casim_moments(params,points,rhoa,Ta,qcf,cficei,0.0_wp,m_s) !get concentration
        n0=m_s(1)
        m_s(1)=0.0
      END IF
!prf
!AH-KID changes 
!#endif
      CALL get_lam_n0(m1, params, lam, n0, params%id)

      IF (l_kfsm) THEN
        ! Code for Kalli's single moment work
        IF (params%id==2) THEN

          ! Rain: Abel and Boutle distribution
          na=0.22
          nb=2.2
          n0=na*params%gam_1_mu*lam**(nb-1.0-params%fix_mu)

        ELSE IF (params%id==5) THEN

          ! Graupel: Swann PSD
          ! JW: Note we could eventually look to using Field et al. (2019)
          !     PSD - JAMC.
          na=5.0e25
          nb=-4.0
          n0=na*params%gam_1_mu*lam**(nb-1.0-params%fix_mu)

        ELSE IF (params%id==3) THEN
!#if DEF_MODEL==UM
          ! Ice: need lsp_moments
          CALL casim_moments(params,points,rhoa,Ta,qcf,cficei,j1,m_s)

          ms=m_s(1)

          CALL casim_moments(params,points,rhoa,Ta,qcf,cficei,j2,m_s)

          IF (m_s(1) < ms*params%a_x/a_i) THEN
            j1=j2
            ms=m_s(1)
            a_s(k)=a_i
            b_s(k)=b_i
          ELSE
            a_s(k) = params%a_x
            b_s(k) = params%b_x
          END IF ! m_s(1)

          ! lams(1) is the slope, while lams(2) is n0
          CALL get_lam_n0(m1, ms, params, lams(1), lams(2), params%id)
!#endif
        END IF ! params%id

      END IF ! l_kfsm

    IF (lam <= 0) THEN
      WRITE(std_msg, *) 'ERROR in lookup', params%id, params%i_2m, m1, num, lam, n0, qcf(1)
      CALL throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_slope_generic_kf

  SUBROUTINE get_mu(m1, m2, m3, p1, p2, p3, mu, mu_g_o, mu_i_o)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_MU'

    REAL(wp), INTENT(IN) :: m1, m2, m3
    REAL(wp), INTENT(IN) :: p1, p2, p3
    REAL(wp), INTENT(IN), OPTIONAL :: mu_g_o(:), mu_i_o(:)
    REAL(wp), INTENT(OUT) :: mu

    REAL(wp) :: G
    REAL(wp) :: muG(nmu), muI(nmu)
    INTEGER  :: i
    REAL(wp) :: k1, k2, pos

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (present(mu_g_o) .AND. present(mu_i_o)) THEN
      muG=mu_g_o
      muI=mu_i_o
    ELSE
      muG=mu_g
      muI=mu_i
    END IF

    k1=p3-p2
    k2=p1-p3

    pos=sign(1.0_wp,k1*k2) ! scaling to +-1 preserves precision for < condition

    G=Hfunc(m1,m2,m3,p1,p2,p3)

    IF (pos*G < pos*muG(1)) THEN
      mu=-1.0e-10  ! set to be small and negative so that it is picked up
      ! in checks later on
    ELSE
      DO i=1, nmu-1
        IF ((muG(i)-G)*(muG(i+1)-G) <= 0.0) EXIT
      END DO

      IF (i==nmu) THEN
        mu=max_mu
      ELSE
        mu=muI(i)+(G-muG(i))/(muG(i+1)-muG(i))*(muI(i+1)-muI(i))
      END IF
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_mu

  ! 3M version
  SUBROUTINE get_lam_n0_3M(m1, m2, m3, params, mu, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_LAM_N0_3M'

    REAL(wp), INTENT(IN) :: m1, m2, m3, mu
                                !! m3 NEEDED for 3rd moment code
    TYPE(hydro_params), INTENT(IN) :: params
    !real(wp), intent(in) :: p1, p2, p3
    REAL(wp), INTENT(OUT) :: lam, n0

    REAL(wp) :: p, m
    REAL(wp) :: l2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !    l2 = 1./(p2-p3)
    !
    !    lam = ((GammaFunc(1.+mu+p2)/GammaFunc(1.+mu+p3)) &
    !         *(m3/m2))**l2
    ! Make sure we use m1 and m2 to calculate lambda, in case we've gone
    ! out of bounds with mu and need to modify m3.
    l2=1.0/(params%p2-params%p1)

    lam=((GammaFunc(1.0+mu+params%p2)/GammaFunc(1.0+mu+params%p1))*(m1/m2))**l2

    m=m2
    p=params%p2

    n0=lam**(p)*m*GammaFunc(1.0+mu)/GammaFunc(1.0+mu+p)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_lam_n0_3M

  ! 2M version
  SUBROUTINE get_lam_n0_2M(m1, m2, params, lam, n0, id)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_LAM_N0_2M'

    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: m1, m2 !, mu
    !real(wp), intent(in) :: p1, p2
    REAL(wp), INTENT(OUT) :: lam, n0
    INTEGER, INTENT(IN) :: id

    REAL(wp) :: p, m

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    IF (l_kfsm .AND. id == 3) THEN 
       ! field formulation of ice psd, p2 can vary, so need to calc the gamma and inverse  
       ! each timestep. Note p2 is set outside the call - nasty!
       lam=((params%gam_1_mu_p1/GammaFunc(1.0+params%fix_mu+params%p2)) &
            *(m2/m1))**(1.0/(params%p1-params%p2))

       m=m2
       p=params%p2

       n0=lam**(p)*m*params%gam_1_mu/GammaFunc(1.0+params%fix_mu+p)
    ELSE
       lam = ((params%gam_1_mu_p1/params%gam_1_mu_p2)* &
            (m2/m1))**(params%inv_p1_p2)
       
 
       n0=lam**(params%p2)*m2*params%gam_1_mu/params%gam_1_mu_p2
    END IF
 
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_lam_n0_2M

  ! 1M version
  SUBROUTINE get_lam_n0_1M(m1, params, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_LAM_N0_1M'

    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: m1, n0 
    REAL(wp), INTENT(OUT) :: lam

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !
    ! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
    ! (c.f. LEM formulation, na, nb)
    ! if we want to fix na and nb, we can do this and lambda is
    ! given by:
    ! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
    !

    lam = (n0*params%gam_1_mu_p1/params%gam_1_mu*m1**(-1.0))**(params%inv_p1)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_lam_n0_1M

  SUBROUTINE get_lam_n0_1M_KF(m1, params, lam, n0, id)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_LAM_N0_1M_KF'

    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: m1, n0
    !real(wp), intent(in) :: p1
    REAL(wp), INTENT(OUT) :: lam
    INTEGER, INTENT(IN) :: id
    REAL(wp) :: na, nb

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !
    ! Fixing Nx is equivalent to having n0=na*lam**(1+mu)
    ! (c.f. LEM formulation, na, nb)
    ! if we want to fix na and nb, we can do this and lambda is
    ! given by:
    ! lam=(na*gamma(1+mu+p1)/m1)**(1/(1+mu+p1-nb))
    !

    IF (l_kfsm .AND. id == 2) THEN
      ! Rain: Abel and Boutle (2012)
      na  = 0.22
      nb  = 2.2
      lam = (na*params%gam_1_mu_p1*m1**(-1.0)) &
           **(1.0/(1+params%fix_mu+params%p1-nb))
    ELSE IF (l_kfsm .AND. id == 5) THEN
      ! Graupel: Swann (LEM); Forbes and Halliwell (2003, internal report)
      na  = 5.0e25
      nb  = -4.0
      lam=(na*params%gam_1_mu_p1*m1**(-1.0)) &
           **(1.0/(1+params%fix_mu+params%p1-nb))
    ELSE
      ! not l_kfsm or species other than rain or graupel
      lam = (n0*params%gam_1_mu_p1/params%gam_1_mu*m1**(-1.0))**(params%inv_p1)
    END IF ! id / l_kfsm

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_lam_n0_1M_KF

  ! Get n0 given a moment and lamda and mu
  SUBROUTINE get_n0(m, p, mu, lam, n0)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_N0'

    REAL(wp), INTENT(IN) :: m, mu, lam
    REAL(wp), INTENT(IN) :: p
    REAL(wp), INTENT(OUT) :: n0

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    n0=m/(GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_n0

  FUNCTION moment(n0,lam,mu,p)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='MOMENT'

    REAL(wp), INTENT(IN) :: n0, lam, mu, p
    REAL(wp) :: moment

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    moment=n0*GammaFunc(1+mu+p)*lam**(-p)/GammaFunc(1+mu)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION moment
END MODULE lookup
