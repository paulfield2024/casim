! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Shipway (2015) revised droplet activation scheme (peak supersaturation approximation).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE shipway_activation_mod
  USE variable_precision, ONLY: wp

  ! A revised activation scheme
  ! See Shipway, B. J.: Revisiting Twomey's approximation for peak supersaturation,
  ! Atmos. Chem. Phys., 15, 3803-3814, doi:10.5194/acp-15-3803-2015, 2015.

  USE shipway_parameters
  USE shipway_constants
  USE shipway_lookup, ONLY: lookup_I, xmax, ymax, xmin, ymin
  USE mphys_die, ONLY: throw_mphys_error, incorrect_opt, mphys_message

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SHIPWAY_ACTIVATION_MOD'

  REAL(wp) :: dv_flag=0 ! Determines method for calculating Dv
  REAL(wp) :: est_smax_time
  REAL(wp) :: psi1 ! placed here so it can be used for estimating smax time
  REAL(wp) :: C1=1.058, C2=1.904 ! See Korolev et al. paper

  REAL(wp) :: wdiag, Tdiag, Ndiag
  INTEGER :: counter=0

  REAL(wp) :: cumulative_xmin=1.e10
  REAL(wp) :: cumulative_xmax=0.0
  INTEGER :: cumulative_brent=0
  INTEGER :: cumulative_quad=0
  INTEGER :: cumulative_secant=0
  INTEGER :: cumulative_mid=0

  REAL :: LHS_hold

  CHARACTER(len=600) :: std_msg

  PUBLIC solve_nccn_household,solve_nccn_brent

CONTAINS

  SUBROUTINE calc_nccn(smax, nccn, nccni_dg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: smax
    REAL(wp), INTENT(OUT) :: nccn
    REAL(wp), INTENT(OUT), OPTIONAL :: nccni_dg(max_nmodes)

    ! Local variables
    REAL(wp) :: nccni, erfarg
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: RoutineName='CALC_NCCN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    nccn=0.0
    DO i=1,nmodes
      IF (use_mode(i)) THEN
        erfarg=log(smax/s0i(i))/(sqrt(2.)*log(sigmas(i)))
        nccni = 0.5*Ndi(i)*(1.0+erf(erfarg))

        nccn=nccn+nccni
        nccni_dg(i)=nccni
      END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE calc_nccn


  SUBROUTINE calc_LHS(w, T, pressure, alpha_c, entrain_fraction, LHSout)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: T ! Temperature (K)
    REAL(wp), INTENT(IN) :: w ! vertical velocity (m/s)
    REAL(wp), INTENT(IN) :: pressure ! Pressure (Pa)
    REAL(wp), INTENT(IN) :: alpha_c ! Condensation coefficient
    REAL(wp), INTENT(IN) :: entrain_fraction ! Entrainment_fraction
    REAL(wp), INTENT(OUT) :: LHSout

    ! Local variables

    REAL(wp) :: LHS ! intermediate value

    REAL(wp) :: Tc ! Temperature (C)
    REAL(wp) :: es ! Saturation vapour pressure

    REAL(wp) :: psi2, G, Dv_here

    REAL(wp) :: Lvt ! Temperature dependent Lv

    CHARACTER(len=*), PARAMETER :: RoutineName='CALC_LHS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Tc = T-273.15

    es = (100.*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))

    LvT = Lv -(4217.4-1870.0)*Tc

    psi1 = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-entrain_fraction)
    psi2 = eps*pressure/es+LvT**2/(Rv*T**2*cp)

    IF (int(dv_flag)==1)THEN
       Dv_here=Dv
    ELSE IF(dv_flag <= 0)THEN
       Dv_here=Dv_mean(T,alpha_c)
    ELSE
       Dv_here=dv_flag
    END IF

    G = 1/(rhow*(Rv*T/(es*Dv_here)+LvT*(LvT/(Rv*T)-1)/(ka*T)))

    LHS = sqrt(2.0*pi)*rho*(psi1*w)**1.5/(4*pi*rhow*psi2*G**1.5)

    LHSout=LHS

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE calc_LHS

  SUBROUTINE set_inputs(s,sp,sm,ds)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Routine to choose sm,s,sp and ensure they sit within (xminx, xmax)
    REAL(wp), INTENT(INOUT) :: s, sp, sm, ds

    ! Local variables
    LOGICAL :: adjust
    REAL(wp) :: s0ratioi_p, s0ratioi_m
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: RoutineName='SET_INPUTS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ds=min(0.001, .5*abs(s))
    sp = s + ds
    sm = s - ds

    adjust=.TRUE.
    DO WHILE (adjust)
      DO i=1,nmodes
        IF (use_mode(i)) THEN
          s0ratioi_p = sp/s0i(i)
          s0ratioi_m = sm/s0i(i)
          IF (s0ratioi_p > xmax)THEN
            WRITE (std_msg, *) 'Out of range X too big: ', i, xmin, s0i(i), &
                               s0ratioi_p, xmax
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            ! reset
            WRITE (std_msg, *) 'old sm,s,sp', sm,s,sp
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            sp = xmax*s0i(i)*.99
            s  = sp - ds
            sm = s - ds
            WRITE (std_msg, *) 'new sm,s,sp', sm,s,sp
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            EXIT
          END IF
          IF (s0ratioi_m < xmin)THEN
            WRITE (std_msg, *)'Out of range X too small: ', i, xmin, s0i(i), &
                              s0ratioi_m, xmax
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            ! reset
            WRITE (std_msg, *)'old sm,s,sp', sm,s,sp
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            sm = xmin*s0i(i)*1.01
            ds=min(0.001, sm)
            s  = sm + ds
            sp = s + ds
            WRITE (std_msg, *)'new sm,s,sp', sm,s,sp
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
            EXIT
          END IF
        END IF
        IF (i==nmodes)adjust=.FALSE.
      END DO
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE set_inputs


  SUBROUTINE set_inputs_safe(sp,sm)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Routine to choose sm,sp and ensure they sit within (xminx, xmax)
    REAL(wp), INTENT(INOUT) :: sp, sm

    ! Local variables
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: RoutineName='SET_INPUTS_SAFE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    sm=0.0
    sp=HUGE(sp)
    DO i=1,nmodes
      IF (use_mode(i)) THEN
        sm = max(xmin*s0i(i)*1.01, sm)
        sp = min(xmax*s0i(i)*0.99, sp)

      END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE set_inputs_safe

  SUBROUTINE get_extent(s, sp, sm)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Routine to choose sm,sp and ensure they sit within (xminx, xmax)
    REAL(wp), INTENT(IN) :: s
    REAL(wp), INTENT(OUT) :: sp, sm

    ! Local variables

    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: RoutineName='GET_EXTENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    sm=0.0
    sp=HUGE(sp)
    DO i=1,nmodes
      IF (use_mode(i)) THEN
        sm = max(s/s0i(i), sm)
        sp = min(s/s0i(i), sp)
      END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE get_extent

  SUBROUTINE calc_RHS(s,RHS)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: s
    REAL(wp), INTENT(OUT) :: RHS

    ! Local variables

    REAL(wp) :: s0ratioi
    INTEGER :: i

    REAL(wp) :: J1

    CHARACTER(len=*), PARAMETER :: RoutineName='CALC_RHS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    RHS=0.0
    DO i=1,nmodes
      IF (use_mode(i)) THEN
        s0ratioi = s/s0i(i)
        IF (logsigmas(i) > ymax .OR. logsigmas(i)<ymin)THEN
          WRITE (std_msg,*) 'Out of range Y: ', ymin, logsigmas(i), ymax, &
                             s0ratioi
          CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                 std_msg)
        END IF
      END IF
    END DO
    DO i=1,nmodes
      IF (use_mode(i)) THEN
        CALL lookup_I(s0ratioi, logsigmas(i), J1)
        RHS = RHS + ai(i)*J1
      END IF
    END DO

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE calc_RHS

  SUBROUTINE calc_KC(T)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: T ! Temperature(K)
    REAL(wp) :: Ak
    INTEGER :: i

    CHARACTER(len=*), PARAMETER :: RoutineName='CALC_KC'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Ak = 2.*Mw*zetasa/(Ru*T*rhow)

    DO i=1,nmodes
      IF (use_mode(i)) THEN
       ! These are defined in Khvorostyanov and Curry (2006)
       s0i(i) = rdi(i)**(-(1.+betai(i)))*sqrt(4*Ak**3/(27*bi(i)))
       sigmas(i) = sigmad(i)**(1.+betai(i))

       ai(i) = Ndi(i)*s0i(i)*s0i(i)/log(sigmas(i))
       logsigmas(i) = log(sigmas(i))
     END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE calc_KC

  SUBROUTINE solve_nccn_household(order,niter,sa_in,  &
     w,T,pressure,alpha_c,ent_fraction,               &
     smax,nccn,nccni)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE
    ! Use householders methods (e.g. Newton, Halley) to find roots.

    INTEGER, INTENT(IN) :: order, niter
    REAL(wp), INTENT(IN)    :: sa_in
    REAL(wp), INTENT(IN)    :: w,T,pressure, alpha_c, ent_fraction
    REAL(wp), INTENT(OUT)   :: smax, nccn
    REAL(wp), INTENT(OUT)   :: nccni(3) ! diagnostic for each mode

    REAL(wp) :: sa
    INTEGER :: it
    INTEGER :: maxiter=10
    INTEGER :: bigiter
    REAL(wp) :: RHSout,RHSout_p1,RHSout_m1, ds
    REAL(wp) :: F, dF, d2F, diff
    REAL(wp) :: LHS, sa_p1, sa_m1, sa_tm1

    REAL(wp) :: tolx=1e-6, toly ! This tolerence only necessary for very high numbers
    CHARACTER(len=*), PARAMETER :: RoutineName='SOLVE_NCCN_HOUSEHOLD'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (.NOT. any(use_mode)) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN  ! If there's no aerosol present in any of the modes
    END IF

    CALL calc_KC(T)
    CALL calc_LHS(w, T, pressure, alpha_c, ent_fraction, LHS)
    LHS_hold = LHS
    wdiag = w
    Tdiag = T
    Ndiag = sum(Ndi)

    !toly=LHS*.1 ! 10 percent error
    toly=LHS*1e5 ! Essentially not used

    sa = sa_in
    sa_tm1=sa

    diff=1000.
    RHSout=LHS+LHS
    counter=counter+1
    bigiter=0
    DO WHILE((abs(diff) > tolx .OR. abs(RHSout-LHS)>toly) .AND. bigiter<maxiter)
      DO it=1,niter
        CALL set_inputs(sa, sa_p1, sa_m1, ds)
        sa_tm1=sa
        CALL calc_RHS(sa,RHSout)
        CALL calc_RHS(sa_p1,RHSout_p1)

        F=RHSout-LHS
        dF=(RHSout_p1-RHSout)/ds

        SELECT CASE(order)
        CASE(1)
          diff=-F/dF
        CASE(2)
          CALL calc_RHS(sa_m1,RHSout_m1)
          d2F=(RHSout_p1-2*RHSout + RHSout_m1)/(ds*ds)
          diff=-2.*F*dF/(2.*dF*dF-F*d2F)
        END SELECT
        sa = sa + diff

      END DO

       bigiter=bigiter+1
    END DO

    smax=sa
    est_smax_time=1./((psi1*w*C1/C2)/smax)

    CALL calc_nccn(smax, nccn, nccni)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

    END SUBROUTINE solve_nccn_household

  SUBROUTINE solve_nccn_brent( &
     w,T,pressure,alpha_c,ent_fraction,               &
     smax,nccn,nccni)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN)    :: w,T,pressure, alpha_c, ent_fraction
    REAL(wp), INTENT(OUT)   :: smax, nccn
    REAL(wp), INTENT(OUT)   :: nccni(3) ! diagnostic for each mode

    REAL(wp) :: sa, sb
    REAL(wp) :: LHS, sa_p1, sa_m1
    INTEGER :: iflag  ! flag to indicate if solution found
    INTEGER :: ibrent ! indicates number of iteration in bre

    REAL(wp) :: tolx=1e-4           ! tolerence for brent
    REAL(wp) :: tolf=1e-1            ! tolerence for brent
    LOGICAL :: verbose=.FALSE.  ! set brent to verbose output

    REAL(wp) :: ds, RHS
!    real(wp) :: J1
    INTEGER :: nscan=101, i, iquad, isecant, imidpoint
    CHARACTER(len=500) :: message

    CHARACTER(len=*), PARAMETER :: RoutineName='SOLVE_NCCN_BRENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (.NOT. any(use_mode)) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN  ! If there's no aerosol present in any of the modes
    END IF

    CALL calc_KC(T)
    CALL calc_LHS(w, T, pressure, alpha_c, ent_fraction, LHS)

    LHS_hold = LHS
    wdiag = w
    Tdiag = T
    Ndiag = sum(Ndi)
    ibrent = 0
    RHS = 0.0

    CALL set_inputs_safe(sa_p1, sa_m1)
    sa = sa_m1
    sb = sa_p1

    CALL brent(smax_eq, sa, sb, smax, iflag, ibrent &
       ,tolx, tolf, verbose,iquad, isecant, imidpoint)

    IF (iflag == -1) THEN
      IF (smax_eq(sb) - LHS < 0)THEN
        ! In this case we've probably exceeded the bounds of the lookup table
        ! but we should be able to extrapolate using the asymptotic limit of the
        ! integral, i.e. lim x->inf I(x,ls) = x**2*I(xmax,ls)/xmax**2
        DO i=1,nmodes
          IF (use_mode(i)) THEN
!            call lookup_I(xmax, logsigmas(i), J1)
!            RHS = RHS + (ai(i)/s0i(i)/s0i(i))*(J1/xmax/xmax)
            ! NB I've left this in a long form for comparison with
            ! the integral form, but this simplifies considerable
            ! to give smax ~ w^3/4 N^-1/2
            RHS = RHS + (ai(i)/s0i(i)/s0i(i))*2*sqrt(pi)*logsigmas(i)
          END IF
        END DO
        smax=sqrt(LHS/RHS)
      ELSE
        IF (verbose)THEN
          DO i=1,nmodes
            WRITE (std_msg, *) i, smax/s0i(i)
            CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
          END DO
          WRITE (std_msg, *) w,T, Ndi(:),logsigmas(:), LHS, smax
          CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
          ds=(sb-sa)/nscan
          DO i=1,nscan-1
            sa=sa+ds
            WRITE(*,'(e18.4,A,e18.4,A)') sa, ',', smax_eq(sa),','
          END DO
        END IF
        WRITE(message,'(A,e18.4,A,e18.4,A,e18.4,A,e18.4,A,e18.4,A,e18.4)') &
           'w=',w,'; p=',pressure,'; T=',T,'; N=',Ndiag,'; sm=',sa_m1,'; sp=',sa_p1
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
             'Beyond scope of numerics of the parametrization:'//trim(message))
      END IF
    END IF

    IF (verbose)THEN
      WRITE (std_msg, *) 'Number of brent calls:', ibrent
      CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
      cumulative_brent=cumulative_brent + ibrent
      cumulative_quad=cumulative_quad + iquad
      cumulative_secant=cumulative_secant + isecant
      cumulative_mid=cumulative_mid + imidpoint
      WRITE (std_msg, *) 'Cumulative brent calls:', cumulative_brent, &
         cumulative_quad, cumulative_secant, cumulative_mid
      CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
      CALL get_extent(smax, sa, sb)
      WRITE (std_msg, *) 'Extrema of ratios:', sa, sb
      CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
      cumulative_xmin=min(sa, cumulative_xmin)
      cumulative_xmax=max(sb, cumulative_xmax)
      WRITE (std_msg, *) 'Cumulative extrema', cumulative_xmin, cumulative_xmax
      CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
    END IF

    est_smax_time=1.0/((psi1*w*C1/C2)/smax)
    CALL calc_nccn(smax, nccn, nccni)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE solve_nccn_brent

  FUNCTION smax_eq(sa)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(INOUT) :: sa
    REAL(wp) :: smax_eq

    REAL(wp) :: RHSout

    CHARACTER(len=*), PARAMETER :: RoutineName='SMAX_EQ'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL calc_RHS(sa,RHSout)

    smax_eq = RHSout - LHS_hold

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION smax_eq


  SUBROUTINE brent(f, a, b, x, iflag, ibrent, tolx, tolf, verbose, &
     iquad, isecant, imidpoint)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp) :: f
    REAL(wp), INTENT(INOUT) :: a, b
    REAL(wp), INTENT(OUT) :: x
    INTEGER, INTENT(OUT) :: &
         iflag  & ! =-1 if f(a)f(b) > 0
                  ! =0 if successful
         , ibrent ! number of calls to function f (this does not
                  ! include calls made for verbose output.)

    REAL(wp), INTENT(INOUT) :: tolx,tolf
    LOGICAL, INTENT(INOUT) :: verbose
    INTEGER, OPTIONAL :: iquad, isecant, imidpoint

    REAL(wp) :: c, d, tmp, s
    LOGICAL :: mflag
    INTEGER :: icount
    INTEGER, PARAMETER :: icountmax=50

    REAL(wp) :: pfs, fa, fs, fb, fc

    CHARACTER(len=*), PARAMETER :: RoutineName='BRENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    iflag=-999
    ibrent=0
    iquad=0
    isecant=0
    imidpoint=0

    fa=f(a)
    ibrent=ibrent+1
    fb=f(b)
    ibrent=ibrent+1

    IF (verbose)THEN
       WRITE (std_msg, *) 'a = ', a, ', b = ', b
       CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
       WRITE (std_msg, *) 'F(a) = ', fa, ', F(b) = ', fb
       CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
    END IF

    IF (fa == 0)THEN
       x=a
       iflag=0
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       RETURN
    END IF

    IF (fb == 0)THEN
       x=b
       iflag=0
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       RETURN
    END IF

    IF ( fa*fb > 0 )THEN
       x=-999
       iflag=-1
       IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
       RETURN
    END IF

    IF (abs(fa) < abs(fb))THEN
       c=b
       fc=fb
       b=a
       fb=fa
       a=c
       fa=fc
    ELSE
       c=a
       fc=fa
    END IF

    mflag=.TRUE.
    d=0
    fs=tolf+1.

    icount=0
    DO WHILE (fb /=0 .AND. abs(b-a) > tolx .AND. abs(fs) > tolf &
         .AND. icount < icountmax )
       icount=icount+1
       IF ((fa /= fc) .AND. (fb /= fc)) THEN
          !inverse quadratic interpolation
          s=a*fb*fc/((fa-fb)*(fa-fc)) &
               + b*fa*fc/((fb-fa)*(fb-fc)) &
               + c*fa*fb/((fc-fa)*(fc-fb))
          IF (verbose)THEN
             pfs=f(s)
             WRITE (std_msg, *) 'inverse quadratic ','s=',s,'f(s)=', pfs
             CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
          END IF
          iquad=iquad+1
       ELSE
          ! secant method
          s=b-fb*(b-a)/(fb-fa)
          IF (verbose)THEN
             pfs=f(s)
             WRITE (std_msg, *) 'secant ','s=',s,'f(s)=', pfs
             CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
          END IF
          isecant=isecant+1
       END IF
       IF ((s < b .AND. s <(3*a+b)/4) &
            .OR. (s > b .AND. s >(3*a+b)/4) &
            .OR. (mflag .AND. abs(s-b)>=abs(b-c)/2.) &
            .OR. ((.NOT.mflag) .AND. (abs(s-b)>=abs(c-d)/2.)))THEN
          ! midpoint
          s=(a+b)/2
          IF (verbose)THEN
             pfs=f(s)
             WRITE (std_msg, *) 'midpoint ','s=',s,'f(s)=', pfs
             CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
          END IF
          imidpoint=imidpoint+1
          mflag=.TRUE.
       ELSE
          mflag=.FALSE.
       END IF
       fs=f(s)
       ibrent=ibrent+1
       d=c
       c=b
       IF (fa*fs<0)THEN
          b=s
          fb=fs
       ELSE
          a=s
          fa=fs
       END IF


       IF (abs(fa) < abs(fb))THEN
          tmp=b
          b=a
          a=tmp
          tmp=fb
          fb=fa
          fa=tmp
       END IF

       IF (verbose)THEN
         WRITE (std_msg, *) '|a-b|=',abs(a-b)
         CALL mphys_message(ModuleName//':'//RoutineName, std_msg)
       END IF
    END DO


    IF (icount==icountmax)THEN
       x=-999
       iflag=-2
    ELSE
       x=b
       iflag=0
    END IF


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    RETURN
  END SUBROUTINE brent

END MODULE shipway_activation_mod
