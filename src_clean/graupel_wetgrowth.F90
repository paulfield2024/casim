! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Graupel wet growth and shedding (wetgrowth).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE graupel_wetgrowth
  USE variable_precision, ONLY: wp
  USE process_routines, ONLY: process_rate, process_name,   &
       i_gshd, i_gacw, i_gacr, i_gaci, i_gacs
  USE passive_fields, ONLY: TdegC,  qws0, rho

  USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,   &
       rain_params, T_hom_freeze, DR_melt
  USE mphys_switches, ONLY: i_qv, l_kfsm, l_gamma_online, l_prf_cfrac, i_cfg
  USE mphys_constants, ONLY: Lv, Lf, Ka, Cwater, Cice, Dv
  USE thresholds, ONLY: thresh_small, cfliq_small

  USE ventfac, ONLY: ventilation_3M, ventilation_1M_2M
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='GRAUPEL_WETGROWTH'

  PUBLIC wetgrowth
CONTAINS
  SUBROUTINE wetgrowth(ixy_inner, nz, l_Tcold, qfields, cffields, procs, l_sigevap)
    !< Subroutine to determine if all liquid accreted by graupel can be
    !< frozen or if there will be some shedding
    !<
    !< CODE TIDYING: Should move efficiencies into parameters


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:)
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    REAL(wp), INTENT(IN) :: cffields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    LOGICAL, INTENT(IN) :: l_sigevap(:) ! logical to determine significant evaporation

    ! Local variables

    TYPE(process_name) :: iproc ! processes selected depending on
    ! which species we're modifying

    REAL(wp) :: qv
    REAL(wp) :: dmass, dnumber
    REAL(wp) :: mass

    REAL(wp) :: n0, lam, mu, V_x
    REAL(wp) :: pgacw, pgacr, pgaci, pgacs
    REAL(wp) :: Eff, sdryfac, idryfac
    REAL(wp) :: pgaci_dry, pgacs_dry, pgdry
    REAL(wp) :: pgwet       !< Amount of liquid that graupel can freeze withouth shedding
    REAL(wp) :: pgacsum     !< Sum of all graupel accretion terms
    REAL(wp) :: cf_graupel

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='WETGROWTH'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (.NOT. l_sigevap(k) .AND. l_Tcold(k)) THEN 
          IF (l_prf_cfrac) THEN
            IF (cffields(k,i_cfg) > cfliq_small) THEN
               cf_graupel=cffields(k,i_cfg)
            ELSE
               cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
            END IF
          ELSE
            cf_graupel=1.0
          END IF

          mass=qfields(k, graupel_params%i_1m) / cf_graupel

          IF (mass > thresh_small(graupel_params%i_1m)) THEN
             pgacw=procs(graupel_params%i_1m, i_gacw%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             pgacr=procs(graupel_params%i_1m, i_gacr%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             pgaci=procs(graupel_params%i_1m, i_gaci%id)%column_data(k)/cf_graupel !convert to ingraupel rates
             
             IF (l_kfsm) THEN
                pgacs = 0.0 !<KF. Single ice prognostics. No process gacs. Set to zero.>
             ELSE
                pgacs=procs(graupel_params%i_1m, i_gacs%id)%column_data(k) /cf_graupel !convert to ingraupel rates
             END IF ! l_kfsm
             
             ! Factors for converting wet collection efficiencies to dry ones
             ! AH - Changed 1.0 to 0.001 to match RA and GA tests. 0.001 prevents a divide by zero but is equivalent 
             !      to a 0.0 eff. 
             Eff=0.001
             sdryfac=1.0_wp ! min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)
             Eff=0.001
             idryfac=1.0_wp !min(1.0_wp, 0.2*exp(0.08*TdegC(k))/Eff)
             
             pgaci_dry=idryfac*pgaci
             pgacs_dry=sdryfac*pgacs
             
             pgdry=pgacw+pgacr+pgaci_dry+pgacs_dry
             
             pgacsum=pgacw+pgacr+pgaci+pgacs
             
             IF (pgacsum > thresh_small(graupel_params%i_1m)) THEN
                qv=qfields(k, i_qv)
                mass=qfields(k, graupel_params%i_1m) / cf_graupel
                
                n0=dist_n0(k,graupel_params%id)
                mu=dist_mu(k,graupel_params%id)
                lam=dist_lambda(k,graupel_params%id)
                
                IF (l_gamma_online) THEN 
                   CALL ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, graupel_params)
                ELSE 
                   CALL ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, graupel_params)
                END IF

                pgwet=(910.0/graupel_params%density)**0.625*(Lv*Dv*(qws0(k,ixy_inner)-qv)- &
                     Ka*TdegC(k,ixy_inner)/rho(k,ixy_inner))/(Lf+Cwater*TdegC(k,ixy_inner))*V_x *cf_graupel ! grid mean
                pgwet=pgwet+(pgaci+pgacs)*(1.0-Cice*TdegC(k,ixy_inner)/(Lf+Cwater*TdegC(k,ixy_inner)))

                IF (pgdry < pgwet .OR. TdegC(k,ixy_inner) < T_hom_freeze) THEN ! Dry growth, so use recalculated gaci, gacs
                   IF (pgaci + pgacs > 0.0) THEN

                      IF ( .NOT. l_kfsm ) THEN
                         ! KF. Single ice prognostics. No process gacs
                         dmass=pgacs_dry * cf_graupel !convert back to gridbox mean
                         procs(graupel_params%i_1m, i_gacs%id)%column_data(k)=dmass
                         procs(snow_params%i_1m, i_gacs%id)%column_data(k)=-dmass
                         IF (snow_params%l_2m) THEN
                            dnumber=procs(snow_params%i_2m, i_gacs%id)%column_data(k)*sdryfac
                            procs(snow_params%i_2m, i_gacs%id)%column_data(k)=dnumber
                         END IF
                      END IF ! l_kfsm
                      
                      dmass=pgaci_dry * cf_graupel !convert back to gridbox mean
                      procs(graupel_params%i_1m, i_gaci%id)%column_data(k)=dmass
                      procs(ice_params%i_1m, i_gaci%id)%column_data(k)=-dmass
                      IF (ice_params%l_2m) THEN
                         dnumber=procs(ice_params%i_2m, i_gaci%id)%column_data(k)*idryfac
                         procs(ice_params%i_2m, i_gaci%id)%column_data(k)=dnumber
                      END IF
                   END IF
                ELSE ! Wet growth mode so recalculate gacr, gshd
                   IF (abs(procs(graupel_params%i_1m, i_gacr%id)%column_data(k)) > &
                        thresh_small(graupel_params%i_1m)) THEN
                      dmass=(pgacr+pgwet-pgacsum) * cf_graupel !convert back to gridbox mean
                      IF (dmass > 0.0) THEN
                         procs(graupel_params%i_1m, i_gacr%id)%column_data(k)=dmass
                         procs(rain_params%i_1m, i_gacr%id)%column_data(k)=-dmass
                      ELSE
                         iproc=i_gshd
                         dmass=min(pgacw * cf_graupel, -1.0*dmass)
                         procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                         IF (rain_params%l_2m) THEN
                            dnumber=dmass/(rain_params%c_x*DR_melt**3)
                            procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                         END IF
                      END IF
                   END IF
                END IF
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE wetgrowth
END MODULE graupel_wetgrowth
