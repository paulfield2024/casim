! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Vapour deposition growth of ice, snow and graupel (idep).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE ice_deposition
  USE variable_precision, ONLY: wp, iwp
  USE passive_fields, ONLY: rho, pressure, exner, TdegK
  USE mphys_switches, ONLY: i_qv, i_th   &
       , i_am6, i_an2  &
       , i_am7, i_an6                   &
       , l_process, l_passivenumbers_ice, l_passivenumbers &
       , i_an12 &
       , i_am8, i_am2, i_an11, l_gamma_online
  USE type_process, ONLY: process_name
  USE process_routines, ONLY: process_rate, i_idep,    &
       i_dsub, i_sdep, i_gdep, i_dssub, i_dgsub &
       , i_isub, i_ssub, i_gsub, i_iacw, i_raci, i_sacw, i_sacr  &
       , i_gacw, i_gacr
  USE mphys_parameters, ONLY: hydro_params, rain_params, cloud_params
  USE mphys_constants, ONLY: Ls, Lf, ka, Dv, Rv
  USE qsat_funs, ONLY: qisaturation
  USE thresholds, ONLY: thresh_small
  USE aerosol_routines, ONLY: aerosol_active

  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
  USE ventfac, ONLY: ventilation_1M_2M, ventilation_3M
! removing line below changes answers
  USE special, ONLY: pi

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ICE_DEPOSITION'

  LOGICAL :: l_latenteffects = .FALSE.

  PUBLIC idep
CONTAINS

  !< Subroutine to determine the deposition/sublimation onto/from
  !< ice, snow and graupel.  There is no source/sink for number
  !< when undergoing deposition, but there is a sink when sublimating.
  SUBROUTINE idep(ixy_inner, dt, nz, l_Tcold, params, qfields, cffields, procs, dustact, aeroice, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfi, i_cfs, i_cfg


    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:) 
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN) :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
    REAL(wp), INTENT(IN) :: cffields(:,:)

    ! aerosol fields
    TYPE(aerosol_active), INTENT(IN) :: dustact(:), aeroice(:)

    ! optional aerosol fields to be processed
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)


    ! Local Variables
    TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    TYPE(process_name) :: i_acw, i_acr ! collection processes with cloud and rain
    REAL(wp) :: dmass, dnumber, dmad, dnumber_a, dnumber_d, dmac

    REAL(wp) :: th
    REAL(wp) :: qv
    REAL(wp) :: num, mass

    REAL(wp) :: qis(nz)
    REAL(wp) :: n0, lam, mu
    REAL(wp) :: V_x, AB
    REAL(wp) :: cf

    LOGICAL :: l_suball

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='IDEP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
      IF (l_Tcold(k)) THEN
        th=qfields(k, i_th)
        qis(k) = qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
      END IF
    END DO

    DO k = 1, nz
      IF (l_Tcold(k)) THEN 
        l_suball=.FALSE. ! do we want to sublimate everything?
          
        mass=qfields(k, params%i_1m)
          
        qv=qfields(k, i_qv)
          
        IF (qv>qis(k)) THEN
          SELECT CASE (params%id)
          CASE (3_iwp) !ice
            iproc=i_idep
            iaproc=i_dsub
            i_acw=i_iacw
            i_acr=i_raci
            cf=cffields(k,i_cfi)
          CASE (4_iwp) !snow
            iproc=i_sdep
            iaproc=i_dssub
            i_acw=i_sacw
            i_acr=i_sacr
            cf=cffields(k,i_cfs)
          CASE (5_iwp) !graupel
            iproc=i_gdep
            iaproc=i_dgsub
            i_acw=i_gacw
            i_acr=i_gacr
            cf=cffields(k,i_cfg)
          END SELECT
        ELSE
          SELECT CASE (params%id)
          CASE (3_iwp) !ice
            iproc=i_isub
            iaproc=i_dsub
            i_acw=i_iacw
            i_acr=i_raci
            cf=cffields(k,i_cfi)
          CASE (4_iwp) !snow
            iproc=i_ssub
            iaproc=i_dssub
            i_acw=i_sacw
            i_acr=i_sacr
            cf=cffields(k,i_cfs)
          CASE (5_iwp) !graupel
            iproc=i_gsub
            iaproc=i_dgsub
            i_acw=i_gacw
            i_acr=i_gacr
            cf=cffields(k,i_cfg)
          END SELECT
        END IF
          
        IF (mass > thresh_small(params%i_1m)) THEN ! if no existing ice, we don't grow/deplete it.
             
          IF (params%l_2m) num=qfields(k, params%i_2m)
             
          n0=dist_n0(k,params%id)
          mu=dist_mu(k,params%id)
          lam=dist_lambda(k,params%id)
            
          IF (l_gamma_online) THEN
            CALL ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, params)
          ELSE
            CALL ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, params)
          END IF
             
          AB=1.0/(Ls*Ls/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))*rho(k,ixy_inner)+1.0/(Dv*qis(k)))
          dmass=(qv/qis(k)-1.0)*V_x*AB *cf ! grid mean
             
          ! Include latent heat effects of collection of rain and cloud
          ! as done in Milbrandt & Yau (2005)
          IF (l_latenteffects) THEN
            dmass=dmass - Lf*Ls/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))                      &
                  *(procs(cloud_params%i_1m, i_acw%id)%column_data(k)          &
                  + procs(rain_params%i_1m, i_acr%id)%column_data(k))
          END IF

          ! Check we haven't become subsaturated and limit if we have (dep only)
          ! NB doesn't account for simultaneous ice/snow growth - checked elsewhere
          IF (dmass > 0.0) dmass=min((qv-qis(k))/dt,dmass)
          ! Check we don't remove too much (sub only)
          IF (dmass < 0.0) dmass=max(-mass/dt,dmass)
             
          IF (params%l_2m) THEN
            dnumber=0.0
            IF (dmass < 0.0) dnumber=dmass*num/mass
          END IF
             
          IF (-dmass*dt >0.98*mass .OR. (params%l_2m .AND.                     &
              -dnumber*dt > 0.98*num)) THEN
            l_suball=.TRUE.
            dmass=-mass/dt
            dnumber=-num/dt
          END IF
             
          procs(i_qv, iproc%id)%column_data(k)=-dmass
          procs(params%i_1m, iproc%id)%column_data(k)=dmass
             
          IF (params%l_2m) procs(params%i_2m,iproc%id)%column_data(k)=dnumber
             
          IF (dmass < 0.0 .AND. l_process) THEN ! Only process aerosol if sublimating
            IF (iaproc%id==i_dsub%id) THEN
              dmad=dnumber*dustact(k)%mact1_mean*dustact(k)%nratio1
              dnumber_d=dnumber*dustact(k)%nratio1
            ELSE IF (iaproc%id==i_dssub%id) THEN
              dmad=dnumber*dustact(k)%mact2_mean*dustact(k)%nratio2
              dnumber_d=dnumber*dustact(k)%nratio2
            ELSE IF (iaproc%id==i_dgsub%id) THEN
              dmad=dnumber*dustact(k)%mact3_mean*dustact(k)%nratio3
              dnumber_d=dnumber*dustact(k)%nratio3
            END IF

            !checking that mass change from
            !active insol in ice is negative
            !and larger in magnitude than epsilon
            !for sublimation of ice
            IF (dmad < -epsilon(dmad)) THEN 
              aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
              aerosol_procs(i_am6, iaproc%id)%column_data(k)=-dmad
              ! <WARNING: putting back in coarse mode
              IF (l_passivenumbers_ice) THEN
                aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d 
              END IF
              aerosol_procs(i_an6, iaproc%id)%column_data(k)=-dnumber_d
              ! <WARNING: putting back in coarse mode
            END IF

            IF (iaproc%id==i_dsub%id) THEN
              dmac=dnumber*aeroice(k)%mact1_mean*aeroice(k)%nratio1
              dmac=min(dmac,aeroice(k)%mact1/dt)
              dnumber_a=dnumber*aeroice(k)%nratio1
            ELSE IF (iaproc%id==i_dssub%id) THEN
              dmac=dnumber*aeroice(k)%mact2_mean*aeroice(k)%nratio2
              dmac=min(dmac,aeroice(k)%mact2/dt)
              dnumber_a=dnumber*aeroice(k)%nratio2
            ELSE IF (iaproc%id==i_dgsub%id) THEN
              dmac=dnumber*aeroice(k)%mact3_mean*aeroice(k)%nratio3
              dmac=min(dmac,aeroice(k)%mact3/dt)
              dnumber_a=dnumber*aeroice(k)%nratio3
            END IF

            !checking that mass change from
            !active sol in ice is negative
            !and larger in magnitude than epsilon
            !for sublimation of ice
            IF (dmac < -epsilon(dmac)) THEN 
              aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
              aerosol_procs(i_am2, iaproc%id)%column_data(k)=-dmac
              ! <WARNING: putting back in accumulation mode                  
              IF (l_passivenumbers) THEN
                aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
              END IF
              aerosol_procs(i_an2, iaproc%id)%column_data(k)=-dnumber_a
              ! <WARNING: putting back in accumulation mode
            END IF
          END IF
        END IF
      END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE idep
END MODULE ice_deposition
