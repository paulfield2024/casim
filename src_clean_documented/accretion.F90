! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Rain-cloud accretion (racw) and related warm-rain accretion process rates.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Accretion (collection) of cloud liquid water by rain drops.
!
! Paper reference:
!   Field et al. (2023) Appendix A.5.1, Eq. (A10)-(A11):
!   Pracw = 67(qw*qr)^1.15 following Khairoutdinov & Kogan (2000);
!   Nracw = Pracw*qw/nw.
!
MODULE accretion
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho
! use mphys_switches, only: i_m3r, l_3mr
  USE mphys_switches, ONLY: i_ql, i_qr, i_nl, l_2mc, &
       l_aacc, i_am4, i_am5, l_process, active_rain, isol, l_preventsmall, &
       l_prf_cfrac, i_cfl, i_cfr, l_kk00
  USE mphys_constants, ONLY: fixed_cloud_number
  USE mphys_parameters, ONLY: hydro_params
! use mphys_parameters, only: p1, p2, p3, rain_params
  USE process_routines, ONLY: process_rate, i_pracw, i_aacw
  USE thresholds, ONLY: ql_small, qr_small, cfliq_small
  USE sweepout_rate, ONLY: sweepout
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
! use m3_incs, only: m3_inc_type2
  USE casim_stph, ONLY: l_rp2_casim, fixed_cloud_number_rp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ACCRETION'

  PRIVATE

  PUBLIC racw
CONTAINS
  !---------------------------------------------------------------------------
  !> @author
  !> Ben Shipway
  !
  !> @brief
  !> This subroutine calculates increments due to the accretion of
  !> cloud water by rain
  !--------------------------------------------------------------------------- !
  ! Pracw/Nracw: Khairoutdinov & Kogan (2000) accretion,
  ! Field et al. (2023) Eq. (A10)-(A11).
  SUBROUTINE racw(ixy_inner, dt, qfields, cffields, aerofields, procs, params, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt !< microphysics time increment (s)   
                         !! dt NEEDED for 3rd moment code
    REAL(wp), INTENT(IN) :: qfields(:,:)     !< hydrometeor fields
    REAL(wp), INTENT(IN) :: cffields(:,:)     ! < cloud fractions
    REAL(wp), INTENT(IN) :: aerofields(:,:)  !< aerosol fields
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)         !< hydrometeor process rates
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:) !< aerosol process rates
    TYPE(hydro_params), INTENT(IN) :: params !< parameters describing hydrometor size distribution/fallspeeds etc.

    ! Local Variables

    REAL(wp) :: dmass, dnumber, damass
!   real(wp) :: m1, m2, m3, dm1, dm2, dm3

    REAL(wp) :: cloud_mass
    REAL(wp) :: cloud_number
    REAL(wp) :: rain_mass
!   real(wp) :: rain_number
!   real(wp) :: rain_m3
    REAL(wp) :: cf_liquid, cf_rain


    REAL(wp) :: mu, n0, lam
    LOGICAL :: l_kk_acw=.TRUE.

    INTEGER :: k ! local index for k

    CHARACTER(len=*), PARAMETER :: RoutineName='RACW'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    IF ( l_rp2_casim ) THEN
        fixed_cloud_number = fixed_cloud_number_rp
    END IF

    DO k = 1, ubound(qfields,1)
       IF (l_prf_cfrac) THEN
          IF (cffields(k,i_cfl) > cfliq_small) THEN
             ! only doing liquid cloud fraction at the moment
             cf_liquid=cffields(k,i_cfl)
          ELSE
             cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
          END IF
          IF (cffields(k,i_cfr) > cfliq_small) THEN
             ! only doing liquid cloud fraction at the moment
             cf_rain=cffields(k,i_cfr)
          ELSE
             cf_rain=cfliq_small !nonzero value - maybe move cf test higher up
          END IF
       ELSE
          cf_liquid=1.0_wp
          cf_rain=1.0_wp
       END IF
       
       cloud_mass = qfields(k, i_ql) / cf_liquid
       rain_mass = qfields(k, i_qr) / cf_rain
       
       IF (l_2mc ) THEN
          cloud_number=qfields(k, i_nl) / cf_liquid
       ELSE
          cloud_number=fixed_cloud_number / cf_liquid !check that fixed_cloud_number is grid mean
       END IF
       
    
       ! if (l_3mr) rain_m3 = qfields(k, i_m3r)
       
       IF (cloud_mass*cf_liquid > ql_small .AND. rain_mass*cf_rain > qr_small) THEN
          IF (l_kk_acw) THEN
             !        dmass=min(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
             IF (l_kk00) THEN
                ! Use KK accretion parametrisation but limit to 90% of cloud mass removal
                dmass = MIN(0.9*cloud_mass, 67.0*(cloud_mass*rain_mass)**1.15)
             ELSE
                ! Use Kogan(2013) accretion parametrisation but limit to 90% 
                ! of cloud mass removal
                dmass = min(0.9*cloud_mass, 8.53*(cloud_mass**1.05)*(rain_mass)**0.98)
             END IF
             
          ELSE
             n0=dist_n0(k,params%id)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)
             dmass=sweepout(n0, lam, mu, params, rho(k,ixy_inner))*cloud_mass
          END IF
          
          IF (l_preventsmall .AND. dmass < qr_small) dmass=0.0
          IF (l_2mc) dnumber=dmass/(cloud_mass/cloud_number)
          
          
          IF (l_prf_cfrac) THEN
             ! convert back to grid mean
             dmass=dmass*min(cf_liquid, cf_rain)
             dnumber=dnumber*min(cf_liquid, cf_rain)
             cloud_mass=cloud_mass*cf_liquid
          END IF
          
          
          procs(i_ql, i_pracw%id)%column_data(k)=-dmass
          procs(i_qr, i_pracw%id)%column_data(k)=dmass
          IF (l_2mc) THEN
             procs(i_nl, i_pracw%id)%column_data(k)=-dnumber
          END IF
       
      ! if (l_3mr) then
      !    m1=rain_mass/rain_params%c_x
      !    m2=rain_number
      !    m3=rain_m3
             
      !    dm1=dt*dmass/rain_params%c_x
      !    dm2=0
      !    call m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)
      !    dm3=dm3/dt
      !    procs(i_m3r, i_pracw%id)%column_data(k) = dm3
      ! end if

    
       END IF
    END DO

    IF (l_aacc .AND. l_process) THEN
       IF (active_rain(isol)) THEN
          DO k = 1, ubound(qfields,1)
             dmass = procs(i_qr, i_pracw%id)%column_data(k)
             damass=dmass/cloud_mass*aerofields(k,i_am4)
             aerosol_procs(i_am4, i_aacw%id)%column_data(k)=-damass
             aerosol_procs(i_am5, i_aacw%id)%column_data(k)=damass
          END DO
       END IF
    END IF
          
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE racw
END MODULE accretion
