! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Cloud-to-rain autoconversion (raut), KK2000/Kogan2013 schemes.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE autoconversion
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho
  USE mphys_switches, ONLY: i_ql, i_qr, i_nl, i_nr, l_2mc, &
       l_2mr, l_aaut, i_am4, i_am5, cloud_params, rain_params, l_process, &
       l_separate_rain, l_preventsmall, l_prf_cfrac, i_cfl, l_kk00
! use mphys_switches, only: m3r, l_3mr
  USE mphys_constants, ONLY: fixed_cloud_number
  USE mphys_parameters, ONLY: rain_params
! use mphys_parameters, only: mu_aut
  USE process_routines, ONLY: process_rate, i_praut, i_aaut
  USE thresholds, ONLY: ql_small, nl_small, qr_small, cfliq_small
! use m3_incs, only: m3_inc_type3
  USE casim_stph, ONLY: l_rp2_casim, fixed_cloud_number_rp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='AUTOCONVERSION'

  PRIVATE

  PUBLIC raut
CONTAINS

  SUBROUTINE raut(ixy_inner, dt, qfields, cffields, aerofields, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:),    cffields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

    ! Local variables
    REAL(wp) :: dmass, dnumber1, dnumber2, damass
!   real(wp) :: dm1,dm2,dm3
    REAL(wp) :: cloud_mass
    REAL(wp) :: cloud_number
!   real(wp) :: p1, p2, p3
!   real(wp) :: k1, k2, k3
    REAL(wp) :: mu_qc ! < cloud shape parameter (currently only used diagnostically here)
    REAL(wp) :: cf_liquid

    INTEGER :: k
    CHARACTER(len=*), PARAMETER :: RoutineName='RAUT'

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

      IF (cffields(k,i_cfl) > cfliq_small) THEN  !only doing liquid cloud fraction at the moment
        cf_liquid=cffields(k,i_cfl)
      ELSE
        cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
      END IF

    ELSE ! l_prf_cfrac

      cf_liquid = 1.0_wp

    END IF

    cloud_mass=qfields(k, i_ql)/cf_liquid

    IF (l_2mc) THEN
      cloud_number=qfields(k, i_nl)/cf_liquid
    ELSE
      cloud_number=fixed_cloud_number
    END IF

    IF (cloud_mass *cf_liquid> ql_small .AND. &
         cloud_number * cf_liquid > nl_small) THEN
      IF (l_kk00) THEN
         dmass = 1350.*cloud_mass**2.47*  &
              (cloud_number/1.e6*rho(k,ixy_inner))**(-1.79)
      ELSE
         ! new method, k13 scheme
         dmass = 7.98e10*cloud_mass**4.22*  &
              (cloud_number/1.e6*rho(k,ixy_inner))**(-3.01)
      END IF

      dmass=min(.25*cloud_mass/dt, dmass)
      IF (l_preventsmall .AND. dmass < qr_small) dmass=0.0
      IF (l_2mc) dnumber1=dmass/(cloud_mass/cloud_number)
      mu_qc=min(15.0_wp, (1000.0E6/cloud_number + 2.0))
      ! The following line is correct for l_kk00 = .true., i.e. use KK2000 for 
      ! autoconversion and accretion, which is the default switch in 
      ! mphys_switches. However, need to confirm this is 
      ! correct for Kogan (k13), i.e. l_kk00 = .false., since if the 50.0E-6
      ! is the notional diameter threshold for autoconversion it needs to 
      ! be changed 
      IF (l_2mr) dnumber2=dmass/(rain_params%c_x*(mu_qc/3.0)*(50.0E-6)**3)

          ! AH - found that at the cloud edges rain number produced 
          !      by autoconversion can be larger than cloud number removed. This
          !      is not physically possible, so limit the dnumber2 to be the 
          !      same as dnumber1
          IF (dnumber2 > dnumber1) THEN
             dnumber2 = dnumber1
          END IF

          ! if (l_3mr) then
          !    dm1=dt*dmass/rain_params%c_x
          !    dm2=dt*dnumber2
          !    p1=rain_params%p1
          !    p2=rain_params%p2
          !    p3=rain_params%p3
          !    call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_aut)
          !    dm3=dm3/dt
          ! end if

!convert back to grid mean
      dmass=dmass*cf_liquid
      dnumber1=dnumber1*cf_liquid
      dnumber2=dnumber2*cf_liquid
      cloud_mass=cloud_mass*cf_liquid !for aerosol processing below

          procs(i_ql, i_praut%id)%column_data(k)=-dmass
          procs(i_qr, i_praut%id)%column_data(k)=dmass

          IF (cloud_params%l_2m) THEN
             procs(i_nl, i_praut%id)%column_data(k)=-dnumber1
          END IF
          IF (rain_params%l_2m) THEN
             procs(i_nr, i_praut%id)%column_data(k)=dnumber2
          END IF
          ! if (rain_params%l_3m) then
          !    procs(i_m3r, i_praut%id)%column_data(k)=dm3
          ! end if
          
          IF (l_separate_rain) THEN
             IF (l_aaut .AND. l_process) THEN
                ! Standard Single soluble mode, 2 activated species
                damass=dmass/cloud_mass*aerofields(k,i_am4)
                aerosol_procs(i_am4, i_aaut%id)%column_data(k)=-damass
                aerosol_procs(i_am5, i_aaut%id)%column_data(k)=damass
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE raut
END MODULE autoconversion
