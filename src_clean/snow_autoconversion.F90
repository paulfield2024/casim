! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Ice-to-snow autoconversion (saut), Bergeron process.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE snow_autoconversion
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho, exner, TdegK, pressure
  USE mphys_switches, ONLY:  &
       i_qi, i_qs, i_ni, i_ns, &
       i_qv, i_th, &
       l_harrington
! use mphys_switches, only:  i_m3s
  USE mphys_constants, ONLY: fixed_ice_number,   &
       Lv,ka, Dv, Rv
! use mphys_parameters, only: mu_saut
  USE mphys_parameters, ONLY: snow_params, ice_params,   &
       DImax, tau_saut, DI2S
  USE process_routines, ONLY: process_rate, i_saut
  USE qsat_funs, ONLY: qisaturation
  USE thresholds, ONLY: thresh_small
  USE special, ONLY: pi
! use m3_incs, only: m3_inc_type3

  USE distributions, ONLY: dist_lambda, dist_mU

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SNOW_AUTOCONVERSION'

  PUBLIC saut
CONTAINS

  SUBROUTINE saut(ixy_inner, dt, nz, l_Tcold, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SAUT'

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:) 
    REAL(wp), INTENT(IN) :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    REAL(wp) :: dmass, dnumber
!   real(wp) :: dm1,dm2,dm3
    REAL(wp) :: ice_lam, ice_mu
    REAL(wp) :: ice_mass
    REAL(wp) :: ice_number
    REAL(wp) :: th
    REAL(wp) :: qv
    REAL(wp) :: lami_min , AB, qis
    INTEGER :: k

    LOGICAL :: l_condition ! logical condition to switch on process

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
      IF (l_Tcold(k)) THEN
          
        ice_mass=qfields(k, i_qi)
        IF (ice_params%l_2m) THEN
          ice_number=qfields(k, i_ni)
        ELSE
          ice_number=fixed_ice_number
        END IF
          
        IF (ice_mass >= thresh_small(i_qi) .AND. ice_number >= thresh_small(i_ni)) THEN
          ice_mu=dist_mu(k,ice_params%id)
          ice_lam=dist_lambda(k,ice_params%id)
             
          lami_min=(1.0 + ice_mu)/DImax
             
          IF (l_harrington) THEN
            qv=qfields(k, i_qv)
            th=qfields(k, i_th)
            qis=qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
            l_condition=qv > qis
          ELSE
            l_condition=ice_lam < lami_min ! LEM autconversion
          END IF
          
          IF (l_condition) THEN
             
            IF (l_harrington) THEN
              !< AB This is used elsewhere, so we should do it more efficiently.
              AB=1.0/(Lv*Lv/(Rv*ka*TdegK(k,ixy_inner)*TdegK(k,ixy_inner))*   &
                 rho(k,ixy_inner)+1.0/(Dv*qis))
              dnumber=4.0/DImax/ice_params%density*(qv-qis)*                 &
                      rho(k,ixy_inner)*ice_number*exp(-ice_lam*DImax)*Dv/AB
              dnumber=min(dnumber,0.9*ice_number/dt)
              dmass=pi/6.0*ice_params%density*DImax**3*dnumber
              dmass=min(dmass, 0.7*ice_mass/dt)
            ELSE ! LEM version
              dmass=((lami_min/ice_lam)**ice_params%d_x - 1.0)*ice_mass/tau_saut
              IF (ice_mass < dmass*dt) THEN
                dmass=.5*ice_mass/dt
              END IF
              dmass=min(dmass, .5*ice_mass/dt)
              dnumber=dmass/(ice_params%c_x*DI2S**ice_params%d_x)
            END IF
             
            IF (ice_number < dnumber*dt) THEN
              dnumber=dmass/ice_mass * ice_number
            END IF
             

            IF (dmass*dt > thresh_small(i_qs)) THEN
              ! if (snow_params%l_3m) then
              !   dm1=dt*dmass/snow_params%c_x
              !   dm2=dt*dnumber
              !   call m3_inc_type3(p1, p2, p3, dm1, dm2, dm3, mu_saut)
              !   dm3=dm3/dt
              ! end if
               
              procs(i_qi, i_saut%id)%column_data(k)=-dmass
              procs(i_qs, i_saut%id)%column_data(k)=dmass
                
              IF (ice_params%l_2m) THEN
                 procs(i_ni, i_saut%id)%column_data(k)=-dnumber
              END IF
              IF (snow_params%l_2m) THEN
                procs(i_ns, i_saut%id)%column_data(k)=dnumber
              END IF
              ! if (snow_params%l_3m) then
              !   procs(i_m3s, i_saut%id)%column_data(k)=dm3
              ! end if
            END IF
            !==============================
            ! No aerosol processing needed
            !==============================
          END IF
        END IF
      END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE saut
END MODULE snow_autoconversion
