! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Rain evaporation and CCN re-activation bookkeeping (revp).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE evaporation
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho, qws, TdegK
! use mphys_switches, only: i_m3r, l_3mr
  USE mphys_switches, ONLY: i_qr, i_nr, i_qv, i_ql, l_2mr, i_am2, &
       i_an2, i_am3, i_an3, i_am4, i_am5, l_process, aero_index, &
       l_separate_rain, i_am6, i_an6, i_am9, l_warm, i_an11, i_an12, l_passivenumbers, &
       l_passivenumbers_ice, l_inhom_revp, l_gamma_online, &
       l_bypass_which_mode, iopt_which_mode
  USE mphys_constants, ONLY: Lv,ka, Dv, Rv
  USE mphys_parameters, ONLY: c_r, rain_params
  USE process_routines, ONLY: process_rate, i_prevp, i_arevp
  USE thresholds, ONLY: qr_small, ss_small, qr_tidy, ql_tidy
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
  USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
  USE ventfac, ONLY: ventilation_1M_2M, ventilation_3M
  USE which_mode_to_use, ONLY : which_mode
  USE mphys_die, ONLY: throw_mphys_error, incorrect_opt, std_msg
  
  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='EVAPORATION'

  PRIVATE

  PUBLIC revp
CONTAINS

  SUBROUTINE revp(ixy_inner, dt, nz, qfields, cffields, aerophys, aerochem, aeroact, dustliq, procs, aerosol_procs, l_sigevap)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfr

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: qfields(:,:)
    INTEGER, INTENT(IN) :: nz
    TYPE(aerosol_phys), INTENT(IN) :: aerophys(:)
    TYPE(aerosol_chem), INTENT(IN) :: aerochem(:)
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
    TYPE(process_rate), INTENT(INOUT) :: aerosol_procs(:,:)
    LOGICAL, INTENT(OUT) :: l_sigevap(:) ! Determines if there is significant evaporation
    REAL(wp), INTENT(IN) :: cffields(:,:)

    ! Local variables
    REAL(wp) :: dmass, dnumber, dnumber_a, dnumber_d
    REAL(wp) :: m1, m2, dm1
!   real(wp) :: m3, dm3

    REAL(wp) :: n0, lam, mu
    REAL(wp) :: V_r, AB

    REAL(wp) :: rain_mass
    REAL(wp) :: rain_number
!   real(wp) :: rain_m3
    REAL(wp) :: qv
    REAL(wp) :: cf

    LOGICAL :: l_rain_test ! conditional test on rain

    REAL(wp) :: dmac, dmac1, dmac2, dnac1, dnac2, dmacd
    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='REVP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    DO k = 1, nz
       
       l_sigevap(k)=.FALSE.
       
       qv=qfields(k, i_qv)
       rain_mass=qfields(k, i_qr)
       IF (l_2mr)rain_number=qfields(k, i_nr)
       ! if (l_3mr)rain_m3=qfields(k, i_m3r)
       
       IF (qv/qws(k,ixy_inner) < 1.0-ss_small .AND. qfields(k, i_ql) < ql_tidy .AND. rain_mass > qr_tidy) THEN
          
          m1=rain_mass/c_r
          IF (l_2mr) m2=rain_number
          ! if (l_3mr) m3=rain_m3

          l_rain_test=.FALSE.
          IF (l_2mr) l_rain_test=rain_number>0
          IF (rain_mass > qr_small .AND. (.NOT. l_2mr .OR. l_rain_test)) THEN
             n0=dist_n0(k,rain_params%id)
             mu=dist_mu(k,rain_params%id)
             lam=dist_lambda(k,rain_params%id)
             cf=cffields(k,i_cfr)

             IF (l_gamma_online) THEN 
                CALL ventilation_3M(ixy_inner, k, V_r, n0, lam, mu, rain_params)
             ELSE
                CALL ventilation_1M_2M(ixy_inner, k, V_r, n0, lam, mu, rain_params)
             END IF

             AB=1.0/(Lv**2/(Rv*ka)*rho(k,ixy_inner)*TdegK(k,ixy_inner)**(-2)+1.0/(Dv*qws(k,ixy_inner)))
             dmass=(1.0-qv/qws(k,ixy_inner))*V_r*AB  *cf !grid mean
             
             dmass=MIN(dmass,rain_mass/dt)
           
          ELSE
             dmass=rain_mass/dt
          END IF

          dm1=dmass/c_r
          IF (l_2mr) THEN
             dnumber=0.0
             IF (l_inhom_revp) dnumber=dm1*m2/m1
          END IF
          ! if (l_3mr) dm3=dm1*m3/m1

          IF (l_2mr) THEN
             IF (dnumber*dt > rain_number .OR. dmass*dt >= rain_mass-qr_tidy) THEN
                dmass=rain_mass/dt
                dnumber=rain_number/dt
               ! if (l_3mr) dm3=m3/dt
             END IF
          END IF

          procs(i_qr, i_prevp%id)%column_data(k)=-dmass
          procs(i_qv, i_prevp%id)%column_data(k)=dmass

          IF (dmass*dt/rain_mass > .8) l_sigevap(k)=.TRUE.

          IF (l_2mr) THEN
             procs(i_nr, i_prevp%id)%column_data(k)=-dnumber
          END IF
          ! if (l_3mr) then
          !    procs(i_m3r, i_prevp%id)%column_data(k)=-dm3
          ! end if
     
      !============================
      ! aerosol processing
      !============================
          IF (l_process .AND. abs(dnumber) >0) THEN

             dmac=dnumber*aeroact(k)%nratio2*aeroact(k)%mact2_mean
             dmac=min(dmac,aeroact(k)%mact2/dt)
             IF (l_separate_rain) THEN
                aerosol_procs(i_am5, i_arevp%id)%column_data(k)=-dmac
             ELSE
                aerosol_procs(i_am4, i_arevp%id)%column_data(k)=-dmac
             END IF
             IF (l_passivenumbers) THEN
                dnumber_a=-dnumber*aeroact(k)%nratio2
                aerosol_procs(i_an11, i_arevp%id)%column_data(k)=dnumber_a
             END IF

             IF (l_passivenumbers_ice) THEN
                dnumber_d=-dnumber*dustliq(k)%nratio2
                aerosol_procs(i_an12, i_arevp%id)%column_data(k)=dnumber_d
             END IF

             ! Return aerosol
             IF (aero_index%i_accum >0 .AND. aero_index%i_coarse >0) THEN
                ! Coarse and accumulation mode being used. Which one to return to?
                IF (l_bypass_which_mode) THEN
                   !Don't use which_mode - just transfer all to either
                   !accum or coarse mode
                  IF (iopt_which_mode==1) THEN !All to accum
                       dmac1=dmac
                       dmac2=0.0
                       dnac1=dnumber_a
                       dnac2=0.0
                  ELSE IF (iopt_which_mode==2) THEN !All to coarse
                       dmac1=0.0
                       dmac2=dmac
                       dnac1=0.0
                       dnac2=dnumber_a
                  ELSE
                       WRITE(std_msg, '(A)') "incorrect iopt_which_mode option selected"
                       CALL throw_mphys_error(incorrect_opt,ModuleName, std_msg)
                  END IF
                ELSE
                  CALL which_mode(dmac, dnumber*aeroact(k)%nratio2, aerophys(k)%rd(aero_index%i_accum), &
                     aerophys(k)%rd(aero_index%i_coarse), aerochem(k)%density(aero_index%i_accum),    &
                     aerophys(k)%sigma(aero_index%i_accum),                                           &
                     dmac1, dmac2, dnac1, dnac2)
                END IF !end of bypass whichmode
                
                aerosol_procs(i_am2, i_arevp%id)%column_data(k)=dmac1
                aerosol_procs(i_an2, i_arevp%id)%column_data(k)=dnac1
                aerosol_procs(i_am3, i_arevp%id)%column_data(k)=dmac2
                aerosol_procs(i_an3, i_arevp%id)%column_data(k)=dnac2
             ELSE
                IF (aero_index%i_accum >0) THEN
                   aerosol_procs(i_am2, i_arevp%id)%column_data(k) = dmac
                   aerosol_procs(i_an2, i_arevp%id)%column_data(k) = dnumber *             &
                        aeroact(k)%nratio2
                END IF
                IF (aero_index%i_coarse >0) THEN
                   aerosol_procs(i_am3, i_arevp%id)%column_data(k) = dmac
                   aerosol_procs(i_an3, i_arevp%id)%column_data(k) = dnumber *             &
                        aeroact(k)%nratio2
                END IF
             END IF

             dmacd=dnumber*dustliq(k)%nratio2*dustliq(k)%mact2_mean
             IF (.NOT. l_warm .AND. dmacd /=0.0) THEN
                aerosol_procs(i_am9, i_arevp%id)%column_data(k)=-dmacd
                aerosol_procs(i_am6, i_arevp%id)%column_data(k)=dmacd
                aerosol_procs(i_an6, i_arevp%id)%column_data(k)=dnumber*dustliq(k)%nratio2
             END IF

          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE revp
END MODULE evaporation
