! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Melting of ice, snow and graupel (melting).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Melting of cloud ice (instantaneous) and thermal-balance
!   melting of snow and graupel into rain.
!
! Paper reference:
!   Field et al. (2023) Appendix A.8.4, Eq. (A30)-(A34): Pimlt
!   = qi/Delta t; Psmlt/Pgmlt from thermal heat balance (including
!   the riming heating terms Psacw/Psacr, Pgacw/Pgacr); Nsmlt/Ngmlt
!   assume mean size is preserved during melting.
!
MODULE ice_melting
  USE variable_precision, ONLY: wp
  USE process_routines, ONLY: process_rate, process_name, i_imlt, i_smlt, i_gmlt, i_sacw, i_sacr, i_gacw, i_gacr, i_gshd, &
       i_dimlt, i_dsmlt, i_dgmlt
  USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
  USE passive_fields, ONLY: TdegC, qws0, rho
  USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params, rain_params, DR_melt, hydro_params, ZERO_REAL_WP
  USE mphys_switches, ONLY: i_qv, i_am4, i_am7, i_am8, i_am9, l_process, l_gamma_online
  USE mphys_constants, ONLY: Lv, Lf, Ka, Cwater, Cice, cp, Dv
  USE thresholds, ONLY: thresh_tidy
  USE m3_incs, ONLY: m3_inc_type2, m3_inc_type3, m3_inc_type4
  USE ventfac, ONLY: ventilation_1M_2M, ventilation_3M
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ICE_MELTING'

CONTAINS

  !> Subroutine to calculate rate of melting of ice species
  !>
  !> OPTIMISATION POSSIBILITIES: Shouldn't have to recalculate all 3m quantities
  !>                             If just rescaling mass conversion for dry mode
  ! Ice/snow/graupel melting, Field et al. (2023)
  ! Eq. (A30)-(A34).
  SUBROUTINE melting(ixy_inner, dt, nz, params, qfields, cffields, procs, l_sigevap, aeroice, dustact, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    USE mphys_switches,       ONLY: i_cfs, i_cfg

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
    REAL(wp), INTENT(IN) :: cffields(:,:)


    ! aerosol fields
    TYPE(aerosol_active), INTENT(IN) :: aeroice(:), dustact(:)

    ! optional aerosol fields to be processed
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

    LOGICAL, INTENT(IN) :: l_sigevap(:) ! logical to determine significant evaporation

    TYPE(process_name) :: iproc, iaproc ! processes selected depending on
    ! which species we're modifying
    TYPE(process_name) :: i_acw, i_acr ! accretion processes
    REAL(wp) :: qv
    REAL(wp) :: dmass, dnumber
!   real(wp) :: dm1, dm2, dm3, dm3_r, m2, m3
    REAL(wp) :: num, mass, m1
    REAL(wp) :: n0, lam, mu, V_x
    REAL(wp) :: acc_correction
    LOGICAL :: l_meltall ! do we melt everything?
    REAL(wp) :: dmac, dmad
    REAL(wp) :: cf

    INTEGER :: k
    
    CHARACTER(len=*), PARAMETER :: RoutineName='MELTING'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!intiliase
    dmac=0.0
    dmad=0.0
    dmass=0.0
    dnumber=0.0
    
    DO k = 1, nz
       IF (.NOT. l_sigevap(k)) THEN 
          l_meltall=.FALSE.
          mass=qfields(k, params%i_1m)
          IF (mass > thresh_tidy(params%i_1m) .AND. TdegC(k,ixy_inner) > 0.0) THEN
             IF (params%l_2m) num=qfields(k, params%i_2m)
             IF (params%id == ice_params%id) THEN ! instantaneous removal
                iaproc=i_dimlt
                iproc=i_imlt
                dmass=mass/dt
                
                IF (params%l_2m)THEN
                   dnumber=num/dt
                   procs(ice_params%i_2m, iproc%id)%column_data(k)=-dnumber
                   IF (rain_params%l_2m) THEN
                      procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                   END IF
                END IF
                
                procs(ice_params%i_1m, iproc%id)%column_data(k)=-dmass
                procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                
                ! 3-moment code retained for future implementation
                ! if (rain_params%l_3m) then
                !   m1=qfields(k, rain_params%i_1m)/rain_params%c_x
                !   m2=qfields(k, rain_params%i_2m)
                !   m3=qfields(k, rain_params%i_3m)
                !   dm1=dt*dmass/rain_params%c_x
                !   dm2=dt*dnumber
                
                !   if (m1 > 0.0) then
                !     call m3_inc_type2(m1, m2, m3, rain_params%p1,      &
                !          rain_params%p2, rain_params%p3, dm1, dm2, dm3_r)
                !   else
                !     call m3_inc_type3(rain_params%p1, rain_params%p2, rain_params%p3,      &
                !          dm1, dm2, dm3_r, rain_params%fix_mu)
                !   end if
                !   dm3_r=dm3_r/dt
                !   procs(rain_params%i_3m, iproc%id)%column_data(k) = dm3_r
                ! end if
             ELSE
                IF (params%id==snow_params%id) THEN
                   i_acw=i_sacw
                   i_acr=i_sacr
                   iproc=i_smlt
                   iaproc=i_dsmlt
                   cf=cffields(k,i_cfs)
                ELSE IF (params%id==graupel_params%id) THEN
                   i_acw=i_gacw
                   i_acr=i_gacr
                   iproc=i_gmlt
                   iaproc=i_dgmlt
                   cf=cffields(k,i_cfg)
                END IF
                acc_correction=0.0
                IF (i_acw%on)acc_correction=procs(params%i_1m, i_acw%id)%column_data(k)
                IF (i_acr%on)acc_correction=acc_correction + procs(params%i_1m, i_acr%id)%column_data(k)
                IF (params%id==graupel_params%id .AND. i_gshd%on) THEN
                   acc_correction=acc_correction + procs(params%i_1m, i_gshd%id)%column_data(k)
                END IF
                
                qv=qfields(k, i_qv)
                
                m1=mass/params%c_x
                IF (params%l_2m) num=qfields(k, params%i_2m)
                ! 3-moment code retained for future implementation
                ! if (params%l_3m) m3=qfields(k, params%i_3m)
                
                n0=dist_n0(k,params%id)
                mu=dist_mu(k,params%id)
                lam=dist_lambda(k,params%id)
                
                IF (l_gamma_online) THEN 
                   CALL ventilation_3M(ixy_inner, k, V_x, n0, lam, mu, params)
                ELSE 
                   CALL ventilation_1M_2M(ixy_inner, k, V_x, n0, lam, mu, params)
                END IF
                
                dmass=(1.0/(rho(k,ixy_inner)*Lf))*(Ka*TdegC(k,ixy_inner) + Lv*Dv*rho(k,ixy_inner)*(qv - qws0(k,ixy_inner))) * V_x&
                     + (Cwater*TdegC(k,ixy_inner)/Lf)*acc_correction
                dmass=dmass*cf ! grid mean


                dmass=max(dmass, ZERO_REAL_WP) ! ensure positive
                
                
                dmass=min(dmass, mass/dt) ! ensure we don't remove too much
                IF (dmass*dt > 0.95*mass) THEN ! we're pretty much removing everything
                   l_meltall=.TRUE.
                   dmass=mass/dt
                END IF
                
                !--------------------------------------------------
                ! Apply spontaneous rain breakup if drops are large
                !--------------------------------------------------
                !< RAIN BREAKUP TO BE ADDED
                
                procs(params%i_1m, iproc%id)%column_data(k)=-dmass
                procs(rain_params%i_1m, iproc%id)%column_data(k)=dmass
                
                IF (params%l_2m) THEN
                   dnumber=dmass*num/mass
                   procs(params%i_2m, iproc%id)%column_data(k)=-dnumber
                   procs(rain_params%i_2m, iproc%id)%column_data(k)=dnumber
                END IF
                ! 3-moment code retained for future implementation
                ! if (params%l_3m) then
                !   if (l_meltall) then
                !     dm3=-m3/dt
                !   else
                !     dm1=-dt*dmass/params%c_x
                !     dm2=-dt*dnumber
                !     m2=num
                !     call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
                !     dm3=dm3/dt
                !   end if
                !   procs(params%i_3m, iproc%id)%column_data(k)=dm3
                ! end if
                
                ! if (rain_params%l_3m) then
                !   if (params%l_3m) then
                !     call m3_inc_type4(dm3, rain_params%c_x, params%c_x, params%p3, dm3_r)
                !   else
                !     m1=qfields(k, rain_params%i_1m)/rain_params%c_x
                !     m2=qfields(k, rain_params%i_2m)
                !     m3=qfields(k, rain_params%i_3m)
                
                !     dm1=dt*dmass/rain_params%c_x
                !     dm2=dt*dnumber
                !     call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3_r, rain_params%fix_mu)
                !     dm3_r=dm3_r/dt
                !   end if
                !   procs(rain_params%i_3m, iproc%id)%column_data(k) = dm3_r
                ! end if
             END IF
             !----------------------
             ! Aerosol processing...
             !----------------------
             
             IF (l_process) THEN
                
                IF (params%id == ice_params%id) THEN
                   dmac=dnumber*aeroice(k)%nratio1*aeroice(k)%mact1_mean
                   dmac=min(dmac, aeroice(k)%mact1 /dt)
                   dmad=dnumber*dustact(k)%nratio1*dustact(k)%mact1_mean
                   dmad=min(dmad, dustact(k)%mact1 /dt)
                ELSE IF (params%id == snow_params%id) THEN
                   dmac=dnumber*aeroice(k)%nratio2*aeroice(k)%mact2_mean
                   dmac=min(dmac, aeroice(k)%mact2 /dt)
                   dmad=dnumber*dustact(k)%nratio2*dustact(k)%mact2_mean
                   dmad=min(dmad, dustact(k)%mact2 /dt)
                ELSE IF (params%id == graupel_params%id) THEN
                   dmac=dnumber*aeroice(k)%nratio3*aeroice(k)%mact3_mean
                   dmac=min(dmac, aeroice(k)%mact3/dt)
                   dmad=dnumber*dustact(k)%nratio3*dustact(k)%mact3_mean
                   dmad=min(dmad, dustact(k)%mact3/dt)
                END IF
                
                aerosol_procs(i_am8, iaproc%id)%column_data(k)=-dmac
                aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
                aerosol_procs(i_am9, iaproc%id)%column_data(k)=dmad
                aerosol_procs(i_am7, iaproc%id)%column_data(k)=-dmad
             END IF
          END IF
       END IF
    END DO
    
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE melting
END MODULE ice_melting
