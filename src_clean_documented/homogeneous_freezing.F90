! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Homogeneous freezing of rain and cloud droplets below -40C (ihom_rain, ihom_droplets).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Heterogeneous immersion freezing of rain drops to graupel
!   (Bigg, 1953) and homogeneous freezing of cloud droplets to
!   ice below the homogeneous-freezing temperature threshold.
!
! Paper reference:
!   Field et al. (2023) Appendix A.9.2, Eq. (A35)-(A37): rain
!   freezing probability follows Bigg (1953); Appendix A.9.3,
!   Eq. (A38): homogeneous droplet freezing Nhomc balances the
!   Squires supersaturation equation for a given updraught,
!   assuming 50 micron ice spheres are formed.
!
MODULE homogeneous
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho, pressure, w, exner
  USE mphys_switches, ONLY: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am6, i_an2, l_2mi, l_2ms, l_2mg &
       , i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_am9 , i_m3r, i_m3g, i_qr, i_qg, i_nr, i_ng, i_nl          &
       , isol, i_am4, i_am8, active_ice, l_process, l_prf_cfrac
  USE process_routines, ONLY: process_rate, i_homr, i_homc, i_dhomc, i_dhomr
  USE mphys_parameters, ONLY: rain_params, graupel_params, cloud_params, ice_params, T_hom_freeze
  USE mphys_constants, ONLY:  Ls, cp, Mw, g, Rd, Ru, Lv, Lf, Dv, Rv
  USE qsat_funs, ONLY: qsaturation, qisaturation
  USE thresholds, ONLY: thresh_small, thresh_tidy
  USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
  USE special, ONLY: GammaFunc, pi
  USE m3_incs, ONLY: m3_inc_type2, m3_inc_type4

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='HOMOGENEOUS'

CONTAINS

  !> Calculates immersion freezing of rain drops
  !> See Bigg 1953
  ! Bigg (1953) immersion freezing of rain to graupel,
  ! Field et al. (2023) Eq. (A35)-(A37).
  SUBROUTINE ihom_rain(ixy_inner, dt, nz, l_Tcold, qfields, l_sigevap, aeroact, dustliq,  &
                      procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:)
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! optional aerosol fields to be processed
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

    LOGICAL, INTENT(IN) :: l_sigevap(:) ! logical to determine significant evaporation

    REAL(wp) :: dmass, dnumber, dmac, coef, dmadl
!    real(wp) :: dm1, dm2, dm3, dm3_g, m1, m2, m3
    REAL(wp) :: n0, lam, mu

    REAL(wp) :: th
    REAL(wp) :: qr, nr

    REAL(wp) :: Tc

    REAL(wp), PARAMETER :: A_bigg = 0.66, B_bigg = 100.0

    LOGICAL :: l_condition, l_freezeall
    LOGICAL :: l_ziegler=.TRUE.

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='IHOM_RAIN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (l_Tcold(k)) THEN 
          th = qfields(k, i_th)
          Tc = th*exner(k,ixy_inner) - 273.15
          qr = qfields(k, i_qr)
          IF (rain_params%l_2m)THEN
             nr = qfields(k, i_nr)
          ELSE
             nr = 1000. ! PRAGMATIC SM HACK
          END IF
          
          l_condition=(Tc < -4.0 .AND. qr > thresh_small(i_qr) .AND. .NOT. l_sigevap(k))
          l_freezeall=.FALSE.
          IF (l_condition) THEN
             
             n0 = dist_n0(k,rain_params%id)
             mu = dist_mu(k,rain_params%id)
             lam = dist_lambda(k,rain_params%id)
             
             IF (l_ziegler) THEN
                dnumber = B_bigg*(exp(-A_bigg*Tc)-1.0)*rho(k,ixy_inner)*qr/rain_params%density
                dnumber = min(dnumber, nr/dt)
                dmass = (qr/nr)*dnumber
             ELSE
                coef = B_bigg*(pi/6.0)*(exp(-A_bigg*Tc)-1.0)/rho(k,ixy_inner)               &
                     * n0 /(lam*lam*lam)/GammaFunc(1.0 + mu)
                
                dmass = coef * rain_params%c_x * lam**(-rain_params%d_x)          &
                     * GammaFunc(4.0 + mu + rain_params%d_x)
                
                dnumber = coef                                                    &
                     * GammaFunc(4.0 + mu)
                
             END IF
             
             ! PRAGMATIC HACK - FIX ME
             ! If most of the drops are frozen, do all of them
             IF (dmass*dt >0.95*qr .OR. dnumber*dt > 0.95*nr) THEN
                dmass=qr/dt
                dnumber=nr/dt
                l_freezeall=.TRUE.
             END IF

             procs(i_qr, i_homr%id)%column_data(k) = -dmass
             procs(i_qg, i_homr%id)%column_data(k) = dmass
             
             IF (rain_params%l_2m) THEN
                procs(i_nr, i_homr%id)%column_data(k) = -dnumber
             END IF
             IF (graupel_params%l_2m) THEN
                procs(i_ng, i_homr%id)%column_data(k) = dnumber
             END IF
             ! 3-moment code is commented for future implementation
             ! if (rain_params%l_3m) then
             !   if (l_freezeall) then
             !     dm3=-qfields(k,i_m3r)/dt
             !   else
             !     m1=qr/rain_params%c_x
             !     m2=qfields(k,i_nr)
             !     m3=qfields(k,i_m3r)
             
             !     dm1=-dt*dmass/rain_params%c_x
             !     dm2=-dt*dnumber
             
             !     call m3_inc_type2(m1, m2, m3, rain_params%p1, rain_params%p2, rain_params%p3, dm1, dm2, dm3)
             !     dm3=dm3/dt
             !   end if
             !   procs(i_m3r, i_homr%id)%column_data(k) = dm3
             ! end if
             
             ! if (graupel_params%l_3m) then
             !   if (rain_params%l_3m) then
             !     call m3_inc_type4(dm3, graupel_params%c_x, rain_params%c_x, rain_params%p3, dm3_g)
             !   else
             !     m1=qfields(k,i_qg)/graupel_params%c_x
             !     m2=qfields(k,i_ng)
             !     m3=qfields(k,i_m3g)
             
             !     dm1=-dt*dmass/graupel_params%c_x
             !     dm2=-dt*dnumber
             !     call m3_inc_type2(m1, m2, m3, graupel_params%p1, graupel_params%p2, graupel_params%p3, dm1, dm2, dm3_g)
             !     dm3_g=dm3_g/dt
             !   end if
             !   procs(i_m3g, i_homr%id)%column_data(k) = dm3_g
             ! end if

             IF (l_process) THEN
                dmac = dnumber*aeroact(k)%mact2_mean

                aerosol_procs(i_am8, i_dhomr%id)%column_data(k) = dmac
                aerosol_procs(i_am4, i_dhomr%id)%column_data(k) = -dmac
                
                ! Dust already in the liquid phase
                dmadl = dnumber*dustliq(k)%mact2_mean*dustliq(k)%nratio2
                IF (dmadl /=0.0) THEN
                   aerosol_procs(i_am9, i_dhomr%id)%column_data(k) = -dmadl
                   aerosol_procs(i_am7, i_dhomr%id)%column_data(k) = dmadl
                END IF
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ihom_rain

  !> Calculates homogeneous freezing of cloud drops
  !> See Wisener 1972
  ! Homogeneous freezing of cloud droplets below the -38C
  ! threshold, Field et al. (2023) Eq. (A38).
  SUBROUTINE ihom_droplets(ixy_inner, dt, nz, l_Tcold, qfields, aeroact, dustliq, procs, aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:) 
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    ! aerosol fields
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustliq(:)

    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

    REAL(wp) :: dmass, dnumber, dmac, dmadl
!    real(wp) :: m1, m2, m3

    REAL(wp) :: th
    REAL(wp) :: qv, ql

    REAL(wp) :: Tc

    REAL(wp) :: Tk, cap, rhoi, Ei, Ew, bm, Ai, B0, Bis, aw, dnumberi, &
                ka, min_homog_ni, dniraw, d0_homog
    LOGICAL :: l_use_critical_w = .TRUE.  !ni controlled by w and environmental conditions
    LOGICAL :: l_use_ni_limit = .FALSE.   !ni limited to max per timestep
    
    
    LOGICAL :: l_condition

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='IHOM_DROPLETS'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (l_Tcold(k)) THEN
          th = qfields(k, i_th)
          Tc = th*exner(k,ixy_inner) - 273.15
          ql = qfields(k, i_ql) 
          
          l_condition=(Tc < T_hom_freeze .AND. ql > thresh_tidy(i_ql))
          IF (l_condition) THEN
             dmass=min(ql, cp*(T_hom_freeze - Tc)/Lf)/dt    
             IF (cloud_params%l_2m) THEN    
                dnumber=dmass*(qfields(k, i_nl))/ql   
                dnumberi=dnumber ! this will be overwritten if l_use_critical_W is used   
                
                !limit production to 1/cc in timestep - do this properly with KL but need to adapt for droplets
                IF (l_use_ni_limit) THEN
                   dnumberi=min(dnumber, 1e2)
                   dmass=qfields(k, i_ql)*dnumberi/qfields(k, i_nl)
                END IF
                
                IF (l_use_critical_w) THEN
                   !! alternative limiter based on w - this is explicit w. Will overwrite dnumber
                   !! Method is based on derivation for sink of vapour to ice defined on Field et al, 2014, 
                   !! Mixed-phase clouds in a turbulent environment. Part 2: Analytic treatment, 
                   !! https://rmets.onlinelibrary.wiley.com/doi/full/10.1002/qj.2175
                   qv = qfields(k, i_qv)
                   Tk = th*exner(k,ixy_inner)
                   ka=((5.69+0.017*(Tk-273.15))*1e-5) * 418.6 
                   
                   cap=1.0 !assume spheres and radius used
                   rhoi=200.0
                   min_homog_ni=1e2 !kg-1
                   d0_homog=50e-6 !m
                   
                   Ei=qisaturation(Tk,pressure(k,ixy_inner)/100.)*pressure(k,ixy_inner)/ &
                        (qisaturation(Tk,pressure(k,ixy_inner)/100.)+0.62198) !vap press over ice [Pa]
                   Ew=qsaturation(Tk,pressure(k,ixy_inner)/100.)*pressure(k,ixy_inner)/ &
                        (qsaturation(Tk,pressure(k,ixy_inner)/100.)+0.62198) !vap press over liq [Pa]
                   
                   bm=1.0/(qv)+Lv*Lf/(cp*Rv*Tk**2)
                   Ai=1.0/(rhoi*Lf**2/(ka*Rv*Tk**2)+rhoi*Rv*Tk/(Ei*Dv))
                   B0=4.0*pi*cap*rhoi*Ai/rho(k,ixy_inner)
                   Bis=bm*B0*(Ew/Ei-1.0)
                   aw=g/(Rd*Tk)*(Lv*Rd/(cp*Rv*Tk)-1.0)
                   
                   dnumberi=max(w(k,ixy_inner),0.0)*(aw/Bis)/d0_homog  / dt  !convert to a rate
                   dniraw=dnumberi
                   !make max just below limit
                   dnumberi=min(dnumber*0.90, max(min_homog_ni, dnumberi))
                END IF
             END IF


             procs(i_ql, i_homc%id)%column_data(k)=-dmass  
             procs(i_qi, i_homc%id)%column_data(k)=dmass   
             
             IF (cloud_params%l_2m) THEN
                procs(i_nl, i_homc%id)%column_data(k)=-dnumber   
                IF (ice_params%l_2m) THEN
                   procs(i_ni, i_homc%id)%column_data(k)=dnumberi   
                   ! if l_use_critical_w or l_limit_ni is true then 
                   ! dnumber not necessarily equal to dnumberi 
                END IF
             END IF
             
             
!!!NB if we keep l_use_critical_w approach need to deal with dnumber=/=dnumberi for processing!!

             IF (l_process) THEN
                dmac=dnumber*aeroact(k)%mact1_mean
                
                aerosol_procs(i_am8, i_dhomc%id)%column_data(k)=dmac
                aerosol_procs(i_am4, i_dhomc%id)%column_data(k)=-dmac
                
                ! Dust already in the liquid phase
                dmadl=dnumber*dustliq(k)%mact1_mean*dustliq(k)%nratio1
                IF (dmadl /=0.0) THEN
                   aerosol_procs(i_am9, i_dhomc%id)%column_data(k)=-dmadl
                   aerosol_procs(i_am7, i_dhomc%id)%column_data(k)=dmadl
                END IF
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ihom_droplets
END MODULE homogeneous
