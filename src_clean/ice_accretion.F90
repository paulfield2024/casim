! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Binary accretion process rates between hydrometeor species (iacc), e.g. snow/graupel accreting cloud or rain.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE ice_accretion
USE variable_precision, ONLY: wp, iwp
USE mphys_die, ONLY: throw_mphys_error, incorrect_opt, std_msg
USE passive_fields, ONLY: rho, TdegC, TdegK
USE mphys_switches, ONLY: l_process,   &
     cloud_params, rain_params,  ice_params, snow_params, &
     i_am4, i_am7, i_am8, i_am9, l_prf_cfrac, &
     l_g, l_raci_g, l_prf_cfrac, &
     i_cfi, i_cfs, i_cfg, i_cfl, i_cfr, mpof, l_srg
USE type_process, ONLY: process_name
USE process_routines, ONLY: process_rate,   &
     i_iacw, i_sacw, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, &
     i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr, i_draci
USE mphys_parameters, ONLY: hydro_params
USE thresholds, ONLY: thresh_small, cfliq_small, thresh_sig 
USE aerosol_routines, ONLY: aerosol_active

USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
USE sweepout_rate, ONLY: sweepout, binary_collection
!  use sweepout_rate, only: binary_collection_1M2M, sweepout_1M2M
  USE casim_stph, ONLY: l_rp2_casim, mpof_casim_rp

IMPLICIT NONE
PRIVATE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ICE_ACCRETION'

PUBLIC iacc
CONTAINS

SUBROUTINE iacc(ixy_inner, dt, nz, l_Tcold, params_X, params_Y, params_Z, qfields, cffields, procs, &
                l_sigevap, aeroact, dustliq, aerosol_procs, params_snow)
    !
    !< CODE TIDYING: Move efficiencies into parameters

USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

! Subroutine arguments
INTEGER, INTENT(IN) :: ixy_inner
REAL(wp), INTENT(IN) :: dt
INTEGER, INTENT(IN) :: nz
LOGICAL, INTENT(IN) :: l_Tcold(:) 
TYPE(hydro_params), INTENT(IN) :: params_X !< parameters for species which does the collecting
TYPE(hydro_params), INTENT(IN) :: params_Y !< parameters for species which is collected
TYPE(hydro_params), INTENT(IN) :: params_Z !< parameters for species to which resulting amalgamation is sent
TYPE(hydro_params), INTENT(IN), OPTIONAL :: params_snow 
                                              !< snow parameters to which resulting 
                                              !< amalgamation may be sent when rain collecting snow, where 
                                              !< rain mass determines whether graupel (params_Z) or snow,
                                              !< and snow collecting rain T determines whether graupel 
                                              !< (params_Z) or snow.
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
REAL(wp), INTENT(IN) :: cffields(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
TYPE(aerosol_active), INTENT(IN) :: aeroact(:)
TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)

LOGICAL, INTENT(IN) :: l_sigevap(:) ! logical to determine significant evaporation

! Local variables
TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
! which species we're depositing on.

REAL(wp) :: dmass_Y, dnumber_Y, dmac, dmad
REAL(wp) :: dmass_X, dmass_Z, dnumber_X, dnumber_Z

REAL(wp) :: number_X, mass_X
REAL(wp) :: mass_Y, number_Y

REAL(wp) :: n0_X, lam_X, mu_X
REAL(wp) :: n0_Y, lam_Y, mu_Y

REAL(wp) :: cf_X, cf_Y, overlap_cf

REAL(wp) :: Eff ! collection efficiencies need to re-evaluate these and put them in properly to mphys_parameters

LOGICAL :: l_condition, l_alternate_Z
LOGICAL :: l_slow ! fall speed is slow compared to interactive species
LOGICAL :: l_aero ! If true then this process will modify aerosol

INTEGER :: k

! local variables for the indexing of the amalgation hydrometeor of X and Y. These 
! are required because when rain collects ice or snow collects rain, the resulting 
! hydrometeor depends on rain mass and T, respectively
INTEGER :: Zid ! id of the hydrometeor
INTEGER :: Zi_1m ! index for mass of hydrometeor
LOGICAL :: l_Zi_2m ! logical to determine if double moment resulting hydrometeor
INTEGER :: Zi_2m ! index for number of hydrometeor

CHARACTER(len=*), PARAMETER :: RoutineName='IACC'

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise increments
dmass_X = 0.0_wp
dnumber_X = 0.0_wp

! Apply RP scheme
IF ( l_rp2_casim ) THEN
    mpof = mpof_casim_rp
END IF

!setup cloud fractions
DO k = 1, nz
  IF (l_Tcold(k)) THEN   
    IF (l_prf_cfrac) THEN
      SELECT CASE (params_X%id)
      CASE(1) 
        cf_X=max(cffields(k,i_cfl), cfliq_small)
      CASE(2)
        cf_X=max(cffields(k,i_cfr), cfliq_small)
      CASE(3) 
        cf_X=max(cffields(k,i_cfi), cfliq_small)
      CASE(4)
        cf_X=max(cffields(k,i_cfs), cfliq_small)
      CASE(5) 
        cf_X=max(cffields(k,i_cfg), cfliq_small)
      END SELECT
      SELECT CASE (params_Y%id)
      CASE(1) 
        cf_Y=max(cffields(k,i_cfl), cfliq_small)
      CASE(2)
        cf_Y=max(cffields(k,i_cfr), cfliq_small)
      CASE(3) 
        cf_Y=max(cffields(k,i_cfi), cfliq_small)
      CASE(4)
        cf_Y=max(cffields(k,i_cfs), cfliq_small)
      CASE(5) 
        cf_Y=max(cffields(k,i_cfg), cfliq_small)
      END SELECT
    ELSE
      cf_X=1.0
      cf_Y=1.0
    END IF

    mass_Y=qfields(k, params_Y%i_1m) / cf_Y  !only do cloud species
    mass_X=qfields(k, params_X%i_1m) / cf_X  !only do cloud species

    !prf
    
    l_condition=mass_X * cf_X > thresh_small(params_X%i_1m) .AND. & 
                mass_Y * cf_Y > thresh_small(params_Y%i_1m)
    
    ! Check for rain in collector or collected and then check for significant evap
    IF ((params_X%id == rain_params%id .OR. params_Y%id == rain_params%id) &
         .AND. l_sigevap(k))  l_condition =.FALSE.

    
    IF (l_condition) THEN
      ! initialize variables which may not have been set
      number_X=0.0
      dnumber_Y=0.0

      Zid = params_Z%id 
      Zi_1m = params_Z%i_1m 
      IF (params_Z%l_2m) THEN 
        l_Zi_2m = params_Z%l_2m
        Zi_2m = params_Z%i_2m
      END IF

      ! determine resulting hydrometeor type from rain collecting ice
      !
      IF (params_X%id == rain_params%id .AND.  params_Y%id == ice_params%id) THEN
        IF (mass_X > thresh_sig(params_X%i_1m) .AND. l_g .AND. l_raci_g) THEN
          ! resulting hydrometeor is graupel
          Zid = params_Z%id
          Zi_1m = params_Z%i_1m
          IF (params_Z%l_2m) THEN
            l_Zi_2m = params_Z%l_2m
            Zi_2m = params_Z%i_2m
          END IF
        ELSE ! resulting hydrometeor is snow
          IF (PRESENT(params_snow)) THEN
            Zid = params_snow%id
            Zi_1m = params_snow%i_1m
            IF (params_snow%l_2m) THEN
              l_Zi_2m = params_snow%l_2m
              Zi_2m = params_snow%i_2m
            END IF
          ELSE
            WRITE(std_msg, '(A)') 'Resultant of iacc, praci is snow but snow_params '//&
                         'opt arg not passed'
            CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                          std_msg)
          END IF
        END IF
      END IF
      ! determine resulting hydrometeor type from snow collecting rain
      !
      IF (params_X%id == snow_params%id .AND.  params_Y%id == rain_params%id) THEN
        IF (TdegK(k,ixy_inner) < 268.15 .AND. l_g .AND. l_srg) THEN
          ! resulting hydrometeor is graupel
          Zid = params_Z%id 
          Zi_1m = params_Z%i_1m 
          IF (params_Z%l_2m) THEN 
            l_Zi_2m = params_Z%l_2m
            Zi_2m = params_Z%i_2m
          END IF
        ELSE ! resulting hydrometeor is snow 
          IF (PRESENT(params_snow)) THEN 
            Zid = params_snow%id 
            Zi_1m = params_snow%i_1m 
            IF (params_snow%l_2m) THEN 
               l_Zi_2m = params_snow%l_2m
               Zi_2m = params_snow%i_2m
            END IF
          ELSE
            WRITE(std_msg, '(A)') 'Resultant of iacc, pracr is snow but snow_params '//&
                                  'opt arg not passed'
            CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, std_msg)
          END IF
        END IF
      END IF

      l_alternate_Z=params_X%id /=  Zid ! params_Z%id
      l_aero=.FALSE.

      SELECT CASE (params_Y%id)
      CASE(1_iwp) ! cloud liquid is collected
        l_slow=.TRUE. ! Y catagory falls slowly compared to X
        SELECT CASE (params_X%id)
        CASE (3_iwp) !ice collects
          iproc=i_iacw
          Eff=0.0 !1.0  !PRF turn off ice collecting cloud drops - too small
          iaproc=i_diacw
          l_aero=.TRUE.
        CASE (4_iwp) !snow collects
          iproc=i_sacw
          Eff=0.5 !make assumption that collecting area is approx half of circle - similar to operational.
          iaproc=i_dsacw
          l_aero=.TRUE.
        CASE (5_iwp) !graupel collects
          iproc=i_gacw
          Eff=1.0
          iaproc=i_dgacw
          l_aero=.TRUE.
        END SELECT
      CASE(2_iwp) ! rain is collected
        l_slow=.FALSE. ! Y catagory falls slowly compared to X
        SELECT CASE (params_X%id)
        CASE (4_iwp) !snow collects
          iproc=i_sacr
          Eff=1.0
          iaproc=i_dsacr
          l_aero=.TRUE.
        CASE (5_iwp) !graupel collects
          iproc=i_gacr
          Eff=1.0
          iaproc=i_dgacr
          l_aero=.TRUE.
        END SELECT
      CASE(3_iwp)! cloud ice is collected
        l_slow=.TRUE. ! Y catagory falls slowly compared to X
        SELECT CASE (params_X%id)
        CASE (2_iwp) !rain collects
          iproc=i_raci
          Eff=1.0          ! Rain collecting ice efficiency to make graupel =~0.5? Lew and Pruppacher 1983
          iaproc=i_draci
          l_aero=.TRUE.
        CASE (4_iwp) !snow collects
          iproc=i_saci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        CASE (5_iwp) !graupel collects
          iproc=i_gaci
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        END SELECT
      CASE(4_iwp) ! snow is collected
        l_slow=.FALSE. ! Y catagory falls slowly compared to X
        SELECT CASE (params_X%id)
        CASE (5_iwp) !graupel collects
          iproc=i_gacs
          !Eff=1.0
          Eff = min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
        END SELECT
      END SELECT

      IF (params_X%l_2m) number_X=qfields(k, params_X%i_2m)  / cf_X

      n0_X=dist_n0(k,params_X%id)
      mu_X=dist_mu(k,params_X%id)
      lam_X=dist_lambda(k,params_X%id)

      IF (l_slow) THEN  ! collected species is approximated to have zero fallspeed
        !if (l_gamma_online) then
        dmass_Y=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k,ixy_inner))*mass_Y
       !else
       !   dmass_Y=-Eff*sweepout_1M2M(n0_X, lam_X, params_X, rho(k))*mass_Y
       !endif
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        IF (params_Y%l_2m) THEN
          number_Y=qfields(k, params_Y%i_2m) / cf_Y
          dnumber_Y=dmass_Y*number_Y/mass_Y
        END IF

        IF (l_alternate_Z) THEN ! We move resulting collision to another species.
          IF (params_Y%l_2m) THEN
            number_Y=qfields(k, params_Y%i_2m)  / cf_Y
          ELSE
            ! we need an else in here
          END IF

          !if (l_gamma_online) then
          dmass_X=-Eff*sweepout(n0_X, lam_X, mu_X, params_X, rho(k,ixy_inner), mass_weight=.TRUE.) * number_Y
          !else
          !   dmass_X=-Eff*sweepout_1M2M(n0_X, lam_X, params_X, rho(k), mass_weight=.true.) * number_Y
          !endif

          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y !dmass_Y*number_X/mass_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        END IF
      ELSE  ! both species have significant fall velocity
        IF (params_Y%l_2m) THEN
          number_Y=qfields(k, params_Y%i_2m) / cf_Y
        END IF

        n0_Y=dist_n0(k,params_Y%id)
        mu_Y=dist_mu(k,params_Y%id)
        lam_Y=dist_lambda(k,params_Y%id)

        !if (l_gamma_online) then
        dmass_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, mu_Y,&
                   params_X, params_Y, rho(k,ixy_inner), mass_weight=.TRUE.)
        !else
        !   dmass_Y=-Eff*binary_collection_1M2M(n0_X, lam_X, n0_Y, lam_Y, params_X, params_Y, rho(k), mass_weight=.true.)
        !endif
        dmass_Y=max(dmass_Y, -mass_Y/dt)

        IF (params_Y%l_2m) THEN
           !if (l_gamma_online) then
          dnumber_Y=-Eff*binary_collection(n0_X, lam_X, mu_X, n0_Y, lam_Y, &
                                mu_Y, params_X, params_Y, rho(k,ixy_inner))
           !else
           !   dnumber_Y=-Eff*binary_collection_1M2M(n0_X, lam_X, n0_Y, lam_Y, params_X, params_Y, rho(k))
           !endif
        END IF

        IF (l_alternate_Z) THEN ! We move resulting collision to another species.
           !if (l_gamma_online) then 
          dmass_X=-Eff*binary_collection(n0_Y, lam_Y, mu_Y, n0_X, lam_X,   &
            mu_X, params_Y, params_X, rho(k,ixy_inner), mass_weight=.TRUE.)
         !else
         !   dmass_X=-Eff*binary_collection_1M2M(n0_Y, lam_Y, n0_X, lam_X, params_Y, params_X, rho(k), mass_weight=.true.)
         !endif
          dmass_X=max(dmass_X, -mass_X/dt)
          dnumber_X=dnumber_Y
          dmass_Z=-1.0*(dmass_X + dmass_Y)
          dnumber_Z=-1.0*dnumber_X
        END IF
      END IF

      ! PRAGMATIC HACK - FIX ME (ALTHOUGH I QUITE LIKE IT)
      ! If most of the collected particles are to be removed then remove all of them
      IF (-dmass_Y*dt >0.95*mass_Y .OR. (params_Y%l_2m .AND. -dnumber_Y*dt > 0.95*number_Y)) THEN
        dmass_Y=-mass_Y/dt
        dnumber_Y=-number_Y/dt
      END IF
      ! If most of the collecting particles are to be removed then remove all of them
      IF (l_alternate_Z .AND. (-dmass_X*dt >0.95*mass_X .OR.                                          &
           (params_X%l_2m .AND. -dnumber_X*dt > 0.95*number_X))) THEN
        dmass_X=-mass_X/dt
        dnumber_X=-number_X/dt
        dmass_Z = -1.0*( dmass_X + dmass_Y )
        dnumber_Z = -1.0*dnumber_X
      END IF

      IF (.NOT. l_alternate_Z) THEN
        dmass_X=-dmass_Y
        dnumber_X=0.0
        dnumber_Z=0.0
        dmass_Z=0.0
      END IF


!convert back to grid mean
     !use mixed-phase overlap function
      IF (((params_Y%id == cloud_params%id) .OR. (params_X%id == cloud_params%id)) &
         .OR. ((params_Y%id == rain_params%id) .OR. (params_X%id == rain_params%id)) ) THEN
        overlap_cf=min(1.0,max(0.0,mpof*min(cf_X, cf_Y) +  max(0.0,(1.0-mpof)*(cf_X+cf_Y-1.0))))
      ELSE
        overlap_cf = min(cf_Y, cf_X)
      END IF
       !!overlap_cf=1.0
     
      dmass_Y = dmass_Y * overlap_cf
      dmass_X = dmass_X * overlap_cf
      IF (params_Y%l_2m) dnumber_Y = dnumber_Y * overlap_cf
      IF (params_X%l_2m) dnumber_X = dnumber_X * overlap_cf
      IF (l_alternate_Z) THEN
        dmass_Z = dmass_Z * overlap_cf
        IF (params_Z%l_2m) dnumber_Z = dnumber_Z * overlap_cf
      END IF

      procs(params_Y%i_1m, iproc%id)%column_data(k)=dmass_Y
      IF (params_Y%l_2m) THEN
        procs(params_Y%i_2m, iproc%id)%column_data(k)=dnumber_Y 
      END IF

      procs(params_X%i_1m, iproc%id)%column_data(k)=dmass_X  
      IF (params_X%l_2m) THEN
        procs(params_X%i_2m, iproc%id)%column_data(k)=dnumber_X  
      END IF

      IF (l_alternate_Z) THEN
        procs(Zi_1m, iproc%id)%column_data(k)=dmass_Z
        IF (l_Zi_2m) THEN
          procs(Zi_2m, iproc%id)%column_data(k)=dnumber_Z
        END IF
      END IF

      !----------------------
      ! Aerosol processing...
      !----------------------

      IF (l_process .AND. l_aero) THEN
        ! We note that all processes result in a source of frozen water
        ! i.e. params_Z is either ice, snow or graupel
     
        IF (params_Y%id == cloud_params%id) THEN
          dmac=abs(dnumber_Y)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_Y)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        ELSE IF (params_Y%id == rain_params%id) THEN
          dmac=abs(dnumber_Y)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_Y)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        ELSE IF (params_X%id == cloud_params%id) THEN ! This is never the case !
          dmac=abs(dnumber_X)*aeroact(k)%mact1_mean*aeroact(k)%nratio1
          dmad=abs(dnumber_X)*dustliq(k)%mact1_mean*dustliq(k)%nratio1
        ELSE IF (params_X%id == rain_params%id) THEN
          dmac=abs(dnumber_X)*aeroact(k)%mact2_mean*aeroact(k)%nratio2
          dmad=abs(dnumber_X)*dustliq(k)%mact2_mean*dustliq(k)%nratio2
        END IF

        aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
        aerosol_procs(i_am4, iaproc%id)%column_data(k)=-dmac
        aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
        aerosol_procs(i_am9, iaproc%id)%column_data(k)=-dmad
      END IF
    END IF
  END IF
END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE iacc
END MODULE ice_accretion
