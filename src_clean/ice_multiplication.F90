! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Secondary ice production: Hallet-Mossop rime splintering, droplet shattering, and ice-ice collisional fragmentation.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE ice_multiplication
  USE variable_precision, ONLY: wp, iwp
  USE process_routines, ONLY: process_rate, process_name, i_gacw, i_sacw,      &
                                 i_ihal, i_idps, i_iics, i_gaci, i_gacs, i_homr
  USE passive_fields, ONLY: TdegC, TdegK
  USE mphys_parameters, ONLY: ice_params, snow_params, graupel_params,         &
                       dN_hallet_mossop, M0_hallet_mossop, dN_droplet_shatter, &
                       P_droplet_shatter, coef_ice_breakup
  USE thresholds, ONLY: thresh_small, cfliq_small
  USE m3_incs, ONLY: m3_inc_type2
  USE mphys_switches, ONLY: l_prf_cfrac, i_cfs, i_cfg, i_cfl, mpof, i_qg,      &
                            i_ng, i_cfi
  USE mphys_constants, ONLY: pi
  USE casim_stph, ONLY: l_rp2_casim, mpof_casim_rp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ICE_MULTIPLICATION'

CONTAINS
  !> Subroutine to determine the ice splintering by Hallet-Mossop
  !> This effect requires prior calculation of the accretion rate of
  !> graupel and snow.
  !> This is a source of ice number and mass and a sink of liquid
  !> (but this is done via the accretion processes already so is
  !> represented here as a sink of snow/graupel)
  !> For triple moment species there is a corresponding change in the
  !> 3rd moment assuming shape parameter is not changed
  !>
  !> Subroutine to determine the droplet shattering
  !> This effect requires prior calculation of the Bigg freezing
  !> rate for raindrops.
  !> This is a source of ice number and mass, and a sink of graupel
  !> number and mass (sink for snow number and mass already done in 
  !> the Bigg's raindrop freezing).
  !>
  !> Subroutine to determine the ice-ice collision
  !> This effect requires prior calculation of collision tendency
  !> rates for graupel-snow accretion and snow collecting snow.
  !> This is a source of snow and ice number concentration, and not 
  !> changing the graupel mass and number.
  !>
  !> AEROSOL: All aerosol sinks/sources are assumed to come from soluble modes
  !
  !> OPTIMISATION POSSIBILITIES:
  SUBROUTINE hallet_mossop(ixy_inner, dt, nz, cffields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN) :: cffields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! Local variables
    REAL(wp) :: gacw, sacw  ! accretion process rates
    REAL(wp) :: dnumber_s, dnumber_g  ! number conversion rate from snow/graupel
    REAL(wp) :: dmass_s, dmass_g      ! mass conversion rate from snow/graupel
    REAL(wp) :: Eff  !< splintering efficiency
    REAL(wp) :: cf_snow, cf_graupel, cf_liquid, overlap_cfsnow, overlap_cfgraupel

    INTEGER :: k

    TYPE(process_name) :: iproc ! processes selected depending on which species we're modifying

    CHARACTER(len=*), PARAMETER :: RoutineName='HALLET_MOSSOP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    IF ( l_rp2_casim ) THEN
        mpof = mpof_casim_rp
    END IF

    IF (.NOT. ice_params%l_2m) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    
    DO k = 1, nz
       IF (TdegC(k,ixy_inner) < 0.0_wp) THEN 
          IF (l_prf_cfrac) THEN
             IF (cffields(k,i_cfl) > cfliq_small) THEN
                cf_liquid=cffields(k,i_cfl)
             ELSE
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfs) > cfliq_small) THEN
                cf_snow=cffields(k,i_cfs)
             ELSE
                cf_snow=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfg) > cfliq_small) THEN
                cf_graupel=cffields(k,i_cfg)
             ELSE
                cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
          ELSE
             cf_snow=1.0
             cf_graupel=1.0
             cf_liquid=1.0
          END IF

          !use mixed-phase overlap function
          overlap_cfsnow=min(1.0,max(0.0,mpof*min(cf_liquid, cf_snow) +         &
               max(0.0,(1.0-mpof)*(cf_liquid+cf_snow-1.0))))
          overlap_cfgraupel=min(1.0,max(0.0,mpof*min(cf_liquid, cf_graupel) +   &
               max(0.0,(1.0-mpof)*(cf_liquid+cf_graupel-1.0))))

          Eff=1.0 - abs(TdegC(k,ixy_inner) + 5.0)/2.5 ! linear increase between -2.5/-7.5 and -5C

          IF (Eff > 0.0) THEN
             sacw=0.0
             gacw=0.0
             !! should use cf_overlap as in ice accretion
             IF (snow_params%i_1m > 0) &
                  sacw=procs(snow_params%i_1m, i_sacw%id)%column_data(k)/overlap_cfsnow  !insnow process rate
             IF (graupel_params%i_1m > 0) &
                  gacw=procs(graupel_params%i_1m, i_gacw%id)%column_data(k)/overlap_cfgraupel ! ingraupel process rate
             
             IF ((sacw*overlap_cfsnow + gacw*overlap_cfgraupel)*dt > thresh_small(snow_params%i_1m)) THEN
                iproc=i_ihal

                dnumber_g=dN_hallet_mossop * Eff * (gacw) ! Number of splinters from graupel
                dnumber_s=dN_hallet_mossop * Eff * (sacw) ! Number of splinters from snow
                
                dnumber_g=min(dnumber_g, 0.5*gacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid
                dnumber_s=min(dnumber_s, 0.5*sacw/M0_hallet_mossop) ! don't remove more than 50% of rimed liquid

                dmass_g=dnumber_g * M0_hallet_mossop * overlap_cfgraupel  ! convert back to grid mean
                dmass_s=dnumber_s * M0_hallet_mossop * overlap_cfsnow ! convert back to grid mean
                
                dnumber_g=dnumber_g * overlap_cfgraupel  ! convert back to grid mean
                dnumber_s=dnumber_s * overlap_cfsnow  ! convert back to grid mean
        

                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_g + dmass_s 
                procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_g + dnumber_s
                
                !-------------------
                ! Sinks for snow...
                !-------------------
                IF (sacw > 0.0) THEN
                   procs(snow_params%i_1m, iproc%id)%column_data(k)=-dmass_s
                   procs(snow_params%i_2m, iproc%id)%column_data(k)=0.0
                END IF
                
                !---------------------
                ! Sinks for graupel...
                !---------------------
                IF (gacw > 0.0) THEN
                   procs(graupel_params%i_1m, iproc%id)%column_data(k)=-dmass_g
                   procs(graupel_params%i_2m, iproc%id)%column_data(k)=0.0
                END IF
                
             END IF
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE hallet_mossop

!--------------------------------------------------------------------------------
!> Subroutine for droplet shattering (Sullivan et al., 2018)
!--------------------------------------------------------------------------------
  SUBROUTINE droplet_shattering(ixy_inner, dt, nz, cffields, qfields, procs)

   USE yomhook, ONLY: lhook, dr_hook
   USE parkind1, ONLY: jprb, jpim

   IMPLICIT NONE

   ! Subroutine arguments
   INTEGER, INTENT(IN) :: ixy_inner
   REAL(wp), INTENT(IN) :: dt
   INTEGER, INTENT(IN) :: nz
   REAL(wp), INTENT(IN) :: cffields(:,:)
   REAL(wp), INTENT(IN) :: qfields(:,:)
   TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

   ! Local variables
   REAL(wp) :: homr_mass, homr_number  ! rate of homogeneous freezing of rain (Bigg, 1953)
   REAL(wp) :: dnumber_i, dnumber_g  ! number conversion rate for ice crystal and graupel
   REAL(wp) :: dmass_i, dmass_g      ! mass conversion rate for ice crystal and graupel
   REAL(wp) :: prob_DS ! temperature-dependent shattering probability
   REAL(wp) :: cf_graupel, graupel_mass, graupel_number

   INTEGER :: k

   TYPE(process_name) :: iproc ! processes selected depending on which species we're modifying

   CHARACTER(len=*), PARAMETER :: RoutineName='DROPLET_SHATTERING'

   INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
   INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
   REAL(KIND=jprb)               :: zhook_handle

   !--------------------------------------------------------------------------
   ! End of header, no more declarations beyond here
   !--------------------------------------------------------------------------
   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

   IF (.NOT. ice_params%l_2m) THEN
     IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
     RETURN
   END IF
   
   DO k = 1, nz
      IF (TdegC(k,ixy_inner) < 0.0_wp) THEN 
         IF (l_prf_cfrac) THEN
            IF (cffields(k,i_cfg) > cfliq_small) THEN
               cf_graupel=cffields(k,i_cfg)
            ELSE
               cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
            END IF
         ELSE
            cf_graupel=1.0
         END IF
         
         ! Shattering Probability:
         ! Normal distribution centred at 258 K with standard deviation of 3 K
         ! probablity of droplet shatter P_Droplet_shatter defaults to 0.2
         ! maximum of distribution is 0.13298 (Sullivan (2018))
         prob_DS= (P_droplet_shatter / 0.13298) * (1 / (SQRT(2 * pi) * 3))          &
                 * EXP((-(TdegK(k,ixy_inner) - 258)**2) / (18))
         
         IF (prob_DS > 0.0) THEN
            homr_mass = 0.0   ! mass tendency of raindrop frozen
            homr_number = 0.0   ! number tendency of raindrop frozen

            graupel_mass=qfields(k,i_qg)
            graupel_number=qfields(k,i_ng)
            
            IF (graupel_params%i_1m > 0) &
                 homr_mass=procs(graupel_params%i_1m, i_homr%id)%column_data(k)/cf_graupel ! ingraupel process rate (kg / kg-1)

            IF (graupel_params%i_2m > 0) &
                 homr_number=procs(graupel_params%i_2m, i_homr%id)%column_data(k)/cf_graupel ! ingraupel process rate (number / kg-1)

            IF ((homr_mass*cf_graupel)*dt > thresh_small(graupel_params%i_1m)  &
                .AND. (homr_number*cf_graupel)*dt > thresh_small(graupel_params%i_2m)) THEN 
               
               dnumber_i=(1 + prob_DS * dN_droplet_shatter) * homr_number ! Number of splinters from graupel
               ! No more than 50% of the graupels created from the frozen raindrops
               dmass_i=min(dnumber_i * M0_hallet_mossop,0.5*homr_mass) 

               dnumber_i=dnumber_i * cf_graupel ! Convert back to grid-box mean
               dmass_i=dmass_i * cf_graupel

               dmass_g=-dmass_i
               dnumber_g=dmass_g * graupel_number / graupel_mass 

               IF (homr_mass > 0.0) THEN
                  iproc = i_idps
                  !-------------------
                  ! Sources for ice...
                  !-------------------
                  procs(ice_params%i_1m, iproc%id)%column_data(k)=dmass_i
                  procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_i
                  
                  !---------------------
                  ! Sinks for graupel...
                  !---------------------
                  procs(graupel_params%i_1m, iproc%id)%column_data(k)=dmass_g
                  procs(graupel_params%i_2m, iproc%id)%column_data(k)=dnumber_g
               END IF
            END IF
         END IF
      END IF
   END DO

   IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

 END SUBROUTINE droplet_shattering

!-----------------------------------------------------------------------------------
!> Subroutine for ice-ice collision

 SUBROUTINE ice_collision(ixy_inner, dt, nz, cffields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner    
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN) :: cffields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! Local variables
    REAL(wp) :: gaci, gacs  ! accretion process rates
    REAL(wp) :: dnumber_i, dnumber_s ! number tendency for the collided hydrometeors
    REAL(wp) :: BR_fragments ! Temperature-dependent fragments from ice-ice collision
    REAL(wp) :: cf_snow, cf_graupel, cf_ice, overlap_cfsg, overlap_cfig

    INTEGER :: k

    TYPE(process_name) :: iproc ! processes selected depending on which species we're modifying

    CHARACTER(len=*), PARAMETER :: RoutineName='ICE_COLLISION'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (.NOT. ice_params%l_2m) THEN
      IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
      RETURN
    END IF
    
    DO k = 1, nz
       IF (TdegC(k,ixy_inner) < 0.0_wp) THEN 
          IF (l_prf_cfrac) THEN
             IF (cffields(k,i_cfs) > cfliq_small) THEN
                cf_snow=cffields(k,i_cfs)
             ELSE
                cf_snow=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfg) > cfliq_small) THEN
                cf_graupel=cffields(k,i_cfg)
             ELSE
                cf_graupel=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfi) > cfliq_small) THEN
               cf_ice=cffields(k,i_cfi)
            ELSE
               cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
            END IF
          ELSE
             cf_snow=1.0
             cf_graupel=1.0
             cf_ice=1.0
          END IF
         
          overlap_cfsg = min(cf_snow, cf_graupel)
          overlap_cfig = min(cf_ice, cf_graupel)
           
          ! Number of fragments generated based on Takahashi et al., (1995)
          BR_fragments=coef_ice_breakup * ((TdegK(k,ixy_inner) - 252) ** 1.2)      &
                      * EXP(-(TdegK(k,ixy_inner) - 252)/5)

          IF (BR_fragments > 0.0) THEN
             gacs=0.0
             gaci=0.0
             
             IF (snow_params%i_2m > 0) &
                  gacs=-procs(snow_params%i_2m, i_gacs%id)%column_data(k)/overlap_cfsg  ! insnow process rate
             IF (graupel_params%i_2m > 0) &
                  gaci=-procs(ice_params%i_2m, i_gaci%id)%column_data(k)/overlap_cfig   ! inice process rate
             
             IF ((gacs*cf_snow)*dt > thresh_small(snow_params%i_2m)) THEN
                iproc=i_iics

                dnumber_s=BR_fragments * (gacs) * overlap_cfsg   ! Number of splinters from graupel and convert back to grid mean
                !-------------------
                ! Sources for snow...
                !-------------------
                procs(snow_params%i_2m, iproc%id)%column_data(k)=dnumber_s                
             END IF

             IF ((gaci*cf_ice)*dt > thresh_small(ice_params%i_2m)) THEN
                iproc=i_iics

                dnumber_i=BR_fragments * (gaci) * overlap_cfig   ! Number of splinters from graupel and convert back to grid mean
                !-------------------
                ! Sources for ice...
                !-------------------
                procs(ice_params%i_2m, iproc%id)%column_data(k)=dnumber_i
             END IF

          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ice_collision

END MODULE ice_multiplication
