! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Condensation/evaporation of cloud water and droplet activation (condevp).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Saturation-adjustment condensation/evaporation of cloud
!   water (used when CASIM's own condensation is active, e.g.
!   in MONC/KiD, rather than the UM cloud-fraction scheme).
!
! Paper reference:
!   Field et al. (2023) Appendix A.6.1 describes the
!   "all-or-nothing" saturation-adjustment scheme used when
!   CASIM's own condensation step is active; no numbered
!   equation is given for it there (Eq. A19-A20 instead cover
!   the UM cloud-fraction coupling, not this module).
!
MODULE condensation
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho, pressure, w, exner
  USE mphys_switches, ONLY: i_qv, i_ql, i_nl, i_th, i_qr, l_warm, &
       i_am4, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am6, i_an6, i_am9, i_an11, i_an12,  &
       cloud_params, l_process, l_passivenumbers,l_passivenumbers_ice, aero_index, &
       l_cfrac_casim_diag_scheme
  USE process_routines, ONLY: process_rate, i_cond, i_aact
  USE mphys_constants, ONLY: Lv, cp
  USE qsat_funs, ONLY: qsaturation, dqwsatdt
  USE thresholds, ONLY: ql_small, ss_small, thresh_tidy
  USE activation, ONLY: activate
  USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active
  USE which_mode_to_use, ONLY : which_mode
  USE casim_runtime, ONLY: casim_time, casim_smax, casim_smax_limit_time
  USE casim_parent_mod, ONLY: casim_parent, parent_um, parent_kid
  USE cloud_frac_scheme, ONLY: cloud_frac_casim_mphys

! #if DEF_MODEL==MODEL_KiD
!   use diagnostics, only: save_dg, i_dgtime, i_here, k_here
!   use runtime, only: time
!   Use namelists, only : smax, smax_limit_time
! #endif

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CONDENSATION'

  PRIVATE

      REAL(wp), ALLOCATABLE :: dnccn_all(:),dmac_all(:)
      REAL(wp), ALLOCATABLE :: dnccnd_all(:),dmad_all(:)

!$OMP THREADPRIVATE(dnccn_all, dmac_all, dnccnd_all, dmad_all)

  PUBLIC condevp_initialise, condevp_finalise, condevp
!PRF
  PUBLIC dnccn_all, dmac_all, dnccnd_all, dmad_all
!PRF

CONTAINS

  SUBROUTINE condevp_initialise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='CONDEVP_INITIALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ALLOCATE(dnccn_all(aero_index%nccn))
    ALLOCATE(dmac_all(aero_index%nccn))
    ALLOCATE(dnccnd_all(aero_index%nin))
    ALLOCATE(dmad_all(aero_index%nin))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE condevp_initialise  

  SUBROUTINE condevp_finalise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='CONDEVP_FINALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DEALLOCATE(dnccn_all)
    DEALLOCATE(dmac_all)
    DEALLOCATE(dnccnd_all)
    DEALLOCATE(dmad_all)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE condevp_finalise  

  ! All-or-nothing saturation-adjustment condensation/
  ! evaporation, Field et al. (2023) Sec. A.6.1.
  SUBROUTINE condevp(ixy_inner, dt, nz, qfields, procs, aerophys, aerochem,   &
       aeroact, dustphys, dustchem, dustliq, aerosol_procs, rhcrit_lev)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
    TYPE(aerosol_phys), INTENT(IN) :: aerophys(:)
    TYPE(aerosol_chem), INTENT(IN) :: aerochem(:)
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:)
    TYPE(aerosol_phys), INTENT(IN) :: dustphys(:)
    TYPE(aerosol_chem), INTENT(IN) :: dustchem(:)
    TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
    TYPE(process_rate), INTENT(INOUT), OPTIONAL, TARGET :: aerosol_procs(:,:)

    REAL(wp), INTENT(IN) :: rhcrit_lev(:)

    ! Local variables
    REAL(wp) :: dmass, dnumber, dmac, dmad, dnumber_a, dnumber_d
    REAL(wp) :: dmac1, dmac2, dnac1, dnac2

    REAL(wp) :: th
    REAL(wp) :: qv
    REAL(wp) :: cloud_mass
    REAL(wp) :: cloud_number

    REAL(wp) :: qs, dqsdt, qsatfac

    REAL(wp) :: tau   ! timescale for adjustment of condensate
    REAL(wp) :: w_act ! vertical velocity to use for activation

    REAL(wp) :: smax,ait_ccn, acc_ccn, tot_ccn, activated_arg, &
         activated_cloud
    ! local variables for diagnostics cloud scheme (if needed)
    REAL(wp) :: cloud_mass_new, abs_liquid_t

    REAL(wp) :: cfrac, cfrac_old

    LOGICAL :: l_docloud  ! do we want to do the calculation of cond/evap
    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='CONDEVP'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    tau=dt ! adjust instantaneously

    ! Initializations
    dmac=0.0
    dmad=0.0
    dnumber=0.0
    dnumber_a=0.0
    dnumber_d=0.0
    
    dnccn_all=0.0
    dmac_all=0.0

    DO k = 1, nz

    cfrac=1.0
    cfrac_old=1.0

    ! Set pointers for convenience
    cloud_mass=qfields(k, i_ql)
    IF (cloud_params%l_2m) cloud_number=qfields(k, i_nl)

    th=qfields(k, i_th)

    IF (casim_parent == parent_um ) THEN

      IF (l_cfrac_casim_diag_scheme ) THEN
        qv=qfields(k, i_qv)+cloud_mass
      ELSE
        qv=qfields(k, i_qv)
      END IF

      IF (l_cfrac_casim_diag_scheme) THEN
        ! work out saturation vapour pressure/mixing ratio based on
        ! liquid water temperature
        abs_liquid_T=(th*exner(k,ixy_inner))-((lv * cloud_mass )/cp)
        qs=qsaturation(abs_liquid_T, pressure(k,ixy_inner)/100.0)
      ELSE
        qs=qsaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)
      END IF
    
    ELSE ! casim_parent /= parent_um
      qv=qfields(k, i_qv)
      qs=qsaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)

    END IF ! casim_parent == parent_um

    l_docloud=.TRUE.
    IF (qs==0.0) l_docloud=.FALSE.

      
    IF (casim_parent == parent_kid) THEN
      IF ((qv/qs > 1.0 - ss_small .OR. cloud_mass > 0.0) .AND. l_docloud) THEN
!AH - following code limits the maximum supersaturation permitted. This 
!     is needed for the KiD-A 2d Sc case  - USE WITH CAUTION!!
        IF ( (((qv/qs)-1.)*100.) > casim_smax .AND. casim_time <= casim_smax_limit_time) THEN
          qs = qv/(1+(casim_smax/100.))
        END IF
      END IF
    END IF ! casim_parent == parent_kid
! AH - set the cloud fraction for new calc of activation
    IF (cloud_mass > epsilon(1.0_wp)) THEN 
       cfrac_old = 1.0_wp
    ELSE 
       cfrac_old = 0.0_wp
    END IF
       
    IF ((qv/qs > 1.0 - ss_small .OR. cloud_mass > 0.0 .OR. l_cfrac_casim_diag_scheme) .AND. l_docloud) THEN
! DPG - allow the cloud scheme to operate even if we are sub-saturated (since
! this is it's purpose!)
      IF (l_cfrac_casim_diag_scheme .AND. casim_parent == parent_um ) THEN

        !Call Smith scheme before setting up microphysics vars, to work out
        ! cloud fraction, which is used to derive in-cloud mass and number
        !
        !IMPORTANT - qv is total water at this stage!
        CALL cloud_frac_casim_mphys(k, pressure(k,ixy_inner), th*exner(k,ixy_inner), abs_liquid_T, rhcrit_lev(k),  &
             qs, qv, cloud_mass, qfields(k,i_qr), cloud_mass_new )

        dmass=max(-cloud_mass, (cloud_mass_new-cloud_mass))/dt
      ELSE
        dqsdt=dqwsatdt(qs, th*exner(k,ixy_inner))
        qsatfac=1.0/(1.0 + Lv/Cp*dqsdt)
        dmass=max(-cloud_mass, (qv-qs)*qsatfac )/dt
      END IF ! l_cfrac_casim_diag_scheme

      IF (dmass > 0.0_wp) THEN ! condensation
        IF (dmass*dt + cloud_mass > ql_small) THEN ! is it worth bothering with?
           ! AH - if dmass > 0.0 there is a change in mass, so assume cloud fraction is 1.0
           ! this assumption is only valid with all-or-nothing scheme and no cloud fraction
           ! scheme
          cfrac = 1.0_wp
          IF (cloud_params%l_2m) THEN
            ! If significant cloud formed then assume minimum velocity of 0.01m/s
            w_act=max(w(k,ixy_inner), 0.01_wp)

            CALL activate(tau, cloud_mass, cloud_number, w_act,         &
                 rho(k,ixy_inner), dnumber, dmac, th*exner(k,ixy_inner), pressure(k,ixy_inner),       &
                 cfrac, cfrac_old, aerophys(k), aerochem(k),            & 
                 aeroact(k), dustphys(k), dustchem(k), dustliq(k),      &
                 dnccn_all, dmac_all, dnumber_d, dmad,                  &
                 dnccnd_all, dmad_all, smax, ait_ccn, acc_ccn,          &
                 tot_ccn,activated_arg,activated_cloud)

            dnumber_a=dnumber
          END IF
        ELSE
          dmass=0.0 ! not worth doing anything
        END IF
      ELSE  ! evaporation
        IF (cloud_mass > thresh_tidy(i_ql)) THEN ! anything significant to remove or just noise?
          IF (dmass*dt + cloud_mass < ql_small) THEN  ! Remove all cloud
            ! Remove small quantities.
            dmass=-cloud_mass/dt
            ! liberate all number and aerosol
            IF (cloud_params%l_2m) THEN
              dnumber=-cloud_number/dt

              !============================
              ! aerosol processing
              !============================
              IF (l_process) THEN
                dmac=-aeroact(k)%mact1/dt
                dmad=-dustliq(k)%mact1/dt

                IF (l_passivenumbers) THEN
                  dnumber_a=-aeroact(k)%nact1/dt
                ELSE
                  dnumber_a=dnumber
                END IF
                IF (l_passivenumbers_ice) THEN
                  dnumber_d=-dustliq(k)%nact1/dt
                ELSE
                  dnumber_d=dnumber
                END IF

                IF (aero_index%nin > 0) THEN 
                   dmad_all(aero_index%i_coarse_dust) = dmad
                   dnccnd_all(aero_index%i_coarse_dust) = dnumber_d
                END IF

                IF (aero_index%i_accum >0 .AND. aero_index%i_coarse >0) THEN
                  ! We have both accumulation and coarse modes
                  IF (dnumber_a*dmac<= 0) THEN
                    dnumber_a=dmac/1.0e-18/dt
                  END IF
                  CALL which_mode(dmac, dnumber_a,                                 &
                       aerophys(k)%rd(aero_index%i_accum), aerophys(k)%rd(aero_index%i_coarse), &
                       aerochem(k)%density(aero_index%i_accum),     &
                       aerophys(k)%sigma(aero_index%i_accum),       &
                       dmac1, dmac2, dnac1, dnac2)

                  dmac_all(aero_index%i_accum)=dmac1  ! put it back into accumulation mode
                  dnccn_all(aero_index%i_accum)=dnac1
                  dmac_all(aero_index%i_coarse)=dmac2  ! put it back into coarse mode
                  dnccn_all(aero_index%i_coarse)=dnac2

                ELSE

                  IF (aero_index % i_accum > 0) THEN
                    dmac_all  ( aero_index % i_accum) = dmac
                    dnccn_all ( aero_index % i_accum) = dnumber_a
                  END IF

                  IF (aero_index % i_coarse > 0) THEN
                    dmac_all  (aero_index % i_coarse) = dmac
                    dnccn_all (aero_index % i_coarse) = dnumber_a
                  END IF

                END IF
              END IF
            END IF
          ELSE ! Still some cloud will be left behind
            dnumber=0.0 ! we assume no change in number during evap
            dnccn_all=0.0 ! we assume no change in number during evap
            dmac=0.0 ! No aerosol processing required
            dmac_all=0.0 ! No aerosol processing required
            dnumber_a=0.0 ! No aerosol processing required
            dnumber_d=0.0 ! No aerosol processing required
            dnccnd_all = 0.0
            dmad_all = 0.0
          END IF
        ELSE  ! Nothing significant here to remove - the tidying routines will deal with this
          dmass=0.0 ! no need to do anything since this is now just numerical noise
          dnumber=0.0 ! we assume no change in number during evap
          dnccn_all=0.0 ! we assume no change in number during evap
          dmac=0.0 ! No aerosol processing required
          dmac_all=0.0 ! No aerosol processing required
          dnumber_a=0.0 ! No aerosol processing required
          dnumber_d=0.0 ! No aerosol processing required
          dnccnd_all = 0.0
          dmad_all = 0.0
        END IF
      END IF

      IF (dmass /= 0.0_wp) THEN

        procs(i_qv, i_cond%id)%column_data(k)=-dmass
        procs(i_ql, i_cond%id)%column_data(k)=dmass

        IF (cloud_params%l_2m) THEN
          procs(i_nl, i_cond%id)%column_data(k)=dnumber
        END IF

        !============================
        ! aerosol processing
        !============================
        IF (l_process) THEN

          aerosol_procs(i_am4, i_aact%id)%column_data(k)=dmac
          IF (l_passivenumbers) aerosol_procs(i_an11, i_aact%id)%column_data(k)=dnumber_a
          IF (l_passivenumbers_ice) aerosol_procs(i_an12, i_aact%id)%column_data(k)=dnumber_d

          IF (aero_index%i_aitken > 0) THEN
            aerosol_procs(i_am1, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_aitken)
            aerosol_procs(i_an1, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_aitken)
          END IF
          IF (aero_index%i_accum > 0) THEN
            aerosol_procs(i_am2, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_accum)
            aerosol_procs(i_an2, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_accum)
          END IF
          IF (aero_index%i_coarse > 0) THEN
            aerosol_procs(i_am3, i_aact%id)%column_data(k)=-dmac_all(aero_index%i_coarse)
            aerosol_procs(i_an3, i_aact%id)%column_data(k)=-dnccn_all(aero_index%i_coarse)
          END IF

          IF (.NOT. l_warm .AND. dmad /=0.0) THEN
            ! We may have some dust in the liquid...
             IF (aero_index%nin > 0 ) THEN 
                aerosol_procs(i_am9, i_aact%id)%column_data(k)= dmad_all(aero_index%i_coarse_dust)
                aerosol_procs(i_am6, i_aact%id)%column_data(k)=-dmad_all(aero_index%i_coarse_dust)! < USING COARSE
                aerosol_procs(i_an6, i_aact%id)%column_data(k)=-dnccnd_all(aero_index%i_coarse_dust) ! < USING COARSE
             END IF
          END IF
        END IF
       END IF
     END IF 
   END DO
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE condevp
END MODULE condensation
