! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Main CASIM entry point (shipway_microphysics): orchestrates per-column microphysics over a full 3D grid.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Top-level microphysics driver (shipway_microphysics /
!   microphysics_common): loops over grid columns/levels,
!   preconditions state, dispatches every enabled physics
!   process in sequence, and gathers diagnostics.
!
! Paper reference:
!   Corresponds to the overall CASIM time-step structure shown
!   schematically in Field et al. (2023) Figure A1 and described
!   throughout Appendix A; the specific dispatch/loop logic here
!   is CASIM software infrastructure, not itself a numbered
!   equation.
!
MODULE micro_main
  USE variable_precision, ONLY: wp
  USE mphys_parameters, ONLY: nz, nq, rain_params, cloud_params, ice_params, &
       snow_params, graupel_params, nspecies, ZERO_REAL_WP, a_s, b_s, &
       nxy_inner
  USE process_routines, ONLY: process_rate, zero_procs, allocate_procs, deallocate_procs, i_cond, i_praut, &
       i_pracw, i_pracr, i_prevp, i_psedr, i_psedl, i_aact, i_aaut, i_aacw, i_aevp, i_asedr, i_asedl, i_arevp, &
       i_tidy2, i_atidy2, i_inuc, i_idep, i_dnuc, i_dsub, i_saut, i_iacw, i_sacw, i_pseds, &
       i_sdep, i_saci, i_raci, i_sacr, i_gacw, i_gacr, i_gaci, i_gacs, i_gdep, i_psedg, i_sagg, &
       i_gshd, i_ihal, i_smlt, i_gmlt, i_psedi, i_homr, i_homc, i_imlt, i_isub, i_ssub, i_gsub, i_sbrk, i_dssub, &
       i_dgsub, i_dsedi, i_dseds, i_dsedg, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, &
       i_dgacr, i_draci, i_dhomr, i_dhomc, i_idps, i_iics
  USE sum_process, ONLY: sum_procs, sum_aprocs, tend_temp, aerosol_tend_temp
  USE aerosol_routines, ONLY: examine_aerosol, aerosol_phys, aerosol_chem, aerosol_active, allocate_aerosol, &
       deallocate_aerosol
  USE mphys_switches, ONLY: hydro_complexity, aero_complexity, i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, i_qi, &
       i_qs, i_qg, i_ni, i_ns, i_ng, i_m3s, i_m3g, i_am1, i_an1, i_am2, i_an2, i_am3, i_an3, i_am4, i_am5, i_am6, &
       i_an6, i_am7, i_am8 , i_am9, i_am10, i_an10, i_an11, i_an12, i_ak1, i_ak2, i_ak3, &
       aerosol_option, l_warm, l_passivenumbers, l_passivenumbers_ice, &
       l_sed, l_idep, aero_index, nq_l, nq_r, nq_i, nq_s, nq_g, &
       l_sg, l_g, l_process, max_sed_length, max_step_length, l_harrington, l_passive, ntotala, ntotalq, &
       l_onlycollect, pswitch, aswitch, l_isub, l_pos1, l_pos2, l_pos3, l_pos4, l_no_pgacs_in_sumprocs, &
       l_pos5, l_pos6, i_hstart, l_tidy_negonly, l_separate_rain,  &
       iopt_act, iopt_shipway_act, l_prf_cfrac, l_kfsm, l_gamma_online, l_subseds_maxv, &
       i_cfl, i_cfr, i_cfi, i_cfs, i_cfg, l_reisner_graupel_embryo
! use mphys_switches, only: l_rain,
  USE passive_fields, ONLY: rexner, min_dz
  USE mphys_constants, ONLY: cp, Lv
  USE distributions, ONLY: query_distributions, initialise_distributions, dist_lambda, dist_mu, dist_n0, dist_lams
  USE passive_fields, ONLY: initialise_passive_fields, set_passive_fields, TdegK, rhcrit_1d
  USE autoconversion, ONLY: raut
  USE evaporation, ONLY: revp
  USE condensation, ONLY: condevp_initialise, condevp_finalise, condevp
  USE accretion, ONLY: racw
  USE aggregation, ONLY: racr, ice_aggregation
  USE sedimentation, ONLY: sedr, sedr_1M_2M, terminal_velocity_CFL
  USE ice_nucleation, ONLY: inuc
  USE ice_deposition, ONLY: idep
  USE ice_accretion, ONLY: iacc
  USE breakup, ONLY: ice_breakup
  USE snow_autoconversion, ONLY: saut
  USE ice_multiplication, ONLY: hallet_mossop, droplet_shattering, ice_collision
  USE graupel_wetgrowth, ONLY: wetgrowth
  USE graupel_embryo, ONLY: graupel_embryos
  USE ice_melting, ONLY: melting
  USE homogeneous, ONLY: ihom_rain, ihom_droplets
  USE adjust_deposition, ONLY: adjust_dep
  USE mphys_constants, ONLY: fixed_aerosol_sigma, fixed_aerosol_density
!AJM  removing line below causes model to fail
   USE lookup, ONLY: get_slope_generic

  USE mphys_tidy, ONLY: initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, &
       ensure_saturated, tidy_qin, tidy_ain, ensure_positive_aerosol
  USE preconditioning, ONLY: precondition, preconditioner

  USE casim_reflec_mod, ONLY: casim_reflec, setup_reflec_constants
  ! For initialization of Shipway (2015) activation scheme
  USE shipway_lookup, ONLY: generate_tables

  USE generic_diagnostic_variables, ONLY: casdiags

  USE casim_runtime, ONLY: casim_time
  USE casim_parent_mod, ONLY: casim_parent, parent_monc

  USE casim_stph, ONLY: l_rp2_casim, snow_a_x_rp, ice_a_x_rp

! #if DEF_MODEL==MODEL_KiD
!   ! Kid modules
!   use diagnostics, only: save_dg, i_dgtime, n_sub, n_subsed
!   use runtime, only: time
!   use parameters, only: nx
!   Use namelists, only : no_precip_time, l_sediment
! #endif


  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='MICRO_MAIN'

  LOGICAL :: l_tendency_loc
  LOGICAL :: l_warm_loc

!$OMP THREADPRIVATE(l_tendency_loc, l_warm_loc)

  INTEGER :: i_start, i_end ! upper and lower i levels which are to be used
  INTEGER :: j_start, j_end ! upper and lower j levels
  INTEGER :: k_start, k_end ! upper and lower k levels

!$OMP THREADPRIVATE(i_start,i_end,j_start,j_end,k_start,k_end)
!  integer :: nxny
  REAL(wp), ALLOCATABLE, SAVE :: precip(:,:) ! diagnostic for surface precip rate

!--Add one more dimention for arrays to be used in ixy_inner loop--!
  REAL(wp), ALLOCATABLE :: dqfields(:,:,:), qfields(:,:,:), tend(:,:,:)
  REAL(wp), ALLOCATABLE :: daerofields(:,:,:), aerofields(:,:,:), aerosol_tend(:,:,:)
  REAL(wp), ALLOCATABLE :: cffields(:,:,:) !cloudfraction fields

!$OMP THREADPRIVATE(precip, dqfields, qfields, cffields, tend,                   &
!$OMP               daerofields, aerofields, aerosol_tend)

  TYPE(process_rate), ALLOCATABLE :: procs(:,:,:)
  TYPE(process_rate), ALLOCATABLE :: aerosol_procs(:,:,:)

!$OMP THREADPRIVATE(procs, aerosol_procs)

  TYPE(aerosol_active), ALLOCATABLE :: aeroact(:)
  TYPE(aerosol_phys), ALLOCATABLE   :: aerophys(:)
  TYPE(aerosol_chem), ALLOCATABLE   :: aerochem(:)

!$OMP THREADPRIVATE(aeroact, aerophys, aerochem)

  TYPE(aerosol_active), ALLOCATABLE :: dustact(:)
  TYPE(aerosol_phys), ALLOCATABLE   :: dustphys(:)
  TYPE(aerosol_chem), ALLOCATABLE   :: dustchem(:)

!$OMP THREADPRIVATE(dustact, dustphys, dustchem)

  TYPE(aerosol_active), ALLOCATABLE :: aeroice(:)  ! Soluble aerosol in ice
  TYPE(aerosol_active), ALLOCATABLE :: dustliq(:)! Insoluble aerosol in liquid

!$OMP THREADPRIVATE(aeroice, dustliq)

  REAL(wp), ALLOCATABLE :: qfields_in(:,:,:)
  REAL(wp), ALLOCATABLE :: qfields_mod(:,:,:)
  REAL(wp), ALLOCATABLE :: aerofields_in(:,:,:)
  REAL(wp), ALLOCATABLE :: aerofields_mod(:,:,:)

!$OMP THREADPRIVATE(qfields_in, qfields_mod, aerofields_in, aerofields_mod)

  LOGICAL, ALLOCATABLE :: l_Tcold(:,:) ! temperature below freezing, i.e. .not. l_Twarm
  LOGICAL, ALLOCATABLE :: l_sigevap(:,:) ! logical to determine significant evaporation

!$OMP THREADPRIVATE(l_Tcold, l_sigevap)

!Arrays needed by inner loop

  INTEGER, ALLOCATABLE :: columns(:)

!$OMP THREADPRIVATE(columns)

  REAL :: DTPUD ! Time step for puddle diagnostic

  PUBLIC initialise_micromain, finalise_micromain, shipway_microphysics, DTPUD

!PRF
  PUBLIC aerophys, aeroact, aerochem, dustact, dustphys, dustchem, aeroice, dustliq
!PRF
CONTAINS

  SUBROUTINE initialise_micromain(il, iu, jl, ju, kl, ku,                 &
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='INITIALISE_MICROMAIN'

    INTEGER, INTENT(IN) :: il, iu ! upper and lower i levels
    INTEGER, INTENT(IN) :: jl, ju ! upper and lower j levels
    INTEGER, INTENT(IN) :: kl, ku ! upper and lower k levels

    INTEGER, INTENT(IN) :: is_in, ie_in ! upper and lower i levels which are to be used
    INTEGER, INTENT(IN) :: js_in, je_in ! upper and lower j levels
    INTEGER, INTENT(IN) :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    LOGICAL, INTENT(IN) :: l_tendency

    ! Local variables

    INTEGER :: k

    INTEGER :: nprocs     ! number of process rates stored
    INTEGER :: naeroprocs ! number of process rates stored
    INTEGER :: naero      ! number of aerosol fields

    REAL(wp) :: beta_init

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    l_warm_loc=l_warm ! Original setting

    i_start=is_in
    i_end=ie_in
    j_start=js_in
    j_end=je_in
    k_start=ks_in
    k_end=ke_in

    !nxy_inner = (i_end-i_start+1)*(j_end-j_start+1) ! test

    ALLOCATE(rhcrit_1d(kl:ku))
    ! Set RHCrit to 1.0 as default; parent model can then overwrite
    ! this if needed
    rhcrit_1d(:) = 1.0

    l_tendency_loc = l_tendency

    nq=sum(hydro_complexity%nmoments)+2 ! also includes vapour and theta
    nz=k_end-k_start+1
    nprocs = hydro_complexity%nprocesses


    ALLOCATE(precondition(nz, nxy_inner))
    precondition=.TRUE. ! Assume all points need to be considered
    ALLOCATE(l_Tcold(nz, nxy_inner))
    l_Tcold =.FALSE. ! Assumes no cold points, this is set at the beginning of microphysics_common
    ALLOCATE(l_sigevap(nz, nxy_inner))
    l_sigevap = .FALSE.
    ALLOCATE(qfields(nz, nq, nxy_inner))
    ALLOCATE(dqfields(nz, nq, nxy_inner))
    qfields=ZERO_REAL_WP
    dqfields=ZERO_REAL_WP
    !allocate(procs(nz, nprocs))
    ALLOCATE(procs(ntotalq, nprocs, nxy_inner))
    ALLOCATE(tend(nz, nq, nxy_inner))
    ALLOCATE(tend_temp(nz,nq))
    ALLOCATE(cffields(nz,5, nxy_inner)) !5 'cloud' fractions
    cffields=ZERO_REAL_WP

    ! Allocate aerosol storage
    IF (aerosol_option > 0) THEN
      naero=ntotala
      naeroprocs=aero_complexity%nprocesses
      ALLOCATE(aerofields(nz, naero, nxy_inner))
      ALLOCATE(daerofields(nz, naero, nxy_inner))
      aerofields=ZERO_REAL_WP
      daerofields=ZERO_REAL_WP
      ! allocate(aerosol_procs(nz, naeroprocs))
      ALLOCATE(aerosol_procs(ntotala, naeroprocs, nxy_inner))
      ALLOCATE(aerosol_tend(nz, naero, nxy_inner))
      ALLOCATE(aerosol_tend_temp(nz, naero))
    ELSE
      ! Dummy arrays required
      ALLOCATE(aerofields(1,1,nxy_inner))
      ALLOCATE(daerofields(1,1,nxy_inner))
      ALLOCATE(aerosol_procs(1,1,nxy_inner))
      ALLOCATE(aerosol_tend(1,1,nxy_inner))
      ALLOCATE(aerosol_tend_temp(1,1))
    END IF

    ALLOCATE(aerophys(nz))
    ALLOCATE(aerochem(nz))
    ALLOCATE(aeroact(nz))

    CALL allocate_aerosol(aerophys, aerochem, aero_index%nccn)
    ALLOCATE(dustphys(nz))
    ALLOCATE(dustchem(nz))
    ALLOCATE(dustact(nz))
    CALL allocate_aerosol(dustphys, dustchem, aero_index%nin)

    ALLOCATE(aeroice(nz))
    ALLOCATE(dustliq(nz))

    ! Preserve initial values for non-Shipway activation
    IF ( iopt_act == iopt_shipway_act ) THEN
      beta_init = 0.5
    ELSE
      beta_init = 1.0
    END IF

    ! Temporary initialization of chem and sigma
    DO k =1,size(aerophys)
      aerophys(k)%sigma(:)=fixed_aerosol_sigma
      aerophys(k)%rpart(:)=0.0
      aerochem(k)%vantHoff(:)=3.0
      aerochem(k)%massMole(:)=132.0e-3
      aerochem(k)%density(:)=fixed_aerosol_density
      aerochem(k)%epsv(:)=1.0
      aerochem(k)%beta(:)=beta_init
    END DO
    DO k =1,size(dustphys)
      dustphys(k)%sigma(:)=fixed_aerosol_sigma
      dustphys(k)%rpart(:)=0.0
      dustchem(k)%vantHoff(:)=3.0
      dustchem(k)%massMole(:)=132.0e-3
      dustchem(k)%density(:)=fixed_aerosol_density
      dustchem(k)%epsv(:)=1.0
      dustchem(k)%beta(:)=beta_init
    END DO

    !allocate space for the process rates
!    do ixy = 1, nxy_inner ! push ixy loop into allocate_procs
    CALL allocate_procs(nxy_inner, procs, nz, nprocs, ntotalq)
    IF (l_process) CALL allocate_procs(nxy_inner, aerosol_procs, nz, naeroprocs, ntotala)
!    end do


    ! allocate diagnostics
    ALLOCATE(precip(il:iu,jl:ju))

    CALL initialise_passive_fields(k_start, k_end)

    ALLOCATE(qfields_in(nz, nq, nxy_inner))
    ALLOCATE(qfields_mod(nz, nq, nxy_inner))
    IF (aerosol_option > 0) THEN
      ALLOCATE(aerofields_in(nz, naero, nxy_inner))
      ALLOCATE(aerofields_mod(nz, naero, nxy_inner))
    END IF
    CALL initialise_distributions(nz, nspecies)
    CALL initialise_mphystidy()
    CALL condevp_initialise()
! Here we initialise some things for the Shipway(2015) aerosol activation.

!$OMP SINGLE
    IF ( iopt_act == iopt_shipway_act ) THEN
      CALL generate_tables()
    END IF

    CALL setup_reflec_constants()
!$OMP END SINGLE

! Array needed for inner loop
    ALLOCATE(columns(nxy_inner))


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE initialise_micromain

  SUBROUTINE finalise_micromain()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='FINALISE_MICROMAIN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!    nxy_inner = 2

    ! deallocate diagnostics
    DEALLOCATE(precip)

    ! deallocate process rates
    CALL deallocate_procs(nxy_inner, procs)

    DEALLOCATE(procs)
    DEALLOCATE(qfields)
    DEALLOCATE(tend)
    DEALLOCATE(tend_temp)
    DEALLOCATE(precondition)
    DEALLOCATE(cffields)

    ! aerosol fields
    IF (l_process) CALL deallocate_procs(nxy_inner, aerosol_procs)

    DEALLOCATE(dustliq)
    DEALLOCATE(aeroice)

    CALL deallocate_aerosol(aerophys, aerochem)
    DEALLOCATE(aerophys)
    DEALLOCATE(aerochem)
    DEALLOCATE(aeroact)
    CALL deallocate_aerosol(dustphys, dustchem)
    DEALLOCATE(dustphys)
    DEALLOCATE(dustchem)
    DEALLOCATE(dustact)
    DEALLOCATE(aerosol_procs)
    DEALLOCATE(aerosol_tend)
    DEALLOCATE(aerosol_tend_temp)
    DEALLOCATE(aerofields)
    DEALLOCATE(daerofields)
    DEALLOCATE(rhcrit_1d)
    DEALLOCATE(qfields_in)
    DEALLOCATE(qfields_mod)
    IF (allocated(aerofields_in)) DEALLOCATE(aerofields_in)
    IF (allocated (aerofields_mod)) DEALLOCATE(aerofields_mod)
    DEALLOCATE(dist_lambda)
    DEALLOCATE(dist_mu)
    DEALLOCATE(dist_n0)
    DEALLOCATE(dist_lams)
    DEALLOCATE(a_s)
    DEALLOCATE(b_s)
    CALL finalise_mphystidy()
    CALL condevp_finalise()

! Array needed for inner loop
    DEALLOCATE(columns)


    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE finalise_micromain

  SUBROUTINE shipway_microphysics(il, iu, jl, ju, kl, ku, dt,               &
       qv, q1, q2, q3, q4, q5, q6, q7, q8, q9, q10, q11, q12, q13,          &
       theta, a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11, a12, a13,       &
       a14, a15, a16, a17, a18, a19, a20,                                   &
       exner, pressure, rho, w, tke, dz,                                    &
       cfliq, cfice, cfsnow, cfrain, cfgr,    &
       dqv, dq1, dq2, dq3, dq4, dq5, dq6, dq7, dq8, dq9, dq10, dq11, dq12,  &
       dq13, dth, da1, da2, da3, da4, da5, da6, da7, da8, da9, da10, da11,  &
       da12, da13, da14, da15, da16, da17,                                  &
       is_in, ie_in, js_in, je_in)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SHIPWAY_MICROPHYSICS'

    INTEGER, INTENT(IN) :: il, iu ! upper and lower i levels
    INTEGER, INTENT(IN) :: jl, ju ! upper and lower j levels
    INTEGER, INTENT(IN) :: kl, ku ! upper and lower k levels

    REAL(wp), INTENT(IN) :: dt    ! parent model timestep (s)

    ! hydro fields in... 1-5 should be warm rain, 6+ are ice
    ! see mphys_casim for details of what is passed in
    REAL(wp), INTENT(IN) :: q1( kl:ku, il:iu, jl:ju ), q2( kl:ku, il:iu, jl:ju )   &
         , q3( kl:ku, il:iu, jl:ju ), q4( kl:ku, il:iu, jl:ju ), q5( kl:ku, il:iu, jl:ju ) &
         , q6( kl:ku, il:iu, jl:ju ), q7( kl:ku, il:iu, jl:ju ), q8( kl:ku, il:iu, jl:ju ) &
         , q9( kl:ku, il:iu, jl:ju ), q10( kl:ku, il:iu, jl:ju ), q11( kl:ku, il:iu, jl:ju ) &
         , q12( kl:ku, il:iu, jl:ju ), q13( kl:ku, il:iu, jl:ju )

    REAL(wp) :: cfliq(kl:ku, il:iu, jl:ju ), cfrain(kl:ku, il:iu, jl:ju ), cfice(kl:ku, il:iu, jl:ju ), &
                cfsnow(kl:ku, il:iu, jl:ju ), cfgr(kl:ku, il:iu, jl:ju )



    REAL(wp), INTENT(IN) :: qv( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: theta( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: exner( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: pressure( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: rho( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: w( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: tke( kl:ku, il:iu, jl:ju )
    REAL(wp), INTENT(IN) :: dz( kl:ku, il:iu, jl:ju )

    ! Aerosol fields in
    REAL(wp), INTENT(IN) :: a1( kl:ku, il:iu, jl:ju ), a2( kl:ku, il:iu, jl:ju )   &
         , a3( kl:ku, il:iu, jl:ju ), a4( kl:ku, il:iu, jl:ju ), a5( kl:ku, il:iu, jl:ju ) &
         , a6( kl:ku, il:iu, jl:ju ), a7( kl:ku, il:iu, jl:ju ), a8( kl:ku, il:iu, jl:ju ) &
         , a9( kl:ku, il:iu, jl:ju ), a10( kl:ku, il:iu, jl:ju ), a11(kl:ku, il:iu, jl:ju ) &
         , a12( kl:ku, il:iu, jl:ju ), a13( kl:ku, il:iu, jl:ju ), a14( kl:ku, il:iu, jl:ju ) &
         , a15( kl:ku, il:iu, jl:ju ), a16( kl:ku, il:iu, jl:ju ), a17(kl:ku,il:iu, jl:ju )  &
         , a18( kl:ku, il:iu, jl:ju ), a19( kl:ku, il:iu, jl:ju ), a20( kl:ku,il:iu, jl:ju )

    ! hydro tendencies in:  from parent model forcing i.e. advection
    ! hydro tendencies out: from microphysics only...
    REAL(wp), INTENT(INOUT) :: dq1( kl:ku, il:iu, jl:ju ), dq2( kl:ku, il:iu, jl:ju ) &
         , dq3( kl:ku, il:iu, jl:ju ), dq4( kl:ku, il:iu, jl:ju ), dq5( kl:ku, il:iu, jl:ju ) &
         , dq6( kl:ku, il:iu, jl:ju ), dq7( kl:ku, il:iu, jl:ju ), dq8( kl:ku, il:iu, jl:ju ) &
         , dq9( kl:ku, il:iu, jl:ju ), dq10( kl:ku, il:iu, jl:ju ), dq11( kl:ku, il:iu, jl:ju ) &
         , dq12( kl:ku, il:iu, jl:ju ), dq13( kl:ku, il:iu, jl:ju )

    ! qv/theta tendencies in:  from parent model forcing i.e. advection
    ! qv/theta tendencies out: from microphysics only
    REAL(wp), INTENT(INOUT) :: dqv( kl:ku, il:iu, jl:ju ), dth( kl:ku, il:iu, jl:ju )

    ! aerosol tendencies in:  from parent model forcing i.e. advection
    ! aerosol tendencies out: from microphysics only
    REAL(wp), INTENT(INOUT) :: da1( kl:ku, il:iu, jl:ju ), da2( kl:ku, il:iu, jl:ju ) &
         , da3( kl:ku, il:iu, jl:ju ), da4( kl:ku, il:iu, jl:ju ), da5( kl:ku, il:iu, jl:ju ) &
         , da6( kl:ku, il:iu, jl:ju ), da7( kl:ku, il:iu, jl:ju ), da8( kl:ku, il:iu, jl:ju ) &
         , da9( kl:ku, il:iu, jl:ju ), da10( kl:ku, il:iu, jl:ju ), da11( kl:ku, il:iu, jl:ju ) &
         , da12( kl:ku, il:iu, jl:ju ), da13( kl:ku, il:iu, jl:ju ), da14( kl:ku, il:iu, jl:ju ) &
         , da15( kl:ku, il:iu, jl:ju ), da16( kl:ku, il:iu, jl:ju ), da17( kl:ku, il:iu, jl:ju )

    INTEGER, INTENT(IN), OPTIONAL :: is_in, ie_in ! upper and lower i levels which are to be used
    INTEGER, INTENT(IN), OPTIONAL :: js_in, je_in ! upper and lower j levels

    ! Local variables

    INTEGER :: k, i, j, ixy
    INTEGER :: nxy_all, ixy_outer, ixy_inner, nxy_inner_loop

    REAL(wp) :: precip_l(nxy_inner)
    REAL(wp) :: precip_r(nxy_inner)
    REAL(wp) :: precip_i(nxy_inner)
    REAL(wp) :: precip_s(nxy_inner)
    REAL(wp) :: precip_g(nxy_inner)

    REAL(wp) :: precip_l1d(nz, nxy_inner)
    REAL(wp) :: precip_r1d(nz, nxy_inner)
    REAL(wp) :: precip_i1d(nz, nxy_inner)
    REAL(wp) :: precip_s1d(nz, nxy_inner)
    REAL(wp) :: precip_so1d(nz, nxy_inner)
    REAL(wp) :: precip_g1d(nz, nxy_inner)

    REAL(wp) :: waterpath

    REAL(wp) :: dbz_tot_c(nz), dbz_g_c(nz), dbz_i_c(nz), &
                dbz_s_c(nz),   dbz_l_c(nz), dbz_r_c(nz)

    INTEGER :: kc ! Casim Z-level

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (casim_parent == parent_monc) casim_time = casim_time + dt
    ! (In the KiD model, casim_time is set to variable 'time')

! #if DEF_MODEL==MODEL_KiD
!     !AH - following code limits stops precip processes and sedimentation for no_precip_time. This
!     !     is needed for the KiD-A 2d Sc case
!     if (time <= no_precip_time) then
!        pswitch%l_praut=.false.
!        pswitch%l_pracw=.false.
!        pswitch%l_pracr=.false.
!        pswitch%l_prevp=.false.
!        l_sed=.false.
!     endif
!     ! if (time > no_precip_time .and. l_rain .and. l_sediment) then
!     !    pswitch%l_praut=.true.
!     !    pswitch%l_pracw=.true.
!     !    pswitch%l_pracr=.true.
!     !    pswitch%l_prevp=.true.
!     !    l_sed=.true.
!     ! endif
! #endif



    nxy_all = (je_in-js_in+1)*(ie_in-is_in+1)
    DO ixy_outer=1, int(ceiling(dble(nxy_all)/nxy_inner))
       nxy_inner_loop = min(nxy_inner, &
                        nxy_all-(ixy_outer-1)*nxy_inner)
       DO ixy_inner=1, nxy_inner_loop
          ixy = (ixy_outer-1)*nxy_inner + ixy_inner
          j = modulo(ixy-1,(je_in-js_in+1))+js_in
          i = (ixy-1)/(je_in-js_in+1)+is_in

          precip_l(ixy_inner) = 0.0
          precip_r(ixy_inner) = 0.0
          precip_i(ixy_inner) = 0.0
          precip_s(ixy_inner) = 0.0
          precip_g(ixy_inner) = 0.0

          precip_l1d(:, ixy_inner)  = 0.0
          precip_r1d(:, ixy_inner)  = 0.0
          precip_i1d(:, ixy_inner)  = 0.0
          precip_s1d(:, ixy_inner)  = 0.0
          precip_g1d(:, ixy_inner)  = 0.0
          precip_so1d(:, ixy_inner) = 0.0
          tend(:,:,ixy_inner)=ZERO_REAL_WP
          CALL zero_procs(procs(:,:,ixy_inner))
          aerosol_tend(:,:,ixy_inner)=ZERO_REAL_WP
          IF (l_process) THEN
             CALL zero_procs(aerosol_procs(:,:,ixy_inner))
          END IF
          !set cloud fraction fields
          cffields(:,i_cfl,ixy_inner)=cfliq(k_start:k_end,i,j)
          cffields(:,i_cfr,ixy_inner)=cfrain(k_start:k_end,i,j)
          cffields(:,i_cfi,ixy_inner)=cfice(k_start:k_end,i,j)
          cffields(:,i_cfs,ixy_inner)=cfsnow(k_start:k_end,i,j)
          cffields(:,i_cfg,ixy_inner)=cfgr(k_start:k_end,i,j)
          ! Set the qfields
          qfields(:, i_qv, ixy_inner)=qv(k_start:k_end,i,j)
          qfields(:, i_th, ixy_inner)=theta(k_start:k_end,i,j)
          IF (nq_l > 0) qfields(:,i_ql,ixy_inner)=q1(k_start:k_end,i,j)
          IF (nq_r > 0) qfields(:,i_qr,ixy_inner)=q2(k_start:k_end,i,j)
          IF (nq_l > 1) qfields(:,i_nl,ixy_inner)=q3(k_start:k_end,i,j)
          IF (nq_r > 1) qfields(:,i_nr,ixy_inner)=q4(k_start:k_end,i,j)
          IF (nq_r > 2) qfields(:,i_m3r,ixy_inner)=q5(k_start:k_end,i,j)
          IF (nq_i > 0) qfields(:,i_qi,ixy_inner)=q6(k_start:k_end,i,j)
          IF (nq_s > 0) qfields(:,i_qs,ixy_inner)=q7(k_start:k_end,i,j)
          IF (nq_g > 0) qfields(:,i_qg,ixy_inner)=q8(k_start:k_end,i,j)
          IF (nq_i > 1) qfields(:,i_ni,ixy_inner)=q9(k_start:k_end,i,j)
          IF (nq_s > 1) qfields(:,i_ns,ixy_inner)=q10(k_start:k_end,i,j)
          IF (nq_g > 1) qfields(:,i_ng,ixy_inner)=q11(k_start:k_end,i,j)
          IF (nq_s > 2) qfields(:,i_m3s,ixy_inner)=q12(k_start:k_end,i,j)
          IF (nq_g > 2) qfields(:,i_m3g,ixy_inner)=q13(k_start:k_end,i,j)
          dqfields(:, i_qv, ixy_inner)=dqv(k_start:k_end,i,j)
          dqfields(:, i_th, ixy_inner)=dth(k_start:k_end,i,j)
          IF (nq_l > 0) dqfields(:,i_ql,ixy_inner)=dq1(k_start:k_end,i,j)
          IF (nq_r > 0) dqfields(:,i_qr,ixy_inner)=dq2(k_start:k_end,i,j)
          IF (nq_l > 1) dqfields(:,i_nl,ixy_inner)=dq3(k_start:k_end,i,j)
          IF (nq_r > 1) dqfields(:,i_nr,ixy_inner)=dq4(k_start:k_end,i,j)
          IF (nq_r > 2) dqfields(:,i_m3r,ixy_inner)=dq5(k_start:k_end,i,j)
          IF (nq_i > 0) dqfields(:,i_qi,ixy_inner)=dq6(k_start:k_end,i,j)
          IF (nq_s > 0) dqfields(:,i_qs,ixy_inner)=dq7(k_start:k_end,i,j)
          IF (nq_g > 0) dqfields(:,i_qg,ixy_inner)=dq8(k_start:k_end,i,j)
          IF (nq_i > 1) dqfields(:,i_ni,ixy_inner)=dq9(k_start:k_end,i,j)
          IF (nq_s > 1) dqfields(:,i_ns,ixy_inner)=dq10(k_start:k_end,i,j)
          IF (nq_g > 1) dqfields(:,i_ng,ixy_inner)=dq11(k_start:k_end,i,j)
          IF (nq_s > 2) dqfields(:,i_m3s,ixy_inner)=dq12(k_start:k_end,i,j)
          IF (nq_g > 2) dqfields(:,i_m3g,ixy_inner)=dq13(k_start:k_end,i,j)
          IF (aerosol_option > 0) THEN
             IF (i_am1 >0) aerofields(:, i_am1, ixy_inner)=a1(k_start:k_end,i,j)
             IF (i_an1 >0) aerofields(:, i_an1, ixy_inner)=a2(k_start:k_end,i,j)
             IF (i_am2 >0) aerofields(:, i_am2, ixy_inner)=a3(k_start:k_end,i,j)
             IF (i_an2 >0) aerofields(:, i_an2, ixy_inner)=a4(k_start:k_end,i,j)
             IF (i_am3 >0) aerofields(:, i_am3, ixy_inner)=a5(k_start:k_end,i,j)
             IF (i_an3 >0) aerofields(:, i_an3, ixy_inner)=a6(k_start:k_end,i,j)
             IF (i_am4 >0) aerofields(:, i_am4, ixy_inner)=a7(k_start:k_end,i,j)
             IF (i_am5 >0) aerofields(:, i_am5, ixy_inner)=a8(k_start:k_end,i,j)
             IF (i_am6 >0) aerofields(:, i_am6, ixy_inner)=a9(k_start:k_end,i,j)
             IF (i_an6 >0) aerofields(:, i_an6, ixy_inner)=a10(k_start:k_end,i,j)
             IF (i_am7 >0) aerofields(:, i_am7, ixy_inner)=a11(k_start:k_end,i,j)
             IF (i_am8 >0) aerofields(:, i_am8, ixy_inner)=a12(k_start:k_end,i,j)
             IF (i_am9 >0) aerofields(:, i_am9, ixy_inner)=a13(k_start:k_end,i,j)
             IF (i_am10 >0) aerofields(:, i_am10, ixy_inner)=a14(k_start:k_end,i,j)
             IF (i_an10 >0) aerofields(:, i_an10, ixy_inner)=a15(k_start:k_end,i,j)
             IF (i_an11 >0) aerofields(:, i_an11, ixy_inner)=a16(k_start:k_end,i,j)
             IF (i_an12 >0) aerofields(:, i_an12, ixy_inner)=a17(k_start:k_end,i,j)
             IF (i_ak1 >0) aerofields(:, i_ak1, ixy_inner)=a18(k_start:k_end,i,j)
             IF (i_ak2 >0) aerofields(:, i_ak2, ixy_inner)=a19(k_start:k_end,i,j)
             IF (i_ak3 >0) aerofields(:, i_ak3, ixy_inner)=a20(k_start:k_end,i,j)
             IF (i_am1 >0) daerofields(:, i_am1, ixy_inner)=da1(k_start:k_end,i,j)
             IF (i_an1 >0) daerofields(:, i_an1, ixy_inner)=da2(k_start:k_end,i,j)
             IF (i_am2 >0) daerofields(:, i_am2, ixy_inner)=da3(k_start:k_end,i,j)
             IF (i_an2 >0) daerofields(:, i_an2, ixy_inner)=da4(k_start:k_end,i,j)
             IF (i_am3 >0) daerofields(:, i_am3, ixy_inner)=da5(k_start:k_end,i,j)
             IF (i_an3 >0) daerofields(:, i_an3, ixy_inner)=da6(k_start:k_end,i,j)
             IF (i_am4 >0) daerofields(:, i_am4, ixy_inner)=da7(k_start:k_end,i,j)
             IF (i_am5 >0) daerofields(:, i_am5, ixy_inner)=da8(k_start:k_end,i,j)
             IF (i_am6 >0) daerofields(:, i_am6, ixy_inner)=da9(k_start:k_end,i,j)
             IF (i_an6 >0) daerofields(:, i_an6, ixy_inner)=da10(k_start:k_end,i,j)
             IF (i_am7 >0) daerofields(:, i_am7, ixy_inner)=da11(k_start:k_end,i,j)
             IF (i_am8 >0) daerofields(:, i_am8, ixy_inner)=da12(k_start:k_end,i,j)
             IF (i_am9 >0) daerofields(:, i_am9, ixy_inner)=da13(k_start:k_end,i,j)
             IF (i_am10 >0) daerofields(:, i_am10, ixy_inner)=da14(k_start:k_end,i,j)
             IF (i_an10 >0) daerofields(:, i_an10, ixy_inner)=da15(k_start:k_end,i,j)
             IF (i_an11 >0) daerofields(:, i_an11, ixy_inner)=da16(k_start:k_end,i,j)
             IF (i_an12 >0) daerofields(:, i_an12, ixy_inner)=da17(k_start:k_end,i,j)
          END IF
       END DO ! ixy_inner

       !inner loop pushed into set_passive_fields
       !--------------------------------------------------
       ! set fields which will not be modified
       !--------------------------------------------------
       CALL set_passive_fields(nxy_inner_loop, ixy_outer, is_in, js_in, je_in, &
            dt, rho,    &
            pressure, exner,    &
            dz,                            &
            w, tke, qfields)

       !--------------------------------------------------
       ! Do the business...
       !--------------------------------------------------
       CALL microphysics_common(&
            nxy_inner_loop, &
            ! Inner loop size
            ixy_outer, is_in, js_in, je_in, &
            ! To calculate i, j
            dt, &
            !i , j,
            ! To be calculated so no need anymore
            qfields, cffields, dqfields, tend, procs &
            !, precip(i,j)
            , precip &
            , precip_l, precip_r, precip_i, precip_s, precip_g       &
            , precip_r1d, precip_s1d, precip_so1d, precip_g1d                     &
            , aerophys, aerochem, aeroact                                         &
            , dustphys, dustchem, dustact                                         &
            , aeroice, dustliq                                                    &
            , aerofields, daerofields, aerosol_tend, aerosol_procs                &
            , rhcrit_1d)

       !end do ! ixy_inner
       !--------------------------------------------------
       ! Relate back tendencies
       ! Check indices in mphys_switches that the appropriate
       ! fields are being passed back to mphys_casim
       !--------------------------------------------------

       DO ixy_inner=1, nxy_inner_loop
           ixy = (ixy_outer-1)*nxy_inner + ixy_inner
           j = modulo(ixy-1,(je_in-js_in+1))+js_in
           i = (ixy-1)/(je_in-js_in+1)+is_in

           dqv(k_start:k_end,i,j)=tend(:,i_qv,ixy_inner)
           dth(k_start:k_end,i,j)=tend(:,i_th,ixy_inner)
           dq1(k_start:k_end,i,j)=tend(:,i_ql,ixy_inner)
           dq2(k_start:k_end,i,j)=tend(:,i_qr,ixy_inner)
           IF (cloud_params%l_2m) dq3(k_start:k_end,i,j)=tend(:,i_nl,ixy_inner)
           IF (rain_params%l_2m) dq4(k_start:k_end,i,j)=tend(:,i_nr,ixy_inner)
           IF (rain_params%l_3m) dq5(k_start:k_end,i,j)=tend(:,i_m3r,ixy_inner)

           IF (.NOT. l_warm) THEN
              IF (ice_params%l_1m) dq6(k_start:k_end,i,j)=tend(:,i_qi,ixy_inner)
              IF (snow_params%l_1m) dq7(k_start:k_end,i,j)=tend(:,i_qs,ixy_inner)
              IF (graupel_params%l_1m) dq8(k_start:k_end,i,j)=tend(:,i_qg,ixy_inner)
              IF (ice_params%l_2m) dq9(k_start:k_end,i,j)=tend(:,i_ni,ixy_inner)
              IF (snow_params%l_2m) dq10(k_start:k_end,i,j)=tend(:,i_ns,ixy_inner)
              IF (graupel_params%l_2m) dq11(k_start:k_end,i,j)=tend(:,i_ng,ixy_inner)
              IF (snow_params%l_3m) dq12(k_start:k_end,i,j)=tend(:,i_m3s,ixy_inner)
              IF (graupel_params%l_3m) dq13(k_start:k_end,i,j)=tend(:,i_m3g,ixy_inner)
           END IF

           IF (l_process) THEN
              IF (i_am1 >0) da1(k_start:k_end,i,j)=aerosol_tend(:,i_am1,ixy_inner)
              IF (i_an1 >0) da2(k_start:k_end,i,j)=aerosol_tend(:,i_an1,ixy_inner)
              IF (i_am2 >0) da3(k_start:k_end,i,j)=aerosol_tend(:,i_am2,ixy_inner)
              IF (i_an2 >0) da4(k_start:k_end,i,j)=aerosol_tend(:,i_an2,ixy_inner)
              IF (i_am3 >0) da5(k_start:k_end,i,j)=aerosol_tend(:,i_am3,ixy_inner)
              IF (i_an3 >0) da6(k_start:k_end,i,j)=aerosol_tend(:,i_an3,ixy_inner)
              IF (i_am4 >0) da7(k_start:k_end,i,j)=aerosol_tend(:,i_am4,ixy_inner)
              IF (i_am5 >0) da8(k_start:k_end,i,j)=aerosol_tend(:,i_am5,ixy_inner)
              IF (i_am6 >0) da9(k_start:k_end,i,j)=aerosol_tend(:,i_am6,ixy_inner)
              IF (i_an6 >0) da10(k_start:k_end,i,j)=aerosol_tend(:,i_an6,ixy_inner)
              IF (i_am7 >0) da11(k_start:k_end,i,j)=aerosol_tend(:,i_am7,ixy_inner)
              IF (i_am8 >0) da12(k_start:k_end,i,j)=aerosol_tend(:,i_am8,ixy_inner)
              IF (i_am9 >0) da13(k_start:k_end,i,j)=aerosol_tend(:,i_am9,ixy_inner)
              IF (i_am10 >0) da14(k_start:k_end,i,j)=aerosol_tend(:,i_am10,ixy_inner)
              IF (i_an10 >0) da15(k_start:k_end,i,j)=aerosol_tend(:,i_an10,ixy_inner)
              IF (i_an11 >0) da16(k_start:k_end,i,j)=aerosol_tend(:,i_an11,ixy_inner)
              IF (i_an12 >0) da17(k_start:k_end,i,j)=aerosol_tend(:,i_an12,ixy_inner)
           ELSE
              da1(k_start:k_end,i,j)=0.0
              da2(k_start:k_end,i,j)=0.0
              da3(k_start:k_end,i,j)=0.0
              da4(k_start:k_end,i,j)=0.0
              da5(k_start:k_end,i,j)=0.0
              da6(k_start:k_end,i,j)=0.0
              da7(k_start:k_end,i,j)=0.0
              da9(k_start:k_end,i,j)=0.0
              da10(k_start:k_end,i,j)=0.0
              da11(k_start:k_end,i,j)=0.0
              da12(k_start:k_end,i,j)=0.0
              da13(k_start:k_end,i,j)=0.0
              da14(k_start:k_end,i,j)=0.0
              da15(k_start:k_end,i,j)=0.0
              da16(k_start:k_end,i,j)=0.0
              da17(k_start:k_end,i,j)=0.0
           END IF

           IF ( l_warm ) THEN
              IF ( casdiags % l_surface_cloud ) casdiags % SurfaceCloudR(i,j) = precip_l(ixy_inner)
              IF ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r(ixy_inner)
              IF ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = 0.0
              IF ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = 0.0
              IF ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,k_start:k_end)  = precip_r1d(:,ixy_inner)
              IF ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,k_start:k_end)  = 0.0
              IF ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,k_start:k_end)  = 0.0
              IF ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,k_start:k_end) = 0.0
           ELSE ! l_warm

              IF ( casdiags % l_surface_rain ) casdiags % SurfaceRainR(i,j)  = precip_r(ixy_inner)
              IF ( casdiags % l_surface_snow ) casdiags % SurfaceSnowR(i,j)  = precip_s(ixy_inner)
              IF ( casdiags % l_surface_graup) casdiags % SurfaceGraupR(i,j) = precip_g(ixy_inner)
              IF ( casdiags % l_rainfall_3d ) casdiags % rainfall_3d(i,j,k_start:k_end)  = precip_r1d(:,ixy_inner)
              IF ( casdiags % l_snowfall_3d ) casdiags % snowfall_3d(i,j,k_start:k_end)  = precip_s1d(:,ixy_inner)
              IF ( casdiags % l_snowonly_3d ) casdiags % snowonly_3d(i,j,k_start:k_end)  = precip_so1d(:,ixy_inner)
              IF ( casdiags % l_graupfall_3d) casdiags % graupfall_3d(i,j,k_start:k_end) = precip_g1d(:,ixy_inner)
           END IF ! l_warm

           IF ( casdiags % l_radar ) THEN

              CALL tidy_qin(ixy_inner, qfields(:,:,ixy_inner))  !check this is conserving. If i do this here do we need a tidy_ain?
              CALL casim_reflec(ixy_inner, nz, nq, rho(k_start:k_end,i,j), qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),  &
                               dbz_tot_c, dbz_g_c, dbz_i_c,                      &
                               dbz_s_c,   dbz_l_c, dbz_r_c  )

              casdiags % dbz_tot(i,j, k_start:k_end) = dbz_tot_c(:)
              casdiags % dbz_g(i,j,   k_start:k_end) = dbz_g_c(:)
              casdiags % dbz_s(i,j,   k_start:k_end) = dbz_s_c(:)
              casdiags % dbz_i(i,j,   k_start:k_end) = dbz_i_c(:)
              casdiags % dbz_l(i,j,   k_start:k_end) = dbz_l_c(:)
              casdiags % dbz_r(i,j,   k_start:k_end) = dbz_r_c(:)

           END IF ! casdiags % l_radar

           IF ( casdiags % l_tendency_dg ) THEN
              DO k = k_start, k_end
                 kc = k - k_start + 1
                 casdiags % dth_cond_evap(i,j,k) = procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc) * &
                                                   Lv/cp * rexner(kc,ixy_inner)
                 casdiags % dqv_cond_evap(i,j,k) = -(procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc))
                 casdiags % dth_total(i,j,k) = tend(kc,i_th,ixy_inner)
                 casdiags % dqv_total(i,j,k) = tend(kc,i_qv,ixy_inner)
                 casdiags % dqc(i,j,k) = tend(kc,i_ql,ixy_inner)
                 casdiags % dqr(i,j,k) = tend(kc,i_qr,ixy_inner)

                 IF (.NOT. l_warm) THEN
                    casdiags % dqi(i,j,k) = tend(kc,i_qi,ixy_inner)
                    casdiags % dqs(i,j,k) = tend(kc,i_qs,ixy_inner)
                    casdiags % dqg(i,j,k) = tend(kc,i_qg,ixy_inner)
                 END IF
              END DO
           END IF

           IF ( casdiags % l_lwp ) THEN
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_ql,ixy_inner) )
              END DO
              casdiags % lwp(i,j)=waterpath
           END IF
           IF ( casdiags % l_rwp ) THEN
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qr,ixy_inner) )
              END DO
              casdiags % rwp(i,j)=waterpath
           END IF
           IF ( casdiags % l_iwp ) THEN
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qi,ixy_inner) )
              END DO
              casdiags % iwp(i,j)=waterpath
           END IF
           IF ( casdiags % l_swp ) THEN
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qs,ixy_inner) )
              END DO
              casdiags % swp(i,j)=waterpath
           END IF
           IF ( casdiags % l_gwp ) THEN
              waterpath=0.0
              DO k = k_start, k_end
                 waterpath = waterpath + (rho(k,i,j)*dz(k,i,j) * qfields(k,i_qg,ixy_inner) )
              END DO
              casdiags % gwp(i,j)=waterpath
           END IF
       !end do ! i
    !end do   ! j
        END DO ! ixy_inner
     END DO ! ixy_outer


! #if DEF_MODEL==MODEL_KiD
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny, 'precip', i_dgtime)
!     call save_dg(sum(casdiags % SurfaceRainR(:, :))/nxny*3600.0, 'surface_precip_mmhr', i_dgtime)
! #endif

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE shipway_microphysics

  SUBROUTINE microphysics_common(&
         nxy_inner_loop &
       , ixy_outer, is_in, js_in, je_in &
       , dt &
       !ix, jy,
       , qfields, cffields, dqfields, tend &
       , procs, precip, precip_l, precip_r, precip_i, precip_s, precip_g      &
       , precip_r1d, precip_s1d, precip_so1d, precip_g1d                      &
       , aerophys, aerochem, aeroact                                          &
       , dustphys, dustchem, dustact                                          &
       , aeroice, dustliq                                                     &
       , aerofields, daerofields, aerosol_tend, aerosol_procs                 &
       , rhcrit_1d)


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE casim_parent_mod, ONLY: casim_parent, parent_um

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nxy_inner_loop
    INTEGER, INTENT(IN) :: ixy_outer
    INTEGER, INTENT(IN) :: is_in, js_in, je_in

    REAL(wp), INTENT(IN) :: dt  ! timestep from parent model
    ! integer, intent(in) :: ix, jy

    REAL(wp), INTENT(IN) :: rhcrit_1d(:)
    REAL(wp), INTENT(INOUT) :: qfields(:,:,:), dqfields(:,:,:), tend(:,:,:)
    REAL(wp), INTENT(IN) :: cffields(:,:,:)

    TYPE(process_rate), INTENT(INOUT) :: procs(:,:,:)
    ! real(wp), intent(out) :: precip
    REAL(wp), INTENT(OUT) :: precip(:,:)
    REAL(wp), INTENT(INOUT) :: precip_l(:)
    REAL(wp), INTENT(INOUT) :: precip_r(:)
    REAL(wp), INTENT(INOUT) :: precip_i(:)
    REAL(wp), INTENT(INOUT) :: precip_s(:)
    REAL(wp), INTENT(INOUT) :: precip_g(:)

    REAL(wp), INTENT(INOUT) :: precip_r1d(:,:)
    REAL(wp), INTENT(INOUT) :: precip_s1d(:,:)
    REAL(wp), INTENT(INOUT) :: precip_so1d(:,:)
    REAL(wp), INTENT(INOUT) :: precip_g1d(:,:)

    REAL(wp) :: mindz

    ! Aerosol fields
    TYPE(aerosol_phys), INTENT(INOUT)   :: aerophys(:)
    TYPE(aerosol_chem), INTENT(IN)      :: aerochem(:)
    TYPE(aerosol_active), INTENT(INOUT) :: aeroact(:)
    TYPE(aerosol_phys), INTENT(INOUT)   :: dustphys(:)
    TYPE(aerosol_chem), INTENT(IN)      :: dustchem(:)
    TYPE(aerosol_active), INTENT(INOUT) :: dustact(:)

    TYPE(aerosol_active), INTENT(INOUT) :: aeroice(:)
    TYPE(aerosol_active), INTENT(INOUT) :: dustliq(:)

    REAL(wp), INTENT(INOUT) :: aerofields(:,:,:), daerofields(:,:,:), aerosol_tend(:,:,:)
    TYPE(process_rate), INTENT(INOUT), OPTIONAL :: aerosol_procs(:,:,:)

    REAL(wp) :: step_length

    REAL(wp) :: sed_length, sed_length_cloud, sed_length_rain, sed_length_ice, sed_length_snow &
      , sed_length_graupel

    !--not input anymore--
    INTEGER :: ix, jy, ixy_inner, ixy, i_column

    INTEGER :: n, k, nsed, iq

    LOGICAL :: l_Twarm   ! temperature above freezing

    INTEGER, PARAMETER :: level1 = 1

    ! Local working precipitation rates
    REAL(wp) :: precip_l_w(nxy_inner) ! Liquid cloud precip
    REAL(wp) :: precip_r_w(nxy_inner) ! Rain precip
    REAL(wp) :: precip_i_w(nxy_inner) ! Ice precip
    REAL(wp) :: precip_g_w(nxy_inner) ! Graupel precip
    REAL(wp) :: precip_s_w(nxy_inner) ! Snow precip

    !AH - note that nz is derived in mphys_init and accounts for the lowest level
    !     not equal to 1
    REAL(wp) :: precip1d(nz,nxy_inner) ! local working precip rate
    REAL(wp) :: precip_l_w1d(nz,nxy_inner) ! Liquid cloud precip 1D
    REAL(wp) :: precip_r_w1d(nz,nxy_inner) ! Rain precip 1D
    REAL(wp) :: precip_i_w1d(nz,nxy_inner) ! Ice precip 1D
    REAL(wp) :: precip_g_w1d(nz,nxy_inner) ! Graupel precip 1D
    REAL(wp) :: precip_s_w1d(nz,nxy_inner) ! Snow precip

    INTEGER :: nsubsteps, nsubseds, n_inner

    REAL :: inv_nsubsteps, inv_nsubseds, inv_allsubs
    ! inverse number of substeps for each hydrometor
    REAL :: inv_nsubseds_cloud, inv_nsubseds_ice, inv_nsubseds_rain,           &
      inv_nsubseds_snow, inv_nsubseds_graupel
    REAL :: inv_allsubs_cloud, inv_allsubs_ice, inv_allsubs_rain,              &
      inv_allsubs_snow, inv_allsubs_graupel

    ! number of substeps for each hydrometor
    INTEGER :: nsubseds_cloud, nsubseds_ice, nsubseds_rain,                    &
      nsubseds_snow, nsubseds_graupel

    CHARACTER(len=*), PARAMETER :: RoutineName='MICROPHYSICS_COMMON'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    IF ( l_rp2_casim ) THEN
      snow_params%a_x = snow_a_x_rp
      ice_params%a_x = ice_a_x_rp
    END IF

! Calculations for substeps do not need to stay inside inner loop
! So moving this part after those calculations
!------
!    do ixy_inner=1, nxy_inner_loop

!       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
!       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
!       ix = (ixy-1)/(je_in-js_in+1)+is_in

!       qfields_in(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Initial values of q
!       qfields_mod(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Modified initial values of q (may be modified if bad values sent in)
!-----
    !! AH - Derive the number of microphysical substeps. The default
    !!      max_step_length = 10.0 s in MONC (default) and 120 s in the UM.
    !!      This is set in mphys_switches


    nsubsteps=max(1, ceiling(dt/max_step_length))
    step_length=dt/nsubsteps
    inv_nsubsteps = 1.0 / REAL(nsubsteps)


    !! AH - Derive the maximum number of sedimentation substeps, which
    !!      are performed every microphysics step. max_sed_length = 2.0 s for MONC
    !!      and 120 s for the UM and is set in mphys_switches.
    !!      If step_length (the microphysical timestep) is longer than max_sed_length
    !!      then the sedimentation will substep

    nsubseds=max(1, ceiling(step_length/max_sed_length))
    sed_length=step_length/nsubseds
    inv_nsubseds = 1.0 / REAL(nsubseds)

    inv_allsubs = 1.0 / REAL( nsubseds * nsubsteps)

    nsubseds_cloud = nsubseds
    sed_length_cloud = sed_length
    inv_nsubseds_cloud = inv_nsubseds
    inv_allsubs_cloud = inv_allsubs
    ! rain
    nsubseds_rain = nsubseds
    sed_length_rain = sed_length
    inv_nsubseds_rain = inv_nsubseds
    inv_allsubs_rain = inv_allsubs
    IF (.NOT. l_warm_loc) THEN
       ! ice
       nsubseds_ice = nsubseds
       sed_length_ice = sed_length
       inv_nsubseds_ice = inv_nsubseds
       inv_allsubs_ice = inv_allsubs
       ! snow
       nsubseds_snow = nsubseds
       sed_length_snow = sed_length
       inv_nsubseds_snow = inv_nsubseds
       inv_allsubs_snow = inv_allsubs
       ! graupel
       nsubseds_graupel = nsubseds
       sed_length_graupel = sed_length
       inv_nsubseds_graupel = inv_nsubseds
       inv_allsubs_graupel = inv_allsubs
    END IF


   ! starting inner loop here --
    DO ixy_inner=1, nxy_inner_loop

       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
       ix = (ixy-1)/(je_in-js_in+1)+is_in

       IF (l_subseds_maxv) THEN

          IF (casim_parent == parent_um) THEN
             mindz = 20.0
          ELSE
             !mindz = min_dz ! derived in passive_fields
             mindz = min_dz(ixy_inner)
          END IF
          ! cloud
          !print *, 'terminal vt called'                    ! be careful of mindz, need to have inner dimension too?
          CALL terminal_velocity_CFL(step_length, cloud_params%maxv, nsubseds_cloud, &
               sed_length_cloud, nsubseds, sed_length, mindz)
               inv_nsubseds_cloud = 1.0 / REAL(nsubseds_cloud)
               inv_allsubs_cloud = 1.0 / REAL(nsubseds_cloud * nsubsteps)
               ! rain
          CALL terminal_velocity_CFL(step_length, rain_params%maxv, nsubseds_rain, &
               sed_length_rain, nsubseds, sed_length, mindz)
               inv_nsubseds_rain = 1.0 / REAL(nsubseds_rain)
               inv_allsubs_rain = 1.0 / REAL(nsubseds_rain * nsubsteps)
          IF (.NOT. l_warm_loc) THEN
             ! ice
             CALL terminal_velocity_CFL(step_length, ice_params%maxv, nsubseds_ice, &
                  sed_length_ice, nsubseds, sed_length, mindz)
             inv_nsubseds_ice = 1.0 / REAL(nsubseds_ice)
             inv_allsubs_ice = 1.0 / REAL(nsubseds_ice * nsubsteps)
             ! snow
             CALL terminal_velocity_CFL(step_length, snow_params%maxv, nsubseds_snow, &
                  sed_length_snow, nsubseds, sed_length, mindz)
             inv_nsubseds_snow = 1.0 / REAL(nsubseds_snow)
             inv_allsubs_snow = 1.0 / REAL(nsubseds_snow * nsubsteps)
             ! graupel
             CALL terminal_velocity_CFL(step_length, graupel_params%maxv, nsubseds_graupel, &
                  sed_length_graupel, nsubseds, sed_length, mindz)
             inv_nsubseds_graupel = 1.0 / REAL(nsubseds_graupel)
             inv_allsubs_graupel = 1.0 / REAL(nsubseds_graupel * nsubsteps)
             !print *, nsubseds_cloud,nsubseds_rain,nsubseds_ice,nsubseds_snow,nsubseds_graupel
             !print *, sed_length_cloud,sed_length_rain,sed_length_ice,sed_length_snow,sed_length_graupel
          END IF
       END IF !! V - Do we need to calculate the subseds in every timestep?

       qfields_in(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Initial values of q
       !V -- Don't need to give qields_mod value two times (?)
       !V qfields_mod(:,:,ixy_inner)=qfields(:,:,ixy_inner) ! Modified initial values of q (may be modified if bad values sent in)

       IF (l_tendency_loc) THEN! Parent model uses tendencies
          qfields_mod(:,:,ixy_inner)=qfields_in(:,:,ixy_inner)+dt*dqfields(:,:,ixy_inner)
       ELSE! Parent model uses increments
          qfields_mod(:,:,ixy_inner)=qfields_in(:,:,ixy_inner)+dqfields(:,:,ixy_inner)
       END IF

       IF (.NOT. l_passive) THEN
          CALL tidy_qin(ixy_inner, qfields_mod(:,:,ixy_inner))
       END IF

       !---------------------------------------------------------------
       ! Determine (and possibly limit) size distribution
       !---------------------------------------------------------------
       CALL query_distributions(ixy_inner, cloud_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       CALL query_distributions(ixy_inner, rain_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       IF (.NOT. l_warm_loc) THEN
          CALL query_distributions(ixy_inner, ice_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
          CALL query_distributions(ixy_inner, snow_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
          CALL query_distributions(ixy_inner, graupel_params, qfields_mod(:,:,ixy_inner), cffields(:,:,ixy_inner))
       END IF

       qfields(:,:,ixy_inner)=qfields_mod(:,:,ixy_inner)

       IF (aerosol_option > 0) THEN
          aerofields_in(:,:,ixy_inner)=aerofields(:,:,ixy_inner) ! Initial values of aerosol
          aerofields_mod(:,:,ixy_inner)=aerofields(:,:,ixy_inner) ! Modified initial values  (may be modified if bad values sent in)

          IF (l_tendency_loc) THEN! Parent model uses tendencies
             aerofields_mod(:,:,ixy_inner)=aerofields_in(:,:,ixy_inner)+dt*daerofields(:,:,ixy_inner)
          ELSE! Parent model uses increments
             aerofields_mod(:,:,ixy_inner)=aerofields_in(:,:,ixy_inner)+daerofields(:,:,ixy_inner)
          END IF

          IF (l_process) CALL tidy_ain(qfields_mod(:,:,ixy_inner), aerofields_mod(:,:,ixy_inner))

          aerofields(:,:,ixy_inner)=aerofields_mod(:,:,ixy_inner)
       END IF

    END DO ! ixy_inner


    ! switch n and ixy_inner loop using an "if" condition
    DO n = 1, nsubsteps

       n_inner = 0 !(how many columns has precondition=true under n-loop)

       DO ixy_inner=1, nxy_inner_loop

          ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
          jy = modulo(ixy-1,(je_in-js_in+1))+js_in
          ix = (ixy-1)/(je_in-js_in+1)+is_in

          !do n=1,nsubsteps

          CALL preconditioner(ixy_inner, qfields(:,:,ixy_inner))

          IF ( casdiags % l_mphys_pts ) THEN
             ! Set microphysics points flag based on precondition
             DO k = 1, nz
                casdiags % mphys_pts(ix, jy, k) = precondition(k,ixy_inner)
             END DO
          END IF ! casdiags % l_mphys_pts

         !!------------------------------------------------------
         !! Early exit if we will have nothing to do.
         !! (i.e. no hydrometeors and subsaturated)
         !!------------------------------------------------------
         !if (.not. any(precondition(:,ixy_inner))) exit
         IF (any(precondition(:,ixy_inner))) THEN

            !-------------------------------
            ! Derive aerosol distribution
            ! parameters
            !-------------------------------
            IF (aerosol_option > 0)                                           &
               CALL examine_aerosol(aerofields(:,:,ixy_inner),                &
                    qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,      &
                    dustphys, dustchem, dustact, aeroice, dustliq, icall=1)

            ! In order to get rid of "if precondition", here collect the columns that
            ! have "if precondition" to be "true"
            !--------------------------------------------------------------------------------
            n_inner = n_inner + 1
            columns(n_inner) = ixy_inner

         END IF
       END DO ! ixy_inner


       ! The microphysics will only do in columns with precondtion = true
       IF (n_inner == 0) EXIT ! nothing to do for all columns

       DO i_column = 1, n_inner

          ixy_inner = columns(i_column)

          ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
          jy = modulo(ixy-1,(je_in-js_in+1))+js_in
          ix = (ixy-1)/(je_in-js_in+1)+is_in

          ! Later on the microphysics will only do in columns with precondition=true

          DO k=1,nz

             l_Twarm=TdegK(k,ixy_inner) > 273.15
             l_Tcold(k,ixy_inner)=.NOT. l_Twarm

          END DO
          !
          !=================================
          !
          ! WARM MICROPHYSICAL PROCESSES....
          !
          !=================================
          !
          !-------------------------------
          ! Do the autoconversion to rain
          !-------------------------------
          IF (pswitch%l_praut) THEN
             CALL raut(ixy_inner, step_length, qfields(:,:,ixy_inner),         &
               cffields(:,:,ixy_inner), aerofields(:,:,ixy_inner),             &
               procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
          END IF

          !-------------------------------
          ! Do the rain accreting cloud
          !-------------------------------
          IF (pswitch%l_pracw) THEN
             CALL racw(ixy_inner, step_length, qfields(:,:,ixy_inner),         &
               cffields(:,:,ixy_inner), aerofields(:,:,ixy_inner),             &
               procs(:,:,ixy_inner), rain_params, aerosol_procs(:,:,ixy_inner))
          END IF

          !-------------------------------
          ! Do the rain self-collection
          !-------------------------------
          IF (pswitch%l_pracr) THEN
              CALL racr(ixy_inner, step_length, qfields(:,:,ixy_inner),        &
                procs(:,:,ixy_inner))
          END IF

          !-------------------------------
          ! Do the evaporation of rain
          !-------------------------------
          IF (pswitch%l_prevp) THEN
             !! initialise l_sigevap to false for all levels so that previous
             !! "trues" are not included when rain_precondition is false
             l_sigevap(:,ixy_inner) = .FALSE.
             CALL revp(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),     &
               cffields(:,:,ixy_inner), aerophys, aerochem, aeroact, dustliq,  &
               procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner),             &
               l_sigevap(:,ixy_inner))
          END IF

          !=================================
          !
          ! ICE MICROPHYSICAL PROCESSES....
          !
          !=================================
          ! Start of all ice processes that occur at T < 0C
          !
          IF (.NOT. l_warm_loc) THEN

             !------------------------------------------------------
             ! Condensation/immersion/contact nucleation of cloud ice
             !------------------------------------------------------
             IF (pswitch%l_pinuc) THEN
                CALL inuc(ixy_inner,step_length, nz, l_Tcold(:,ixy_inner),     &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), dustphys, aeroact, dustliq,            &
                  aerosol_procs(:,:,ixy_inner))
             END IF

             !------------------------------------------------------
             ! Autoconverion to snow
             !------------------------------------------------------
             IF (pswitch%l_psaut .AND. .NOT. l_kfsm) THEN
                CALL saut(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),    &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
             END IF
             !------------------------------------------------------
             ! Accretion processes
             !------------------------------------------------------
             ! Ice -> Cloud -> Ice
             IF (pswitch%l_piacw) THEN
                CALL iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),    &
                  ice_params, cloud_params, ice_params, qfields(:,:,ixy_inner),&
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner),               &
                  l_sigevap(:,ixy_inner), aeroact, dustliq,                    &
                  aerosol_procs(:,:,ixy_inner))
             END IF
             ! Snow -> Cloud -> Snow
             IF (l_sg) THEN
                IF (pswitch%l_psacw .AND. .NOT. l_kfsm) THEN
                   CALL iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner), &
                     snow_params, cloud_params, snow_params,                   &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner))
                END IF
                !
                ! Snow -> Ice -> Snow
                IF (pswitch%l_psaci .AND. .NOT. l_kfsm) THEN
                   CALL iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner), &
                     snow_params, ice_params, snow_params,                     &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner))
                END IF
                !
                IF (pswitch%l_praci) THEN
                   ! Rain -> Ice -> Graupel AND Rain -> Ice -> snow, decision made in iacc
                   CALL iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),&
                     rain_params, ice_params, graupel_params,                  &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner), snow_params)
                   ! only one call needed and the decision of graupel or snow is made within iacc
                   ! NOTE: this will break kfsm!!

                END IF
                IF (pswitch%l_psacr .AND. .NOT. l_kfsm) THEN
                   ! Snow -> Rain -> Graupel AND Snow -> Rain -> Snow, decision made in iacc
                   CALL iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),&
                     snow_params, rain_params, graupel_params,                 &
                     qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),          &
                     procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact,    &
                     dustliq, aerosol_procs(:,:,ixy_inner), snow_params)
                END IF
                IF (l_g) THEN
                   ! Graupel -> Cloud -> Graupel
                   IF (pswitch%l_pgacw) THEN
                      CALL iacc(ixy_inner, step_length, nz,                    &
                        l_Tcold(:,ixy_inner), graupel_params, cloud_params,    &
                        graupel_params, qfields(:,:,ixy_inner),                &
                        cffields(:,:,ixy_inner), procs(:,:,ixy_inner),         &
                        l_sigevap(:,ixy_inner), aeroact, dustliq,              &
                        aerosol_procs(:,:,ixy_inner))
                   END IF
                   ! Graupel -> Rain -> Graupel
                   IF (pswitch%l_pgacr) THEN
                      CALL iacc(ixy_inner, step_length, nz,                    &
                        l_Tcold(:,ixy_inner), graupel_params, rain_params,     &
                        graupel_params, qfields(:,:,ixy_inner),                &
                        cffields(:,:,ixy_inner), procs(:,:,ixy_inner),         &
                        l_sigevap(:,ixy_inner), aeroact, dustliq,              &
                        aerosol_procs(:,:,ixy_inner))
                   END IF
                   ! Graupel -> Ice -> Graupel
                   !                   if(pswitch%l_gsaci)call iacc(step_length, k, graupel_params, ice_params, graupel_params, qfields, &
                   !                       procs, aeroact, dustliq, aerosol_procs)
                   ! Graupel -> Snow -> Graupel
                   !                   if(pswitch%l_gsacs)call iacc(step_length, k, graupel_params, snow_params, graupel_params, qfields, &
                   !                       procs, aeroact, dustliq, aerosol_procs)

                ! Graupel -> Ice -> Graupel
                IF (pswitch%l_pgaci) THEN
                   CALL iacc(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner), graupel_params, ice_params, &
                        graupel_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), & 
                        procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact, dustliq, aerosol_procs(:,:,ixy_inner))
                END IF

                ! Graupel -> Snow -> Graupel
                IF (pswitch%l_pgacs) THEN
                   CALL iacc(ixy_inner, step_length, nz, l_Tcold(:,ixy_inner),  graupel_params, snow_params, &
                        graupel_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), &
                        procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroact, dustliq, aerosol_procs(:,:,ixy_inner))
                END IF
                END IF
             END IF

             !------------------------------------------------------
             ! Small snow accreting cloud should be sent to graupel
             ! (Ikawa & Saito 1991)
             !------------------------------------------------------
             IF (.NOT. l_kfsm .AND. l_reisner_graupel_embryo) THEN
                ! Only do this process when Kalli's single moment code is
                ! not in use and l_reisner_graupel_embryo is true; otherwise we ignore it.
                IF (l_g .AND. .NOT. l_onlycollect) THEN
                CALL graupel_embryos(ixy_inner, step_length, nz,               &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner),                &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner))
                END IF ! l_g
             END IF ! not l_kfsm and precondition

             !------------------------------------------------------
             ! Wet deposition/shedding (resulting from graupel
             ! accretion processes)
             ! NB This alters some of the accretion processes, so
             ! must come after their calculation and before they
             ! are used/rescaled elsewhere
             !------------------------------------------------------
             IF (l_g .AND. .NOT. l_onlycollect) THEN
                CALL wetgrowth(ixy_inner, nz, l_Tcold(:,ixy_inner),            &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner))
             END IF

             !------------------------------------------------------
             ! Aggregation (self-collection)
             !------------------------------------------------------
             IF (pswitch%l_psagg .AND. .NOT. l_kfsm) THEN
                CALL ice_aggregation(ixy_inner, step_length, nz,               &
                  l_Tcold(:,ixy_inner), snow_params, qfields(:,:,ixy_inner),   &
                  procs(:,:,ixy_inner))
             END IF

             !------------------------------------------------------
             ! Break up (snow only)
             !------------------------------------------------------
             IF (pswitch%l_psbrk .AND. .NOT. l_kfsm) THEN
                CALL ice_breakup(nz, l_Tcold(:,ixy_inner), snow_params,        &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
             END IF

             !------------------------------------------------------
             ! Ice multiplication (Hallet-mossop)
             !------------------------------------------------------
             IF (pswitch%l_pihal .AND. .NOT. l_kfsm) THEN
                CALL hallet_mossop(ixy_inner, step_length, nz,                 &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner))
             END IF

             !------------------------------------------------------
             ! Homogeneous freezing (rain and cloud)
             !------------------------------------------------------
             IF (pswitch%l_phomr) THEN
                CALL ihom_rain(ixy_inner, step_length, nz,                     &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner),                &
                  l_sigevap(:,ixy_inner), aeroact, dustliq,                    &
                  procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
             END IF

             IF (pswitch%l_phomc) THEN
                CALL ihom_droplets(ixy_inner, step_length, nz,                 &
                  l_Tcold(:,ixy_inner), qfields(:,:,ixy_inner), aeroact,       &
                  dustliq, procs(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner))
             END IF
          
          !------------------------------------------------------
          ! Droplet shattering
          !------------------------------------------------------
          IF (pswitch%l_pidps .AND. .NOT. l_kfsm) THEN
             CALL droplet_shattering(ixy_inner, step_length, nz, cffields(:,:,ixy_inner), &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner))
          END IF

          !------------------------------------------------------
          ! Ice-ice collision (breakup)
          IF (pswitch%l_piics .AND. .NOT. l_kfsm) THEN
             CALL ice_collision(ixy_inner, step_length, nz, cffields(:,:,ixy_inner), &
                  procs(:,:,ixy_inner))
          END IF

             !------------------------------------------------------
             ! Deposition/sublimation of ice/snow/graupel
             !------------------------------------------------------
             IF (pswitch%l_pidep) THEN
                CALL idep(ixy_inner, step_length, nz,  l_Tcold(:,ixy_inner),   &
                  ice_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner), &
                  procs(:,:,ixy_inner), dustact, aeroice,                      &
                  aerosol_procs(:,:,ixy_inner))
             END IF

             IF (pswitch%l_psdep .AND. .NOT. l_kfsm ) THEN
                CALL idep(ixy_inner,step_length, nz, l_Tcold(:,ixy_inner),     &
                  snow_params, qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),&
                  procs(:,:,ixy_inner), dustact, aeroice,                      &
                  aerosol_procs(:,:,ixy_inner))
             END IF

             IF (pswitch%l_pgdep) THEN
                CALL idep(ixy_inner,step_length, nz,  l_Tcold(:,ixy_inner),    &
                  graupel_params, qfields(:,:,ixy_inner),                      &
                  cffields(:,:,ixy_inner), procs(:,:,ixy_inner), dustact,      &
                  aeroice, aerosol_procs(:,:,ixy_inner))
             END IF

             IF (l_harrington .AND. .NOT. l_onlycollect) THEN
                CALL adjust_dep(nz, l_Tcold(:,ixy_inner), procs(:,:,ixy_inner))
             END IF

             !-----------------------------------------------------------
             ! Make sure we don't remove more than saturation allows
             !-----------------------------------------------------------
             IF (l_idep) THEN
                CALL ensure_saturated(ixy_inner, nz, l_Tcold(:,ixy_inner),     &
                  step_length, qfields(:,:,ixy_inner), procs(:,:,ixy_inner),   &
                  (/i_idep, i_sdep, i_gdep/))
             END IF
             !-----------------------------------------------------------
             ! Make sure we don't put back more than saturation allows
             !-----------------------------------------------------------
             IF (l_isub) THEN
                CALL ensure_saturated(ixy_inner, nz, l_Tcold(:,ixy_inner),     &
                  step_length, qfields(:,:,ixy_inner), procs(:,:,ixy_inner),   &
                  (/i_isub, i_ssub, i_gsub/))
             END IF
             ! END all processes at T < 0C
             !
             ! start all ice processes that occur T > 0C, i.e. melting
             !------------------------------------------------------
             ! Melting of ice/snow/graupel
             !------------------------------------------------------
             IF (pswitch%l_psmlt .AND. .NOT. l_kfsm) THEN
                CALL melting(ixy_inner, step_length, nz, snow_params,          &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             END IF
             IF (pswitch%l_pgmlt) THEN
                CALL melting(ixy_inner, step_length, nz, graupel_params,       &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             END IF
             IF (pswitch%l_pimlt) THEN
                CALL melting(ixy_inner, step_length, nz, ice_params,           &
                  qfields(:,:,ixy_inner), cffields(:,:,ixy_inner),             &
                  procs(:,:,ixy_inner), l_sigevap(:,ixy_inner), aeroice,       &
                  dustact, aerosol_procs(:,:,ixy_inner))
             END IF
         END IF ! end if .not. warm_loc
         !-----------------------------------------------------------
         ! Make sure we don't remove more than we have to start with
         !-----------------------------------------------------------
         IF (.NOT. l_warm_loc) THEN
            IF (l_pos1) CALL ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), cloud_params,      &
              (/i_praut, i_pracw, i_iacw, i_sacw, i_gacw, i_homc, i_inuc/),    &
              aeroprocs=aerosol_procs(:,:,ixy_inner),                          &
              iprocs_dependent=(/i_aaut, i_aacw/))

            IF (l_pos2) CALL ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), ice_params,        &
              (/i_raci, i_saci, i_gaci, i_saut, i_isub, i_imlt/),              &
              (/i_ihal, i_idps, i_iics, i_gshd, i_inuc, i_homc, i_iacw, i_idep/))

            IF (l_pos3) CALL ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), rain_params,       &
              (/i_prevp, i_sacr, i_gacr, i_homr/),                             &
              (/i_praut, i_pracw, i_raci, i_gshd, i_smlt, i_gmlt/),            &
              aeroprocs=aerosol_procs(:,:,ixy_inner),                          &
              iprocs_dependent=(/i_arevp/))

            IF (l_pos4) CALL ensure_positive(nz, step_length,                  &
              qfields(:,:,ixy_inner), procs(:,:,ixy_inner), snow_params,       &
              (/i_gacs, i_smlt, i_sacr, i_ssub /),                             &
              (/i_sdep, i_sacw, i_saut, i_saci, i_raci, i_gshd, i_ihal, i_iics/)) 
         ELSE
            IF (pswitch%l_praut .AND. pswitch%l_pracw) THEN
                IF (l_pos5) CALL ensure_positive(nz, step_length,              &
                  qfields(:,:,ixy_inner), procs(:,:,ixy_inner), cloud_params,  &
                  (/i_praut, i_pracw/),                                        &
                  aeroprocs=aerosol_procs(:,:,ixy_inner),                      &
                  iprocs_dependent=(/i_aaut, i_aacw/))
            END IF

            IF (pswitch%l_prevp) THEN
               IF (l_pos6) CALL ensure_positive(nz, step_length,               &
                 qfields(:,:,ixy_inner), procs(:,:,ixy_inner),                 &
                 rain_params, (/i_prevp/), (/i_praut, i_pracw/),               &
                 aerosol_procs(:,:,ixy_inner), (/i_arevp/))
            END IF

         END IF

         !-------------------------------
         ! Collect terms we have so far
         !-------------------------------

         CALL sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),      &
           tend(:,:,ixy_inner), (/i_praut, i_pracw, i_pracr, i_prevp/),        &
           l_thermalexchange=.TRUE., qfields=qfields(:,:,ixy_inner),           &
           l_passive=l_passive)


         IF (.NOT. l_warm_loc) THEN
           IF (.NOT. l_no_pgacs_in_sumprocs) THEN
             CALL sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
                i_gacw, i_gacr, i_gaci, i_gacs, i_ihal, i_iics, i_idps,  i_gshd, i_sbrk,&
                i_saut, i_sagg, i_isub, i_ssub, i_gsub/),        &
                l_thermalexchange=.TRUE., qfields=qfields(:,:,ixy_inner),&
                l_passive=l_passive, i_thirdmoment=2)
           ELSE
             CALL sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                (/i_idep, i_sdep, i_gdep, i_iacw, i_sacr, i_sacw, i_saci, i_raci,&
                i_gacw, i_gacr, i_gaci, i_ihal, i_iics, i_idps, i_gshd, i_sbrk,&
                i_saut, i_sagg, i_isub, i_ssub, i_gsub/),        &
                l_thermalexchange=.TRUE., qfields=qfields(:,:,ixy_inner),&
                l_passive=l_passive, i_thirdmoment=2)
           END IF

            CALL sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),   &
                 tend(:,:,ixy_inner),                                          &
                 (/i_inuc, i_imlt, i_smlt, i_gmlt, i_homr, i_homc/),           &
                 l_thermalexchange=.TRUE., qfields=qfields(:,:,ixy_inner),     &
                 l_passive=l_passive)
         END IF

         CALL update_q(qfields_mod(:,:,ixy_inner), qfields(:,:,ixy_inner),     &
                                          tend(:,:,ixy_inner), l_fixneg=.TRUE.)

         IF (l_process) THEN
             CALL ensure_positive_aerosol(nz, step_length,                     &
                  aerofields(:,:,ixy_inner), aerosol_procs(:,:,ixy_inner),     &
                  (/i_aaut, i_aacw, i_aevp, i_arevp, i_dnuc, i_dsub, i_dssub,  &
                  i_dgsub, i_dimlt, i_dsmlt, i_dgmlt, i_diacw, i_dsacw,        &
                  i_dgacw, i_dsacr, i_dgacr, i_draci, i_dhomr, i_dhomc /) )

            CALL sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),     &
                 aerosol_tend(:,:,ixy_inner),                                  &
                 (/i_aaut, i_aacw, i_aevp, i_arevp/) )

            IF (.NOT. l_warm) THEN
               CALL sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),  &
                    aerosol_tend(:,:,ixy_inner),                               &
                    (/i_dnuc, i_dsub, i_dssub, i_dgsub, i_dimlt, i_dsmlt,      &
                    i_dgmlt, i_diacw, i_dsacw, i_dgacw, i_dsacr, i_dgacr,      &
                    i_draci, i_dhomr, i_dhomc /) )
            END IF

            CALL update_q(aerofields_mod(:,:,ixy_inner),                       &
                 aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),       &
                 l_aerosol=.TRUE.,l_fixneg=.TRUE.)

            !-------------------------------
            ! Re-Derive aerosol distribution
            ! parameters
            !-------------------------------
            CALL examine_aerosol(aerofields(:,:,ixy_inner),                    &
                 qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,          &
                 dustphys, dustchem, dustact, aeroice, dustliq, icall=2)
         END IF

         !-------------------------------
         ! Do the condensation/evaporation
         ! of cloud and activation of new
         ! drops
         !-------------------------------

         IF (casim_parent == parent_um .AND. l_prf_cfrac) THEN
            ! In um and using Paul's cloud fraction so just do update of fields

            !call update_q(qfields_mod, qfields, tend, l_fixneg=.true.)
            CALL query_distributions(ixy_inner, cloud_params,                  &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            CALL query_distributions(ixy_inner, rain_params,                   &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            IF (.NOT. l_warm_loc) THEN
               CALL query_distributions(ixy_inner, ice_params,                 &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, snow_params,                &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, graupel_params,             &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            END IF

         ELSE ! Not in UM or not using Paul's cloud fraction in UM
            IF (pswitch%l_pcond)THEN
               CALL condevp(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),&
                    procs(:,:,ixy_inner), aerophys, aerochem, aeroact,         &
                    dustphys, dustchem, dustliq,                               &
                    aerosol_procs(:,:,ixy_inner), rhcrit_1d)
               !-------------------------------
               ! Collect terms we have so far
               !-------------------------------
               CALL sum_procs(ixy_inner, step_length, nz, procs(:,:,ixy_inner),&
                    tend(:,:,ixy_inner), (/i_cond/),                           &
                    l_thermalexchange=.TRUE., qfields=qfields(:,:,ixy_inner),  &
                    l_passive=l_passive)

               CALL update_q(qfields_mod(:,:,ixy_inner),                       &
                    qfields(:,:,ixy_inner), tend(:,:,ixy_inner),               &
                    l_fixneg=.TRUE.)
            END IF ! pswitch%l_pcond

            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
            CALL query_distributions(ixy_inner, cloud_params,                  &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            CALL query_distributions(ixy_inner, rain_params,                   &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            IF (.NOT. l_warm_loc) THEN
               CALL query_distributions(ixy_inner, ice_params,                 &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, snow_params,                &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, graupel_params,             &
                      qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            END IF

            IF (l_process) THEN
               CALL sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),  &
                    aerosol_tend(:,:,ixy_inner), (/i_aact /) )

               CALL update_q(aerofields_mod(:,:,ixy_inner),                    &
                    aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),    &
                    l_aerosol=.TRUE.)

               !-------------------------------
               ! Re-Derive aerosol distribution
               ! parameters
               !-------------------------------
               CALL examine_aerosol(aerofields(:,:,ixy_inner),                 &
                    qfields(:,:,ixy_inner), aerophys, aerochem, aeroact,       &
                    dustphys, dustchem, dustact, aeroice, dustliq, icall=3)
            END IF

         END IF ! (casim_parent == parent_um .and. l_prf_cfrac)

!!PRF

         !-------------------------------
         ! Do the sedimentation
         !-------------------------------

         precip_l_w(ixy_inner) = 0.0
         precip_r_w(ixy_inner) = 0.0
         precip_i_w(ixy_inner) = 0.0
         precip_g_w(ixy_inner) = 0.0
         precip_s_w(ixy_inner) = 0.0

         DO k = 1, nz
            precip_l_w1d(k, ixy_inner) = 0.0
            precip_r_w1d(k, ixy_inner) = 0.0
            precip_i_w1d(k, ixy_inner) = 0.0
            precip_g_w1d(k, ixy_inner) = 0.0
            precip_s_w1d(k, ixy_inner) = 0.0
         END DO

         IF ( casdiags % l_process_rates ) THEN
            CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=0)
         END IF

         IF (l_sed) THEN
            IF (.NOT. l_subseds_maxv) THEN ! need to add check for 3M code

            !! AH - Following block of code performs sedimentation using the standard (original)
            !!      method, where number of substeps for all hydrometeors are derived using
            !!      max_sed_length and this substep is applied to all hydrometeors
            !!
               DO nsed=1,nsubseds

                  IF (nsed > 1) THEN
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     !call zero_procs_exp(procs)
                     CALL zero_procs(procs(:,:,ixy_inner))
                     IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner))
                     !---------------------------------------------------------------
                     ! Re-Determine (and possibly limit) size distribution
                     !---------------------------------------------------------------
                     CALL query_distributions(ixy_inner, cloud_params,         &
                          qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     CALL query_distributions(ixy_inner, rain_params,          &
                          qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     IF (.NOT. l_warm_loc) THEN
                        CALL query_distributions(ixy_inner, ice_params,        &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        CALL query_distributions(ixy_inner, snow_params,       &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        CALL query_distributions(ixy_inner, graupel_params,    &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                     END IF

                     !-------------------------------
                     ! Re-Derive aerosol distribution
                     ! parameters
                     !-------------------------------
                     IF (l_process) CALL examine_aerosol(                      &
                             aerofields(:,:,ixy_inner), qfields(:,:,ixy_inner),&
                             aerophys, aerochem, aeroact, dustphys, dustchem,  &
                             dustact, aeroice, dustliq, icall=2)
                  END IF ! nsed > 1

                  ! NOTE: if l_gamma_online is true then the original CASIM sedimentation will be used,
                  !       which will calculate the gamma function every timestep. This has to be done
                  !       when 3M or diagnostic shape is used, as shape will change. For single and
                  !       and double moment simulations, it is recommended that l_gamma_online is false.
                  !       This is computationally much more efficient (on CPU and GPU).
                  IF (pswitch%l_psedl) THEN
                     IF (l_gamma_online) THEN
                        CALL sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, cloud_params, procs(:,:,ixy_inner),      &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     ELSE
                        CALL sedr_1M_2M(ixy_inner, sed_length,                 &
                              qfields(:,:,ixy_inner), aeroact, dustliq,        &
                              cloud_params, procs(:,:,ixy_inner),              &
                              aerosol_procs(:,:,ixy_inner),                    &
                              precip1d(:,ixy_inner), l_process)
                     END IF

                     precip_l_w(ixy_inner) = precip_l_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     DO k = 1, nz
                        precip_l_w1d(k,ixy_inner) = precip_l_w1d(k,ixy_inner)  &
                                                        + precip1d(k,ixy_inner)
                     END DO

                     CALL sum_procs(ixy_inner, sed_length, nz,                 &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedl/), qfields=qfields(:,:,ixy_inner))

                  END IF

                  IF (pswitch%l_psedr) THEN
                     IF (l_gamma_online) THEN
                        CALL sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, rain_params, procs(:,:,ixy_inner),       &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     ELSE
                        CALL sedr_1M_2M(ixy_inner, sed_length,                 &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             rain_params, procs(:,:,ixy_inner),                &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     END IF

                     precip_r_w(ixy_inner) = precip_r_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     DO k = 1, nz
                        precip_r_w1d(k,ixy_inner) = precip_r_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     END DO

                     CALL sum_procs(ixy_inner, sed_length, nz,                 &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedr/), qfields=qfields(:,:,ixy_inner))

                  END IF

                  IF (.NOT. l_warm_loc) THEN

                     IF (pswitch%l_psedi) THEN
                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, ice_params,                  &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                ice_params, procs(:,:,ixy_inner),              &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF
                        precip_i_w(ixy_inner) = precip_i_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_i_w1d(k,ixy_inner) =                         &
                              precip_i_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        CALL sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedi/), qfields=qfields(:,:,ixy_inner))

                     END IF

                     IF (pswitch%l_pseds .AND. .NOT. l_kfsm) THEN
                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, snow_params,                 &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                snow_params, procs(:,:,ixy_inner),             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF
                        precip_s_w(ixy_inner) = precip_s_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_s_w1d(k,ixy_inner) =                         &
                              precip_s_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        CALL sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/ i_pseds /), qfields=qfields(:,:,ixy_inner))

                     END IF

                     IF (pswitch%l_psedg) THEN
                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, graupel_params,              &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length,              &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                graupel_params, procs(:,:,ixy_inner),          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF
                        precip_g_w(ixy_inner) = precip_g_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_g_w1d(k,ixy_inner) =                         &
                              precip_g_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        CALL sum_procs(ixy_inner, sed_length, nz,              &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedg/), qfields=qfields(:,:,ixy_inner))

                     END IF

                     !!call sum_procs(sed_length, nz, procs, tend, (/i_psedi, i_pseds, i_psedg/), qfields=qfields)
                  END IF !.not. l_warm_loc

                  CALL update_q(qfields_mod(:,:,ixy_inner),                    &
                       qfields(:,:,ixy_inner), tend(:,:,ixy_inner),            &
                       l_fixneg=.TRUE.)

                  IF (l_process) THEN
                     IF (l_warm) THEN
                        CALL ensure_positive_aerosol(nz, step_length,          &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             (/i_asedr, i_asedl/) )
                        CALL sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedr, i_asedl/))
                        CALL update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                    ELSE ! not l_warm - includes ice procs
                        CALL ensure_positive_aerosol(nz, step_length,          &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             (/i_asedr, i_asedl,i_dsedi, i_dseds, i_dsedg/) )
                        CALL sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner),                      &
                             (/i_asedr, i_asedl,i_dsedi, i_dseds, i_dsedg/) )
                        CALL update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                    END IF ! l_warm
                  END IF ! l_process

                  IF ( casdiags % l_process_rates ) THEN
                     CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end, ncall=1)
                  END IF

               END DO ! nsed

              ELSE ! l_subseds_maxv = true
              !! AH - Following block of code performs sedimentation using a CFL based on the
              !!      prescribed max terminal velocity (mphys_params) for each hydrometeor This
              !!      creates a number of substeps for each  hydrometeor and hence a loop for
              !!      each hydrometeor
              !!
               IF (pswitch%l_psedl) THEN
                  DO nsed=1,nsubseds_cloud
                     IF (nsed > 1) THEN
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        CALL query_distributions(ixy_inner, cloud_params,      &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        IF (l_process) CALL examine_aerosol(                   &
                                            aerofields(:,:,ixy_inner),         &
                                            qfields(:,:,ixy_inner), aerophys,  &
                                            aerochem, aeroact , dustphys,      &
                                            dustchem, dustact, aeroice,        &
                                            dustliq, icall=2)

                     END IF

                     IF (l_gamma_online) THEN
                        CALL sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, cloud_params, procs(:,:,ixy_inner),      &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     ELSE
                        CALL sedr_1M_2M(ixy_inner, sed_length_cloud,           &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             cloud_params, procs(:,:,ixy_inner),               &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     END IF

                     precip_l_w(ixy_inner) = precip_l_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     DO k = 1, nz
                        precip_l_w1d(k,ixy_inner) = precip_l_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     END DO

                     IF ( casdiags % l_process_rates ) THEN
                        CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                     END IF

                     CALL sum_procs(ixy_inner, sed_length_cloud, nz,           &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedl/), qfields=qfields(:,:,ixy_inner))

                     CALL update_q(qfields_mod(:,:,ixy_inner),                 &
                          qfields(:,:,ixy_inner), tend(:,:,ixy_inner),         &
                          l_fixneg=.TRUE.)

                     IF (l_process) THEN
                        CALL sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedl/) )
                        CALL update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                     END IF
                     ! !-------------------------------
                     ! ! Reset process rates if they
                     ! ! are to be re-used
                     ! !-------------------------------
                     CALL zero_procs(procs(:,:,ixy_inner), (/i_psedl/))
                     IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner), (/i_asedl/))
                  END DO ! k
               END IF ! pswitch%l_psedl

               IF (pswitch%l_psedr) THEN
                  DO nsed=1,nsubseds_rain
                     IF (nsed > 1) THEN
                        !---------------------------------------------------------------
                        ! Re-Determine (and possibly limit) size distribution
                        !---------------------------------------------------------------
                        CALL query_distributions(ixy_inner, rain_params,       &
                             qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                        !-------------------------------
                        ! Re-Derive aerosol distribution
                        ! parameters
                        !-------------------------------
                        IF (l_process) CALL examine_aerosol(                   &
                             aerofields(:,:,ixy_inner), qfields(:,:,ixy_inner),&
                             aerophys, aerochem, aeroact, dustphys, dustchem,  &
                             dustact, aeroice, dustliq, icall=2)
                     END IF
                     IF (l_gamma_online) THEN
                        CALL sedr(ixy_inner, qfields(:,:,ixy_inner), aeroact,  &
                             dustliq, rain_params, procs(:,:,ixy_inner),       &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     ELSE
                        CALL sedr_1M_2M(ixy_inner, sed_length_rain,            &
                             qfields(:,:,ixy_inner), aeroact, dustliq,         &
                             rain_params, procs(:,:,ixy_inner),                &
                             aerosol_procs(:,:,ixy_inner),                     &
                             precip1d(:,ixy_inner), l_process)
                     END IF

                     precip_r_w(ixy_inner) = precip_r_w(ixy_inner) +           &
                                                     precip1d(level1,ixy_inner)

                     DO k = 1, nz
                        precip_r_w1d(k,ixy_inner) = precip_r_w1d(k,ixy_inner) +&
                                                          precip1d(k,ixy_inner)
                     END DO

                     IF ( casdiags % l_process_rates ) THEN
                        CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                     END IF

                     CALL sum_procs(ixy_inner, sed_length_rain, nz,            &
                          procs(:,:,ixy_inner), tend(:,:,ixy_inner),           &
                          (/i_psedr/), qfields=qfields(:,:,ixy_inner))
                     CALL update_q(qfields_mod(:,:,ixy_inner),                 &
                          qfields(:,:,ixy_inner), tend(:,:,ixy_inner),         &
                          l_fixneg=.TRUE.)
                     IF (l_process) THEN
                        CALL sum_aprocs(sed_length, nz,                        &
                             aerosol_procs(:,:,ixy_inner),                     &
                             aerosol_tend(:,:,ixy_inner), (/i_asedr, i_asedl/))
                        CALL update_q(aerofields_mod(:,:,ixy_inner),           &
                             aerofields(:,:,ixy_inner),                        &
                             aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                     END IF
                     !-------------------------------
                     ! Reset process rates if they
                     ! are to be re-used
                     !-------------------------------
                     CALL zero_procs(procs(:,:,ixy_inner), (/i_psedr/))
                     IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner), (/i_asedr/))

                  END DO ! nsed
               END IF ! pswitch%l_psedr

               IF (.NOT. l_warm_loc) THEN

                  IF (pswitch%l_psedi) THEN
                     DO nsed=1,nsubseds_ice
                        IF (nsed > 1) THEN
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           CALL query_distributions(ixy_inner, ice_params,     &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           !-------------------------------
                           ! Re-Derive aerosol distribution
                           ! parameters
                           !-------------------------------
                           IF (l_process) CALL examine_aerosol(                &
                                        aerofields(:,:,ixy_inner),             &
                                        qfields(:,:,ixy_inner), aerophys,      &
                                        aerochem, aeroact , dustphys, dustchem,&
                                        dustact, aeroice, dustliq, icall=2)
                        END IF
                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, ice_params,                  &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length_ice,          &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                ice_params, procs(:,:,ixy_inner),              &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF

                        precip_i_w(ixy_inner) = precip_i_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_i_w1d(k,ixy_inner) =                         &
                              precip_i_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        IF ( casdiags % l_process_rates ) THEN
                           CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        END IF

                        CALL sum_procs(ixy_inner, sed_length_ice, nz,          &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedi/), qfields=qfields(:,:,ixy_inner))
                        CALL update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.TRUE.)
                        IF (l_process) THEN
                           CALL sum_aprocs(sed_length, nz,                     &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dsedi/))
                           CALL update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                        END IF
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        CALL zero_procs(procs(:,:,ixy_inner), (/i_psedi/))
                        IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dsedi/))
                     END DO ! nsed
                  END IF !pswitch%l_psedi

                  IF (pswitch%l_pseds) THEN
                     DO nsed=1,nsubseds_snow
                        IF (nsed > 1) THEN
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           CALL query_distributions(ixy_inner, snow_params,    &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           ! !-------------------------------
                           ! ! Re-Derive aerosol distribution
                           ! ! parameters
                           ! !-------------------------------
                           IF (l_process) CALL examine_aerosol(                &
                                    aerofields(:,:,ixy_inner),                 &
                                    qfields(:,:,ixy_inner), aerophys, aerochem,&
                                    aeroact, dustphys, dustchem, dustact,      &
                                    aeroice, dustliq, icall=2)
                        END IF

                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, snow_params,                 &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length_snow,         &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                snow_params, procs(:,:,ixy_inner),             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF

                        precip_s_w(ixy_inner) = precip_s_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_s_w1d(k,ixy_inner) =                         &
                              precip_s_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        IF ( casdiags % l_process_rates ) THEN
                           CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        END IF

                        CALL sum_procs(ixy_inner, sed_length_snow, nz,         &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_pseds/), qfields=qfields(:,:,ixy_inner))
                        CALL update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.TRUE.)
                        IF (l_process) THEN
                           CALL sum_aprocs(sed_length, nz,                     &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dseds/) )
                           CALL update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                        END IF
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        CALL zero_procs(procs(:,:,ixy_inner), (/i_pseds/))
                        IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dseds/))
                     END DO ! nsed
                  END IF ! pswitch%l_pseds

                  IF (pswitch%l_psedg) THEN
                     DO nsed=1,nsubseds_graupel
                        IF (nsed > 1) THEN
                           !---------------------------------------------------------------
                           ! Re-Determine (and possibly limit) size distribution
                           !---------------------------------------------------------------
                           CALL query_distributions(ixy_inner, graupel_params, &
                               qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
                           !-------------------------------
                           ! Re-Derive aerosol distribution
                           ! parameters
                           !-------------------------------
                           IF (l_process) CALL examine_aerosol(                &
                                   aerofields(:,:,ixy_inner),                  &
                                   qfields(:,:,ixy_inner), aerophys, aerochem, &
                                   aeroact, dustphys, dustchem, dustact,       &
                                   aeroice, dustliq, icall=2)
                        END IF
                        IF (l_gamma_online) THEN
                           CALL sedr(ixy_inner, qfields(:,:,ixy_inner),        &
                                aeroice, dustact, graupel_params,              &
                                procs(:,:,ixy_inner),                          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        ELSE
                           CALL sedr_1M_2M(ixy_inner, sed_length_graupel,      &
                                qfields(:,:,ixy_inner), aeroice, dustact,      &
                                graupel_params, procs(:,:,ixy_inner),          &
                                aerosol_procs(:,:,ixy_inner),                  &
                                precip1d(:,ixy_inner), l_process)
                        END IF

                        precip_g_w(ixy_inner) = precip_g_w(ixy_inner) +        &
                                                     precip1d(level1,ixy_inner)

                        DO k = 1, nz
                           precip_g_w1d(k,ixy_inner) =                         &
                              precip_g_w1d(k,ixy_inner) + precip1d(k,ixy_inner)
                        END DO

                        IF ( casdiags % l_process_rates ) THEN
                           CALL gather_process_diagnostics(ixy_inner, ix, jy, k_start, k_end,ncall=1)
                        END IF

                        CALL sum_procs(ixy_inner, sed_length_graupel, nz,      &
                             procs(:,:,ixy_inner), tend(:,:,ixy_inner),        &
                             (/i_psedg/), qfields=qfields(:,:,ixy_inner))
                        CALL update_q(qfields_mod(:,:,ixy_inner),              &
                             qfields(:,:,ixy_inner), tend(:,:,ixy_inner),      &
                             l_fixneg=.TRUE.)
                        IF (l_process) THEN
                           CALL sum_aprocs(sed_length_graupel, nz,             &
                                aerosol_procs(:,:,ixy_inner),                  &
                                aerosol_tend(:,:,ixy_inner), (/i_dsedg/))

                           CALL update_q(aerofields_mod(:,:,ixy_inner),        &
                                aerofields(:,:,ixy_inner),                     &
                                aerosol_tend(:,:,ixy_inner), l_aerosol=.TRUE.)
                        END IF
                        !-------------------------------
                        ! Reset process rates if they
                        ! are to be re-used
                        !-------------------------------
                        CALL zero_procs(procs(:,:,ixy_inner), (/i_psedg/))
                        IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner), (/i_dsedg/))
                     END DO ! nsed
                  END IF ! pswitch%l_psedg
               END IF  ! l_warm_loc
              END IF !.not. l_subseds_maxv
         END IF ! l_sed

         precip_l(ixy_inner) = precip_l(ixy_inner) + precip_l_w(ixy_inner)
         ! For diagnostic purposes, set precip_r, precip_s and precip to pass out
         ! For the UM, rainfall rate is assumed as sum of all liquid components
         ! (so includes sedimentation of rain and liquid cloud)
         precip_r(ixy_inner) = precip_r(ixy_inner) + precip_r_w(ixy_inner)

         ! For the UM, snowfall rate is assumed to be a sum of all solid components
         ! (so includes ice, snow and graupel)
         precip_s(ixy_inner) = precip_s(ixy_inner) + precip_s_w(ixy_inner)
         precip_i(ixy_inner) = precip_i(ixy_inner) + precip_i_w(ixy_inner)
         ! For the UM, graupel rate is just itself
         precip_g(ixy_inner) = precip_g(ixy_inner) + precip_g_w(ixy_inner)

         DO k = 1, nz
            precip_r1d(k,ixy_inner)  = precip_r1d(k,ixy_inner)  +              &
                          precip_l_w1d(k,ixy_inner) + precip_r_w1d(k,ixy_inner)
            precip_g1d(k,ixy_inner)  = precip_g1d(k,ixy_inner)  +              &
                                                      precip_g_w1d(k,ixy_inner)
            precip_s1d(k,ixy_inner)  = precip_s1d(k,ixy_inner)  +              &
                       precip_s_w1d(k,ixy_inner) + precip_i_w1d(k,ixy_inner) + &
                       precip_g_w1d(k,ixy_inner)
            precip_so1d(k,ixy_inner) = precip_so1d(k,ixy_inner) +              &
                          precip_s_w1d(k,ixy_inner) + precip_i_w1d(k,ixy_inner)
         END DO

         IF (nsubsteps>1)THEN
            !-------------------------------
            ! Reset process rates if they
            ! are to be re-used
            !-------------------------------
            !call zero_procs_exp(procs)
            CALL zero_procs(procs(:,:,ixy_inner))
            IF (l_process) CALL zero_procs(aerosol_procs(:,:,ixy_inner))
            !---------------------------------------------------------------
            ! Re-Determine (and possibly limit) size distribution
            !---------------------------------------------------------------
            CALL query_distributions(ixy_inner, cloud_params,                  &
                 qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            CALL query_distributions(ixy_inner, rain_params,                   &
                 qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            IF (.NOT. l_warm_loc) THEN
               CALL query_distributions(ixy_inner, ice_params,                 &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, snow_params,                &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
               CALL query_distributions(ixy_inner, graupel_params,             &
                    qfields(:,:,ixy_inner), cffields(:,:,ixy_inner))
            END IF
         END IF
       !end if ! precondition
    END DO ! i_column
    END DO ! nsubsteps

    DO ixy_inner=1, nxy_inner_loop

       ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
       jy = modulo(ixy-1,(je_in-js_in+1))+js_in
       ix = (ixy-1)/(je_in-js_in+1)+is_in

       ! We want the mean precipitation over the parent timestep - so divide
       ! by the total number of substeps - multiply by inv_allsubs is quicker
       precip_l(ixy_inner) = precip_l(ixy_inner) * inv_allsubs_cloud
       precip_r(ixy_inner) = precip_r(ixy_inner) * inv_allsubs_rain
       precip_i(ixy_inner) = precip_i(ixy_inner) * inv_allsubs_ice
       precip_s(ixy_inner) = precip_s(ixy_inner) * inv_allsubs_snow
       precip_g(ixy_inner) = precip_g(ixy_inner) * inv_allsubs_graupel

       ! UM precip rates are
       precip_r(ixy_inner) = precip_l(ixy_inner) + precip_r(ixy_inner)
       precip_s(ixy_inner) = precip_i(ixy_inner) + precip_s(ixy_inner) + precip_g(ixy_inner)

       DO k = 1, nz
          precip_r1d(k,ixy_inner)  = precip_r1d(k,ixy_inner)  * inv_allsubs
          precip_s1d(k,ixy_inner)  = precip_s1d(k,ixy_inner)  * inv_allsubs
          precip_so1d(k,ixy_inner) = precip_so1d(k,ixy_inner) * inv_allsubs
          precip_g1d(k,ixy_inner)  = precip_g1d(k,ixy_inner)  * inv_allsubs
       END DO ! k

       ! Precip is a sum of everything, so just add rain and snow together which
       ! has all components added.
       ! Do not add precip_g, otherwise graupel contributions will be double-counted
       precip(ix,jy) = precip_r(ixy_inner)   + precip_s(ixy_inner)

       !--------------------------------------------------
       ! Tidy up any small/negative numbers
       ! we may have generated.
       !--------------------------------------------------
       IF (pswitch%l_tidy2) THEN

          CALL qtidy(ixy_inner, step_length, nz, qfields(:,:,ixy_inner),       &
               procs(:,:,ixy_inner), aerofields(:,:,ixy_inner), aeroact,       &
               dustact, aeroice, dustliq , aerosol_procs(:,:,ixy_inner),       &
               i_tidy2, i_atidy2, l_negonly=l_tidy_negonly)

          CALL sum_procs(ixy_inner, step_length, nz,                           &
               procs(:,:,ixy_inner), tend(:,:,ixy_inner), (/i_tidy2/),         &
               qfields=qfields(:,:,ixy_inner), l_passive=l_passive)

          CALL update_q(qfields_mod(:,:,ixy_inner), qfields(:,:,ixy_inner), tend(:,:,ixy_inner))

          IF (l_process) THEN
             CALL sum_aprocs(step_length, nz, aerosol_procs(:,:,ixy_inner),    &
                  aerosol_tend(:,:,ixy_inner), (/i_atidy2/) )
             CALL update_q(aerofields_mod(:,:,ixy_inner),                      &
                  aerofields(:,:,ixy_inner), aerosol_tend(:,:,ixy_inner),      &
                  l_aerosol=.TRUE.)
          END IF
       END IF

       !
       ! Add on initial adjustments that may have been made
       !

       IF (l_tendency_loc) THEN! Convert back from cumulative value to tendency
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)+qfields_mod(:,:,ixy_inner)-  &
                           qfields_in(:,:,ixy_inner)-dqfields(:,:,ixy_inner)*dt
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)/dt
       ELSE
          tend(:,:,ixy_inner)=tend(:,:,ixy_inner)+qfields_mod(:,:,ixy_inner)-  &
                              qfields_in(:,:,ixy_inner)-dqfields(:,:,ixy_inner)
          ! prevent negative values
          DO iq=i_hstart,ntotalq
             DO k=1,nz
                tend(k,iq,ixy_inner)=max(tend(k,iq,ixy_inner),                 &
                        -(qfields_in(k,iq,ixy_inner)-dqfields(k,iq,ixy_inner)))
            END DO
          END DO
       END IF

       IF (aerosol_option > 0) THEN
          ! processing
          IF (l_process) THEN
             IF (l_tendency_loc) THEN! Convert back from cumulative value to tendency
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)+       &
                    aerofields_mod(:,:,ixy_inner)-aerofields_in(:,:,ixy_inner)-&
                    daerofields(:,:,ixy_inner)*dt
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)/dt
             ELSE
                aerosol_tend(:,:,ixy_inner)=aerosol_tend(:,:,ixy_inner)+       &
                    aerofields_mod(:,:,ixy_inner)-aerofields_in(:,:,ixy_inner)-&
                    daerofields(:,:,ixy_inner)
                ! prevent negative values
                DO iq=1,ntotala
                   DO k=1,nz
                      aerosol_tend(k,iq,ixy_inner)=                            &
                           max(aerosol_tend(k,iq,ixy_inner),                   &
                           -(aerofields_in(k,iq,ixy_inner)-daerofields(k,iq,ixy_inner)))
                   END DO
                END DO
             END IF
          END IF
       END IF

    END DO ! ixy_inner

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE microphysics_common

  SUBROUTINE update_q(qfields_in, qfields, tend, l_aerosol, l_fixneg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: qfields_in(:,:)
    REAL(wp), INTENT(INOUT) :: qfields(:,:)
    REAL(wp), INTENT(IN) :: tend(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: l_aerosol ! flag to indicate updating of aerosol
    LOGICAL, INTENT(IN), OPTIONAL :: l_fixneg  ! Flag to use cludge to bypass negative/zero numbers
    INTEGER :: k, iqx
    LOGICAL :: l_fix

    CHARACTER(len=*), PARAMETER :: RoutineName='UPDATE_Q'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_fix=.FALSE.
    IF (present(l_fixneg)) l_fix=l_fixneg

    DO iqx=1, ubound(tend,2)
       DO k=lbound(tend,1), ubound(tend,1)
          qfields(k,iqx)=qfields_in(k,iqx)+tend(k,iqx)
       END DO
    END DO

     IF (.NOT. present(l_aerosol) .AND. l_fix) THEN
      !quick lem fixes  - this code should never be used ?
        DO iqx=1, ubound(tend,2)
          DO k=lbound(tend,1), ubound(tend,1)
             IF (iqx==i_ni .AND. qfields(k,iqx)<=0.0) THEN
                qfields(k,iqx)=0.0
                qfields(k,i_qi)=0.0
             END IF
             IF (iqx==i_nr .AND. qfields(k,iqx)<=0.0) THEN
                qfields(k,iqx)=0.0
                qfields(k,i_qr)=0.0
                IF (i_m3r/=0) qfields(k,i_m3r)=0.0
             END IF
             IF (iqx==i_nl .AND. qfields(k,iqx)<=0.0) THEN
                qfields(k,iqx)=0.0
                qfields(k,i_ql)=0.0
             END IF
             IF (iqx==i_ns .AND. qfields(k,iqx)<=0.0) THEN
                qfields(k,iqx)=0.0
                qfields(k,i_qs)=0.0
                IF (i_m3s/=0) qfields(k,i_m3s)=0.0
             END IF
             IF (iqx==i_ng .AND. qfields(k,iqx)<=0.0) THEN
                qfields(k,iqx)=0.0
                qfields(k,i_qg)=0.0
                IF (i_m3g/=0) qfields(k,i_m3g)=0.0
             END IF
          END DO
       END DO
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE update_q

  SUBROUTINE gather_process_diagnostics(ixy_inner, i, j, k_start, k_end,ncall)

    ! Gathers all process rate diagnostics if in use and outputs them to the
    ! CASIM generic diagnostic fields, ready for use in any model.

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Indices of this particular grid square
    INTEGER, INTENT(IN) :: ixy_inner, i, j,ncall
    INTEGER, INTENT(IN) :: k_start, k_end ! Start/end points of grid

    ! Local variables

    CHARACTER(len=*), PARAMETER :: RoutineName='GATHER_PROCESS_DIAGNOSTICS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    INTEGER :: kc ! Casim Z-level
    INTEGER :: k  ! Loop counter in z-direction

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Based on code from RGS: fill in process rates:

    IF (ncall==0) THEN
    IF (casdiags % l_phomc) THEN
      IF (pswitch%l_phomc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % phomc(i,j,k) = procs(ice_params%i_1m,i_homc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % phomc(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_phomc

    IF (casdiags % l_nhomc) THEN
      IF ((pswitch%l_phomc) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nhomc(i,j,k) = procs(ice_params%i_2m,i_homc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomc(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_phomc

    IF (casdiags % l_pinuc) THEN
      IF (pswitch%l_pinuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pinuc(i,j,k) = procs(ice_params%i_1m,i_inuc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pinuc(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_ninuc) THEN
      IF ((pswitch%l_pinuc) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ninuc(i,j,k) = procs(ice_params%i_2m,i_inuc%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ninuc (i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! casdiags % l_pinuc

    IF (casdiags % l_pidep) THEN
      IF (pswitch%l_pidep) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pidep(i,j,k) = procs(ice_params%i_1m,i_idep%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pidep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psdep) THEN
      IF (pswitch%l_psdep) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psdep(i,j,k) = procs(snow_params%i_1m,i_sdep%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psdep(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_piacw) THEN
      IF (pswitch%l_piacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % piacw(i,j,k) = procs(ice_params%i_1m,i_iacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % piacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niacw) THEN
      IF ((pswitch%l_piacw) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niacw(i,j,k) = procs(cloud_params%i_2m,i_iacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF  

    IF (casdiags % l_psacw) THEN
      IF (pswitch%l_psacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psacw(i,j,k) = procs(snow_params%i_1m,i_sacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsacw) THEN
      IF ((pswitch%l_psacw) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsacw(i,j,k) = procs(cloud_params%i_2m,i_sacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psacr) THEN
      IF (pswitch%l_psacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psacr(i,j,k) = procs(snow_params%i_1m,i_sacr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psacr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsacr) THEN
      IF ((pswitch%l_psacr) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsacr(i,j,k) = procs(rain_params%i_2m,i_sacr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsacr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pisub) THEN
      IF (pswitch%l_pisub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pisub(i,j,k) = -1.0 * procs(ice_params%i_1m,i_isub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pisub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nisub) THEN
      IF ((pswitch%l_pisub) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nisub(i,j,k) = -1.0 * procs(ice_params%i_2m,i_isub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nisub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pssub) THEN
      IF (pswitch%l_pssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pssub(i,j,k) = -1.0 * procs(snow_params%i_1m,i_ssub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pssub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nssub) THEN
      IF ((pswitch%l_pssub) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nssub(i,j,k) = -1.0 * procs(snow_params%i_2m,i_ssub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nssub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pimlt) THEN
      IF (pswitch%l_pimlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pimlt(i,j,k) = procs(rain_params%i_1m,i_imlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pimlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nimlt) THEN
      IF ((pswitch%l_pimlt) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nimlt(i,j,k) = procs(rain_params%i_2m,i_imlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nimlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psmlt) THEN
      IF (pswitch%l_psmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psmlt(i,j,k) = procs(rain_params%i_1m,i_smlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsmlt) THEN
      IF ((pswitch%l_psmlt) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsmlt(i,j,k) = procs(rain_params%i_2m,i_smlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaut) THEN
      IF (pswitch%l_psaut) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psaut(i,j,k) = procs(snow_params%i_1m,i_saut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsaut) THEN
      IF ((pswitch%l_psaut) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsaut(i,j,k) = procs(snow_params%i_2m,i_saut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psaci) THEN
      IF (pswitch%l_psaci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psaci(i,j,k) = procs(snow_params%i_1m,i_saci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsaci) THEN
      IF ((pswitch%l_psaci) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsaci(i,j,k) = procs(ice_params%i_2m,i_saci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praut) THEN
      IF (pswitch%l_praut) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praut(i,j,k) = procs(rain_params%i_1m,i_praut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraut) THEN
      IF ((pswitch%l_praut) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraut(i,j,k) = procs(rain_params%i_2m,i_praut%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pracw) THEN
      IF (pswitch%l_pracw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pracw(i,j,k) = procs(rain_params%i_1m,i_pracw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pracw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nracw) THEN
      IF ((pswitch%l_pracw) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nracw(i,j,k) = procs(cloud_params%i_2m,i_pracw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nracw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nracr) THEN
      IF ((pswitch%l_pracr) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nracr(i,j,k) = procs(rain_params%i_2m,i_pracr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nracr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_prevp) THEN
      IF (pswitch%l_prevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % prevp(i,j,k) = -1.0 * procs(rain_params%i_1m,i_prevp%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % prevp(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nrevp) THEN
      IF ((pswitch%l_prevp) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nrevp(i,j,k) = -1.0 * procs(rain_params%i_2m,i_prevp%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nrevp(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

! #if DEF_MODEL==MODEL_KiD
!      IF (casdiags % l_praut) THEN
!        IF (pswitch%l_praut) THEN
!           DO k = k_start, k_end
!              kc = k - k_start + 1
!              call save_dg(k, casdiags % praut(i,j,k) , 'praut', i_dgtime)
!           END DO
!        ENDIF
!     ENDIF

!     IF (casdiags % l_pracw) THEN
!         IF (pswitch%l_pracw) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, casdiags % pracw(i,j,k) , 'pracw', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_prevp) THEN
!         IF (pswitch%l_prevp) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, casdiags % prevp(i,j,k) , 'prevp', i_dgtime)
!             END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedr) THEN
!         IF (pswitch%l_psedr) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc) , 'psedr', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_psedl) THEN
!         IF (pswitch%l_psedl) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc) , 'psedl', i_dgtime)
!            END DO
!         END IF
!      END IF

!      IF (casdiags % l_pracr) THEN
!         IF (pswitch%l_pracr) THEN
!            DO k = k_start, k_end
!               kc = k - k_start + 1
!               call save_dg(k, procs(rain_params%i_1m,i_pracr%id,ixy_inner)%column_data(kc), 'pracr', i_dgtime)
!            END DO
!         END IF
!      END IF
! #endif
    IF (casdiags % l_pgacw) THEN
      IF (pswitch%l_pgacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgacw(i,j,k) = procs(graupel_params%i_1m, i_gacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngacw) THEN
      IF ((pswitch%l_pgacw) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngacw(i,j,k) = procs(cloud_params%i_2m, i_gacw%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgacs) THEN
      IF (pswitch%l_pgacs) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgacs(i,j,k) = procs(graupel_params%i_1m, i_gacs%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgacs(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngacs) THEN
      IF ((pswitch%l_pgacs) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngacs(i,j,k) = procs(snow_params%i_2m, i_gacs%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngacs(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgmlt) THEN
      IF (pswitch%l_pgmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgmlt(i,j,k) = procs(rain_params%i_1m, i_gmlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngmlt) THEN
      IF ((pswitch%l_pgmlt) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngmlt(i,j,k) = procs(rain_params%i_2m, i_gmlt%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngmlt(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgsub) THEN
      IF (pswitch%l_pgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgsub(i,j,k) = -1.0 * procs(graupel_params%i_1m,i_gsub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgsub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngsub) THEN
      IF ((pswitch%l_pgsub) .AND. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngsub(i,j,k) = -1.0 * procs(graupel_params%i_2m, i_gsub%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngsub(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedi(i,j,k) = procs(ice_params%i_1m, i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedi(i,j,k) = procs(ice_params%i_2m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedi(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pseds(i,j,k) = procs(snow_params%i_1m, i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nseds(i,j,k) = procs(snow_params%i_2m, i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nseds(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedr(i,j,k) = procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedr) THEN
      IF ((pswitch%l_psedr) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedr(i,j,k) = procs(rain_params%i_2m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedg(i,j,k) = procs(graupel_params%i_1m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .AND. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedg(i,j,k) = procs(graupel_params%i_2m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedg(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedl(i,j,k) = procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % psedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nsedl) THEN
      IF ((pswitch%l_psedl) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedl(i,j,k) = procs(cloud_params%i_2m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pcond) THEN
      IF (pswitch%l_pcond) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pcond(i,j,k) = procs(cloud_params%i_1m,i_cond%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pcond(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_phomr) THEN
      IF (pswitch%l_phomr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % phomr(i,j,k) = procs(graupel_params%i_1m,i_homr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % phomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

   IF (casdiags % l_pihal) THEN
      IF (pswitch%l_pihal) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pihal(i,j,k) = procs(ice_params%i_1m,i_ihal%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pihal(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nihal) THEN
      IF ((pswitch%l_pihal) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nihal(i,j,k) = procs(ice_params%i_2m,i_ihal%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nihal(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nhomr) THEN
      IF ((pswitch%l_phomr) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nhomr(i,j,k) = procs(graupel_params%i_2m,i_homr%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nhomr(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_g) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_g(i,j,k) = procs(graupel_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_g(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_r) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_r(i,j,k) = procs(rain_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_r(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_praci_i) THEN
      IF (pswitch%l_praci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % praci_i(i,j,k) = procs(ice_params%i_1m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % praci_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_g) THEN
      IF ((pswitch%l_praci) .AND. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_g(i,j,k) = procs(graupel_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_g(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_r) THEN
      IF ((pswitch%l_praci) .AND. (rain_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_r(i,j,k) = procs(rain_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_r(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nraci_i) THEN
      IF ((pswitch%l_praci) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nraci_i(i,j,k) = procs(ice_params%i_2m,i_raci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nraci_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pidps) THEN
      IF ((pswitch%l_pidps) .AND. (ice_params%l_1m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pidps(i,j,k) = procs(ice_params%i_1m,i_idps%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pidps(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_nidps) THEN
      IF ((pswitch%l_pidps) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nidps(i,j,k) = procs(ice_params%i_2m,i_idps%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nidps(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_pgaci) THEN
      IF ((pswitch%l_pgaci) .AND. (ice_params%l_1m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pgaci(i,j,k) = procs(ice_params%i_1m,i_gaci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % pgaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_ngaci) THEN
      IF ((pswitch%l_pgaci) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % ngaci(i,j,k) = procs(ice_params%i_2m,i_gaci%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % ngaci(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niics_s) THEN
      IF ((pswitch%l_piics) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niics_s(i,j,k) = procs(snow_params%i_2m,i_iics%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niics_s(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    IF (casdiags % l_niics_i) THEN
      IF ((pswitch%l_piics) .AND. (ice_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % niics_i(i,j,k) = procs(ice_params%i_2m,i_iics%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % niics_i(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    !-----------------------------------------------------
    !  aerosol stash
    !-----------------------------------------------------

    IF (casdiags % l_aact_am1) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_am1(i,j,k) = aerosol_procs(i_am1, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_am1(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 601

    IF (casdiags % l_aact_an1) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_an1(i,j,k) = aerosol_procs(i_an1, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_an1(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 602

    IF (casdiags % l_aact_am2) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_am2(i,j,k) = aerosol_procs(i_am2, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_am2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 603

    IF (casdiags % l_aact_an2) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_an2(i,j,k) = aerosol_procs(i_an2, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_an2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 604

    IF (casdiags % l_aact_am3) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_am3(i,j,k) = aerosol_procs(i_am3, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_am3(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 605

    IF (casdiags % l_aact_an3) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_an3(i,j,k) = aerosol_procs(i_an3, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_an3(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 606

    IF (casdiags % l_aact_am9) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_am9(i,j,k) = aerosol_procs(i_am9, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 607

    IF (casdiags % l_aact_an6) THEN
      IF (aswitch%l_aact) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aact_an6(i,j,k) = aerosol_procs(i_an6, i_aact%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aact_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 608

    IF (casdiags % l_aaut) THEN
      IF ((aswitch%l_aaut) .AND. (l_separate_rain)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aaut(i,j,k) = aerosol_procs(i_am5, i_aaut%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aaut(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 609

    IF (casdiags % l_aacw) THEN
      IF ((aswitch%l_aacw) .AND. (l_separate_rain)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % aacw(i,j,k) = aerosol_procs(i_am5, i_aacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % aacw(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 610

    IF (casdiags % l_arevp_am2) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_am2(i,j,k) = aerosol_procs(i_am2, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_am2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 614

    IF (casdiags % l_arevp_an2) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_an2(i,j,k) = aerosol_procs(i_an2, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_an2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 615

    IF (casdiags % l_arevp_am3) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_am3(i,j,k) = aerosol_procs(i_am3, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_am3(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 616

    IF (casdiags % l_arevp_an3) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_an3(i,j,k) = aerosol_procs(i_an3, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_an3(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 617

    IF (casdiags % l_arevp_am4) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_am4(i,j,k) = aerosol_procs(i_am4, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_am4(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 618

    IF (casdiags % l_arevp_am5) THEN
      IF ((aswitch%l_arevp) .AND. (l_separate_rain)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_am5(i,j,k) = aerosol_procs(i_am5, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_am5(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 619

    IF (casdiags % l_arevp_am6) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_am6(i,j,k) = aerosol_procs(i_am6, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_am6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 620

    IF (casdiags % l_arevp_an6) THEN
      IF (aswitch%l_arevp) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % arevp_an6(i,j,k) = aerosol_procs(i_an6, i_arevp%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % arevp_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 621

    IF (casdiags % l_dnuc_am8) THEN
      IF (aswitch%l_dnuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dnuc_am8(i,j,k) = aerosol_procs(i_am8, i_dnuc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dnuc_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 625

    IF (casdiags % l_dnuc_am6) THEN
      IF (aswitch%l_dnuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dnuc_am6(i,j,k) = aerosol_procs(i_am6, i_dnuc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dnuc_am6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 626

    IF (casdiags % l_dnuc_am9) THEN
      IF (aswitch%l_dnuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dnuc_am9(i,j,k) = aerosol_procs(i_am9, i_dnuc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dnuc_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 627

    IF (casdiags % l_dnuc_an6) THEN
      IF (aswitch%l_dnuc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dnuc_an6(i,j,k) = aerosol_procs(i_an6, i_dnuc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dnuc_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 628

    IF (casdiags % l_dsub_am2) THEN
      IF (aswitch%l_dsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsub_am2(i,j,k) = aerosol_procs(i_am2, i_dsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsub_am2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 629

    IF (casdiags % l_dsub_an2) THEN
      IF (aswitch%l_dsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsub_an2(i,j,k) = aerosol_procs(i_an2, i_dsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsub_an2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 630

    IF (casdiags % l_dsub_am6) THEN
      IF (aswitch%l_dsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsub_am6(i,j,k) = aerosol_procs(i_am6, i_dsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsub_am6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 631

    IF (casdiags % l_dsub_an6) THEN
      IF (aswitch%l_dsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsub_an6(i,j,k) = aerosol_procs(i_an6, i_dsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsub_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 632

    IF (casdiags % l_dssub_am2) THEN
      IF (aswitch%l_dssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dssub_am2(i,j,k) = aerosol_procs(i_am2, i_dssub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dssub_am2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 645

    IF (casdiags % l_dssub_an2) THEN
      IF (aswitch%l_dssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dssub_an2(i,j,k) = aerosol_procs(i_an2, i_dssub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dssub_an2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 646

    IF (casdiags % l_dssub_am6) THEN
      IF (aswitch%l_dssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dssub_am6(i,j,k) = aerosol_procs(i_am6, i_dssub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dssub_am6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 647

    IF (casdiags % l_dssub_an6) THEN
      IF (aswitch%l_dssub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dssub_an6(i,j,k) = aerosol_procs(i_an6, i_dssub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dssub_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 648

    IF (casdiags % l_dgsub_am2) THEN
      IF (aswitch%l_dgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgsub_am2(i,j,k) = aerosol_procs(i_am2, i_dgsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgsub_am2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 649

    IF (casdiags % l_dgsub_an2) THEN
      IF (aswitch%l_dgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgsub_an2(i,j,k) = aerosol_procs(i_an2, i_dgsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgsub_an2(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 650

    IF (casdiags % l_dgsub_am6) THEN
      IF (aswitch%l_dgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgsub_am6(i,j,k) = aerosol_procs(i_am6, i_dgsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgsub_am6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 651

    IF (casdiags % l_dgsub_an6) THEN
      IF (aswitch%l_dgsub) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgsub_an6(i,j,k) = aerosol_procs(i_an6, i_dgsub%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgsub_an6(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 652

    IF (casdiags % l_dhomc_am8) THEN
      IF (aswitch%l_dhomc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dhomc_am8(i,j,k) = aerosol_procs(i_am8, i_dhomc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dhomc_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 653

    IF (casdiags % l_dhomc_am7) THEN
      IF (aswitch%l_dhomc) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dhomc_am7(i,j,k) = aerosol_procs(i_am7, i_dhomc%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dhomc_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 654

    IF (casdiags % l_dhomr_am8) THEN
      IF (aswitch%l_dhomr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dhomr_am8(i,j,k) = aerosol_procs(i_am8, i_dhomr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dhomr_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 655

    IF (casdiags % l_dhomr_am7) THEN
      IF (aswitch%l_dhomr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dhomr_am7(i,j,k) = aerosol_procs(i_am7, i_dhomr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dhomr_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 656

    IF (casdiags % l_dimlt_am4) THEN
      IF (aswitch%l_dimlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dimlt_am4(i,j,k) = aerosol_procs(i_am4, i_dimlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dimlt_am4(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 657

    IF (casdiags % l_dimlt_am9) THEN
      IF (aswitch%l_dimlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dimlt_am9(i,j,k) = aerosol_procs(i_am9, i_dimlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dimlt_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 658

    IF (casdiags % l_dsmlt_am4) THEN
      IF (aswitch%l_dsmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsmlt_am4(i,j,k) = aerosol_procs(i_am4, i_dsmlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsmlt_am4(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 659

    IF (casdiags % l_dsmlt_am9) THEN
      IF (aswitch%l_dsmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsmlt_am9(i,j,k) = aerosol_procs(i_am9, i_dsmlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsmlt_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 660

    IF (casdiags % l_dgmlt_am4) THEN
      IF (aswitch%l_dgmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgmlt_am4(i,j,k) = aerosol_procs(i_am4, i_dgmlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgmlt_am4(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 661

    IF (casdiags % l_dgmlt_am9) THEN
      IF (aswitch%l_dgmlt) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgmlt_am9(i,j,k) = aerosol_procs(i_am9, i_dgmlt%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgmlt_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 662

    IF (casdiags % l_diacw_am8) THEN
      IF (aswitch%l_diacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % diacw_am8(i,j,k) = aerosol_procs(i_am8, i_diacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % diacw_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 663

    IF (casdiags % l_diacw_am7) THEN
      IF (aswitch%l_diacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % diacw_am7(i,j,k) = aerosol_procs(i_am7, i_diacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % diacw_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 664

    IF (casdiags % l_dsacw_am8) THEN
      IF (aswitch%l_dsacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsacw_am8(i,j,k) = aerosol_procs(i_am8, i_dsacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsacw_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 665

    IF (casdiags % l_dsacw_am7) THEN
      IF (aswitch%l_dsacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsacw_am7(i,j,k) = aerosol_procs(i_am7, i_dsacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsacw_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 666

    IF (casdiags % l_dgacw_am8) THEN
      IF (aswitch%l_dgacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgacw_am8(i,j,k) = aerosol_procs(i_am8, i_dgacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgacw_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 667

    IF (casdiags % l_dgacw_am7) THEN
      IF (aswitch%l_dgacw) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgacw_am7(i,j,k) = aerosol_procs(i_am7, i_dgacw%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgacw_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 668

    IF (casdiags % l_dsacr_am8) THEN
      IF (aswitch%l_dsacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsacr_am8(i,j,k) = aerosol_procs(i_am8, i_dsacr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsacr_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 669

    IF (casdiags % l_dsacr_am7) THEN
      IF (aswitch%l_dsacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsacr_am7(i,j,k) = aerosol_procs(i_am7, i_dsacr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsacr_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 670

    IF (casdiags % l_dgacr_am8) THEN
      IF (aswitch%l_dgacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgacr_am8(i,j,k) = aerosol_procs(i_am8, i_dgacr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgacr_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 671

    IF (casdiags % l_dgacr_am7) THEN
      IF (aswitch%l_dgacr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dgacr_am7(i,j,k) = aerosol_procs(i_am7, i_dgacr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dgacr_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 672

    IF (casdiags % l_draci_am8) THEN
      IF (aswitch%l_draci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % draci_am8(i,j,k) = aerosol_procs(i_am8, i_draci%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % draci_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 673

    IF (casdiags % l_draci_am7) THEN
      IF (aswitch%l_draci) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % draci_am7(i,j,k) = aerosol_procs(i_am7, i_draci%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % draci_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 674

    !---------------------------------------------------------------

    ELSE !ncall > 0

    IF (casdiags % l_psedi) THEN
      IF (pswitch%l_psedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedi(i,j,k) = casdiags % psedi(i,j,k)+                   &
                    procs(ice_params%i_1m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedi) THEN
      IF ((pswitch%l_psedi) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedi(i,j,k) = casdiags % nsedi(i,j,k)+                  &
                   procs(ice_params%i_2m,i_psedi%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_pseds) THEN
      IF (pswitch%l_pseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % pseds(i,j,k) = casdiags % pseds(i,j,k)+                  &
                  procs(snow_params%i_1m,i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nseds) THEN
      IF ((pswitch%l_pseds) .AND. (snow_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nseds(i,j,k) = casdiags % nseds(i,j,k)+                   &
                   procs(snow_params%i_2m,i_pseds%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedr) THEN
      IF (pswitch%l_psedr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedr(i,j,k) = casdiags % psedr(i,j,k)+                   &
                   procs(rain_params%i_1m,i_psedr%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedg) THEN
      IF (pswitch%l_psedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedg(i,j,k) = casdiags % psedg(i,j,k)+                   &
                procs(graupel_params%i_1m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedg) THEN
      IF ((pswitch%l_psedg) .AND. (graupel_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedg(i,j,k) = casdiags % nsedg(i,j,k)+                   &
                procs(graupel_params%i_2m,i_psedg%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_psedl) THEN
      IF (pswitch%l_psedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % psedl(i,j,k) = casdiags % psedl(i,j,k)+                   &
                  procs(cloud_params%i_1m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      END IF
    END IF

    IF (casdiags % l_nsedl) THEN
      IF ((pswitch%l_psedl) .AND. (cloud_params%l_2m)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % nsedl(i,j,k) = procs(cloud_params%i_2m,i_psedl%id,ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % nsedl(i,j,:) = ZERO_REAL_WP
      END IF
    END IF

    !---------------------------------------------
    ! aerosol stash
    !---------------------------------------------

    IF (casdiags % l_asedr_am) THEN
      IF (aswitch%l_asedr) THEN
        IF (l_separate_rain) THEN
          DO k = k_start, k_end
            kc = k - k_start + 1
            casdiags % asedr_am(i,j,k) = aerosol_procs(i_am5, i_asedr%id, ixy_inner)%column_data(kc)
          END DO
        ELSE
          DO k = k_start, k_end
            kc = k - k_start + 1
            casdiags % asedr_am(i,j,k) = aerosol_procs(i_am4, i_asedr%id, ixy_inner)%column_data(kc)
          END DO
        END IF ! separate rain aerosol
      ELSE
        casdiags % asedr_am(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 611

    IF (casdiags % l_asedr_an11) THEN
      IF ((aswitch%l_asedr) .AND. (l_passivenumbers)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedr_an11(i,j,k) = aerosol_procs(i_an11, i_asedr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedr_an11(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 612

    IF (casdiags % l_asedr_an12) THEN
      IF ((aswitch%l_asedr) .AND. (l_passivenumbers_ice)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedr_an12(i,j,k) = aerosol_procs(i_an12, i_asedr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedr_an12(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 613


    IF (casdiags % l_asedl_am4) THEN
      IF (aswitch%l_asedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedl_am4(i,j,k) = aerosol_procs(i_am4, i_asedl%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedl_am4(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 622

    IF (casdiags % l_asedl_an11) THEN
      IF ((aswitch%l_asedl) .AND. (l_passivenumbers)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedl_an11(i,j,k) = aerosol_procs(i_an11, i_asedl%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedl_an11(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 623

    IF (casdiags % l_asedl_an12) THEN
      IF ((aswitch%l_asedl) .AND. (l_passivenumbers_ice)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedl_an12(i,j,k) = aerosol_procs(i_an12, i_asedl%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedl_an12(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 624

    IF (casdiags % l_dsedi_am7) THEN
      IF (aswitch%l_dsedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedi_am7(i,j,k) = aerosol_procs(i_am7, i_dsedi%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedi_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 633

    IF (casdiags % l_dsedi_am8) THEN
      IF (aswitch%l_dsedi) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedi_am8(i,j,k) = aerosol_procs(i_am8, i_dsedi%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedi_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 634

    IF (casdiags % l_dsedi_an11) THEN
      IF ((aswitch%l_dsedi) .AND. (l_passivenumbers)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedi_an11(i,j,k) = aerosol_procs(i_an11, i_dsedi%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedi_an11(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 635

    IF (casdiags % l_dsedi_an12) THEN
      IF ((aswitch%l_dsedi) .AND. (l_passivenumbers_ice)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedi_an12(i,j,k) = aerosol_procs(i_an12, i_dsedi%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedi_an12(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 636

    IF (casdiags % l_dseds_am7) THEN
      IF (aswitch%l_dseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dseds_am7(i,j,k) = aerosol_procs(i_am7, i_dseds%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dseds_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 637

    IF (casdiags % l_dseds_am8) THEN
      IF (aswitch%l_dseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dseds_am8(i,j,k) = aerosol_procs(i_am8, i_dseds%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dseds_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 638

    IF (casdiags % l_dseds_an11) THEN
      IF ((aswitch%l_dseds) .AND. (l_passivenumbers)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dseds_an11(i,j,k) = aerosol_procs(i_an11, i_dseds%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dseds_an11(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 639

    IF ((casdiags % l_dseds_an12) .AND. (l_passivenumbers_ice)) THEN
      IF (aswitch%l_dseds) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dseds_an12(i,j,k) = aerosol_procs(i_an12, i_dseds%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dseds_an12(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 640

    IF (casdiags % l_dsedg_am7) THEN
      IF (aswitch%l_dsedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedg_am7(i,j,k) = aerosol_procs(i_am7, i_dsedg%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedg_am7(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 641

    IF (casdiags % l_dsedg_am8) THEN
      IF (aswitch%l_dsedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedg_am8(i,j,k) = aerosol_procs(i_am8, i_dsedg%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedg_am8(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 642

    IF (casdiags % l_dsedg_an11) THEN
      IF ((aswitch%l_dsedg) .AND. (l_passivenumbers)) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedg_an11(i,j,k) = aerosol_procs(i_an11, i_dsedg%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedg_an11(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 643

    IF ((casdiags % l_dsedg_an12) .AND. (l_passivenumbers_ice)) THEN
      IF (aswitch%l_dsedg) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % dsedg_an12(i,j,k) = aerosol_procs(i_an12, i_dsedg%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % dsedg_an12(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 644

    IF (casdiags % l_asedl_am9) THEN
      IF (aswitch%l_asedl) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedl_am9(i,j,k) = aerosol_procs(i_am9, i_asedl%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedl_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 675

    IF (casdiags % l_asedr_am9) THEN
      IF (aswitch%l_asedr) THEN
        DO k = k_start, k_end
          kc = k - k_start + 1
          casdiags % asedr_am9(i,j,k) = aerosol_procs(i_am9, i_asedr%id, ixy_inner)%column_data(kc)
        END DO
      ELSE
        casdiags % asedr_am9(i,j,:) = ZERO_REAL_WP
      END IF
    END IF ! stash 676

    !-----------------------------------------------

    END IF ! ncall

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE gather_process_diagnostics
END MODULE micro_main
