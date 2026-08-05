! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Top-level CASIM initialisation/finalisation (mphys_init, mphys_finalise): validates options, sets constants and lookup tables.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   One-off initialisation entry points: sets microphysics
!   switches, initialises look-up tables (gamma function, mu
!   look-up, sedimentation, aerosol) before the first timestep.
!
! Paper reference:
!   CASIM software infrastructure supporting the physics
!   described throughout Field et al. (2023) Appendix A; not
!   itself the subject of a numbered equation.
!
MODULE initialize
  USE mphys_die, ONLY: throw_mphys_error, incorrect_opt, std_msg
  USE variable_precision, ONLY: wp
  USE lookup, ONLY: set_mu_lookup, mu_i, mu_g, mu_i_sed, mu_g_sed, nmu
  USE derived_constants, ONLY: set_constants
  USE mphys_parameters, ONLY: p1, p2, p3, sp1, sp2, sp3, snow_params
  USE mphys_switches, ONLY: aerosol_option, l_warm, cloud_params, &
       rain_params, ice_params, snow_params, graupel_params,&
       iopt_act, l_g, l_sg, l_override_checks, i_am10, i_an10, &
       isol, iinsol, active_rain, active_cloud, aero_index, active_number, process_level, iopt_act, l_process
  USE gauss_casim_micro, ONLY: gaussfunclookup
  USE micro_main, ONLY : initialise_micromain, finalise_micromain
  USE sedimentation, ONLY : initialise_sedr, finalise_sedr

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='INITIALIZE'

  PUBLIC mphys_init, mphys_finalise

CONTAINS 

  SUBROUTINE mphys_init(il, iu, jl, ju, kl, ku,                 &       
       is_in, ie_in, js_in, je_in, ks_in, ke_in, l_tendency)


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: il, iu ! upper and lower i levels
    INTEGER, INTENT(IN) :: jl, ju ! upper and lower j levels
    INTEGER, INTENT(IN) :: kl, ku ! upper and lower k levels

    INTEGER, INTENT(IN), OPTIONAL :: is_in, ie_in ! upper and lower i levels which are to be used
    INTEGER, INTENT(IN), OPTIONAL :: js_in, je_in ! upper and lower j levels
    INTEGER, INTENT(IN), OPTIONAL :: ks_in, ke_in ! upper and lower k levels

    ! New optional l_tendency logical added...
    ! if true then a tendency is returned (i.e. units/s)
    ! if false then an increment is returned (i.e. units/timestep)
    LOGICAL, INTENT(IN), OPTIONAL :: l_tendency

    ! Local variables
    REAL(wp) :: tmp

    INTEGER :: i_start, i_end ! upper and lower i levels which are to be used
    INTEGER :: j_start, j_end ! upper and lower j levels
    INTEGER :: k_start, k_end ! upper and lower k levels
    LOGICAL :: l_tendency_loc

    CHARACTER(len=*), PARAMETER :: RoutineName='MPHYS_INIT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP SINGLE

    CALL check_options()
    CALL set_constants()

!$OMP END SINGLE

    ! Set grid extents to operate on
    IF (present(is_in)) i_start=is_in
    IF (present(ie_in)) i_end=ie_in
    IF (present(js_in)) j_start=js_in
    IF (present(je_in)) j_end=je_in
    IF (present(ks_in)) k_start=ks_in
    IF (present(ke_in)) k_end=ke_in

    ! if not passed in, then default to full grid
    IF (.NOT. present(is_in)) i_start=il
    IF (.NOT. present(ie_in)) i_end=iu
    IF (.NOT. present(js_in)) j_start=jl
    IF (.NOT. present(je_in)) j_end=ju
    IF (.NOT. present(ks_in)) k_start=kl
    IF (.NOT. present(ke_in)) k_end=ku
    
    IF (present(l_tendency)) THEN
      l_tendency_loc=l_tendency
    ELSE
      l_tendency_loc=.TRUE.
    END IF

    CALL initialise_micromain(il, iu, jl, ju, kl, ku, i_start, i_end, j_start, j_end, k_start, k_end, l_tendency_loc)

    CALL initialise_sedr()

!$OMP SINGLE
    CALL initialise_lookup_tables()
!$OMP END SINGLE

    !Gaussfunc
    CALL gaussfunclookup(snow_params%id, tmp, a=snow_params%fix_mu, b=snow_params%b_x)

!$OMP SINGLE
    !Initialise the gamma function so not calced on every timestep
    CALL gamma_initialize()
!$OMP END SINGLE

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE mphys_init

  SUBROUTINE mphys_finalise()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='MPHYS_FINALISE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    CALL finalise_sedr()
    CALL finalise_micromain()

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE mphys_finalise

  SUBROUTINE initialise_lookup_tables()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='INITIALISE_LOOKUP_TABLES'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !------------------------------------------------------
    ! Look up tables... (These need to be extended)
    !------------------------------------------------------
    ! mu lookup

    ALLOCATE(mu_i(nmu))
    ALLOCATE(mu_g(nmu))
    CALL set_mu_lookup(p1, p2, p3, mu_i, mu_g)
    ALLOCATE(mu_i_sed(nmu))
    ALLOCATE(mu_g_sed(nmu))
    CALL set_mu_lookup(sp1, sp2, sp3, mu_i_sed, mu_g_sed)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE initialise_lookup_tables

  ! Check that the options that have been selected are consitent
  SUBROUTINE check_options()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='CHECK_OPTIONS'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (.NOT. l_override_checks) THEN
      IF (l_warm .AND. l_process) THEN
        WRITE(std_msg, '(A)') 'processing does not currently work with l_warm=.true.'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      IF (.NOT. ice_params%l_1m .AND. .NOT. l_warm) THEN
        WRITE(std_msg, '(A)') 'l_warm must be true if not using ice'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      IF (cloud_params%l_2m .AND. aerosol_option == 0     &
           .AND. (iopt_act== 1 .OR. iopt_act==3)) THEN
        WRITE(std_msg, '(A)') 'for double moment cloud you must have aerosol_option>0'// &
             'or else activation should be independent of aerosol'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
             
      END IF

      IF (.NOT. cloud_params%l_2m .AND. aerosol_option > 0) THEN
        WRITE(std_msg, '(A)') 'aerosol_option must be 0 if not using double moment microphysics'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      IF (l_g .AND. .NOT. l_sg) THEN
        WRITE(std_msg, '(A)') 'Cannot run with graupel but not with snow'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
        l_g=.FALSE.
      END IF

      ! Some options are not yet working or well tested so don't let anyone use these

      IF (i_am10 > 0 .OR. i_an10 > 0) THEN
        WRITE(std_msg, '(A)') 'Accumulation mode dust is not yet used.'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      ! Aerosol consistency
      IF (active_rain(isol) .AND. .NOT. active_cloud(isol)) THEN
        WRITE(std_msg, '(A)') 'active_rain(isol) .and. .not. active_cloud(isol)'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      ! Aerosol consistency
      IF (aero_index%i_accum == 0 .AND. aero_index%i_coarse==0 .AND. &
         (iopt_act==1 .OR. iopt_act==3) )                            THEN
        WRITE(std_msg, '(A)') 'Must have accumulation or coarse mode aerosol '//&
                              'for chosen activation option'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      ! Aerosol consistency
      IF (aero_index%i_accum == 0 .AND. aero_index%nccn /= 1 .AND. &
          active_number(isol))                                     THEN
      
        WRITE(std_msg, '(A)') 'Soluble modes must only be accumulation mode '//&
                              'with soluble active_number'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF

      ! Aerosol consistency
      IF (aero_index%nin > 0) THEN 
         IF (aero_index%i_accum_dust == 0 .AND. aero_index%nin /= 1 .AND. &
              active_number(iinsol))                                  THEN
            WRITE(std_msg, '(A)') 'Dust modes must only be accumulation mode '//&
                 'with insoluble active_number'
            CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                 std_msg)
         END IF
      END IF

      ! Aerosol processing consistency
      IF (process_level > 0 .AND. iopt_act < 3) THEN
        WRITE(std_msg, '(A)') 'If processing aerosol, must use higher '//&
                              'level activation code, i.e check iopt_act'
        CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                               std_msg)
      END IF
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE check_options
  
  SUBROUTINE gamma_initialize()
    ! Purpose: When CASIM is used in 1M or 2M mode with a fixed shape parameter, 
    !          all the gamma functions can be calculated at the beginning of the job.
    !          This routine does this calculation so that results can be used in 
    !          sedimentation, lookup...
    
    !use mphys_parameters, only: cloud_params, rain_params, ice_params, snow_params, &
    !     graupel_params
    USE special, ONLY: Gammafunc

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='GAMMA_INITIALIZE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! moment gamma functions

    ! Cloud
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    cloud_params%gam_1_mu_p1 = GammaFunc(1.0_wp+cloud_params%fix_mu+cloud_params%p1)
    cloud_params%gam_1_mu_p2 = GammaFunc(1.0_wp+cloud_params%fix_mu+cloud_params%p2)
    cloud_params%gam_1_mu = GammaFunc(1.0_wp+cloud_params%fix_mu)
    ! Moment exponent for lamba
    cloud_params%inv_p1_p2 = (1.0_wp)/(cloud_params%p1 - cloud_params%p2)
    cloud_params%inv_p1 = (1.0_wp)/(cloud_params%p1)

    ! Rain
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    rain_params%gam_1_mu_p1 = GammaFunc(1.0_wp+rain_params%fix_mu+rain_params%p1)
    rain_params%gam_1_mu_p2 = GammaFunc(1.0_wp+rain_params%fix_mu+rain_params%p2)
    rain_params%gam_1_mu = GammaFunc(1.+rain_params%fix_mu)
    ! Moment exponent for lamba
    rain_params%inv_p1_p2 = (1.0_wp)/(rain_params%p1 - rain_params%p2)
    rain_params%inv_p1 = (1.0_wp)/(rain_params%p1)

    ! Ice
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    ice_params%gam_1_mu_p1 = GammaFunc(1.0_wp+ice_params%fix_mu+ice_params%p1)
    ice_params%gam_1_mu_p2 = GammaFunc(1.0_wp+ice_params%fix_mu+ice_params%p2)
    ice_params%gam_1_mu = GammaFunc(1.0_wp+ice_params%fix_mu)
    ! Moment exponent for lamba
    ice_params%inv_p1_p2 = (1.0_wp)/(ice_params%p1 - ice_params%p2)
    ice_params%inv_p1 = (1.0_wp)/(ice_params%p1)

    ! Snow
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    snow_params%gam_1_mu_p1 = GammaFunc(1.0_wp+snow_params%fix_mu+snow_params%p1)
    snow_params%gam_1_mu_p2 = GammaFunc(1.0_wp+snow_params%fix_mu+snow_params%p2)
    snow_params%gam_1_mu = GammaFunc(1.0_wp+snow_params%fix_mu)
    ! Moment exponent for lamba
    snow_params%inv_p1_p2 = (1.0_wp)/(snow_params%p1 - snow_params%p2)
    snow_params%inv_p1 = (1.0_wp)/(snow_params%p1)

    ! Graupel
    ! gamma function constants (used in Lookup, get_lam_n0_1M and 2M)
    graupel_params%gam_1_mu_p1 = GammaFunc(1.0_wp+graupel_params%fix_mu+graupel_params%p1)
    graupel_params%gam_1_mu_p2 = GammaFunc(1.0_wp+graupel_params%fix_mu+graupel_params%p2)
    graupel_params%gam_1_mu = GammaFunc(1.0_wp+graupel_params%fix_mu)
    ! Moment exponent for lamba
    graupel_params%inv_p1_p2 = (1.0_wp)/(graupel_params%p1 - graupel_params%p2)
    graupel_params%inv_p1 = (1.0_wp)/(graupel_params%p1)   

    ! gamma functions for sedimentation moments and ventilation
    
    ! Cloud
    ! gamma function constants (used in sedimentation)
    cloud_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp1)
    cloud_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp1+cloud_params%b_x)
    cloud_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp2)
    cloud_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+cloud_params%fix_mu+cloud_params%sp2+cloud_params%b_x)
    ! exponent of lambda for sedimentation
    cloud_params%exp_1_mu_sp1 = 1.0+cloud_params%fix_mu+cloud_params%sp1
    cloud_params%exp_1_mu_sp1_bx = 1.0 + cloud_params%fix_mu + cloud_params%sp1 + cloud_params%b_x
    cloud_params%exp_1_mu_sp2 = 1.0+cloud_params%fix_mu+cloud_params%sp2
    cloud_params%exp_1_mu_sp2_bx = 1.0 + cloud_params%fix_mu + cloud_params%sp2 + cloud_params%b_x
    ! gamma function used in ventilation (not really needed for cloud, just added for completeness)
    cloud_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*cloud_params%b_x+cloud_params%fix_mu+2.5)

    ! Rain
    ! gamma function constants (used in sedimentation)
    rain_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1)
    rain_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1+rain_params%b_x)
    rain_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2)
    rain_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2+rain_params%b_x)
    ! gamma function constants for Abel and shipway (used in sedimentation)
    rain_params%gam_1_mu_sp1_b2x = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp1+rain_params%b2_x)
    rain_params%gam_1_mu_sp2_b2x = & 
         Gammafunc(1.0_wp+rain_params%fix_mu+rain_params%sp2+rain_params%b2_x)
    ! exponent of lambda for sedimentation
    rain_params%exp_1_mu_sp1 = 1.0+rain_params%fix_mu+rain_params%sp1
    rain_params%exp_1_mu_sp1_bx = 1.0 + rain_params%fix_mu + rain_params%sp1 + rain_params%b_x
    rain_params%exp_1_mu_sp2 = 1.0+rain_params%fix_mu+rain_params%sp2
    rain_params%exp_1_mu_sp2_bx = 1.0 + rain_params%fix_mu + rain_params%sp2 + rain_params%b_x
    ! exponent of lambda for Abel and Shipway sedimentation
    rain_params%exp_1_mu_sp1_b2x = &
         1.0 + rain_params%fix_mu + rain_params%sp1 + rain_params%b2_x
    rain_params%exp_1_mu_sp2_b2x = &
         1.0 + rain_params%fix_mu + rain_params%sp2 + rain_params%b2_x
    ! gamma function used in ventilation
    rain_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*rain_params%b_x+rain_params%fix_mu+2.5)

    ! Ice
    ! gamma function constants (used in sedimentation)
    ice_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp1)
    ice_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp1+ice_params%b_x)
    ice_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp2)
    ice_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+ice_params%fix_mu+ice_params%sp2+ice_params%b_x)
    ! exponent of lambda for sedimentation
    ice_params%exp_1_mu_sp1 = 1.0+ice_params%fix_mu+ice_params%sp1
    ice_params%exp_1_mu_sp1_bx = 1.0 + ice_params%fix_mu + ice_params%sp1 + ice_params%b_x
    ice_params%exp_1_mu_sp2 = 1.0+ice_params%fix_mu+ice_params%sp2
    ice_params%exp_1_mu_sp2_bx = 1.0 + ice_params%fix_mu + ice_params%sp2 + ice_params%b_x
    ! gamma function used in ventilation
    ice_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*ice_params%b_x+ice_params%fix_mu+2.5)

    
    ! Snow
    ! gamma function constants (used in sedimentation)
    snow_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp1)
    snow_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp1+snow_params%b_x)
    snow_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp2)
    snow_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+snow_params%fix_mu+snow_params%sp2+snow_params%b_x)
    ! exponent of lambda for sedimentation
    snow_params%exp_1_mu_sp1 = 1.0+snow_params%fix_mu+snow_params%sp1
    snow_params%exp_1_mu_sp1_bx = 1.0 + snow_params%fix_mu + snow_params%sp1 + snow_params%b_x
    snow_params%exp_1_mu_sp2 = 1.0+snow_params%fix_mu+snow_params%sp2
    snow_params%exp_1_mu_sp2_bx = 1.0 + snow_params%fix_mu + snow_params%sp2 + snow_params%b_x
    ! gamma function used in ventilation
    snow_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*snow_params%b_x+snow_params%fix_mu+2.5)
    
    
    ! Graupel
    ! gamma function constants (used in sedimentation)
    graupel_params%gam_1_mu_sp1 = Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp1)
    graupel_params%gam_1_mu_sp1_bx = & 
         Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp1+graupel_params%b_x)
    graupel_params%gam_1_mu_sp2 = Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp2)
    graupel_params%gam_1_mu_sp2_bx = & 
         Gammafunc(1.0_wp+graupel_params%fix_mu+graupel_params%sp2+graupel_params%b_x)
    ! exponent of lambda for sedimentation
    graupel_params%exp_1_mu_sp1 = 1.0+graupel_params%fix_mu+graupel_params%sp1
    graupel_params%exp_1_mu_sp1_bx = 1.0 + graupel_params%fix_mu + graupel_params%sp1 + graupel_params%b_x
    graupel_params%exp_1_mu_sp2 = 1.0+graupel_params%fix_mu+graupel_params%sp2
    graupel_params%exp_1_mu_sp2_bx = 1.0 + graupel_params%fix_mu + graupel_params%sp2 + graupel_params%b_x
    ! gamma function used in ventilation
    graupel_params%gam_0p5bx_mu_2p5 = GammaFunc(.5*graupel_params%b_x+graupel_params%fix_mu+2.5)

! gamma functions for the calculation of the ice_accretion using the functions sweepout and binary collection

    ! Cloud 
    ! gamma function (used in sweepout)
    cloud_params%gam_3_bx_mu = GammaFunc(3.0+cloud_params%b_x+cloud_params%fix_mu)
    cloud_params%gam_3_bx_mu_dx = GammaFunc(3.0+cloud_params%b_x+cloud_params%fix_mu+cloud_params%d_x)
    ! gamma function constants (used in binary_collection)
    cloud_params%gam_2_mu = GammaFunc(2.0_wp+cloud_params%fix_mu)
    cloud_params%gam_3_mu = GammaFunc(3.0_wp+cloud_params%fix_mu)
    cloud_params%gam_2_mu_dx = GammaFunc(2.0_wp+cloud_params%fix_mu+cloud_params%d_x)
    cloud_params%gam_3_mu_dx = GammaFunc(3.0_wp+cloud_params%fix_mu+cloud_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    cloud_params%gam_1_mu_dx_bx = GammaFunc(1.0 + cloud_params%fix_mu + cloud_params%d_x + cloud_params%b_x)
    cloud_params%gam_1_mu_dx = GammaFunc(1.0 + cloud_params%fix_mu + cloud_params%d_x) 
    
    ! Rain
    ! gamma function (used in sweepout)
    rain_params%gam_3_bx_mu = GammaFunc(3.0+rain_params%b_x+rain_params%fix_mu)
    rain_params%gam_3_bx_mu_dx = GammaFunc(3.0+rain_params%b_x+rain_params%fix_mu+rain_params%d_x)
    ! gamma function constants (used in binary_collection)
    rain_params%gam_2_mu = GammaFunc(2.0_wp+rain_params%fix_mu)
    rain_params%gam_3_mu = GammaFunc(3.0_wp+rain_params%fix_mu)
    rain_params%gam_2_mu_dx = GammaFunc(2.0_wp+rain_params%fix_mu+rain_params%d_x)
    rain_params%gam_3_mu_dx = GammaFunc(3.0_wp+rain_params%fix_mu+rain_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    rain_params%gam_1_mu_dx_bx = GammaFunc(1.0 + rain_params%fix_mu + rain_params%d_x + rain_params%b_x)
    rain_params%gam_1_mu_dx = GammaFunc(1.0 + rain_params%fix_mu + rain_params%d_x) 
    
    ! Ice
    ! gamma function (used in sweepout)
    ice_params%gam_3_bx_mu = GammaFunc(3.0+ice_params%b_x+ice_params%fix_mu)
    ice_params%gam_3_bx_mu_dx = GammaFunc(3.0+ice_params%b_x+ice_params%fix_mu+ice_params%d_x)
    ! gamma function constants (used in binary_collection)
    ice_params%gam_2_mu = GammaFunc(2.0_wp+ice_params%fix_mu)
    ice_params%gam_3_mu = GammaFunc(3.0_wp+ice_params%fix_mu)
    ice_params%gam_2_mu_dx = GammaFunc(2.0_wp+ice_params%fix_mu+ice_params%d_x)
    ice_params%gam_3_mu_dx = GammaFunc(3.0_wp+ice_params%fix_mu+ice_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    ice_params%gam_1_mu_dx_bx = GammaFunc(1.0 + ice_params%fix_mu + ice_params%d_x + ice_params%b_x)
    ice_params%gam_1_mu_dx = GammaFunc(1.0 + ice_params%fix_mu + ice_params%d_x) 
    
    ! Snow
    ! gamma function (used in sweepout)
    snow_params%gam_3_bx_mu = GammaFunc(3.0+snow_params%b_x+snow_params%fix_mu)
    snow_params%gam_3_bx_mu_dx = GammaFunc(3.0+snow_params%b_x+snow_params%fix_mu+snow_params%d_x)
    ! gamma function constants (used in binary_collection)
    snow_params%gam_2_mu = GammaFunc(2.0_wp+snow_params%fix_mu)
    snow_params%gam_3_mu = GammaFunc(3.0_wp+snow_params%fix_mu)
    snow_params%gam_2_mu_dx = GammaFunc(2.0_wp+snow_params%fix_mu+snow_params%d_x)
    snow_params%gam_3_mu_dx = GammaFunc(3.0_wp+snow_params%fix_mu+snow_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    snow_params%gam_1_mu_dx_bx = GammaFunc(1.0 + snow_params%fix_mu + snow_params%d_x + snow_params%b_x)
    snow_params%gam_1_mu_dx = GammaFunc(1.0 + snow_params%fix_mu + snow_params%d_x) 

    ! Graupel
    ! gamma function (used in sweepout)
    graupel_params%gam_3_bx_mu = GammaFunc(3.0+graupel_params%b_x+graupel_params%fix_mu)
    graupel_params%gam_3_bx_mu_dx = GammaFunc(3.0+graupel_params%b_x+graupel_params%fix_mu+graupel_params%d_x)
    ! gamma function constants (used in binary_collection)
    graupel_params%gam_2_mu = GammaFunc(2.0_wp+graupel_params%fix_mu)
    graupel_params%gam_3_mu = GammaFunc(3.0_wp+graupel_params%fix_mu)
    graupel_params%gam_2_mu_dx = GammaFunc(2.0_wp+graupel_params%fix_mu+graupel_params%d_x)
    graupel_params%gam_3_mu_dx = GammaFunc(3.0_wp+graupel_params%fix_mu+graupel_params%d_x)
    ! gamma function used in Vx and Vy (binary_collection)
    graupel_params%gam_1_mu_dx_bx = GammaFunc(1.0 + graupel_params%fix_mu + graupel_params%d_x + graupel_params%b_x)
    graupel_params%gam_1_mu_dx = GammaFunc(1.0 + graupel_params%fix_mu + graupel_params%d_x) 
    

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
    
  END SUBROUTINE gamma_initialize

END MODULE initialize
