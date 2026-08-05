! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Aerosol mode partitioning and CCN/dust activation (Abdul-Razzak & Ghan 2000), used when aerosol_option>0.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Shared helper routines for transferring aerosol mass/number
!   between the free-aerosol and in-cloud/in-hydrometeor aerosol
!   prognostics during activation, nucleation and collection.
!
! Paper reference:
!   Field et al. (2023) Sec. 1 describes activated aerosol being
!   carried as an in-cloud prognostic transported by dynamics
!   and sedimentation (citing Ghan & Easter, 2006); the specific
!   bookkeeping routines here are not given numbered equations.
!
MODULE aerosol_routines
  USE mphys_die, ONLY: throw_mphys_error, warn, std_msg
  USE type_aerosol, ONLY: aerosol_phys, aerosol_chem, aerosol_active
  USE variable_precision, ONLY: wp
  USE mphys_constants, ONLY: Mw, zetasa, Ru, rhow, g, Lv, mp_eps, cp, Rd, Rv, ka, Dv, pi
  USE special, ONLY: casim_erfc, erfinv
  USE mphys_switches, ONLY: i_am4, i_nl, i_nr, &
       i_am5, i_am7, i_am8, i_am9, i_an11, i_an12, i_ni, i_ns, i_ng, i_ql, &
       i_qr, i_qi, i_qs, i_qg, l_active_inarg2000, aero_index, l_warm,  &
       l_process, l_passivenumbers,              &
       l_passivenumbers_ice, l_separate_rain, l_ukca_casim, l_bypass_which_mode
  USE thresholds, ONLY: nr_tidy, nl_tidy, ni_tidy, ccn_tidy, qr_tidy, aeromass_small, aeronumber_small
  USE mphys_parameters, ONLY: sigma_arc, nz
  USE lognormal_funcs, ONLY: MNtoRm ! DPG - added this for MNtoRm since was
                                    ! causing circular conflicts as wanted
                                    ! to use it in which_mode_to_use.F90

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='AEROSOL_ROUTINES'

  LOGICAL :: l_warned_of_prag_hack = .FALSE.

  PRIVATE

  PUBLIC aerosol_active, aerosol_phys, aerosol_chem, abdulRazzakGhan2000, invert_partial_moment, upperpartial_moment_logn, &
       invert_partial_moment_approx, invert_partial_moment_betterapprox, examine_aerosol, allocate_aerosol, &
       deallocate_aerosol, abdulRazzakGhan2000_dust
CONTAINS
  !
  ! Allocate space for aerosol_chem and aerosol_phys
  !
  SUBROUTINE allocate_aerosol(aerophys, aerochem, nmodes, initphys, initchem)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    TYPE(aerosol_phys), INTENT(INOUT) :: aerophys(:)
    TYPE(aerosol_chem), INTENT(INOUT) :: aerochem(:)
    INTEGER, INTENT(IN) :: nmodes    ! number of modes
    TYPE(aerosol_phys), INTENT(IN), OPTIONAL :: initphys
    TYPE(aerosol_chem), INTENT(IN), OPTIONAL :: initchem

    ! Local variables

    INTEGER :: k
    CHARACTER(len=*), PARAMETER :: RoutineName='ALLOCATE_AEROSOL'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !------------------------
    ! First allocate aerophys
    !------------------------
    DO k=1, size(aerophys)
      aerophys(k)%nmodes=nmodes
      ALLOCATE(aerophys(k)%N(max(1,nmodes)))
      ALLOCATE(aerophys(k)%M(max(1,nmodes)))
      ALLOCATE(aerophys(k)%rd(max(1,nmodes)))
      ALLOCATE(aerophys(k)%sigma(max(1,nmodes)))
      ALLOCATE(aerophys(k)%rpart(max(1,nmodes)))

      IF (present(initphys)) aerophys(k)=initphys
    END DO

    !------------------------
    ! Then allocate aerochem
    !------------------------
    DO k=1, size(aerochem)
      aerochem(k)%nmodes=nmodes
      ALLOCATE(aerochem(k)%vantHoff(max(1,nmodes)))
      ALLOCATE(aerochem(k)%massMole(max(1,nmodes)))
      ALLOCATE(aerochem(k)%density(max(1,nmodes)))
      ALLOCATE(aerochem(k)%epsv(max(1,nmodes)))
      ALLOCATE(aerochem(k)%beta(max(1,nmodes)))
      ALLOCATE(aerochem(k)%bk(max(1,nmodes)))

      IF (present(initphys)) aerochem(k)=initchem
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE allocate_aerosol

  !
  ! Deallocate space for aerosol_chem and aerosol_phys
  !
  SUBROUTINE deallocate_aerosol(aerophys, aerochem)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE


    ! Subroutine arguments

    TYPE(aerosol_phys), INTENT(INOUT) :: aerophys(:)
    TYPE(aerosol_chem), INTENT(INOUT) :: aerochem(:)

    ! Local variables

    INTEGER :: k
    CHARACTER(len=*), PARAMETER :: RoutineName='DEALLOCATE_AEROSOL'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !------------------------
    ! First deallocate aerophys
    !------------------------
    DO k=1, size(aerophys)
      DEALLOCATE(aerophys(k)%N)
      DEALLOCATE(aerophys(k)%M)
      DEALLOCATE(aerophys(k)%rd)
      DEALLOCATE(aerophys(k)%sigma)
      DEALLOCATE(aerophys(k)%rpart)
    END DO

    !------------------------
    ! Then deallocate aerochem
    !------------------------
    DO k=1, size(aerochem)
      DEALLOCATE(aerochem(k)%vantHoff)
      DEALLOCATE(aerochem(k)%massMole)
      DEALLOCATE(aerochem(k)%density)
      DEALLOCATE(aerochem(k)%epsv)
      DEALLOCATE(aerochem(k)%beta)
      DEALLOCATE(aerochem(k)%bk)
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE deallocate_aerosol

  !
  !< This routine examines the aerosol input and derives physical and
  !< chemical properties in
  !<     - aerophys: information about physical composition of aerosol
  !<     - aerochem: information about chemical composition of aerosol (currently not derived)
  !<     - aeroact:  information about physical composition activated aerosol
  !<     - dustphys: information about physical composition of dust
  !<     - dustchem: information about chemical composition of dust (currently not derived)
  !<     - dustact:  information about physical composition activated dust
  !<     - aeroice:  information about physical composition activated aerosol in ice
  !<     - dustliq:  information about physical composition activated dust in liquid water
  !
  SUBROUTINE examine_aerosol(aerofields, qfields, aerophys, aerochem, aeroact, &
       dustphys, dustchem, dustact, aeroice, dustliq, icall)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN), TARGET :: aerofields(:,:)
    TYPE(aerosol_phys), INTENT(INOUT), TARGET :: aerophys(:)
    TYPE(aerosol_chem), INTENT(IN), TARGET :: aerochem(:)
    TYPE(aerosol_active), INTENT(INOUT) :: aeroact(:)
    TYPE(aerosol_phys), INTENT(INOUT), TARGET, OPTIONAL :: dustphys(:)
    TYPE(aerosol_chem), INTENT(IN), TARGET, OPTIONAL :: dustchem(:)
    TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: dustact(:)
    TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: aeroice(:)
    TYPE(aerosol_active), INTENT(INOUT), OPTIONAL :: dustliq(:)
    REAL(wp), INTENT(INOUT), TARGET :: qfields(:,:)
    INTEGER, INTENT(IN), OPTIONAL :: icall

    ! Local variables

    REAL(wp) :: mac, mar, mad, maai, madl, cloud_number, rain_number
    REAL(wp) :: cloud_mass, rain_mass
    REAL(wp) :: ice_number, snow_number, graupel_number
    REAL(wp) :: ice_mass, snow_mass, graupel_mass
    REAL(wp) :: mact2
    REAL(wp) :: density
    INTEGER :: k

    REAL(wp) :: rm_arc
    CHARACTER(2) :: chcall
    REAL(wp) :: ratio_l, ratio_r, nratio_l, nratio_r
    REAL(wp) :: ratio_i, ratio_s, ratio_g, nratio_i, nratio_s, nratio_g
!   real(wp) :: mratio_i, mratio_s, mratio_g
    REAL(wp) :: ratio_dil, ratio_dli, ratio_ali, ratio_ail

    INTEGER :: imode
    REAL(wp) :: mode_N, mode_M, mode_K
    REAL(wp) :: nitot, ntot, nhtot, mitot
    REAL(wp) :: mact, rcrit2

    LOGICAL :: l_condition, l_condition_r

    CHARACTER(len=*), PARAMETER :: RoutineName='EXAMINE_AEROSOL'

    REAL, PARAMETER :: aero_mact_mean_max=2.3e-15 !kg 3xsigma+mean(0.3micron) density 1777
    REAL, PARAMETER :: dust_mact_mean_max=3.1e-13 !kg 3xsigma+mean(1micron) density 1777
    ! When l_bypass_whichmode is True (or no larger modes exist to move aerosol to)
    ! aerosol mean sizes in hydrometeors can become so large that they lead to instability
    ! during hydrometeor sedimentation. Limiting the maximum mean mass stops these edge
    ! effects dominating the aerosol and eventually hydrometeor fields.

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Some diagnostic strings
    IF (present(icall)) THEN
      WRITE(chcall,'(A1,i1)') '_', icall
    ELSE
      chcall=''
    END IF

    DO k=1, nz
        cloud_number=qfields(k,i_nl)
        rain_number=qfields(k,i_nr)
        cloud_mass=qfields(k,i_ql)
        rain_mass=qfields(k,i_qr)
        density=aerochem(k)%density(1)

        DO imode=1, aero_index%nccn
          mode_N=aerofields(k,aero_index%ccn_n(imode))
          mode_M=aerofields(k,aero_index%ccn_m(imode))
          mode_K=aerofields(k,aero_index%ccn_k(imode))
          IF (mode_N > ccn_tidy .AND. mode_M > ccn_tidy*epsilon(1.0_wp)) THEN ! FiX This
            aerophys(k)%N(imode)=mode_N
            aerophys(k)%M(imode)=mode_M
            aerophys(k)%rd(imode)=MNtoRm(mode_M,mode_N,density,aerophys(k)%sigma(imode))
            aerochem(k)%bk(imode)=mode_K
          ELSE
            aerophys(k)%N(imode)=0.0
            aerophys(k)%M(imode)=0.0
            aerophys(k)%rd(imode)=0.0
            aerochem(k)%bk(imode)=0.0
          END IF
        END DO

        aeroact(k)%nact=0.0
        aeroact(k)%mact=0.0
        aeroact(k)%rcrit=999.0
        aeroact(k)%mact_mean=0.0
        aeroact(k)%nact2=0.0
        aeroact(k)%nact1=0.0
        aeroact(k)%mact1=0.0
        aeroact(k)%rcrit1=999.0
        aeroact(k)%nact2=0.0
        aeroact(k)%rcrit2=999.0
        aeroact(k)%mact2=0.0
        aeroact(k)%mact2_mean=0.0
        aeroact(k)%mact1_mean=0.0

        aeroact(k)%nratio1=0.0
        aeroact(k)%nratio2=0.0

        IF (l_process) THEN
          ! Examine activated aerosol
          mac=aerofields(k,i_am4)

          IF (l_separate_rain) THEN
            mar=aerofields(k,i_am5)
            mact=mac+mar
            l_condition=mac > epsilon(1.0_wp) .AND. cloud_number > nl_tidy
            l_condition_r=mar > epsilon(1.0_wp) .AND. rain_number > nr_tidy .AND. rain_mass > qr_tidy
          ELSE
            mact=mac
            l_condition=mac > epsilon(1.0_wp) .AND. cloud_number + rain_number > nl_tidy
            l_condition_r=mac > epsilon(1.0_wp) .AND. rain_number > nr_tidy .AND. rain_mass > qr_tidy
          END IF

          IF (l_warm)THEN
            maai=0.0
          ELSE
            maai=aerofields(k,i_am8)
          END IF

          ratio_ali=mact/(mact+maai+epsilon(1.0_wp))

          nhtot=cloud_number+rain_number

          IF (l_passivenumbers) THEN
            ntot=aerofields(k,i_an11)*ratio_ali
          ELSE
            ntot=nhtot
          END IF

          l_condition=l_condition .AND. ntot > epsilon(1.0_wp)

          nhtot=nhtot+epsilon(1.0_wp) ! prevent possible later division by zero
          nratio_l=cloud_number/nhtot
          nratio_r=rain_number/nhtot

          IF (l_condition) THEN
            aeroact(k)%nact=ntot
            aeroact(k)%mact=mac
            aeroact(k)%rcrit=0.0
            aeroact(k)%mact_mean=aeroact(k)%mact /(aeroact(k)%nact + epsilon(ntot))
            aeroact(k)%mact_mean=min(aero_mact_mean_max, aeroact(k)%mact_mean)
            ! Get mean radius of distribution
            rm_arc=MNtoRm(aeroact(k)%mact,aeroact(k)%nact,density,sigma_arc)

            aeroact(k)%rd=rm_arc
            aeroact(k)%sigma=sigma_arc

          END IF

          IF (l_condition_r) THEN
            IF (l_separate_rain) THEN! Separate rain category
              rcrit2=0
              mact2=mar
            ELSE IF (.NOT. l_condition) THEN ! No cloud here
              rcrit2=0
              mact2=mac
            ELSE ! Diagnostic partitioning of aerosol between cloud and rain.
              rcrit2=0.0
              mact2=aeroact(k)%mact*rain_mass/(cloud_mass + rain_mass + epsilon(nhtot))
            END IF

            aeroact(k)%nact2=ntot*nratio_r
            aeroact(k)%rcrit2=rcrit2
            aeroact(k)%mact2=mact2
            aeroact(k)%mact2_mean=aeroact(k)%mact2/(aeroact(k)%nact2 + epsilon(ntot))
            aeroact(k)%mact2_mean=min(aero_mact_mean_max, aeroact(k)%mact2_mean)

          END IF

          aeroact(k)%nact1=max(0.0_wp, aeroact(k)%nact-aeroact(k)%nact2)
          aeroact(k)%mact1=max(0.0_wp, aeroact(k)%mact-aeroact(k)%mact2)
          aeroact(k)%rcrit1=0.0
          aeroact(k)%mact1_mean=aeroact(k)%mact1/(aeroact(k)%nact1+epsilon(mar))
          aeroact(k)%mact1_mean=min(aero_mact_mean_max, aeroact(k)%mact1_mean)

          IF (cloud_number > epsilon(1.0_wp)) THEN
            aeroact(k)%nratio1=max(0.0,min(1.0,aeroact(k)%nact1/cloud_number))
          END IF
          IF (rain_number > epsilon(1.0_wp)) THEN
            aeroact(k)%nratio2=max(0.0,min(1.0,aeroact(k)%nact2/rain_number))
          END IF
          
        END IF
    END DO

    IF (.NOT. l_warm) THEN
      ! Activated dust
      DO k=1,nz
        density=dustchem(k)%density(1)
        ice_number=qfields(k,i_ni)
        snow_number=qfields(k,i_ns)
        nhtot=ice_number+snow_number
        IF (i_ng > 0) THEN
          graupel_number=qfields(k,i_ng)
          nhtot=nhtot+graupel_number
        END IF

        IF (l_process) THEN
          mad=aerofields(k,i_am7)
          madl=aerofields(k,i_am9)
          ratio_dil=mad/(madl+mad+epsilon(mad))
        END IF

        IF (l_passivenumbers_ice) THEN
          nitot=aerofields(k,i_an12)*ratio_dil
        ELSE
          nitot=nhtot
        END IF

        nhtot=nhtot+epsilon(nhtot) ! prevent possible later division by zero
        nratio_i=ice_number/nhtot
        nratio_s=snow_number/nhtot
        nratio_g=graupel_number/nhtot

        !        mratio_i=ice_mass/mitot
        !        mratio_s=snow_mass/mitot
        !        mratio_g=graupel_mass/mitot

        ! Use nratios...
        ratio_i=nratio_i
        ratio_s=nratio_s
        ratio_g=nratio_g

        ! Examine interstitial dust
        DO imode=1, aero_index%nin
          mode_N=aerofields(k,aero_index%in_n(imode))
          mode_M=aerofields(k,aero_index%in_m(imode))
          IF (mode_m < 0 .OR. mode_n < 0 ) THEN

            IF ( .NOT. l_warned_of_prag_hack ) THEN

              WRITE(std_msg,*)   'Problem! Using pragmatic hack to continue.  '//       &
                                 'If you see this message, report it to Ben Shipway!!!'

              CALL throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)
              l_warned_of_prag_hack = .TRUE.
            END IF

            IF (mode_m < 0) mode_m=aeromass_small
            IF (mode_n < 0) mode_n=aeronumber_small
          END IF

          dustphys(k)%N(imode)=mode_N
          dustphys(k)%M(imode)=mode_M

          dustphys(k)%rd(imode)=MNtoRm(mode_M,mode_N,density,dustphys(k)%sigma(imode))
        END DO

        ! Examine activated dust
        ! Initialize to zero/defaults
        dustact(k)%nact=0.0
        dustact(k)%mact=0.0
        dustact(k)%rcrit=999.0
        dustact(k)%mact_mean=0.0
        dustact(k)%nact1=0.0
        dustact(k)%nratio1=0.0
        dustact(k)%mact1=0.0
        dustact(k)%rcrit1=999.0
        dustact(k)%mact1_mean=0.0
        dustact(k)%mact2=0.0
        dustact(k)%nact2=0.0
        dustact(k)%nratio2=0.0
        dustact(k)%rcrit2=999.0
        dustact(k)%mact2_mean=0.0
        dustact(k)%mact3=0.0
        dustact(k)%nact3=0.0
        dustact(k)%nratio3=0.0
        dustact(k)%rcrit3=999.0
        dustact(k)%mact3_mean=0.0

        IF (l_process) THEN
          IF (mad > epsilon(1.0_wp) .AND. nitot > ni_tidy) THEN

            ! Equal partitioning of dust distribution across all ice species
            ! This will cause problems, e.g. when ice splinters, there will be a
            ! spurious increase in the distributed dust mass in the ice category
            ! (This isn't a problem if there is no snow or graupel)

            dustact(k)%nact=nitot
            dustact(k)%mact=mad
            dustact(k)%rcrit=0.0
            dustact(k)%mact_mean=dustact(k)%mact/(dustact(k)%nact+epsilon(mad))
            dustact(k)%mact_mean=min(dust_mact_mean_max, dustact(k)%mact_mean)

            IF (snow_number > epsilon(1.0_wp)) THEN
              dustact(k)%nact2=nitot*ratio_s
              dustact(k)%rcrit2=0.0
              dustact(k)%mact2=mad *ratio_s
              dustact(k)%mact2_mean=dustact(k)%mact2/(dustact(k)%nact2+epsilon(nhtot))
              dustact(k)%mact2_mean=min(dust_mact_mean_max, dustact(k)%mact2_mean)
              dustact(k)%nratio2=max(0.0,min(1.0,dustact(k)%nact2/(snow_number+epsilon(nhtot)) ))
            END IF

            IF (i_ng > 0 .AND. graupel_number > ni_tidy) THEN
              dustact(k)%nact3=nitot*ratio_g
              dustact(k)%rcrit3=0.0
              dustact(k)%mact3=mad *ratio_g
              dustact(k)%mact3_mean=dustact(k)%mact3/(dustact(k)%nact3+epsilon(nhtot))
              dustact(k)%mact3_mean=min(dust_mact_mean_max, dustact(k)%mact3_mean)
              dustact(k)%nratio3=max(0.0,min(1.0,dustact(k)%nact3/(graupel_number+epsilon(nhtot)) ))
            END IF
            
            IF (ice_number > epsilon(1.0_wp)) THEN
               dustact(k)%nact1=max(0.0_wp, dustact(k)%nact-dustact(k)%nact2-dustact(k)%nact3)
               dustact(k)%mact1=max(0.0_wp, dustact(k)%mact-dustact(k)%mact2-dustact(k)%mact3)
               dustact(k)%rcrit1=0.0
               dustact(k)%mact1_mean=dustact(k)%mact1/(dustact(k)%nact1+epsilon(mar))
               dustact(k)%mact1_mean=min(dust_mact_mean_max, dustact(k)%mact1_mean)
               dustact(k)%nratio1=max(0.0,min(1.0,dustact(k)%nact1/(ice_number+epsilon(nhtot)) ))
            END IF

          END IF
        END IF
      END DO

      IF (present(aeroice)) THEN
        IF (l_process) THEN
          ! Activated aerosol in ice
          DO k=1, nz
            density=aerochem(k)%density(1)
            maai=aerofields(k,i_am8)
            mac=aerofields(k,i_am4)
            ratio_ail = maai/(maai+mac+epsilon(1.0_wp))

            ice_number=qfields(k,i_ni)
            snow_number=qfields(k,i_ns)
            ice_mass=qfields(k,i_qi)
            snow_mass=qfields(k,i_qs)
            graupel_mass=qfields(k,i_qg)

            nhtot=ice_number+snow_number
            mitot=ice_mass+snow_mass
            IF (i_ng > 0) THEN
              graupel_number=qfields(k,i_ng)
              nhtot=nhtot+graupel_number
              mitot=mitot+graupel_mass
            END IF

            IF (l_passivenumbers) THEN
              nitot=aerofields(k,i_an11)*ratio_ail
            ELSE
              nitot=nhtot
            END IF

            nhtot=nhtot+epsilon(1.0_wp) ! prevent possible later division by zero
            nratio_i=ice_number/nhtot
            nratio_s=snow_number/nhtot
            nratio_g=graupel_number/nhtot
            !            mratio_i=ice_mass/mitot
            !            mratio_s=snow_mass/mitot
            !            mratio_g=graupel_mass/mitot

            ! Use nratios...
            ratio_i=nratio_i
            ratio_s=nratio_s
            ratio_g=nratio_g

            ! Initialize to zero/defaults
            aeroice(k)%nact=0.0
            aeroice(k)%mact=0.0
            aeroice(k)%rcrit=999.0
            aeroice(k)%mact_mean=0.0
            aeroice(k)%nact1=0.0
            aeroice(k)%nratio1=0.0
            aeroice(k)%mact1=0.0
            aeroice(k)%rcrit1=999.0
            aeroice(k)%mact1_mean=0.0
            aeroice(k)%mact2=0.0
            aeroice(k)%nact2=0.0
            aeroice(k)%nratio2=0.0
            aeroice(k)%rcrit2=999.0
            aeroice(k)%mact2_mean=0.0
            aeroice(k)%mact3=0.0
            aeroice(k)%nact3=0.0
            aeroice(k)%nratio3=0.0
            aeroice(k)%rcrit3=999.0
            aeroice(k)%mact3_mean=0.0

            IF (maai > epsilon(1.0_wp) .AND. nitot > ni_tidy) THEN
              ! Equal partitioning of dust distribution across all ice species
              aeroice(k)%nact=nitot
              aeroice(k)%mact=maai
              aeroice(k)%rcrit=0.0
              aeroice(k)%mact_mean=aeroice(k)%mact/(aeroice(k)%nact+epsilon(nhtot))
              aeroice(k)%mact_mean=min(aero_mact_mean_max, aeroice(k)%mact_mean)

              IF (ratio_s > epsilon(1.0_wp)) THEN
                aeroice(k)%nact2=nitot*ratio_s
                aeroice(k)%rcrit2=0.0
                aeroice(k)%mact2=maai*ratio_s
                aeroice(k)%mact2_mean=aeroice(k)%mact2/(aeroice(k)%nact2+epsilon(nhtot))
                aeroice(k)%mact2_mean=min(aero_mact_mean_max, aeroice(k)%mact2_mean)
                aeroice(k)%nratio2=max(0.0, min(1.0,aeroice(k)%nact2/(snow_number+epsilon(nhtot)) ))
              END IF

              IF (ratio_g > epsilon(1.0_wp)) THEN
                aeroice(k)%nact3=nitot*ratio_g
                aeroice(k)%rcrit3=0.0
                aeroice(k)%mact3=maai*ratio_g
                aeroice(k)%mact3_mean=aeroice(k)%mact3/(aeroice(k)%nact3+epsilon(nhtot))
                aeroice(k)%mact3_mean=min(aero_mact_mean_max, aeroice(k)%mact3_mean)
                aeroice(k)%nratio3=max(0.0, min(1.0, aeroice(k)%nact3/(graupel_number+epsilon(nhtot)) ))
              END IF
              
              IF (ice_number > epsilon(1.0_wp)) THEN
                aeroice(k)%nact1=max(0.0_wp, aeroice(k)%nact-aeroice(k)%nact2-aeroice(k)%nact3)
                aeroice(k)%mact1=max(0.0_wp, aeroice(k)%mact-aeroice(k)%mact2-aeroice(k)%mact3)
                aeroice(k)%rcrit1=0.0
                aeroice(k)%mact1_mean=aeroice(k)%mact1/(aeroice(k)%nact1+epsilon(mar))
                aeroice(k)%mact1_mean=min(aero_mact_mean_max, aeroice(k)%mact1_mean)
                aeroice(k)%nratio1=max(0.0, min(1.0,aeroice(k)%nact1/(ice_number+epsilon(nhtot)) ))
              END IF

            END IF
          END DO
        END IF
      END IF

      IF (present(dustliq)) THEN
        IF (l_process) THEN
          ! Activated aerosol in ice
          DO k=1, nz
            density=dustchem(k)%density(1)
            madl=aerofields(k,i_am9)
            mad=aerofields(k,i_am7)
            ratio_dli = madl/(mad+madl+epsilon(1.0_wp))

            cloud_number=qfields(k,i_nl)
            rain_number=qfields(k,i_nr)

            nratio_l=min(1.0,max(0.0,cloud_number/(cloud_number+rain_number+epsilon(cloud_number))))
            nratio_r=max(0.0,min(1.0,1.0-nratio_l))

            ! Use nratios...
            ratio_l=nratio_l
            ratio_r=nratio_r

            IF (l_passivenumbers_ice) THEN
              ntot=aerofields(k,i_an12)*ratio_dli
            ELSE
              ntot=cloud_number + rain_number
            END IF

            ! Initialize to zero/defaults
            dustliq(k)%nact=0.0
            dustliq(k)%mact=0.0
            dustliq(k)%rcrit=999.0
            dustliq(k)%mact_mean=0.0
            dustliq(k)%nact1=0.0
            dustliq(k)%nratio1=0.0
            dustliq(k)%mact1=0.0
            dustliq(k)%rcrit1=999.0
            dustliq(k)%mact1_mean=0.0
            dustliq(k)%mact2=0.0
            dustliq(k)%nact2=0.0
            dustliq(k)%nratio2=0.0
            dustliq(k)%rcrit2=999.0
            dustliq(k)%mact2_mean=0.0
            dustliq(k)%mact3=0.0
            dustliq(k)%nact3=0.0
            dustliq(k)%nratio3=0.0
            dustliq(k)%rcrit3=999.0
            dustliq(k)%mact3_mean=0.0

            IF (madl > epsilon(1.0_wp) .AND. ntot > nr_tidy) THEN
              ! Equal partitioning of aerosol distribution across all liquid species

              dustliq(k)%nact=ntot
              dustliq(k)%mact=madl
              dustliq(k)%rcrit=0.0
              dustliq(k)%mact_mean=dustliq(k)%mact/(dustliq(k)%nact+epsilon(nhtot))
              dustliq(k)%mact_mean=min(dust_mact_mean_max, dustliq(k)%mact_mean)

              IF (ratio_r > epsilon(1.0_wp)) THEN
                dustliq(k)%nact2=ntot*ratio_r
                dustliq(k)%rcrit2=0.0
                dustliq(k)%mact2=madl*ratio_r
                dustliq(k)%mact2_mean=dustliq(k)%mact2/(dustliq(k)%nact2+epsilon(nhtot))
                dustliq(k)%mact2_mean=min(dust_mact_mean_max, dustliq(k)%mact2_mean)
                dustliq(k)%nratio2=max(0.0,min(1.0,dustliq(k)%nact2/(rain_number+epsilon(nhtot)) ))
              END IF

              IF (ratio_l > epsilon(1.0_wp)) THEN
                dustliq(k)%nact1=max(0.0_wp, dustliq(k)%nact-dustliq(k)%nact2)
                dustliq(k)%mact1=max(0.0_wp, dustliq(k)%mact-dustliq(k)%mact2)
                dustliq(k)%rcrit1=0.0
                dustliq(k)%mact1_mean=dustliq(k)%mact1/(dustliq(k)%nact1+epsilon(mar))
                dustliq(k)%mact1_mean=min(dust_mact_mean_max, dustliq(k)%mact1_mean)
                dustliq(k)%nratio1=max(0.0,min(1.0, dustliq(k)%nact1/(cloud_number+epsilon(nhtot)) ))
              END IF

            END IF
          END DO
        END IF
      END IF
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE examine_aerosol

  ! Note: despite appearences this is only coded up
  ! for a single mode so far...
  SUBROUTINE AbdulRazzakGhan2000(w, p, T, phys, chem, nccn, Smax, active_phys, nccn_active, l_useactive )

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN) :: w ! vertical velocity (ms-1)
    REAL(wp), INTENT(IN) :: p ! pressure (Pa)
    REAL(wp), INTENT(IN) :: T ! temperature (K)
    TYPE(aerosol_phys), INTENT(IN) :: phys
    TYPE(aerosol_chem), INTENT(IN) :: chem
    TYPE(aerosol_active), INTENT(IN), OPTIONAL :: active_phys
    REAL(wp), INTENT(OUT) :: nccn_active ! notional nccn that would be generated by activated aerosol
    REAL(wp), INTENT(OUT) :: Nccn(:) ! number of activated aerosol for each mode
    REAL(wp), INTENT(OUT) :: smax ! peak supersaturation
    LOGICAL, INTENT(OUT) :: l_useactive ! Do we use consider already activated aerosol

    ! Local variables

    REAL(wp), ALLOCATABLE :: s_cr(:)
    REAL(wp) :: s_cr_active
    REAL(wp) :: Tc, Ak, Bk, eta, alpha, gamma_var, es, bigG
    REAL(wp) :: f1, f2, zeta, error_func
    REAL(wp) :: rsmax2  ! 1/(smax*smax)
    REAL(wp) :: diff
    INTEGER :: i
    REAL(wp) :: total

    CHARACTER(len=*), PARAMETER :: RoutineName='ABDULRAZZAKGHAN2000'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_useactive=.FALSE.
    IF (present(active_phys)) l_useactive=active_phys%nact > ccn_tidy .AND. l_active_inarg2000

    nccn_active=0.0

    ALLOCATE(s_cr(phys%nmodes))

    Tc=T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*TC)) * 1e-3

    Ak=2.0*Mw*zetasa/(Ru*T*rhow)
    alpha=g*(Lv/(mp_eps*cp*T)-1)/(T*Rd)
    es=(100.0*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))
    gamma_var=mp_eps*p/es+Lv**2/(Rv*T**2*cp)
    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4
    bigG=1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))
    zeta=(2.0/3.0)*Ak*(w*alpha/BigG)**0.5
    rsmax2=0.0

    DO i=1, phys%nmodes
      IF (phys%N(i) > ccn_tidy) THEN
        IF(l_ukca_casim) THEN
          Bk =chem%bk(i)
        ELSE
          Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        END IF
        s_cr(i)=(2.0/sqrt(Bk))*(Ak/(3.0*phys%rd(i)))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*phys%N(i))
        f1=0.5*exp(2.5*(log(phys%sigma(i)))**2)
        f2=1.0+0.25*log(phys%sigma(i))
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr(i)*s_cr(i)/(eta+3.0*zeta))**.75)/(s_cr(i)*s_cr(i))
      END IF
    END DO

    IF (l_useactive) THEN  ! We should include all the activated aerosol
      IF (active_phys%nact > ccn_tidy) THEN
        i=aero_index%i_accum ! use accumulation mode chem (should do this properly)
        IF(l_ukca_casim) THEN
          Bk =chem%bk(i)
        ELSE
          Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        END IF
        s_cr_active=(2.0/sqrt(Bk))*(Ak/(3.0*active_phys%rd))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*active_phys%nact)
        f1=0.5*exp(2.5*(log(active_phys%sigma))**2)
        f2=1.0+0.25*log(active_phys%sigma)
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr_active*s_cr_active/(eta+3.0*zeta))**.75)/(s_cr_active*s_cr_active)
      END IF
    END IF
    smax=sqrt(1.0/rsmax2)
    nccn=0.0

    IF (l_useactive) THEN  ! We should include all the activated aerosol
      IF (active_phys%nact > ccn_tidy) THEN
        error_func=1.0-erf(2.0*log(s_cr_active/smax)/(3.0*sqrt(2.0)*log(active_phys%sigma)))
        nccn_active=0.5*active_phys%nact*error_func
      END IF
    END IF

    DO i=1, phys%nmodes
      IF (phys%N(i) > ccn_tidy) THEN
        error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(phys%sigma(i))))
        nccn(i)=0.5*phys%N(i)*error_func

        ! Make sure we don't activate too many...
        nccn(i)=min(nccn(i), .999*phys%N(i))
        ! And don't bother if we don't do much
        IF (nccn(i) < ccn_tidy) nccn(i)=0.0
      END IF
    END DO

    IF (l_useactive) THEN  ! Remove already activated aerosol
      ! from newly activated
      ! Start from smallest mode
      diff =active_phys%nact-nccn_active
      total=0.0
      DO i=1, size(nccn)
        total=nccn(i)+total
        nccn(i)=min(nccn(i), max(0.0_wp, total - diff))
      END DO
      !nccn(1) = max(0.0, nccn(1) - diff)
      !nccn(2) = min(nccn(2), max(0.0, nccn(2) + nccn(1) - diff))
      !nccn(3) = min(nccn(3), max(0.0, nccn(3) + nccn(2) + nccn(1) - diff))
    END IF
    nccn=.99*nccn ! Don't allow all aerosol to be removed.
    DEALLOCATE(s_cr)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE AbdulRazzakGhan2000

  SUBROUTINE AbdulRazzakGhan2000_dust(w, p, T, phys, chem, nccn, Smax, active_phys, nccn_active, &
      nccn_dactive, dphys, dchem, active_dphys, dnccn, l_useactive )

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ABDULRAZZAKGHAN2000_DUST'

    REAL(wp), INTENT(IN) :: w ! vertical velocity (ms-1)
    REAL(wp), INTENT(IN) :: p ! pressure (Pa)
    REAL(wp), INTENT(IN) :: T ! temperature (K)
    TYPE(aerosol_phys), INTENT(IN) :: phys
    TYPE(aerosol_chem), INTENT(IN) :: chem
    TYPE(aerosol_active), INTENT(IN), OPTIONAL :: active_phys
    REAL(wp), INTENT(OUT) :: nccn_active ! notional nccn that would be generated by activated aerosol
    REAL(wp), INTENT(OUT) :: Nccn(:) ! number of activated aerosol for each mode
    REAL(wp), INTENT(OUT) :: smax ! peak supersaturation
    LOGICAL, INTENT(OUT) :: l_useactive ! Do we use consider already activated aerosol

    REAL(wp), ALLOCATABLE :: s_cr(:)
    REAL(wp) :: s_cr_active
    REAL(wp) :: Tc, Ak, Bk, eta, alpha, gamma_var, es, bigG
    REAL(wp) :: f1, f2, zeta, error_func
    REAL(wp) :: rsmax2  ! 1/(smax*smax)
    REAL(wp) :: diff
    INTEGER :: i
    REAL(wp) :: total

    TYPE(aerosol_phys), INTENT(IN) :: dphys
    TYPE(aerosol_chem), INTENT(IN) :: dchem
    TYPE(aerosol_active), INTENT(IN), OPTIONAL :: active_dphys
    REAL(wp), INTENT(OUT) :: nccn_dactive ! notional nccn that would be generated by activated aerosol
    REAL(wp), INTENT(OUT) :: dnccn(:) ! number of activated aerosol
    REAL(wp) :: s_cr_active_d

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_useactive=.FALSE.
    IF (present(active_phys)) l_useactive=active_phys%nact > ccn_tidy .AND. l_active_inarg2000

    nccn_active=0.0

    ALLOCATE(s_cr(phys%nmodes+dphys%nmodes))

    Tc=T-273.15 ! Temperature (C)

    ! surface tension of water, not solute (units are N/m)
    !zetasa = (76.1-(0.155*TC)) * 1e-3

    Ak=2.0*Mw*zetasa/(Ru*T*rhow)
    alpha=g*(Lv/(mp_eps*cp*T)-1)/(T*Rd)
    es=(100.0*6.1121)*exp((18.678-Tc/(234.5))*Tc/(257.14+Tc))
    gamma_var=mp_eps*p/es+Lv**2/(Rv*T**2*cp)
    !Dv = (0.211*((T/273.15)**(1.94))*((100000.)/(p))) / 1.e4
    bigG=1.0/(rhow*(Rv*T/(es*Dv)+Lv*(Lv/(Rv*T)-1)/(ka*T)))
    zeta=(2.0/3.0)*Ak*(w*alpha/BigG)**0.5
    rsmax2=0.0

    DO i=1, phys%nmodes
      IF (phys%N(i) > ccn_tidy) THEN
        IF(l_ukca_casim) THEN
          Bk =chem%bk(i)
        ELSE
          Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        END IF
        s_cr(i)=(2.0/sqrt(Bk))*(Ak/(3.0*phys%rd(i)))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*phys%N(i))
        f1=0.5*exp(2.5*(log(phys%sigma(i)))**2)
        f2=1.0+0.25*log(phys%sigma(i))
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr(i)*s_cr(i)/(eta+3.0*zeta))**.75)/(s_cr(i)*s_cr(i))
      END IF
    END DO

    DO i=1,dphys%nmodes
      IF (dphys%N(i) > ni_tidy) THEN
        Bk = dchem%vantHoff(i)*Mw*dchem%density(i)/(dchem%massMole(i)*rhow)
        s_cr(phys%nmodes+i) = (2.0/SQRT(Bk))*(Ak/(3.0*dphys%rd(i)))**1.5
        eta = (w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*dphys%N(i))
        f1 = 0.5*EXP(2.5*(LOG(dphys%sigma(i)))**2)
        f2 = 1.0 + 0.25*LOG(dphys%sigma(i))
        rsmax2 = rsmax2+(f1*(zeta/eta)**1.5+ f2*(s_cr(phys%nmodes+i)*s_cr(phys%nmodes+i)/(eta+3.0*zeta))**.75)/ &
             (s_cr(phys%nmodes+i)*s_cr(phys%nmodes+i))
      END IF
    END DO

    IF (l_useactive) THEN  ! We should include all the activated aerosol
      IF (active_phys%nact > ccn_tidy) THEN
        i=aero_index%i_accum ! use accumulation mode chem (should do this properly)
        IF(l_ukca_casim) THEN
          Bk =chem%bk(i)
        ELSE
          Bk=chem%vantHoff(i)*Mw*chem%density(i)/(chem%massMole(i)*rhow)
        END IF
        s_cr_active=(2.0/sqrt(Bk))*(Ak/(3.0*active_phys%rd))**1.5
        eta=(w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*active_phys%nact)
        f1=0.5*exp(2.5*(log(active_phys%sigma))**2)
        f2=1.0+0.25*log(active_phys%sigma)
        rsmax2=rsmax2+(f1*(zeta/eta)**1.5+f2*(s_cr_active*s_cr_active/(eta+3.0*zeta))**.75)/(s_cr_active*s_cr_active)
      END IF

      IF (active_dphys%nact > ni_tidy) THEN
        i = aero_index%i_coarse_dust
        Bk = dchem%vantHoff(i)*Mw*dchem%density(i)/(dchem%massMole(i)*rhow)
        s_cr_active_d = (2.0/SQRT(Bk))*(Ak/(3.0*active_dphys%rd))**1.5
        eta = (w*alpha/BigG)**1.5/(2.0*pi*rhow*gamma_var*active_dphys%nact)
        f1 = 0.5*EXP(2.5*(LOG(active_dphys%sigma))**2)
        f2 = 1.0+0.25*LOG(active_dphys%sigma)
        rsmax2 = rsmax2 + (f1*(zeta/eta)**1.5     &
           +f2*(s_cr_active_d*s_cr_active_d/(eta+3.0*zeta))**.75)/(s_cr_active_d*s_cr_active_d)
      END IF
    END IF
    smax=sqrt(1.0/rsmax2)
    nccn=0.0
    dnccn = 0.0

    IF (l_useactive) THEN  ! We should include all the activated aerosol
      IF (active_phys%nact > ccn_tidy) THEN
        error_func=1.0-erf(2.0*log(s_cr_active/smax)/(3.0*sqrt(2.0)*log(active_phys%sigma)))
        nccn_active=0.5*active_phys%nact*error_func
      END IF

      IF (active_dphys%nact > ni_tidy) THEN
        error_func = 1.0-erf(2.0*LOG(s_cr_active_d/smax)/(3.0*SQRT(2.0)*LOG(active_dphys%sigma)))
        nccn_dactive = 0.5*active_dphys%nact*error_func
      END IF
    END IF

    DO i=1, phys%nmodes
      IF (phys%N(i) > ccn_tidy) THEN
        error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(phys%sigma(i))))
        nccn(i)=0.5*phys%N(i)*error_func

        ! Make sure we don't activate too many...
        nccn(i)=min(nccn(i), .999*phys%N(i))
        ! And don't bother if we don't do much
        IF (nccn(i) < ccn_tidy) nccn(i)=0.0
      END IF
    END DO

    DO i=1,dphys%nmodes
      IF (dphys%N(i) > ni_tidy) THEN
        error_func = 1.0 - erf(2.0*LOG(s_cr(phys%nmodes+i)/smax)/(3.0*SQRT(2.0)*LOG(dphys%sigma(i))))
        dnccn(i) = 0.5*dphys%N(i)*error_func
        dnccn(i) = MIN(dnccn(i),.999*dphys%N(i))
        IF (dnccn(i)<ni_tidy) dnccn(i)=0.0
      END IF
    END DO

    IF (l_useactive) THEN  ! Remove already activated aerosol
      ! from newly activated
      ! Start from smallest mode
      diff =active_phys%nact-nccn_active
      total=0.0
      DO i=1, size(nccn)
        total=nccn(i)+total
        nccn(i)=min(nccn(i), max(0.0_wp, total - diff))
      END DO
      !nccn(1) = max(0.0, nccn(1) - diff)
      !nccn(2) = min(nccn(2), max(0.0, nccn(2) + nccn(1) - diff))
      !nccn(3) = min(nccn(3), max(0.0, nccn(3) + nccn(2) + nccn(1) - diff))

      diff = active_dphys%nact - nccn_dactive
      total = 0.0
      DO i=1,SIZE(dnccn)
        total = dnccn(i) + total
        dnccn(i) = MIN(dnccn(i),MAX(0.0_wp, total-diff))
      END DO
    END IF
    nccn=.99*nccn ! Don't allow all aerosol to be removed.
    dnccn=.99*dnccn
    DEALLOCATE(s_cr)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
 
 END SUBROUTINE AbdulRazzakGhan2000_dust

  !
  ! Calculate the moments of a lognormal distribution
  !
  ! int_0^infinty r^p*n(r)dr
  !
  ! where
  !
  !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
  !            *exp(-.5*(log(r/rm)/log(sigma))^2)
  !
  ! rm = median radius
  ! log(rm) = mean of log(n(r)/N)
  ! log(sigma) = s.d. of log(n(r)/N)
  !
  FUNCTION moment_logn(N, rm, sigma, p)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN) :: N, rm, sigma
    REAL(wp), INTENT(IN) :: p ! calculate pth moment


    ! Local variables
    REAL(wp) :: moment_logn
    CHARACTER(len=*), PARAMETER :: RoutineName='MOMENT_LOGN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle
    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    moment_logn=N*rm**p*exp(.5*p*p*log(sigma)**2)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION moment_logn

  !
  ! Calculate the upper partial moment of a
  ! lognormal distribution
  !
  ! int_rcrit^infinty r^p*n(r)dr
  !
  ! where
  !
  !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
  !            *exp(-.5*(log(r/rm)/log(sigma))^2)
  !
  ! rm = median radius
  ! log(rm) = mean of log(n(r)/N)
  ! log(sigma) = s.d. of log(n(r)/N)
  !
  FUNCTION upperpartial_moment_logn(N, rm, sigma, p, rcrit)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN) :: N, rm, sigma
    REAL(wp), INTENT(IN) :: p ! calculate pth moment
    REAL(wp), INTENT(IN) :: rcrit ! lower threshold for partial moment
    REAL(wp) :: upperpartial_moment_logn

    ! Local variables

    CHARACTER(len=*), PARAMETER :: RoutineName='UPPERPARTIAL_MOMENT_LOGN'
    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    IF (rcrit==0.0) THEN
      upperpartial_moment_logn=moment_logn(N, rm, sigma, p)
    ELSE
      upperpartial_moment_logn=N*rm**p*exp(.5*p*p*log(sigma)**2)     &
           * .5*casim_erfc((log(rcrit/rm)/log(sigma) - p*log(sigma))/sqrt(2.0))
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION upperpartial_moment_logn

  !
  ! Calculate the lower bound of integral given
  ! a complete and upper partial moment
  !
  ! m = int_0^infinty r^p*n(r)dr
  ! mup = int_x^infinty r^p*n(r)dr
  !
  ! where
  !
  !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
  !            *exp(-.5*(log(r/rm)/log(sigma))^2)
  !
  ! rm = median radius
  ! log(rm) = mean of log(n(r)/N)
  ! log(sigma) = s.d. of log(n(r)/N)
  !
  FUNCTION invert_partial_moment(m, mup, p, rm, sigma)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(wp), INTENT(IN) :: m, mup
    REAL(wp), INTENT(IN) :: p ! pth moment
    REAL(wp), INTENT(IN) :: rm, sigma
    REAL(wp) :: x, invert_partial_moment

    ! Local variables
    CHARACTER(len=*), PARAMETER :: RoutineName='INVERT_PARTIAL_MOMENT'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    x=rm*exp(p*log(sigma)**2+sqrt(2.0)*log(sigma)*erfinv(1.0-2.0*mup/m))
    invert_partial_moment=x

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION invert_partial_moment

  !
  ! Calculate the lower bound of integral given
  ! a complete and upper partial moment
  !
  ! mup = int_x^infinty r^p*n(r)dr
  !
  ! where
  !
  !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
  !            *exp(-.5*(log(r/rm)/log(sigma))^2)
  !
  ! rm = median radius
  ! log(rm) = mean of log(n(r)/N)
  ! log(sigma) = s.d. of log(n(r)/N)
  !
  ! This uses an approximation for erfc(log(x)) to achieve this
  FUNCTION invert_partial_moment_approx(mup, p, rm, sigma)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN) :: mup, rm, sigma, p ! pth moment

    ! Local variables
    REAL(wp) :: invert_partial_moment_approx
    REAL(wp) :: beta, c, lsig, mbeta
    REAL(wp) :: x

    CHARACTER(len=*), PARAMETER :: RoutineName='INVERT_PARTIAL_MOMENT_APPROX'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    c=-0.1021 ! constant in approximation
    lsig=sqrt(2.0)*log(sigma)

    beta=mup*exp(-0.25*p*p*lsig*lsig)/rm**p
    mbeta=1.0-beta

    IF (beta < 1.0e-4) THEN
      x=rm *(-c)**(-lsig)
    ELSE IF (beta >.9999) THEN
      x=0.0
    ELSE
      x=rm * (sigma**(-0.5*p)*(((mbeta)*c+sqrt((mbeta)**2*c*c + 4.0*beta*(mbeta)))/(2.0*beta)))**(lsig)
    END IF
    invert_partial_moment_approx=x

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION invert_partial_moment_approx

  !
  ! Calculate the lower bound of integral given
  ! a complete and upper partial moment
  !
  ! mup = int_x^infinty r^p*n(r)dr
  !
  ! where
  !
  !    n(r) = N/(r*sqrt(2*pi)*log(sigma)) &
  !            *exp(-.5*(log(r/rm)/log(sigma))^2)
  !
  ! rm = median radius
  ! log(rm) = mean of log(n(r)/N)
  ! log(sigma) = s.d. of log(n(r)/N)
  !
  ! This uses an approximation 26.2.23 from A&S
  FUNCTION invert_partial_moment_betterapprox(mup, p, N, rm, sigma)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE


    ! Subroutine arguments
    REAL(wp), INTENT(IN) :: mup, rm, sigma, N, p ! pth moment

    ! Local variables
    REAL(wp) :: invert_partial_moment_betterapprox
    REAL(wp) :: frac, x
    REAL(wp) :: small_frac=1e-6

    CHARACTER(len=*), PARAMETER :: &
    RoutineName='INVERT_PARTIAL_MOMENT_BETTERAPPROX'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (mup==0.0) THEN
      frac=epsilon(x)
    ELSE
      frac = mup/(N*rm**p)*exp(-.5*p*p*log(sigma)*log(sigma))
    END IF

    IF (frac>=1.0 - epsilon(x)) frac=.999
    IF (frac > small_frac) THEN
      x=normal_quantile(frac)
      invert_partial_moment_betterapprox=rm*exp(x*log(sigma)+p*log(sigma)*log(sigma))
    ELSE
      invert_partial_moment_betterapprox=0.0
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION invert_partial_moment_betterapprox

  ! This uses an approximation 26.2.23 from A&S
  FUNCTION normal_quantile(p)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    REAL(wp), INTENT(IN) :: p

    ! Local variables
    REAL(wp) :: normal_quantile
    REAL(wp) :: pp, t,x
    REAL(wp) :: c0,c1,c2,d1,d2,d3
    CHARACTER(len=*), PARAMETER :: RoutineName='NORMAL_QUANTILE'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    c0=2.515517
    c1=0.802853
    c2=0.010328
    d1=1.432788
    d2=0.189269
    d3=0.001308

    pp=p

    IF (p<=epsilon(x)) THEN
      x=huge(x)
    ELSE IF (p >=1.0 - epsilon(x)) THEN
      x=-huge(x)
    ELSE
      IF (p>0.5) pp=1.0-p
      t=sqrt(log(1.0/(pp*pp)))
      x=t-(c0+c1*t+c2*t*t)/(1.0+d1*t+d2*t*t+d3*t*t*t)
      IF (p>0.5)x=-x
    END IF
    normal_quantile=x

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION normal_quantile

END MODULE aerosol_routines
