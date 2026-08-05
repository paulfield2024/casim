! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Cloud droplet activation from aerosol via a fixed-CCN/simple activation scheme.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Droplet activation from aerosol: generates increments in
!   cloud droplet number (and mass/aerosol number where coupled
!   to an aerosol scheme) from a given supersaturation/updraught.
!
! Paper reference:
!   Field et al. (2023) Appendix A.6.2: droplet number increment
!   may only increase unless cloud fraction decreases, following
!   Stevens et al. (1996); the aerosol-based calculation itself
!   is provided by src/shipway_activation/shipway_activation_mod.F90
!   (Gordon et al., 2020, adapting Abdul-Razzak & Ghan, 2000).
!
MODULE activation
  USE variable_precision, ONLY: wp
  USE mphys_constants, ONLY: fixed_cloud_number
  USE mphys_parameters, ONLY: C1, K1, cloud_params, zero_real_wp
  USE mphys_switches, ONLY: iopt_act, iopt_inuc, aero_index,   &
       l_warm,              &
       iopt_shipway_act,l_ukca_casim,      &
       activate_in_cloud
  USE aerosol_routines, ONLY: aerosol_active, aerosol_phys, aerosol_chem, &
       abdulRazzakGhan2000, upperpartial_moment_logn, &
       invert_partial_moment_betterapprox, &
       AbdulRazzakGhan2000_dust
  USE special, ONLY: pi,GammaFunc
  USE thresholds, ONLY: w_small, nl_tidy, ni_tidy, ccn_tidy, ql_small
  USE shipway_parameters, ONLY: max_nmodes, nmodes, Ndi, &
     rdi, sigmad, bi, betai, use_mode, nd_min
  USE shipway_constants, ONLY: Mw, rhow, eps, Rd, Dv, Lv, cp, &
      Dv_mean, alpha_c, zetasa, Ru
  USE qsat_funs, ONLY: qsaturation,dqwsatdt
  USE shipway_activation_mod, ONLY: solve_nccn_household, solve_nccn_brent
  USE casim_stph, ONLY: l_rp2_casim, fixed_cloud_number_rp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ACTIVATION'

  PRIVATE

  PUBLIC activate

  ! Variables used in the call to the Shipway (2015) activation scheme

  REAL(wp) :: ent_fraction=0.
  INTEGER  :: order=2
  INTEGER  :: niter=8
  REAL(wp) :: smax0=0.001
  !real(wp) :: alpha_c=0.05 !kinetic parameter
  REAL(wp) :: nccni(max_nmodes) ! maximum number of modes that can be used (usually=3)

CONTAINS

  ! Droplet number/mass activation increment, Field et al.
  ! (2023) Sec. A.6.2, after Stevens et al. (1996).
  SUBROUTINE activate(dt, cloud_mass, cloud_number, w, rho, dnumber, dmac, T, p,  &
       cfliq,cfliq_old, aerophys, aerochem, aeroact, dustphys, dustchem,   &
       dustliq, dnccn_all, dmac_all, dnumber_d, dmass_d, dnccnd_all,dmad_all,  &
       smax,ait_cdnc,accum_cdnc, tot_cdnc,activated_arg,activated_cloud)
!PRF NB extra cfliq argument missing from activate call in condensation for non-um


    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE thresholds, ONLY: cfliq_small
    IMPLICIT NONE

    ! Subroutine arguments

    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: cloud_mass, cloud_number, w, rho, T, p, &
         cfliq, cfliq_old
    REAL(wp), INTENT(OUT) ::  dnumber, dmac
    TYPE(aerosol_phys), INTENT(IN) :: aerophys
    TYPE(aerosol_chem), INTENT(IN) :: aerochem
    TYPE(aerosol_active), INTENT(IN) :: aeroact
    TYPE(aerosol_phys), INTENT(IN) :: dustphys
    TYPE(aerosol_chem), INTENT(IN) :: dustchem
    TYPE(aerosol_active), INTENT(IN) :: dustliq
    REAL(wp), INTENT(OUT) :: dnccn_all(:),dmac_all(:)
    REAL(wp), INTENT(OUT) :: dnccnd_all(:),dmad_all(:)
    REAL(wp), INTENT(OUT) :: dnumber_d, dmass_d ! activated dust number and mass
    REAL(wp), INTENT(OUT) :: smax,ait_cdnc, accum_cdnc,tot_cdnc
    REAL(wp), INTENT(OUT) :: activated_cloud, activated_arg
    ! Local Variables

    REAL(wp) :: cloud_number_work, cloud_number_work_old, cloud_mass_work
    REAL(wp) :: cloud_radius_work, smax_cloud, smax_act
    REAL(wp) :: active, rcrit, nccn_active, dactive,nccn_dactive
    REAL(wp) :: Nd, rm, sigma, density
    REAL(wp) :: dnccn_cloud(aero_index%nccn)
    REAL(wp) :: cf_liquid, cf_thresh, cf_liquid_old

    INTEGER :: imode
    LOGICAL :: l_useactive

    REAL(wp) :: LvT, alpha, lam, tau, qs, &
         Dv_here, erfarg
    REAL(wp) :: Ak, bigGthermal, bigGdiffusion, gammaL,     &
         gammaR, gammastar, bigG
    REAL(wp) :: s0i(max_nmodes), sigmas(max_nmodes)
    REAL(wp) :: m1, m2, j1
    REAL(wp) :: kwdqsdz, dqsdt
    REAL(wp), PARAMETER :: smax_act_min = 0.02

    INTEGER, PARAMETER :: solve_household = 1
    INTEGER, PARAMETER :: solve_brent = 2

    REAL(wp) :: dv_flag=0

    INTEGER :: solve_select

    CHARACTER(len=*), PARAMETER :: RoutineName='ACTIVATE'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ! Apply RP scheme
    IF ( l_rp2_casim ) THEN
        fixed_cloud_number = fixed_cloud_number_rp
    END IF

    cf_thresh = cfliq_small ! should be a small number, shouldn't matter too much
    IF (cfliq > cfliq_small) THEN  !only doing liquid cloud fraction at the moment
      cf_liquid=cfliq
    ELSE
      cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
    END IF
    IF (cfliq_old > cfliq_small) THEN  !only doing liquid cloud fraction at the moment
      cf_liquid_old=cfliq_old
    ELSE
      cf_liquid_old=cfliq_small !nonzero value - maybe move cf test higher up
    END IF

!make in-cloud number conc
! Should this be the new cloud fraction or the old cloud fraction?
! I think it should be the cloud fraction after advection but before
! the new condensation/cloud scheme call, to be consistent with the CDNC used
! here.

    cloud_number_work_old = cloud_number / cf_liquid_old
    cloud_number_work = cloud_number / cf_liquid
    ! threshold cloud number conc, so that cloud number
    ! is larger 1 droplet per m^3. This ensures stability of the
    ! tau calculation.
    IF (cloud_number < 1.0e-6) THEN
       cloud_number_work_old = 1.0e-6/cf_liquid_old
       cloud_number_work = 1.0e-6/cf_liquid
    END IF

! This is cloud_mass_preap2 which is the OLD cloud mass.
    cloud_mass_work = cloud_mass / cf_liquid_old
    l_useactive=.FALSE.
    dmac=0.0
    dnumber_d=0.0
    dmass_d=0.0
    dnumber=0.0
    dmac_all=0.0
    dnccn_all=0.0
    dmad_all=0.0
    dnccnd_all=0.0
    dnccn_cloud(:)=0.0
    smax = 0.0
    ait_cdnc=0.0
    accum_cdnc=0.0
    tot_cdnc=0.0
    activated_arg=0.0
    activated_cloud=0.0
    solve_select = solve_brent
    gammastar=0.0
    bigG=0
    IF (int(dv_flag)==1)THEN
       Dv_here=Dv
    ELSE IF(dv_flag <= 0)THEN
       Dv_here=Dv_mean(T,1.0_wp)
    ELSE
      Dv_here=dv_flag
    END IF
    m2=cloud_number_work_old
    m1=cloud_mass_work/cloud_params%c_x
    j1=1.0/(cloud_params%p1-cloud_params%p2)
    qs=(0.029/0.018)*qsaturation(T, p/100.0)
    dqsdt = dqwsatdt (qs,T)

    LvT = Lv -(4217.4-1870.0)*(T-273.15)
    IF (cloud_mass_work > ql_small .AND. cf_liquid > cf_thresh) THEN
      ! The following calculations are based on Gordon et al, 2020,
      ! "Development of aerosol activation in the double-moment
      ! Unified Model and evaluation with CLARIFY measurements",
      ! Atmos. Chem. Phys., 20, 10997-11024,
      ! https://doi.org/10.5194/acp-20-10997-2020, 2020
      lam=((GammaFunc(1.0+cloud_params%fix_mu+cloud_params%p1) / &
           GammaFunc(1.0+cloud_params%fix_mu+cloud_params%p2))*(m2/m1))**(j1)
      cloud_radius_work = 0.5*GammaFunc(2.0+cloud_params%fix_mu)/ &
           (lam*GammaFunc(1.0+cloud_params%fix_mu))
      bigGdiffusion = rhow*Ru*T/(p*qs*Dv_here*0.018)
      bigGthermal = (LvT*rhow/(0.024*T))*(LvT*0.018/(Ru*T) -1)
      bigG = 1.0/(bigGdiffusion+bigGthermal)
      gammaR= 0.018*LvT*LvT/(cp*Ru*T*T)
      gammaL = 0.029/(qs*0.018)
      gammastar = 4*pi*rhow*(gammaL+gammaR)/rho
      tau = 1.0/(gammastar*bigG*cloud_number_work_old*rho*cloud_radius_work)
    ELSE
      tau =1000.0  ! This ensures activation will always happen. Activation is
      ! only called when water is condensing
    END IF
    IF (tau >= 1000.0 .OR. tau <= 0.0) THEN
      tau=1000.0
    END IF

    ! This agrees well with the version in Ghan et al
    alpha = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-ent_fraction)

    smax_cloud= alpha*w*tau

    kwdqsdz = w*5.3e5*15*dqsdt*0.006 ! constant*w*time-threshold*dqsat/dT*dT/dz

    SELECT CASE(iopt_act)
    CASE default
      ! fixed number
      active=fixed_cloud_number
    CASE(1)
      ! activate 100% aerosol
      active=sum(aerophys%N(:))
    CASE(2)
      ! simple Twomey law Cs^k expressed as
      ! a function of w (Rogers and Yau 1989)
      active=0.88*C1**(2.0/(K1+2.0))*(7.0E-2*(w*100.0)**1.5)**(K1/(K1+2.0))*1.0e6/rho
    CASE(3)
      ! Use scheme of Abdul-Razzak and Ghan
      ! setup Shipway parameters to ensure calc_nccn works correctly
      nmodes=aero_index%nccn
      Ak = 2.*Mw*zetasa/(Ru*T*rhow)

      DO imode=1,aero_index%nccn
        IF(l_ukca_casim) THEN
          bi(imode) =aerochem%bk(imode)*aerochem%epsv(imode)
        ELSE
          bi(imode) =aerochem%vantHoff(imode)*aerochem%epsv(imode)* &
             aerochem%density(imode)*Mw/(rhow*aerochem%massMole(imode))
        END IF
        Ndi(imode)=aerophys%N(imode)
        rdi(imode)=aerophys%rd(imode)
        sigmad(imode)=aerophys%sigma(imode)
        betai(imode)=0.5      !aerochem%beta(imode) This is set to 0.5 for
                              !    Shipway not for ARG
        !print *, 'imode, bi', bi(imode), imode, Ak,betai(imode),rdi(imode)
        IF (rdi(imode) > epsilon(1.0_wp)) THEN
           s0i(imode) = rdi(imode)**(-(1.0+betai(imode))) * &
                sqrt(4.0*Ak**3.0/(27.0*bi(imode)))
        ELSE
           s0i(imode) = 0.0
        END IF
        sigmas(imode) = sigmad(imode)**(1.+betai(imode))
        ! only use the mode if there's significant number
        use_mode(imode) = aerophys%N(imode) > Nd_min
      END DO ! loop over modes

      IF (w > w_small .AND. sum(aerophys%N(:)) > ccn_tidy) THEN
         ! The following derivations of nccn is based on an
         ! adaptation of Abdul-Razzak and Ghan scheme, which is
         ! described in Gordon et al, 2020,
         ! "Development of aerosol activation in the double-moment
         ! Unified Model and evaluation with CLARIFY measurements",
         ! Atmos. Chem. Phys., 20, 10997-11024,
         ! https://doi.org/10.5194/acp-20-10997-2020, 2020
        CALL AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax_act, aeroact, &
               nccn_active, l_useactive)

        activated_arg = sum(dnccn_all(:))

        IF (Smax_act > smax_act_min .AND. .NOT. l_warm) THEN

          IF (iopt_inuc < 4) THEN
            ! For lower-order ice nucleation options, need to initialise
            ! dactive and dmass_d
            dactive = zero_real_wp
            dmass_d = zero_real_wp
          ELSE
            ! For higher-order ice nucleation schemes, dactive and dmass_d
            ! can be based on dustphys
            dactive = 0.01*dustphys%N(1)

            IF ( dustphys%N(1) > ni_tidy ) THEN
              dmass_d = dactive * dustphys%M(1) / dustphys%N(1)
            ELSE
              ! Prevent divide by zero generating nonsense.
              dmass_d = zero_real_wp
            END IF ! dustphys%N(1) > ni_tidy

          END IF ! iopt_inuc

        END IF ! Smax_act > 0.02 and .not. l_warm

        IF(activate_in_cloud ==2 .AND. cf_liquid > cf_thresh .AND. &
            smax_cloud <= smax_act) THEN
          ! In this case use a weighted sum of ARG-activation outside cloud and equilibrium inside cloud
          ! should really use the equilibrium*old-cloud-frac+arg*(new-cloud-frac-old-cloud-frac)
          activated_cloud=0
          dnccn_cloud(:)=0
          DO imode=1,aero_index%nccn
             IF (use_mode(imode)) THEN
                  ! In ARG2000, error_func=1.0-erf(2.0*log(s_cr(i)/smax)/(3.0*sqrt(2.0)*log(phys%sigma(i))))
                  ! so insert the factor of 2/3 for consistency
                erfarg=2.0*log(smax_cloud/s0i(imode))/(3.0*sqrt(2.)*log(sigmas(imode)))
                dnccn_cloud(imode) = 0.5*Ndi(imode)*(1.0+erf(erfarg))
                activated_cloud=activated_cloud+dnccn_cloud(imode)
             END IF
          END DO
          ! New CDNC is a sum of CDNC created in old and new cloud
          active = activated_cloud*cf_liquid_old+(cf_liquid-cf_liquid_old)*activated_arg
          active = active/cf_liquid ! This is just because we multiply by CF
          dnccn_all(:) = (dnccn_cloud(:)*cf_liquid_old+(cf_liquid-cf_liquid_old)*dnccn_all(:))/cf_liquid
          ! at the end of the subroutine, and we need to avoid multiplying by
          ! CF twice
        ELSE IF(activate_in_cloud==1 .OR. cf_liquid <= cf_thresh        &
          .OR. smax_cloud > smax_act ) THEN
          ! In this case use Abdul-Razzak & Ghan
          smax = smax_act
          active=sum(dnccn_all(:))
        ELSE
          dnccn_all(:)=0.0
          active=0.0
        END IF ! activate_in_cloud == 2 etc

      ELSE
        active=0.0
        dactive=0.0
      END IF

    CASE(4)
      ! Use scheme of Abdul-Razzak and Ghan (including for insoluble aerosol by
      ! assuming small amount of soluble material on it)
      IF (w > w_small .AND. (sum(aerophys%N(:))+sum(dustphys%N(:))) > ccn_tidy) THEN
        IF (l_warm) THEN
          CALL AbdulRazzakGhan2000(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact, &
             nccn_active, l_useactive)
          active=sum(dnccn_all(:))
        END IF
        IF (.NOT. l_warm) THEN
           CALL AbdulRazzakGhan2000_dust(w, p, T, aerophys, aerochem, dnccn_all, Smax, aeroact, &
                nccn_active, nccn_dactive, dustphys, dustchem, dustliq, dnccnd_all, l_useactive)
           active   = sum(dnccn_all(:))
           dactive  = dnccnd_all(aero_index%i_coarse_dust) !SUM(dnccnd_all(:)) i_accum_dust currently not used!
        END IF
      ELSE
        active=0.0
        dactive=0.0
      END IF

    CASE(iopt_shipway_act)
      ! Use scheme of Shipway 2015
      ! This is a bit clunk and could be harmonized
      nmodes=aero_index%nccn
      DO imode=1,aero_index%nccn
        IF(l_ukca_casim) THEN
          bi(imode) =aerochem%bk(imode)*aerochem%epsv(imode)
        ELSE
          bi(imode) =aerochem%vantHoff(imode)*aerochem%epsv(imode)* &
             aerochem%density(imode)*Mw/(rhow*aerochem%massMole(imode))
        END IF
        Ndi(imode)=aerophys%N(imode)
        rdi(imode)=aerophys%rd(imode)
        sigmad(imode)=aerophys%sigma(imode)
        betai(imode)=aerochem%beta(imode)
        ! only use the mode if there's significant number
        use_mode(imode) = aerophys%N(imode) > Nd_min
      END DO

      IF (any(use_mode) .AND. cloud_number_work < sum(Ndi))THEN
        IF (tau > kwdqsdz .OR. activate_in_cloud==1) THEN
          SELECT CASE (solve_select)
            CASE (solve_household)
              CALL solve_nccn_household( order, niter, smax0, w, T, p, alpha_c, &
                                         ent_fraction, smax, active, nccni     )
            CASE (solve_brent)
              CALL solve_nccn_brent(w, T, p, alpha_c, ent_fraction,      &
                                  smax, active, nccni)
          END SELECT
        ELSE
           LvT = Lv -(4217.4-1870.0)*(T-273.15)
           alpha = 9.8*(LvT/(eps*cp*T)-1.0)/(T*Rd)*(1-ent_fraction)
           smax= alpha*w*tau
           !call calc_nccn(smax,active,nccni)
        END IF

        dnccn_all(1:aero_index%nccn) = nccni(1:aero_index%nccn)
      ELSE
        active=0.0
      END IF
    END SELECT

    IF(activate_in_cloud==2) THEN
      smax = smax_act
      ait_cdnc = tau
      accum_cdnc = smax_cloud
    ELSE IF(activate_in_cloud==0) THEN
      ait_cdnc = tau
      accum_cdnc = dnccn_all(2)
    ELSE
      ait_cdnc = aerophys%N(1)
      accum_cdnc = aerophys%N(2)
    END IF
    tot_cdnc = active

    SELECT CASE(iopt_act)
    CASE default
      ! fixed number, so no need to calculate aerosol changes
      dnumber=max(0.0_wp,(active-cloud_number_work))
    CASE (1:5)
      IF (active > nl_tidy) THEN
        IF (l_useactive) THEN
          dnumber=active
          dnumber_d=dactive
        ELSE
          dnumber=max(0.0_wp,(active+dactive-cloud_number_work))
          ! Rescale to ensure total removal of aerosol number=creation of cloud number
          dnccn_all = dnccn_all*(dnumber/(sum(dnccn_all) + sum(dnccnd_all) + tiny(dnumber)))
          dnccnd_all = dnccnd_all*(dnumber/(sum(dnccn_all) + sum(dnccnd_all) + tiny(dnumber)))
          dnumber = sum(dnccn_all)
          dnumber_d = sum(dnccnd_all)
        END IF
        ! Need to make this consistent with all aerosol_options
        DO imode = 1, aero_index%nccn
          Nd=aerophys%N(imode)
          IF (Nd > ccn_tidy) THEN
            rm=aerophys%rd(imode)
            sigma=aerophys%sigma(imode)
            density=aerochem%density(imode)

            rcrit=invert_partial_moment_betterapprox(dnccn_all(imode), 0.0_wp, Nd, rm, sigma)

            dmac_all(imode)=(4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
            dmac_all(imode)=min(dmac_all(imode),0.999*aerophys%M(imode)) ! Don't remove more than 99.9%
            dmac=dmac+dmac_all(imode)
          END IF
        END DO
        ! for the insoluble mode
        IF (iopt_act == 4) THEN
          dmass_d=zero_real_wp
          dmad_all(:)=zero_real_wp
          DO imode = 1, aero_index%nin
            Nd=dustphys%N(imode)
            IF (Nd > ccn_tidy .AND. dnumber_d > ccn_tidy) THEN
              rm=dustphys%rd(imode)
              sigma=dustphys%sigma(imode)
              density=dustchem%density(imode)

              rcrit=invert_partial_moment_betterapprox(dnccnd_all(imode), 0.0_wp, Nd, rm, sigma)

              dmad_all(imode)=(4.0*pi*density/3.0)*(upperpartial_moment_logn(Nd, rm, sigma, 3.0_wp, rcrit))
              dmad_all(imode)=min(dmad_all(imode),0.999*dustphys%M(imode)) ! Don't remove more than 99.9%
              dmass_d=dmass_d+dmad_all(imode)
            END IF
          END DO
        END IF
      END IF
    END SELECT

    ! Convert to rates rather than increments .... and back to gridbox means
    dmac=dmac/dt * cf_liquid
    dmac_all=dmac_all/dt * cf_liquid
    dnccn_all=dnccn_all/dt * cf_liquid
    dnumber=dnumber/dt * cf_liquid
    dmass_d = dmass_d/dt * cf_liquid
    dmad_all = dmad_all/dt * cf_liquid
    dnccnd_all=dnccnd_all/dt * cf_liquid
    dnumber_d=dnumber_d/dt * cf_liquid

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE activate
END MODULE activation
