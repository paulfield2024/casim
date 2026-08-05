! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Heterogeneous ice nucleation (inuc), five selectable schemes.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Primary heterogeneous ice nucleation: production of cloud ice
!   number/mass from water vapour, either by nudging to a
!   temperature-dependent concentration or via aerosol/dust-based
!   ice-nucleating-particle parametrisations.
!
! Paper reference:
!   Field et al. (2023) Appendix A.9.1: default nudges cloud ice
!   number to Cooper (1987) at water saturation and T < -8C;
!   aerosol-aware options use DeMott et al. (2010, 2015), Niemand
!   et al. (2012), Atkinson et al. (2013) or Tobo et al. (2013),
!   as summarised by Miltenberger et al. (2020). No single
!   numbered equation is given for these parametrisations.
!
MODULE ice_nucleation
  USE mphys_die, ONLY: throw_mphys_error, bad_values, std_msg
  USE variable_precision, ONLY: wp
  USE passive_fields, ONLY: rho, pressure, w, exner
  USE mphys_switches, ONLY: i_qv, i_ql, i_qi, i_ni, i_th , hydro_complexity, i_am4, i_am6, i_an2, l_2mi, l_2ms, l_2mg, &
       i_am8, i_am9, aerosol_option, i_nl, i_ns, i_ng, iopt_inuc, i_am7, i_an6, i_an12, l_process, l_passivenumbers, &
       l_passivenumbers_ice, active_number, active_ice, isol, iinsol, l_itotsg, contact_efficiency, immersion_efficiency, &
       aero_index, l_prf_cfrac, iopt_act, i_cfl, i_cfi, l_nudge_to_cooper
  USE process_routines, ONLY: process_rate, i_inuc, i_dnuc
  USE mphys_parameters, ONLY: nucleated_ice_mass, cloud_params, ice_params
  USE mphys_constants, ONLY: Ls, cp, pi, m3_to_cm3
  USE qsat_funs, ONLY: qsaturation, qisaturation
  USE thresholds, ONLY: ql_small, w_small, ni_tidy, nl_tidy, cfliq_small
  USE aerosol_routines, ONLY: aerosol_phys, aerosol_chem, aerosol_active

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='ICE_NUCLEATION'

CONTAINS

  !> Currently this routine considers heterogeneous nucleation
  !> notionally as a combination of deposition and/or condensation freezing.
  !> Immersion (i.e. freezing through preeixisting resident IN within a cloud drop) and
  !> Contact freezing (i.e. collision between cloud drop and IN) are not
  !> yet properly concidered.  Such freezing mechanisms should consider the
  !> processing of the aerosol in different ways.
  ! Heterogeneous ice nucleation (Cooper 1987 default, or
  ! aerosol-based schemes), Field et al. (2023) Sec. A.9.1.
  SUBROUTINE inuc(ixy_inner, dt, nz, l_Tcold, qfields, cffields, procs, dustphys, aeroact, dustliq, &
       aerosol_procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:) 
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    REAL(wp), INTENT(IN) :: cffields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! aerosol fields
    TYPE(aerosol_phys), INTENT(IN) :: dustphys(:)
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:)
    TYPE(aerosol_active), INTENT(IN) :: dustliq(:)

    ! optional aerosol fields to be processed
    TYPE(process_rate), INTENT(INOUT), OPTIONAL, TARGET :: aerosol_procs(:,:)


    ! Local variables
    REAL(wp) :: dmass, dnumber, dmad, dmac, dmadl

    ! Liquid water and ice saturation for Meyers equation
    REAL(wp) :: lws_meyers, is_meyers

    !coefficients for Demott parametrization
    REAL(wp) :: a_demott, b_demott, c_demott, d_demott, cf
    REAL(wp) :: Tp01 ! 273.16-Tk (or 0.01 - Tc)
    REAL(wp) :: Tc ! local temperature in C
    REAL(wp) :: Tk ! local temperature in K
    REAL(wp) :: th
    REAL(wp) :: qv
    REAL(wp) :: ice_number
    REAL(wp) :: cloud_number, cloud_mass
    REAL(wp) :: qs, qis, dN_imm, dN_contact, ql
    REAL(wp) :: Si(nz), Sw(nz)
    
    REAL(wp) :: cf_liquid, cf_ice

    ! parameters for Meyers et al (1992)
    ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
    ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708–721
    REAL(wp), PARAMETER :: meyers_a = -0.639 ! Meyers eq 2.4 coeff a
    REAL(wp), PARAMETER :: meyers_b = 0.1296 ! Meyers eq 2.4 coeff b

    ! parameters for Tobo et al. (2013)
    REAL(wp) :: a_tobo, b_tobo, c_tobo, d_tobo

    ! variables for surface site based parameterisations
    REAL(wp) :: n_sites, surf_area
    
    INTEGER :: k

    LOGICAL :: l_condition

    CHARACTER(len=*), PARAMETER :: RoutineName='INUC'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    DO k = 1, nz
      IF (l_Tcold(k)) THEN
        th=qfields(k, i_th)
        Tk=th*exner(k,ixy_inner)
        qs=qsaturation(Tk, pressure(k,ixy_inner)/100.0)
        qis=qisaturation(Tk, pressure(k,ixy_inner)/100.0)
        qv=qfields(k, i_qv)
         
        !!    if (qs==0.0 .or. qis==0.0) then
        !!      write(std_msg, '(A)') 'Error in saturation calculation - qs or qis is zero'
        !!      call throw_mphys_error(bad_values,  ModuleName//':'//RoutineName, std_msg)
        !!    end if
        IF (qs > 0.0) THEN
          Sw(k)=qv/qs - 1.0
        ELSE
          Sw(k)=-999.0
        END IF
        IF (qis > 0.0) THEN
          Si(k)=qv/qis - 1.0
        ELSE
          Si(k)=-999.0
        END IF
      END IF
    END DO

    DO k = 1, nz
       IF (l_Tcold(k)) THEN
          IF (l_prf_cfrac) THEN
             IF (cffields(k,i_cfl) > cfliq_small) THEN
                cf_liquid=cffields(k,i_cfl)
             ELSE
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfi) > cfliq_small) THEN
                cf_ice=cffields(k,i_cfi)
             ELSE
                cf_ice=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
          ELSE
             cf_liquid=1.0
             cf_ice=1.0
          END IF

          qv=qfields(k, i_qv)
          th=qfields(k, i_th)
          
          ql=qfields(k, i_ql)
          Tk=th*exner(k,ixy_inner)
          cloud_mass=qfields(k, i_ql) / cf_liquid
          IF (cloud_params%l_2m) THEN
             cloud_number=qfields(k, i_nl) / cf_liquid
          ELSE
             cloud_number=cloud_params%fix_N0
          END IF
          
          Tc=Tk - 273.15

          ! What's the condition for ice nucleation...?
          ! This condition needs to be consistent with the mechanisms we're
          ! parametrized
          !l_condition=(( Sw(k) >= -1.e-8 .and. TdegC(k) < -8)) .or. Si(k) > 0.25

          SELECT CASE(iopt_inuc)
          CASE default
             l_condition=(( Sw(k) >= -0.001 .AND. Tc < -8 .AND. Tc > -38) .OR. &
                            Si(k) >= 0.08)
          CASE (2)
             ! Meyers, same condition as DeMott
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
          CASE (4)
             ! DeMott Depletion of dust (contact and immersion)
             l_condition=( Sw(k) >= -0.001  .AND. Tc < 0 .AND. Tc > -38)
          CASE (6)
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
          CASE (7)
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
          CASE (8)
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
          CASE (9)
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
          CASE (10)
             l_condition=( cloud_number >= nl_tidy .AND. Tc < 0)
             
          END SELECT
          
          IF (l_condition) THEN
             
             IF (ice_params%l_2m)THEN
                ice_number=qfields(k, i_ni) / cf_ice
             ELSE
                ice_number=1.e3 ! PRAGMATIC SM HACK
             END IF
             
             dN_contact=0.0
             dN_imm=0.0
             SELECT CASE(iopt_inuc)
             CASE default
                ! Cooper
                ! Cooper WA (1986) Ice Initiation in Natural Clouds. Precipitation Enhancement -
                ! A Scientific Challenge. Meteor Monogr, (Am Meteor Soc, Boston, MA), 21, pp 29-32.
                
                IF (ql * cf_liquid > ql_small) THEN  !only make ice when liquid present
                   dN_imm=5.0*exp(-0.304*Tc)/rho(k,ixy_inner)
                   IF (iopt_act == 0) THEN !for fixed number concs adjust the rate
                                             !to nudge back to climatology- can be negative -
                                             !this just represents a nudging incr for cooper
                      IF (l_nudge_to_cooper) THEN 
                        dN_imm =  (dN_imm-ice_number)*0.8 
                      ELSE
                         dN_imm = MAX( dN_imm-ice_number, 0.0 )
                      END IF 
                   END IF
                END IF
             CASE (2)
                ! Meyers
                ! Meyers MP, DeMott PJ, Cotton WR (1992) New primary ice-nucleation
                ! parameterizations in an explicit cloud model. J Appl Meteorol 31:708-721
                lws_meyers = 6.112 * exp(17.62*Tc/(243.12 + Tc))
                is_meyers  = 6.112 * exp(22.46*Tc/(272.62 + Tc))
                dN_imm     = 1.0e3 * exp(meyers_a + meyers_b *(100.0*(lws_meyers/is_meyers-1.0)))/rho(k,ixy_inner)
                dN_imm     = MAX( dN_imm-ice_number, 0.0 )
                ! Applied just for water saturation, deposition freezing ignored
                
             CASE (3)
                ! Fletcher NH (1962) The Physics of Rain Clouds (Cambridge Univ Press, Cambridge, UK)
                dN_imm=0.01*exp(-0.6*Tc)/rho(k,ixy_inner)
             CASE (4)
                ! DeMott Depletion of dust
                ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
                ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
                a_demott=5.94e-5
                b_demott=3.33
                c_demott=0.0264
                d_demott=0.0033
                Tp01=0.01-Tc
                
                IF (dustphys(k)%N(1) > ni_tidy) THEN
                   dN_contact=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                                &
                        (rho(k,ixy_inner) * m3_to_cm3 * contact_efficiency*dustphys(k)%N(1))**(c_demott*Tp01+d_demott)
                   dN_contact=min(.9*dustphys(k)%N(1), dN_contact)
                END IF
                
                IF (dustliq(k)%nact1 > ni_tidy) THEN
                   dN_imm=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                                    &
                        (rho(k,ixy_inner) * m3_to_cm3 * dustliq(k)%nact1)**(c_demott*Tp01+d_demott)
                   dN_imm=immersion_efficiency*dN_imm
                   dN_imm=min(dustliq(k)%nact1, dN_imm)
                END IF
                
             CASE (5)
                dN_imm=0.0
                dN_contact=max(0.0_wp, (dustphys(k)%N(1)-ice_number))
                
             CASE (6)
                ! DeMott Depletion of dust (2015)
                ! 'Integrating laborator and field data to quantify the immersion freezing ice nucleation
                !  activity of mineral dust particles', Atmos. Chem. Phys., 15, 393-409, doi:10.5194/acp-15-393-2015
                a_demott = 0.0
                b_demott = 1.25
                c_demott = 0.46
                d_demott = -11.6
                cf = 1.0     ! cf is the default callibration factor from Demott (eq 2, figures 5 and 6)
                Tp01 = 0.01 - Tc
                
                IF (dustphys(k)%N(1) > ni_tidy) THEN
                   dN_contact=1.0e3/rho(k,ixy_inner)*cf*                                                        &
                        (rho(k,ixy_inner)*m3_to_cm3*contact_efficiency*dustphys(k)%N(1))**(a_demott*(273.16-Tk)+b_demott)*  &
                        exp(c_demott*(273.16-Tk)+d_demott)
                   dN_contact=min(.9*dustphys(k)%N(1), dN_contact)
                END IF
                
                IF (dustliq(k)%nact1 > ni_tidy) THEN
                   dN_imm=1.0e3/rho(k,ixy_inner)*cf*                                                            &
                        (rho(k,ixy_inner)*dustliq(k)%nact1*m3_to_cm3)**(a_demott*(273.16-Tk)+b_demott)*       &
                        exp(c_demott*(273.16-Tk)+d_demott)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                END IF
                
             CASE (7)
                ! Niemand et al. (2012) - using only insoluble in liquid!
                ! 'A particle-surface-area-based parameterization of immersion freezing on desert dust particles',
                ! J. Atmos. Sci., 69, 3077-3092, doi:10.1175/JAS-D-11-0249.1
                IF (dustliq(k)%nact1 > ni_tidy) THEN
                   surf_area = 4*pi*dustliq(k)%nact1*rho(k,ixy_inner)*dustphys(k)%rd(aero_index%i_coarse_dust)**2* &
                        EXP(2*dustphys(k)%sigma(aero_index%i_coarse_dust)**2)      ! m2/m3
                   n_sites = EXP(-0.517*Tc+8.934)    ! 1/m2
                   dN_imm = n_sites*surf_area/rho(k,ixy_inner) ! 1/kg
                   ! AKM: this approximation is only valid for small particles with Sae,j*ns <<1 !
                   ! according to the paper for monodisperse aerosol this is ok for d<3mum and T < -30degC
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                END IF
                
             CASE (8)
                ! Atkinson et al. (2013) - using only insoluble in liquid!
                ! 'The importance of feldspar for ice nucleation by mineral dust in mixed-phase clouds',
                ! Nature, 498, 355-358, doi:10.1038/nature12278
                IF (dustliq(k)%nact1 > ni_tidy) THEN
                   surf_area=0.35*4*pi*dustliq(k)%nact1*rho(k,ixy_inner)*(dustphys(k)%rd(aero_index%i_coarse_dust))**2* &
                        EXP(2*dustphys(k)%sigma(aero_index%i_coarse_dust)**2) ! cm2/m3
                   ! AKM: assuming fraction of K-feldspar in insoluble dust is 0.35
                   n_sites = EXP(-1.038*Tk+275.26) !1/cm2
                   dN_imm = n_sites*surf_area/rho(k,ixy_inner)      ! 1/kg
                   ! AKM: this approximation is only valid for small particles ! (s. comment for case (7))
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                END IF
                
             CASE (9)
                ! Tobo et al. (2013) - using only insoluble in liquid!
                ! 'Biological aerosol particles as a key determinant of ice nuclei populations in a forest
                !  ecosystem', J. Geophys. Res., 118, 10100-10110, doi:10.1002/jgrd.50801
                a_tobo = -0.074
                b_tobo = 3.8
                c_tobo = 0.414
                d_tobo = -9.671
                IF (dustliq(k)%nact1 > ni_tidy) THEN
                   dN_imm=1.0e3/rho(k,ixy_inner)                                                    &
                        *(rho(k,ixy_inner)*m3_to_cm3*dustliq(k)%nact1)**(a_tobo*(273.16-Tk)+b_tobo) &
                        *EXP(c_tobo*(273.16-Tk)+d_tobo)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                END IF
                
             CASE (10)
                ! DeMott Depletion of dust - not distinguishing between insoluble in
                ! liquid and interstitial aerosol
                ! 'Predicting global atmospheric ice nuclei distributions and their impacts on climate',
                ! Proc. Natnl. Acad. Sci., 107 (25), 11217-11222, 2010, doi:10.1073/pnas.0910818107
                a_demott = 5.94e-5
                b_demott = 3.33
                c_demott = 0.0264
                d_demott = 0.0033
                Tp01 = 0.01 - Tc
                
                IF ((dustliq(k)%nact1 > ni_tidy) .OR. (dustphys(k)%N(1) > ni_tidy)) THEN
                   dN_imm=1.0e3/rho(k,ixy_inner)*a_demott*(Tp01)**b_demott*                               &
                        (rho(k,ixy_inner)*m3_to_cm3*(dustliq(k)%nact1+dustphys(k)%N(1)))**(c_demott*Tp01+d_demott)
                   dN_imm=MAX(dN_imm-ice_number,0.0)
                   ! distribute INP between interstital and activated dust (for budgeting
                   ! simulations with l_process > 0)
                   IF ((dustliq(k)%nact1 > ni_tidy) .AND. (dustphys(k)%N(1) > ni_tidy)) THEN
                      dN_contact = dN_imm - dustliq(k)%nact1
                      dN_imm = dN_imm - dN_contact
                      dN_contact=MIN(0.9*dustphys(k)%N(1), dN_contact)
                   ELSE IF (dustliq(k)%nact1 > ni_tidy) THEN
                      dN_imm=MIN(dustliq(k)%nact1, dN_imm)
                      dN_contact=0.0
                   ELSE IF (dustphys(k)%N(1) > ni_tidy) THEN
                      dN_contact=MIN(0.9*dustphys(k)%N(1), dN_imm)
                      dN_imm=0.0
                   END IF
                END IF
                
             END SELECT
             
             IF (cloud_params%l_2m) dN_imm=min(dN_imm, cloud_number)
             
             dN_imm=dN_imm/dt
             dN_contact=dN_contact/dt
             dnumber=dN_imm + dN_contact
             dnumber=min(dnumber, cloud_number/dt)  !this is limit for condensation/imm fzg
             
             !convert back to gridbox mean
             dnumber=dnumber*cf_liquid
             cloud_number=cloud_number*cf_liquid
             cloud_mass=cloud_mass*cf_liquid
             
             IF (dnumber > ni_tidy) THEN
                dmass=cloud_mass*dnumber/cloud_number
                procs(i_qi, i_inuc%id)%column_data(k)=dmass
                
                IF (l_2mi) THEN
                   procs(i_ni, i_inuc%id)%column_data(k)=dnumber
                END IF
                IF (cloud_params%l_2m) THEN
                   procs(i_nl, i_inuc%id)%column_data(k)=-dnumber
                END IF
                procs(i_ql, i_inuc%id)%column_data(k)=-dmass
                
                IF (l_process) THEN
                   
                   ! New ice nuclei
                   dmad=dN_contact*dustphys(k)%M(1)/dustphys(k)%N(1)
                   
                   ! Frozen soluble aerosol
                   dmac=dnumber*aeroact(k)%mact1_mean*aeroact(k)%nratio1
                   
                   ! Dust already in the liquid phase
                   dmadl=dN_imm*dustliq(k)%mact1_mean*dustliq(k)%nratio1
                   
                   aerosol_procs(i_am8, i_dnuc%id)%column_data(k)=dmac
                   aerosol_procs(i_am4, i_dnuc%id)%column_data(k)=-dmac
                   aerosol_procs(i_am9, i_dnuc%id)%column_data(k)=-dmadl
                   
                   aerosol_procs(i_am7, i_dnuc%id)%column_data(k)=dmad+dmadl
                   aerosol_procs(i_am6, i_dnuc%id)%column_data(k)=-dmad    ! <WARNING: using coarse mode
                   aerosol_procs(i_an6, i_dnuc%id)%column_data(k)=-dN_contact ! <WARNING: using coarse mode
                   
                   IF (l_passivenumbers_ice) THEN
                      ! we retain information on what'd been nucleated
                      aerosol_procs(i_an12, i_dnuc%id)%column_data(k)=dN_contact
                   END IF
                END IF
             END IF
          END IF
       END IF ! l_Tcold
    END DO  ! k loop

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE inuc
END MODULE ice_nucleation
