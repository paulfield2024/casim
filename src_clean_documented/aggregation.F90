! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Rain self-collection (racr) and snow self-aggregation (ice_aggregation) process rates.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Self-collection (aggregation) processes: rain-rain and
!   snow-snow self collection, which change number concentration
!   only (no mass change), per hydrometeor type.
!
! Paper reference:
!   Field et al. (2023) Appendix A.5.3, Eq. (A17)-(A18): rain
!   self-collection follows Beheng (1994); snow self-collection
!   follows the double-moment Large Eddy Model formulation of
!   Gray et al. (2001).
!
MODULE aggregation
  USE variable_precision, ONLY: wp, iwp
  USE passive_fields, ONLY: rho, TdegC
  USE mphys_switches, ONLY: i_qr, i_nr, l_2mr, i_ns
! use mphys_switches, only: i_m3r, l_3mr
  USE process_routines, ONLY: process_rate,  process_name, i_pracr, i_iagg, i_sagg, i_gagg
  USE mphys_parameters, ONLY: d_r, hydro_params
! use mphys_parameters, only: p1, p2, p3, rain_params
  USE mphys_constants, ONLY: rhow, rho0
  USE thresholds, ONLY: qr_small, thresh_small, nr_small
! use m3_incs, only: m3_inc_type2
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
  USE special, ONLY: pi
  USE gauss_casim_micro, ONLY: gaussfunclookup, gaussfunclookup_2d
! use gauss_casim_micro, only: gauss_casim_func

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='AGGREGATION'

  PRIVATE

  PUBLIC racr, ice_aggregation
CONTAINS

  ! Rain-rain self-collection number tendency Nracr,
  ! Field et al. (2023) Eq. (A17), after Beheng (1994).
  SUBROUTINE racr(ixy_inner, dt, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt           !! dt NEEDED for 3rd moment code
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! Local variables
    REAL(wp) :: dnumber, Dr, Eff
    REAL(wp) :: rain_mass
    REAL(wp) :: rain_number
!! variables below NEEDED for 3rd moment code
!   real(wp) :: m1, m2, m3, dm1, dm2, dm3
!   real(wp) :: rain_m3
    LOGICAL :: l_beheng=.TRUE.

    INTEGER :: k ! k index for looping over column

    CHARACTER(len=*), PARAMETER :: RoutineName='RACR'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dnumber=0.0
    IF (l_2mr) THEN
       DO k = 1, ubound(qfields,1)
          rain_mass=qfields(k, i_qr)
          rain_number=qfields(k, i_nr)
          ! if (l_3mr) rain_m3=qfields(k, i_m3r)

          IF (rain_mass > qr_small .AND. rain_number>nr_small) THEN
             IF (l_beheng) THEN
                ! This commented code block imposes a limit on the size of 
                ! agg. This is based on a pragmatic choice that drops large
                ! than 600 microns will breakup. This code effectively negates
                ! excessive numerical size sorting in kinematic cases
                 Dr=(.75/pi)*(rain_mass/rain_number/rhow)**(1.0/d_r)
                 IF ( Dr < 600.0e-6) THEN
                    ! Modified from original
                    Eff=.5
                 ELSE
                    Eff=0.0
                 END IF
                !
                ! For RA3 testing Beheng 1994 Atmos. Res. (equ 10) was used (below)
                ! Beheng uses cgs units, ie. g m-3, CASIM uses SI. Hence, the dnumber
                ! equation below is multiplied by rho to convert from cgs to SI 
                ! (also the 10e3 factor from Beheng is not included)
                !Eff=1.0 !Beheng 1994 Atmos. Res.
                dnumber=Eff * 8.0 * rain_number * rain_mass * rho(k,ixy_inner)   ! #/kg/s
                 !CFL limitation 
                dnumber=rain_number * min( dnumber*dt/rain_number, 0.5  )/dt                            
             END IF
             procs(i_nr, i_pracr%id)%column_data(k)=-dnumber

            ! if (l_3mr) then
            !    m1=rain_mass/rain_params%c_x
            !    m2=rain_number
            !    m3=rain_m3
            !
            !    dm1=0.0
            !    dm2=-Dt*dnumber
            !    call m3_inc_type2(m1, m2, m3, p1, p2, p3, dm1, dm2, dm3)
            !    dm3=dm3/dt
            !    procs(i_m3r, i_pracr%id)%column_data(k) = dm3
            ! end if
          END IF
       END DO
    END IF



    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE racr

  !< Subroutine to determine the aggregation of
  !< ice, snow and graupel.  Aggregation causes a sink in number
  !< and for triple moment species a corresponding change in the
  !< 3rd moment assuming shape parameter is not changed
  !< NB: Aerosol mass is not modified by this process
  !
  !< CODE TIDYING: Move efficiencies into parameters
  ! Snow-snow self-collection number tendency Nsacs,
  ! Field et al. (2023) Eq. (A18), after Gray et al. (2001).
  SUBROUTINE ice_aggregation(ixy_inner, dt, nz, l_Tcold, params, qfields, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim
    USE special, ONLY: Gammafunc

    IMPLICIT NONE

    ! Subroutine arguments
    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dt
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:)   
    TYPE(hydro_params), INTENT(IN) :: params
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)

    ! Local variables
    TYPE(process_name) :: iproc ! processes selected depending on which species we're modifying
    REAL(wp) :: dnumber
    REAL(wp) :: Eff ! collection efficiencies need to re-evaluate these and put them in properly to mphys_parameters
    REAL(wp) :: mass, gaussterm
    REAL(wp) :: n0, lam, mu

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='ICE_AGGREGATION'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k = 1, nz
       IF (l_Tcold(k)) THEN 
          Eff=min(1.0_wp, 0.2*exp(0.08*TdegC(k,ixy_inner)))
          SELECT CASE (params%id)
          CASE (3_iwp) !ice
             iproc=i_iagg
          CASE (4_iwp) !snow
             iproc=i_sagg
             Eff=MIN(1.0_wp, 0.1*exp(0.08*TdegC(k,ixy_inner)))
          CASE (5_iwp) !graupel
             iproc=i_gagg
          END SELECT
          
          mass=qfields(k, params%i_1m)

          IF (mass > thresh_small(params%i_1m) .AND. params%l_2m) THEN ! if no significant ice, we don't bother

             n0=dist_n0(k,params%id)
             mu=dist_mu(k,params%id)
             lam=dist_lambda(k,params%id)

             n0=n0*(lam**(mu+1.0))/(GammaFunc(1.0+mu))

             IF (params%l_3m) THEN
                !gaussterm = gauss_casim_func(mu, params%b_x)
                CALL gaussfunclookup_2d(params%id, gaussterm, mu, params%b_x)
             ELSE
                CALL gaussfunclookup(params%id, gaussterm)
             END IF
             
             dnumber=-1.0*gaussterm * pi*0.125*                                &
                 (rho0/rho(k,ixy_inner))**params%g_x*params%a_x*n0*n0*         &
                 Eff*lam**(-(4.0 + 2.0*mu + params%b_x))
      
      
             dnumber=-1.0 * qfields(k, i_ns) *                                 &
                 min( -dnumber*dt/qfields(k, i_ns), 0.5  )/dt            

      
             procs(params%i_2m, iproc%id)%column_data(k)=dnumber
          END IF
       END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ice_aggregation
END MODULE aggregation
