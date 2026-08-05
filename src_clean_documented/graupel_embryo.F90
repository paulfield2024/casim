! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Small-snow-plus-cloud conversion into graupel embryos (graupel_embryos).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Formation of new graupel embryos from snow riming past a
!   critical rate (currently disabled in the UM configuration).
!
! Paper reference:
!   Field et al. (2023) Appendix A.11.1, Eq. (A44)-(A45), after
!   Reisner et al. (1998); the paper notes this process is
!   currently disabled in this UM implementation.
!
MODULE graupel_embryo
  USE variable_precision, ONLY: wp
  USE process_routines, ONLY: process_rate, i_sacw
  USE passive_fields, ONLY: rho
  USE mphys_parameters, ONLY: snow_params, graupel_params, cloud_params
  USE mphys_constants, ONLY: rho0
  USE thresholds, ONLY: thresh_sig, cfliq_small
  USE special, ONLY: pi, Gammafunc
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0
  USE mphys_switches, ONLY: l_prf_cfrac, i_cfs, i_cfl, mpof
  USE casim_stph, ONLY: l_rp2_casim, mpof_casim_rp, snow_a_x_rp


  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='GRAUPEL_EMBRYO'

  PUBLIC graupel_embryos
CONTAINS

  ! Snow-to-graupel embryo formation, Field et al. (2023)
  ! Eq. (A44)-(A45), after Reisner et al. (1998). Disabled by
  ! default in the UM configuration described in the paper.
  SUBROUTINE graupel_embryos(ixy_inner, dt, nz, l_Tcold, qfields, cffields, procs)

    !< Subroutine to convert some of the small rimed snow to graupel
    !< (Ikawa & Saito 1991)
    !<    !
    !< OPTIMISATION POSSIBILITIES: See gamma functions and distribution calculations
    !<
    !< AEROSOL: NOT DONE YET - internal category transfer

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

    ! Local variables
    REAL(wp) :: cloud_mass, snow_mass, snow_number
    REAL(wp) :: dmass, dnumber
    REAL(wp) :: snow_n0, snow_lam, snow_mu ! distribution parameters

    REAL(wp) :: dnembryo ! rate of embryo creation
    REAL(wp) :: embryo_mass = 1.6e-10 ! mass(kg) of a new graupel embryo (should be in parameters)

    REAL(wp) :: pgsacw  ! Rate of mass transfer to graupel
    REAL(wp) :: Eff

    REAL(wp) :: alpha = 4.0 ! Tunable parameter see Reisner 1998, Ikawa+Saito 1991
    REAL(wp) :: rhogms      ! Difference between graupel density and snow density

    REAL(wp) :: Garg ! Argument for Gamma function
    INTEGER  :: pid  ! Process id

    REAL(wp) :: cf_snow, cf_liquid, overlap_cf

    INTEGER :: k

    CHARACTER(len=*), PARAMETER :: RoutineName='GRAUPEL_EMBRYOS'

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
        snow_params%a_x = snow_a_x_rp
    END IF

    DO k = 1, nz
       IF (l_Tcold(k)) THEN 
          IF (l_prf_cfrac) THEN
             IF (cffields(k,i_cfs) > cfliq_small) THEN
                cf_snow=cffields(k,i_cfs)
             ELSE
                cf_snow=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
             IF (cffields(k,i_cfl) > cfliq_small) THEN
                cf_liquid=cffields(k,i_cfl)
             ELSE
                cf_liquid=cfliq_small !nonzero value - maybe move cf test higher up
             END IF
          ELSE
             cf_snow=1.0
             cf_liquid=1.0
          END IF


          snow_mass=qfields(k, snow_params%i_1m)   / cf_snow
          cloud_mass=qfields(k, cloud_params%i_1m)  / cf_liquid
          
          IF (snow_mass * cf_snow > thresh_sig(snow_params%i_1m) .AND. &
               cloud_mass * cf_liquid > thresh_sig(cloud_params%i_1m)) THEN
             
             IF (snow_params%l_2m) snow_number=qfields(k, snow_params%i_2m) / cf_snow

             snow_n0=dist_n0(k,snow_params%id)
             snow_mu=dist_mu(k,snow_params%id)
             snow_lam=dist_lambda(k,snow_params%id)
             !convert to reisner N0, based on Reisner et al, QJMRS, 1998 
             !( https://doi.org/10.1002/qj.49712454804 )
             snow_n0=snow_n0*(snow_lam**(snow_mu+1.0))/(GammaFunc(1.0+snow_mu))    
             
             !< This efficiency is the same as used in sacw calculation, 
             !< so that we make the assumption that collecting area is approx half 
             !< of circle - similar to operational.
             !< This efficiency should be defined consistently in the parameters
             Eff=0.5_wp
             rhogms=graupel_params%density-snow_params%density  !this is inconsistent for mass~D**2
             
             Garg=2.0+2*snow_params%b_x + snow_mu
            pgsacw = (0.75*alpha*dt*pi/rhogms)*Eff*Eff*rho(k,ixy_inner)*rho(k,ixy_inner)*cloud_mass*cloud_mass   &
           *snow_params%a_x*snow_params%a_x*snow_n0*GammaFunc(Garg)*(2*snow_params%f_x + 2*snow_lam)**(-Garg) &
           !the 2* lambda etc doesnt look the same as in reisner)
           *(rho0/rho(k,ixy_inner))**(2*snow_params%g_x) !in-graupel rate
             
             dnembryo=max(snow_params%density*pgsacw/rhogms/embryo_mass/rho(k,ixy_inner), 0.0_wp)
             dnumber=min(dnembryo, 0.95*snow_number/dt)
             
             !use mixed-phase overlap function
             overlap_cf=min(1.0,max(0.0,mpof*min(cf_snow, cf_liquid)                      &
                  + max(0.0,(1.0-mpof)*(cf_snow+cf_liquid-1.0))))
     
             pid=i_sacw%id
             pgsacw=pgsacw*overlap_cf !convert back to grid box mean
             dmass=procs(snow_params%i_1m, pid)%column_data(k)- pgsacw !grid box mean
             dnumber=dnumber*overlap_cf !convert back to grid box mean
             
             procs(snow_params%i_1m,pid)%column_data(k)=dmass           
             procs(graupel_params%i_1m,pid)%column_data(k)=pgsacw   
             IF (graupel_params%l_2m) THEN
                procs(graupel_params%i_2m,pid)%column_data(k)=dnumber
             END IF
             IF (snow_params%l_2m) THEN
                procs(snow_params%i_2m,pid)%column_data(k)=-dnumber
             END IF

          END IF
       END IF
    END DO
      
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE graupel_embryos
END MODULE graupel_embryo
