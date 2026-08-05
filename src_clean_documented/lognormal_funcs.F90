! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Inverts lognormal aerosol-mode mass/number to a mean radius (MNtoRm).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Converts between lognormal aerosol-mode moments (mass M,
!   number N) and the modal (number-mean) radius Rm.
!
! Paper reference:
!   Supports the MURK/ARCL aerosol-to-droplet-number coupling
!   described qualitatively in Field et al. (2023) Sec. A.7-A.8;
!   no numbered equation is given there for this conversion.
!
MODULE lognormal_funcs

USE variable_precision, ONLY: wp
USE thresholds, ONLY: ccn_tidy, aeromass_small, aeronumber_small
USE mphys_constants, ONLY: pi

IMPLICIT NONE

CONTAINS

  !
  ! Calculate mean radius of lognormal distribution
  ! given Mass and number
  !
  ! Lognormal mass/number -> modal radius conversion, used
  ! by the MURK/ARCL aerosol coupling (Field et al., 2023,
  ! Sec. A.7-A.8).
  FUNCTION MNtoRm(M, N, density, sigma)   

    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: M, N, density, sigma
    REAL(wp) :: MNtoRm

!    if (N==0 .or. M==0) then ! shouldn't really be here - DPG_bug_checks - M and N can be negative since gets called for increments from evaporation (removal of aerosol from activated mode).
    IF (abs(N)<aeronumber_small .OR. abs(M)<aeromass_small) THEN ! shouldn't really be here. Maybe should make thresholds smaller since will be divided by the timestep in some instances.
      MNtoRm=0.0
    ELSE
      MNtoRm=( 3.0*M*exp(-4.5*log(sigma)**2)/(4.0*N*pi*density) )**(1.0/3.0)
    END IF
  END FUNCTION MNtoRm

END MODULE lognormal_funcs
