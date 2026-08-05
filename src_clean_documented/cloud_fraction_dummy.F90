! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Dummy cloud-fraction scheme (module cloud_frac_scheme) used so MONC and KiD builds link.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Stand-in cloud-fraction scheme used when CASIM is run
!   without a host cloud-fraction scheme (e.g. standalone/
!   column-model use), providing a trivial 0/1 cloud fraction.
!
! Paper reference:
!   Field et al. (2023) Appendix A.6.1 describes the
!   "all-or-nothing" (cloud fraction of 1 or 0) saturation-
!   adjustment assumption used by CASIM's own condensation when
!   no host cloud-fraction scheme (e.g. the UM bimodal scheme,
!   Van Weverberg et al., 2021) is coupled in.
!
! Dummy module for cloud fraction scheme to ensure that MONC and KiD build

MODULE cloud_frac_scheme

IMPLICIT NONE

CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CLOUD_FRAC_SCHEME'

CONTAINS

SUBROUTINE cloud_frac_casim_mphys(k, pressure, T_in, T_l, rhcrit_lev, &
                                  qs, qv, cloud_mass, qfields, cloud_mass_new)

USE variable_precision, ONLY: wp

USE yomhook,  ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim

IMPLICIT NONE

INTEGER,  INTENT(IN)  :: k
REAL(wp), INTENT(IN)  :: pressure, T_in, T_l, rhcrit_lev, qs, qv, cloud_mass, qfields
REAL(wp), INTENT(OUT) :: cloud_mass_new

CHARACTER(LEN=*),   PARAMETER :: RoutineName='CLOUD_FRAC_CASIM_MPHYS'
INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cloud_mass_new = cloud_mass

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE cloud_frac_casim_mphys

END MODULE cloud_frac_scheme
