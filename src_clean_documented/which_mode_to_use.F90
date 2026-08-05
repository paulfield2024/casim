! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Partitions evaporated aerosol mass/number between CCN modes (which_mode).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Decides, for aerosol processing, whether a droplet/particle
!   should be treated via the simple mass-threshold method or a
!   more complete aerosol-mode-tracking method.
!
! Paper reference:
!   Supports the aerosol-cloud coupling discussed qualitatively
!   in Field et al. (2023) Sec. 1 and Appendix A.6.2-A.8; this
!   internal dispatch logic is not described by a numbered
!   equation in the paper.
!
! Routine decides which modes to put back
! re-evaporated (or potentially any other) aerosol
MODULE which_mode_to_use
  USE variable_precision, ONLY: wp
  USE mphys_switches, ONLY: l_aeroproc_midway
  USE lognormal_funcs, ONLY: MNtoRm
  USE mphys_constants, ONLY: pi

  IMPLICIT NONE

  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='WHICH_MODE_TO_USE'

  INTEGER, PARAMETER :: imethod = 2 ! method to use
  INTEGER, PARAMETER :: iold_method = 1
  INTEGER, PARAMETER :: isimple_method = 2 ! simple method
  REAL(wp), PARAMETER :: r_thresh_fixed = 0.5e-6 ![m] for simple method
  !set threshold to arithmetically half way between accum and coarse sizes

  REAL(wp) :: max_accumulation_mean_radius = 0.25e-6
  REAL(wp) :: min_coarse_mean_radius = 1.0e-6

  PUBLIC which_mode
CONTAINS

  ! Subroutine calculates how much of the total increments
  ! to mass and number should be sent to each of two modes.
  SUBROUTINE which_mode(dm, dn, r1_in, r2_in, density, sigma, dm1, dm2, dn1, dn2)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='WHICH_MODE'

    REAL(wp), INTENT(IN) :: dm  ! total mass increment
    REAL(wp), INTENT(IN) :: dn  ! total number increment

    REAL(wp), INTENT(IN) :: r1_in  ! mean radius of mode 1
    REAL(wp), INTENT(IN) :: r2_in  ! mean radius of mode 2

    REAL(wp), INTENT(IN) :: density ! density of aerosol, assumed the same 
                                    ! across all modes

    REAL(wp), INTENT(IN) :: sigma ! sigma of the aerosol that is being 
                                  ! evaporated/added

    REAL(wp), INTENT(OUT) :: dm1  ! mass increment to mode 1
    REAL(wp), INTENT(OUT) :: dn1  ! number increment to mode 1
    REAL(wp), INTENT(OUT) :: dm2  ! mass increment to mode 2
    REAL(wp), INTENT(OUT) :: dn2  ! number increment to mode 2

    ! local variables
    REAL(wp) :: rm    ! mean radius of increment
    REAL(wp) :: gamma_var ! convenience variable

    REAL(wp) :: r1, r2     ! r1 and r2 (possibly modified) 
    REAL(wp) :: r1_3, r2_3 ! r1**3 and r2**3

    REAL(wp) :: ftpi ! 4/3*pi
    REAL(wp) :: rftpi ! 1./(4/3*pi)

    REAL(wp) :: r_thresh

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dm1=0.0
    dm2=0.0
    dn1=0.0
    dn2=0.0

    IF (dm*dn > 0.0) THEN 
      ! dm and dn should be positive and of the same sign
      r1=min(max_accumulation_mean_radius, r1_in)
      r2=max(min_coarse_mean_radius, r2_in)

      IF (l_aeroproc_midway) THEN
        r_thresh = (r1+r2)/2.0
      ELSE
        r_thresh = r2
      END IF

      ftpi = pi *4./3.
      rftpi=1./ftpi

      IF (l_aeroproc_midway) THEN
        rm = MNtoRm(dm,dn,density,sigma) !DPG_bug_changes
      ELSE
         ! Keep original code to preserve answers when switched off
        rm = (rftpi*dm/dn/density)**(1.0/3.0)
      END IF

      SELECT CASE(imethod)
      CASE default
        IF (rm >= r_thresh) THEN
          dm1=0.0
          dn1=0.0
          dm2=dm
          dn2=dn
        ELSE IF (rm < r1 .OR. l_aeroproc_midway) THEN
          dm1=dm
          dn1=dn
          dm2=0.0
          dn2=0.0
        ELSE
          r1_3=r1*r1*r1
          r2_3=r2*r2*r2
          gamma_var=(rm*rm*rm-r1_3)/(r2_3-r1_3)
          dn2=dn*(gamma_var)
          dn1=dn - dn2
          IF (l_aeroproc_midway) THEN
            dm2=ftpi*density*r2_3*dn2*exp(4.5*log(sigma)**2) !DPG_bug_changes
          ELSE
            dm2=FTPI*density*r2_3*dn2
          END IF

          dm1=dm-dm2
        END IF
      CASE(isimple_method)
        IF (rm >= r_thresh_fixed) THEN
          dm1=0.0
          dn1=0.0
          dm2=dm
          dn2=dn
        ELSE
          dm1=dm
          dn1=dn
          dm2=0.0
          dn2=0.0
        END IF
      END SELECT
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE which_mode
END MODULE which_mode_to_use
