! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Stochastic-physics (random parameter) settings for CASIM, updated by the parent model and applied to fall-speed/cloud-number parameters.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Module for random parameters scheme
! The parameters are updated in the UM and
! passed through to CASIM via this module
MODULE casim_stph

  USE variable_precision, ONLY: wp

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='CASIM_STPH'

  ! Stochastic physics settings
  LOGICAL :: l_rp2_casim = .FALSE. ! Switch for RP scheme
  REAL(wp) :: snow_a_x_rp = 12.0 ! Snow fallspeed
  REAL(wp) :: ice_a_x_rp = 6000000.0 ! Ice fallspeed
  REAL(wp) :: mpof_casim_rp = 0.5 ! mixed-phase overlap factor
  REAL(wp) :: fixed_cloud_number_rp = 150.0*1.0e6 ! fixed cloud number

END MODULE casim_stph
