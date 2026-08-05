! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Defines identifiers for which parent model (KiD/MONC/UM) is driving CASIM.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Records which host model (UM, MONC, KiD, standalone) is
!   driving CASIM, so that host-specific code paths (e.g. UM
!   cloud-fraction coupling vs CASIM's own saturation-adjustment
!   condensation) can be selected.
!
! Paper reference:
!   Field et al. (2023) Appendix A.6.1 explains that CASIM's own
!   condensation/activation step is disabled when running inside
!   the UM (condensation and cloud fraction are instead supplied
!   by the UM cloud scheme) but active in MONC/KiD; this module
!   implements that host-selection switch.
!
! Module to contain CASIM parent information

MODULE casim_parent_mod

IMPLICIT NONE

! No RoutineName as no subroutines or functions in this module

! Parameters for each model CASIM is driven from
INTEGER, PARAMETER :: parent_unset = 0
INTEGER, PARAMETER :: parent_kid   = 1
INTEGER, PARAMETER :: parent_monc  = 2
INTEGER, PARAMETER :: parent_um    = 3

! Variable containing casim parent - initialised
! as unset for now, but should be set by each
! parent model
INTEGER :: casim_parent = parent_unset

END MODULE casim_parent_mod
