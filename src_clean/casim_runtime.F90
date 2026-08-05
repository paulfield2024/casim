! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Module-level runtime state shared across CASIM: current column/level indices and elapsed simulation time.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE casim_runtime

USE variable_precision, ONLY: wp
USE mphys_parameters,   ONLY: zero_real_wp

IMPLICIT NONE
SAVE

INTEGER :: i_here = 0
INTEGER :: j_here = 0
INTEGER :: k_here = 0

! Time in casim in seconds from start. Will be updated with parent model
! timestep as run progresses.
REAL(wp) :: casim_time = zero_real_wp
REAL(wp) :: casim_smax = 1.5 ! A reasonable guess at maximum supersaturation
                             ! overwritten by the KiD namelist
REAL(wp) :: casim_smax_limit_time = 1.0e20 ! A very long time which is 
                                           ! overwritten by the KiD namelist
                                           

END MODULE casim_runtime
