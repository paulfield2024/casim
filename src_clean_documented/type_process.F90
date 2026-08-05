! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Declares the process_rate and process_name derived types used for process-rate bookkeeping.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Defines the process_name derived type used to identify
!   which process-rate index a given calculation should
!   accumulate into.
!
! Paper reference:
!   Software infrastructure supporting the process-rate naming
!   convention in Field et al. (2023) Table A2; not itself a
!   numbered equation.
!
! declares derived types used by processes
MODULE type_process
  USE variable_precision, ONLY: wp

  IMPLICIT NONE

  TYPE :: process_rate
     REAL(wp), ALLOCATABLE :: column_data(:)
  END TYPE process_rate

  TYPE :: process_name
     INTEGER :: id          ! Id for array indexing
     INTEGER :: unique_id   ! Unique id for diagnostic identification
     CHARACTER(20) :: p_name  ! Process name
     LOGICAL :: on          ! is the process going to be used, i.e. on=.true.
  END TYPE process_name
END MODULE type_process
