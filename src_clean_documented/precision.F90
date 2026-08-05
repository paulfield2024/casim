! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Determines default real precision (defp) from the host compiler.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Defines additional numeric kind parameters used where a
!   specific (not necessarily working) precision is required.
!
! Paper reference:
!   CASIM software infrastructure; not discussed in Field et al.
!   (2023).
!
MODULE PRECISION

  IMPLICIT NONE

  REAL :: dummy

  INTEGER, PARAMETER :: defp = KIND(dummy) ! default precision

  INTEGER, PARAMETER :: OneByteInt = selected_int_kind(2), TwoByteInt = selected_int_kind(4), &
       FourByteInt = selected_int_kind(9), EightByteInt = selected_int_kind(18)

  INTEGER, PARAMETER :: FourByteReal = selected_real_kind(P =  6, R =  37)   &
       ,EightByteReal = selected_real_kind(P = 13, R =  307)

  INTEGER, PARAMETER :: wp=EightByteReal                  & ! real working precision
       ,iwp=EightByteInt                & ! integer working precision
       ,ncdfp=FourByteReal              & ! netcdf real precision
       ,incdfp=FourByteInt                ! netcdf integer precision

  INTEGER, PARAMETER :: sp=FourByteReal, dp=EightByteReal
END MODULE PRECISION
