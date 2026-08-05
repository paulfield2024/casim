! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Hard-coded placeholder dust/soluble-aerosol number/mass concentrations used as a stand-in aerosol field.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE dust_hack

! No ModuleName required as no subroutines are present

  IMPLICIT NONE
  ! Hacks to use 3rd moments as aerosol fields
  ! These are plumbed through in mphys_casim_um
  !                      0-2km        2-4km          4-6km
  ! dust /cc               0.0466311    0.0506125  0.000694496
  ! soluble /cc              543.755      176.647      185.030
  ! soluble mass kg/m3   1.27056e-09  5.96604e-10  1.01946e-10

  REAL :: n_dust_1=0.0466311*1.0e6
  REAL :: n_dust_2=0.0506125*1.0e6
  REAL :: n_dust_3=0.000694496*1.0e6

  REAL :: n_sol_1=543.755*1.0e6
  REAL :: n_sol_2=176.647*1.0e6
  REAL :: n_sol_3=185.030*1.0e6

  REAL :: m_sol_1=1.27056e-09
  REAL :: m_sol_2=5.96604e-10
  REAL :: m_sol_3=1.01946e-10

  REAL :: dust_factor = 1.0
  REAL :: sol_factor = 1.0

END MODULE dust_hack
