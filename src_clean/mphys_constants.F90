! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Physical constants used throughout CASIM (e.g. density of water, reference density).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE mphys_constants
  USE variable_precision, ONLY: wp

  IMPLICIT NONE

  REAL, PARAMETER :: rhow = 997.0

  REAL :: rho0 = 1.22 ! reference density

  REAL :: fixed_cloud_number = 150.0*1.0e6
  REAL :: fixed_rain_number = 5.0e3
  REAL :: fixed_rain_mu = 0.0

  REAL :: fixed_ice_number = 1000.0

  REAL :: fixed_snow_number = 500.0
  REAL :: fixed_snow_mu = 0.0

  REAL :: fixed_graupel_number = 500.0
  REAL :: fixed_graupel_mu = 0.0

  REAL :: fixed_aerosol_number  = 50.0*1.0e6
  REAL :: fixed_aerosol_rm      = 0.05*1.0e-6
  REAL :: fixed_aerosol_sigma   = 1.5
  REAL :: fixed_aerosol_density = 1777.0

  REAL :: cp = 1005.0

  !< Some of these should really be functions of temperature...
  !  AH - constants set to match UM.
  REAL(wp) :: visair = 1.44E-5  ! kinematic viscosity of air
  REAL(wp) :: Cwater = 4180.0  ! specific heat capacity of liquid water
  REAL(wp) :: Cice = 2100.0    ! specific heat capacity of ice

  REAL(wp) :: Mw = 0.18015e-1 ! Molecular weight of water  [kg mol-1].
  REAL(wp) :: zetasa = 0.8e-1 ! Surface tension at solution-air
  ! interface
  REAL(wp) :: Ru = 8.314472   ! Universal gas constant
  REAL(wp) :: Rd = 287.05     ! gas constant for dry air
  REAL(wp) :: Rv = 461.5100164      ! gas constant for water vapour (matched the UM)
  ! Do not use 'eps' below - this is a Fortran intrinsic!
  REAL(wp) :: mp_eps = 1.6077    ! (Rv/Rd)
  REAL(wp) :: repsilon = 0.62198 ! 1.0 / mp_eps (matched the UM)
  REAL(wp) :: Dv = 0.226e-4   ! diffusivity of water vapour in air
  REAL(wp), PARAMETER :: Lv = 0.2501e7   ! Latent heat of vapourization
  REAL(wp), PARAMETER :: Ls = 0.2835e7   ! Latent heat of sublimation
  REAL(wp), PARAMETER :: Lf = Ls - Lv    ! Latent heat of fusion
  REAL(wp) :: ka = 0.243e-1   ! thermal conductivity of air
  REAL(wp) :: g = 9.8         ! gravitational acceleration ms-2
  REAL(wp) :: SIGLV=8.0e-02   ! Liquid water-air surface tension [N m-1].

  REAL(wp), PARAMETER :: pi = 3.14159265358979323846 ! pi (as used in UM)

  REAL(wp) :: m3_to_cm3 = 1.0e-6 ! Conversion factor from [m-3] to [cm-3]

END MODULE mphys_constants
