! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Aerosol-mode parameters (number, radius, solubility) used by the Shipway activation scheme.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Aerosol-mode parameters (number concentration, dry radius,
!   geometric standard deviation, solubility/hygroscopicity per
!   mode) used as input to the Shipway activation scheme.
!
! Paper reference:
!   Field et al. (2023) Sec. A.7 describes, for the MURK
!   coupling, a lognormal ammonium-sulphate mode with mode size
!   9.5e-8 m, geometric standard deviation 1.4, density
!   1769 kg/m3 and hygroscopicity B=0.4 supplied via parameters
!   of this kind.
!
MODULE shipway_parameters
  USE variable_precision, ONLY: wp

  IMPLICIT NONE

  !-------------------------------------------------------------
  ! The following variables describe the different 
  ! modes of dry aerosol
  !------------------------------------------------------------
  INTEGER, PARAMETER :: max_nmodes=3
  INTEGER :: nmodes ! number of modes in aerosol distribution.
                    ! This is the size for the following arrays
  INTEGER :: imode  ! current mode being considered
  REAL(wp) :: Ndi(max_nmodes)      ! Number concentration (m-3)
  REAL(wp) :: rdi(max_nmodes)      ! Geometric mean radius (m)
  REAL(wp) :: sigmad(max_nmodes)   ! standard deviation
  REAL(wp) :: bi(max_nmodes)       ! Solubility parameters (see K&C)
  REAL(wp) :: betai(max_nmodes)    ! Solubility parameters (see K&C)
  LOGICAL :: use_mode(max_nmodes)  ! A flag that is set to true if there is significant aerosol number in a given mode
  LOGICAL :: l_aerosol_set     ! Has the aerosol been set?

  !-------------------------------------------------------------
  ! The following variables are the terms used in the K&C
  ! formulation of differential activation
  !-------------------------------------------------------------

  REAL(wp) :: sigmas(max_nmodes)   ! supersaturation dispersion
  REAL(wp) :: s0i(max_nmodes)      ! geometric mean supersaturation

  REAL(wp) :: ai(max_nmodes)        ! Coefficient for lookup method
  REAL(wp) :: logsigmas(max_nmodes) ! Log of sigmas

  !-------------------------------------------------------------
  ! Some additions which we are convenient to include here for
  ! Coupling to Dan's driver
  !-------------------------------------------------------------
  REAL(wp) :: vantHoff(max_nmodes)  ! van't Hoff factor
  REAL(wp) :: epsv(max_nmodes)      ! fraction of soluble mass
  REAL(wp) :: massMole(max_nmodes)  ! Molar mass
  REAL(wp) :: density(max_nmodes)  ! Molar mass

  ! Numerical issues
  REAL(wp), PARAMETER :: Nd_min=1.e0 ! minimum number concentration use to set use_mode

END MODULE shipway_parameters
