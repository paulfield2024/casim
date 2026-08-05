! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Physical constants used by the Shipway activation scheme.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Physical constants (gas constant, molecular weight of
!   water, surface tension, densities) specific to the Shipway
!   activation scheme.
!
! Paper reference:
!   Provides constants used by the Gordon et al. (2020)
!   activation calculation cited in Field et al. (2023)
!   Sec. A.6.2; the constants themselves are not tabulated in
!   the paper.
!
MODULE shipway_constants
  USE variable_precision, ONLY: wp
  USE mphys_constants, ONLY: pi

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SHIPWAY_CONSTANTS'

  REAL(wp), PARAMETER :: &
       Ru = 8.314472    & ! Universal gas constant
       ,Mw = 0.18015e-1 & ! Molecular weight of water 
       ,zetasa = 0.8e-1 & ! Surface tension at solution-air interface
       ,rhow = 1000.    & ! Density of water
       ,rho = 1           ! Density of  air 
  
  REAL(wp), PARAMETER :: &
         eps = 1.608   &  !(Rv/Rd)
         ,Rd = 287.05   &
         ,Rv = 461.5    &
         ,Dv = 0.226e-4 &
         ,Lv = 0.2501e7 &
         ,cp = 1005.    &
         ,ka = 0.243e-1 

  REAL(wp) :: alpha_c = 0.05 ! set do default value which we may be changed elsewhere
                         ! kinetic parameter, expressing probability
                         ! of water vapour molecules being incorporated
                         ! into droplet upon collision + dissolution 
                         ! kinetics (see fountoukis & nenes 2007 and 
                         ! Asa-Awuku and Nenes,2007)
  
  REAL(wp) :: &
       Dp_big  & ! droplet size upper bound (m)
       ,Dp_low & ! droplet size lower bound (m)
       ,delDp    ! difference between upper and lower bound
  

CONTAINS

  REAL(wp) FUNCTION Dv_mean(T, alpha_c)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    REAL(wp), INTENT(IN):: T ! temperature (Kelvin)
    REAL(wp), INTENT(IN):: alpha_c ! condensation coefficient
                               ! (Sassen & Dodd 1988 use 0.05)

    ! Local Variables

    REAL(wp) :: B ! a combination of variables

    CHARACTER(len=*), PARAMETER :: RoutineName='DV_MEAN'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    Dp_big=5.e-6
    Dp_low=min(0.207683*alpha_c**(-0.33048)*1e-6,Dp_big-1e-15)
    delDp=Dp_big-Dp_low

    B=(2.*Dv/alpha_c)*sqrt(2*pi*Mw/(Ru*T))
    Dv_mean=Dv*(delDp-B*log((Dp_big+B)/(Dp_low+B)))/delDp

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END FUNCTION Dv_mean
       
END MODULE shipway_constants
