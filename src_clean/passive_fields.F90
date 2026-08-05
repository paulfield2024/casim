! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Column fields that are read but not updated during a microphysics timestep (rho, pressure, z, exner, w, tke, ...).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
!Fields which aren't updated during the microphysics calculation
MODULE passive_fields
  USE mphys_die, ONLY: throw_mphys_error, bad_values, std_msg
  USE variable_precision, ONLY: wp
  USE qsat_funs, ONLY: qsaturation
  USE mphys_switches, ONLY: i_th
  USE mphys_parameters, ONLY: nxy_inner

  IMPLICIT NONE
  PRIVATE

  REAL(wp), ALLOCATABLE :: rho(:,:), pressure(:,:), z(:), exner(:,:), rexner(:,:)
  REAL(wp), ALLOCATABLE :: dz(:,:)
  REAL(wp), ALLOCATABLE :: qws(:,:), qws0(:,:), TdegC(:,:), TdegK(:,:), w(:,:), tke(:)

!$OMP THREADPRIVATE(rho, pressure, z, exner, rexner, dz, &
!$OMP               qws, qws0, TdegC, TdegK, w, tke)

  REAL(wp), ALLOCATABLE :: rhcrit_1d(:)
!$OMP THREADPRIVATE(rhcrit_1d)

  REAL(wp), ALLOCATABLE :: rdz_on_rho(:,:)
!$OMP THREADPRIVATE(rdz_on_rho)

  REAL(wp) :: dt
  REAL(wp), ALLOCATABLE :: min_dz(:) ! minimum vertical resolution
  INTEGER :: kl, ku, nz

!$OMP THREADPRIVATE(dt,min_dz,nz,kl,ku)

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='PASSIVE_FIELDS'

  PUBLIC set_passive_fields, rho, pressure, initialise_passive_fields, z, qws, &
       exner, rexner, dz, qws0, &
       TdegC, TdegK, w, tke, rhcrit_1d, rdz_on_rho, min_dz
CONTAINS

  SUBROUTINE initialise_passive_fields(kl_arg, ku_arg)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='INITIALISE_PASSIVE_FIELDS'

    INTEGER, INTENT(IN) :: kl_arg, ku_arg

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


    kl=kl_arg
    ku=ku_arg
    nz=ku-kl+1
    ALLOCATE(rho(nz, nxy_inner))
    ALLOCATE(pressure(nz, nxy_inner))
    ALLOCATE(exner(nz, nxy_inner))
    ALLOCATE(rexner(nz, nxy_inner))
    ALLOCATE(dz(nz, nxy_inner))
    ALLOCATE(w(nz, nxy_inner))
    ALLOCATE(tke(nz)) !! tke(:) is not used in other subroutines. So no need to have inner dimension
    ALLOCATE(rdz_on_rho(nz, nxy_inner))
    ALLOCATE(qws(nz, nxy_inner))
    ALLOCATE(qws0(nz, nxy_inner))
    ALLOCATE(TdegK(nz, nxy_inner))
    ALLOCATE(TdegC(nz, nxy_inner))
    ALLOCATE(min_dz(nxy_inner))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE initialise_passive_fields

  SUBROUTINE set_passive_fields(nxy_inner_loop, ixy_outer, is_in, js_in, je_in, &
                                dt_in, rho_in, p_in, exner_in,   &
                                dz_in, w_in, tke_in, qfields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SET_PASSIVE_FIELDS'

    INTEGER, INTENT(IN) :: nxy_inner_loop
    INTEGER, INTENT(IN) :: ixy_outer
    INTEGER, INTENT(IN) :: is_in, js_in, je_in

    REAL(wp), INTENT(IN) :: dt_in
    REAL(wp), INTENT(IN) :: rho_in(:,:,:), p_in(:,:,:), exner_in(:,:,:)
    REAL(wp), INTENT(IN) :: dz_in(:,:,:)
    REAL(wp), INTENT(IN) :: w_in(:,:,:), tke_in(:,:,:)
    REAL(wp), INTENT(IN), TARGET :: qfields(:,:,:)
    INTEGER :: k, ixy_inner, ixy, jy, ix

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    dt=dt_in

    DO ixy_inner=1, nxy_inner_loop

      ixy  = (ixy_outer-1)*nxy_inner + ixy_inner
      jy = modulo(ixy-1,(je_in-js_in+1))+js_in
      ix = (ixy-1)/(je_in-js_in+1)+is_in


      DO k=1,nz
        rho(k,ixy_inner)=rho_in(k,ix,jy)
        pressure(k, ixy_inner)=p_in(k,ix,jy)
        exner(k, ixy_inner)=exner_in(k,ix,jy)
        rexner(k, ixy_inner)=1.0/exner(k, ixy_inner)
        w(k, ixy_inner)=w_in(k,ix,jy)
        !tke(k)=tke_in(kl:ku)
        dz(k, ixy_inner)=dz_in(k,ix,jy)
        rdz_on_rho(k, ixy_inner)=1.0/(dz_in(k,ix,jy)*rho_in(k,ix,jy))
        !  do k=1,nz
        TdegK(k, ixy_inner)=qfields(k,i_th,ixy_inner)*exner(k, ixy_inner)
        TdegC(k, ixy_inner)=TdegK(k, ixy_inner)-273.15
        qws(k, ixy_inner)=qsaturation(TdegK(k, ixy_inner), pressure(k, ixy_inner)/100.0)
        qws0(k, ixy_inner)=qsaturation(273.15_wp, pressure(k, ixy_inner)/100.0)
      END DO ! k

      qws(nz, ixy_inner)=1.0e-8

      min_dz(ixy_inner) = minval(dz(:,ixy_inner))

      IF (any(qws(:,ixy_inner)==0.0)) THEN
        WRITE(std_msg, '(A)') 'Error in saturation calculation - qws is zero'
        CALL throw_mphys_error( bad_values, ModuleName//':'//RoutineName,     &
                                std_msg )

      END IF
    END DO ! ixy_inner
    !nullify(theta)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE set_passive_fields
END MODULE passive_fields
