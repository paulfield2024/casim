! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Sedimentation of hydrometeors (Eulerian flux-form, sedr/sedr_1M_2M) and CFL substep selection.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
MODULE sedimentation
  USE mphys_die, ONLY: throw_mphys_error, incorrect_opt, std_msg
  USE variable_precision, ONLY: wp, iwp
  USE mphys_parameters, ONLY: nz, hydro_params, cloud_params, rain_params   &
       , ice_params, snow_params, graupel_params, a_s, b_s
  USE passive_fields, ONLY: rho, rdz_on_rho, dz
  USE type_process, ONLY: process_name
  USE mphys_switches, ONLY: l_abelshipway, l_sed_3mdiff, &
       i_am4, l_ased, i_am5, i_am7, i_am8, i_am9, l_passivenumbers, l_passivenumbers_ice,   &
       i_an11, i_an12, &
       l_separate_rain, l_warm, i_aerosed_method, l_sed_icecloud_as_1m, l_sed_rain_1m, l_sed_snow_1m, &
       l_sed_graupel_1m, l_kfsm, & 
       cfl_vt_max, l_sed_eulexp
  USE mphys_constants, ONLY: rho0
  USE process_routines, ONLY: process_rate, i_psedr, i_asedr, i_asedl, i_psedl, &
       i_pseds, i_psedi, i_psedg, i_dsedi, i_dseds, i_dsedg
  USE thresholds, ONLY: thresh_small
  USE special, ONLY: Gammafunc

  USE lookup, ONLY: moment
  USE distributions, ONLY: dist_lambda, dist_mu, dist_n0, dist_lams
  USE aerosol_routines, ONLY: aerosol_active

  USE lsp_sedim_eulexp_mod, ONLY: lsp_sedim_eulexp
  ! lsp_sedim_eulexp_mod is a UM module, so can not be stored in CASIM or MONC repo. 
  ! if required please contact Adrian Hill. With MONC this uses a dummy routine in the CASIM
  ! component

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SEDIMENTATION'

  REAL(wp), ALLOCATABLE :: flux_n1(:)
  REAL(wp), ALLOCATABLE :: flux_n2(:)
  REAL(wp), ALLOCATABLE :: flux_n3(:)
  REAL(wp), ALLOCATABLE :: Grho(:)

!$OMP THREADPRIVATE(flux_n1, flux_n2, flux_n3, Grho)

  PUBLIC sedr, initialise_sedr, finalise_sedr, sedr_1M_2M, terminal_velocity_CFL
CONTAINS

  SUBROUTINE initialise_sedr()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='INITIALISE_SEDR'


    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ALLOCATE(flux_n1(nz), flux_n2(nz), flux_n3(nz), Grho(nz))

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE initialise_sedr

  SUBROUTINE finalise_sedr()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='FINALISE_SEDR'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DEALLOCATE(flux_n1, flux_n2, flux_n3, Grho)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE finalise_sedr

  SUBROUTINE sedr(ixy_inner, qfields, aeroact, dustact,   &
       params, procs, aerosol_procs, precip1d, l_doaerosol)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SEDR'

    INTEGER, INTENT(IN) :: ixy_inner

    REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
    TYPE(hydro_params), INTENT(IN) :: params
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
    TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)
    REAL(wp), INTENT(OUT) :: precip1d(nz)
    LOGICAL, OPTIONAL, INTENT(IN) :: l_doaerosol

    REAL(wp) :: dm1, dm2, dm3
    REAL(wp) :: dn1, dn2, dn3
    REAL(wp) :: m1, m2, m3
    REAL(wp) :: n1, n2, n3
    REAL(wp) :: hydro_mass
    REAL(wp) :: n0, lam, mu, u1r, u2r, u1r2, u2r2
!   real(wp) :: u3r, u3r2
    INTEGER :: k

    REAL(wp) :: p1, p2, p3
    REAL(wp) :: sp1, sp2, sp3
    REAL(wp) :: a_x, b_x, f_x, c_x
    REAL(wp) :: a2_x, b2_x, f2_x
    LOGICAL :: l_fluxin, l_fluxout

    TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
    ! which species we're depositing on.

    REAL(wp) :: dmac
    REAL(wp) :: dmad
    REAL(wp) :: dnumber_a, dnumber_d
    LOGICAL :: l_sedim_generic
    LOGICAL :: l_da_local  ! local tranfer of l_doaerosol
    ! If this is used, you can't trust the results.

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_da_local=.FALSE.
    IF (present(l_doaerosol)) l_da_local=l_doaerosol

    IF (l_kfsm) THEN
      l_sedim_generic = (params%id==ice_params%id)
    ELSE
      l_sedim_generic = .FALSE.
    END IF

    ! precip diag
    DO k = 1, nz
      precip1d(k) = 0.0
      flux_n1(k)  = 0.0
      Grho(k)=(rho0/rho(k,ixy_inner))**params%g_x
    END DO

    m1=0.0
    m2=0.0
    m3=0.0

    p1=params%p1
    p2=params%p2
    p3=params%p3

    IF (l_sed_3mdiff) THEN
      sp1=params%sp1
      sp2=params%sp2
      sp3=params%sp3
    ELSE
      sp1=params%p1
      sp2=params%p2
      sp3=params%p3
    END IF

    ! we don't want flexible approach in this version....
    IF (p1/=sp1 .OR. p2/=sp2 .OR. p3/=sp3) THEN
      WRITE(std_msg, '(A)') 'Cannot have flexible sedimentation options '//&
                            'with CASIM aerosol'
      CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                             std_msg)
    END IF

    IF (.NOT. l_kfsm) THEN
      a_x=params%a_x
      b_x=params%b_x
    END IF

    f_x=params%f_x
    c_x=params%c_x
    a2_x=params%a2_x
    b2_x=params%b2_x
    f2_x=params%f2_x

    IF (params%l_2m) flux_n2(:)=0.0     
!!$    if (params%l_3m) flux_n3=0.0      

    SELECT CASE (params%id)
    CASE (1_iwp) !cloud
      iproc=i_psedl
      iaproc=i_asedl
    CASE (2_iwp) !rain
      iproc=i_psedr
      iaproc=i_asedr
    CASE (3_iwp) !ice
      iproc=i_psedi
      iaproc=i_dsedi
    CASE (4_iwp) !snow
      iproc=i_pseds
      iaproc=i_dseds
    CASE (5_iwp) !graupel
      iproc=i_psedg
      iaproc=i_dsedg
    END SELECT

    DO k=nz-1, 1, -1

      ! initialize to zero
      dm1=0.0
      dm2=0.0
      dm3=0.0
      dn1=0.0
      dn2=0.0
      dn3=0.0
      n1=0.0
      n2=0.0
      n3=0.0
      m1=0.0
      m2=0.0
      m3=0.0

      hydro_mass=qfields(k, params%i_1m)
      IF (params%l_2m) m2=qfields(k, params%i_2m)
!!$      if (params%l_3m) m3=qfields(k, params%i_3m)

      l_fluxin=.FALSE.
      l_fluxout=.FALSE.
      IF (hydro_mass > thresh_small(params%i_1m)) l_fluxout=.TRUE.
      IF (qfields(k+1, params%i_1m) > thresh_small(params%i_1m)) l_fluxin=.TRUE.

      IF (l_sedim_generic) THEN
        a_x=a_s(k)
        b_x=b_s(k)
      ELSE
        a_x=params%a_x
        b_x=params%b_x
      END IF

      IF (l_fluxout) THEN
        m1=(hydro_mass/c_x)

        n0=dist_n0(k,params%id)
        mu=dist_mu(k,params%id)

        IF (l_sedim_generic) THEN
          lam=dist_lams(k,params%id,1)
        ELSE
          lam=dist_lambda(k,params%id)
        END IF

        IF (l_sed_3mdiff) THEN
          ! Moment transfer
          n1=moment(n0, lam, mu, sp1)
          n2=moment(n0, lam, mu, sp2)
          n3=moment(n0, lam, mu, sp3)
        ELSE
          n1=m1*rho(k,ixy_inner)
          n2=m2*rho(k,ixy_inner)
          n3=m3*rho(k,ixy_inner)
        END IF

        u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
             *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))

        IF (params%l_2m)     &
             u2r=a_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f_x)**(-(1.0+mu+sp2+b_x))) &
             *(Gammafunc(1.0+mu+sp2+b_x)/Gammafunc(1.0+mu+sp2))

!!$        if (params%l_3m)     &
!!$             u3r=a_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f_x)**(-(1.0+mu+sp3+b_x))) &
!!$             *(Gammafunc(1.0+mu+sp3+b_x)/Gammafunc(1.0+mu+sp3))

        IF (l_abelshipway .AND. params%id==rain_params%id) THEN ! rain can use abel and shipway formulation
          u1r2=a2_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f2_x)**(-(1.0+mu+sp1+b2_x)))   &
               *(Gammafunc(1.0+mu+sp1+b2_x)/Gammafunc(1.0+mu+sp1))
          u1r=u1r+u1r2

          IF (params%l_2m) THEN
            u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(1.0+mu+sp2+b2_x))) &
                 *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
            u2r=u2r+u2r2
          END IF

!!$          if (params%l_3m) then
!!$            u3r2=a2_x*Grho(k)*(lam**(1.0+mu+sp3)*(lam+f2_x)**(-(1.0+mu+sp3+b2_x))) &
!!$                 *(Gammafunc(1.0+mu+sp3+b2_x)/Gammafunc(1.0+mu+sp3))
!!$            u3r=u3r+u3r2
!!$          end if
        END IF

        ! fall speeds shouldn't get too big...
        u1r=min(u1r,params%maxv)

        ! For clouds and ice, we only use a 1M representation of sedimentation
        ! so that spurious size sorting doesn't lead to overactive autoconversion
        IF (params%id==cloud_params%id .OR. params%id==ice_params%id) THEN
          IF (l_sed_icecloud_as_1m)THEN
            IF (params%l_2m)u2r=u1r
          ELSE
            IF (params%l_2m)u2r=min(u2r,params%maxv)
          END IF
        ELSE
          IF (params%l_2m)u2r=min(u2r,params%maxv)
!!$          if (params%l_3m)u3r=min(u3r,params%maxv)
        END IF

        ! fall speeds shouldn't be negative (can happen with original AS formulation)
        u1r=max(u1r,0.0_wp)
        IF (params%l_2m)u2r=max(u2r,0.0_wp)
!!$        if (params%l_3m)u3r=max(u3r,0.0_wp)

         flux_n1(k)=n1*u1r
           
         IF (params%l_2m) flux_n2(k)=n2*u2r
        
        !if (params%l_3m) flux_n3(k)=n3*u3r

        precip1d(k) = flux_n1(k)*c_x
        
     END IF

      dmac=0.0
      dmad=0.0
      dnumber_a=0.0
      dnumber_d=0.0

      IF (l_fluxout) THEN !flux out (flux(k+1) will be zero if no flux in)
        dn1=(flux_n1(k+1)-flux_n1(k))*rdz_on_rho(k,ixy_inner)
        IF (params%l_2m) dn2=(flux_n2(k+1)-flux_n2(k))*rdz_on_rho(k,ixy_inner)
        !if (params%l_3m) dn3=(flux_n3(k+1)-flux_n3(k))*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        IF (l_ased .AND. l_da_local) THEN
          IF (params%id == cloud_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1 - flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k,ixy_inner)
            END IF
            IF (.NOT. l_warm) THEN
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-  &
                   flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
              IF (l_passivenumbers_ice .AND. dustact(k)%mact_mean > 0.0) THEN
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
              END IF
            END IF
          ELSE IF (params%id == rain_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers .AND. aeroact(k)%mact_mean > 0.0) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (.NOT. l_warm) THEN
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-  &
                   flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
              IF (l_passivenumbers_ice) THEN
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
              END IF
            END IF
          END IF

          IF (params%id == ice_params%id) THEN
            dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean-    &
                 flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-    &
                 flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1-flux_n2(k)*aeroact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
            END IF
          ELSE IF (params%id == snow_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-    &
                 flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-    &
                 flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
          ELSE IF (params%id == graupel_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean-    &
                 flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean-    &
                 flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)

            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3-flux_n2(k)*aeroact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3-flux_n2(k)*dustact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
            END IF
          END IF
        END IF
      ELSE IF (l_fluxin) THEN !flux in, but not out
        dn1=flux_n1(k+1)*rdz_on_rho(k,ixy_inner)
        IF (params%l_2m) dn2=flux_n2(k+1)*rdz_on_rho(k,ixy_inner)
        !if (params%l_3m) dn3=flux_n3(k+1)*rdz_on_rho(k)

        !============================
        ! aerosol processing
        !============================
        IF (l_ased .AND. l_da_local) THEN
          IF (params%id == cloud_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (.NOT. l_warm) THEN
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
              IF (l_passivenumbers_ice) THEN
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
              END IF
            END IF
          ELSE IF (params%id == rain_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (.NOT. l_warm) THEN
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
              IF (l_passivenumbers_ice) THEN
                dnumber_d=flux_n2(k+1)*dustact(k+1)%nratio2*rdz_on_rho(k,ixy_inner)
              END IF
            END IF
          END IF

          IF (params%id == ice_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
            END IF
          ELSE IF (params%id == snow_params%id) THEN
            dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
            IF (l_passivenumbers) THEN
              dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
            IF (l_passivenumbers_ice) THEN
              dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
            END IF
          ELSE IF (params%id == graupel_params%id) THEN
            IF (i_aerosed_method==1) THEN
              dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
              dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
              IF (l_passivenumbers) THEN
                dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
              END IF
              IF (l_passivenumbers_ice) THEN
                dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
              END IF
            ELSE
              WRITE(std_msg, '(A)') 'ERROR: GET RID OF i_aerosed_method variable!'
              CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                                     std_msg )
            END IF
          END IF
        END IF
      END IF

      ! Store the aerosol process terms...
      IF (l_da_local) THEN
        IF (params%id == cloud_params%id .OR. params%id == rain_params%id) THEN
          !liquid phase
          IF (l_separate_rain .AND. params%id == rain_params%id) THEN
            aerosol_procs(i_am5, iaproc%id)%column_data(k)=dmac
          ELSE
            aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
          END IF
          IF (.NOT. l_warm) aerosol_procs(i_am9, iaproc%id)%column_data(k) =dmad
          IF (l_passivenumbers) THEN
            aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
          END IF
          IF (l_passivenumbers_ice) THEN
            aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
          END IF

        ELSE
          !ice phase
          aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
          aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
          IF (l_passivenumbers) THEN
            aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
          END IF
          IF (l_passivenumbers_ice) THEN
            aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
          END IF
        END IF
      END IF

      dm1=dn1
      dm2=dn2
      !dm3=dn3

      procs(params%i_1m, iproc%id)%column_data(k)=c_x*dm1

      IF (params%l_2m) procs(params%i_2m, iproc%id)%column_data(k)=dm2

      !if (params%l_3m) procs(params%i_3m, iproc%id)%column_data(k)=dm3
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE sedr

SUBROUTINE sedr_1M_2M(ixy_inner, step_length, qfields, aeroact, dustact,   &
       params, procs, aerosol_procs, precip1d, l_doaerosol)

!!! This routine has the same functionality as sedr, except it will only work with 
!!! single moment or double moment settings, since the fallspeeds are calced using 
!!! pre-calculated gamma functions with fixed shape. This routine is default when 
!!! the third prognostic moment is switched off. 

! UM dummy routines
USE yomhook, ONLY: lhook, dr_hook
USE parkind1, ONLY: jprb, jpim
USE um_types, ONLY: real_lsprec

USE casim_parent_mod, ONLY: casim_parent, parent_um

IMPLICIT NONE

CHARACTER(len=*), PARAMETER :: RoutineName='SEDR_1M_2M'

INTEGER, INTENT(IN) :: ixy_inner

REAL(wp), INTENT(IN) :: step_length
REAL(wp), INTENT(IN), TARGET :: qfields(:,:)
TYPE(hydro_params), INTENT(IN) :: params
TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:)
TYPE(process_rate), INTENT(INOUT), TARGET :: procs(:,:)
TYPE(process_rate), INTENT(INOUT), TARGET :: aerosol_procs(:,:)
REAL(wp), INTENT(OUT) :: precip1d(nz)
LOGICAL, OPTIONAL, INTENT(IN) :: l_doaerosol

REAL(wp) :: dn1, dn2
REAL(wp) :: m1, m2
REAL(wp) :: n1, n2
REAL(wp) :: hydro_mass
REAL(wp) :: n0, lam, mu, u1r, u2r, u1r2, u2r2
INTEGER :: k

INTEGER :: i_1m, i_2m
REAL(wp) :: p1, p2 
REAL(wp) :: sp1, sp2 
REAL :: a_x, b_x, f_x, g_x, c_x, d_x
REAL(wp) :: a2_x, b2_x, f2_x
LOGICAL :: l_fluxin, l_fluxout

! AH: variables needed for lsp_sedim_eulexp. All declared using real_lsprec from the UM or 
!     the dummy real_lsprec from MONC casim component, which is wp
REAL(real_lsprec) :: lsp_sedim_c_x
REAL(real_lsprec) :: u1r_above, u2r_above
!   real(real_lsprec) :: u1w, u2w
!   real(real_lsprec) :: mixingratio_fromabove, mixingratio
!   real(real_lsprec) :: numberconc, numberconc_fromabove

REAL(real_lsprec) :: m0 
INTEGER :: points
REAL(real_lsprec) :: dhi(1), dhir(1),rhor(1),rhoin(1),mixratio_thislayer(1), &
     flux_fromabove(1),fallspeed_fromabove(1), fallspeed_thislayer(1),total_flux_out(1)
! End declaration of variables for lsp_sedim_eulexp

TYPE(process_name) :: iproc, iaproc  ! processes selected depending on
! which species we're depositing on.

REAL(wp) :: dmac
REAL(wp) :: dmad
REAL(wp) :: dnumber_a, dnumber_d
LOGICAL :: l_da_local  ! local tranfer of l_doaerosol
LOGICAL :: l_sedim_generic
! If this is used, you can't trust the results.

INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
REAL(KIND=jprb)               :: zhook_handle

!--------------------------------------------------------------------------
! End of header, no more declarations beyond here
!--------------------------------------------------------------------------
IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

IF (present(l_doaerosol)) THEN
  l_da_local=l_doaerosol
ELSE
  l_da_local=.FALSE.
END IF

IF (l_kfsm) THEN
  l_sedim_generic = (params%id==ice_params%id)
ELSE
  l_sedim_generic = .FALSE.
  a_x=params%a_x
  b_x=params%b_x
END IF
    
! precip diag
DO k = 1, nz
  precip1d(k) = 0.0
  flux_n1(k)  = 0.0
  Grho(k)=(rho0/rho(k,ixy_inner))**params%g_x
END DO

i_1m=params%i_1m
i_2m=params%i_2m
p1=params%p1
p2=params%p2

sp1=params%p1
sp2=params%p2

! we don't want flexible approach in this version....
IF (p1/=sp1 .OR. p2/=sp2) THEN
  WRITE(std_msg, '(A)') 'Cannot have flexible sedimentation options '//&
                        'with CASIM aerosol'
  CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
                         std_msg)
END IF

f_x=params%f_x
g_x=params%g_x
c_x=params%c_x
lsp_sedim_c_x = c_x
d_x=params%d_x
a2_x=params%a2_x
b2_x=params%b2_x
f2_x=params%f2_x

IF (params%l_2m) flux_n2(:)=0.0          

SELECT CASE (params%id)
CASE (1_iwp) !cloud
  iproc=i_psedl
  iaproc=i_asedl
CASE (2_iwp) !rain
  iproc=i_psedr
  iaproc=i_asedr
CASE (3_iwp) !ice
  iproc=i_psedi
  iaproc=i_dsedi
CASE (4_iwp) !snow
  iproc=i_pseds
  iaproc=i_dseds
CASE (5_iwp) !graupel
  iproc=i_psedg
  iaproc=i_dsedg
END SELECT

! AH : initialise terminal velocity for level+1 before loop starts
!     
u1r_above=0.0
u2r_above=0.0

DO k=nz-1, 1, -1

  ! initialize to zero
  dn1=0.0
  dn2=0.0
  n1=0.0
  n2=0.0
  m2=0.0

  hydro_mass=qfields(k, params%i_1m)
  IF (params%l_2m) m2=qfields(k, params%i_2m)

  IF (casim_parent == parent_um .AND. l_sed_eulexp ) THEN
    ! AH - the code below is required to run in the UM (GA) and 
    !      it over-rides the above conditions. This code should be wrapped in 
    !      in a condition that relates to l_sed_eulexp
    !      NOTE: at this stage it is not clear to me why this has to 
    !            be set to true when sed_eulexp is true!
    l_fluxin=.TRUE. 
    l_fluxout=.TRUE.
  ELSE
    l_fluxin=.FALSE.
    l_fluxout=.FALSE.
    IF (hydro_mass > thresh_small(params%i_1m)) l_fluxout=.TRUE.
    IF (qfields(k+1, params%i_1m) > thresh_small(params%i_1m)) l_fluxin=.TRUE.
  END IF

  IF (l_sedim_generic) THEN
    a_x=a_s(k)
    b_x=b_s(k)
  ELSE
    a_x=params%a_x
    b_x=params%b_x
  END IF

  IF (l_fluxout) THEN
    m1=(hydro_mass/c_x)

    n0=dist_n0(k,params%id)
    mu=dist_mu(k,params%id)
    lam=dist_lambda(k,params%id)

    n1=m1*rho(k,ixy_inner)
    n2=m2*rho(k,ixy_inner)

    IF (lam > 0.0_wp) THEN 
      IF (l_kfsm) THEN 
        u1r=a_x*Grho(k)*(lam**(1.0+mu+sp1)*(lam+f_x)**(-(1.0+mu+sp1+b_x)))     &
            *(Gammafunc(1.0+mu+sp1+b_x)/Gammafunc(1.0+mu+sp1))
      ELSE 
        u1r=a_x*Grho(k)*(lam**(params%exp_1_mu_sp1)*(lam+f_x)**(-(params%exp_1_mu_sp1_bx)))     &
            *(params%gam_1_mu_sp1_bx/params%gam_1_mu_sp1)
        IF (params%l_2m)     &
            u2r=a_x*Grho(k)*(lam**(params%exp_1_mu_sp2)*(lam+f_x)**(-(params%exp_1_mu_sp2_bx))) &
                *(params%gam_1_mu_sp2_bx/params%gam_1_mu_sp2)
      END IF
        
      IF (l_abelshipway .AND. params%id==rain_params%id) THEN ! rain can use abel and shipway formulation
        u1r2=a2_x*Grho(k)*(lam**(params%exp_1_mu_sp1)*(lam+f2_x)**(-(params%exp_1_mu_sp1_b2x)))   &
             *(params%gam_1_mu_sp1_b2x/params%gam_1_mu_sp1)
        u1r=u1r+u1r2
           
        IF (params%l_2m) THEN
          !u2r2=a2_x*Grho(k)*(lam**(1.0+mu+sp2)*(lam+f2_x)**(-(params%exp_1_mu_sp1_b2x))) &
          !     *(Gammafunc(1.0+mu+sp2+b2_x)/Gammafunc(1.0+mu+sp2))
          u2r2=a2_x*Grho(k)*(lam**(params%exp_1_mu_sp2)*(lam+f2_x)**(-(params%exp_1_mu_sp2_b2x))) &
               *(params%gam_1_mu_sp2_b2x/params%gam_1_mu_sp2)
          u2r=u2r+u2r2
        END IF
      END IF

      ! fall speeds shouldn't get too big...
      u1r=min(u1r,params%maxv)

      ! For clouds and ice, we only use a 1M representation of sedimentation
      ! so that spurious size sorting doesn't lead to overactive autoconversion
      IF (params%id==cloud_params%id .OR. params%id==ice_params%id) THEN
        IF (l_sed_icecloud_as_1m)THEN
          IF (params%l_2m)u2r=u1r
        ELSE
          IF (params%l_2m)u2r=min(u2r,params%maxv)
        END IF
      ELSE
        IF (params%l_2m)u2r=min(u2r,params%maxv)
      END IF

      IF ( (params%id==rain_params%id .AND. l_sed_rain_1m) .OR. &
           (params%id==snow_params%id .AND. l_sed_snow_1m) .OR. &
           (params%id==graupel_params%id .AND. l_sed_graupel_1m) )  THEN
        IF (params%l_2m)u2r=u1r
      END IF

      ! fall speeds shouldn't be negative (can happen with original AS formulation)
      u1r=max(u1r,0.0_wp)
      IF (params%l_2m)u2r=max(u2r,0.0_wp)
    ELSE 
      u1r=0.0_wp
      IF (params%l_2m)u2r=0.0_wp
    END IF
    flux_n1(k)=n1*u1r
        
    IF (params%l_2m) flux_n2(k)=n2*u2r
           

    ! AH : derive the sedimentation flux
    IF (casim_parent == parent_um .AND. l_sed_eulexp) THEN 
            ! use the method, which is based on UM and is most appropriate 
            ! (and stable) for  long timesteps
      !PRF use WB method -overwrite everything for now - put logicals in
      !print *, 'sed_eulexp called', k, hydro_mass
      IF (hydro_mass < thresh_small(params%i_1m)) THEN
        u1r=0.0
        n1=0.0
        IF (params%l_2m) THEN
          u2r=0.0
          n2=0.0
        END IF
      END IF
           
      points=1
      m0=1e-10/lsp_sedim_c_x  !for the mass - does not apply to number - 
                              !fix zero in one but not the other below)
      dhi(1)=step_length/dz(k,ixy_inner)
      dhir(1)=dz(k,ixy_inner)/step_length
      rhoin(1)=rho(k,ixy_inner)
      rhor(1)=1.0/rho(k,ixy_inner)
      mixratio_thislayer(1)=lsp_sedim_c_x*n1   /rho(k,ixy_inner)
      flux_fromabove(1)=lsp_sedim_c_x*flux_n1(k+1)
      fallspeed_fromabove(1)=u1r_above
      fallspeed_thislayer(1)=u1r
      total_flux_out(1)=0.0

      CALL lsp_sedim_eulexp(points,m0,dhi,dhir,rhoin,rhor,                     &
                            flux_fromabove, fallspeed_thislayer,               &
                            mixratio_thislayer, fallspeed_fromabove,           &
                            total_flux_out)
      flux_n1(k)=total_flux_out(1)/lsp_sedim_c_x
      u1r_above=fallspeed_thislayer(1)
           
      IF (params%l_2m) THEN
        points=1
        m0=1e-10/lsp_sedim_c_x  !for the mass - does not apply to number - 
                                !fix zero in one but not the other below)
        dhi(1)=step_length/dz(k,ixy_inner)
        dhir(1)=dz(k,ixy_inner)/step_length
        rhoin(1)=rho(k,ixy_inner)
        rhor(1)=1.0/rho(k,ixy_inner)
        mixratio_thislayer(1)=n2/rho(k,ixy_inner)
        flux_fromabove(1)=flux_n2(k+1)
        fallspeed_fromabove(1)=u2r_above
        fallspeed_thislayer(1)=u2r
        total_flux_out(1)=0.0
        
        CALL lsp_sedim_eulexp(points,m0,dhi,dhir,rhoin,rhor,                   &
                              flux_fromabove, fallspeed_thislayer,             &
                              mixratio_thislayer, fallspeed_fromabove,         &
                              total_flux_out)
        flux_n2(k)=total_flux_out(1)
        u2r_above=fallspeed_thislayer(1)
      END IF
           
      IF (flux_n1(k) < epsilon(1.0_wp) ) THEN
        flux_n1(k)=0.0
        u1r_above=0.0
      END IF

      ! if the fallspeed is within a typical number value of the precision
      ! used in lsp_sedim_eulexp, there is the chance a floating point
      ! error will occur on the next level down, hence stop sedimentation     
      IF (params%l_2m .AND. (flux_n2(k) < epsilon(1.0_wp) &
           .OR. u2r_above < 1.0e5*tiny(1.0_real_lsprec))) THEN
        flux_n2(k)=0.0
        u2r_above=0.0
      END IF    
        
!!PRF
    END IF
!         else if (parent_model
!                       ! ! Ported code that produces very different results to standard sed
!             flux_fromabove(1) = c_x * flux_n1(k+1)
!             mixingratio_fromabove = (flux_fromabove(1)*step_length/dz(k)*(1/rho(k)))
!             mixingratio = c_x * n1 * (1.0/rho(k))
           
!             if (mixingratio_fromabove + mixingratio > epsilon(1.0_wp)) then 
!                u1w=(u1r*mixingratio + &
!                     u1r_above*mixingratio_fromabove)/ & 
!                     (mixingratio_fromabove + mixingratio) 
!             else

!               u1w = 0.0_wp
              
!            endif

!            if (u1w > 0.0_wp) then 

!               flux_n1(k)= flux_fromabove(1) + (dz(k)/step_length) * &
!                    (rho(k)*mixingratio - flux_fromabove(1)/u1w) * &
!                    (1.0_wp - exp(-u1w*(step_length/dz(k))))
              
!               flux_n1(k) = flux_n1(k)/c_x
              
!            else
              
!               flux_n1(k) = 0.0
             
!            endif
              
!            u1r_above = u1w

!            if (flux_n1(k) < epsilon(1.0_wp)) then
!               flux_n1(k)=0.0
!               u1r_above=0.0
!            endif
           
!            if  (params%l_2m) Then
              
!               numberconc = n2 * (1.0/rho(k))
!               numberconc_fromabove = (flux_n2(k+1)*step_length/dz(k)*(1/rho(k)))

!               if (numberconc_fromabove + numberconc > epsilon(1.0_wp)) then 
!                  u2w=(u2r*numberconc + & 
!                       u2r_above*numberconc_fromabove)/ &
!                       (numberconc+numberconc_fromabove) ! PF:number also mass wtd 
!                                                         ! AH: no it is not, it is number weighted
!               else
                 
!                  u2w = 0.0_wp
                 
!               endif
              
!               if ( u2w > 0.0_wp ) then 
                 
!                  flux_n2(k)=flux_n2(k+1) + (dz(k)/step_length) * &
!                       (rho(k)*numberconc - flux_n2(k+1)/u2w) * &
!                       (1.0-exp(-u2w*(step_length/dz(k))))
!               else
                 
!                  flux_n2(k) = 0.0
                 
!               endif
                 
!               u2r_above = u2w

!               if (flux_n2(k) < epsilon(1.0_wp)) then
!                  flux_n2(k)=0.0
!                  u2r_above=0.0
!               endif
!            endif
! #endif
!        endif ! L_sed_eulexp

    precip1d(k) = flux_n1(k)*c_x
        ! diagnostic for precip

  END IF

  dmac=0.0
  dmad=0.0
  dnumber_a=0.0
  dnumber_d=0.0

  IF (l_fluxout) THEN !flux out (flux(k+1) will be zero if no flux in)
    dn1=(flux_n1(k+1)-flux_n1(k))*rdz_on_rho(k,ixy_inner)
    IF (params%l_2m) dn2=(flux_n2(k+1)-flux_n2(k))*rdz_on_rho(k,ixy_inner)

    !============================
    ! aerosol processing
    !============================
    IF (l_ased .AND. l_da_local) THEN
      IF (params%id == cloud_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean -      &
             flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)* rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1 -                       &
                    flux_n2(k)*aeroact(k)%nratio1)* rdz_on_rho(k,ixy_inner)
        END IF
        IF (.NOT. l_warm) THEN
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-     &
               flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
          IF (l_passivenumbers_ice .AND. dustact(k)%mact_mean > 0.0) THEN
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-                      &
                      flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
          END IF
        END IF
      ELSE IF (params%id == rain_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-       &
             flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers .AND. aeroact(k)%mact_mean > 0.0) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-                        &
                    flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (.NOT. l_warm) THEN
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-     &
               flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
          IF (l_passivenumbers_ice) THEN
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-                      &
                      flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
          END IF
        END IF
      END IF

      IF (params%id == ice_params%id) THEN
        dmac = (flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean-     &
             flux_n2(k)*aeroact(k)%nratio1*aeroact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        dmad = (flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean-     &
             flux_n2(k)*dustact(k)%nratio1*dustact(k)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1-                        &
                    flux_n2(k)*aeroact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1-                        &
                    flux_n2(k)*dustact(k)%nratio1)*rdz_on_rho(k,ixy_inner)
        END IF
      ELSE IF (params%id == snow_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean-       &
             flux_n2(k)*aeroact(k)%nratio2*aeroact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean-       &
             flux_n2(k)*dustact(k)%nratio2*dustact(k)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2-                        &
                    flux_n2(k)*aeroact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2-                        &
                    flux_n2(k)*dustact(k)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
      ELSE IF (params%id == graupel_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean-       &
             flux_n2(k)*aeroact(k)%nratio3*aeroact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean-       &
             flux_n2(k)*dustact(k)%nratio3*dustact(k)%mact3_mean)*rdz_on_rho(k,ixy_inner)

        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3-                        &
                    flux_n2(k)*aeroact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3-                        &
                    flux_n2(k)*dustact(k)%nratio3)*rdz_on_rho(k,ixy_inner)
        END IF
      END IF
    END IF
  ELSE IF (l_fluxin) THEN !flux in, but not out
    dn1=flux_n1(k+1)*rdz_on_rho(k,ixy_inner)
    IF (params%l_2m) dn2=flux_n2(k+1)*rdz_on_rho(k,ixy_inner)

    !============================
    ! aerosol processing
    !============================
    IF (l_ased .AND. l_da_local) THEN
      IF (params%id == cloud_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (.NOT. l_warm) THEN
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
          IF (l_passivenumbers_ice) THEN
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
          END IF
        END IF
      ELSE IF (params%id == rain_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (.NOT. l_warm) THEN
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
          IF (l_passivenumbers_ice) THEN
            dnumber_d=flux_n2(k+1)*dustact(k+1)%nratio2*rdz_on_rho(k,ixy_inner)
          END IF
        END IF
      END IF

      IF (params%id == ice_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio1*aeroact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio1*dustact(k+1)%mact1_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio1)*rdz_on_rho(k,ixy_inner)
        END IF
      ELSE IF (params%id == snow_params%id) THEN
        dmac=(flux_n2(k+1)*aeroact(k+1)%nratio2*aeroact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        dmad=(flux_n2(k+1)*dustact(k+1)%nratio2*dustact(k+1)%mact2_mean)*rdz_on_rho(k,ixy_inner)
        IF (l_passivenumbers) THEN
          dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
        IF (l_passivenumbers_ice) THEN
          dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio2)*rdz_on_rho(k,ixy_inner)
        END IF
      ELSE IF (params%id == graupel_params%id) THEN
        IF (i_aerosed_method==1) THEN
          dmac=(flux_n2(k+1)*aeroact(k+1)%nratio3*aeroact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
          dmad=(flux_n2(k+1)*dustact(k+1)%nratio3*dustact(k+1)%mact3_mean)*rdz_on_rho(k,ixy_inner)
          IF (l_passivenumbers) THEN
            dnumber_a=(flux_n2(k+1)*aeroact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
          END IF
          IF (l_passivenumbers_ice) THEN
            dnumber_d=(flux_n2(k+1)*dustact(k+1)%nratio3)*rdz_on_rho(k,ixy_inner)
          END IF
        ELSE
          WRITE(std_msg, '(A)') 'ERROR: GET RID OF i_aerosed_method variable!'
          CALL throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, std_msg )
        END IF
      END IF
    END IF
  END IF

  ! Store the aerosol process terms...
  IF (l_da_local) THEN
    IF (params%id == cloud_params%id .OR. params%id == rain_params%id) THEN
      !liquid phase
      IF (l_separate_rain .AND. params%id == rain_params%id) THEN
        aerosol_procs(i_am5, iaproc%id)%column_data(k)=dmac
      ELSE
        aerosol_procs(i_am4, iaproc%id)%column_data(k)=dmac
      END IF
      IF (.NOT. l_warm) aerosol_procs(i_am9, iaproc%id)%column_data(k) =dmad
      IF (l_passivenumbers) THEN
        aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
      END IF
      IF (l_passivenumbers_ice) THEN
        aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
      END IF

    ELSE
      !ice phase
      aerosol_procs(i_am7, iaproc%id)%column_data(k)=dmad
      aerosol_procs(i_am8, iaproc%id)%column_data(k)=dmac
      IF (l_passivenumbers) THEN
        aerosol_procs(i_an11, iaproc%id)%column_data(k)=dnumber_a
      END IF
      IF (l_passivenumbers_ice) THEN
        aerosol_procs(i_an12, iaproc%id)%column_data(k)=dnumber_d
      END IF
    END IF
  END IF

  procs(params%i_1m, iproc%id)%column_data(k)=c_x*dn1

  IF (params%l_2m) procs(params%i_2m, iproc%id)%column_data(k)=dn2

END DO

IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

END SUBROUTINE sedr_1M_2M


  SUBROUTINE terminal_velocity_CFL(dt, vt, nsubseds_hydro, sed_length_hydro, nsubseds_max, & 
    sed_length_max, mindz)

    ! work out the number of substeps for each hydrometeor
    ! using the minimum dz and maxv to determine a CFL limit
    ! NOTE: if the parent_model has a fixed timestep this calc
    !       should be done in initialise_sedr. If there is a dynamic
    !       timestep, this calc has to be done on every timestep
    
    IMPLICIT NONE

    REAL(wp), INTENT(IN) :: dt  ! timestep from parent model
    REAL(wp), INTENT(IN) :: vt  ! max hydrometeor terminal velocity for a hydrometeor
                                ! (comes from mphys_params)
    ! Both of the following varaibles are derived from mphys_switches
    INTEGER, INTENT(IN)  ::  nsubseds_max  ! number of sedimentation substeps for sedimentation
    REAL(wp), INTENT(IN) :: sed_length_max ! sedimentation substep length sedimentation
    ! cfl_vt_max is declared and set in mphys_switches
    REAL(wp), INTENT(IN) :: mindz  ! pass it in 
    

    INTEGER, INTENT(OUT)  ::  nsubseds_hydro  ! number of sedimentation substeps for a hydrometeor
    REAL(wp), INTENT(OUT) :: sed_length_hydro ! sedimentation substep length for a hydormeteor

    ! local variables
    REAL(wp) :: cfl_vt ! CFL value for the terminal velocity argument
    

    
    cfl_vt =( dt * vt ) / mindz

    !!write(6,*) 'PRF1s ',cfl_vt, cfl_vt_max

    IF ( cfl_vt > cfl_vt_max ) THEN
       nsubseds_hydro = max(1, ceiling(cfl_vt / cfl_vt_max))
       sed_length_hydro = dt / (REAL(nsubseds_hydro, KIND=wp))
       
       IF (sed_length_hydro< sed_length_max) THEN
          nsubseds_hydro = nsubseds_max
          sed_length_hydro = sed_length_max
       END IF

       
    ELSE    
       nsubseds_hydro = 1
       sed_length_hydro = dt
    END IF

   !! write(6,*) 'PRF2s ',dt,nsubseds_hydro,sed_length_hydro,nsubseds_max,sed_length_max
  
  END SUBROUTINE terminal_velocity_CFL

END MODULE sedimentation
