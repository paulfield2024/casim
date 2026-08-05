! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Enforces positivity/consistency of hydrometeor and aerosol moments after each process update.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Post-process clean-up: removes negligibly small residual
!   mass/number left in a hydrometeor after the process rates
!   have been applied, to avoid numerical noise.
!
! Paper reference:
!   Related to the threshold/limiting philosophy discussed
!   qualitatively in Field et al. (2023) Appendix A.13, but the
!   specific clean-up thresholds here are not given a numbered
!   equation in the paper.
!
MODULE mphys_tidy
  USE variable_precision, ONLY: wp
  USE process_routines, ONLY: process_rate,  process_name
  USE aerosol_routines, ONLY: aerosol_active
  USE thresholds, ONLY: thresh_tidy, thresh_atidy
  USE passive_fields, ONLY: exner, pressure
  USE mphys_switches, ONLY:                       &
       i_qv, i_ql, i_nl, i_qr, i_nr, i_m3r, i_th, &
       i_qi, i_ni, i_qs, i_ns, i_m3s,             &
       i_qg, i_ng, i_m3g,                         &
!       l_3mr, l_3mg, l_3ms,                      &
       l_2mc, l_2mr,                              &
       l_2mi, l_2ms, l_2mg,                       &
       i_an2, i_am2, i_am4, i_am5, l_warm,        &
       i_an6, i_am6, i_am7, i_am8, i_am9,         &
       i_an11, i_an12,                            &
       l_process, ntotalq, ntotala,               &
       i_qstart, i_nstart, i_m3start,             &
       l_separate_rain, l_tidy_conserve_E, l_tidy_conserve_q, &
       l_passivenumbers, l_passivenumbers_ice
  USE mphys_constants, ONLY: Lv, Ls, cp
  USE qsat_funs, ONLY: qisaturation
  USE mphys_parameters, ONLY: hydro_params
  USE mphys_die, ONLY: throw_mphys_error, bad_values, warn, std_msg

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_TIDY'

  LOGICAL, PARAMETER :: l_rescale_on_number = .FALSE.
! logical :: l_tidym3 = .false.  ! Don't tidy based on m3 values

  REAL(wp), ALLOCATABLE :: thresh(:), athresh(:), qin_thresh(:)
!$OMP THREADPRIVATE(thresh, athresh, qin_thresh)
  
  LOGICAL :: current_l_negonly, current_qin_l_negonly
!$OMP THREADPRIVATE(current_l_negonly, current_qin_l_negonly)

  PUBLIC initialise_mphystidy, finalise_mphystidy, qtidy, ensure_positive, ensure_saturated, tidy_qin, &
       tidy_ain, ensure_positive_aerosol
CONTAINS

  SUBROUTINE initialise_mphystidy()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='INITIALISE_MPHYSTIDY'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    ALLOCATE(thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)), qin_thresh(lbound(thresh_tidy,1):ubound(thresh_tidy,1)))

    IF (l_process) THEN
      ALLOCATE(athresh(lbound(thresh_atidy,1):ubound(thresh_atidy,1)))
    END IF

    current_l_negonly=.TRUE.
    CALL recompute_constants(.TRUE.)
    current_qin_l_negonly=.TRUE.
    CALL recompute_qin_constants(.TRUE.)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE initialise_mphystidy

  SUBROUTINE recompute_qin_constants(l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='RECOMPUTE_QIN_CONSTANTS'

    LOGICAL, INTENT(IN) :: l_negonly

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    qin_thresh=thresh_tidy
    IF (l_negonly) THEN
      qin_thresh=0.0*qin_thresh
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE recompute_qin_constants

  SUBROUTINE recompute_constants(l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='RECOMPUTE_CONSTANTS'

    LOGICAL, INTENT(IN) :: l_negonly

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (l_negonly) THEN
      thresh=0.0
    ELSE
      thresh=thresh_tidy
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE recompute_constants

  SUBROUTINE finalise_mphystidy()

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='FINALISE_MPHYSTIDY'

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (l_process) THEN
      DEALLOCATE(athresh)
    END IF
    DEALLOCATE(thresh, qin_thresh)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE finalise_mphystidy

  SUBROUTINE qtidy(ixy_inner, dt, nz, qfields, procs, aerofields, aeroact, dustact, aeroice, dustliq, &
       aeroprocs, i_proc, i_aproc, l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='QTIDY'

    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: qfields(:,:), aerofields(:,:)
    TYPE(aerosol_active), INTENT(IN) :: aeroact(:), dustact(:), aeroice(:), dustliq(:)
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
    TYPE(process_rate), INTENT(INOUT) :: aeroprocs(:,:)
    TYPE(process_name), INTENT(IN) :: i_proc, i_aproc
    LOGICAL, INTENT(IN), OPTIONAL :: l_negonly

    LOGICAL :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
    LOGICAL :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
    LOGICAL :: qg_reset, ng_reset, m3g_reset
    LOGICAL :: am4_reset, am5_reset, am7_reset, am8_reset, am9_reset
    LOGICAL :: an11_reset, an12_reset

    REAL(wp) :: dmass, dnumber

    LOGICAL :: l_qsig(0:ntotalq), l_qpos, l_qsmall, l_qsneg(0:ntotalq)
    !    l_qpos: q variable is positive
    !    l_qsmall:   q variable is positive, but below tidy threshold
    !    l_qsneg:    q variable is small or negative

    LOGICAL :: l_qice, l_qliquid

    LOGICAL :: l_apos, l_asmall, l_asneg(ntotala), l_asig(ntotala)

    INTEGER :: iq, k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    l_qsig(0)=.FALSE.
    l_qsneg(0)=.FALSE.
    
    DO k = 1, nz
    nr_reset=.FALSE.
    m3r_reset=.FALSE.
    qi_reset=.FALSE.
    qs_reset=.FALSE.
    ns_reset=.FALSE.
    m3s_reset=.FALSE.
    qg_reset=.FALSE.
    ng_reset=.FALSE.
    m3g_reset=.FALSE.

    am4_reset=.FALSE.
    am5_reset=.FALSE.
    am7_reset=.FALSE.
    am8_reset=.FALSE.
    am9_reset=.FALSE.
    an11_reset=.FALSE.
    an12_reset=.FALSE.

    IF (present(l_negonly)) THEN
      IF (l_negonly .NEQV. current_l_negonly) THEN
        current_l_negonly=l_negonly
        CALL recompute_constants(l_negonly)
      END IF
    ELSE IF (current_l_negonly) THEN
      current_l_negonly=.FALSE.
      CALL recompute_constants(.FALSE.)
    END IF

    DO iq=1, ntotalq
      l_qsig(iq)=qfields(k, iq) > thresh(iq)

      l_qpos=qfields(k, iq) > 0.0

      l_qsmall = (.NOT. l_qsig(iq)) .AND. l_qpos

      l_qsneg(iq)=qfields(k, iq) < 0.0 .OR. l_qsmall
    END DO

    IF (l_process) THEN
      athresh=thresh_atidy
      IF (present(l_negonly)) THEN
        IF (l_negonly) athresh=0.0
      END IF

      DO iq=1, ntotala
        l_asig(iq)=aerofields(k, iq) > athresh(iq)

        l_apos=aerofields(k, iq) > 0.0

        l_asmall=(.NOT. l_asig(iq)) .AND. l_apos

        l_asneg(iq)=aerofields(k, iq) < 0.0 .OR. l_asmall
      END DO
    END IF

    l_qliquid=l_qsig(i_ql) .OR. l_qsig(i_qr)
    l_qice=l_qsig(i_qi) .OR. l_qsig(i_qs) .OR. l_qsig(i_qg)

    ! Tidying of small and negative numbers and/or incompatible numbers (e.g.nl>0 and ql=0)
    ! - Mass and energy conserving...
    !==============================
    ! What should be reset?
    !==============================
    ql_reset=l_qsneg(i_ql)
    IF (l_2mc) THEN
      nl_reset=l_qsneg(i_nl) .OR. (l_qsig(i_nl) .AND. ql_reset)
      ql_reset=ql_reset .OR. (l_qsig(i_ql) .AND. nl_reset)
    END IF
    qr_reset=l_qsneg(i_qr)
    IF (l_2mr) THEN
      nr_reset=l_qsneg(i_nr) .OR. (l_qsig(i_nr) .AND. qr_reset)
      qr_reset=qr_reset .OR. (l_qsig(i_qr) .AND. nr_reset)
    END IF
    ! if (l_3mr .and. l_tidym3) then
    !   m3r_reset=l_qsneg(i_m3r) .or. (l_qsig(i_m3r) .and. (qr_reset .or. nr_reset))
    !   nr_reset=nr_reset .or. (l_qsig(i_nr) .and. m3r_reset)
    !   qr_reset=qr_reset .or. (l_qsig(i_qr) .and. m3r_reset)
    ! end if

    IF (.NOT. l_warm) THEN
      qi_reset=l_qsneg(i_qi)
      IF (l_2mi) THEN
        ni_reset=l_qsneg(i_ni) .OR. (l_qsig(i_ni) .AND. qi_reset)
        qi_reset=qi_reset .OR. (l_qsig(i_qi) .AND. ni_reset)
      END IF

      qs_reset=l_qsneg(i_qs)
      IF (l_2ms) THEN
        ns_reset=l_qsneg(i_ns) .OR. (l_qsig(i_ns) .AND. qs_reset)
        qs_reset=qs_reset .OR. (l_qsig(i_qs) .AND. ns_reset)
      END IF
      ! if (l_3ms .and. l_tidym3) then
      !   m3s_reset=l_qsneg(i_m3s) .or. (l_qsig(i_m3s) .and. (qs_reset .or. ns_reset))
      !   ns_reset=ns_reset .or. (l_qsig(i_ns) .and. m3s_reset)
      !   qs_reset=qs_reset .or. (l_qsig(i_qs) .and. m3s_reset)
      ! end if

      qg_reset=l_qsneg(i_qg)
      IF (l_2mg) THEN
        ng_reset=l_qsneg(i_ng) .OR. (l_qsig(i_ng) .AND. qg_reset)
        qg_reset=qg_reset .OR. (l_qsig(i_qg) .AND. ng_reset)
      END IF
      ! if (l_3mg .and. l_tidym3) then
      !   m3g_reset=l_qsneg(i_m3g) .or. (l_qsig(i_m3g) .and. (qg_reset .or. ng_reset))
      !   ng_reset=ng_reset .or. (l_qsig(i_ng) .and. m3g_reset)
      !   qg_reset=qg_reset .or. (l_qsig(i_qg) .and. m3g_reset)
      ! end if
    END IF

    !===========================================================
    ! Aerosol tests...
    !===========================================================
    IF (l_process) THEN
      ! Aerosols in liquid water

      ! If small/neg values...
      IF (l_asneg(i_am4))am4_reset=.TRUE.
      IF (l_passivenumbers) THEN
        IF (l_asneg(i_an11))an11_reset=.TRUE.  ! can be in liq or ice
      END IF
      IF (l_passivenumbers_ice) THEN
        IF (l_asneg(i_an12))an12_reset=.TRUE.  ! can be in liq or ice
      END IF

      IF (l_separate_rain) THEN
        IF (l_asneg(i_am5))am5_reset=.TRUE.
      END IF
      IF (.NOT. l_Warm) THEN
        IF (l_asneg(i_am9))am9_reset=.TRUE.
      END IF
      ! If no hydrometeors...
      IF ((ql_reset .AND. qr_reset) .OR. .NOT. l_qliquid) THEN
        IF (l_asig(i_am4))am4_reset=.TRUE.
        IF (.NOT.l_warm) THEN
          IF (l_asig(i_am9))am9_reset=.TRUE.
        END IF
        IF (l_separate_rain) THEN
          IF (l_asig(i_am5)) am5_reset=.TRUE.
        END IF
      END IF

      ! If no active aerosol, then we shouldn't have any hydrometeor...(what about SIP?)
      ql_reset=ql_reset .OR. (am4_reset .AND. am9_reset .AND. l_qsig(i_ql))
      qr_reset=qr_reset .OR. (am4_reset .AND. am9_reset .AND. l_qsig(i_qr))
      qr_reset=qr_reset .OR. (am5_reset .AND. l_qsig(i_qr))

      ! Aerosols in ice
      ! If small/neg values...
      IF (.NOT. l_Warm)THEN
        IF (l_asneg(i_am7)) am7_reset=.TRUE.
        IF (l_asneg(i_am8)) THEN
          am8_reset=.TRUE.
        END IF
        ! If no hydrometeors...
        IF ((qi_reset .AND. qs_reset .AND. qg_reset) .OR. .NOT. l_qice) THEN
          IF (l_asig(i_am7)) am7_reset=.TRUE.
          IF (l_asig(i_am8)) THEN
            am8_reset=.TRUE.
          END IF
        END IF
      END IF

      ! If no active aerosol, then we shouldn't have any hydrometeor...
      qi_reset=qi_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qi))
      qs_reset=qs_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qs))
      qg_reset=qg_reset .OR. (am7_reset .AND. am8_reset .AND. l_qsig(i_qg))
    END IF

    !===========================================================
    ! Consistency following aerosol
    !===========================================================
    nl_reset=ql_reset .AND. l_qsig(i_nl)
    nr_reset=nr_reset .AND. l_qsig(i_nr)
    m3r_reset=m3r_reset .AND. l_qsig(i_m3r)
    ni_reset=qi_reset .AND. l_qsig(i_ni)
    ns_reset=ns_reset .AND. l_qsig(i_ns)
    m3s_reset=m3s_reset .AND. l_qsig(i_m3s)
    ng_reset=ng_reset .AND. l_qsig(i_ng)
    m3g_reset=m3g_reset .AND. l_qsig(i_m3g)

    !==============================
    ! Now reset things...
    !==============================
    IF (ql_reset .OR. nl_reset) THEN
      dmass=qfields(k, i_ql)/dt
      procs(i_ql,i_proc%id)%column_data(k)=-dmass
      IF (l_tidy_conserve_q)procs(i_qv,i_proc%id)%column_data(k)=dmass
      IF (l_tidy_conserve_E)procs(i_th,i_proc%id)%column_data(k)=-Lv*dmass/cp/exner(k,ixy_inner)
      IF (l_2mc) THEN
        dnumber=qfields(k, i_nl)/dt
        procs(i_nl,i_proc%id)%column_data(k)=-dnumber
      END IF
      !--------------------------------------------------
      ! aerosol - not reset, but adjusted acordingly
      !--------------------------------------------------
      IF (l_process) THEN
        IF (.NOT. am4_reset) THEN
          dmass=aeroact(k)%mact1/dt
          dnumber=aeroact(k)%nact1/dt
          IF (dmass > 0.0) THEN
            aeroprocs(i_am4,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=dnumber
          END IF
        END IF
        IF (.NOT. am9_reset) THEN
          dmass=dustliq(k)%mact1/dt
          dnumber=dustliq(k)%nact1/dt
          IF (dmass > 0.0) THEN
            aeroprocs(i_am9,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=dnumber
          END IF
        END IF
      END IF
    END IF

    IF (qr_reset .OR. nr_reset .OR. m3r_reset) THEN
      dmass=qfields(k, i_qr)/dt
      procs(i_qr,i_proc%id)%column_data(k)=-dmass
      IF (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=             &
                                   procs(i_qv,i_proc%id)%column_data(k) + dmass
      IF (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=             &
          procs(i_th,i_proc%id)%column_data(k) - Lv*dmass/cp/exner(k,ixy_inner)
      IF (l_2mr) THEN
        dnumber=qfields(k, i_nr)/dt
        procs(i_nr,i_proc%id)%column_data(k)=-dnumber
      END IF
      ! if (l_3mr) then
      !   procs(i_m3r,i_proc%id)%column_data(k)=-qfields(k, i_m3r)/dt
      ! end if
      !--------------------------------------------------
      ! aerosol
      !--------------------------------------------------
      IF (l_process) THEN
        IF (.NOT. am4_reset) THEN
          dmass=aeroact(k)%mact2/dt
          dnumber=aeroact(k)%nact2/dt
          aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
          aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
          aeroprocs(i_am4,i_aproc%id)%column_data(k)=aeroprocs(i_am4,i_aproc%id)%column_data(k)-dmass
        END IF
        IF (.NOT. am5_reset) THEN
          IF (l_separate_rain) THEN
            dmass=aeroact(k)%mact2/dt
            dnumber=aeroact(k)%nact2/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am5,i_aproc%id)%column_data(k)=aeroprocs(i_am5,i_aproc%id)%column_data(k)-dmass
          END IF
        END IF
        IF (.NOT. am9_reset .AND. .NOT. l_warm) THEN
          dmass=dustliq(k)%mact2/dt
          dnumber=dustliq(k)%nact2/dt
          IF (dmass>0.0) THEN
            aeroprocs(i_am9,i_aproc%id)%column_data(k)=-dmass
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=dnumber
          END IF
        END IF
      END IF
    END IF

    IF (.NOT. l_warm) THEN
      IF (qi_reset .OR. ni_reset) THEN
        dmass=qfields(k, i_qi)/dt
        procs(i_qi,i_proc%id)%column_data(k)=-dmass
        IF (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=           &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        IF (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=           &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        IF (l_2mi) THEN
          dnumber=qfields(k, i_ni)/dt
          procs(i_ni,i_proc%id)%column_data(k)=-dnumber
        END IF
        !--------------------------------------------------
        !aerosol
        !--------------------------------------------------
        IF (l_process) THEN
          IF (.NOT. am7_reset) THEN
            dmass=dustact(k)%mact1/dt
            dnumber=dustact(k)%nact1/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          END IF
          IF (.NOT. am8_reset) THEN
            dmass=aeroice(k)%mact1/dt
            dnumber=aeroice(k)%nact1/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          END IF
        END IF
      END IF

      IF (qs_reset .OR. ns_reset .OR. m3s_reset) THEN
        dmass=qfields(k, i_qs)/dt
        procs(i_qs,i_proc%id)%column_data(k)=-dmass
        IF (l_tidy_conserve_q) procs(i_qv,i_proc%id)%column_data(k)=           &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        IF (l_tidy_conserve_E) procs(i_th,i_proc%id)%column_data(k)=           &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        IF (l_2ms) THEN
          dnumber=qfields(k, i_ns)/dt
          procs(i_ns,i_proc%id)%column_data(k)=-dnumber
        END IF
        ! if (l_3ms) then
        !   procs(i_m3s,i_proc%id)%column_data(k)=-qfields(k, i_m3s)/dt
        ! end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        IF (l_process) THEN
          IF (.NOT. am7_reset) THEN
            dmass=dustact(k)%mact2/dt
            dnumber=dustact(k)%nact2/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          END IF
          IF (.NOT. am8_reset) THEN
            dmass=aeroice(k)%mact2/dt
            dnumber=aeroice(k)%nact2/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          END IF
        END IF
      END IF

      IF (qg_reset .OR. ng_reset .OR. m3g_reset) THEN
        dmass=qfields(k, i_qg)/dt
        procs(i_qg,i_proc%id)%column_data(k)=-dmass
        IF (l_tidy_conserve_q)procs(i_qv,i_proc%id)%column_data(k)=            &
                                     procs(i_qv,i_proc%id)%column_data(k)+dmass
        IF (l_tidy_conserve_E)procs(i_th,i_proc%id)%column_data(k)=            &
            procs(i_th,i_proc%id)%column_data(k)-Ls*dmass/cp/exner(k,ixy_inner)
        IF (l_2mg) THEN
          dnumber=qfields(k, i_ng)/dt
          procs(i_ng,i_proc%id)%column_data(k)=-dnumber
        END IF
        ! if (l_3mg) then
        !   procs(i_m3g,i_proc%id)%column_data(k)=-qfields(k, i_m3g)/dt
        ! end if
        !--------------------------------------------------
        ! aerosol
        !--------------------------------------------------
        IF (l_process) THEN
          IF (.NOT. am7_reset) THEN
            dmass=dustact(k)%mact3/dt
            dnumber=dustact(k)%nact3/dt
            aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aeroprocs(i_an6,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          END IF
          IF (.NOT. am8_reset) THEN
            dmass=aeroice(k)%mact3/dt
            dnumber=aeroice(k)%nact3/dt
            aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aeroprocs(i_an2,i_aproc%id)%column_data(k)+dnumber
            aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          END IF
        END IF
      END IF
    END IF
    !==============================
    ! Now reset aerosol...
    !==============================

    IF (l_process) THEN
      IF (am4_reset) THEN
        dmass=aerofields(k, i_am4)/dt
        aeroprocs(i_am4,i_aproc%id)%column_data(k)=aeroprocs(i_am4,i_aproc%id)%column_data(k)-dmass
        aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
      END IF

      IF (am5_reset .AND. l_separate_rain) THEN
        dmass=aerofields(k, i_am5)/dt
        aeroprocs(i_am5,i_aproc%id)%column_data(k)=aeroprocs(i_am5,i_aproc%id)%column_data(k)-dmass
        aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
      END IF

      IF (.NOT. l_warm) THEN
        IF (am7_reset) THEN
          dmass=aerofields(k, i_am7)/dt
          aeroprocs(i_am7,i_aproc%id)%column_data(k)=aeroprocs(i_am7,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
        END IF

        IF (am8_reset) THEN
          dmass=aerofields(k, i_am8)/dt
          aeroprocs(i_am8,i_aproc%id)%column_data(k)=aeroprocs(i_am8,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am2,i_aproc%id)%column_data(k)=aeroprocs(i_am2,i_aproc%id)%column_data(k)+dmass
        END IF

        IF (am9_reset) THEN
          dmass=aerofields(k, i_am9)/dt
          aeroprocs(i_am9,i_aproc%id)%column_data(k)=aeroprocs(i_am9,i_aproc%id)%column_data(k)-dmass
          aeroprocs(i_am6,i_aproc%id)%column_data(k)=aeroprocs(i_am6,i_aproc%id)%column_data(k)+dmass
        END IF

        !only do this is passive numbers are used!
        IF ( l_passivenumbers ) THEN        
          IF (an11_reset) THEN  !just reset it and put number into accum sol
            aeroprocs(i_an11,i_aproc%id)%column_data(k)=-aerofields(k,i_an11)/dt
            aeroprocs(i_an2,i_aproc%id)%column_data(k)=aerofields(k,i_an11)/dt
          END IF
        END IF
        IF ( l_passivenumbers_ice ) THEN        
          IF (an12_reset) THEN  !just reset it and put number into coarse insol
            aeroprocs(i_an12,i_aproc%id)%column_data(k)=-aerofields(k,i_an12)/dt
            aeroprocs(i_an6,i_aproc%id)%column_data(k)=aerofields(k,i_an12)/dt
          END IF
        END IF

      END IF
    END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE qtidy

  SUBROUTINE tidy_qin(ixy_inner, qfields, l_negonly)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='TIDY_QIN'

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(INOUT) :: qfields(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: l_negonly

    LOGICAL :: ql_reset, nl_reset, qr_reset, nr_reset, m3r_reset
    LOGICAL :: qi_reset, ni_reset, qs_reset, ns_reset, m3s_reset
    LOGICAL :: qg_reset, ng_reset, m3g_reset
    INTEGER :: k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    IF (present(l_negonly)) THEN
      IF (l_negonly .NEQV. current_qin_l_negonly) THEN
        current_qin_l_negonly=l_negonly
        CALL recompute_qin_constants(l_negonly)
      END IF
    ELSE IF (current_qin_l_negonly) THEN
      current_qin_l_negonly=.FALSE.
      CALL recompute_qin_constants(.FALSE.)
    END IF

    DO k=1, ubound(qfields,1)
      nl_reset=.FALSE.
      nr_reset=.FALSE.
      m3r_reset=.FALSE.
      qi_reset=.FALSE.
      ni_reset=.FALSE.
      qs_reset=.FALSE.
      ns_reset=.FALSE.
      m3s_reset=.FALSE.
      qg_reset=.FALSE.
      ng_reset=.FALSE.
      m3g_reset=.FALSE.

      ql_reset=qfields(k, i_ql) < 0.0 .OR. (qfields(k, i_ql) < qin_thresh(i_ql) .AND. qfields(k, i_ql) >0)
      IF (l_2mc) THEN
        nl_reset=qfields(k, i_nl) < 0.0 .OR. (qfields(k, i_nl) < qin_thresh(i_nl) .AND. qfields(k, i_nl) >0) .OR. &
             (qfields(k, i_nl) > 0.0 .AND. qfields(k, i_ql) <= 0.0)
        ql_reset=ql_reset .OR. (qfields(k, i_ql) > 0.0 .AND. qfields(k, i_nl) <= 0.0)
      END IF
      qr_reset=qfields(k, i_qr) < 0.0 .OR. (qfields(k, i_qr) < qin_thresh(i_qr) .AND. qfields(k, i_qr) > 0.0)

      IF (l_2mr) THEN
        nr_reset=qfields(k, i_nr) < 0.0 .OR. (qfields(k, i_nr) < qin_thresh(i_nr) .AND. qfields(k, i_nr) > 0.0) .OR. &
             (qfields(k, i_nr) > 0.0 .AND. qfields(k, i_qr) <= 0.0)
        qr_reset=qr_reset .OR. (qfields(k, i_qr) > 0.0 .AND. qfields(k, i_nr) <= 0.0)
      END IF

      ! if (l_3mr .and. l_tidym3) then
      !   m3r_reset=qfields(k, i_m3r) < 0.0 .or. (qfields(k, i_m3r) < qin_thresh(i_m3r) .and. qfields(k, i_m3r) >0) .or. &
      !        (qfields(k, i_m3r) > 0.0 .and. (qfields(k, i_qr) <=0.0 .or. qfields(k, i_nr) <=0.0))
      !   qr_reset=qr_reset .or. (qfields(k, i_qr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
      !   nr_reset=nr_reset .or. (qfields(k, i_nr) > 0.0 .and. qfields(k, i_m3r) <= 0.0)
      ! end if

      IF (.NOT. l_warm) THEN
        qi_reset=qfields(k, i_qi) < 0.0 .OR. (qfields(k, i_qi) < qin_thresh(i_qi) .AND. qfields(k, i_qi) > 0.0)
        IF (l_2mi) THEN
          ni_reset=qfields(k, i_ni) < 0.0 .OR. (qfields(k, i_ni) < qin_thresh(i_ni) .AND. qfields(k, i_ni) > 0.0) .OR. &
               (qfields(k, i_ni) > 0.0 .AND. qfields(k, i_qi) <= 0.0)
          qi_reset=qi_reset .OR. (qfields(k, i_qi) > 0.0 .AND. qfields(k, i_ni) <= 0.0)
        END IF

        qs_reset=qfields(k, i_qs) < 0.0 .OR. (qfields(k, i_qs) < qin_thresh(i_qs) .AND. qfields(k, i_qs) > 0.0)
        IF (l_2ms) THEN
          ns_reset=qfields(k, i_ns) < 0.0 .OR. (qfields(k, i_ns) < qin_thresh(i_ns) .AND. qfields(k, i_ns) > 0.0) .OR. &
               (qfields(k, i_ns) > 0.0 .AND. qfields(k, i_qs) <= 0.0)
          qs_reset=qs_reset .OR. (qfields(k, i_qs) > 0.0 .AND. qfields(k, i_ns) <= 0.0)
        END IF
        ! if (l_3ms .and. l_tidym3) then
        !   m3s_reset=qfields(k, i_m3s) < 0.0 .or. (qfields(k, i_m3s) < qin_thresh(i_m3s) .and. qfields(k, i_m3s) >0) .or.&
        !        (qfields(k, i_m3s) > 0.0 .and. (qfields(k, i_qs) <=0.0 .or. qfields(k, i_ns) <=0.0))
        !   qs_reset=qs_reset .or. (qfields(k, i_qs) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
        !   ns_reset=ns_reset .or. (qfields(k, i_ns) > 0.0 .and. qfields(k, i_m3s) <= 0.0)
        ! end if

        qg_reset=qfields(k, i_qg) < 0.0 .OR. (qfields(k, i_qg) < qin_thresh(i_qg) .AND. qfields(k, i_qg) > 0.0)
        IF (l_2mg) THEN
          ng_reset=qfields(k, i_ng) < 0.0 .OR. (qfields(k, i_ng) < qin_thresh(i_ng) .AND. qfields(k, i_ng) > 0.0) .OR. &
               (qfields(k, i_ng) > 0.0 .AND. qfields(k, i_qg) <= 0.0)
          qg_reset=qg_reset .OR. (qfields(k, i_qg) > 0.0 .AND. qfields(k, i_ng) <= 0.0)
        END IF
        ! if (l_3mg .and. l_tidym3) then
        !   m3g_reset=qfields(k, i_m3g) < 0.0 .or. (qfields(k, i_m3g) < qin_thresh(i_m3g) .and. qfields(k, i_m3g) >0) .or.&
        !        (qfields(k, i_m3g) > 0.0 .and. (qfields(k, i_qg) <=0.0 .or. qfields(k, i_ng) <=0.0))
        !   qg_reset=qg_reset .or. (qfields(k, i_qg) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
        !   ng_reset=ng_reset .or. (qfields(k, i_ng) > 0.0 .and. qfields(k, i_m3g) <= 0.0)
        ! end if
      END IF

      !==============================
      ! Now reset things...
      !==============================
      IF (ql_reset .OR. nl_reset) THEN
        IF (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_ql)/exner(k,ixy_inner)
        IF (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_ql)
        qfields(k,i_ql)=0.0
        IF (l_2mc) THEN
          qfields(k,i_nl)=0.0
        END IF
      END IF

      IF (qr_reset .OR. nr_reset .OR. m3r_reset) THEN
        IF (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Lv/cp*qfields(k,i_qr)/exner(k,ixy_inner)
        IF (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qr)
        qfields(k,i_qr)=0.0
        IF (l_2mr) THEN
          qfields(k,i_nr)=0.0
        END IF
        ! if (l_3mr) then
        !   qfields(k,i_m3r)=0.0
        ! end if
      END IF

      IF (qi_reset .OR. ni_reset) THEN
        IF (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qi)/exner(k,ixy_inner)
        IF (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qi)
        qfields(k,i_qi)=0.0
        IF (l_2mi) THEN
          qfields(k,i_ni)=0.0
        END IF
      END IF

      IF (qs_reset .OR. ns_reset .OR. m3s_reset) THEN
        IF (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qs)/exner(k,ixy_inner)
        IF (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qs)
        qfields(k,i_qs)=0.0
        IF (l_2ms) THEN
          qfields(k,i_ns)=0.0
        END IF
        ! if (l_3ms) then
        !   qfields(k,i_m3s)=0.0
        ! end if
      END IF

      IF (qg_reset .OR. ng_reset .OR. m3g_reset) THEN
        IF (l_tidy_conserve_E) qfields(k,i_th)=qfields(k,i_th)-Ls/cp*qfields(k,i_qg)/exner(k,ixy_inner)
        IF (l_tidy_conserve_q) qfields(k,i_qv)=qfields(k,i_qv)+qfields(k,i_qg)
        qfields(k,i_qg)=0.0
        IF (l_2mg) THEN
          qfields(k,i_ng)=0.0
        END IF
        ! if (l_3mg) then
        !   qfields(k,i_m3g)=0.0
        ! end if
      END IF

    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE tidy_qin

  SUBROUTINE tidy_ain(qfields, aerofields)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='TIDY_AIN'

    REAL(wp), INTENT(IN) :: qfields(:,:)
    REAL(wp), INTENT(INOUT) :: aerofields(:,:)

    INTEGER :: k

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k=1, ubound(qfields,1)
      IF ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am4)>0.0) .OR. aerofields(k,i_am4) < 0.0) THEN

        aerofields(k,i_am4)=0.0
      END IF
      IF (i_am9 > 0)THEN
        IF ((qfields(k, i_ql)+qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am9)>0.0) .OR. aerofields(k,i_am9) < 0.0) THEN
          aerofields(k,i_am9)=0.0
        END IF
      END IF

      IF (i_am5 > 0) THEN
        IF (((qfields(k,i_qr) <=0.0 .AND. aerofields(k,i_am5)>0.0) .OR. aerofields(k,i_am5) < 0.0)) THEN
          aerofields(k,i_am5)=0.0
        END IF
      END IF
      IF (i_am7 > 0) THEN
        IF (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .AND. aerofields(k,i_am7)>0.0) &
             .OR. aerofields(k,i_am7) < 0.0)) THEN
          aerofields(k,i_am7)=0.0
        END IF
      END IF
      IF (i_am8 > 0) THEN
        IF (((qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) <=0.0 .AND. aerofields(k,i_am8)>0.0) &
             .OR. aerofields(k,i_am8) < 0.0)) THEN
          aerofields(k,i_am8)=0.0
        END IF
      END IF
      IF (i_an12 > 0) THEN
        IF (((qfields(k, i_ql) + qfields(k,i_qr) + qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) &
                                                                 <=0.0 .AND. aerofields(k,i_an12)>0.0) &
                                                            .OR. aerofields(k,i_an12) < 0.0)) THEN
          aerofields(k,i_an12)=0.0
        END IF
      END IF
      IF (i_an11 > 0) THEN
        IF (((qfields(k, i_ql) + qfields(k,i_qr) + qfields(k,i_qi) + qfields(k,i_qs) + qfields(k,i_qg) &
                                                                 <=0.0 .AND. aerofields(k,i_an12)>0.0) &
                                                            .OR. aerofields(k,i_an11) < 0.0)) THEN
          aerofields(k,i_an11)=0.0
        END IF
      END IF

    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE tidy_ain

  ! Subroutine to ensure parallel processes don't remove more
    ! mass than is available and then rescales all processes
    ! (including number and other terms)
  SUBROUTINE ensure_positive(nz, dt, qfields, procs, params, iprocs_scalable, &
                             iprocs_nonscalable, aeroprocs, iprocs_dependent, &
                             iprocs_dependent_ns)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ENSURE_POSITIVE'

    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:)         ! microphysical process rates

    TYPE(hydro_params), INTENT(IN) :: params        ! parameters from hydrometeor variable to test
    TYPE(process_name), INTENT(IN) :: iprocs_scalable(:)    ! list of processes to rescale
    TYPE(process_name), INTENT(IN), OPTIONAL ::      &
         iprocs_nonscalable(:) ! list of other processes which
    ! provide source or sink, but
    ! which we don't want to rescale
    TYPE(process_rate), INTENT(INOUT), OPTIONAL ::   &
         aeroprocs(:,:)        ! associated aerosol process rates
    TYPE(process_name), INTENT(IN), OPTIONAL ::      &
         iprocs_dependent(:)   ! list of aerosol processes which
    ! are dependent on rescaled processes and
    ! so should be rescaled themselves
    TYPE(process_name), INTENT(IN), OPTIONAL ::      &
         iprocs_dependent_ns(:)   ! list of aerosol processes which
    ! are dependent on rescaled processes but
    ! we don't want to rescale
    INTEGER :: iproc, id, iq, k
    REAL(wp) :: delta_scalable, delta_nonscalable, ratio, maxratio
    INTEGER :: i_1m, i_2m, i_3m
    LOGICAL :: l_rescaled

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
    
    DO k=1,nz
    i_1m=params%i_1m
    IF (params%l_2m) i_2m=params%i_2m
    IF (params%l_3m) i_3m=params%i_3m
    l_rescaled=.FALSE.

    IF (qfields(k, i_1m) > 0.0) THEN
      ! First calculate the unscaled increments
      delta_scalable=0.0
      delta_nonscalable=0.0
      DO iproc=1,size(iprocs_scalable)
        IF (iprocs_scalable(iproc)%on) THEN
          id=iprocs_scalable(iproc)%id
          delta_scalable=delta_scalable+procs(i_1m, id)%column_data(k)
        END IF
      END DO
      delta_scalable=delta_scalable*dt
      IF (present(iprocs_nonscalable)) THEN
        DO iproc=1, size(iprocs_nonscalable)
          IF (iprocs_nonscalable(iproc)%on) THEN
            id=iprocs_nonscalable(iproc)%id
            delta_nonscalable=delta_nonscalable+procs(i_1m, id)%column_data(k)
          END IF
        END DO
      END IF
      delta_nonscalable=delta_nonscalable*dt

      ! Test to see if we need to do anything
      IF (delta_scalable+delta_nonscalable+qfields(k, i_1m) < spacing(qfields(k, i_1m)) &
           .AND. abs(delta_scalable) > epsilon(delta_scalable)) THEN
        ratio=-(qfields(k, i_1m)+delta_nonscalable)/delta_scalable
        IF (ratio > 1.0) THEN
          IF (present(iprocs_nonscalable)) THEN
            DO iproc=1, size(iprocs_nonscalable)
              IF (iprocs_nonscalable(iproc)%on) THEN
                id=iprocs_nonscalable(iproc)%id
              END IF
            END DO
          END IF
          DO iproc=1, size(iprocs_scalable)
            IF (iprocs_scalable(iproc)%on) THEN
              id=iprocs_scalable(iproc)%id
            END IF
          END DO

          DO iproc=1, size(iprocs_scalable)
            IF (iprocs_scalable(iproc)%on) THEN
              id=iprocs_scalable(iproc)%id
            END IF
          END DO
          DO iproc=1, size(iprocs_nonscalable)
            IF (iprocs_nonscalable(iproc)%on) THEN
              id=iprocs_nonscalable(iproc)%id
            END IF
          END DO

          WRITE(std_msg, '(A, F7.4)') 'Problem with ratio > 1.0: ratio = ', ratio
          ! Tell mphys_error that this is due to bad values
          CALL throw_mphys_error(bad_values, ModuleName//':'//RoutineName, std_msg)

        END IF
        ! Now rescale the scalable processes
        DO iproc=1, size(iprocs_scalable)
          IF (iprocs_scalable(iproc)%on) THEN
            id=iprocs_scalable(iproc)%id
            DO iq = i_qstart, i_nstart-1
               procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
            END DO
          END IF
        END DO
        ! Set flag to indicate a rescaling was performed
        l_rescaled=.TRUE.
      END IF

      ! Now we need to rescale additional moments
      IF (l_rescaled) THEN
        IF (params%l_2m) THEN ! second moment
          ! First calculate the unscaled increments
          delta_scalable=0.0
          delta_nonscalable=0.0
          DO iproc=1,size(iprocs_scalable)
            IF (iprocs_scalable(iproc)%on) THEN
              id=iprocs_scalable(iproc)%id
              delta_scalable=delta_scalable+procs(i_2m, id)%column_data(k)
            END IF
          END DO
          delta_scalable=delta_scalable*dt
          IF (abs(delta_scalable) > epsilon(delta_scalable)) THEN
            IF (present(iprocs_nonscalable)) THEN
              DO iproc=1, size(iprocs_nonscalable)
                IF (iprocs_nonscalable(iproc)%on) THEN
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(i_2m, id)%column_data(k)
                END IF
              END DO
            END IF
            delta_nonscalable=delta_nonscalable*dt
            ! ratio may now be greater than 1
            ratio=-(qfields(k, i_2m)+delta_nonscalable)/delta_scalable
            ! Now rescale the scalable processes
            DO iproc=1, size(iprocs_scalable)
              IF (iprocs_scalable(iproc)%on) THEN
                id=iprocs_scalable(iproc)%id
                DO iq = i_nstart, i_m3start-1
                   procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
                END DO
              END IF
            END DO
          END IF
        END IF
        ! ! if (params%l_3m) then ! third moment
        ! !   ! First calculate the unscaled increments
        ! !   delta_scalable=0.0
        ! !   delta_nonscalable=0.0
        ! !   do iproc=1,size(iprocs_scalable)
        ! !     if (iprocs_scalable(iproc)%on) then
        ! !       id=iprocs_scalable(iproc)%id
        ! !       delta_scalable=delta_scalable+procs(i_3m, id)%column_data(k)
        ! !     end if
        ! !   end do
        ! !   delta_scalable=delta_scalable*dt
        ! !   if (abs(delta_scalable) > epsilon(delta_scalable)) then
        ! !     if (present(iprocs_nonscalable)) then
        ! !       do iproc=1, size(iprocs_nonscalable)
        ! !         if (iprocs_nonscalable(iproc)%on) then
        ! !           id=iprocs_nonscalable(iproc)%id
        ! !           delta_nonscalable=delta_nonscalable+procs(i_3m, id)%column_data(k)
        ! !         end if
        ! !       end do
        ! !     end if
        ! !     delta_nonscalable=delta_nonscalable*dt
        ! !     ! ratio may now be greater than 1
        ! !     ratio=-(qfields(k, i_3m)+delta_nonscalable)/(delta_scalable + epsilon(1.0))
        ! !     ! Now rescale the scalable processes
        ! !     do iproc=1, size(iprocs_scalable)
        ! !       if (iprocs_scalable(iproc)%on) then
        ! !         id=iprocs_scalable(iproc)%id
        ! !         do iq = i_m3start, ubound(procs,1)
        ! !            procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
        ! !         enddo
        ! !       end if
        ! !     end do
        ! !   end if
        ! ! end if

        !Now rescale the increments to aerosol
        ! How do we do this?????
        !        print*, 'WARNING: Should be rescaling aerosol?, but not done!'

      ELSE ! What if we haven't rescaled mass, but number is now not conserved?
        IF (l_rescale_on_number) THEN
          IF (params%l_2m) THEN ! second moment
            ! First calculate the unscaled increments
            delta_scalable=0.0
            delta_nonscalable=0.0
            DO iproc=1, size(iprocs_scalable)
              IF (iprocs_scalable(iproc)%on) THEN
                id=iprocs_scalable(iproc)%id
                delta_scalable=delta_scalable+procs(i_2m, id)%column_data(k)
              END IF
            END DO
            delta_scalable=delta_scalable*dt
            IF (present(iprocs_nonscalable)) THEN
              DO iproc=1, size(iprocs_nonscalable)
                IF (iprocs_nonscalable(iproc)%on) THEN
                  id=iprocs_nonscalable(iproc)%id
                  delta_nonscalable=delta_nonscalable+procs(i_2m, id)%column_data(k)
                END IF
              END DO
            END IF
            delta_nonscalable=delta_nonscalable*dt
            ! Test to see if we need to do anything
            IF (delta_scalable+delta_nonscalable+qfields(k, i_2m) < spacing(qfields(k, i_2m)) &
                 .AND. abs(delta_scalable) > epsilon(delta_scalable) ) THEN
              maxratio=(1.0-spacing(delta_scalable))
              IF (abs(delta_scalable) < spacing(delta_scalable)) THEN
                IF (present(iprocs_nonscalable)) THEN
                  DO iproc=1, size(iprocs_nonscalable)
                    IF (iprocs_nonscalable(iproc)%on) THEN
                      id=iprocs_nonscalable(iproc)%id
                    END IF
                  END DO
                END IF
                DO iproc=1, size(iprocs_scalable)
                  IF (iprocs_scalable(iproc)%on) THEN
                    id=iprocs_scalable(iproc)%id
                  END IF
                END DO
              END IF
              ratio=-(maxratio*qfields(k, i_2m)+delta_nonscalable)/delta_scalable
              ratio=max(ratio, 0.0_wp)
              IF (ratio==0.0_wp) THEN
                IF (present(iprocs_nonscalable)) THEN
                  DO iproc=1, size(iprocs_nonscalable)
                    IF (iprocs_nonscalable(iproc)%on) id=iprocs_nonscalable(iproc)%id
                  END DO
                END IF
                DO iproc=1, size(iprocs_scalable)
                  IF (iprocs_scalable(iproc)%on) id=iprocs_scalable(iproc)%id
                END DO
              END IF

              IF (ratio<0.95) THEN
                ! Some warnings for testing

                WRITE(std_msg, *) 'WARNING: Significantly rescaled number, but not sure ' // &
                                   'what to do with other moments. id, ratio, bad',           &
                                    params%id, ratio, delta_scalable + delta_nonscalable +    &
                                    qfields(k, i_2m), spacing(qfields(k, i_2m)),              &
                                    (qfields(k, i_2m) + delta_nonscalable), delta_scalable,   &
                                    'qfields', qfields(k, i_1m), qfields(k, i_2m)

                CALL throw_mphys_error(warn, ModuleName//':'//RoutineName, std_msg)


                IF (present(iprocs_nonscalable)) THEN
                  DO iproc=1, size(iprocs_nonscalable)
                    IF (iprocs_nonscalable(iproc)%on) THEN
                      id=iprocs_nonscalable(iproc)%id
                    END IF
                  END DO
                END IF
                DO iproc=1, size(iprocs_scalable)
                  IF (iprocs_scalable(iproc)%on) THEN
                    id=iprocs_scalable(iproc)%id
                  END IF
                END DO
              END IF

              ! Now rescale the scalable processes
              DO iproc=1, size(iprocs_scalable)
                IF (iprocs_scalable(iproc)%on) THEN
                  id=iprocs_scalable(iproc)%id
                  IF (i_2m > 0) THEN
                     DO iq = i_nstart, i_m3start-1
                        procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
                     END DO
                  END IF
                END IF
              END DO

              ! Set flag to indicate a rescaling was performed
              l_rescaled=.TRUE.
            END IF
          END IF
        END IF
      END IF
    END IF
    END DO
  
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ensure_positive

  ! Subroutine to ensure parallel aerosol processes don't remove more
  ! mass than is available and then rescales all processes
  ! (NB this follows any rescaling due to the parent microphysical processes and
  ! we might lose consistency between number and mass here)
  SUBROUTINE ensure_positive_aerosol(nz, dt, aerofields, aerosol_procs, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ENSURE_POSITIVE_AEROSOL'

    INTEGER, INTENT(IN) :: nz
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: aerofields(:,:)
    TYPE(process_rate), INTENT(INOUT) :: aerosol_procs(:,:)  ! aerosol process rates
    TYPE(process_name), INTENT(IN) :: iprocs(:)    ! list of processes to rescale

    INTEGER :: iq, iproc, id, k
    REAL(wp) :: ratio, delta_scalable

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k=1,nz
      DO iq=1, ntotala
        delta_scalable=0.0
        DO iproc=1, size(iprocs)
          IF (iprocs(iproc)%on) THEN
            id=iprocs(iproc)%id
            delta_scalable=delta_scalable + aerosol_procs(iq, id)%column_data(k)
          END IF
        END DO
        delta_scalable=delta_scalable*dt
        IF (delta_scalable + aerofields(k, iq) < spacing(aerofields(k, iq))    &
             .AND. abs(delta_scalable) > spacing(aerofields(k, iq))) THEN
          ratio=(spacing(aerofields(k, iq))-aerofields(k, iq))/(delta_scalable)

          DO iproc=1 ,size(iprocs)
            IF (iprocs(iproc)%on) THEN
              id=iprocs(iproc)%id
              aerosol_procs(iq,id)%column_data(k)=aerosol_procs(iq,id)%column_data(k)*ratio
            END IF
          END DO
        END IF
      END DO
    END DO
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ensure_positive_aerosol

  ! Subroutine to ensure parallel ice processes don't remove more
    ! vapour than is available (i.e. so become subsaturated)
    ! and then rescales processes

    ! Modified so it can also prevent sublimation processes from putting
    ! back too much vapour and so become supersaturated
  SUBROUTINE ensure_saturated(ixy_inner, nz, l_Tcold, dt, qfields, procs, iprocs_scalable)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ENSURE_SATURATED'

    INTEGER, INTENT(IN) :: ixy_inner
    INTEGER, INTENT(IN) :: nz
    LOGICAL, INTENT(IN) :: l_Tcold(:) 
    REAL(wp), INTENT(IN) :: dt
    REAL(wp), INTENT(IN) :: qfields(:,:)
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:)         ! microphysical process rates
    TYPE(process_name), INTENT(IN) :: iprocs_scalable(:)    ! list of processes to rescale

    INTEGER :: iproc, id, iq, k
    REAL(wp) :: delta_scalable, ratio, delta_sat
    REAL(wp) :: th, qis

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO k=1, nz
    IF (l_Tcold(k)) THEN
   
    delta_scalable=0.0
    DO iproc=1, size(iprocs_scalable)
      IF (iprocs_scalable(iproc)%on) THEN
        id=iprocs_scalable(iproc)%id
        delta_scalable=delta_scalable+procs(i_qv, id)%column_data(k)*dt
      END IF
    END DO

    delta_scalable=abs(delta_scalable)

    IF (delta_scalable > spacing(delta_scalable)) THEN
      th=qfields(k, i_th)

      qis=qisaturation(th*exner(k,ixy_inner), pressure(k,ixy_inner)/100.0)

      delta_sat=abs(qis-qfields(k, i_qv))

      IF (delta_scalable > delta_sat) THEN
        ratio=delta_sat/delta_scalable
        DO iproc=1, size(iprocs_scalable)
          IF (iprocs_scalable(iproc)%on) THEN
            id=iprocs_scalable(iproc)%id
            DO iq = i_qstart, i_nstart-1
               procs(iq, id)%column_data(k)=procs(iq, id)%column_data(k)*ratio
            END DO
          END IF
        END DO
      END IF
    END IF
    END IF
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE ensure_saturated

END MODULE mphys_tidy
