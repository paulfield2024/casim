! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Accumulates hydrometeor and aerosol process rates into tendencies, including latent heating (sum_procs, sum_aprocs).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Sums the individual process-rate contributions (module
!   sum_process) into net mass/number tendencies per
!   hydrometeor, applied to update the prognostic fields each
!   (sub)step.
!
! Paper reference:
!   Implements the right-hand-side summation of the tendency
!   equations Dqv/Dt ... Dng/Dt given in Field et al. (2023)
!   Eq. (A3); the summation code itself is CASIM infrastructure.
!
MODULE sum_process
  USE variable_precision, ONLY: wp
! use mphys_die, only: throw_mphys_error, incorrect_opt, std_msg
  USE type_process, ONLY: process_name, process_rate
  USE mphys_switches, ONLY: i_th, i_ql, i_qr, i_qs, i_qi, i_qg, l_warm, ntotalq, ntotala
  USE passive_fields, ONLY: rexner
! use passive_fields, only: rho
  USE mphys_constants, ONLY: cp, Lv, Ls
  USE mphys_parameters, ONLY: ZERO_REAL_WP
! use mphys_parameters, only: hydro_params, snow_params, rain_params, graupel_params
! use m3_incs, only: m3_inc_type2
  USE process_routines, ONLY: process_name
  
  IMPLICIT NONE
  PRIVATE

  ! allocated and deallocated in micro_main (initialise and finalise_micromain)
  REAL(wp), ALLOCATABLE :: tend_temp(:,:) ! Temporary storage for accumulated tendendies
  REAL(wp), ALLOCATABLE :: aerosol_tend_temp(:,:) ! Temporary storage for accumulated aerosol tendendies

!$OMP THREADPRIVATE(tend_temp, aerosol_tend_temp)

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='SUM_PROCESS'

  PUBLIC sum_aprocs, sum_procs,  tend_temp, aerosol_tend_temp
CONTAINS

  SUBROUTINE sum_aprocs(dst, nz, procs, tend, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SUM_APROCS'

    REAL(wp), INTENT(IN) :: dst  ! step length (s)
    INTEGER, INTENT(IN) :: nz ! number of points in a column
    TYPE(process_rate), INTENT(IN) :: procs(:,:)
    TYPE(process_name), INTENT(IN) :: iprocs(:)
    REAL(wp), INTENT(INOUT) :: tend(:,:)

    !real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies

    INTEGER :: k, iq, iproc, i
    INTEGER :: nproc

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    !allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    aerosol_tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)

    DO i=1, nproc
      IF (iprocs(i)%on) THEN
        iproc=iprocs(i)%id
        DO iq=1, ntotala
           DO k=1,nz
         ! if (.not. all(procs(k,iproc)%source(:)==ZERO_REAL_WP)) then
              aerosol_tend_temp(k, iq)=aerosol_tend_temp(k, iq) + &
                   procs(iq,iproc)%column_data(k)*dst
            END DO
         ! end if
        END DO
      END IF
    END DO

    ! Add on tendencies to those already passed in.
    tend=tend+aerosol_tend_temp
    !deallocate(tend_temp)

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE sum_aprocs

  SUBROUTINE sum_procs(ixy_inner, dst, nz, procs, tend, iprocs, l_thermalexchange, i_thirdmoment, qfields, l_passive )

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='SUM_PROCS'

    INTEGER, INTENT(IN) :: ixy_inner
    REAL(wp), INTENT(IN) :: dst  ! step length (s)
    INTEGER, INTENT(IN) :: nz ! number of points in a column
    TYPE(process_rate), INTENT(IN) :: procs(:,:)
    TYPE(process_name), INTENT(IN) :: iprocs(:)
    REAL(wp), INTENT(INOUT) :: tend(:,:)
    LOGICAL, INTENT(IN), OPTIONAL :: l_thermalexchange  ! Calculate the thermal exchange terms
    INTEGER, INTENT(IN), OPTIONAL :: i_thirdmoment  ! Calculate the tendency of the third moment
    REAL(wp), INTENT(IN), OPTIONAL :: qfields(:,:) ! Required for debugging or with i_thirdmoment
    LOGICAL, INTENT(IN), OPTIONAL :: l_passive ! If true don't apply final tendency (testing diagnostics with pure sedimentation)

    !real(wp), allocatable :: tend_temp(:,:) ! Temporary storage for accumulated tendendies
!   type(hydro_params) :: params
!   real(wp) :: dm1,dm2,dm3,m1,m2,m3
    INTEGER :: k, iq, iproc, i
    INTEGER :: nproc

    LOGICAL :: do_thermal
!   logical :: do_third  ! Currently not plumbed in, so explicitly calculated for each process
    LOGICAL :: do_update ! update the tendency
!   integer :: third_type

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    do_thermal=.FALSE.
    IF (present(l_thermalexchange)) do_thermal=l_thermalexchange

    ! do_third=.false.
    ! third_type = 0 ! Set up a default value
    ! if (present(i_thirdmoment)) then
    !   do_third=.true.
    !   third_type=i_thirdmoment
    ! end if

    do_update=.TRUE.
    IF (present(l_passive)) do_update= .NOT. l_passive

    ! this allocation is based on the shape of tend, which is (k, proc) originally.
    ! switch tend_temp to be the same. Should this allocation be done here? 
    !allocate(tend_temp(lbound(tend,1):ubound(tend,1), lbound(tend,2):ubound(tend,2)))
    tend_temp=ZERO_REAL_WP

    nproc=size(iprocs)

    DO i=1, nproc
       iproc=iprocs(i)%id
        DO iq=1, ntotalq
          DO k = 1, nz
             tend_temp(k, iq)=tend_temp(k, iq)+procs(iq,iproc)%column_data(k)*dst
          END DO
       END DO
    END DO
  
    ! if (do_third) then
    !    ! calculate increment to third moment based on collected increments from
    !    ! q and n NB This overwrites any previously calculated values
    !    ! Rain
    !    params=rain_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k, params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'rain i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg )
    !                end select
    !                tend_temp(k, params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if
       
    !    ! Snow
    !    params=snow_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k, params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k, params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3> 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'snow i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg)
    !                end select
    !                tend_temp(k, params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if
       
    !    ! Graupel
    !    params=graupel_params
    !    if (params%l_3m) then
    !       do k = 1, nz
    !          m1=qfields(k, params%i_1m)*rho(k)/params%c_x
    !          m2=qfields(k, params%i_2m)
    !          m3=qfields(k, params%i_3m)
    !          dm1=tend_temp(k,params%i_1m)*rho(k)/params%c_x
    !          dm2=tend_temp(k,params%i_2m)
    !          if (dm1 < -.99*m1 .or. dm2 < -.99*m2) then
    !             dm1=-m1
    !             dm2=-m2
    !             dm3=-m3
    !          else
    !             if (m3 > 0.0 .and. m1 > 0.0 .and. m2 > 0.0 .and. (abs(dm1) > 0.0 .or. abs(dm2) > 0.0)) then
    !                select case (third_type)
    !                case (2)
    !                   call m3_inc_type2(m1, m2, m3, params%p1, params%p2, params%p3, dm1, dm2, dm3)
    !                case default
    !                   write(std_msg, '(A)') 'graupel i_thirdmoment incorrectly set'
    !                   call throw_mphys_error(incorrect_opt, ModuleName//':'//RoutineName, &
    !                        std_msg)
    !                end select
    !                tend_temp(k,params%i_3m)=dm3
    !             end if
    !          end if
    !       enddo
    !    end if

    ! endif ! endif for l_dothird
    
    ! Calculate the thermal exchange values
    ! (this overwrites anything that was already stored in the theta tendency)
    IF (do_thermal) THEN
       DO k=1,nz
          tend_temp(k, i_th)=(tend_temp(k, i_ql)+tend_temp(k,i_qr))*Lv/cp * rexner(k,ixy_inner)
       END DO
       IF (.NOT. l_warm) THEN 
          DO k=1, nz
             tend_temp(k, i_th)=tend_temp(k,i_th)+ &
               (tend_temp(k, i_qi)+tend_temp(k, i_qs)+tend_temp(k,i_qg))*Ls/cp *rexner(k,ixy_inner)
          END DO
       END IF
    END IF

    IF (do_update) THEN
      tend = tend+tend_temp
    END IF

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE sum_procs
END MODULE sum_process
