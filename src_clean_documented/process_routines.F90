! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Process-rate bookkeeping: allocation/deallocation/zeroing of the procs and aerosol_procs arrays.
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Method:
!   Defines the process_rate derived type and the master list
!   of process-rate indices (i_praut, i_pracw, ...) shared by
!   every physics-process module to accumulate tendencies.
!
! Paper reference:
!   Provides the software infrastructure for the process-rate
!   subscript naming convention described in Field et al. (2023)
!   Table A2; it is CASIM infrastructure, not itself the subject
!   of a numbered equation.
!
MODULE process_routines
  USE variable_precision, ONLY: wp
  USE type_process, ONLY: process_name, process_rate
  USE mphys_parameters, ONLY: hydro_params, ZERO_REAL_WP

  IMPLICIT NONE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='PROCESS_ROUTINES'

  ! NB These ids are now overwritten in mphys_switches to only allocate those
  ! that are needed given namelist switches
  TYPE(process_name) :: i_cond  = process_name(0, 1, 'pcond', on=.FALSE.)
  TYPE(process_name) :: i_praut = process_name(0, 2, 'praut', on=.FALSE.)
  TYPE(process_name) :: i_pracw = process_name(0, 3, 'pracw', on=.FALSE.)
  TYPE(process_name) :: i_pracr = process_name(0, 4, 'pracr', on=.FALSE.)
  TYPE(process_name) :: i_prevp = process_name(0, 5, 'prevp', on=.FALSE.)
  TYPE(process_name) :: i_psedl = process_name(0, 6, 'psedl', on=.FALSE.)
  TYPE(process_name) :: i_psedr = process_name(0, 7, 'psedr', on=.FALSE.)
  TYPE(process_name) :: i_tidy  = process_name(0, 8, 'ptidy', on=.FALSE.)
  TYPE(process_name) :: i_tidy2 = process_name(0, 9, 'ptidy2', on=.FALSE.)
  TYPE(process_name) :: i_inuc  = process_name(0, 10, 'pinuc', on=.FALSE.)
  TYPE(process_name) :: i_idep  = process_name(0, 11, 'pidep', on=.FALSE.)
  TYPE(process_name) :: i_iacw  = process_name(0, 12, 'piacw', on=.FALSE.)
  TYPE(process_name) :: i_saut  = process_name(0, 13, 'psaut', on=.FALSE.)
  TYPE(process_name) :: i_sdep  = process_name(0, 14, 'psdep', on=.FALSE.)
  TYPE(process_name) :: i_sacw  = process_name(0, 15, 'psacw', on=.FALSE.)
  TYPE(process_name) :: i_gdep  = process_name(0, 16, 'pgdep', on=.FALSE.)
  TYPE(process_name) :: i_pseds = process_name(0, 17, 'pseds', on=.FALSE.)
  TYPE(process_name) :: i_psedi = process_name(0, 18, 'psedi', on=.FALSE.)
  TYPE(process_name) :: i_psedg = process_name(0, 19, 'psedg', on=.FALSE.)
  TYPE(process_name) :: i_saci  = process_name(0, 20, 'psaci', on=.FALSE.)
  TYPE(process_name) :: i_raci  = process_name(0, 21, 'praci', on=.FALSE.)
  TYPE(process_name) :: i_sacr  = process_name(0, 22, 'psacr', on=.FALSE.)
  TYPE(process_name) :: i_gacr  = process_name(0, 23, 'pgacr', on=.FALSE.)
  TYPE(process_name) :: i_gacw  = process_name(0, 24, 'pgacw', on=.FALSE.)
  TYPE(process_name) :: i_gaci  = process_name(0, 25, 'pgaci', on=.FALSE.)
  TYPE(process_name) :: i_gacs  = process_name(0, 26, 'pgacs', on=.FALSE.)
  TYPE(process_name) :: i_iagg  = process_name(0, 27, 'piagg', on=.FALSE.)
  TYPE(process_name) :: i_sagg  = process_name(0, 28, 'psagg', on=.FALSE.)
  TYPE(process_name) :: i_gagg  = process_name(0, 29, 'pgagg', on=.FALSE.)
  TYPE(process_name) :: i_sbrk  = process_name(0, 30, 'psbrk', on=.FALSE.)
  TYPE(process_name) :: i_gshd  = process_name(0, 31, 'pgshd', on=.FALSE.)
  TYPE(process_name) :: i_ihal  = process_name(0, 32, 'pihal', on=.FALSE.)
  TYPE(process_name) :: i_smlt  = process_name(0, 33, 'psmlt', on=.FALSE.)
  TYPE(process_name) :: i_gmlt  = process_name(0, 34, 'pgmlt', on=.FALSE.)
  TYPE(process_name) :: i_homr  = process_name(0, 35, 'phomr', on=.FALSE.)
  TYPE(process_name) :: i_homc  = process_name(0, 36, 'phomc', on=.FALSE.)
  TYPE(process_name) :: i_ssub  = process_name(0, 37, 'pssub', on=.FALSE.)
  TYPE(process_name) :: i_gsub  = process_name(0, 38, 'pgsub', on=.FALSE.)
  TYPE(process_name) :: i_isub  = process_name(0, 39, 'pisub', on=.FALSE.)
  TYPE(process_name) :: i_imlt  = process_name(0, 40, 'pimlt', on=.FALSE.)
  TYPE(process_name) :: i_iics  = process_name(0, 41, 'piics', on=.FALSE.)
  TYPE(process_name) :: i_idps  = process_name(0, 42, 'pidps', on=.FALSE.)
  ! aerosol processes
  TYPE(process_name)  :: i_aact  = process_name(0, 101, 'aact', on=.FALSE.)
  TYPE(process_name)  :: i_aaut  = process_name(0, 102, 'aaut', on=.FALSE.)
  TYPE(process_name)  :: i_aacw  = process_name(0, 103, 'aacw', on=.FALSE.)
  TYPE(process_name)  :: i_aevp  = process_name(0, 104, 'aevp', on=.FALSE.)
  TYPE(process_name)  :: i_asedr = process_name(0, 105, 'asedr', on=.FALSE.)
  TYPE(process_name)  :: i_arevp = process_name(0, 106, 'arevp', on=.FALSE.)
  TYPE(process_name)  :: i_asedl = process_name(0, 107, 'asedl', on=.FALSE.)
  !... additional tidying processes (Need to sort out location for these)
  TYPE(process_name)  :: i_atidy = process_name(0, 108, 'atidy', on=.FALSE.)
  TYPE(process_name)  :: i_atidy2 = process_name(0, 109, 'atidy2', on=.FALSE.)
  !... ice related processes
  TYPE(process_name)  :: i_dnuc  = process_name(0, 110, 'dnuc', on=.FALSE.)
  TYPE(process_name)  :: i_dsub  = process_name(0, 111, 'dsub', on=.FALSE.)
  TYPE(process_name)  :: i_dsedi = process_name(0, 112, 'dsedi', on=.FALSE.)
  TYPE(process_name)  :: i_dseds = process_name(0, 113, 'dseds', on=.FALSE.)
  TYPE(process_name)  :: i_dsedg = process_name(0, 114, 'dsedg', on=.FALSE.)
  TYPE(process_name)  :: i_dssub  = process_name(0, 115, 'dssub', on=.FALSE.)
  TYPE(process_name)  :: i_dgsub  = process_name(0, 116, 'dgsub', on=.FALSE.)
  TYPE(process_name)  :: i_dhomc  = process_name(0, 117, 'dhomc', on=.FALSE.)
  TYPE(process_name)  :: i_dhomr  = process_name(0, 118, 'dhomr', on=.FALSE.)
  TYPE(process_name)  :: i_dimlt  = process_name(0, 119, 'dimlt', on=.FALSE.)
  TYPE(process_name)  :: i_dsmlt  = process_name(0, 120, 'dsmlt', on=.FALSE.)
  TYPE(process_name)  :: i_dgmlt  = process_name(0, 121, 'dgmlt', on=.FALSE.)
  TYPE(process_name)  :: i_diacw  = process_name(0, 122, 'diacw', on=.FALSE.)
  TYPE(process_name)  :: i_dsacw  = process_name(0, 123, 'dsacw', on=.FALSE.)
  TYPE(process_name)  :: i_dgacw  = process_name(0, 124, 'dgacw', on=.FALSE.)
  TYPE(process_name)  :: i_dsacr  = process_name(0, 125, 'dsacr', on=.FALSE.)
  TYPE(process_name)  :: i_dgacr  = process_name(0, 126, 'dgacr', on=.FALSE.)
  TYPE(process_name)  :: i_draci  = process_name(0, 127, 'draci', on=.FALSE.)

CONTAINS

  ! Allocate space to store the microphysical process rates
  SUBROUTINE allocate_procs(nxy_inner, procs, nz, nprocs, ntotalq)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ALLOCATE_PROCS'

    INTEGER, INTENT(IN) :: nxy_inner
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:,:)
    INTEGER, INTENT(IN) :: nz      ! number of height levels
    INTEGER, INTENT(IN) :: nprocs  ! number of physical processes
    INTEGER, INTENT(IN) :: ntotalq ! number of q or aerosol fields

    INTEGER :: iproc, iq, ixy

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO ixy=1, nxy_inner
     DO iproc=1, nprocs
       DO iq=1,ntotalq
         ALLOCATE(procs(iq,iproc,ixy)%column_data(nz))
       END DO
       END DO
     END DO
    
    DO ixy=1, nxy_inner
       CALL zero_procs(procs(:,:,ixy))
    END DO

    !  call zero_procs(procs)

     IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

   END SUBROUTINE allocate_procs

  SUBROUTINE zero_procs(procs, iprocs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='ZERO_PROCS'

    TYPE(process_rate), INTENT(INOUT) :: procs(:,:)
    TYPE(process_name), INTENT(IN), OPTIONAL :: iprocs(:)

    INTEGER :: iproc, nproc, iq, i
    INTEGER :: lb1, lb2, ub1, ub2

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    lb1=lbound(procs,1)
    ub1=ubound(procs,1)
    
    IF(present(iprocs)) THEN 
       nproc = size(iprocs)
       DO i = 1, nproc
          iproc = iprocs(i)%id
          DO iq=lb1,ub1
             procs(iq,iproc)%column_data(:)=0.0
          END DO
       END DO
    ELSE
       lb2=lbound(procs,2)
       ub2=ubound(procs,2)
       DO iproc=lb2, ub2
          DO iq=lb1,ub1
             procs(iq,iproc)%column_data(:)=0.0
          END DO
       END DO
    END IF
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE zero_procs
 

  SUBROUTINE deallocate_procs(nxy_inner, procs)

    USE yomhook, ONLY: lhook, dr_hook
    USE parkind1, ONLY: jprb, jpim

    IMPLICIT NONE

    CHARACTER(len=*), PARAMETER :: RoutineName='DEALLOCATE_PROCS'

    INTEGER, INTENT(IN) :: nxy_inner
    TYPE(process_rate), INTENT(INOUT) :: procs(:,:,:)

    INTEGER :: iproc, iq, ixy

    INTEGER(KIND=jpim), PARAMETER :: zhook_in  = 0
    INTEGER(KIND=jpim), PARAMETER :: zhook_out = 1
    REAL(KIND=jprb)               :: zhook_handle

    !--------------------------------------------------------------------------
    ! End of header, no more declarations beyond here
    !--------------------------------------------------------------------------
    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

    DO ixy=1, nxy_inner
    DO iproc=lbound(procs,2), ubound(procs,2)
      DO iq=lbound(procs,1), ubound(procs,1)
             IF (allocated(procs(iq,iproc,ixy)%column_data)) DEALLOCATE(procs(iq,iproc,ixy)%column_data)
          END DO
      END DO
    END DO

    IF (lhook) CALL dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

  END SUBROUTINE deallocate_procs
END MODULE process_routines
