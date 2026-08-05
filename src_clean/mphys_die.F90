! *****************************COPYRIGHT*******************************
! (C) British Crown Copyright 2014-2017, Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Error/warning reporting and termination for CASIM (throw_mphys_error, mphys_message).
!
! Code Owner: Please refer to the CASIM COPYRIGHT.txt and README.md
! This file belongs in section: CASIM microphysics
!
! Code description:
!   Language: Fortran 2008.
!   This code follows the Met Office UMDP3 style guide.
!
! Routine to exit the microphysics following an error there
MODULE mphys_die

  ! Module for the KiD and MONC models - UM has its own version

! #if DEF_MODEL==MODEL_KiD
!   use runtime, only: time
! #endif

  IMPLICIT NONE
  PRIVATE

  CHARACTER(len=*), PARAMETER, PRIVATE :: ModuleName='MPHYS_DIE'

  INTEGER, PARAMETER :: incorrect_opt = 1
  INTEGER, PARAMETER :: bad_values    = 2
  INTEGER, PARAMETER :: warn          = -1
  INTEGER, PARAMETER :: std_msg_len   = 400 ! Length of a standard message

  CHARACTER(len=std_msg_len) :: std_msg = ''

  PUBLIC throw_mphys_error, incorrect_opt, bad_values, warn, std_msg, mphys_message
  
CONTAINS

  SUBROUTINE throw_mphys_error(itype, routine, info)

    ! If modifying the subroutine or argument list, ensure that the 
    ! UM version of mphys_die is also modified to give the same answers

    IMPLICIT NONE

    INTEGER,      INTENT(IN) :: itype  ! type of error 1 = Incorrect specification of options
                                       !               2 = Bad values found
                                       !               3 = Unknown error
                                       !              <0 = Warning, code will continue
    CHARACTER(*), INTENT(IN) :: routine
    CHARACTER(*), INTENT(IN) :: info ! error information

    INTEGER, PARAMETER :: nstandard_types=3
    CHARACTER(100) :: stdinfo(nstandard_types) =     &
         (/ 'Incorrect specification of options      ' &
         ,  'Bad values found                        ' &
         ,  'Unknown error                           ' &
         /)
    CHARACTER(1000) :: str
    REAL :: minus_one=-1.

    CHARACTER(len=*), PARAMETER :: RoutineName='THROW_MPHYS_ERROR'

    IF (itype > 0) THEN

      !------------------------------------------------------------
      ! Produce error message
      !------------------------------------------------------------

      str='Error in CASIM microphysics: '

      IF ( itype <= 3 ) THEN
        str=trim(str)//trim(stdinfo(itype))
      ELSE
        str=trim(str)//trim(stdinfo(3))
      END IF

      str=trim(str)//' Additional information: '//trim(info)
! #if DEF_MODEL==MODEL_KiD
!       print*, 'Runtime is:' , time
! #endif
      PRINT*, routine,':', trim(str)
      PRINT*, (minus_one)**0.5
      STOP

    ELSE IF ( itype < 0 ) THEN

      !------------------------------------------------------------
      ! Produce warning message
      !------------------------------------------------------------

      str='Warning from CASIM microphysics! '
      str=trim(str)//' Message: '//trim(info)
! #if DEF_MODEL==MODEL_KiD
!       print*, 'Runtime is:' , time
! #endif
      PRINT*, routine,':', trim(str)

    END IF

  END SUBROUTINE throw_mphys_error


  SUBROUTINE mphys_message(routine, msg)

    IMPLICIT NONE

    CHARACTER(*), INTENT(IN) :: routine ! Routine providing the message
    CHARACTER(*), INTENT(IN) :: msg ! Message

    CHARACTER(len=*), PARAMETER :: RoutineName='MPHYS_MESSAGE'

    CHARACTER(2000) :: str

    str = '| Message from CASIM microphysics | Routine:' 
    str = trim(str)//trim(routine)//' | Message: '
    str = trim(str)//trim(msg)//' |'

    PRINT *, trim(str)

  END SUBROUTINE mphys_message

END MODULE mphys_die
