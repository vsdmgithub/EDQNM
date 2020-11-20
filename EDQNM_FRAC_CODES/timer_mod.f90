! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! CODE BY:
! --------   |         |   ---------        /\        |\      |
! |          |         |  |                /  \       | \     |
! |          |         |  |               /    \      |  \    |
! --------   |         |  |   ------|    /------\     |   \   |
!         |  |         |  |         |   /        \    |    \  |
!         |  |         |  |         |  /          \   |     \ |
! ---------   ----------  ----------  /            \  |      \|
! --------------------------------------------------------------------------------------------------------------------------------------------                                                                           
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! -_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_
! #########################
! MODULE: timer_mod
! LAST MODIFIED: 16 November 2020
! #########################
MODULE timer_mod
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This module keeps track of time for the simulation.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

	IMPLICIT  NONE

    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::time_start,time_end,tim_hours
    CHARACTER(LEN=30)::c,date_sim
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
  
    CONTAINS
    
	SUBROUTINE start_timer
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO: START THE TIMER
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT NONE
    
        CALL CPU_TIME(time_start)
        CALL FDATE(date_sim)
        
        WRITE(*,'(A60)'),TRIM(ADJUSTL(' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'))
        WRITE(*,'(2A30)'),TRIM(ADJUSTL(' PROGRAM STARTED ON :')),date_sim
        WRITE(*,'(A60)'),TRIM(ADJUSTL('-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'))
        ! Stores time when the program is initiated
          
	END

	SUBROUTINE finish_timer
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO: FINISH THE TIMER
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT NONE

        CALL FDATE(date_sim)
        CALL CPU_TIME(time_end)    
        tim_hours   =   (time_end-time_start) / 3600.0D0

        WRITE(*,'(A60)'),TRIM(ADJUSTL('-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-'))
        WRITE(*,'(2A30)'),TRIM(ADJUSTL(' PROGRAM ENDED ON :')),date_sim
        WRITE(*,'(A60)'),TRIM(ADJUSTL('TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'))
        WRITE(*,'(A30,F20.4)'),TRIM(ADJUSTL('TOTAL RUN TIME (HOURS) :  ')),tim_hours
        WRITE(*,'(A60)'),TRIM(ADJUSTL(' HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'))
              
	END

END MODULE timer_mod
