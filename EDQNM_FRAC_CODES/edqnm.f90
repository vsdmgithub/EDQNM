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
! PROGRAM : EDQNM
! LAST MODIFIED: 16 November 2020
! _____________________________________
! LIST OF MODULES USED :
!       1. main_run
!       2. solver
!       3. initial_condition
!       4. global_variables
!       5. system_parameters
!       6. constants
!       7. output
!       8. timer_mod
! _____________________________________
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! PROGRAM FOR SOLVING EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
PROGRAM EDQNM
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This program solves the EDQNM equation 
    ! All the work is done in the modules. Calling a few would finish the code. 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE main_run
    USE timer_mod
    
	IMPLICIT NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::visc_ind
    DOUBLE PRECISION,DIMENSION(6)::viscosity_array
    ! Sequence of simulations with different viscosity
    
    CALL start_timer 
    viscosity_array=(/ 12.0D0, 13.0D0, 14.0D0, 15.0D0, 16.0D0, 17.0D0 /)

    DO visc_ind =  1, 6

    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  I  N  I  T  I  A  L  I  Z  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL read_input
    CALL init_global_variables
    CALL init_global_arrays

    viscosity   =   viscosity_array( visc_ind )
    PRINT*,'++++++++++++ VISCOSITY ++++++++++++++++'
    WRITE(*,'(I2,ES10.3)'),visc_ind,viscosity
    PRINT*,'========================================'
    
    CALL init_system_parameters
    CALL init_system_arrays
    ! We get all the variables, and arrays ready to be allocated

    c='y'
    ! Easy way to stop the evolution with only initiation

    IF (c .EQ. 'y') THEN

        CALL pre_analysis
        ! Does time_step check, initial condition and writing details of simulation
        ! Allocating the evolution arrays, if everything is set, 'all_set' will be 1.

        IF (all_set .EQ. 1) THEN

            PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'
            PRINT*,'   S  I  M  U  L  A  T  I  O  N        S  T  A  R  T  S '
            PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'

            CALL time_evolution
            ! Solve the EDQNM equation, in discrete time

            CALL post_analysis
            ! Does the post-analysis, deallocating

            PRINT*,'IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII'      
            PRINT*,'   S  I  M  U  L  A  T  I  O  N        E  N  D  S '
            PRINT*,'TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT'

        END IF

    END IF
   
    END DO
    CALL finish_timer
    
END PROGRAM EDQNM
