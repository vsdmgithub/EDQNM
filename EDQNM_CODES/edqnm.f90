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
    INTEGER(KIND=4)::visc_ind,no_of_visc
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::viscosity_array,dt_array,k_kol_array

    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    ! SEQUENCE OF SIMULATIONS WITH DIFFERENT VISCOSITY
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CALL start_timer
    
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !  I  N  I  T  I  A  L  I  Z  A  T  I  O  N
    !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    CALL read_input

    no_of_visc  =   6
    ALLOCATE( viscosity_array( no_of_visc ), dt_array( no_of_visc ), k_kol_array( no_of_visc) )
    
    viscosity_array = (/ 5.0D0, 8.0D0, 10.0D0, 12.0D0, 16.0D0, 20.0D0 /)
!    viscosity_array=(/ 10.0D0 /)
    ! ARRAY OF VISCOSITIES

    dt_array        = (/ 10.0D0, 10.0D0, 10.0D0, 10.0D0, 10.0D0, 10.0D0 /)
!    dt_array=(/ 2.0D0 /)
    ! ARRAY OF TIMESTEPS FOR EACH VISCOSITY RESPECTIVELY

    k_kol_array     = (/ 511.0D0, 359.0D0, 304.0D0, 265.0D0, 213.0D0, 181.0D0 /)
    ! ARRAY OF K-KOLMO SCALES
    
    DO visc_ind =   1,  no_of_visc
    
    viscosity   =   viscosity_array( visc_ind ) * ( 10.0D0 ** ( - 5.0D0 ) )
    ! Correcting the order of viscosity

    dt          =   dt_array( visc_ind ) * ( 10.0D0 ** ( - 5.0D0 ) )
!    dt          =   dt_array( visc_ind ) * ( 10.0D0 ** ( - 3.0D0 ) )
    ! Correcting the order of time step
    
    PRINT*,'========================================'
    WRITE(*,'(I2,A15,ES10.3)'),visc_ind,'VISCOSITY = ',viscosity
    PRINT*,'========================================'

    CALL init_global_variables
    CALL init_global_arrays   
    ! System is getting ready.

    CALL init_system_parameters
    CALL init_system_arrays
    ! We get all the variables, and arrays ready to be allocated

    CALL pre_analysis
    ! Does time_step check, initial condition and writing details of simulation
    ! Allocating the evolution arrays, if everything is set, 'all_set' will be 1.

    c='y'
    ! Easy way to stop the evolution with only initiation
    
    IF (c .EQ. 'n') THEN

        PRINT*,'========================================'
        WRITE(*,'(A30,ES10.2)'),'RECOMMENDED TIME STEP = ',dt_ref
        WRITE(*,'(A30,ES10.2)'),'GIVEN TIME STEP = ',dt
        PRINT*,'========================================'
           
    END IF

    IF (c .EQ. 'y') THEN

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

    CALL array_deallocation
    ! To deallocate all the used arrays.

    END DO
    
    CALL finish_timer
    
END PROGRAM EDQNM
