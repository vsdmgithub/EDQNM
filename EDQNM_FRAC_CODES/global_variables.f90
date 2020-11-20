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
! MODULE: global_variables
! LAST MODIFIED: 16 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! GLOBAL VARIABLES FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE global_variables
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the global variables and arrays for the simulation space are declared and given values
! here, wheras temporary (IF necessary) are declared within the subroutines.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE constants
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
	IMPLICIT  NONE
    ! _________________________
    ! SPACE VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::N
	INTEGER (KIND=4)::k_ind,q_ind,p_ind    
    ! ---------------------------------------------------------
    DOUBLE PRECISION::mom_base
    DOUBLE PRECISION::log_lambda,lambda
    ! _________________________
    ! TIME (SIMULATION) VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::t_step,t_step_total
	INTEGER (KIND=4)::t_step_save
	INTEGER (KIND=4)::t_step_debug
	INTEGER (KIND=4)::save_no,save_total    
    ! ---------------------------------------------------------
    DOUBLE PRECISION::time_total,time_now,time_save
    DOUBLE PRECISION::dt
    ! _________________________
    ! OUTPUT AND INPUT FILES AND DIR
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER (KIND=4)::state_sim
    ! _________________________
    DOUBLE PRECISION::rd_no
    ! _________________________
    CHARACTER(LEN=60)::input_file
    CHARACTER(LEN=80)::name_sim
    CHARACTER(LEN=60)::N_char
    CHARACTER(LEN=60)::path_dir
    CHARACTER(LEN=60)::name_dir
    ! _________________________
    ! GLOBAL ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::mom, mom_band
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::t_axis
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::laplacian_k
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
    SUBROUTINE read_input
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Read simulation parameters from a file 'input_file'
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
       ! _________________________
       ! LOCAL VARIABLES
       ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        input_file='simulation_parameters.dat'
        ! This file contains all major Input parameters to be fed
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! R  E  A  D  I  N  G       I  N  P  U  T       F  I  L  E
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        OPEN(unit=1001,file=TRIM(ADJUSTL(input_file)))

        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_i6,ADVANCE='yes'),N
        ! No of momentum shells

        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_d12p6,ADVANCE='yes'),dt
        ! Time step size

        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_d8p4,ADVANCE='yes'),time_total
        ! Total time to simulate

        READ(1001,f_d8p4,ADVANCE='yes')
        READ(1001,f_i4,ADVANCE='yes'),save_total
        ! No of saves for the velocity data. 

        CLOSE(1001)
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END

	SUBROUTINE init_global_variables
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the global variables that are used through out the code.
    ! PREREQUISITE: SUBROUTINE 'read_input' has to be called before initializing  these variables.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE

        state_sim=0
        ! This being the first variable to start the simulation. At last, it will be set to '1'

        CALL get_simulation_name( name_sim )
        ! Creating dated and timed name for the simulation
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! S  P  E  C  T  R  U  M        A  N  D         T  I  M  E  
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

!        N       =   61
        ! No of discrete wavenumbers, which will form shells in 3d space with radial width 'dk_i'
       
        lambda      =   two ** ( 0.25D0 )
        ! Ratio of consecutive shells

        log_lambda  =   DLOG( lambda )
        ! Logarithm of Lambda
        
        mom_base    =   two ** ( - thr )
        ! Base momentum

        CALL time_to_step_convert(time_total,t_step_total)
        ! Converts time to steps
        
    	t_step_save   =     t_step_total / save_total
        ! Determines how many time steps after the save has to be made.

        CALL step_to_time_convert(t_step_save,time_save)
        ! Converts steps to time
        
    	t_step_debug =      t_step_total / 10
        ! No of times, the data will be checked for 'NaN' during the simul

        WRITE (N_char, f_i8),N
        ! converting resolution value to CHARACTER

        path_dir    =   '../EDQNM_data/'
        ! path of the main directory relative to this file.

        name_dir    =   'N'//TRIM(ADJUSTL(N_char))//'/'
        ! name of the main directory
        
	END

	SUBROUTINE init_global_arrays
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the global arrays that are declared/allotted here.
    ! PREREQUISITE: SUBROUTINE 'init_global_variables' has to be called before initializing  these arrays.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        
        IMPLICIT  NONE

        ALLOCATE ( mom( N ) , mom_band ( N ) )
        ALLOCATE ( t_axis(0 : t_step_total) )
        ALLOCATE ( laplacian_k( N ) ) 

        ! Defining the band for each discrete k_i as k^{u}_i,k^{l}_i  for above and below
        ! EQN:- $\Delta k_i=k_i \ln(\lambda)$
        ! EQN:- $k_i^{l}=\frac{\Delta k_i}{\lambda -1}
        ! EQN:- $k_i^{u}=\frac{\lambda \Delta k_i}{\lambda -1}
        
        DO k_ind = 1, N

            mom( k_ind )        =   mom_base * ( lambda ** (k_ind - 1) )
            ! Value of momentum for the corresponding index.
            
            mom_band( k_ind )   =   mom( k_ind ) * log_lambda
            ! Gap around the momentum
            
            laplacian_k( k_ind )=   mom( k_ind ) ** two
            ! Laplacian in spectral space
            
        END DO

        DO t_step = 0 , t_step_total

            t_axis ( t_step )   =    dt *  t_step
            ! Time axis declaration
            
        END DO

    END

    SUBROUTINE step_to_time_convert(step,time)
	! CALL this to convert time step into actual time of simulation

    	IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(IN)::step
		DOUBLE PRECISION,INTENT(OUT)::time
        time    =  DBLE(step) * dt

    END

    SUBROUTINE time_to_step_convert(time,step)
	! CALL this to convert time step into actual time of simulation

    	IMPLICIT  NONE
		INTEGER (KIND=4),INTENT(OUT)::step
		DOUBLE PRECISION,INTENT(IN)::time
        step    =   CEILING( time / dt)

    END

END MODULE global_variables
