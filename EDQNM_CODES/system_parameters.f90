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
! MODULE: system_parameters
! LAST MODIFIED: 16 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SYSTEM PARAMETERS FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE system_parameters
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the system parameters, corresponding to the different variant of EDQNM equation we are studying
! is defined here. Note global variables is common amongst these different variants.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE global_variables
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
	IMPLICIT  NONE
    ! _________________________
    ! SYSTEM VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    INTEGER(KIND=4)::ind
    INTEGER(KIND=4)::ind_integral
    INTEGER(KIND=4)::t_step_forcing
    ! _________________________
    DOUBLE PRECISION::viscosity
    DOUBLE PRECISION::forcing
    DOUBLE PRECISION::initial_en
    DOUBLE PRECISION::energy,enstrophy,dissipation_rate
    DOUBLE PRECISION::eddy_constant
    DOUBLE PRECISION::time_visc,time_spec
    ! _________________________
    CHARACTER(LEN=30)::name_sys
    ! _________________________
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec 
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::geom_fac
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::forcer,forcer_template 
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_time,es_time,ds_time
    ! _________________________
    INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE::kqp_triangle_status
    INTEGER(KIND=4),DIMENSION(:,:),ALLOCATABLE::p_ind_min,p_ind_max
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
	SUBROUTINE init_system_parameters
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !       Initialize all the system related parameters.
    ! PREREQUISITE: SUBROUTINE 'init_global_arrays' (in global_variables module)
    !  has to be called before initializing  these variables.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        ! NOTES:
        ! 1. This is forced viscous EDQNM, with forcing given in the first
        !    few shells matching the dissipation rate with slight fluctuations
        !    to keep it random. The time averaged net energy remains constant.
        ! 2. Viscosity levels for resolutions
        !     N45 - Minimum of 0.0005.
        ! 3. Eddy constant is generally not changed.
        ! 4. Two timescales are derived, one from net energy, other from viscosity
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        name_sys    =   'v_08_'

        viscosity   =   8.0D0 * ( 10.0D0 ** ( - 4.0D0 ) )

        initial_en  =   one

        ind_integral     =   FLOOR( DBLE( N ) / 10.0D0 )
        ! Index (position) of integral scale
      
        time_spec        =   one / DSQRT( initial_en * mom( N ) * mom( N ) )
        ! Time scale from energy and largest momentum

        time_visc        =   one / ( viscosity * ( mom( N ) ** two ) )
        ! Time scale from viscosity and largest momentum

        eddy_constant    =   0.54D0
        ! Eddy constant in its expression

        CALL time_to_step_convert(time_visc,t_step_forcing)
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          
	END

    SUBROUTINE init_system_arrays
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to initialize all in-built arrays that need not change during the evolution of system.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION:: mom_p_min,mom_p_max
        DOUBLE PRECISION:: z_f, x_f, y_f
        
        ALLOCATE( p_ind_min( N, N ), p_ind_max( N, N ) ) 
        ALLOCATE( kqp_triangle_status( N, N, N ), geom_fac( N, N, N ) )

        kqp_triangle_status =  0
        geom_fac            =  zero
        ! Reseting to zero for safety
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  L I M I T S   O F    I N T E G R A T I O N 
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        DO k_ind = 1 , N
        DO q_ind = 1 , N
        IF ( q_ind .NE. k_ind ) THEN

            mom_p_min   =  DABS( mom( k_ind ) - mom( q_ind ) )
            mom_p_max   =  DABS( mom( k_ind ) + mom( q_ind ) )

            IF ( mom_p_min .LT. mom( 1 ) ) THEN
                p_ind_min( k_ind, q_ind )   =   1
            ELSE
                p_ind_min( k_ind, q_ind )   =   find_index_ceiling( mom_p_min )
            END IF

            IF ( mom_p_max .GT. mom( N ) ) THEN
                p_ind_max( k_ind, q_ind )   =   N
            ELSE
                p_ind_max( k_ind, q_ind )   =   find_index_floor( mom_p_max )
            END IF

        END IF
        END DO
        END DO
        
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  T R I A N G L E     C H E C K   &   C O S I N E S   O F    I T.
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        DO k_ind = 1 , N
        DO q_ind = 1 , N
        IF ( q_ind .NE. k_ind ) THEN

            DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind ) 
            
                IF ( ( ( mom( q_ind ) + mom( p_ind ) ) .GT. mom( k_ind ) ) &
                .AND. ( DABS(mom( q_ind ) - mom( p_ind )) .LT. mom( k_ind ) ) ) THEN

                    kqp_triangle_status(k_ind, q_ind, p_ind)  =   1
                    
                    CALL find_cosine( k_ind, q_ind, p_ind, z_f )
                    CALL find_cosine( k_ind, p_ind, q_ind, x_f )
                    CALL find_cosine( q_ind, p_ind, k_ind, y_f )
                    ! Finding cosines for all three sides once it is approved to participate in the triad interaction

                    geom_fac( k_ind, q_ind, p_ind )  =  (z_f ** thr ) - x_f * y_f
                    
                END IF      

            END DO

        END IF
        END DO
        END DO

    END
    
    SUBROUTINE find_cosine( i1, i2, i3 , cosine)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this calculate the cosine of angle opposite to side with index 'i3'
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION,INTENT(OUT)::cosine
        INTEGER(KIND=4),INTENT(IN)::i1, i2, i3

        cosine  = ( mom( i1 ) ** two ) + ( mom( i2 ) ** two  ) - ( mom( i3 ) ** two )

        cosine  = cosine / ( two * mom( i1 ) * mom( i2 ) )
        
    END

    INTEGER FUNCTION find_index_floor( mom0 )
    ! ------------
    ! FUNCTION TO: Calculate the largest index whose momentum is smaller than the given momentum 'mom0'
    ! -------------
        DOUBLE PRECISION::mom0
        
        find_index_floor  =   FLOOR ( DLOG( mom0 / mom_base) / log_lambda )

    RETURN
    END

    INTEGER FUNCTION find_index_ceiling( mom0 )
    ! ------------
    ! FUNCTION TO: Calculate the smallest index whose momentum is larger than the given momentum 'mom0'
    ! -------------
        DOUBLE PRECISION::mom0
        
        find_index_ceiling  =   CEILING ( DLOG( mom0 / mom_base) / log_lambda )

    RETURN
    END

END MODULE system_parameters
