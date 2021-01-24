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
! LAST MODIFIED: 05 January 2021
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
    INTEGER(KIND=4)::no_of_triads
    ! _________________________
    DOUBLE PRECISION::viscosity
    DOUBLE PRECISION::viscosity0
    DOUBLE PRECISION::forcing
    DOUBLE PRECISION::mom_kol
    DOUBLE PRECISION::frac_index
    DOUBLE PRECISION::initial_en
    DOUBLE PRECISION::energy,enstrophy
    DOUBLE PRECISION::dissipation_rate,skewness
    DOUBLE PRECISION::eddy_constant
    DOUBLE PRECISION::time_visc,time_spec
    DOUBLE PRECISION::dt_ref
    DOUBLE PRECISION::localness_cutoff_ratio
    ! _________________________
    CHARACTER(LEN=30)::name_sys
    ! _________________________
    DOUBLE PRECISION,DIMENSION(3)::triad_sides
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec 
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::frac_laplacian_k
    DOUBLE PRECISION,DIMENSION(:,:,:),ALLOCATABLE::geom_fac
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::forcer,forcer_template 
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::en_time,es_time
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::ds_time,sk_time
    ! _________________________
    INTEGER(KIND=4),DIMENSION(:,:,:),ALLOCATABLE::kqp_status
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
        ! 2. Viscosity0 levels for resolutions
        !     N45 - Minimum of 0.0005.
        ! 3. Eddy constant is generally not changed.
        ! 4. Two timescales are derived, one from net energy, other from viscosity
        ! 5. This is a fractional viscous model, s varies from (0,1). 
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        name_sys    =   's_zp9_v_'

        frac_index  =   0.9 
        ! Fractional laplacian index, 1 means original laplacian.

        viscosity   =   ( mom_kol ** ( two - two * frac_index ) )
        viscosity   =   viscosity * viscosity0
        ! Modified viscosity, such that at k_kol, the dissipative coefficient is 
        ! same as that of s=1 case.

!       viscosity   =   viscosity0
        ! COMMENT this line to keep Modified viscosity

        initial_en  =   one
        ! Initial energy
        
        ind_integral     =   FLOOR( DBLE( N ) / 10.0D0 )
        ! Index (position) of integral scale
      
        time_spec        =   one / DSQRT( initial_en * mom( N ) * mom( N ) )
        ! Time scale from energy and largest momentum

        time_visc        =   one / ( viscosity * ( mom( N ) ** ( two * frac_index ) ) )
        ! Time scale from viscosity and largest momentum

        eddy_constant    =   0.54D0
        ! Eddy constant in its expression

        localness_cutoff_ratio  =   0.4
        ! Ratio of min to max triad sides, to say it is a nonlocal triad interactions
        ! This is userdefined . Has to be << 1 is must.

        CALL time_to_step_convert(time_visc,t_step_forcing)

        CALL fix_time_step( 0.5D0 * MIN( time_spec, time_visc), dt_ref )
                
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

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  A  R  R  A  Y     A  L  L  O  C  A  T  I  O  N  
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
        ALLOCATE( p_ind_min( N, N ), p_ind_max( N, N ) ) 
        ALLOCATE( kqp_status( N, N, N ), geom_fac( N, N, N ) )
        ALLOCATE( frac_laplacian_k( N ) )
        
        kqp_status          =  0
        geom_fac            =  zero
        ! Reseting to zero for safety
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  L I M I T S   O F    I N T E G R A T I O N 
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        DO k_ind = 1 , N
        DO q_ind = 1 , N

            mom_p_min   =  DABS( mom( k_ind ) - mom( q_ind ) )
            mom_p_max   =  DABS( mom( k_ind ) + mom( q_ind ) )

            IF ( mom_p_min .LT. mom( 1 ) ) THEN
                p_ind_min( k_ind, q_ind )   =   1
            ELSE IF (find_index_floor( mom_p_min ) .LE. 1) THEN
                p_ind_min( k_ind, q_ind )   =    1
            ELSE
                p_ind_min( k_ind, q_ind )   =   find_index_floor( mom_p_min )
            END IF
            ! FINDING THE MINIMUM OF 'p_ind' FOR EVERY PAIR OF 'k_ind,q_ind'
            
            IF ( mom_p_max .GT. mom( N-1 ) ) THEN
                p_ind_max( k_ind, q_ind )   =   N
                
            ELSE IF (find_index_ceiling( mom_p_max ) .GE. N ) THEN
                p_ind_max( k_ind, q_ind )   =   N
                
            ELSE
                p_ind_max( k_ind, q_ind )   =   find_index_ceiling( mom_p_max )+1
            END IF
            ! FINDING THE MAXIMUM OF 'p_ind' FOR EVERY PAIR OF 'k_ind,q_ind'

        END DO
        END DO
        
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  T R I A N G L E     C H E C K   &   C O S I N E S   O F    I T.
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        no_of_triads =  0
        ! Reset for no of triads in the system.
        
        DO k_ind = 1 , N
        DO q_ind = 1 , N

            DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind ) 
            
                IF ( find_triangle_compatibility(k_ind, q_ind, p_ind) .EQ. 1) THEN

                    kqp_status( k_ind, q_ind, p_ind )  =   1
                    ! MEANING, THE GIVEN THREE MOMENTUM 'p,k,q' CAN FORM A TRIANGLE.

                    CALL find_cosine( k_ind, q_ind, p_ind, z_f )
                    CALL find_cosine( k_ind, p_ind, q_ind, y_f )
                    CALL find_cosine( q_ind, p_ind, k_ind, x_f )
                    ! FINDING COSINES FOR ALL THREE SIDES ONCE IT IS APPROVED TO PARTICIPATE IN THE TRIAD INTERACTION

                    geom_fac( k_ind, q_ind, p_ind )  =  (z_f ** thr ) - x_f * y_f
                    ! GEOMETRIC FACTOR IN THE E.D.Q.N.M

                    no_of_triads = no_of_triads + 1
                    ! Counting the triad

                END IF      

            END DO

        END DO
        END DO

        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   T  R  I  A  D      D  E  B  U  G
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        DO k_ind = 1, N
        DO q_ind = 1, N
        DO p_ind = 1, N

            IF (p_ind_max(k_ind,q_ind) .GT. N) THEN
                PRINT*,'ERROR IN UPPER BOUNDARY',p_ind_max(k_ind,q_ind)
            END IF

            IF (p_ind_min(k_ind,q_ind) .LT. 1) THEN
                PRINT*,'ERROR IN LOWER BOUNDARY',p_ind_min(k_ind,q_ind),k_ind,q_ind
            END IF

            IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( k_ind, p_ind, q_ind) ) THEN
                PRINT*,'ERROR IN TRIAD at k=',k_ind,' q=',q_ind,' p=',p_ind
            END IF
            
            IF ( kqp_status( k_ind, q_ind, p_ind ) .NE. kqp_status( q_ind, k_ind, p_ind) ) THEN
                PRINT*,'ERROR IN TRIAD at k=',k_ind,' q=',q_ind,' p=',p_ind
            END IF
        
        END DO
        END DO
        END DO
        
        !XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  F  R  A  C  T  I  O  N  A  L     L  A  P  L  A  C  I  A  N
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        DO k_ind = 1, N

            frac_laplacian_k( k_ind )   =   mom( k_ind ) ** ( two * frac_index )
            ! FRACTIONAL LAPLACIAN - TO BE USED IN WHEREVER VISCOUS TERMS APPEAR
            
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

        ! FINDING COSINE OF ANGLE BETWEEN SIDES 'i1,i2'
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

    INTEGER FUNCTION find_triangle_compatibility( i1, i2, i3 )
    ! ------------
    ! FUNCTION TO: Calculate the compatibility that three momentum can form a triangle, 1 is yes, 0 means no.
    ! -------------
        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER(KIND=4),INTENT(IN)::i1, i2, i3
        
        find_triangle_compatibility =  0

        IF ( mom(i1) + mom(i2) .GT. mom(i3) ) THEN
        IF ( mom(i3) + mom(i2) .GT. mom(i1) ) THEN
        IF ( mom(i1) + mom(i3) .GT. mom(i2) ) THEN

            find_triangle_compatibility =  1

        END IF
        END IF
        END IF

    RETURN
    
    END

END MODULE system_parameters
