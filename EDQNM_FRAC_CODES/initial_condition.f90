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
! MODULE: initial_condition
! LAST MODIFIED: 05 January 2021
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! INITIAL CONDITION FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE initial_condition
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! In this module, initial conditions are provided
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE system_parameters
    
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    IMPLICIT  NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!	INTEGER (KIND=4):: 
    ! ---------------------------------------------------------
!    DOUBLE PRECISION::
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS

    SUBROUTINE make_initial_condition(en0,spec0)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Initialize initial condition 
    ! INPUT : Energy 
    ! OUTPUT : Complex array "spectrum" 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION::A0,s_exp

        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION,INTENT(IN)::en0
        DOUBLE PRECISION,DIMENSION(N),INTENT(OUT)::spec0
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        A0      =   one
        ! Normalization constant

        s_exp   =   two
        ! Integral scale spectrum exponent
        
        DO k_ind = 1, N

            spec0( k_ind )   =   A0 * ( mom( k_ind ) ** s_exp )
            spec0( k_ind )   =   spec0( k_ind ) * DEXP( - hf * ( mom( k_ind ) / mom( ind_integral ) ) ** two )
            
        END DO

        A0  =   en0 / SUM( spec0 * mom_band )
        ! Adjustng the normalization constant.

        DO k_ind = 1, N

            spec0( k_ind )   =   A0 * ( mom( k_ind ) ** s_exp )
            spec0( k_ind )   =   spec0( k_ind ) * DEXP( - hf * ( mom( k_ind ) / mom( ind_integral ) ) ** two )
            
        END DO

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        CALL make_forcing_template

	END

    SUBROUTINE read_initial_condition(en0,spec0)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Initialize initial condition 
    ! INPUT : Energy 
    ! OUTPUT : Complex array "spectrum" 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION::A0,dum
        CHARACTER(LEN=100)::spec0_address
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION,INTENT(IN)::en0
        DOUBLE PRECISION,DIMENSION(N),INTENT(OUT)::spec0
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  I   N   I   T   I   A   L              C    O    N    D    I    T    I    O     N
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        spec0_address=TRIM(ADJUSTL(path_dir)) // 'N' // TRIM(ADJUSTL(N_char)) 
        ! where initial condition is stored

        OPEN( UNIT = 2001 ,FILE = TRIM(ADJUSTL(spec0_address))//'.dat' )        

        DO k_ind = 1, N

            READ( 2001, '(F12.6)',ADVANCE='NO'),dum
            READ( 2001, '(F32.17)',ADVANCE='YES'),spec0( k_ind )
            
        END DO

        CLOSE( 2001 )
        
        A0  =   en0 / SUM( spec0 * mom_band )
        ! Adjustng the normalization constant.

        DO k_ind = 1, N

            spec0( k_ind )   =   A0 * spec0( k_ind )
            
        END DO

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        CALL make_forcing_template

	END

    SUBROUTINE make_forcing_template
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    ! Initialize a template of forcing similar to the initial condition.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 

        IMPLICIT  NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION::A0,s_exp

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !  F  O  R  C  I  N  G       T  E  M  P  L  A  T  E
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        A0      =   zero
        ! Normalization constant

        s_exp   =   two
        ! Integral scale spectrum exponent
        
        DO k_ind = 1, N

            A0 = A0 + ( mom( k_ind ) ** s_exp ) * DEXP( - hf * ( mom( k_ind ) / &
                        mom( ind_integral ) ) ** two ) * mom_band( k_ind )
            
        END DO

        DO k_ind = 1, N

            forcer_template( k_ind )   =   ( mom( k_ind ) ** s_exp )
            forcer_template( k_ind )   =   forcer_template( k_ind ) * DEXP(&
                                         - hf * ( mom( k_ind ) / mom( ind_integral ) ) ** two )
            
        END DO

        forcer_template     =    forcer_template    /   A0

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

	END



END MODULE initial_condition
