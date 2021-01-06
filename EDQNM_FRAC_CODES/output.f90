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
! MODULE: output
! LAST MODIFIED: 05 January 2021
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! OUTPUT FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE output
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the outputs from the simulation are done in this module.
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE constants
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
	IMPLICIT  NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4)::M,T
	INTEGER (KIND=4)::j
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    SUBROUTINE simulation_data_import(N0,T0)
    ! This stores data for writing files later

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        INTEGER(KIND=4),INTENT(IN)::N0,T0

        M   =   N0
        T    =  T0

    END

    SUBROUTINE write_spectrum( nam, k, f_k)
    ! CALL this to write real space data

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=*),INTENT(IN)::nam
        DOUBLE PRECISION,DIMENSION( M ),INTENT(IN):: k,f_k

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  P  R  I  N   T          O  U  T  P  U  T 
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        OPEN(unit = 888,  file = nam)
        DO j = 1, M
             WRITE(888,f_d12p6,ADVANCE = 'no'),k(j)
             WRITE(888,f_d32p17,ADVANCE = 'yes'),f_k(j)
        END DO
        CLOSE(888)
        
    END

    SUBROUTINE write_temporal( nam, x, f_x)
    ! CALL this to write temporal data

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=*),INTENT(IN)::nam
        DOUBLE PRECISION,DIMENSION(0:T),INTENT(IN):: x,f_x

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  P  R  I  N   T          O  U  T  P  U  T 
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        OPEN(unit = 777,  file = nam)
        DO j = 0, T
             WRITE(777,f_d12p6,ADVANCE = 'no'),x(j)
             WRITE(777,f_d32p17,ADVANCE = 'yes'),f_x(j)
        END DO
        CLOSE(777)
        
    END

END MODULE output
