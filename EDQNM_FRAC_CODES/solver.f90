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
! MODULE: solver
! LAST MODIFIED: 16 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! SOLVER FOR  EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE solver
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! Takes the spectral velocity and updates it by a step, using the subroutines
! 1. rk4_algorithm
! 2. time_derivative
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
	USE initial_condition
	
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    IMPLICIT NONE
    ! _________________________
    ! SOLVER VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION::viscous_freq,eddy_freq 
    DOUBLE PRECISION::eddy_k,eddy_q,eddy_p
    DOUBLE PRECISION::integrand
    DOUBLE PRECISION::damping
   ! _________________________
    ! SOLVER ARRAYS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::d_spec1, d_spec2, d_spec3, d_spec4
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::spec_temp,eddy_array
    DOUBLE PRECISION,DIMENSION(:),ALLOCATABLE::transfer_spec,flux
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH

    CONTAINS
    
	SUBROUTINE rk4_algorithm
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to USE RK4 algorithm to move one step forward in time for the matrix 'spec(k,t)-> spec(k,t+1)'
    ! Alg: - Runga kutta 4th order
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! First store the spectral velocity into a temporary matrix, as steps of RK4 algorithm will manipulate 'spec(k)''

        spec_temp   =   spec
        CALL time_derivative(d_spec1) ! This call provides the time derivative for the spectral 

    	spec   =   spec_temp  +   hf * d_spec1
        CALL time_derivative(d_spec2)

        spec   =   spec_temp  +   hf * d_spec2
        CALL time_derivative(d_spec3)
        
		spec   =   spec_temp  +   d_spec3
        CALL time_derivative(d_spec4)

        ! Final increment for 'v(k)'
		spec   =   spec_temp  +   ( d_spec1 + two * d_spec2 + two * d_spec3 + d_spec4 ) / six

    END
    
    SUBROUTINE time_derivative(d_spec)
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this to get the time derivative matrix for matrix 'spec(k_ind)'
    ! This is the EDQNM EQUATION implemented for numerical computation
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        DOUBLE PRECISION ,DIMENSION( N ),INTENT(OUT)::d_spec

        CALL transfer_term
        
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   E   D   Q   N   M          E   Q   N.
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        d_spec  =   transfer_spec 
        ! The transfer term 
        
        d_spec  =   d_spec - two * viscosity * frac_laplacian_k * spec
        ! The dissipative term

        d_spec  =   d_spec + forcer
        ! The forcing term

        d_spec  =   dt * d_spec
        ! Increment in spectrum

	END

    SUBROUTINE transfer_term
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL this 
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   T   R   A   N   S   F   E   R       T   E   R   M
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        transfer_spec   =   zero
        ! Reseting the transfer term
        
        DO k_ind = 1 , N
        DO q_ind = 1 , N
        IF ( q_ind .NE. k_ind ) THEN

            DO p_ind = p_ind_min( k_ind, q_ind ), p_ind_max( k_ind, q_ind )
            ! The third momentum runs through p_min to p_max, determined from the 'k,q'

                IF ( kqp_triangle_status( k_ind, q_ind, p_ind )  .EQ. 1 ) THEN
                ! If three momentum can form a triangle

                    CALL transfer_term_integrand
                    
                    transfer_spec( k_ind ) = transfer_spec( k_ind ) + &
                    integrand * mom_band( q_ind ) * mom_band( p_ind )

                END IF
                
            END DO

        END IF
        END DO
        END DO

    END
    
    SUBROUTINE transfer_term_integrand
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to .
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
      
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   I  N  T  E  G  R  A  N  D
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        CALL damping_factor

        integrand  =  (mom( k_ind )**(two)) * spec( q_ind ) - ( mom( q_ind )**(two) ) * spec( k_ind ) 
        
        integrand  =  integrand * spec( p_ind ) / mom( p_ind )
        
        integrand  =  integrand * geom_fac( k_ind, q_ind, p_ind )
         
        integrand  =  integrand * damping  
          
    END

    SUBROUTINE damping_factor
        ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
        ! ------------
        ! CALL this to .
        ! -------------
        ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   E  D  D  Y            F  R  E  Q  U  E  N  C  Y
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

        eddy_array  =   laplacian_k * spec * mom_band
        
        eddy_k   =   SUM( eddy_array( : k_ind) ) ** hf
        eddy_q   =   SUM( eddy_array( : q_ind) ) ** hf
        eddy_p   =   SUM( eddy_array( : p_ind) ) ** hf
         
        eddy_freq   = eddy_constant * ( eddy_k + eddy_q + eddy_p )      

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        !   D   A   M   P   I   N   G       F  A  C  T  O  R 
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        
        viscous_freq  =  viscosity * ( laplacian_k( k_ind ) + &
                        laplacian_k( q_ind ) + laplacian_k( p_ind ) )

        damping        =  one - DEXP( -( eddy_freq + viscous_freq ) * time_now )

        damping        =  damping / ( eddy_freq + viscous_freq )
        
    END
    
END MODULE solver
