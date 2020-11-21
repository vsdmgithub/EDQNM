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
! MODULE: main
! LAST MODIFIED: 16 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! MAIN MODULE FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE main_run
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! This is the main module. All the other modules are sub-modules to this.
! Then nesting of all sub-modules is as follows
! MAIN MODULE
!   |
!   ∟ ---> OUTPUT MODULE ---> CONSTANTS MODULE
!   |
!   ∟ ---> SOLVER MODULE
!               |
!               ∟--> INITIAL CONDITION MODULE
!                        |
!                        ∟--> SYSTEM PARAMETERS MODULE
!                                |
!                                ∟-->  GLOBAL VARIABLES MODULE
!                                     |  
!                                     ∟--> CONSTANTS MODULE
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! [[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[[
    !  SUB-MODULES
    !  ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]
    USE solver
    USE output
    
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
	IMPLICIT  NONE
    ! _________________________
    !  VARIABLES
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	INTEGER (KIND=4):: all_set   
    ! ---------------------------------------------------------
!    DOUBLE PRECISION:: energy, enstrophy 
    ! ---------------------------------------------------------
    CHARACTER(LEN=100)::file_location
    CHARACTER(LEN=100)::file_address
    CHARACTER(LEN=20)::file_time
    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
    
    CONTAINS
    SUBROUTINE pre_analysis
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Call this to check the validity of parameters for the simulation, like time step, etc.,
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !       T    I    M     E              S    T    E    P              C   H    E   C   K
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF ( ( dt .LT.  time_spec ) .AND. ( dt .LT. time_visc) ) THEN

            all_set =  1
            
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  A  R  R  A  Y        D  E  A  L  L  O  C  A  T  I  O  N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            ALLOCATE( spec( N ) )
            ALLOCATE( forcer( N ), forcer_template( N ) )
            ALLOCATE( en_time (0 : t_step_total) )
            ALLOCATE( es_time (0 : t_step_total) )
            ALLOCATE( ds_time (0 : t_step_total) )
            ALLOCATE( d_spec1( N ), d_spec2( N ), d_spec3( N ), d_spec4( N ))
            ALLOCATE( spec_temp( N ), transfer_spec( N ), eddy_array( N ) )
            ALLOCATE( flux( N ) )

            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  I  N  I  T  I  A  L        C  O  N  D  I  T  I  O  N
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
            CALL make_initial_condition(initial_en, spec)
            ! Calls the subroutine to make a initial condition with normalized energy as 'initial_en'

!            CALL read_initial_condition(initial_en, spec)
            ! Calls the subroutine to read a initial condition from file with normalized energy as 'initial_en'

            CALL write_details
            ! Writes the parameters corresponding to the simulation for reference

            CALL simulation_data_import(N,t_step_total)
            ! Copies the basic data of simulation to a subroutine in output module for easy saving

            file_location=TRIM(ADJUSTL(path_dir))//TRIM(ADJUSTL(name_dir))&
            //TRIM(ADJUSTL(name_sys))//TRIM(ADJUSTL(name_sim))//'/'
            
        ELSE
        
            all_set =   0

            WRITE(*,'(A50)'),'ERROR: TIME STEP TOO LARGE'
            WRITE(*,'(A50)'),'----------------------------------------------------------------------'
            WRITE(*,'(A50,F10.6)'),' RESET THE TIME STEP (AT MAX) AS :',MIN(time_spec,time_visc)

         END IF

    END
    
	SUBROUTINE write_details
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! CALL THIS SUBROUTINE TO:
    !   Write the details of the simulation
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< 
        IMPLICIT  NONE

        CALL SYSTEM('mkdir '//TRIM(ADJUSTL(path_dir))//TRIM(ADJUSTL(name_dir)))
        
        CALL SYSTEM('mkdir '//TRIM(ADJUSTL(path_dir))//TRIM(ADJUSTL(name_dir))&
        //TRIM(ADJUSTL(name_sys))//TRIM(ADJUSTL(name_sim))//'/')
        ! Command to create the main directory and sub directory (name_sim) in the desired path
 
        file_address=TRIM(ADJUSTL(path_dir))    //  TRIM(ADJUSTL(name_dir)) //  &
        'details_'//TRIM(ADJUSTL(name_sys))//TRIM(ADJUSTL(name_sim))

        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
        OPEN(UNIT=233,FILE=TRIM(ADJUSTL(file_address))//'.dat')

        WRITE(233,"(A40)"),TRIM(ADJUSTL('--------------------------------------------------------------------'))
        WRITE(233,"(A40)"),TRIM(ADJUSTL('------  EDQNM (FRAC) EQUATION----------------------'))
        WRITE(233,"(A40)"),TRIM(ADJUSTL(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -'))
        WRITE(233,"(A40)"),TRIM(ADJUSTL('-----------PARAMETERS OF SIMULATION------------'))
        WRITE(233,"(A40)"),TRIM(ADJUSTL('--------------------------------------------------------------------'))
        WRITE(233,"(A2,A20,A2,I5)"),'1.','No of modes    ','= ',N
        WRITE(233,"(A2,A20,A2,ES12.5)")'2.','Time step   ','= ',dt
        WRITE(233,"(A2,A20,A2,I8)")'3.',' Total time steps   ','= ',t_step_total
        WRITE(233,"(A2,A20,A2,F5.2)")'4.','Total time ','= ',time_total
        WRITE(233,"(A2,A20,A2,F8.6)")'5.',' Viscosity  ','= ',viscosity
        WRITE(233,"(A2,A20,A2,F6.3)")'6.',' Fractional index  ','= ',frac_index
        WRITE(233,"(A2,A20,A2,I5)")'7.',' No of saves   ','= ',save_total
        WRITE(233,"(A2,A20,A2,F6.3)")'8.',' Initial energy ','= ',initial_en

        CLOSE(233)
        ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
          
	END
    
    SUBROUTINE time_evolution
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! Loop of time steps, where at each step the spectral velocities
    ! are updated through any of the algoritm. Meanwhile, inter_analysis and
    ! outputs are printed respectively.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE

        ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
        !      B   U   R   G   E   R   S      E   V   O   L   U   T   I   O   N
        ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !             S        T         A         R        T
        ! 8888888888888888888888888888888888888888888888888888888888888888

       WRITE(*,'(A40)'),'-----------------------------------------------------------'
       WRITE(*,'(A43)'),' |   TIME    |    ENERGY     | DISSIP RATE    |'
       WRITE(*,'(A42)'),'-----------------------------------------------------------'
       
        DO t_step = 0,  t_step_total

            CALL inter_analysis
                    
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            !  P  S  E  U  D  O  -  S  P  E  C  T  R  A  L     A  L   G  O  R  I  T  H  M
            !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            CALL rk4_algorithm
            ! Updates velocity field as per EDQNM equation for next time step

            IF (all_set .NE. 1) THEN
                EXIT
                ! Meaning 'NaN' is encountered during the Debug
            END IF
            
        END DO              
        PRINT*,'-----------------------------------------------------------'
        ! |||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
        !                    E     N     D 
        ! 8888888888888888888888888888888888888888888888888888888888888888
        
        state_sim=1
        ! Stating that the simulation has ended.

    END

    SUBROUTINE inter_analysis
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This does all the inter_analysis, making calls to write output during the evolution, debug and statistics part.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  S  A  V  I  N  G    D  A  T  A
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

        IF (MOD(t_step,t_step_save) .EQ. 0) THEN

            CALL step_to_time_convert(t_step, time_now)
            ! Converts the 't_step' to actual time 'time_now'

            WRITE (file_time,f_d8p4),time_now
            ! Writes 'time_now' as a CHARACTER

            file_address  =   TRIM(ADJUSTL(file_location))  //  'spectrum_t_'   //  TRIM(ADJUSTL(file_time))//'.dat'
            CALL write_spectrum(file_address,mom,spec)
            !  ENERGY FILE

            file_address  =   TRIM(ADJUSTL(file_location))  //  'dissipation_t_'   //  TRIM(ADJUSTL(file_time))//'.dat'
            CALL write_spectrum(file_address,mom,two * viscosity * frac_laplacian_k * spec)
            !  DISSIPATION FILE

            CALL transfer_term
            file_address  =   TRIM(ADJUSTL(file_location))  //  'transfer_t_'   //  TRIM(ADJUSTL(file_time))//'.dat'
            CALL write_spectrum(file_address,mom,transfer_spec)
            !  ENERGY TRANSFER FILE
            
            DO ind = 1, N
                flux( ind ) = - SUM( transfer_spec( : ind) * mom_band( : ind) )
            END DO
            file_address  =   TRIM(ADJUSTL(file_location))  //  'flux_t_'   //  TRIM(ADJUSTL(file_time))//'.dat'
            CALL write_spectrum(file_address,mom,flux)
            !  FLUX  FILE
       
        END IF

        energy           =   SUM( spec * mom_band )
        en_time(t_step)  =   energy

        enstrophy        =   SUM( laplacian_k * spec * mom_band )
        es_time(t_step)  =   enstrophy
        
        dissipation_rate =   two * viscosity * SUM( frac_laplacian_k * spec * mom_band )
        ds_time(t_step)  =   dissipation_rate 

        !  ENERGY,ENSTROPHY, AND DISSIPATION VS TIME 

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  F  O  R  C  I  N  G      I  N  I  T  I  A  T  I  N  G
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (MOD(t_step,t_step_forcing) .EQ. 0) THEN

            CALL RANDOM_NUMBER(rd_no)
            rd_no   =   ( rd_no - hf ) ** thr  +  one
            rd_no   =   DSQRT( rd_no )
            forcer  =   rd_no * dissipation_rate * forcer_template
            ! FORCING

        END IF
        
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  D  E  B  U  G         F  O  R          N  a   N
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        IF (MOD(t_step,t_step_debug) .EQ. 0) then
            
            CALL step_to_time_convert(t_step, time_now)
            ! Converts the 't_step' to actual time 'time_now'

            DO ind = 1, N
                IF ( spec (ind) .NE. spec (ind) ) THEN

                    all_set =  0

                    PRINT*,"NaN encountered before t=",time_now
                    EXIT
                    ! IF any NaN is encountered, the loop is exited without any further continuation.

                END IF
           END DO

           WRITE(*,'(A4,F8.4,A4,F12.8,A4,F12.8,A4)'),' | ',time_now,' | '&
           ,energy,' | ',dissipation_rate,' | ' 

        END IF
       
     END

    SUBROUTINE post_analysis
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This does all the post analysis, making calls to write output after the evolution, debug and statistics part.
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        IMPLICIT NONE

        file_address  =  TRIM(ADJUSTL(file_location))   //     'energy_vs_time.dat'
        CALL write_temporal(file_address,t_axis,en_time)
        ! ENERGY VS TIME FILE

        file_address  =  TRIM(ADJUSTL(file_location))   //     'enstrophy_vs_time.dat'
        CALL write_temporal(file_address,t_axis,es_time)
        ! ENSTROPHY VS TIME FILE

        file_address  =  TRIM(ADJUSTL(file_location))   //     'viscous_loss_vs_time.dat'
        CALL write_temporal(file_address,t_axis,ds_time)
        ! DISSIPATION VS TIME FILE

        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !  A  R  R  A  Y        D  E  A  L  L  O  C  A  T  I  O  N
        !  ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        DEALLOCATE(spec)
        DEALLOCATE(frac_laplacian_k)
        DEALLOCATE(geom_fac)
        DEALLOCATE(forcer,forcer_template)
        DEALLOCATE(en_time,es_time,ds_time)
        DEALLOCATE(d_spec1, d_spec2, d_spec3, d_spec4)
        DEALLOCATE(spec_temp, transfer_spec, eddy_array )
        DEALLOCATE(flux )

        DEALLOCATE(p_ind_max,p_ind_min)
        DEALLOCATE(kqp_triangle_status)
        
        DEALLOCATE(mom,mom_band)
        DEALLOCATE(t_axis)
        DEALLOCATE(laplacian_k)
       
     END

END MODULE main_run
