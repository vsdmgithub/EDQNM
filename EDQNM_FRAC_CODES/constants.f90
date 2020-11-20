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
! MODULE: constants
! LAST MODIFIED: 16 November 2020
! #########################
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! CONSTANTS FOR EDQNM EQUATION
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
MODULE constants
! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
! ------------
! All the user defined constants (that can't be changed) are declared here. And refered in other modules. 
! -------------
! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

    ! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
	IMPLICIT  NONE

    ! _________________________
    ! CONSTANTS
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DOUBLE COMPLEX,PARAMETER::i=DCMPLX(0.0D0,1.0D0),c0=DCMPLX(0.0D0,0.0D0)
	DOUBLE PRECISION,PARAMETER::two_pi=DATAN(1.0D0)*8.0D0
	DOUBLE PRECISION,PARAMETER::zero=0.0D0
    DOUBLE PRECISION,PARAMETER::qtr=0.25D0,hf=0.5D0,one=1.0D0
    DOUBLE PRECISION,PARAMETER::two=2.0D0,thr=3.0D0,fiv=5.0D0,six=6.0D0
    DOUBLE PRECISION,PARAMETER::tol_double=0.00000000000000001
    ! _________________________
    ! FORMATS 
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    CHARACTER(LEN=*),PARAMETER::f_i2='(I2)',f_i4='(I4)',f_i6='(I6)'
    CHARACTER(LEN=*),PARAMETER::f_i8='(I8)',f_i12='(I12)',f_i16='(I16)'
    CHARACTER(LEN=*),PARAMETER::f_d8p4='(F8.4)',f_d12p6='(F12.6)',f_d16p8='(F16.8)',f_d5p2='(F5.2)'
    CHARACTER(LEN=*),PARAMETER::f_d12p2='(F12.2)',f_d32p17='(F32.17)',f_d12p8='(F12.8)'
    CHARACTER(LEN=*),PARAMETER::f_e5p2='(ES6.2)',f_e10p4='(ES12.4)'
    CHARACTER(LEN=*),PARAMETER::f_c32p17='(F32.17,F32.17)'

    ! _________________________
    ! USER DEFINED DATA-TYPE
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    CONTAINS

    SUBROUTINE get_simulation_name( sim_char )
    ! INFO - START  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    ! ------------
    ! This provides a string  of format
    !     'run_d121120_t104022' for run dated 12/11/2020 timed 10:40:22 
    !  
    ! -------------
    ! INFO - END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

        IMPLICIT NONE
        ! _________________________
        ! LOCAL  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=8)::year_char,month_char,date_char
        CHARACTER(LEN=8)::hour_char,min_char,sec_char
        CHARACTER(LEN=8)::dummy_char
        INTEGER,DIMENSION(8)::values
        ! _________________________
        ! TRANSFER  VARIABLES
        ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        CHARACTER(LEN=*),INTENT(OUT)::sim_char
 
        CALL DATE_AND_TIME(VALUES=values)
        ! Gets the date and time as integer array
        
        values(1)   =   MOD(values(1),2000)
        ! writing only last two digits of year
        
        WRITE(year_char,f_i4), values(1)

        IF ( values(2) .LT. 10 ) THEN
            WRITE(dummy_char,f_i2), values(2)
            month_char   =   '0'//TRIM(ADJUSTL(dummy_char))
        ELSE 
            WRITE(month_char,f_i2), values(2)
        END IF

        IF ( values(3) .LT. 10 ) THEN
            WRITE(dummy_char,f_i2), values(3)
            date_char   =   '0'//TRIM(ADJUSTL(dummy_char))
        ELSE 
            WRITE(date_char,f_i2), values(3)
        END IF

        IF ( values(5) .LT. 10 ) THEN
            WRITE(dummy_char,f_i2), values(5)
            hour_char   =   '0'//TRIM(ADJUSTL(dummy_char))
        ELSE 
            WRITE(hour_char,f_i2), values(5)
        END IF

        IF ( values(6) .LT. 10 ) THEN
            WRITE(dummy_char,f_i2), values(6)
            min_char   =   '0'//TRIM(ADJUSTL(dummy_char))
        ELSE 
            WRITE(min_char,f_i2), values(6)
        END IF

        IF ( values(7) .LT. 10 ) THEN
            WRITE(dummy_char,f_i2), values(7)
            sec_char   =   '0'//TRIM(ADJUSTL(dummy_char))
        ELSE 
            WRITE(sec_char,f_i2), values(7)
        END IF
        
        sim_char    =  'RUN_D'//TRIM(ADJUSTL(date_char))//TRIM(ADJUSTL(month_char))//&
        TRIM(ADJUSTL(year_char))//'_T'//TRIM(ADJUSTL(hour_char))//&
        TRIM(ADJUSTL(min_char))//TRIM(ADJUSTL(sec_char))

    END

END MODULE constants
