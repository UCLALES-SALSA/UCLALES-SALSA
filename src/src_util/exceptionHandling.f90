MODULE exceptionHandling
  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "exceptionHandling"

  ! 
  ! -------------------------------------------------------
  ! Someday this module will perhaps contain beautifully
  ! formulated functions and classes to check consistency
  ! during model runtime, recovery attemps from exceptions
  ! and error message manager (which is the starting point 
  ! for this module).
  !
  ! Juha Tonttila, FMI, 2017
  ! -------------------------------------------------------
  ! 

  CONTAINS

    SUBROUTINE errorMessage(i_global_name,i_name,i_message)
      IMPLICIT NONE
      CHARACTER(len=*), INTENT(in) :: i_global_name, i_name ! Name of the module/class and name of the procedure triggering the message
      CHARACTER(len=*), INTENT(in) :: i_message ! Error message
      
      CHARACTER(len=50), PARAMETER :: name = "errorMessage"

      WRITE(*,*) ""
      WRITE(*,*) ">>>"//TRIM(global_name)//": "//TRIM(name)//":"  
      WRITE(*,*) "---------------------------------------------"
      WRITE(*,*) TRIM(i_global_name)//": "//TRIM(i_name)//": "//TRIM(i_message)
      WRITE(*,*) "---------------------------------------------"
      WRITE(*,*) ""

    END SUBROUTINE errorMessage



END MODULE exceptionHandling
