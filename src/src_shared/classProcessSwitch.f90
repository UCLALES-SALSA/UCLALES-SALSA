
MODULE classProcessSwitch
  IMPLICIT NONE

  TYPE ProcessSwitch
     LOGICAL :: switch = .FALSE. ! True or False
     REAL    :: delay = 0.   ! delay time for the process to switch
     
     LOGICAL :: state = .FALSE.  ! This gives the current state of the switch and this should be used
                        ! in the code
     INTEGER :: mode = 1         ! Mode switch, e.g. selecting between different implementation of the same physical process
     

  END TYPE ProcessSwitch

  INTERFACE ProcessSwitch
     PROCEDURE cnstr
  END INTERFACE ProcessSwitch

  CONTAINS

    FUNCTION cnstr()
      IMPLICIT NONE
      TYPE(ProcessSwitch) :: cnstr
      
      cnstr%switch = .FALSE.
      cnstr%delay = 0.
      
      cnstr%state = .FALSE.
      
      cnstr%mode = 1


    END FUNCTION cnstr

END MODULE classProcessSwitch
