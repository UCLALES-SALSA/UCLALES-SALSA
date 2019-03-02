MODULE check_state
  USE mpi_interface, ONLY : appl_abort
  USE classFieldArray, ONLY : FieldArray
  IMPLICIT NONE

  LOGICAL :: breakUndefOutput = .FALSE. ! Namelist variable. TRUE == stop program if output varlists contain undefined variable names

  ! Contains routines to check the workflow of model and output variables

  CONTAINS

    SUBROUTINE checkOutputs(user,field)
      CHARACTER(len=*), INTENT(in) :: user(:)
      TYPE(FieldArray), INTENT(in) :: field

      CHARACTER(len=50), ALLOCATABLE :: fieldvars
      
      INTEGER :: i
      
      ! get the list of variable names defined in the field array instance
      ALLOCATE(fieldvars(field%count))
      DO i = 1,field%count
         fieldvars = field%list(i)%name
      END DO
      
      DO i = 1,SIZE(user)
         IF ( ALL(user(i) /= fieldvars(:)) ) THEN
            IF (breakUndefOutput) THEN
               ! Write and error message and stop program
               WRITE(*,*) 'ERROR: user defined outputlist contains variables which are currently undefined: ',user(i)
               CALL appl_abort(0)
            ELSE
               ! WRite a warning
               WRITE(*,*) 'WARNING: user defined outputlist contains variables which are currently undefined: ',user(i)
            END IF               
         END IF
      END DO
      
    END SUBROUTINE checkOutputs
  
END MODULE check_state
