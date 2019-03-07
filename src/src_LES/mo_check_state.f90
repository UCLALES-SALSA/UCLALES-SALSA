MODULE mo_check_state
  USE mpi_interface, ONLY : appl_abort
  USE classFieldArray, ONLY : FieldArray
  IMPLICIT NONE

  LOGICAL :: breakUndefOutput = .FALSE. ! Namelist variable. TRUE == stop program if output varlists contain undefined variable names

  ! Contains routines to check the workflow of model and output variables

  CONTAINS

    SUBROUTINE checkOutputs(N_given,exist,user,field)
      INTEGER, INTENT(in) :: N_given
      LOGICAL, INTENT(inout) :: exist(:)
      CHARACTER(len=*), INTENT(in) :: user(:)
      TYPE(FieldArray), INTENT(in) :: field

      CHARACTER(len=50), ALLOCATABLE :: fieldvars(:)

      INTEGER :: i
      
      IF (SIZE(exist) /= SIZE(user)) THEN
         WRITE(*,*) 'mo_check_state ERROR: the output logical mask should', &
                    'be the same size as the input variable list'
         CALL appl_abort(11)
      END IF
         
      ! get the list of variable names defined in the field array instance
      ALLOCATE(fieldvars(field%count))
      DO i = 1,field%count
         fieldvars(i) = field%list(i)%name
      END DO

      DO i = 1,N_given
            exist(i) = exist(i) .OR. ANY(user(i) == fieldvars(:))            
      END DO
        
      DEALLOCATE(fieldvars)
      
    END SUBROUTINE checkOutputs
  
END MODULE mo_check_state
