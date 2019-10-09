MODULE mo_stats_finder
  USE mo_structured_datatypes
  USE mo_field_state
  IMPLICIT NONE

  CONTAINS

    SUBROUTINE stats_get3d(name,fvar)
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(inout) :: fvar(:,:,:)
      
      TYPE(FloatArray3d), POINTER :: var => NULL()

      fvar = 0.
      
      IF ( Prog%Exist(name) ) THEN  ! Check if the variable name exists in FieldArrays
         CALL Prog%getData(1,var,name=name)
      ELSE IF ( Diag%Exist(name) ) THEN
         CALL Diag%getData(1,var,name=name)
      ELSE IF ( Derived%Exist(name) ) THEN
         CALL Derived%getData(1,var,name=name)
      ELSE         
         RETURN  ! No source data for output variables found. Result will be a zero field
      END IF

      IF (ASSOCIATED(var%onDemand)) THEN
         CALL var%onDemand(name,fvar)
      ELSE
         fvar = var%d
      END IF
      
      var => NULL()
      
    END SUBROUTINE stats_get3d

  
END MODULE mo_stats_finder
