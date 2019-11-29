MODULE mo_stats_finder
  USE mo_structured_datatypes
  USE mo_field_state
  IMPLICIT NONE

  ! Contains subroutines to find source data for domain mean profile
  ! or domain mean timeseries fields. These are supposed to be used
  ! for the case of averaging an EXISTING 3d field already found in
  ! the FieldArray lists Prog, Diag or Derived. The routines will check
  ! if 

  
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
         ! On-demand variable
         CALL var%onDemand(name,fvar)
      ELSE
         ! Variable stored in memory
         fvar = var%d
      END IF
      
      var => NULL()
      
    END SUBROUTINE stats_get3d

    ! ---------------------

    SUBROUTINE stats_get3d_binned(name,fvar,nstr,nend)
      CHARACTER(len=*), INTENT(in) :: name
      INTEGER, INTENT(in) :: nstr,nend
      REAL, INTENT(inout) :: fvar(:,:,:,:)

      TYPE(FloatArray4d), POINTER :: var => NULL()

      fvar = 0.

      IF ( Prog%Exist(name) ) THEN  ! Check if the variable name exists in full-field FieldArrays
         CALL Prog%getData(1,var,name=name)
      ELSE IF ( Diag%Exist(name) ) THEN
         CALL Diag%getData(1,var,name=name)
      ELSE IF ( Derived%Exist(name) ) THEN
         CALL Derived%getData(1,var,name=name)
      ELSE
         RETURN  ! No source data for output variables found. Result will be a zero field
      END IF

      IF (ASSOCIATED(var%onDemand)) THEN
         ! On-demand variable
         CALL var%onDemand(name,fvar,nstr,nend)
      ELSE
         ! Variable stored in memory
         fvar(:,:,:,:) = var%d(:,:,:,:)
      END IF
      
      var => NULL()      
      
    END SUBROUTINE stats_get3d_binned

    
END MODULE mo_stats_finder
