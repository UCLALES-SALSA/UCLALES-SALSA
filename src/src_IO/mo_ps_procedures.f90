MODULE mo_ps_procedures
  USE util, ONLY : get_avg3
  USE mo_structured_datatypes, ONLY : FloatArray3d
  USE grid, ONLY : nzp,nxp,nyp
  USE mo_field_types, ONLY : Prog, Diag, Derived
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: globalAvgProfile
  
  CONTAINS

    SUBROUTINE globalAvgProfile(name,output)
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output(nzp)
      
      TYPE(FloatArray3d), POINTER :: var => NULL()
      REAL :: fvar(nzp,nxp,nyp)

      ! THIS EXPECTS THAT THE CORRESPONDING (FULL FIELD VS AVERAGE) VARIABLES
      ! HAVE IDENTICAL NAMES IN THE FIELDARRAY LISTS. THE POINTER VARIABLE NAMES
      ! CAN BE WHATEVER, THOUGH.

      output = 0.
      
      IF ( Prog%Exist(name) ) THEN  ! Check if the variable name exists in FieldArrays
         CALL Prog%getData(1,var,name=name)
         fvar = var%d
      ELSE IF ( Diag%Exist(name) ) THEN
         CALL Diag%getData(1,var,name=name)
         fvar = var%d
      ELSE IF ( Derived%Exist(name) ) THEN
         CALL Derived%getData(1,var,name=name)
         CALL var%onDemand(name,fvar)
      ELSE         
         RETURN  ! No source data for output variables found ????
      END IF

      ! cond?
      CALL get_avg3(nzp,nxp,nyp,fvar,output)      
      
    END SUBROUTINE globalAvgProfile
    

END MODULE mo_ps_procedures
