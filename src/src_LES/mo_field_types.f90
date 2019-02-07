MODULE mo_field_types
  USE classFieldArray
  IMPLICIT NONE  

  SAVE
  
  ! Field arrays for organizing the prognostic and diagnostic variables and their attributes and output status
  ! ------------------------------------------------------------------------------------------------------------
  TYPE(FieldArray) :: Prog
  TYPE(FieldArray) :: Diag
  TYPE(FieldArray) :: Vector
  ! ------------------------------------------------------------------------------------------------------------

  ! Auxiliary FieldArray instances for pre-selected groups
  TYPE(FieldArray) :: SALSA_tracers_4d  ! 4d SALSA tracers (size distributions and compositions)
  TYPE(FieldArray) :: outProg           ! Contains variables from Prog assigned for output
  TYPE(FieldArray) :: outVector         ! Same for Vector variables
  TYPE(FieldArray) :: outDiag           ! Same for Diag

  TYPE(FieldArray) :: Derived
  TYPE(FieldArray) :: outDerived        ! Derived variables - only for output, the data is only calculated using method onDemand
  TYPE(FieldArray) :: outPS             ! Variables for profile statistic outputs. Data is only calculated with onDemand

END MODULE mo_field_types
