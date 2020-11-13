MODULE mo_structured_datatypes
   !
   ! Written for UCLALES-SALSA (U-S).
   !
   ! Juha Tonttila, FMI, 2017
   !
   ! Juha: Extended the procedure pointer interface
   !       for improved statistical output management.
   !       FMI, 2020
   !  

   IMPLICIT NONE
   
   ! -------------------------------------------------------------
   ! Holds derived types for simple float arrays which are needed
   ! in order to use the unlimited polymorphic fields in the
   ! FieldArray and ArrayElement classes.


   ! This is the base type for the main datatypes
   TYPE, ABSTRACT :: FloatArray
      CHARACTER(len=50) :: shortName = ''               ! Short name of the variable
      
      CHARACTER(len=50) :: srcName = ''                 ! Holds the short name of the grid-variable used as the
                                                        ! source to calculate e.g. corresponding statistical moments
                                                        ! (mean, variance etc.). For other than statistical ouput
                                                        ! variables this does not need to be defined. This is neither
                                                        ! mandatory for all statistical outputs, e.g more complex
                                                        ! parameters depending on several grid variables. 
   END TYPE FloatArray
   
   ! ------------------------------------------
   
   TYPE, EXTENDS(FloatArray) :: FloatArray0d
      REAL, POINTER :: d => NULL()
      LOGICAL :: Initialized = .FALSE.
      PROCEDURE(sproc0d), POINTER :: onDemand => NULL() ! Procedure pointer to a subroutine
                                                        ! for calculating parameters on-demand.
   END TYPE FloatArray0d
   INTERFACE FloatArray0d
      PROCEDURE :: FloatArray0d_constructor
   END INTERFACE FloatArray0d

   ! ------------------------------------------
   
   TYPE, EXTENDS(FloatArray) :: FloatArray1d
      REAL, POINTER :: d(:) => NULL() ! This is used as the generic access name for the data
      LOGICAL :: Initialized = .FALSE.
      PROCEDURE(sproc1d), POINTER :: onDemand => NULL() ! Procedure pointer to a subroutine
                                                        ! for calculating parameters on-demand.      
   END TYPE FloatArray1d
   INTERFACE FloatArray1d
      PROCEDURE :: FloatArray1d_constructor
   END INTERFACE FloatArray1d

   ABSTRACT INTERFACE 
      SUBROUTINE sproc0d(SELF,output)
        IMPORT FloatArray0d
        CLASS(FloatArray0d), INTENT(in) :: SELF
        REAL, INTENT(out) :: output
      END SUBROUTINE sproc0d
      SUBROUTINE sproc1d(SELF,output)
        IMPORT FloatArray1d
        CLASS(FloatArray1d), INTENT(in) :: SELF
        REAL, INTENT(out) :: output(:)
      END SUBROUTINE sproc1d
   END INTERFACE

   
   ! ------------------------------------------
   
   TYPE FloatArray2d
      REAL, POINTER :: d(:,:) => NULL() ! This is used as the generic access name for the data
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
      LOGICAL :: Initialized = .FALSE.
   END TYPE FloatArray2d
   INTERFACE FloatArray2d
      PROCEDURE :: FloatArray2d_constructor
   END INTERFACE FloatArray2d

   ! -------------------------------------------
   
   TYPE FloatArray3d
      REAL, POINTER :: d(:,:,:) => NULL() ! This is used as the generic access name for the data
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
      LOGICAL :: Initialized = .FALSE.
   END TYPE FloatArray3d
   INTERFACE FloatArray3d
      PROCEDURE :: FloatArray3d_constructor
   END INTERFACE FloatArray3d

   ! -------------------------------------------
   
   TYPE FloatArray4d
      REAL, POINTER :: d(:,:,:,:) => NULL() ! This is used as the generic access name for the data
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
      LOGICAL :: Initialized = .FALSE.
   END TYPE FloatArray4d
   INTERFACE FloatArray4d
      PROCEDURE :: FloatArray4d_constructor
   END INTERFACE FloatArray4d


   
 CONTAINS

   ! ----------------------------------------------------
   FUNCTION FloatArray0d_constructor(name,srcname,trgt)
     IMPLICIT NONE
     TYPE(FloatArray0d), TARGET         :: FloatArray0d_constructor
     CHARACTER(len=*), INTENT(in) :: name
     CHARACTER(len=*), INTENT(in), OPTIONAL :: srcname
     REAL, INTENT(in), TARGET, OPTIONAL :: trgt

     FloatArray0d_constructor%shortName = name
     IF (PRESENT(srcname)) THEN
        FloatArray0d_constructor%srcName = srcname
     ELSE
        ! If srcname is not provided, assume that srcName is the same as shortName,
        ! or that the srcName is not used at all.
        FloatArray0d_constructor%srcName = name
     END IF
     
     IF (PRESENT(trgt)) &
          FloatArray0d_constructor%d => trgt

     FloatArray0d_constructor%Initialized = .TRUE.
          
   END FUNCTION FloatArray0d_constructor

   !---------------------------------------------------------------------

   FUNCTION FloatArray1d_constructor(name,srcname,trgt)
     IMPLICIT NONE
     TYPE(FloatArray1d), TARGET         :: FloatArray1d_constructor
     CHARACTER(len=*), INTENT(in)       :: name
     CHARACTER(len=*), INTENT(in), OPTIONAL :: srcname
     REAL, INTENT(in), TARGET, OPTIONAL :: trgt(:)

     FloatArray1d_constructor%shortName = name
     IF (PRESENT(srcname)) THEN
        FloatArray1d_constructor%srcName = srcname
     ELSE
        ! If srcname is not provided, assume that srcName is the same as shortName,
        ! or that the srcName is not used at all.
        FloatArray1d_constructor%srcName = name       
     END IF

     IF (PRESENT(trgt)) &
          FloatArray1d_constructor%d => trgt     

     FloatArray1d_constructor%Initialized = .TRUE.
     
   END FUNCTION FloatArray1d_constructor

   ! ------------------------------------------------------------------
   
   FUNCTION FloatArray2d_constructor(trgt)
     IMPLICIT NONE
     TYPE(FloatArray2d), TARGET         :: FloatArray2d_constructor
     REAL, INTENT(in), TARGET, OPTIONAL :: trgt(:,:)

     IF (PRESENT(trgt)) &
          FloatArray2d_constructor%d => trgt
     FloatArray2d_constructor%Initialized = .TRUE.
     
   END FUNCTION FloatArray2d_constructor

   ! ----------------------------------------------------------------------------
   
   FUNCTION FloatArray3d_constructor(trgt)
     IMPLICIT NONE
     TYPE(FloatArray3d), TARGET         :: FloatArray3d_constructor
     REAL, INTENT(in), TARGET, OPTIONAL :: trgt(:,:,:)

     IF (PRESENT(trgt)) &
          FloatArray3d_constructor%d => trgt     
     FloatArray3d_constructor%Initialized = .TRUE.
     
   END FUNCTION FloatArray3d_constructor

   ! ------------------------------------------------------------------------
   
   FUNCTION FloatArray4d_constructor(trgt)
     IMPLICIT NONE
     TYPE(FloatArray4d), TARGET         :: FloatArray4d_constructor
     REAL, INTENT(in), TARGET, OPTIONAL :: trgt(:,:,:,:)

     IF (PRESENT(trgt)) &
          FloatArray4d_constructor%d => trgt      
     FloatArray4d_constructor%Initialized = .TRUE.
     
   END FUNCTION FloatArray4d_constructor
   
END MODULE mo_structured_datatypes
