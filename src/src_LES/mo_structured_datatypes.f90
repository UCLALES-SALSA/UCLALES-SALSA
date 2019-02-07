MODULE mo_structured_datatypes
   !
   ! Written for UCLALES-SALSA (U-S).
   !
   ! Juha Tonttila, FMI, 2017
   !

   IMPLICIT NONE

   ! -------------------------------------------------------------
   ! Holds derived types for simple float arrays which are needed
   ! in order to use the unlimited polymorphic fields in the
   ! FieldArray and ArrayElement classes.

   ! The FloatArrayXX classes contain three elements. First, the pointer "d", which
   ! provides the access point to the data. "Alloc" is only for internal storage purposes
   ! and should not be directly accessed. The procedure pointer "onDemand" provides the
   ! optional possibility to link a specific subroutine for a specific instance of the
   ! FloatArrayXX class. Typical use example would be: for the instance containing bulk aerosol number
   ! concentration, one would link a subroutine calculating the bulk concentration from
   ! binned concentrations, and this subroutine would be invoked automatically, when writing
   ! the bulk concentration output.

   ! The constructors take 3 arguments:
   !   -  trgt: This is the TARGET data array for the FloatArray instance
   !   -  store: If TRUE, the field "alloc" will be allocated the memory required by "trgt"
   !             and the data is saved in "alloc" and the pointer "d" is associated with alloc.
   !             If False, the data in "trgt" is actually saved in an external TARGET array, and
   !             "d" is associated directly with "trgt".
   !             Thus, in the former case, the FieldArray instance is the storage of the data,
   !             and in the latter case, the FieldArray instance points to data stored somewhere
   !             else. In both casesm the data must be accessed by "d", i.e. myInstance%d(:,:,...)
   !   -  sub: OPTIONAL: the constructor can be given a subroutine name, which "onDemand" will be
   !           associated with. Repeating the example above, for FieldArray instance for bulk aerosol
   !           number concentration, a subroutine calculating the bulk number from binned number concentrations
   !           would be passed, and then the CALL myBulkNumber%onDemand(args) would calculate the bulk value.
   !           this ofcourse forces the external subroutines, especially their interfaces to be pretty generic
   !           and basically all the necessary data has to come from imports from other modules. For examples,
   !           see mo_derived_procedures.
   !
   !           The sub thing may be extended later to cover other than bulk values!

   
   TYPE FloatArray0d
      REAL, POINTER :: d => NULL()
      REAL :: alloc
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
   END TYPE FloatArray0d
   INTERFACE FloatArray0d
      PROCEDURE :: FloatArray0d_constructor
   END INTERFACE
   
   TYPE FloatArray1d
      REAL, POINTER :: d(:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
   END TYPE FloatArray1d
   INTERFACE FloatArray1d
      PROCEDURE :: FloatArray1d_constructor
   END INTERFACE FloatArray1d

   TYPE FloatArray2d
      REAL, POINTER :: d(:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
   END TYPE FloatArray2d
   INTERFACE FloatArray2d
      PROCEDURE :: FloatArray2d_constructor
   END INTERFACE FloatArray2d

   TYPE FloatArray3d
      REAL, POINTER :: d(:,:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
   END TYPE FloatArray3d
   INTERFACE FloatArray3d
      PROCEDURE :: FloatArray3d_constructor
   END INTERFACE FloatArray3d

   TYPE FloatArray4d
      REAL, POINTER :: d(:,:,:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:,:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
      PROCEDURE(), NOPASS, POINTER :: onDemand => NULL()
   END TYPE FloatArray4d
   INTERFACE FloatArray4d
      PROCEDURE :: FloatArray4d_constructor
   END INTERFACE FloatArray4d
   
 CONTAINS

   ! ----------------------------------------------------
   FUNCTION FloatArray0d_constructor(trgt,store)
     IMPLICIT NONE
     TYPE(FloatArray0d), TARGET    :: FloatArray0d_constructor
     REAL, INTENT(in), TARGET      :: trgt
     LOGICAL, INTENT(in), OPTIONAL :: store
     
     LOGICAL :: present_and_true

     IF (PRESENT(store)) THEN
        IF (.NOT. store) THEN
           FloatArray0d_constructor%d => trgt
        ELSE
           present_and_true = .TRUE.
        END IF
     END IF
     
     IF ( .NOT. PRESENT(store) .OR. present_and_true) THEN
        FloatArray0d_constructor%alloc = trgt
        FloatArray0d_constructor%d => FloatArray0d_constructor%alloc
     END IF
          
   END FUNCTION FloatArray0d_constructor

   !---------------------------------------------------------------------

   FUNCTION FloatArray1d_constructor(trgt,store)
     IMPLICIT NONE
     TYPE(FloatArray1d), TARGET         :: FloatArray1d_constructor
     REAL, INTENT(in), TARGET           :: trgt(:)
     LOGICAL, INTENT(in), OPTIONAL      :: store
     
     LOGICAL :: present_and_true
     INTEGER :: dims(1)
     
     present_and_true = .FALSE.
     
     IF (PRESENT(store)) THEN
        IF (.NOT. store) THEN
           FloatArray1d_constructor%d => trgt
        ELSE
           present_and_true = .TRUE.
        END IF
     END IF
     
     IF ( .NOT. PRESENT(store) .OR. present_and_true ) THEN
        dims = SHAPE(trgt)
        ALLOCATE( FloatArray1d_constructor%alloc(dims(1)) )
        FloatArray1d_constructor%alloc = trgt
        ! Associate "d" with the allocated array
        FloatArray1d_constructor%d => FloatArray1d_constructor%alloc
     END IF
     
   END FUNCTION FloatArray1d_constructor

   ! ------------------------------------------------------------------
   
   FUNCTION FloatArray2d_constructor(trgt,store)
     IMPLICIT NONE
     TYPE(FloatArray2d), TARGET         :: FloatArray2d_constructor
     REAL, INTENT(in), TARGET           :: trgt(:,:)
     LOGICAL, INTENT(in), OPTIONAL      :: store
     
     LOGICAL :: present_and_true
     INTEGER :: dims(2)
     
     present_and_true = .FALSE.
     
     IF (PRESENT(store)) THEN
        IF (.NOT. store) THEN
           FloatArray2d_constructor%d => trgt
        ELSE
           present_and_true = .TRUE.
        END IF
     END IF
     
     IF ( .NOT. PRESENT(store) .OR. present_and_true ) THEN
        dims = SHAPE(trgt)
        ALLOCATE( FloatArray2d_constructor%alloc(dims(1),  &
                                                 dims(2))  )
        FloatArray2d_constructor%alloc = trgt
        ! Associate "d" with the allocated array
        FloatArray2d_constructor%d => FloatArray2d_constructor%alloc
     END IF
     
   END FUNCTION FloatArray2d_constructor

   ! ----------------------------------------------------------------------------
   
   FUNCTION FloatArray3d_constructor(trgt,store)
     IMPLICIT NONE
     TYPE(FloatArray3d), TARGET         :: FloatArray3d_constructor
     REAL, INTENT(in), TARGET           :: trgt(:,:,:)
     LOGICAL, INTENT(in), OPTIONAL      :: store
     
     LOGICAL :: present_and_true
     INTEGER :: dims(3)
     
     present_and_true = .FALSE.
     
     IF (PRESENT(store)) THEN
        IF (.NOT. store) THEN
           FloatArray3d_constructor%d => trgt
        ELSE
           present_and_true = .TRUE.
        END IF
     END IF
     
     IF ( .NOT. PRESENT(store) .OR. present_and_true ) THEN
        dims = SHAPE(trgt)
        ALLOCATE( FloatArray3d_constructor%alloc(dims(1),  &
                                                 dims(2),  &
                                                 dims(3))  )
        FloatArray3d_constructor%alloc = trgt
        ! Associate "d" with the allocated array
        FloatArray3d_constructor%d => FloatArray3d_constructor%alloc
     END IF
     
   END FUNCTION FloatArray3d_constructor

   ! ------------------------------------------------------------------------
   
   FUNCTION FloatArray4d_constructor(trgt,store)
     IMPLICIT NONE
     TYPE(FloatArray4d), TARGET         :: FloatArray4d_constructor
     REAL, INTENT(in), TARGET           :: trgt(:,:,:,:)
     LOGICAL, INTENT(in), OPTIONAL      :: store
     
     LOGICAL :: present_and_true
     INTEGER :: dims(4)
     
     present_and_true = .FALSE.
     
     IF (PRESENT(store)) THEN
        IF (.NOT. store) THEN
           FloatArray4d_constructor%d => trgt
        ELSE
           present_and_true = .TRUE.
        END IF
     END IF
     
     IF ( .NOT. PRESENT(store) .OR. present_and_true ) THEN
        dims = SHAPE(trgt)
        ALLOCATE( FloatArray4d_constructor%alloc(dims(1),  &
                                                 dims(2),  &
                                                 dims(3),  &
                                                 dims(4))  )
        FloatArray4d_constructor%alloc = trgt
        ! Associate "d" with the allocated array
        FloatArray4d_constructor%d => FloatArray4d_constructor%alloc
     END IF
     
   END FUNCTION FloatArray4d_constructor
   
END MODULE mo_structured_datatypes
