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

   ! The FloatArrayXX instances are intialized by providing either
   ! a TARGET array, or a regular array to the constructor. With
   ! the first, the type-bound pointer "d" points to memory
   ! allocated to some external storage (such as to UCLALES-SALSA
   ! master scalar arrays). The second allows to store data within
   ! the FloatArrayXX instances, which is again accessed through the
   ! pointer "data".

   TYPE FloatArray1d
      REAL, POINTER :: d(:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
   END TYPE FloatArray1d
   INTERFACE FloatArray1d
      PROCEDURE :: FloatArray1d_constructor
   END INTERFACE FloatArray1d

   TYPE FloatArray2d
      REAL, POINTER :: d(:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
   END TYPE FloatArray2d
   INTERFACE FloatArray2d
      PROCEDURE :: FloatArray2d_constructor
   END INTERFACE FloatArray2d

   TYPE FloatArray3d
      REAL, POINTER :: d(:,:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
   END TYPE FloatArray3d
   INTERFACE FloatArray3d
      PROCEDURE :: FloatArray3d_constructor
   END INTERFACE FloatArray3d

   TYPE FloatArray4d
      REAL, POINTER :: d(:,:,:,:) => NULL() ! This is used as the generic access name for the data
      REAL, ALLOCATABLE :: alloc(:,:,:,:)      ! This is used to store data to a FloatArrayXX instance. The pointer data will still be used to access it.
   END TYPE FloatArray4d
   INTERFACE FloatArray4d
      PROCEDURE :: FloatArray4d_constructor
   END INTERFACE FloatArray4d

   CONTAINS

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
