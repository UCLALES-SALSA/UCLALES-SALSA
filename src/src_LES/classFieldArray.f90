MODULE classFieldArray
   !
   ! Written for UCLALES-SALSA (U-S).
   !
   ! Juha Tonttila, FMI, 2017
   !

   USE mo_structured_datatypes, ONLY : FloatArray0d, FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
   IMPLICIT NONE

   TYPE ArrayElement
      CHARACTER(len=10)  :: name          ! Short name: for output variable name and for fetching data using getData
      CHARACTER(len=100) :: long_name     ! Long name, mainly for output attributes
      CHARACTER(len=10)  :: unit          ! Unit of the variable, e.g. "kg/kg"; used mainly for output attributes
      CHARACTER(len=10)  :: dimension     ! String that gives the dimension environment for output (see ncio.f90)
      CHARACTER(len=10)  :: group         ! A group tag that can be used to fetch a list of certain type variables
      LOGICAL            :: outputstatus  ! TRUE: write this variable to an output file. FALSE: don't.

      ! Below are pointers to scalar variable arrays. In U-S, scalars have two definitions:
      ! the previous value (p) and the tendency (t). These are given here as unlimited
      ! polymorphic pointers. This makes it possible to associate them with any Fortran intrinsic
      ! or derived data type, without duplicating the code. However, for arrays of different ranks,
      ! derived datatypes must be used, which are available for 1-4d arrays (easy to add more...)
      ! in mo_sturctured_datatypes.f90

      ! For diagnostic variables, where tendency is not needed, just use "p". Providing
      ! the tendency arrays e.g. for the constructor routine is optional.

      ! Note, that the polymorphic variables "p" and "t" can be also used to hold more elaborate datatypes,
      ! than FloatArrays, such as aerosol size distributions with binned concentration and composition data
      ! (under development).

      CLASS(*), POINTER :: p => NULL()
      CLASS(*), POINTER :: t => NULL()

      CONTAINS
       
         PROCEDURE :: get_p
         PROCEDURE :: get_t

   END TYPE ArrayElement

   !-----------------------
   INTERFACE ArrayElement
      PROCEDURE :: ArrayElement_constructor
   END INTERFACE ArrayElement

   !----------------------------------------------------------------------------
   !----------------------------------------------------------------------------

   TYPE FieldArray
      TYPE(ArrayElement), ALLOCATABLE  :: list(:)     ! Each element holds data and attributes to one variable given by the class ArrayElement
      INTEGER                          :: count       = 0       ! Number of entries, initialized as 0
      LOGICAL                          :: Initialized = .FALSE. ! Initialized as .FALSE. will be TRUE after the first entry is made

      CONTAINS

         PROCEDURE :: Extend_FieldArray
         GENERIC   :: Extend => Extend_FieldArray

         PROCEDURE :: ExN_FieldArray
         GENERIC   :: ExN => ExN_FieldArray

         PROCEDURE :: NewField

         PROCEDURE :: getField

         PROCEDURE :: getByGroup
         PROCEDURE :: getByOutputstatus
         
         PROCEDURE :: getVarInst_ind, getVarInst_name
         GENERIC   :: getVarInst => getVarInst_ind, getVarInst_name

         PROCEDURE :: getData_0d, getData_1d, getData_2d, getData_3d, getData_4d
         GENERIC   :: getData => getData_0d, getData_1d, getData_2d, getData_3d, getData_4d

         PROCEDURE :: getFieldIndex
       
         PROCEDURE :: destroy_FieldArray
         GENERIC   :: destroy => destroy_FieldArray


   END TYPE FieldArray
   !------------------------
   INTERFACE FieldArray
      PROCEDURE :: FieldArray_constructor
   END INTERFACE FieldArray

   ! ==========================
   CONTAINS
   ! ==========================

   !
   ! ---------------------------
   ! CONSTRUCTORS
   !
   FUNCTION ArrayElement_constructor(in_name, in_long_name, in_unit, in_dimension, in_outputstatus, in_p_data, in_t_data, in_group)
      !
      ! --------------------------------
      ! Instantiate a new ArrayElement
      !
      IMPLICIT NONE
      TYPE(ArrayElement)                        :: ArrayElement_constructor
      CHARACTER(len=*), INTENT(in)              :: in_name        ! Variable name
      CHARACTER(len=*), INTENT(in)              :: in_long_name     ! Long name, mainly for output attributes
      CHARACTER(len=*), INTENT(in)              :: in_unit          ! Unit of the variable, e.g. "kg/kg"; used mainly for output attributes
      CHARACTER(len=*), INTENT(in)              :: in_dimension     ! String that gives the dimension environment for output (see ncio.f90)
      LOGICAL, INTENT(in)                       :: in_outputstatus
      CLASS(*), INTENT(in), POINTER             :: in_p_data      ! Polymorphic pointer to data (values)
      CLASS(*), INTENT(in), POINTER, OPTIONAL   :: in_t_data      ! - '' - (tendencies)
      CHARACTER(len=*), INTENT(in), OPTIONAL    :: in_group

      ArrayElement_constructor%name = in_name
      ArrayElement_constructor%long_name = in_long_name
      ArrayElement_constructor%unit = in_unit
      ArrayElement_constructor%dimension = in_dimension
      ArrayElement_constructor%outputstatus = in_outputstatus
      ArrayElement_constructor%p => in_p_data

      ArrayElement_constructor%group = "default"
      IF (PRESENT(in_group)) ArrayElement_constructor%group = in_group
      IF (PRESENT(in_t_data)) ArrayElement_constructor%t => in_t_data

   END FUNCTION ArrayElement_constructor

   ! ------------------------------------------

   FUNCTION FieldArray_constructor()
      !
      ! ---------------------------------------------------------------
      ! Initialize a FieldArray instance: just puts count to zero and
      ! sets up some switches with initial values
      !
      IMPLICIT NONE
      TYPE(FieldArray) :: FieldArray_constructor

      FieldArray_constructor%count = 0
      FieldArray_constructor%Initialized = .FALSE.
    
   END FUNCTION FieldArray_constructor
  
   !
   ! --------------------------------------------
   ! PROCEDURES BOUND TO FieldArray
   !
   SUBROUTINE newField(SELF, in_name, in_long_name, in_unit, in_dimension, in_outputstatus, in_p_data, in_t_data, in_group)
      !
      ! ------------------------------------------------------------
      ! Create a new variable in the FieldArray list
      !
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(inout)         :: SELF
      CHARACTER(len=*), INTENT(in)             :: in_name     ! Variable name
      CHARACTER(len=*), INTENT(in)             :: in_long_name     ! Long name, mainly for output attributes
      CHARACTER(len=*), INTENT(in)             :: in_unit          ! Unit of the variable, e.g. "kg/kg"; used mainly for output attributes
      CHARACTER(len=*), INTENT(in)             :: in_dimension     ! String that gives the dimension environment for output (see ncio.f90)
      LOGICAL, INTENT(in)                      :: in_outputstatus
      CLASS(*), INTENT(in), POINTER            :: in_p_data   ! Polymorphic pointer to data (values)
      CLASS(*), INTENT(in), POINTER, OPTIONAL  :: in_t_data   ! - '' - (tendencies)
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: in_group
 
      ! Extend the variable list allocation in FieldArray
      CALL SELF%Extend()
      ! Pass the input data and parameters to ArrayElement constructor
      SELF%list(SELF%count) = ArrayElement(in_name,in_long_name,in_unit,in_dimension,in_outputstatus,  &
                                           in_p_data,in_t_data,in_group)

   END SUBROUTINE newField

   ! ------------------------------------------------------------

   SUBROUTINE Extend_FieldArray(SELF)
      !
      ! ----------------------------------------------------------
      ! Extend the memory allocation of the "list" in FieldArray.
      ! Mainly intended to be used by newField procedure
      !
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(inout) :: SELF

      TYPE(ArrayElement), ALLOCATABLE :: tmp(:)
      INTEGER :: ss

      IF (SELF%Initialized) THEN

         ALLOCATE(tmp(SELF%count+1))
         tmp(1:SELF%count) = SELF%list
         DEALLOCATE(SELF%list)
         CALL MOVE_ALLOC(tmp,SELF%list)
       
      ELSE
         ALLOCATE(SELF%list(SELF%count+1))
      END IF

      SELF%count = SELF%count + 1
      SELF%Initialized = .TRUE.

   END SUBROUTINE Extend_FieldArray

   ! -------------------------------------------------------------

   SUBROUTINE ExN_FieldArray(SELF,n)
      !
      ! ----------------------------------------------------------
      ! Extend the memory allocation of the "list" in FieldArray.
      ! Mainly intended to be used by newField procedure
      !
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(inout) :: SELF
      INTEGER, INTENT(in) :: n

      TYPE(ArrayElement), ALLOCATABLE :: tmp(:)
      INTEGER :: ss

      IF (SELF%Initialized) THEN

         ALLOCATE(tmp(SELF%count+n))
         tmp(1:SELF%count) = SELF%list
         DEALLOCATE(SELF%list)
         CALL MOVE_ALLOC(tmp,SELF%list)

      ELSE
         ALLOCATE(SELF%list(SELF%count+n))
      END IF

      SELF%count = SELF%count + n
      SELF%Initialized = .TRUE.

   END SUBROUTINE ExN_FieldArray


   SUBROUTINE getField(SELF,ind,out)
      !
      ! ---------------------------------------------
      ! Returns the ArrayElement instance for given
      ! index "ind" of the FieldArray list.
      !
      IMPLICIT NONE
      CLASS(FieldArray), TARGET, INTENT(in) :: SELF
      INTEGER, INTENT(in) :: ind
      TYPE(ArrayElement), POINTER :: out

      out => SELF%list(ind)

   END SUBROUTINE getField

   ! -----------------------------------------------------------------------------
   ! getGroup returns a FieldArray instance, which contains a subset of variables
   ! from the parent FieldArray instance that belong to the inquired group. Whether
   ! any variables meeting the criterion are found can be checked by the field
   ! FieldArray%Initialized.
   ! 
   SUBROUTINE getByGroup(SELF,groupname,FAout)
     IMPLICIT NONE
     CLASS(FieldArray), INTENT(in)  :: SELF
     CHARACTER(len=*), INTENT(in)   :: groupname
     TYPE(FieldArray), INTENT(out) :: FAout

     LOGICAL :: groupmask(SELF%count)
     INTEGER :: i
     
     FAout = FieldArray()
     
     groupmask = .FALSE.
     groupmask(:) = (SELF%list(:)%group == groupname)
     
     FAout%count = COUNT(groupmask)
     IF (FAout%count > 0) THEN
        ALLOCATE(FAout%list(FAout%count))
        FAout%list(:) = PACK(SELF%list(:),groupmask)     
        FAout%Initialized = .TRUE.
     END IF
        
   END SUBROUTINE getByGroup

   SUBROUTINE getByOutputstatus(SELF,FAout)
     IMPLICIT NONE
     CLASS(FieldArray), INTENT(in) :: SELF
     TYPE(FieldArray), INTENT(out) :: FAout

     LOGICAL :: mask(SELF%count)

     FAout = FieldArray()

     mask = .FALSE.
     mask(:) = (SELF%list(:)%outputstatus)

     FAout%count = COUNT(mask)
     IF (FAout%count > 0) THEN
        ALLOCATE(FAout%list(FAout%count))
        FAout%list(:) = PACK(SELF%list(:),mask)     
        FAout%Initialized = .TRUE.
     END IF
        
   END SUBROUTINE getByOutputstatus
   
   !
   ! ----------------------------------------------------------------------------
   ! getVarInst returns the polymorphic instance for the value (tend=1) or tendency (tend=2)
   ! within ArrayElement referred to by the variable name, or directly by its index.
   ! getVarInst overloads the procedures getVarInst_name and getVarInst_ind bound to
   ! FieldArray.
   !
   SUBROUTINE getVarInst_name(SELF,in_name,out,tend)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)    :: SELF
      CHARACTER(len=*), INTENT(in)     :: in_name  ! Variable name
      CLASS(*), INTENT(out), POINTER   :: out      ! Variable instance
      INTEGER, INTENT(in)              :: tend     ! Get the values (1; default) or tendencies (2)

      TYPE(ArrayElement), POINTER :: Element       ! Pointer to the FieldArray list element
      INTEGER :: ind

      CALL SELF%getFieldIndex(in_name,ind)
      CALL SELF%getField(ind,Element)
      IF (tend==1) THEN
         CALL Element%get_p(out)
      ELSE IF (tend==2) THEN
         CALL Element%get_t(out)
      END IF

      Element => NULL()

   END SUBROUTINE getVarInst_name
   ! -------
   ! -------
   SUBROUTINE getVarInst_ind(SELF,ind,out,tend)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)    :: SELF
      INTEGER, INTENT(in)              :: ind     ! Index of the variable in the FieldArray list
      CLASS(*), INTENT(out), POINTER   :: out     ! Variable instance
      INTEGER, INTENT(in)              :: tend    ! Get the values (1; default) or tendencies (2)

      TYPE(ArrayElement), POINTER :: Element      ! Pointer to the FieldArray list element

      CALL SELF%getField(ind,Element)
      IF (tend==1) THEN
         CALL Element%get_p(out)
      ELSE IF (tend==2) THEN
         CALL Element%get_t(out)
      END IF

      Element => NULL()

   END SUBROUTINE getVarInst_ind

   ! ---------------------------------------------------------------

   !
   ! --------------------------------------------------------------------
   ! getData provides the main user interface for getting the numerical
   ! data (value - tend=1; tendency - tend=2) for a variable corresponding
   ! to the given index or variable name. getData uses the getVarInst
   ! procedure and then essentially collapses the fetched polymorphic
   ! instance to the correct datatype inferred by the output variable "out".
   ! If the datatype of "out" does not correspond to existing variables
   ! given by "index" or "name", an unassociated pointer with error status
   ! is returned. getData currently overloads the procedures getData_1d,
   ! getData_2d, getData_3d and getData_4d bound to FieldArray.
   !
   SUBROUTINE getData_0d(SELF,tend,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray0d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tend   ! Get value (1; default) or tendency (2)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tend)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tend)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, return an error status LISAA TAMA
      SELECT TYPE(pp)
         TYPE IS (FloatArray0d)
            out => pp
      END SELECT

   END SUBROUTINE getData_0d
   
   ! ---------
   ! ---------
   SUBROUTINE getData_1d(SELF,tend,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray1d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tend   ! Get value (1; default) or tendency (2)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tend)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tend)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, return an error status LISAA TAMA
      SELECT TYPE(pp)
         TYPE IS (FloatArray1d)
            out => pp
      END SELECT

   END SUBROUTINE getData_1d
   ! ---------
   ! ---------
   SUBROUTINE getData_2d(SELF,tend,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray2d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tend   ! Get value (1; default) or tendency (2)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tend)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tend)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, return an error status LISAA TAMA
      SELECT TYPE(pp)
         TYPE IS (FloatArray2d)
            out => pp
      END SELECT

   END SUBROUTINE getData_2d
   ! ---------
   ! ---------
   SUBROUTINE getData_3d(SELF,tend,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray3d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tend   ! Get value (1; default) or tendency (2)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tend)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tend)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, return an error status LISAA TAMA
      SELECT TYPE(pp)
         TYPE IS (FloatArray3d)
            out => pp
      END SELECT

   END SUBROUTINE getData_3d
   ! ---------
   ! ---------
   SUBROUTINE getData_4d(SELF,tend,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray4d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tend   ! Get value (1; default) or tendency (2)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tend)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tend)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, return an error status LISAA TAMA
      SELECT TYPE(pp)
         TYPE IS (FloatArray4d)
            out => pp
      END SELECT

   END SUBROUTINE getData_4d

   ! ---------------------------------------------------------------

   SUBROUTINE getFieldIndex(SELF,in_name,ind)
      !
      ! ---------------------------------------------------------
      ! Returns the FieldArray index with name=in_name
      !
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in)  :: in_name
      INTEGER, INTENT(out)          :: ind

      INTEGER :: i

      i = 1
      DO
         IF (i > SELF%count) THEN
            i = 999
            EXIT
         END IF
         IF (SELF%list(i)%name == in_name) THEN
            EXIT
         END IF
         i = i + 1
      END DO
      ind = i
      
   END SUBROUTINE getFieldIndex

   ! ---------------------------------------------------------------

   SUBROUTINE destroy_FieldArray(SELF)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(inout) :: SELF
    
      DEALLOCATE(SELF%list)

   END SUBROUTINE destroy_FieldArray


   !
   ! -------------------------------------------------------
   ! PROCEDURES BOUND TO ArrayElement
   ! --------------------------------
   !
   !

  SUBROUTINE get_p(SELF,out)
      IMPLICIT NONE
      CLASS(ArrayElement), INTENT(in) :: SELF
      CLASS(*), INTENT(out), POINTER  :: out

      out => SELF%p
    
   END SUBROUTINE get_p

   ! --------------------------------------------------------------

   SUBROUTINE get_t(SELF,out)
      IMPLICIT NONE
      CLASS(ArrayElement), INTENT(in) :: SELF
      CLASS(*), INTENT(out), POINTER  :: out

      out => SELF%t

   END SUBROUTINE get_t


END MODULE classFieldArray
