MODULE classFieldArray
   !
   ! Written for UCLALES-SALSA (U-S).
   !
   ! Juha Tonttila, FMI, 2017
   !
   USE mo_structured_datatypes, ONLY : FloatArray0d, FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d

   IMPLICIT NONE

   INTEGER, PARAMETER :: undefined = 999 

   
   TYPE ArrayElement
      CHARACTER(len=50)  :: name          ! Short name: for output variable name and for fetching data using getData
      CHARACTER(len=100) :: long_name     ! Long name, mainly for output attributes
      CHARACTER(len=50)  :: unit          ! Unit of the variable, e.g. "kg/kg"; used mainly for output attributes
      CHARACTER(len=50)  :: dimension     ! String that gives the dimension environment for output (see ncio.f90)
      CHARACTER(len=50), ALLOCATABLE  :: group(:)         ! A group tag that can be used to fetch a list of certain type variables
      LOGICAL            :: outputstatus  ! TRUE: write this variable to an output file. FALSE: don't.

      ! Below are pointers to scalar variable arrays. In U-S, scalars have two definitions:
      ! the previous value (p) and the tendency (t). For vector variables the is also the
      ! current (c) state. These are given here as unlimited
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
      CLASS(*), POINTER :: c => NULL()
      
      CONTAINS
       
         PROCEDURE :: get_p
         PROCEDURE :: get_t
         PROCEDURE :: get_c
         
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

         PROCEDURE :: ExN_FieldArray         ! WHY IS THIS HERE??
         GENERIC   :: ExN => ExN_FieldArray

         PROCEDURE :: NewField

         PROCEDURE :: getField

         PROCEDURE :: getFieldIndex
         
         PROCEDURE :: getByGroup
         PROCEDURE :: getByOutputstatus
         
         PROCEDURE :: getVarInst_ind, getVarInst_name
         GENERIC   :: getVarInst => getVarInst_ind, getVarInst_name

         PROCEDURE :: getData_0d, getData_1d, getData_2d, getData_3d, getData_4d
         GENERIC   :: getData => getData_0d, getData_1d, getData_2d, getData_3d, getData_4d

         PROCEDURE :: getDimension
         
         PROCEDURE :: Exist
       
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
     FUNCTION ArrayElement_constructor(in_name,in_long_name,in_unit,in_dimension,in_outputstatus,   &
                                       in_p_data,in_t_data,in_c_data,in_group)
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
      CLASS(*), INTENT(in), POINTER, OPTIONAL   :: in_c_data      ! - '' - (current; for vectors)
      CHARACTER(len=*), INTENT(in), OPTIONAL    :: in_group(:)

      ArrayElement_constructor%name = in_name
      ArrayElement_constructor%long_name = in_long_name
      ArrayElement_constructor%unit = in_unit
      ArrayElement_constructor%dimension = in_dimension
      ArrayElement_constructor%outputstatus = in_outputstatus
      ArrayElement_constructor%p => in_p_data

      IF (PRESENT(in_group)) THEN
         ALLOCATE(ArrayElement_constructor%group(SIZE(in_group)))
         ArrayElement_constructor%group = in_group
      ELSE
         ! Set group as "default" if nothing given
         ALLOCATE(ArrayElement_constructor%group(1))
         ArrayElement_constructor%group = ['default']
      END IF
      IF (PRESENT(in_t_data)) ArrayElement_constructor%t => in_t_data
      IF (PRESENT(in_c_data)) ArrayElement_constructor%c => in_c_data
      
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
   SUBROUTINE newField(SELF,in_name,in_long_name,in_unit,in_dimension,in_outputstatus,   &
                       in_p_data,in_t_data,in_c_data,in_group)
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
      CLASS(*), INTENT(in), POINTER, OPTIONAL  :: in_c_data   ! - '' - (current)
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: in_group(:)
 
      ! Extend the variable list allocation in FieldArray
      CALL SELF%Extend()
      ! Pass the input data and parameters to ArrayElement constructor
      SELF%list(SELF%count) = ArrayElement(in_name,in_long_name,in_unit,in_dimension,in_outputstatus,  &
                                           in_p_data,in_t_data,in_c_data,in_group)

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

   !
   ! ---------------------------------------------------------
   ! Returns the FieldArray index with name=in_name
   !
   INTEGER FUNCTION getFieldIndex(SELF,in_name)
      CLASS(FieldArray), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in)  :: in_name
      INTEGER :: i

      i = 1
      DO
         IF (i > SELF%count) THEN
            i = undefined
            EXIT
         END IF
         IF (SELF%list(i)%name == in_name) THEN
            EXIT
         END IF
         i = i + 1
      END DO
      getFieldIndex = i
      
   END FUNCTION getFieldIndex
   !
   ! ---------------------------------------------
   ! Returns the ArrayElement instance for given
   ! index "ind" of the FieldArray list.
   !
   SUBROUTINE getField(SELF,out,ind,name)

      IMPLICIT NONE
      CLASS(FieldArray), TARGET, INTENT(in) :: SELF
      INTEGER, INTENT(in), OPTIONAL :: ind
      CHARACTER(len=*), INTENT(in), OPTIONAL :: name
      TYPE(ArrayElement), INTENT(out), POINTER :: out
      INTEGER :: lind
      lind=0
      IF (PRESENT(ind)) THEN
         lind = ind
      ELSE IF (PRESENT(name)) THEN
         lind = SELF%getFieldIndex(name)
      END IF
      IF (lind==0) &
           WRITE(*,*) "classFieldArray:getField -- WARNING: Variable not found ",name
      out => SELF%list(lind)

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
     DO i = 1,SELF%count
        groupmask(i) = (ANY(SELF%list(i)%group(:) == groupname))
     END DO
        
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

     IF (.NOT. SELF%Initialized) RETURN
     
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
   ! THESE SHOULD BE MADE PRIVATE, SINCE THE POLYMORPHIC RESULT REQUIRES MORE WORK BEFORE IT IS USEFUL.
   ! MAINLY THESE ARE NEEDED FOR THE GETDATA-PROCEDURES, WHICH ARE THE ONES THAT SHOULD BE USED.
   !
   SUBROUTINE getVarInst_name(SELF,in_name,out,tlev)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)    :: SELF
      CHARACTER(len=*), INTENT(in)     :: in_name  ! Variable name
      CLASS(*), INTENT(out), POINTER   :: out      ! Variable instance
      INTEGER, INTENT(in)              :: tlev     ! Get the values (1; default) or tendencies (2) or current (for vectors; 3)

      TYPE(ArrayElement), POINTER :: Element      ! Pointer to the FieldArray list element
      INTEGER :: ind

      out => NULL()
      Element => NULL()
      ind = SELF%getFieldIndex(in_name)
      IF (ind /= undefined ) &              ! If the variable is not defined, Element will remain unassociated
           CALL SELF%getField(Element,ind=ind)
      IF (tlev==1 .AND. ASSOCIATED(Element)) THEN ! If Element is unassociated, out will remain unassociated as well
         CALL Element%get_p(out)
      ELSE IF (tlev==2 .AND. ASSOCIATED(Element)) THEN
         CALL Element%get_t(out)
      ELSE IF (tlev==3 .AND. ASSOCIATED(Element)) THEN
         CALL Element%get_c(out)
      END IF

      Element => NULL()

   END SUBROUTINE getVarInst_name
   ! -------
   ! -------
   SUBROUTINE getVarInst_ind(SELF,ind,out,tlev)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)    :: SELF
      INTEGER, INTENT(in)              :: ind     ! Index of the variable in the FieldArray list
      CLASS(*), INTENT(out), POINTER   :: out     ! Variable instance
      INTEGER, INTENT(in)              :: tlev    ! Get the values (1; default) or tendencies (2) or current (for vectors; 3)

      TYPE(ArrayElement), POINTER :: Element      ! Pointer to the FieldArray list element

      Element => NULL()
      out => NULL()
      IF (ind <= SELF%count .AND. ind >= 1 .AND. ind /= undefined) & ! If index is bad, Element and out will remain unassociated
           CALL SELF%getField(Element,ind=ind)      
      IF (tlev==1 .AND. ASSOCIATED(Element)) THEN
         CALL Element%get_p(out)
      ELSE IF (tlev==2 .AND. ASSOCIATED(Element)) THEN
         CALL Element%get_t(out)
      ELSE IF (tlev==3 .AND. ASSOCIATED(Element)) THEN
         CALL Element%get_c(out)
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
   SUBROUTINE getData_0d(SELF,tlev,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray0d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tlev   ! Get value (1; default) or tendency (2) or current (for vectors; 3)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      out => NULL()
      
      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tlev)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tlev)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, out will remain unassociated
      IF (ASSOCIATED(pp)) THEN
         SELECT TYPE(pp)
         TYPE IS (FloatArray0d)
            out => pp
         END SELECT
      END IF

   END SUBROUTINE getData_0d
   
   ! ---------
   ! ---------
   SUBROUTINE getData_1d(SELF,tlev,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray1d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tlev   ! Get value (1; default) or tendency (2) or current (for vectors; 3)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      out => NULL()
      
      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tlev)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tlev)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, out will remain unassociated
      IF (ASSOCIATED(pp)) THEN
         SELECT TYPE(pp)
         TYPE IS (FloatArray1d)
            out => pp
         END SELECT
      END IF
         
   END SUBROUTINE getData_1d
   ! ---------
   ! ---------
   SUBROUTINE getData_2d(SELF,tlev,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray2d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tlev   ! Get value (1; default) or tendency (2) or current (for vectors; 3)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      out => NULL()
      
      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tlev)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tlev)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, out will remain unassociated
      IF (ASSOCIATED(pp)) THEN
         SELECT TYPE(pp)
         TYPE IS (FloatArray2d)
            out => pp
         END SELECT
      END IF
         
   END SUBROUTINE getData_2d
   ! ---------
   ! ---------
   SUBROUTINE getData_3d(SELF,tlev,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray3d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tlev   ! Get value (1; default) or tendency (2) or current (for vectors; 3)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      out => NULL()
      
      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tlev)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tlev)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, out will remain unassociated
      IF (ASSOCIATED(pp)) THEN
         SELECT TYPE(pp)
         TYPE IS (FloatArray3d)
            out => pp
         END SELECT
      END IF
         
   END SUBROUTINE getData_3d
   ! ---------
   ! ---------
   SUBROUTINE getData_4d(SELF,tlev,out,index,name)
      IMPLICIT NONE
      CLASS(FieldArray), INTENT(in)            :: SELF
      TYPE(FloatArray4d), INTENT(out), POINTER :: out    ! Data container instance
      INTEGER, INTENT(in), OPTIONAL            :: index  ! Index of the variable in the FieldArray list
      CHARACTER(len=*), INTENT(in), OPTIONAL   :: name   ! Name of the variable
      INTEGER, INTENT(in)                      :: tlev   ! Get value (1; default) or tendency (2) or current (for vectors; 3)

      CLASS(*), POINTER :: pp  ! Polymorphic pointer to the variable instance in ArrayElement

      out => NULL()
      
      IF (PRESENT(index)) THEN
         CALL SELF%getVarInst(index,pp,tlev)
      ELSE IF (PRESENT(name)) THEN
         CALL SELF%getVarInst(name,pp,tlev)
      END IF

      ! Collapse to inquired datatype. If the in_name or ind does not point to
      ! variable of this datatype, out will remain unassociated
      IF (ASSOCIATED(pp)) THEN
         SELECT TYPE(pp)
         TYPE IS (FloatArray4d)
            out => pp
         END SELECT
      END IF
         
   END SUBROUTINE getData_4d

   ! -----------------------------------------------------------------

   CHARACTER(len=50) FUNCTION getDimension(SELF,iname)
     CLASS(FieldArray), INTENT(in) :: SELF
     CHARACTER(len=*), INTENT(in) :: iname
     INTEGER :: ind
     
     IF ( SELF%Exist(iname) ) THEN
        ind = SELF%getFieldIndex(iname)
        getDimension = SELF%list(ind)%dimension
     ELSE
        getDimension = ''
     END IF
     
   END FUNCTION getDimension
   
   ! -----------------------------------------------------------------
   
   LOGICAL FUNCTION Exist(SELF,iname)
     CLASS(FieldArray), INTENT(in) :: SELF
     CHARACTER(len=*), INTENT(in) :: iname

     INTEGER :: i
     
     ! Check if variable with iname exists
     Exist = .FALSE.
     DO i = 1,SELF%count
        IF (SELF%list(i)%name == iname) THEN
           Exist = .TRUE.
           EXIT
        END IF
     END DO     
   END FUNCTION Exist
   
   
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

   ! --------------------------------------------------------------

   SUBROUTINE get_c(SELF,out)
     IMPLICIT NONE
     CLASS(ArrayElement), INTENT(in) :: SELF
     CLASS(*), INTENT(out), POINTER  :: out

     out => SELF%c

   END SUBROUTINE get_c

   

END MODULE classFieldArray
