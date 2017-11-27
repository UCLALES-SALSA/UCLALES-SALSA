MODULE classComponentIndex
  IMPLICIT NONE

  !
  ! ---------------------------
  ! Main class ComponentIndex
  !
  TYPE ComponentIndex
     INTEGER :: ncomp
     INTEGER, ALLOCATABLE :: ind(:)
     CHARACTER(len=3), ALLOCATABLE :: comp(:)
     
     CONTAINS

       PROCEDURE :: getIndex
       PROCEDURE :: getNComp
       PROCEDURE :: isUsed

  END TYPE ComponentIndex
  
  INTERFACE ComponentIndex
     PROCEDURE :: ComponentIndex_constructor
  END INTERFACE ComponentIndex

  ! --------------------------

  CONTAINS
  
    ! 
    ! --------------------------
    ! CONSTRUCTOR
    !
    FUNCTION ComponentIndex_constructor(ncomp, nlist, listcomp)
      IMPLICIT NONE
      TYPE(ComponentIndex) :: ComponentIndex_constructor 
      INTEGER, INTENT(in) :: ncomp, nlist
      CHARACTER(len=3), INTENT(in) :: listcomp(nlist)

      INTEGER :: i, jj

      ComponentIndex_constructor%ncomp = ncomp
      ALLOCATE(ComponentIndex_constructor%ind(ncomp),   &
               ComponentIndex_constructor%comp(ncomp)   )

      DO i = 1,ncomp
         ComponentIndex_constructor%ind(i) = i
      END DO

      jj = 1
      DO i = 1,nlist
         IF (listcomp(i) == '') CYCLE
         ComponentIndex_constructor%comp(jj) = listcomp(i)
         jj = jj+1
      END DO

    END FUNCTION ComponentIndex_constructor

    ! ---------------------------------

    INTEGER FUNCTION getIndex(SELF,incomp)
      IMPLICIT NONE
      CLASS(ComponentIndex), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      INTEGER :: i

      IF ( ANY(SELF%comp == incomp) ) THEN

         i = 1
         DO WHILE ( (SELF%comp(i) /= incomp) )
            i = i + 1
         END DO
         GetIndex = i
      ELSE IF ( incomp == 'H2O' ) THEN
         GetIndex = SELF%ncomp + 1
      ELSE
         STOP 'getIndex: FAILED, no such component - '
      END IF

      RETURN

    END FUNCTION getIndex

    ! -------------------------------------

    INTEGER FUNCTION getNcomp(SELF)
      CLASS(ComponentIndex), INTENT(in) :: SELF

      GetNcomp = SELF%ncomp+1

      RETURN

    END FUNCTION getNcomp
      
    ! -------------------------------------

    LOGICAL FUNCTION isUsed(SELF,incomp)
      IMPLICIT NONE
      CLASS(ComponentIndex), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      IF ( ANY(SELF%comp == incomp) ) THEN
         IsUsed = .TRUE.
      ELSE
         IsUsed = .FALSE.
      END IF

      RETURN

    END FUNCTION isUsed
      




END MODULE classComponentIndex
