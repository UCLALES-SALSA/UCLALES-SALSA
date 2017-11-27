MODULE classSpecies
  ! THIS WILL REPLACE CLASS_COMPONENTINDEX
  ! THIS IS NOT YET USED ANYWHERE
  ! THIS IS NOT YET TESTED
  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "classSpecies"

  INTEGER, PARAMETER :: maxspec = 7  ! Maximum number of aerosol species, excluding water
  ! Properties for all possible species. The order of the parameters in 
  ! the vectors below is consistent
  CHARACTER(len=3), TARGET   :: allNames(maxspec+1) = ['SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O']
  REAL, TARGET               :: allMM(maxspec+1) = [98.08e-3, 150.e-3, 12.e-3, 100.e-3, 58.44e-3, 62.01e-3, 18.04e-3, 18.016e-3] ! Molecular mass 
  REAL, TARGET               :: allRho(maxspec+1) = [1830., 2000., 2000., 2650., 2165., 1479., 1530., 1000.]                     ! Density
  REAL, TARGET               :: allDiss(maxspec+1) = [3., 1., 0., 0., 2., 1., 1., 1.]                                            ! Dissociation factor

  !
  ! ----------------------------------------------------------------------
  ! Type Species holds a collection of names and properties for the currently 
  ! available aerosol species (including water). Densities, molar masses etc.
  ! in mo_submctl will become deprecated.
  !
  TYPE Species


     ! Logical vector telling which species are used
     LOGICAL            :: used(maxspec+1) = [.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.]
     INTEGER            :: Nused = 0    ! Number of used species, includes water
     INTEGER            :: allInd(maxspec+1) = [0,0,0,0,0,0,0,0]
     ! These are corresponding truncated lists of properties for the species that are used.
     ! They are generated at initialization using logical used and PACK
     CHARACTER(len=3), ALLOCATABLE :: names(:)
     INTEGER         , ALLOCATABLE :: ind(:)
     REAL            , ALLOCATABLE :: MM(:)
     REAL            , ALLOCATABLE :: rho(:)
     REAL            , ALLOCATABLE :: diss(:)

     ! Some additional values
     REAL    :: mas = 132.14e-3  ! Molar mass of ammonium sulphate
     REAL    :: rhoas = 1770.    ! Density of ammonium sulphate
     REAL    :: rhoic = 917.,  & ! Assumed densities of ice and snow
                rhosn = 300.

     ! Pointers to individual names for specific use. Note: is you use these, no
     ! automatic checks are performed whether the species is used. This can be done
     ! manually with isUsed-method.
     CHARACTER(len=3), POINTER :: nsu => NULL(), noc, nbc, ndu, nss, nno, nnh, nwa
     
     REAL, POINTER :: msu => NULL(), moc, mbc, mdu, mss, mno, mnh, mwa
     REAL, POINTER :: rhosu, rhooc, rhobc, rhodu, rhoss, rhono, rhonh, rhowa
     
     CONTAINS

       PROCEDURE :: getIndex
       PROCEDURE :: getNSpec
       PROCEDURE :: isUsed
       PROCEDURE :: getRhoByName, getRhoByIndex
       GENERIC   :: getRho => getRhoByName, getRhoByIndex
       
  END TYPE Species

  INTERFACE Species
       PROCEDURE :: cnstr
  END INTERFACE Species

  CONTAINS

    FUNCTION cnstr(nlist,listcomp)
      IMPLICIT NONE
      TYPE(Species) :: cnstr
      INTEGER, INTENT(in) :: nlist                     ! Number of aerosol species to be used
      CHARACTER(len=3), INTENT(in) :: listcomp(maxspec)  ! Names of the aerosol species to be used

      INTEGER :: i,c

      ! +1 for water
      cnstr%Nused = nlist+1

      DO i = 1, nlist

         DO c = 1,maxspec
            IF ( allNames(c) == listcomp(i) ) THEN
               cnstr%allInd(c) = i
               cnstr%used(c) = .TRUE.
            END IF
         END DO
         
      END DO
      ! For water
      cnstr%allInd(maxspec+1) = nlist+1
      cnstr%used(maxspec+1) = .TRUE.

      ALLOCATE(cnstr%names(cnstr%Nused), cnstr%ind(cnstr%Nused),  &
               cnstr%MM(cnstr%Nused), cnstr%rho(cnstr%Nused),     &
               cnstr%diss(cnstr%Nused))

      cnstr%names = PACK(allNames, cnstr%used)
      cnstr%ind = PACK(cnstr%allInd, cnstr%used)
      cnstr%MM = PACK(allMM, cnstr%used)
      cnstr%rho = PACK(allRho, cnstr%used)
      cnstr%diss = PACK(allDiss, cnstr%used)

      cnstr%nsu  => allNames(1)
      cnstr%noc  => allNames(2)
      cnstr%nbc  => allNames(3)
      cnstr%ndu  => allNames(4)
      cnstr%nss  => allNames(5)
      cnstr%nno  => allNames(6)
      cnstr%nnh  => allNames(7)
      cnstr%nwa => allNames(8)

      cnstr%msu  => allMM(1)
      cnstr%moc  => allMM(2)
      cnstr%mbc  => allMM(3)
      cnstr%mdu  => allMM(4)
      cnstr%mss  => allMM(5)
      cnstr%mno  => allMM(6)
      cnstr%mnh  => allMM(7)
      cnstr%mwa => allMM(8)

      cnstr%rhosu  => allRho(1)
      cnstr%rhooc  => allRho(2)
      cnstr%rhobc  => allRho(3)
      cnstr%rhodu  => allRho(4)
      cnstr%rhoss  => allRho(5)
      cnstr%rhono  => allRho(6)
      cnstr%rhonh  => allRho(7)
      cnstr%rhowa => allRho(8)


    END FUNCTION cnstr

    ! -----------------------------------------

    INTEGER FUNCTION getIndex(SELF,incomp)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      INTEGER :: i

      IF ( SELF%isUsed(incomp) ) THEN

         i = 1
         DO WHILE ( (SELF%names(i) /= incomp) )
            i = i + 1
         END DO
         GetIndex = i
      ELSE
         STOP 'classSpecies: getIndex: FAILED, no such component - '
      END IF

      RETURN

    END FUNCTION getIndex

    ! ---------------------------------------------

    INTEGER FUNCTION getNSpec(SELF)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      
      ! +1 for water
      getNSpec = SELF%Nused

    END FUNCTION getNSpec

    ! ------------------------------------------------

    LOGICAL FUNCTION isUsed(SELF,incomp)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp

      isUsed = .FALSE.

      IF ( ANY(SELF%names == incomp) ) isUsed = .TRUE.

    END FUNCTION isUsed

    ! -------------------------------------------------

    REAL FUNCTION getRhoByIndex(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      INTEGER, INTENT(in) :: nn
      
      getRhoByIndex = SELF%rho(nn)

    END FUNCTION getRhoByIndex
    ! ---------------------------------------
    REAL FUNCTION getRhoByName(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: nn

      INTEGER :: ii

      ii = SELF%getIndex(nn)
      getRhoByName = SELF%getRhoByIndex(ii)

    END FUNCTION getRhoByName

    ! -----------------------------------------------


END MODULE classSpecies
