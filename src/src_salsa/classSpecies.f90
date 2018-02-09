MODULE classSpecies

  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER, PRIVATE :: global_name = "classSpecies"

  INTEGER, PARAMETER, PUBLIC :: maxspec = 7  ! Maximum number of aerosol species, excluding water
  ! Properties for all possible species. The order of the parameters in 
  ! the vectors below is consistent.
  ! The last three: liq,ice,snow
  CHARACTER(len=3), TARGET, PRIVATE   :: allNames(maxspec+1) = ['SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O']
  REAL, TARGET, PRIVATE               :: allMM(maxspec+1) = [98.08e-3, 150.e-3, 12.e-3, 100.e-3,    &
                                                             58.44e-3, 62.01e-3, 18.04e-3, 18.016e-3] ! Molecular mass 
  REAL, TARGET, PRIVATE               :: allRho(maxspec+1) = [1830., 2000., 2000., 2650., 2165., 1479., 1530., 1000.]                   ! Density
  REAL, TARGET, PRIVATE               :: allDiss(maxspec+1) = [3., 1., 0., 0., 2., 1., 1., 1.]                                              ! Dissociation factor
  REAL, TARGET, PRIVATE    :: auxrhoic = 917.     ! ice
  REAL, TARGET, PRIVATE    :: auxrhosn = 300.     ! snow

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
     ! They are generated at initialization using logical used and PACK. The lists with suffix "Fixed"
     ! are given for the compounds in the same order as in the global list "allNames" above. However, user
     ! can define the active compounds in namelist_salsa in any order, which will also be the order in
     ! which the compounds will be in the UCLALES-SALSA mass mixing ratio arrays (except for water, which will
     ! *always* be the last one). Thus, the lists below *without* the suffix "Fixed" are sorted to give the 
     ! properties in the same order as the compounds are in given in the namelist_salsa and in the mass 
     ! mixing ratio arrays, and these are the ones that should be generally used.
     ! Should/can the Fixed-lists be PRIVATE?
     CHARACTER(len=3), ALLOCATABLE :: namesFixed(:), names(:)
     INTEGER         , ALLOCATABLE :: indFixed(:), ind(:)
     REAL            , ALLOCATABLE :: MMFixed(:), MM(:)
     REAL            , ALLOCATABLE :: rholiqFixed(:), rhoiceFixed(:), rhosnowFixed(:),   &
                                      rholiq(:), rhoice(:), rhosnow(:)
     REAL            , ALLOCATABLE :: dissFixed(:), diss(:)

     ! Some additional values
     REAL    :: mas = 132.14e-3  ! Molar mass of ammonium sulphate
     REAL    :: rhoas = 1770.    ! Density of ammonium sulphate


     ! Pointers to individual names for specific use. Note: if you use these, no
     ! automatic checks are performed whether the species is used. This can be done
     ! manually with isUsed-method.
     CHARACTER(len=3), POINTER :: nsu => NULL(), &
                                  noc => NULL(), &
                                  nbc => NULL(), &
                                  ndu => NULL(), &
                                  nss => NULL(), &
                                  nno => NULL(), &
                                  nnh => NULL(), &
                                  nwa => NULL()

     
     REAL, POINTER :: msu => NULL(),   &
                      moc => NULL(),   &
                      mbc => NULL(),   &
                      mdu => NULL(),   &
                      mss => NULL(),   &
                      mno => NULL(),   &
                      mnh => NULL(),   &
                      mwa => NULL()

     REAL, POINTER :: rhosu => NULL(), &
                      rhooc => NULL(), &
                      rhobc => NULL(), &
                      rhodu => NULL(), &
                      rhoss => NULL(), &
                      rhono => NULL(), &
                      rhonh => NULL(), &
                      rhowa => NULL(), &
                      rhoic => NULL(), &
                      rhosn => NULL()
     
     CONTAINS

       PROCEDURE :: getIndex
       PROCEDURE :: getNSpec
       PROCEDURE :: isUsed
       PROCEDURE :: getRhoByName, getRhoByIndex
       GENERIC   :: getRho => getRhoByName, getRhoByIndex
       PROCEDURE :: getDissByName, getDissByIndex
       GENERIC   :: getDiss => getDissByName, getDissByIndex
       PROCEDURE :: getMMByName, getMMByIndex
       GENERIC   :: getMM => getMMByName, getMMByIndex
       
       PROCEDURE, PRIVATE :: sortProperties

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
      
      INTEGER :: nalloc

      INTEGER :: i,c

      ! Define which compounds are used based on namelist-salsa information
      ! and store their corresponding indices in the mass arrays
      DO i = 1, nlist
         DO c = 1,maxspec
            IF ( allNames(c) == listcomp(i) ) THEN
               cnstr%allInd(c) = i
               cnstr%used(c) = .TRUE.
            END IF
         END DO
      END DO

      ! Number of active compounds.
      ! +1 for water. Even though there are 3 water "species" (liq,ice,snow), they are in different arrays 
      ! and each has only one of the water phases.
      cnstr%Nused = nlist+1
         
      ! Add this stuff for water, which is always used (and thus not specified in the namelist)
      cnstr%allInd(maxspec+1) = nlist+1
      cnstr%used(maxspec+1) = .TRUE.
 
      ! Allocate the truncated property lists
      ALLOCATE(cnstr%namesFixed(cnstr%Nused), cnstr%names(cnstr%Nused),          &
               cnstr%indFixed(cnstr%Nused), cnstr%ind(cnstr%Nused),              &
               cnstr%MMFixed(cnstr%Nused), cnstr%MM(cnstr%Nused),                &
               cnstr%rholiqFixed(cnstr%Nused), cnstr%rholiq(cnstr%Nused),        &
               cnstr%rhoiceFixed(cnstr%Nused), cnstr%rhoice(cnstr%Nused),        &
               cnstr%rhosnowFixed(cnstr%Nused), cnstr%rhosnow(cnstr%Nused),      &
               cnstr%dissFixed(cnstr%Nused), cnstr%diss(cnstr%Nused)             )

      ! Truncate the property lists in the order given by the global field "allNames"
      ! in the beginning of this class.
      cnstr%namesFixed = PACK(allNames, cnstr%used)
      cnstr%indFixed = PACK(cnstr%allInd, cnstr%used)
      cnstr%MMFixed = PACK(allMM, cnstr%used)
      cnstr%dissFixed = PACK(allDiss, cnstr%used)
      cnstr%rholiqFixed = PACK(allRho, cnstr%used)
      ! First, use the same arrays for densities in ice and snow arrays
      cnstr%rhoiceFixed = cnstr%rholiqFixed
      cnstr%rhosnowFixed = cnstr%rholiqFixed
      ! Second, replace the water densities with appropriate values
      cnstr%rhoiceFixed(cnstr%Nused) = auxrhoic
      cnstr%rhosnowFixed(cnstr%Nused) = auxrhosn
      
      ! Make another set of truncated property lists, where the values are sorted to the same order
      ! as the compounds appear in the UCLALES-SALSA mass arrays
      CALL cnstr%sortProperties()

      ! Make separate named pointers to properties for use in some special cases (mainly some hard-coded things)
      cnstr%nsu  => allNames(1)
      cnstr%noc  => allNames(2)
      cnstr%nbc  => allNames(3)
      cnstr%ndu  => allNames(4)
      cnstr%nss  => allNames(5)
      cnstr%nno  => allNames(6)
      cnstr%nnh  => allNames(7)
      cnstr%nwa  => allNames(8)

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
      cnstr%rhowa  => allRho(8)
      cnstr%rhoic  => auxrhoic
      cnstr%rhosn  => auxrhosn

    END FUNCTION cnstr

    ! -----------------------------------------

    SUBROUTINE sortProperties(SELF)
      IMPLICIT NONE
      ! -------------------------------------------------------
      ! PRIVATE
      ! Return the vectors of properties sorted to the same order 
      ! in which the compounds are initialized in the UCLALES-SALSA 
      ! mass arrays.
      ! ---------------------------------------------
      CLASS(Species), INTENT(inout) :: SELF

      INTEGER :: ii,jj

      SELF%rholiq(:) = 0.
      SELF%rhoice(:) = 0.
      SELF%rhosnow(:) = 0.
      SELF%names(:) = '   '
      SELF%diss(:) = 0.
      SELF%MM(:) = 0.
      DO ii = 1, SELF%Nused
         jj = SELF%indFixed(ii)  ! This is the index the compound number ii has in the mass arrays
         SELF%rholiq(jj) = SELF%rholiqFixed(ii)
         SELF%rhoice(jj) = SELF%rhoiceFixed(ii)
         SELF%rhosnow(jj) = SELF%rhosnowFixed(ii)
         SELF%names(jj) = SELF%namesFixed(ii)
         SELF%diss(jj) = SELF%dissFixed(ii)
         SELF%MM(jj) = SELF%MMFixed(ii)
         SELF%ind(jj) = SELF%indFixed(ii) ! This is of course trivial but can be used to check everything's ok
      END DO

    END SUBROUTINE sortProperties

    ! ---------------------------------------------------------------------------

    INTEGER FUNCTION getIndex(SELF,incomp,notFoundValue)
      ! ------------------------------------------------------------
      ! Return the SALSA mass array index of the compound by name
      ! ------------------------------------------------------------
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: incomp
      INTEGER, OPTIONAL, INTENT(in) :: notFoundValue ! by default 0

      INTEGER :: i

      IF ( SELF%isUsed(incomp) ) THEN

         i = 1
         DO WHILE ( (SELF%names(i) /= incomp) ) ! Note: using the sorted list of names! (-> SALSA mass array index)
            i = i + 1
         END DO
         getIndex = i
      ELSE
         IF (PRESENT(notFoundValue)) THEN
            getIndex = notFoundValue
         ELSE
            getIndex = 0
         END IF
      END IF

      RETURN

    END FUNCTION getIndex

    ! ---------------------------------------------

    INTEGER FUNCTION getNSpec(SELF,type)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=3), INTENT(in), OPTIONAL :: type

      LOGICAL :: switch

      IF (PRESENT(type)) THEN
         IF (type == 'wet') THEN
            switch = .TRUE.
         ELSE IF (type == 'dry') THEN
            switch = .FALSE.
         END IF
      ELSE
         switch = .TRUE.
      END IF

      IF (switch) THEN
         ! include water
         getNSpec = SELF%Nused
      ELSE
         ! Include only aerosol
         getNSpec = SELF%Nused - 1
      END IF
        
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

    REAL FUNCTION getRhoByIndex(SELF,nn,wat)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      INTEGER, INTENT(in) :: nn
      INTEGER, INTENT(in), OPTIONAL :: wat ! Density of water differes according to phase: 
                                           ! 1: liquid, 2: ice, 3:snow. 1 is the default      
      IF ( PRESENT(wat) ) THEN
         IF (wat == 1) THEN
            getRhoByIndex = SELF%rholiq(nn)
         ELSE IF (wat == 2) THEN
            getRhoByIndex = SELF%rhoice(nn)
         ELSE IF (wat == 3 ) THEN
            getRhoByIndex = SELF%rhosnow(nn)
         END IF
      ELSE IF ( .NOT. PRESENT(wat) ) THEN
         ! By default use liquid array
         getRhoByIndex = SELF%rholiq(nn)
      END IF

    END FUNCTION getRhoByIndex
    ! ---------------------------------------
    REAL FUNCTION getRhoByName(SELF,nn,wat)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: nn
      INTEGER, INTENT(in), OPTIONAL :: wat ! Density of water differes according to phase: 
                                           ! 1: liquid, 2: ice, 3:snow. 1 is the default
      INTEGER :: ii

      ii = SELF%getIndex(nn)
      getRhoByName = SELF%getRhoByIndex(ii,wat)

    END FUNCTION getRhoByName

    ! ---------------------------------------

    REAL FUNCTION getDissByIndex(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      INTEGER, INTENT(in) :: nn

      getDissByIndex = SELF%diss(nn)

    END FUNCTION getDissByIndex
    ! ------------------------------------------------
    REAL FUNCTION getDissByName(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: nn
      
      INTEGER :: ii

      ii = SELF%getIndex(nn)
      getDissByName = SELF%getDissByIndex(ii)

    END FUNCTION getDissByName

    ! ------------------------------------------------------------

    REAL FUNCTION getMMByIndex(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      INTEGER, INTENT(in) :: nn
      
      getMMByIndex = SELF%MM(nn)

    END FUNCTION getMMByIndex 
    ! -------------------------------------------
    REAL FUNCTION getMMByName(SELF,nn)
      IMPLICIT NONE
      CLASS(Species), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: nn

      INTEGER :: ii

      ii = SELF%getIndex(nn)
      getMMByName = SELF%getMMByIndex(ii)

    END FUNCTION getMMByName

    ! -----------------------------------------------

    CHARACTER(len=3) FUNCTION getName(SELF,nn)
      IMPLICIT NONE
      ! -----------------------------------------------------------
      ! Gets species name by index from the subset of used species
      ! -----------------------------------------------------------
      CLASS(species), INTENT(in) :: SELF
      INTEGER, INTENT(in) :: nn

      getName = SELF%names(nn)

    END FUNCTION getName


END MODULE classSpecies
