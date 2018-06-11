MODULE classSpecies

  ! -------------------------------------------------------------------------------------------------------------------------------------
  ! Class Species
  ! --------------
  ! This class contains and organizes physical parameter, array indices and active state of
  ! the chemical compounds in SALSA. It is set up and initialized in mo_salsa_init after 
  ! reading the SALSA namelist. 
  ! 
  ! Basically, the class contains collocated arrays for the compound names 
  ! ('SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O'), logicals for whether a compound is 
  ! configured active, indices in UCLALES and SALSA mass arrays, molecular masses, densities (separate
  ! arrays for liquid, ice and snow particle types solely because of differing H2O density) and
  ! the dissolution factors (0 for insoluble).
  !
  !
  ! GLOBAL FIELDS
  ! --------------
  ! First, there are global (PRIVATE) Master arrays that contain the above information for all possible
  ! compounds in SALSA. These are used just to store the information as a source for the class bound
  ! variables, that are the ones to be used. 
  ! ***
  ! NOTE that the order in which the compounds are given in the PRIVATE master arrays is NOT the order 
  ! in which they will appear in UCLALES and SALSA mass arrays!! 
  ! ***
  ! The order of the compounds in the mass arrays will be the same in which the active compounds are defined 
  ! in the SALSA namelist!! The class bound variables will give the indexing as well as ready-made sorted lists 
  ! of the compound properties (also truncated to contain only the active compounds). Examples of how this works
  ! are given below.
  !
  !
  ! CLASS BOUND FIELDS
  ! --------------------
  ! The first class bound variable is the logical active state "used", which corresponds to "allNames".
  ! Compounds configured active in SALSA namelist get the value TRUE. "Nused" stores the number of 
  ! active compounds. "allInd" (despite the prefix "all" which is otherwise reserved for the global master arrays)
  ! will contain the indices of the species listed in "allNames" based on the order the active compounds
  ! are listed in the SALSA namelist. Nonactive compounds will get the index 0. allInd should be considered
  ! a master array, but it is defined anyway as class bound, since it depends on the model configuration.
  ! 
  ! Next, there are the trucated class bound arrays for compound names, indices, molecular masses, densities
  ! and dissolution factors. As mentioned, trucated means that they contain only the active compounds. Each of 
  ! these have two versions: ones with suffix "Native" and ones without. The "Native" versions give the information
  ! in the order in which the global array allNames is defined. The Non-"Native" version gives the information 
  ! in the order in which the compounds are given in the SALSA namelist and in the model mass arrays (thus, these
  ! are the most important ones).
  !
  ! For example: the SALSA namelist defines SO4, DU and SS as active, so that (note that I mix the order a bit)
  !     listspec = 'DU','SS','SO4','','','',''
  ! Then, the class variable "used" will be (H2O is prescribed always active and will always be the last compound in ALL arrays)
  !   >>  used = [.TRUE.,.FALSE.,.FALSE.,.TRUE.,.TRUE.,.FALSE.,.FALSE.,'TRUE']
  ! and "allInd" will be
  !   >>  allInd = [3,0,0,1,2,0,0,4]
  ! Further,
  !   >> namesNative = ['SO4','DU','SS','H2O']
  ! and
  !   >> indNative = [3,1,2,4]
  ! Then, the arrays sorted to the order of the SALSA namelist should be
  !   >> names = ['DU','SS','SO4','H2O']
  ! and (trivial, but good for checking)
  !   >> ind = [1,2,3,4]
  !
  ! For many (most?) situations, the easiest (and the fastest) way is to use direct access to the class variables.
  ! It is presumably also best to pretty much always use the sorted (Non-"Native") versions of the arrays.
  !
  ! Continuing with the class definition, there are additional property arrays for subsets of compounds, at the 
  ! moment for soluble and insoluble compounds, and variables giving the number of the compounds in both categories.
  ! If no (in)soluble compounds are specified active, the pointers are in nullified state. It will be easiest to
  ! check for this by requiring the number of compounds in each category to be  > 0. These are defined only according
  ! to the namelist sorted Non-"Native" form. 
  !
  ! Last, there are named variables for specific compound names, molar masses and densities, which are just used for
  ! some hardcoded parts. These are defined as pointers just so that there's no unnecessary copying of data. They are 
  ! associated with the global master arrays.
  !
  ! In addition, there is a selection of class bound procedures, that can be used to access the data, if its most convenient.
  ! As saíd, it is always fastest and often easiest to directly access the class variables, but sometimes the code will stay cleaner and
  ! more readable using the class procedures. The currently available list of procedures is definitely not exhaustive, and there's 
  ! plenty of room for extension.
  !
  ! The routines are pretty self explanatory, but will add a documentation about those also later...
  ! ---------------------------------------------------------------------------------------------------------------------------------



  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER, PRIVATE :: global_name = "classSpecies"

  INTEGER,           PARAMETER, PUBLIC  :: maxspec = 7  ! Maximum number of aerosol species, excluding water
  INTEGER,           PARAMETER, PRIVATE :: maxins = 2   ! Maximum number of insoluble compounds
  INTEGER,           PARAMETER, PRIVATE :: maxsol = 5   ! Maximum number of soluble compounds
  
  ! -------------------------------------------------------------------------------------------------------------------------
  ! Master arrays. 
  ! The order of the parameters in the vectors below is consistent. The last three: liq,ice,snow
  CHARACTER(len=3), TARGET, PRIVATE   :: allNames(maxspec+1)  = ['SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O']      ! Names of all possible compounds
  CHARACTER(len=3), TARGET, PRIVATE   :: allNamesInsoluble(maxins) = ['BC ','DU ']                                     ! Names of all possible insoluble compounds
  CHARACTER(len=3), TARGET, PRIVATE   :: allNamesSoluble(maxsol)   = ['SO4','OC ','SS ','NO ','NH ']                   ! Names of all possible soluble compounds
  REAL, TARGET, PRIVATE               :: allMM(maxspec+1)   = [98.08e-3, 150.e-3, 12.e-3, 100.e-3,    &                 
                                                               58.44e-3, 62.01e-3, 18.04e-3, 18.016e-3]                ! Molecular masses 
  REAL, TARGET, PRIVATE               :: allRho(maxspec+1)  = [1830., 2000., 2000., 2650., 2165., 1479., 1530., 1000.] ! Densities
  REAL, TARGET, PRIVATE               :: allDiss(maxspec+1) = [3., 1., 0., 0., 2., 1., 1., 1.]                         ! Dissociation factors
  REAL, TARGET, PRIVATE               :: auxrhoic = 917.                                                               ! (Reference) density of ice                   
  REAL, TARGET, PRIVATE               :: auxrhosn = 300.                                                               ! (Reference) density of snow 
  ! -------------------------------------------------------------------------------------------------------------------------

  !
  ! ----------------------------------------------------------------------
  ! Type Species holds a collection of names and properties for the currently 
  ! available aerosol species (including water). Densities, molar masses etc.
  ! in mo_submctl are therefore deprecated. 
  !
  TYPE Species

     ! The below arrays are compiled according to the configuration given in the SALSA namelist. Check out the constructor method for details
     LOGICAL            :: used(maxspec+1) = [.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.FALSE.,.TRUE.]     ! Logical mask of active compounds
     INTEGER            :: Nused = 0                                                                              ! Number of active compounds
     INTEGER            :: allInd(maxspec+1) = [0,0,0,0,0,0,0,0]                                                  ! Indices of compounds in model mass arrays (0 if not active)

     ! The suffix "Native" for arrays following the order of allNames. Non-"Native" follow the order in which the active compounds are given in SALSA namelist and therefore in model mass arrays.
     ! The lists are trucated to include only the active compound properties using the logical mask "used" and the Fortran intrinsic PACK.
     CHARACTER(len=3), ALLOCATABLE :: namesNative(:), names(:)                                     ! Names of active compounds
     INTEGER         , ALLOCATABLE :: indNative(:), ind(:)                                         ! Indices in model mass arrays
     REAL            , ALLOCATABLE :: MMNative(:), MM(:)                                           ! Molar masses
     REAL            , ALLOCATABLE :: rholiqNative(:), rhoiceNative(:), rhosnowNative(:),   &      ! densities
                                      rholiq(:), rhoice(:), rhosnow(:)
     REAL            , ALLOCATABLE :: rhos(:,:)                                                    ! 2d density array holding all the variants (shape (3,Nused)) 
     REAL            , ALLOCATABLE :: dissNative(:), diss(:)                                       ! Dissociation factors

     ! Subset arrays for soluble and insoluble compounds. If none are used of either category, length 1 is allocated for the arrays and zeros are used to initialize. 
     ! Thus be carefull when using these, check with N(in)soluble!
     INTEGER                       :: Nsoluble = 0, Ninsoluble = 0                                 ! Number of active soluble and insolubl compounds            
     CHARACTER(len=3), ALLOCATABLE :: names_soluble(:), names_insoluble(:)                         ! Names of active soluble and insoluble compounds (Length 1 empty string if non used)
     INTEGER,          ALLOCATABLE :: ind_soluble(:), ind_insoluble(:)                             ! Indices of soluble/insoluble
     REAL,             ALLOCATABLE :: MM_soluble(:), MM_insoluble(:)                               ! Molar masses
     REAL,             ALLOCATABLE :: rholiq_soluble(:),   &                                       ! Densities
                                      rhoice_soluble(:),   &
                                      rhosnow_soluble(:),  &
                                      rholiq_insoluble(:), &
                                      rhoice_insoluble(:), &
                                      rhosnow_insoluble(:)
     REAL,             ALLOCATABLE :: diss_soluble(:), diss_insoluble(:)                           ! Dissociation factors

     REAL    :: mas = 132.14e-3  ! Molar mass of ammonium sulphate
     REAL    :: rhoas = 1770.    ! Density of ammonium sulphate

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
       PROCEDURE, PRIVATE :: defineSubsets

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
      ALLOCATE( cnstr%namesNative(cnstr%Nused),   cnstr%names(cnstr%Nused),     &
                cnstr%indNative(cnstr%Nused),     cnstr%ind(cnstr%Nused),       &
                cnstr%MMNative(cnstr%Nused),      cnstr%MM(cnstr%Nused),        &
                cnstr%rholiqNative(cnstr%Nused),  cnstr%rholiq(cnstr%Nused),    &
                cnstr%rhoiceNative(cnstr%Nused),  cnstr%rhoice(cnstr%Nused),    &
                cnstr%rhosnowNative(cnstr%Nused), cnstr%rhosnow(cnstr%Nused),   &
                cnstr%rhos(3,cnstr%Nused),                                       &
                cnstr%dissNative(cnstr%Nused),    cnstr%diss(cnstr%Nused)       )

      ! Truncate the property lists in the order given by the global field "allNames"
      ! in the beginning of this class.
      cnstr%namesNative = PACK(allNames, cnstr%used)
      cnstr%indNative = PACK(cnstr%allInd, cnstr%used)
      cnstr%MMNative = PACK(allMM, cnstr%used)
      cnstr%dissNative = PACK(allDiss, cnstr%used)
      cnstr%rholiqNative = PACK(allRho, cnstr%used)
      ! First, use the same arrays for densities in ice and snow arrays
      cnstr%rhoiceNative = cnstr%rholiqNative
      cnstr%rhosnowNative = cnstr%rholiqNative
      ! Second, replace the water densities with appropriate values
      cnstr%rhoiceNative(cnstr%Nused) = auxrhoic
      cnstr%rhosnowNative(cnstr%Nused) = auxrhosn
      
      ! Make another set of truncated property lists, where the values are sorted to the same order
      ! as the compounds appear in the UCLALES-SALSA mass arrays
      CALL cnstr%sortProperties()

      ! Build some subset arrays
      CALL cnstr%defineSubsets()

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
      SELF%rhos(:,:) = 0.
      SELF%names(:) = '   '
      SELF%diss(:) = 0.
      SELF%MM(:) = 0.
      DO ii = 1, SELF%Nused
         jj = SELF%indNative(ii)  ! This is the index the compound number ii has in the mass arrays
         SELF%rholiq(jj) = SELF%rholiqNative(ii)
         SELF%rhoice(jj) = SELF%rhoiceNative(ii)
         SELF%rhosnow(jj) = SELF%rhosnowNative(ii)
         SELF%names(jj) = SELF%namesNative(ii)
         SELF%diss(jj) = SELF%dissNative(ii)
         SELF%MM(jj) = SELF%MMNative(ii)
         SELF%ind(jj) = SELF%indNative(ii) ! This is of course trivial but can be used to check everything's ok
      END DO
      
      ! Collect all the density variants in a unified array
      SELF%rhos(1,:) = SELF%rholiq(:)
      SELF%rhos(2,:) = SELF%rhoice(:)
      SELF%rhos(3,:) = SELF%rhosnow(:)


    END SUBROUTINE sortProperties

    ! ---------------------------------------------------------------------------

    SUBROUTINE defineSubsets(SELF)
      IMPLICIT NONE
      CLASS(Species), INTENT(inout) :: SELF
      INTEGER :: i,n,ins,sol
      LOGICAL, ALLOCATABLE :: issol(:), isins(:)

      ALLOCATE(issol(SELF%Nused),isins(SELF%Nused))
      issol(:) = .FALSE.
      isins(:) = .FALSE.

      ! Number of active soluble and isoluble compounds
      n = 0
      DO i = 1,maxsol
         IF ( SELF%isUsed(allNamesSoluble(i)) ) n = n + 1
      END DO      
      SELF%Nsoluble = n

      n = 0
      DO i = 1,maxins
         IF ( SELF%isUsed(allNamesInsoluble(i)) ) n = n + 1
      END DO
      SELF%Ninsoluble = n

      ! Allocate subset arrays. If the determined size is 0, allocate at least 1 (and put zeros/empty strings)
      ! to avoid segmentation faults. Nsoluble/Ninsoluble should always be used to check if there is anyting active...
      ins = MIN(1,SELF%Ninsoluble)
      sol = MIN(1,SELF%Nsoluble)
      ALLOCATE( SELF%names_soluble(sol), SELF%ind_soluble(sol), SELF%MM_soluble(sol),           &
                SELF%rholiq_soluble(sol), SELF%rhoice_soluble(sol), SELF%rhosnow_soluble(sol),  &
                SELF%diss_soluble(sol)                                                          )
      ALLOCATE( SELF%names_insoluble(ins), SELF%ind_insoluble(ins), SELF%MM_insoluble(ins),           &
                SELF%rholiq_insoluble(ins), SELF%rhoice_insoluble(ins), SELF%rhosnow_insoluble(ins),  &
                SELF%diss_insoluble(ins)                                                              )

      SELF%names_soluble = 'x'; SELF%names_insoluble = 'x'
      SELF%ind_soluble = 0; SELF%ind_insoluble = 0
      SELF%MM_soluble = 0.; SELF%MM_insoluble = 0.
      SELF%rholiq_soluble = 0.; SELF%rholiq_insoluble = 0.
      SELF%rhoice_soluble = 0.; SELF%rhoice_insoluble = 0.
      SELF%rhosnow_soluble = 0.; SELF%rhosnow_insoluble = 0.
      SELF%diss_soluble = 0.; SELF%diss_insoluble = 0.

      IF ( SELF%Nsoluble > 0 ) THEN
         DO i = 1,maxsol
            WHERE ( SELF%names == allNamesSoluble(i) ) issol = .TRUE.
         END DO  
      END IF

      IF ( SELF%Ninsoluble > 0 ) THEN
         DO i = 1,maxins
            WHERE ( SELF%names == allNamesInsoluble(i) ) isins = .TRUE.
         END DO
      END IF
         
      SELF%names_soluble = PACK(SELF%names,issol)
      SELF%ind_soluble = PACK(SELF%ind,issol)
      SELF%MM_soluble = PACK(SELF%MM,issol)
      SELF%rholiq_soluble = PACK(SELF%rholiq,issol)
      SELF%rhoice_soluble = PACK(SELF%rhoice,issol)
      SELF%rhosnow_soluble = PACK(SELF%rhosnow,issol)
      SELF%diss_soluble = PACK(SELF%diss,issol)

      SELF%names_insoluble = PACK(SELF%names,isins)
      SELF%ind_insoluble = PACK(SELF%ind,isins)
      SELF%MM_insoluble = PACK(SELF%MM,isins)
      SELF%rholiq_insoluble = PACK(SELF%rholiq,isins)
      SELF%rhoice_insoluble = PACK(SELF%rhoice,isins)
      SELF%rhosnow_insoluble = PACK(SELF%rhosnow,isins)
      SELF%diss_insoluble = PACK(SELF%diss,isins)

      DEALLOCATE(issol,isins)

    END SUBROUTINE defineSubsets

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

      IF (PRESENT(type)) THEN
         IF (type == 'wet') THEN
            getNSpec = SELF%Nused
         ELSE IF (type == 'dry') THEN
            getNSpec = SELF%Nused - 1
         ELSE
            STOP "getNSpec: Bad type option"
         END IF
      ELSE
         getNSpec = SELF%Nused
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
      getRhoByIndex = 0.
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
