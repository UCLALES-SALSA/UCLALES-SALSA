MODULE emission_types
  IMPLICIT NONE

  !
  ! ----------------------------------------------
  ! Contains arrays for emitted size distribution
  !
  TYPE EmitSizeDist
     REAL, ALLOCATABLE :: numc(:)
     REAL, ALLOCATABLE :: mass(:)
  END TYPE EmitSizeDist

  !
  ! --------------------------------------------
  ! Type for emission configuration options
  !
  TYPE  EmitConfig
     INTEGER          :: emitType = 1                 ! 1: Natural seasalt emissions, 2: custom artificial emissions
                                                      ! Ali, addition of emission type 3
                                                      ! 3: artificial emission given by a map (moving source of airborne emission)
                                                      ! Ali, addition of emission type 3
     INTEGER          :: regime = 1                   ! Destination bin regime for emitted aerosol. 1: A, 2: B
     REAL             :: start_time = 0.,  &          ! Start time for emission (s)
                         end_time = 86400.            ! End time for emission (s)

     ! Parameters below valid for emitType > 1
     CHARACTER(len=3) :: species = 'SS '              ! Which aerosol species to emit (must conform to names in src_salsa/classSpecies)
     REAL             :: emitHeightMin = -999.,  &    ! Min height of airborne emissions (m)
                         emitHeightMax = -999.        ! Max height (m)
     INTEGER          :: emitLevMin = -999            ! Integer levels corresponding to the heights; If both heights and lev indices are specified,
     INTEGER          :: emitLevMax = -999            ! the level indices are preferred (so use one or the other...).
     REAL             :: emitLatmin = -1.e6           ! Min and max X/"lon" offsets in meters for the emission area
     REAL             :: emitLatmax = 1.e6
     REAL             :: emitLonmin = -1.e6           ! Min and max Y/"lat" offsets in meters for the emission area
     REAL             :: emitLonmax = 1.e6

     INTEGER          :: emitSizeDistType = 1         ! 1: Monochromatic aerosol, 2: modal size disribution (lognormal)
     REAL             :: emitDiam = 10.e-6,    &      ! Assumed (dry )diameter of the particles (mode diameter for emitType=2).
                         emitNum  = 10000.            ! Number consentration of particles emitted per second #/m3/s (mode concentration for emitType=2)
     REAL             :: emitSigma = 2.0              ! Geometric standard deviation for emitSizeDist=2
     ! Ali, addition of emission type 3
     CHARACTER(len=40):: emitMap = ''                 ! Name of the file providing all location of emission (only for emitType = 3)
     REAL             :: scS = 60.                    ! Source speed (m/s) (only for emitType = 3)
     INTEGER          :: z_expan_up = 0               ! Epands the emission map to adjacent cells above the given map
     INTEGER          :: z_expan_dw = 0               ! Epands the emission map to adjacent cells down the given map      
     ! /
  END TYPE EmitConfig

  ! Ali, addition of emission type 3; this basically extends EmitConfig?
   TYPE EmitType3Config
     INTEGER, ALLOCATABLE :: ix(:)  ! x Index of source location, calculated based on 'emitMap' and computational grid
     INTEGER, ALLOCATABLE :: iy(:)  ! y Index of source location, calculated based on 'emitMap' and computational grid
     INTEGER, ALLOCATABLE :: iz(:)  ! z Index of source location, calculated based on 'emitMap' and computational grid
     REAL, ALLOCATABLE :: t_in(:)   ! Enterance times of emission source in each subdomain/processor
     REAL, ALLOCATABLE :: t_out(:)  ! Exit times of emission source from each subdomain/processor
     REAL, ALLOCATABLE :: t(:)      ! Enterance times of emission source in each cell in each subdomain/processor
     REAL :: t_trac                 ! Tracks the time for the emission source in each subdomain/processor
     INTEGER :: np                  ! Number of points of emission trajectory intersecting with cell boundaries for each subdomain/processor
   END type EmitType3Config

   INTEGER, PARAMETER :: maxEmissionModes = 5                ! Max number of emission modes

   ! NAMELIST variables
   LOGICAL :: emitPristineIN = .TRUE.                        ! TRUE: when aerosol emissions active, IN active particles 
                                                             !       improve the IN efficiency of the target population.
			                                     ! FALSE: IN efficiency of target population remains intact.    
   INTEGER :: nEmissionModes = 0                             ! Number of emission modes (max == maxemissionModes)
   TYPE(EmitConfig), TARGET :: emitModes(MaxEmissionModes)   ! Configuration instances for emission modes
   TYPE(EmitSizeDist), TARGET :: emitData(MaxEmissionModes)  ! Emission size distribution data for each setup/mode.
   TYPE(EmitType3Config), TARGET :: emitType3(MaxEmissionModes)   ! Configuration instances for emission type 3
   
END MODULE emission_types
