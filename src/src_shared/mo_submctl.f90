MODULE mo_submctl
  USE classSpecies, ONLY : Species, maxspec
  USE classProcessSwitch, ONLY : ProcessSwitch
  USE classBinLayout, ONLY : BinLayout
  IMPLICIT NONE

  SAVE

  ! I'd say nothing here needs to be PRIVATE so removed explicit PRIVATE and PUBLIC attributes (PUBLIC is default).
  ! -Juha
  
  TYPE t_parallelbin
     ! Map bin indices between parallel size distributions
     INTEGER :: cur  ! Index for current distribution
     INTEGER :: par  ! Index for corresponding parallel distribution
  END TYPE t_parallelbin

  !
  ! Master switches for SALSA aerosol microphysical processes
  ! --------------------------------------------------------------
  INTEGER, PARAMETER :: Nmaster = 7
  TYPE(ProcessSwitch), TARGET :: lsmaster(Nmaster)     ! Array for master switches. The specific master switches are pointers to this array
  TYPE(ProcessSwitch), POINTER :: lscoag => NULL()     ! Coagulation,      mode = 1: calculate kernels every timestep,
                                                       !                   mode = 2: use reduced update freq
  TYPE(ProcessSwitch), POINTER :: lscnd => NULL()      ! Condensation
  TYPE(ProcessSwitch), POINTER :: lsauto => NULL()     ! Autoconversion,   mode = 1: parameterized simple autoconversion,
                                                       !                   mode = 2: coagulation based precip formation
  TYPE(ProcessSwitch), POINTER :: lsactiv => NULL()    ! Cloud activation, mode = 1: aerosol growth based activation,
                                                       !                   mode = 2: parameterized cloud base activation
  TYPE(ProcessSwitch), POINTER :: lsicenucl => NULL()  ! Ice nucleation,   mode = 1: With %switch=TRUE, %state=FALSE until %delay time is
                                                       !                             reached and ice nucleation is not called, afterwards
                                                       !                             %state=TRUE and ice nucleation is called as usual.
                                                       !                   mode = 2: With %switch=TRUE, %state=FALSE until %delay time is reached,
                                                       !                             but ice nucleation is called anyway. However, it will only
                                                       !                             remove aerosol and not produce ice. After %delay, %state=TRUE
                                                       !                             and ice nucleation will continue normal operation.
  TYPE(ProcessSwitch), POINTER :: lsicemelt => NULL()  ! Melting of ice
  TYPE(ProcessSwitch), POINTER :: lssecice => NULL()


  !
  ! Subprocess switches for microphysical processes; Note that master switch set to FALSE overrides any setting
  ! in the subprocess switches.
  ! ---------------------------------------------------------------------------------------------------------------
  
  ! Coagulation/collision subprocesses
  LOGICAL :: lscgaa  = .TRUE.  ! Coagulation between aerosols
  LOGICAL :: lscgcc  = .TRUE.  ! Collision-coalescence between cloud droplets
  LOGICAL :: lscgca  = .TRUE.  ! Cloud collection of aerosols
  LOGICAL :: lscgpc  = .TRUE.  ! Collection of cloud droplets by rain
  LOGICAL :: lscgpa  = .TRUE.  ! Collection of aerosols by rain
  LOGICAL :: lscgpp  = .TRUE.  ! Collision between rain drops
  LOGICAL :: lscgia  = .TRUE.  ! Ice collection of aerosols
  LOGICAL :: lscgic  = .TRUE.  ! Collection of cloud droplet by ice particles
  LOGICAL :: lscgii  = .TRUE.  ! Collision-coalescence between ice particles
  LOGICAL :: lscgip  = .TRUE.  ! Collection of precipitation by ice particles

  ! Condensation subprocesses
  LOGICAL :: lscndgas   = .FALSE. ! Condensation of precursor gases
  LOGICAL :: lscndh2ocl = .TRUE.  ! Condensation of water vapour on clouds and precipitation
  LOGICAL :: lscndh2oae = .TRUE.  ! Condensation of water vapour on aerosol particles (FALSE -> equilibrium calc.)
  LOGICAL :: lscndh2oic = .TRUE.  ! Condensation of water vapour on ice and snow

  ! Ice nucleation subprocesses
  LOGICAL :: lsicehom = .FALSE.        ! Homogeneous freezing
  LOGICAL :: lsiceimm = .FALSE.        ! Immersion freezing
  LOGICAL :: lsicedep = .FALSE.        ! Deposition freezing

  ! Secondary ice subprocesses
  LOGICAL :: lssiprimespln = .FALSE.                       ! Rime splintering; Hallet-Mossop
  TYPE(ProcessSwitch), POINTER :: lssipdropfrac => NULL()  ! Drop fracturing SIP, %mode = 1: Lawson et al 2015,
                                                           !                      %mode = 2: Phillips et al 2018

  INTEGER, PARAMETER :: Nsub = 1
  TYPE(ProcessSwitch), TARGET :: lssub(Nsub)  ! Holder for type ProcessSwitch subprocess switches (for most just the simple
                                              ! logical switch is required)

  
  ! ---------------------------------------------------------------------------------------------------------------
  ! Reduced coagulation kernel update frequency
  REAL :: cgintvl = 10.             ! If lscoag%mode = 2, gives the coagulation kernel update interval in seconds.
  
  LOGICAL :: lcgupdt = .FALSE.      ! Switch for updating kernels with lscoag%mode = 2; Note: This is determined
                                    ! during runtime -> NOT a NAMELIST parameter!
  ! Aggregation efficiency range for ice. Assumed to be inversely proportional to rime fraction
  ! within this range
  REAL :: Eiagg_max = 0.3
  REAL :: Eiagg_min = 0.05

  ! Contact angle distribution for ice nucleation:
  ! Use contact angle distributions for heterogeneous nucleation
  ! processes according to Savre and Ekman (2015).
  ! initMinTheta is used to specify a minimum contact angle for the entire IN population
  ! during initial stages of the simulation. The period when initMinTheta is applied
  ! is from model initialization to ice_theta_dist%delay (in seconds)
  LOGICAL :: ice_theta_dist = .TRUE.   ! Use contact angle integration in ice nucleation
  TYPE(ProcessSwitch) :: lsFreeTheta   ! Switch: Use freely evolving lower limit in contact angle integration
                                       ! delay: Time until which the low limit contact for contact angle integration
                                       !        is fixed as initMinTheta.
  REAL :: initMinTheta = 0.            ! Fixed lower limit for contact angle integration

  ! Contact angle distribution parameters for immersion freezing
  REAL :: mean_theta_imm = 132. ! mean contact angle
  REAL :: sigma_theta_imm = 20.  ! STD of the contact angle distribution.
  ! Contact angle distribution parameters for deposition freezing
  REAL :: mean_theta_dep = 15.5
  REAL :: sigma_theta_dep = 1.4
  
  LOGICAL :: lsdistupdate = .TRUE.    ! Perform the size distribution update
  LOGICAL :: lscheckarrays = .FALSE.  ! Do some primitive error checking in the SALSA main program

  ! RH Limit:
  ! If lsfreeRH=FALSE, use RH constrained by *rhlim* for SALSA processes. Otherwise use the predicted value.
  ! If lsfreeTH=TRUE and lsfreeRH%delay > 0, the constrain is active until that time, and predicted rh is used after the delay time.
  TYPE(ProcessSwitch) :: lsfreeRH   
  REAL :: rhlim = 1.20

  ! 1) Switches for aerosol microphysical processes ------------------------
  INTEGER, PARAMETER :: nmod = 7

  INTEGER :: nsnucl = 0             ! Choice of the nucleation scheme:
                                    ! 0 = off
                                    ! 1 = binary nucleation
                                    ! 2 = activation type nucleation
                                    ! 3 = kinetic nucleation
                                    ! 4 = ternary nucleation
                                    ! 5 = nucleation with ORGANICs
                                    ! 6 = activation type of nucleation with H2SO4+ORG
                                    ! 7 = heteromolecular nucleation with H2SO4*ORG
                                    ! 8 = homomolecular nucleation of  H2SO4 +
                                    !           heteromolecular nucleation with H2SO4*ORG
                                    ! 9 = homomolecular nucleation of  H2SO4 and ORG +
                                    !           heteromolecular nucleation with H2SO4*ORG
  
  ! 1) Switches for SALSA aerosol microphysical processes ------------------------

  INTEGER ::                    & ! J3 parametrization
      nj3 = 1                     ! 1 = condensational sink (Kerminen&Kulmala, 2002)
                                  ! 2 = coagulational sink (Lehtinen et al. 2007)
                                  ! 3 = coagS+self-coagulation (Anttila et al. 2010)
  REAL :: act_coeff = 1.e-7  ! activation coefficient

  ! SALSA size distribution definitions
  ! ================================================================================================
  ! Define which aerosol species used and initial size distributions
  TYPE(Species), TARGET :: spec  ! Must be initialized in mo_salsa_init (pointer associations). Holds aerosol species indices and properties
  INTEGER :: nspec_dry = 1
  CHARACTER(len=3) :: listspec(maxspec) = (/'SO4','   ','   ','   ','   ','   ','   '/)
  
  ! Volume fractions between aerosol species for A and B-bins
  REAL :: volDistA(maxspec) = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  REAL :: volDistB(maxspec) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)

  INTEGER :: isdtyp = 0  ! Type of input aerosol size distribution: 0 - Uniform
                         !                                          1 - Read vertical profile of the mode
                         !                                              number concentrations from an input file

  ! Geometric standard deviation, mode mean diameter (um) and mode concentration in #/mg (~ #/cm3)
  ! Separate for A and B regimes
  REAL :: sigmagA(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),               & 
          dpgA(nmod) = (/0.03, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2/),          & 
          nA(nmod) = (/1600.,640.,0.,0.,0.,0.,0./)                  
  REAL :: sigmagB(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),                &
          dpgB(nmod) = (/0.03,0.15,0.2,0.2,0.2,0.2,0.2/),                &
          nB(nmod) = (/0.,0.,0.,0.,0.,0.,0./)
  
  ! ==================================================================================================
  
  
  ! Options for ice nucleation (when master switch nlicenucl = .TRUE,)
  ! a) Constant ice number concentration (fixinc > 0 #/kg) is maintained by converting cloud droplets to ice
  REAL :: fixinc = -1.0 ! Default = disabled
  ! b) Modelled ice nucleation || Juha: These could be placed with the other switches. Also, would 
  !                               do the fixed thing with the condition nlicenucl = false and fixinc > 0.,
  !                               otherwise physical nucleation.


  INTEGER, PARAMETER ::            &
       nreg = 2                          ! number of main size regimes

  REAL ::                       &
       reglim(nreg+2) =                            & ! low/high diameter limits
       (/ 3.e-9, 4.e-8, 7.e-7, 1.e-5 /) ! of main size regimes [m]

  INTEGER :: &
       nbin(nreg) = (/ 3, 7 /)   ! number of bins in each main regime
    
  INTEGER ::      & ! start index for bin regimes
       in1a,          & ! regime 1a
       in2a,          & ! regime 2a
       in2b,          & ! regime 2b
       ! last index for bin regimes
       fn1a,          & ! regime 1a
       fn2a,          & ! regime 2a
       fn2b,          & ! regime 2b
       nbins            ! total number of size bins

  ! Juha: Cloud and rain bins:
  TYPE(t_parallelbin) ::   ica, & ! cloud droplets (first, regime a)
                           fca, & ! cloud droplets (last, regime a)
                           icb, & ! cloud droplets (first, regime b)
                           fcb    ! cloud droplets (last, regime b)
  INTEGER             ::   ira,fra! Rain/drizzle bin indices
  INTEGER             ::   ncld   ! Total number of cloud bins
  INTEGER             ::   nprc

  TYPE(BinLayout)     :: bloPrc   ! Precipitation bin definitions
  TYPE(BinLayout)     :: bloIce   ! Ice bin definitions
  
  ! Jaakko: ice bins:
  INTEGER ::   iia, & ! ice particles, first
               fia ! ice particles, last
  INTEGER ::   nice   ! Total number of ice bins
    
  INTEGER :: ntotal  ! Total number of bins across all active particle types
  INTEGER :: nliquid ! Total number of bins across liquid particle types
  INTEGER :: nfrozen ! Total number of bins across frozen particle types

  REAL, ALLOCATABLE :: aerobins(:),  &  ! These are just to deliver information about the bin diameters if the
                       cloudbins(:), &  ! host model needs it (lower limits).
                       precpbins(:), &
                       icebins(:)


  ! Diameter limits for particles in meters. These should be used when creating
  ! classSection instances
  REAL, PARAMETER :: dlaero = 30.e-6,   &
                     dlcloud = 100.e-6, &
                     dlprecp = 5.e-3,   &
                     dlice   = 10.e-3
    
  REAL, PARAMETER ::     &
   avog   = 6.0221e+23,   & ! Avogadro number (#/mol)
   boltz  = 1.3807e-23,   & ! Boltzmann constant (J/K)
   planck = 6.626070040e-34, & ! Planck constant (J*s)
   grav   = 9.81,         & ! gravitational acceleration (m/s^2)
   pstand = 1.01325e+5,   & ! standard pressure (Pa)
   rg     = 8.314,        & ! molar gas constant (J/(mol K))
   pi     = 3.1415927,    & ! self explanatory
   pi6    = 0.5235988,    & ! pi/6
   cpa    = 1005.,        & ! specific heat of dry air, constant P (J/kg/K)
   cwa    = 4200.,        & ! Specific heat capacity of lqiuid water (J/kg/K)
   mair   = 28.97e-3,     & ! molar mass of air (mol/kg)
   deltav = 1.096e-7,     & ! vapor jump length (m)
   deltaT = 2.16e-7,      & ! thermal jump length (m)
   alphaT = 0.96,         & ! thermal accomodation coefficient
   alphac = 1.0,          & ! condensation coefficient
   eps    = epsilon(1.0)       ! epsilon

  REAL, PARAMETER ::   &
   rd    = 287.04,     & ! gas constant for dry air (J/K/kg)
   rv    = 461.5,      & ! gas constant for water vapour (J/K/kg)
   alv    = 2.5e6,     & ! latent heat of vaporisation (J/kg)
   als    = 2.834e6,   & ! latent heat of sublimation (J/kg)
   alf    = 3.3e5        ! Latent heat of freezing (J/kg)

  REAL, PARAMETER ::     & ! molar mass [kg/mol]
       msu = 98.08e-3,        & ! sulphate
       mno = 62.01e-3,        & ! HNO3
       mnh = 18.04e-3,        & ! NH3
       moc = 150.e-3,         & ! organic carbon
       mbc = 12.e-3,          & ! black carbon
       mss = 58.44e-3,        & ! sea salt (NaCl)
       mdu = 100.e-3,         & ! mineral dust
       mwa = 18.016e-3,       & ! water
       mas = 132.14e-3,       & ! ammoniums sulphate ((NH4)2SO4)
       !
       ! densities [kg/m3]
       rhosu = 1830.,         & ! sulphate
       rhono = 1479.,         & ! HNO3
       rhonh = 1530.,         & ! NH3
       rhooc = 2000.,         & ! organic carbon
       rhobc = 2000.,         & ! black carbon
       rhoss = 2165.,         & ! sea salt (NaCl)
       rhodu = 2650.,         & ! mineral dust
       rhowa = 1000.,         & ! water
       rhoic = 917.,          & ! ice
       !
       ! volume of molecule [kg/#]
       mvsu = msu/avog/rhosu,    & ! sulphate
       mvno = mno/avog/rhono,    & ! HNO3
       mvnh = mnh/avog/rhonh,    & ! NH3
       mvoc = moc/avog/rhooc,    & ! organic carbon
       mvss = mss/avog/rhoss,    & ! sea salt
       mvwa = mwa/avog/rhowa,    &
       !
       volratio =                & ! ratio of molecular volumes for
       (msu*rhoss)/(rhosu*mss), & ! sulphate and sea salt
       !
       n3 = 158.79               ! number of H2SO4 molecules in 3 nm cluster
  !  assuming d_sa = 5.54 ???
  
  !-- 4.3) Properties of condensing vapours
  
  REAL, PARAMETER :: & ! diameter of condensing molecule [m]
       d_sa  = 5.539376964394570e-10,  &
       d_oc  = 6.195906936656752e-10,  &
       d_h2o = 3.851565216195334e-10
  
  REAL, PARAMETER :: &
       ions   = 3.0,       & ! van't Hoff factor (ions produced upon dissociation)
       surfw0 = 0.073,     & ! surface tension of pure water @ ~ 293 K [J/m2]
       surfi0 = 0.105,     & ! surface tension of ice
       epsoc  = 0.15         ! water uptake of organic material
  
  !-- 7) Parameters for cloud activation
  
  REAL, ALLOCATABLE :: massacc(:)
  
  REAL, PARAMETER :: &
   nlim = 1.,  & ! Number conc. limit (#/kg) for aerosol and cloud droplets 
   prlim = 1.e-6 ! The same for precipitation and ice species for which concentrations are normally much lower [#/kg]
  
END MODULE mo_submctl
