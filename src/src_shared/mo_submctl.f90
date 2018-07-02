MODULE mo_submctl
  USE classSpecies, ONLY : Species, maxspec
  USE classProcessSwitch, ONLY : ProcessSwitch
  USE classBinLayout, ONLY : BinLayout
  IMPLICIT NONE

  ! I'd say nothing here needs to be PRIVATE so removed explicit PRIVATE and PUBLIC attributes (PUBLIC is default).
  ! -Juha
  
  TYPE t_parallelbin
     ! Map bin indices between parallel size distributions
     INTEGER :: cur  ! Index for current distribution
     INTEGER :: par  ! Index for corresponding parallel distribution
  END TYPE t_parallelbin
       
  !Switches for SALSA aerosol microphysical processes

  INTEGER, PARAMETER :: Nmaster = 7
  TYPE(ProcessSwitch), TARGET :: lsmaster(Nmaster)  ! Array for master switches. The specific master switches are pointers to this array
  TYPE(ProcessSwitch), POINTER :: lscoag => NULL()     ! Coagulation
  TYPE(ProcessSwitch), POINTER :: lscnd => NULL()      ! Condensation
  TYPE(ProcessSwitch), POINTER :: lsauto => NULL()     ! Autoconversion, mode = 1: parameterized simple autoconversion, mode = 2: coagulation based precip formation
  TYPE(ProcessSwitch), POINTER :: lsautosnow => NULL() ! Autoconversion of snow (need to revise)
  TYPE(ProcessSwitch), POINTER :: lsactiv => NULL()    ! Cloud activation, mode = 1: aerosol growth based activation, mode = 2: parameterized cloud base activation
  TYPE(ProcessSwitch), POINTER :: lsicenucl => NULL()  ! Ice nucleation
  TYPE(ProcessSwitch), POINTER :: lsicemelt => NULL()  ! Melting of ice

  ! Collision subprocesses
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
  LOGICAL :: lscgsa  = .TRUE.  ! Collection of aerosols by snow              
  LOGICAL :: lscgsc  = .TRUE.  ! Collection of cloud droplets by snow
  LOGICAL :: lscgsi  = .TRUE.  ! Collection of ice by snow
  LOGICAL :: lscgsp  = .TRUE.  ! Collection of precipitation by snow
  LOGICAL :: lscgss  = .TRUE.  ! Collision-coalescence between snow particles

  ! Condensation subprocesses
  LOGICAL :: lscndgas   = .FALSE. ! Condensation of precursor gases
  LOGICAL :: lscndh2ocl = .TRUE.  ! Condensation of water vapour on clouds and precipitation
  LOGICAL :: lscndh2oae = .TRUE.  ! Condensation of water vapour on aerosol particles (FALSE -> equilibrium calc.)
  LOGICAL :: lscndh2oic = .TRUE.  ! Condensation of water vapour on ice and snow

  LOGICAL :: lsdistupdate = .TRUE.  ! Perform the size distribution update

  LOGICAL :: lscheckarrays = .FALSE.

  TYPE(ProcessSwitch) :: lsfreeRH   ! If FALSE, use RH constrained by *rhlim* for SALSA processes. Otherwise use the predicted value.
                                    ! If lsfreeRH%delay > 0, the constrain is active until that time.
  ! RH Limit: used for initialization and spinup within SALSA to limit the water vapour mixing ratio.
  ! Prevents unrealistically high RH in cloud activation and condensation procedures that is often assigned
  ! in the LES input files to immediately generate cloud. Given in %/100.
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
  
  ! Define which aerosol species used and initial size distributions
  TYPE(Species), TARGET :: spec  ! Must be initialized in mo_salsa_init (pointer associations). Holds aerosol species indices and properties
  INTEGER :: nspec_dry = 1
  CHARACTER(len=3) :: listspec(maxspec) = (/'SO4','   ','   ','   ','   ','   ','   '/)
  
  ! Volume fractions between aerosol species for A and B-bins
  REAL :: volDistA(maxspec) = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  REAL :: volDistB(maxspec) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  ! Number fraction allocated to a-bins in regime 2 (b-bins will get 1-nf2a)
  REAL :: nf2a = 1.0

  ! Options for ice nucleation (when master switch nlicenucl = .TRUE,)
  ! a) Constant ice number concentration (fixinc > 0 #/kg) is maintained by converting cloud droplets to ice
  REAL :: fixinc = -1.0 ! Default = disabled
  ! b) Modelled ice nucleation || Juha: These could be placed with the other switches. Also, would 
  !                               do the fixed thing with the condition nlicenucl = false and fixinc > 0.,
  !                               otherwise physical nucleation.
  LOGICAL :: ice_hom = .FALSE., ice_imm=.FALSE., ice_dep=.FALSE. ! Available ice nucleation modes


  INTEGER :: isdtyp = 0  ! Type of input aerosol size distribution: 0 - Uniform
                         !                                          1 - Read vertical profile of the mode
                         !                                              number concentrations from an input file

  REAL :: sigmag(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),               & ! Stdev
          dpg(nmod) = (/0.03, 0.15, 0.2, 0.2, 0.2, 0.2, 0.2/),    & ! Mode diam in um
          n(nmod) = (/1600.,640.,0.,0.,0.,0.,0./)                   ! #/mg ~ #/cm3

  INTEGER, PARAMETER ::            &
       nreg = 2                          ! number of main size regimes

  REAL ::                       &
       reglim(nreg+2) =                            & ! low/high diameter limits
       (/ 3.e-9, 5.e-8, 7.e-7, 1.e-5 /) ! of main size regimes [m]

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
  
  ! Jaakko: ice bins:
  TYPE(t_parallelbin) ::   iia, & ! ice particles (first, regime a)
                           fia, & ! ice particles (last, regime a)
                           iib, & ! ice particles (first, regime b)
                           fib    ! ice particles (last, regime b)
  INTEGER             ::   isa,fsa! snow bin indices
  INTEGER             ::   nice   ! Total number of ice bins
  INTEGER             ::   nsnw   ! Total number of snow bins
    
  INTEGER :: ntotal ! Total number of bins accross all active particle types

  REAL, ALLOCATABLE :: aerobins(:),  &  ! These are just to deliver information about the bin diameters if the
                       cloudbins(:), &  ! host model needs it (lower limits).
                       precpbins(:), &
                       icebins(:), &
                       snowbins(:)

  ! Diameter limits for particles in meters. These should be used when creating
  ! classSection instances
  REAL, PARAMETER :: dlaero = 30.e-6,   &
                     dlcloud = 100.e-6, &
                     dlprecp = 2.e-3,   &
                     dlice   = 2.e-3,   &
                     dlsnow  = 10.e-3

    
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
   mair   = 28.97e-3,     & ! molar mass of air (mol/kg)
   deltav = 1.096e-7,     & ! vapor jump length (m)
   deltaT = 2.16e-7,      & ! thermal jump length (m)
   alphaT = 0.96,         & ! thermal accomodation coefficient
   alphac = 1.0,          & ! condensation coefficient
   eps    = epsilon(1.0)       ! epsilon

  REAL, PARAMETER ::   &
   rd    = 287.04,     & ! gas constant for dry air (J/K/kg)
   rv    = 461.5,      & ! gas constant for water vapour (J/K/kg)
   alv    = 2.5e6,   & ! latent heat for vaporisation (J/kg)
   als    = 2.834e6      ! latent heat for sublimation (J/kg)

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
       rhosn = 300.,          & ! snow
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
   prlim = 1.e-6 ! The same for precipitation and ice species for which concentrations are normally much lower [#/m3]
  
END MODULE mo_submctl
