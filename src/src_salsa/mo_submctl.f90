MODULE mo_submctl
  USE classSpecies, ONLY : Species, maxspec
  USE classSection
  USE classProcessSwitch, ONLY : ProcessSwitch
  IMPLICIT NONE

  ! I'd say nothing here needs to be PRIVATE so remoed PRIVATE and PUBLIC attributes (PUBLIC is default).
  ! -Juha
  
  ! Datatype used to store information about the binned size distributions of aerosols,cloud,drizzle and ice
  ! ---------------------------------------------------------------------------------------------------------
  
  TYPE t_parallelbin
     ! Map bin indices between parallel size distributions
     INTEGER :: cur  ! Index for current distribution
     INTEGER :: par  ! Index for corresponding parallel distribution
  END TYPE t_parallelbin
  

  ! Particle type specific pointers to "allSALSA" master array defined in mo_salsa_driver.
  ! Pointer association is done in mo_salsa_init. These should be accessed by importin mo_submctl,
  ! not by dummy arguments.
  TYPE(Section), POINTER :: aero(:,:,:)  => NULL(),   &
                            cloud(:,:,:) => NULL(),  &
                            precp(:,:,:) => NULL(),  &
                            ice(:,:,:)   => NULL(),    &
                            snow(:,:,:)  => NULL(),   &
                            liquid(:,:,:) => NULL(),  &
                            frozen(:,:,:) => NULL()
  
  !Switches for SALSA aerosol microphysical processes

  INTEGER, PARAMETER :: Nmaster = 7
  TYPE(ProcessSwitch), TARGET :: lsmaster(Nmaster)  ! Array for master switches. The specific master switches are pointers to this array
  TYPE(ProcessSwitch), POINTER :: lscoag => NULL()     ! Coagulation
  TYPE(ProcessSwitch), POINTER :: lscnd => NULL()      ! Condensation
  TYPE(ProcessSwitch), POINTER :: lsauto => NULL()     ! Autoconversion
  TYPE(ProcessSwitch), POINTER :: lsautosnow => NULL() ! Autoconversion of snow (need to revise)
  TYPE(ProcessSwitch), POINTER :: lsactiv => NULL()    ! Cloud activation
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

  ! Activation schemes
  LOGICAL :: lsactintst  = .TRUE.   ! Switch for interstitial activation: use particle wet size determined by
  LOGICAL :: lsactbase   = .FALSE.  ! Switch for cloud base activation: use the regular parameterized method

  LOGICAL :: lsdistupdate = .TRUE.  ! Perform the size distribution update

  LOGICAL :: lscheckarrays = .FALSE.

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

  ! RH Limit: used for initialization and spinup within SALSA to limit the water vapour mixing ratio.
  ! Prevents unrealistically high RH in cloud activation and condensation procedures that is often assigned
  ! in the LES input files to immediately generate cloud. Given in %/100.
  REAL :: rhlim = 1.20
  
  ! Define which aerosol species used and initial size distributions
  TYPE(Species), TARGET :: spec  ! Must be initialized in mo_salsa_init (pointer associations). Holds aerosol species indices and properties
  INTEGER :: nspec = 1
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
  INTEGER             ::   nprc   ! Total number of precipitation bins
  
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


CONTAINS
  !********************************************************************
  ! Function for calculating terminal velocities for different particle types and size ranges.
  !     Tomi Raatikainen (2.5.2017)
  REAL FUNCTION terminal_vel(radius,rhop,rhoa,visc,beta,flag)
    IMPLICIT NONE
    REAL, INTENT(in) :: radius, rhop ! Particle radius and density
    REAL, INTENT(in) :: rhoa, visc, beta ! Air density, viscocity and Cunningham correction factor
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    ! Constants
    REAL, PARAMETER :: rhoa_ref = 1.225 ! reference air density (kg/m^3)
    
    IF (flag==4) THEN   ! Ice
       ! Ice crystal terminal fall speed from Ovchinnikov et al. (2014)
       !       Dimension D = 2*radius
       terminal_vel = 12.0*sqrt(2.0*radius)
    ELSE IF (flag==5) THEN   ! Snow
       ! The same for snow
       !       Dimension D = 2*radius
       terminal_vel = 12.0*sqrt(2.0*radius)
    ELSE
       ! Aerosol and cloud and rain droplets
       IF (radius<40.0e-6) THEN
          ! Stokes law with Cunningham slip correction factor
          terminal_vel = (4.*radius**2)*(rhop-rhoa)*grav*beta/(18.*visc) ![m s-1]
       ELSE IF (radius<0.6e-3) THEN
          ! Droplets from 40 um to 0.6 mm: linear dependence on particle radius and a correction for reduced pressure
          !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
          terminal_vel = 8.e3*radius*sqrt(rhoa_ref/rhoa)
       ELSE
          ! Droplets larger than 0.6 mm: square root dependence on particle radius and a correction for reduced pressure
          !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
          ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
          terminal_vel = 2.01e2*sqrt( min(radius,2.0e-3)*rhoa_ref/rhoa)
       END IF
    END IF
  END FUNCTION terminal_vel
  
  !********************************************************************
  ! Function for calculating dimension (or wet diameter) for any particle type
  ! - Aerosol, cloud and rain are spherical
  ! - Snow and ice can be irregular and their densities can be size-dependent
  !
  ! Edit this function when needed also note calc_eff_radius in grid.f90
  !
  ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
  ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
  !
  SUBROUTINE CalcDimension(n,ppart,lim,dia,flag)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    TYPE(Section), INTENT(in) :: ppart(n)
    REAL, INTENT(IN) :: lim
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    REAL, INTENT(OUT) :: dia(n)
    INTEGER i
    
    dia(:) = 2.e-10
    DO i=1,n
       IF (ppart(i)%numc>lim) &
            dia(i)=(SUM(ppart(i)%volc(:))/ppart(i)%numc/pi6)**(1./3.)
    END DO
    
  END SUBROUTINE CalcDimension
  
  !********************************************************************
  ! Function for calculating equilibrium water saturation ratio at droplet surface based on Köhler theory
  !
  REAL FUNCTION calc_Sw_eq(part,T)
    TYPE(Section), INTENT(in) :: part ! Any particle
    REAL, INTENT(IN) :: T ! Absolute temperature (K)
    REAL :: dwet
    
    REAL :: znw,zns ! Moles of water and soluble material
    REAL :: zvw, zvs, zvtot ! Volume concentrations of water and soluble material and total dry
    INTEGER :: iwa, ndry ! Index for water, number of "dry" species
    INTEGER :: i
    
    iwa = spec%getIndex("H2O")
    ndry = spec%getNSpec(type="dry")
    
    ! Wet diameter
    dwet = (SUM(part%volc(:))/part%numc/pi6)**(1./3.)
    
    ! Equilibrium saturation ratio = xw*exp(4*sigma*v_w/(R*T*Dwet))
    
    znw = part%volc(iwa)*spec%rhowa/spec%mwa
    zvw = part%volc(iwa)*spec%rhowa/spec%mwa
    zns = 0.
    zvs = 0.
    zvtot = 0.
    DO i = 1,ndry
       zns = zns + spec%diss(i)*part%volc(i)*spec%rholiq(i)/spec%MM(i)
       zvs = zvs + MIN(1.,spec%diss(i)) * part%volc(i) ! Use "diss" here just to select the soluble species
       zvtot = zvtot + part%volc(i)
    END DO
    
    ! Combine the two cases from original code since they're exactly the same??
    IF (zvw > 1.e-28*part%numc .OR. zvs > 1.e-28*part%numc) THEN	
       ! Aqueous droplet OR dry partially soluble particle
       calc_Sw_eq = (znw/(zns+znw)) * exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
    ELSE IF (zvtot-zvs > 1.e-28*part%numc) THEN
       ! Dry insoluble particle
       calc_Sw_eq = exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
    ELSE
       ! Just add eps to avoid divide by zero
       calc_Sw_eq = (znw/(eps+zns+znw)) * exp(4.*surfw0*spec%mwa/(rg*T*spec%rhowa*dwet))
    END IF
    
  END FUNCTIOn calc_Sw_eq
  
  
  ! Function for calculating Pearson's correlation coefficient for two vectors
  REAL FUNCTION calc_correlation(x,y,n)
    INTEGER :: n
    REAL :: x(n), y(n)
    REAL :: sx, sy, sx2, sy2, sxy
    INTEGER :: i
    sx=0.; sy=0.; sx2=0.; sy2=0.; sxy=0.
    DO i=1,n
       sx=sx+x(i)
       sy=sy+y(i)
       sx2=sx2+x(i)**2
       sy2=sy2+y(i)**2
       sxy=x(i)*y(i)
    ENDDO
    IF (sx2*n-sx**2<eps .OR. sy2*n-sy**2<eps) THEN
       calc_correlation = 0.
    ELSE
       calc_correlation = ( sxy*n-sx*sy )/( SQRT(sx2*n-sx**2)*SQRT(sy2*n-sy**2) )
    ENDIF
  END FUNCTION calc_correlation
  
END MODULE mo_submctl
