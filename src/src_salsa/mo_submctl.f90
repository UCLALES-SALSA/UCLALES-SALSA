MODULE mo_submctl
  IMPLICIT NONE

  ! ---------------------------------------------------------------------------------------------------------
  ! 1) Data type definitions

  ! Data type used to store information about aerosol, cloud, rain, ice, and snow size bins
  INTEGER, PARAMETER :: maxnspec=11
  TYPE t_section
     REAL :: vhilim,     & ! bin volume at the high limit
                 vlolim,     & ! - '' - at the low limit
                 vratiohi,   & ! volume ratio between the center and high limit
                 vratiolo,   & ! - '' - and the low limit
                 dmid,       & ! bin middle diameter
                 !******************************************************
                 ! ^ Do NOT change the stuff above after initialization !
                 !******************************************************
                 dwet,       & ! Wet diameter or mean droplet diameter
                 volc(maxnspec),  & ! Volume concentrations of water + aerosol species (water always the first)
                 numc          ! Number concentration of particles/droplets
  END TYPE t_section


  ! ---------------------------------------------------------------------------------------------------------
  ! 2) Switches for SALSA aerosol microphysical processes

  ! Process switches nl* are read from the NAMELIST and NOT changed during runtime.
  ! ls* is the switch actually used and will get the value of nl* except for special
  ! circumstances such as spinup period etc.

  LOGICAL :: nlcoag  = .TRUE., lscoag ! Coagulation master switch
  LOGICAL :: nlcgaa  = .TRUE., lscgaa ! Coagulation between aerosols
  LOGICAL :: nlcgcc  = .TRUE., lscgcc ! Collision-coalescence between cloud droplets
  LOGICAL :: nlcgca  = .TRUE., lscgca ! Cloud collection of aerosols
  LOGICAL :: nlcgpc  = .TRUE., lscgpc ! Collection of cloud droplets by rain
  LOGICAL :: nlcgpa  = .TRUE., lscgpa ! Collection of aerosols by rain
  LOGICAL :: nlcgpp  = .TRUE., lscgpp ! Collision between rain drops
  LOGICAL :: nlcgia  = .TRUE., lscgia ! Ice collection of aerosols
  LOGICAL :: nlcgic  = .TRUE., lscgic ! Collection of cloud droplets by ice particles
  LOGICAL :: nlcgii  = .TRUE., lscgii ! Collision-coalescence between ice particles
  LOGICAL :: nlcgip  = .TRUE., lscgip ! Collection of precipitation by ice particles
  LOGICAL :: nlcgsa  = .TRUE., lscgsa ! Collection of aerosols by snow
  LOGICAL :: nlcgsc  = .TRUE., lscgsc ! Collection of cloud droplets by snow
  LOGICAL :: nlcgsi  = .TRUE., lscgsi ! Collection of ice by snow
  LOGICAL :: nlcgsp  = .TRUE., lscgsp ! Collection of precipitation by snow
  LOGICAL :: nlcgss  = .TRUE., lscgss ! Collision-coalescence between snow particles

  LOGICAL :: nlcnd      = .TRUE.,  lscnd      ! Condensation master switch
  LOGICAL :: nlcndh2ocl = .TRUE.,  lscndh2ocl ! Condensation of water vapour on clouds and precipitation
  LOGICAL :: nlcndh2oae = .TRUE.,  lscndh2oae ! Condensation of water vapour on aerosol particles (FALSE -> equilibrium calc.)
  LOGICAL :: nlcndh2oic = .TRUE.,  lscndh2oic ! Condensation of water vapour on ice and snow
  LOGICAL :: nlcndgas   = .FALSE., lscndgas   ! Condensation of H2SO4 and organic vapors

  LOGICAL :: nlauto     = .TRUE.,  lsauto     ! Autoconversion of cloud droplets (needs activation)
  LOGICAL :: nlautosnow = .TRUE.,  lsautosnow ! Autoconversion of ice particles to snow (needs activation)

  LOGICAL :: nlactiv    = .TRUE.,  lsactiv    ! Cloud droplet activation master switch
  LOGICAL :: nlactintst = .TRUE.,  lsactintst ! Switch for interstitial activation
  LOGICAL :: nlactbase  = .FALSE., lsactbase  ! Switch for cloud base activation

  LOGICAL :: nlicenucl  = .FALSE., lsicenucl  ! ice nucleation master switch
  LOGICAL :: nlicmelt   = .FALSE., lsicmelt   ! ice melting

  ! Other switches
  LOGICAL :: lsdistupdate = .TRUE.  ! Perform the size distribution update


  ! ---------------------------------------------------------------------------------------------------------
  ! 3) Parameters and options for SALSA microphysics

  ! New particle formation or nucleation
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


  ! Autoconversion
  !   Cloud to rain
  REAL :: autoc_rain_zd0 = 50.e-6 ! Cloud-rain diameter limit
  REAL :: autoc_rain_sigmag = 1.2 ! Assumed log-normal cloud drop size distribution width
  LOGICAL :: auto_sb = .FALSE. ! Use the Seifert & Beheng (2001) autoconversion method
  !   Ice to snow
  REAL :: autoc_snow_zd0 = 250.e-6 ! Ice-snow diameter limit
  REAL :: autoc_snow_sigmag = 1.2 ! Assumed log-normal ice particle size distribution width


  ! Options for ice nucleation (when master switch nlicenucl = .TRUE,)
  ! a) Constant ice number concentration (fixinc > 0 #/kg) is maintained by converting cloud droplets to ice/snow
  REAL :: fixinc = -1.0 ! Default = disabled
  ! Cloud freezing order: >0: start from the largest bin, 0: all bins evenly, <0: start from the smallest bin
  INTEGER :: ice_source_opt = 1 ! Default = start from the largest bin
  ! b) Modelled ice nucleation
  LOGICAL :: ice_hom = .FALSE., ice_imm=.FALSE., ice_dep=.FALSE. ! Available ice nucleation modes
  ! c) Start time (s) for ice formation
  REAL :: icenucl_tstart = 0. ! Default: right after initialization
  ! d) Where to put new ice/snow: <0: parallel ice bin, 0: find matching snow bin, >0 snow bin specified by ice_target_opt
  INTEGER :: ice_target_opt = -1 ! Default = parallel ice bins


  ! Gas phase parameters
  INTEGER, PARAMETER :: maxngas=15
  ! Total number of prognostic gas phase species
  INTEGER :: ngases = 0
  ! Names of the prognostic gas phase species
  CHARACTER(len=3) :: zgas(maxngas) = '   '
  ! Molecular weights (kg/mol)
  REAL :: mws_gas(maxngas) = 0.1
  ! Total number of diagnostic gas phase species
  INTEGER :: ngases_diag = 0
  ! Concentration array for the diagnostic gas phase species (mol/kg)
  REAL :: zgas_diag(maxngas) = 0.

  ! Simple sulfate (H2SO4(g) => SO4) and non-volatile organic vapor (LVOA(g)  => OC) partitioning
  REAL :: conc_h2so4=-1., conc_ocnv=-1.            ! Initial concentrations
  LOGICAL :: part_h2so4=.FALSE., part_ocnv=.FALSE. ! Calculated when non-negative input concentrations
  INTEGER :: isog=1, iocg=1                        ! Index to gas phase

  ! Detailed SOA formation including gas phase oxidants, VOCs and VBS and aqSOA species
  INTEGER :: nvbs_setup = -1   ! VBS setup option (hard-coded schemes)
  LOGICAL :: laqsoa = .FALSE.  ! Enable aqSOA formation
  INTEGER :: nvocs=0, nvbs=0, naqsoa=0 ! Number of VOCs, VBS species and aqSOA species
  REAL :: conc_voc(maxngas)=0., conc_vbsg(maxngas)=0., conc_aqsoag(maxngas)=0. ! Initial concetrations (kg/kg)
  REAL :: zvbs_k_OH = 0., zvbs_Eact_p_OH =0. ! Rate coefficient for VBS(g) aging [cm3/#/s]
  ! VOC oxidants
  LOGICAL :: ox_prescribed = .TRUE.              ! Prescribed or prognostic concentration fields
  INTEGER :: id_oh=-1, id_o3=-1, id_no3=-1       ! Indexes to gas phase
  REAL :: conc_oh=-1., conc_o3=-1., conc_no3=-1. ! Initial concentrations (number mixing ratio)
  REAL :: zdayfac_oh=1., zdayfac_o3=1., znightfac_no3=1. ! Scaling factors for diurnal oxidant concentrations
  REAL :: zphotofac_aqsoa=0.                     ! The same for aqSOA photodissociation
  INTEGER :: ox_conc_flag=0, aqsoa_photo_flag=0  ! Flags related to the scaling factors
  ! Communication between LES and SALSA/VBS
  REAL :: model_lat=31.5, & ! Mean latitude from the model domain [deg]
      start_doy=-1.    ! Decimal day of year [-]


  ! ---------------------------------------------------------------------------------------------------------
  ! 4) Other inputs and parameters

  ! RH Limit: used for initialization and spinup within SALSA to limit the water vapour mixing ratio.
  ! Prevents unrealistically high RH in cloud activation and condensation procedures that is often assigned
  ! in the LES input files to immediately generate cloud. Given in %/100.
  REAL :: rhlim = 1.20

  ! Eddy dissipation rate - parameter for the turbulent component in coagulation
  ! Negative value means take from LES, otherwise constant
  real :: eddy_dis_rt = -10.e-4

  ! Define which aerosol species are used and their initial size distributions
  ! Initial aerosol species
  INTEGER :: nspec = 1 ! Does not include water
  INTEGER, PARAMETER :: maxspec = 8
  CHARACTER(len=3) :: listspec(maxspec) = (/'SO4','   ','   ','   ','   ','   ','   ','   '/)

  ! Volume fractions between aerosol species for A and B-bins
  REAL :: volDistA(maxspec) = (/1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  REAL :: volDistB(maxspec) = (/0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0/)
  ! Limit 1a composition to OC and/or SO4
  LOGICAL :: salsa1a_SO4_OC = .TRUE.
  ! Number fraction allocated to a-bins in regime 2 (b-bins will get 1-nf2a)
  REAL :: nf2a = 1.0

  ! Type of the input aerosol size distribution
  !     0 - Uniform, log-normal size distribution parameters given in the NAMELIST
  !     1 - Read vertical profile of those from an input file
  INTEGER :: isdtyp = 0
  ! For isdtyp = 0
  INTEGER, PARAMETER :: nmod = 7
  REAL :: sigmag(nmod) = (/2.0,2.0,2.0,2.0,2.0,2.0,2.0/),   & ! Stdev
             dpg(nmod) = (/0.03,0.15,0.2,0.2,0.2,0.2,0.2/), & ! Mode diam in um
               n(nmod) = (/1600.,640.,0.,0.,0.,0.,0./)        ! 1e6#/kg ~ #/cm3

  ! Aerosol, cloud and ice bin limits (based on dry size)
  INTEGER, PARAMETER :: maxnreg = 5 ! maximum number of subregimes (the first is region 1 and the rest are for region 2)
  REAL :: reglim(maxnreg+1) = (/ 3.e-9, 5.e-8, 7.e-7, 1.e-5, 0., 0. /) ! low/high diameter limits of main size regimes [m]
  INTEGER :: nbin(maxnreg) = (/ 3, 4, 3, 0, 0 /)   ! number of bins in each main regime

  ! Rain and snow bin limits
  !  For example, this line in the NAMELIST would set the default rain bins:
  !    rainbinlim(1:8)=50.,55.,65.,100.,200.,500.,1000.,2000.
  REAL :: rainbinlim(100) = -1.
  REAL :: snowbinlim(100) = -1.

  ! Cloud and ice bin limits for histogram outputs (profiles only)
  REAL :: cldbinlim(100) = -1., icebinlim(100) = -1.
  INTEGER :: nout_cld=0, nout_ice=0


  ! ---------------------------------------------------------------------------------------------------------
  ! 5) Global SALSA parameters

  ! Aerosol bins
  INTEGER ::      & ! start and end indexes for aerosol bin regimes
   in1a, fn1a,    & ! regime 1a
   in2a, fn2a,    & ! regime 2a
   in2b, fn2b,    & ! regime 2b
   nbins            ! total number of aerosol size bins

   ! Cloud, rain, ice and snow bins
   !    cloud & ice: no 1a bins
   !    rain & snow: bins based on wet size
  INTEGER :: inp2a, fnp2a, inp2b, fnp2b ! Bin indexes for cloud and ice (no 1a, just 2a and 2b)
  INTEGER :: ncld, nprc, nice, nsnw ! Total number of cloud, rain, ice and snow bins

  ! These are just to deliver information about the bin limits (radius) to the host model.
  REAL, ALLOCATABLE :: aerobins(:), precpbins(:), snowbins(:)

  ! Indexes linking named specied to the active species (-1 = not used)
  INTEGER :: iso=-1, ioc=-1, ibc=-1, idu=-1, iss=-1, inh=-1, ino=-1, ih2o=-1
  ! Chemical and physics properties alls of aerosol (cloud, ice, ...) species
  REAL :: dens(maxnspec)=1000., mws(maxnspec)=0.1, diss(maxnspec)=0.
  ! Names, e.g. 'SO4'
  CHARACTER(len=3) :: zspec(maxnspec)='   '


  ! ---------------------------------------------------------------------------------------------------------
  ! 6) Constants

  REAL, PARAMETER :: &
   nlim   = 1.e-3,        & ! Number conc. limit (#/m3) for aerosol and cloud droplets
   prlim  = 1.e-6,        & ! The same for precipitation and ice species for which concentrations are normally much lower [#/m3]
   eps    = epsilon(1.0)    ! epsilon

  REAL, PARAMETER :: &
   pi     = 3.1415927,    & ! pi
   pi6    = 0.5235988,    & ! pi/6
   avog   = 6.0221e+23,   & ! Avogadro number (#/mol)
   boltz  = 1.3807e-23,   & ! Boltzmann constant (J/K)
   planck = 6.626070040e-34, & ! Planck constant (J*s)
   grav   = 9.81,         & ! gravitational acceleration (m/s^2)
   pstand = 1.01325e+5,   & ! standard pressure (Pa)
   rg     = 8.314,        & ! molar gas constant (J/(mol K))
   rda    = 287.04,       & ! gas constant for dry air (J/K/kg)
   cpa    = 1005.,        & ! specific heat of dry air, constant P (J/kg/K)
   alv    = 2.5e6,        & ! latent heat for vaporisation (J/kg)
   als    = 2.834e6,      & ! latent heat for sublimation (J/kg)
   mair   = 28.967e-3,    & ! molar mass of air (mol/kg)
   surfw0 = 0.073,        & ! surface tension of pure water @ ~ 293 K [J/m2]
   surfi0 = 0.105           ! surface tension of ice

  REAL, SAVE :: & ! molar mass [kg/mol], dissociation constant [-], and density [kg/m3]
   msu = 98.08e-3,  disssu = 3.0, rhosu = 1830., & ! sulphate
   mno = 62.01e-3,  dissno = 1.0, rhono = 1479., & ! HNO3
   mnh = 18.04e-3,  dissnh = 1.0, rhonh = 1530., & ! NH3
   moc = 150.e-3,   dissoc = 1.0, rhooc = 2000., & ! organic carbon
   mbc = 12.e-3,    dissbc = 0.0, rhobc = 2000., & ! black carbon
   mss = 58.44e-3,  dissss = 2.0, rhoss = 2165., & ! sea salt (NaCl)
   mdu = 100.e-3,   dissdu = 0.0, rhodu = 2650., & ! mineral dust
   mwa = 18.016e-3, disswa = 1.0, rhowa = 1000., & ! water
   rhoic = 917.,    rhosn = 300. ! densitities of ice and snow

  REAL, PARAMETER :: & ! diameter of condensing molecule [m]
   d_sa = 5.539376964394570e-10, & ! H2SO4
   d_oc = 6.195906936656752e-10    ! Non-volatile organics

  REAL, PARAMETER :: n3 = 158.79  ! number of H2SO4 molecules in 3 nm cluster assuming d_sa = 5.54 ???

  REAL, ALLOCATABLE :: massacc(:)

contains
  !********************************************************************
  ! Function for calculating terminal velocities for different particle types and size ranges.
  !     Tomi Raatikainen (2.5.2017)
  REAL FUNCTION terminal_vel(radius,rhop,rhoa,visc,beta,flag)
    implicit none
    REAL, intent(in) :: radius, rhop ! Particle radius and density
    REAL, intent(in) :: rhoa, visc, beta ! Air density, viscocity and Cunningham correction factor
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    ! Constants
    real, parameter :: rhoa_ref = 1.225 ! reference air density (kg/m^3)

    IF (flag==4) THEN   ! Ice
        ! Ice crystal terminal fall speed from Ovchinnikov et al. (2014)
        !       Dimension D = 2*radius
        terminal_vel = 12.0*sqrt(2.0*radius)
    ELSEIF (flag==5) THEN   ! Snow
        ! The same for snow
        !       Dimension D = 2*radius
        terminal_vel = 12.0*sqrt(2.0*radius)
    ELSE
        ! Aerosol and cloud and rain droplets
        IF (radius<40.0e-6) THEN
            ! Stokes law with Cunningham slip correction factor
            terminal_vel = (4.*radius**2)*(rhop-rhoa)*grav*beta/(18.*visc) ![m s-1]
        ELSEIF (radius<0.6e-3) THEN
            ! Droplets from 40 um to 0.6 mm: linear dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            terminal_vel = 8.e3*radius*sqrt(rhoa_ref/rhoa)
        ELSE
            ! Droplets larger than 0.6 mm: square root dependence on particle radius and a correction for reduced pressure
            !   R.R. Rogers: A Short Course in Cloud Physics, Pergamon Press Ltd., 1979.
            ! Note: this is valid up to 2 mm or 9 m/s (at 1000 mbar), where droplets start to break
            terminal_vel = 2.01e2*sqrt( min(radius,2.0e-3)*rhoa_ref/rhoa)
        ENDIF
    ENDIF
  END FUNCTION terminal_vel

  !********************************************************************
  ! Functions for calculating dimension (or wet diameter) for any particle type
  ! - Aerosol, cloud and rain are spherical
  ! - Snow and ice can be irregular and their densities can be size-dependent
  !
  ! Edit these functions when needed
  !
  ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition
  ! and coagulation) and capacitance (condensation). Otherwise compact spherical structure can be expected.
  !
  ! This function is for SALSA t_section arrays and assumes volume-based concentration units
  SUBROUTINE CalcDimension(n,ppart,lim,flag)
    IMPLICIT NONE
    INTEGER, INTENT(in) :: n
    TYPE(t_section), INTENT(inout) :: ppart(n)
    REAL, INTENT(IN) :: lim
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    INTEGER i

    ppart(:)%dwet = 2.e-10
    DO i=1,n
        IF (ppart(i)%numc>lim) &
            ppart(i)%dwet=(SUM(ppart(i)%volc(:))/ppart(i)%numc/pi6)**(1./3.)
    ENDDO

  END SUBROUTINE CalcDimension
  !
  ! This function is for SALSA t_section arrays and assumes volume-based concentration units
  SUBROUTINE CalcMass(mass,n,ppart,lim,flag)
    IMPLICIT NONE
    REAL, INTENT(OUT) :: mass(n)
    INTEGER, INTENT(in) :: n
    TYPE(t_section), INTENT(in) :: ppart(n)
    REAL, INTENT(IN) :: lim
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)
    INTEGER i

    mass(:)=1e-30
    DO i=1,n
        IF (ppart(i)%numc<lim) THEN
            ! No particles
        ELSEIF (flag==4) THEN   ! Ice
            mass(i)=(ppart(i)%volc(1)*rhoic+SUM(ppart(i)%volc(2:)*dens(2:)))/ppart(i)%numc
        ELSEIF (flag==5) THEN   ! Snow
            mass(i)=(ppart(i)%volc(1)*rhosn+SUM(ppart(i)%volc(2:)*dens(2:)))/ppart(i)%numc
        ELSE
            mass(i)=SUM(ppart(i)%volc(:)*dens(:))/ppart(i)%numc
        ENDIF
    ENDDO

  END SUBROUTINE CalcMass
  !
  ! This function is for single LES size bin and assumes that concentration is given as mass per particle
  REAL FUNCTION calc_eff_radius(n,mass,flag)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n ! Number of chemical species
    REAL, INTENT(IN) :: mass(n) ! Mass (kg) per particle
    INTEGER, INTENT(IN) :: flag ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)

    IF (flag==4) THEN   ! Ice
        ! Spherical ice
        calc_eff_radius=0.5*( (mass(1)/rhoic+SUM(mass(2:)/dens(2:n)))/pi6)**(1./3.)
    ELSEIF (flag==5) THEN   ! Snow
        ! Spherical snow
        calc_eff_radius=0.5*( (mass(1)/rhosn+SUM(mass(2:)/dens(2:n)))/pi6)**(1./3.)
    ELSE
        ! Radius from total volume of a spherical particle or aqueous droplet
        calc_eff_radius=0.5*( SUM(mass(:)/dens(1:n))/pi6)**(1./3.)
    ENDIF

  END FUNCTION calc_eff_radius

  !********************************************************************
  ! Function for calculating equilibrium water saturation ratio at droplet surface based on Kohler theory
  !
  REAL FUNCTION calc_Sw_eq(part,T)
    TYPE(t_section), INTENT(in) :: part ! Any particle
    REAL, INTENT(IN) :: T ! Absolute temperature (K)
    REAL :: dwet

    !   #  Name Diss
    !   1   SO4   3
    !   2   OC     1
    !   3   BC     0
    !   4   DU    0
    !   5   SS     2
    !   6   NO    1
    !   7   NH    1
    !   8   H2O

    ! Wet diameter
    dwet=(SUM(part%volc(:))/part%numc/pi6)**(1./3.)

    ! Equilibrium saturation ratio = xw*exp(4*sigma*v_w/(R*T*Dwet))
    IF (part%volc(1)>1e-28*part%numc) THEN
        ! An aqueous droplet
        calc_Sw_eq=part%volc(1)*rhowa/mwa/SUM(part%volc(:)*diss(:)*dens(:)/mws(:))* &
            exp(4.*surfw0*mwa/(rg*T*rhowa*dwet))
    ELSEIF (SUM(part%volc(2:)*diss(2:))>1e-28*part%numc) THEN
        ! Dry particle with soluble substances: allow complete dissolution (the same equation as for aqueous droplets)
        calc_Sw_eq=part%volc(1)*rhowa/mwa/SUM(part%volc(:)*diss(:)*dens(:)/mws(:))* &
            exp(4.*surfw0*mwa/(rg*T*rhowa*dwet))
    ELSEIF (SUM(part%volc(2:))>1e-28*part%numc) THEN
        ! Dry insoluble particle: xw = 1 even with trace amounts of water
        calc_Sw_eq=exp(4.*surfw0*mwa/(rg*T*rhowa*dwet))
    ELSE
        ! Just add eps to avoid divide by zero
        calc_Sw_eq=part%volc(1)*rhowa/mwa/(eps+SUM(part%volc(:)*diss(:)*dens(:)/mws(:)))* &
            exp(4.*surfw0*mwa/(rg*T*rhowa*dwet))
    ENDIF

  END FUNCTIOn calc_Sw_eq

  ! Function for calculating Pearson's correlation coefficient for two vectors
  REAL FUNCTION calc_correlation(x,y,n)
    INTEGER :: n
    REAL :: x(n), y(n)
    REAL :: sx, sy, sx2, sy2, sxy
    INTEGER :: i
    IF (n<=1) THEN
        calc_correlation = 0.
    ELSE
        calc_correlation = (SUM(x*y)*n-SUM(x)*SUM(y))/( SQRT(SUM(x**2)*n-SUM(x)**2)*SQRT(SUM(y**2)*n-SUM(y)**2) )
    ENDIF
    RETURN
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

  ! Function for finding aerosol and gas phase species indices based on their names.
  ! Zero means that species is not found.
  INTEGER FUNCTION find_spec_id(spec)
    CHARACTER(len=*) :: spec
    !
    IF (nspec>0 .AND. LEN_TRIM(spec)>0) THEN
        DO  find_spec_id=1,nspec+1
            IF (zspec(find_spec_id)==spec) RETURN
        ENDDO
    ENDIF
    find_spec_id=0
  END FUNCTION find_spec_id
  !
  INTEGER FUNCTION find_gas_id(spec)
    CHARACTER(len=*) :: spec
    !
    IF (ngases>0 .AND. LEN_TRIM(spec)>0) THEN
        DO  find_gas_id=1,ngases
            IF (zgas(find_gas_id)==spec) RETURN
        ENDDO
    ENDIF
    find_gas_id=0
  END FUNCTION find_gas_id

  ! Function for setting scaling factors for oxidant concentrations and aqSOA photodissociation
  SUBROUTINE set_vbs_diag(etime)
    REAL, INTENT(IN) :: etime ! Elapsed time (s)
    !
    ! Oxidants - only when diagnostic
    IF (ox_prescribed) THEN
        ! Calculate zdayfac_oh, zdayfac_o3 and znightfac_no3 based on ox_conc_flag
        IF (ox_conc_flag==1) THEN
            ! Diurnal scaling for OH and NOx
            CALL calc_day_night_fac(etime,zdayfac_oh,znightfac_no3)
        ELSEIF (ox_conc_flag==2) THEN
            ! Diurnal scaling for all
            CALL calc_day_night_fac(etime,zdayfac_oh,znightfac_no3)
            zdayfac_o3=zdayfac_oh
        ELSE
            ! Default: constant as is
        ENDIF
    ENDIF
    !
    ! aqSOA photodissociation
    IF (naqsoa>0) THEN
        IF (aqsoa_photo_flag==1) THEN
            ! Scaling based on cosine of the solar zenith angle; from 1 to 0
            zphotofac_aqsoa=MAX(0., zenith(0., start_doy+etime/86400.) )
        ELSEIF (aqsoa_photo_flag==2) THEN
            ! Scaling based on cosine of the solar zenith angle; limits based on latitude
            zphotofac_aqsoa=MAX(0., zenith(model_lat, start_doy+etime/86400.) )
        ELSE
            ! Default: constant as is
        ENDIF
    ENDIF
    !
    CONTAINS
        ! Scaling based on 6.073e-5 * u0**1.743 * exp(-0.474 / u0),
        ! where u0 is cosine of the solar zenith angle
        SUBROUTINE calc_day_night_fac(etime,zdayfac,znightfac)
            REAL, INTENT(IN) :: etime ! Elapsed time (s)
            REAL, INTENT(OUT) :: zdayfac, znightfac
            ! Local parameters
            REAL, SAVE :: rate_ave=-999., maxdayfac=-999.
            ! Local variables
            REAL :: u0, rate
            INTEGER :: i
            !
            ! Initilization
            IF (rate_ave<-1.) THEN
                rate_ave = 0.
                maxdayfac = 0.
                do i = 0, 23
                    u0=zenith(model_lat, start_doy+i/24.)
                    IF (u0>0.)then
                        rate = 6.073e-5 * u0**1.743 * exp(-0.474 / u0)
                        if(rate > maxdayfac) maxdayfac = rate
                        rate_ave = rate_ave + rate
                    endif
                enddo
                rate_ave = rate_ave / 24.
                maxdayfac = maxdayfac / rate_ave
            ENDIF
            !
            ! Cosine of the solar zenith angle
            u0=zenith(model_lat, start_doy+etime/86400.)
            IF (u0<=0.)then
                zdayfac = 0.
            else
                zdayfac = 6.073e-5 * u0**1.743 * exp(-0.474 / u0)
                zdayfac = zdayfac / rate_ave
            endif
            znightfac = (maxdayfac - zdayfac) / (maxdayfac - 1.)
        END SUBROUTINE calc_day_night_fac
        !
        ! From LES/rad_driver.f90
        ! ---------------------------------------------------------------------------
        ! Return the cosine of the solar zenith angle give the decimal day and
        ! the latitude
        !
        real function zenith(alat,time)

            real, intent (in)  :: alat, time

            real, parameter :: pi     = 3.14159265358979323846264338327
            real :: lamda, d, sig, del, h, day

            day    = floor(time)
            lamda  = alat*pi/180.
            d      = 2.*pi*int(time)/365.
            sig    = d + pi/180.*(279.9340 + 1.914827*sin(d) - 0.7952*cos(d) &
                 &                      + 0.019938*sin(2.*d) - 0.00162*cos(2.*d))
            del    = asin(sin(23.4439*pi/180.)*sin(sig))
            h      = 2.*pi*((time-day)-0.5)
            zenith = sin(lamda)*sin(del) + cos(lamda)*cos(del)*cos(h)

        end function zenith
  END SUBROUTINE set_vbs_diag

END MODULE mo_submctl
