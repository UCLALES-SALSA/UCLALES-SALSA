!****************************************************************
!*                                                              *
!*   module MO_SALSA_INIT                                       *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to initialize the particle size grid and aerosol           *
!*   processes                                                  *
!*                                                              *
!****************************************************************
MODULE mo_salsa_init

  IMPLICIT NONE

CONTAINS

  ! fxm: when dpmid is used for calculating coagulation coefficients
  ! (only?), would it make more sense to use approximate wet radii
  ! e.g. for sea salt particles?
  !********************************************************************
  !
  ! subroutine SET_SIZEBINS()
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Initializes particle size distribution grid by
  ! calculating size bin limits and mid-size for
  ! *dry* particles in each bin
  !
  !
  ! Method:
  ! -------
  ! Size distribution described using
  !   1) moving center method (regimes 1 and 2)
  !   (Jacobson, Atmos. Env., 31, 131-144, 1997)
  !   2) fixed sectional method (regime 3)
  !
  ! Size bins in each regime are spaced logarithmically
  ! based on given regime size limits and bin number.
  !
  !
  ! Interface:
  ! ----------
  ! Called from model driver
  ! (only at the beginning of simulation)
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005
  ! Harri Kokkola (FMI) 2006
  !
  ! Bug fixes for box model + updated for the new aerosol datatype:
  ! ----------
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------

  SUBROUTINE set_sizebins()

    USE mo_submctl, ONLY : &
         pi6,         & ! pi/6
         reglim,      & ! diameter limits for size regimes [m]
         nbin,        & ! number of size bins in each (sub)regime
         in1a, fn1a,  & ! size regime bin indices: 1a
         in2a, fn2a,  & !     - " -       2a
         in2b, fn2b,  & !     - " -       2b
         nbins,  &
         aerobins

    USE mo_salsa_driver, ONLY : &
         kbdim,klev,  &
         aero

    IMPLICIT NONE

    !-- local variables ----------
    INTEGER :: ii, jj,cc,dd,vv ! loop indices
    INTEGER :: nbin2, nbin3
    REAL ::  ratio ! ratio of regime upper and lower diameter

    nbin2 = 4
    nbin3 = nbin(2) - nbin2

    ! -------------------------
    ! Allocate aerosol tracers
    !--------------------------
    ALLOCATE(aero(kbdim,klev,nbins))

    ! Juha: Corrected some bugs with the section volume definitions
    !-------------------------------------------------------------------------------

    DO jj = 1,klev
       DO ii = 1,kbdim

          !-- 1) size regime 1: --------------------------------------
          !  - minimum & maximum *dry* volumes [fxm]
          !  - bin mid *dry* diameter [m]
          ratio = reglim(2)/reglim(1)   ! section spacing

          DO cc = in1a,fn1a
             aero(ii,jj,cc)%vlolim = pi6*(reglim(1)*ratio**(real(cc-1)/nbin(1)))**3
             aero(ii,jj,cc)%vhilim = pi6*(reglim(1)*ratio**(real(cc)/nbin(1)))**3
             aero(ii,jj,cc)%dmid = ( (aero(ii,jj,cc)%vhilim + aero(ii,jj,cc)%vlolim) /  &
                                     (2.*pi6) )**(1./3.)
             aero(ii,jj,cc)%vratiohi = aero(ii,jj,cc)%vhilim/(pi6*aero(ii,jj,cc)%dmid**3)
             aero(ii,jj,cc)%vratiolo = aero(ii,jj,cc)%vlolim/(pi6*aero(ii,jj,cc)%dmid**3)
          END DO

          !-- 2) size regime 2: --------------------------------------
          !  - minimum & maximum *dry* volumes [fxm]
          !  - bin mid *dry* diameter [m]

          !-- 2.1) first for subregime 2a
          ratio = reglim(3)/reglim(2)   ! section spacing

          DO dd = in2a,in2a+nbin2-1
             cc = dd - in2a
             aero(ii,jj,dd)%vlolim = pi6*(reglim(2)*ratio**(real(cc)/nbin2))**3
             aero(ii,jj,dd)%vhilim = pi6*(reglim(2)*ratio**(real(cc+1)/nbin2))**3
             aero(ii,jj,dd)%dmid = ( (aero(ii,jj,dd)%vhilim + aero(ii,jj,dd)%vlolim) /  &
                                     (2.*pi6) )**(1./3.)
             aero(ii,jj,dd)%vratiohi = aero(ii,jj,dd)%vhilim/(pi6*aero(ii,jj,dd)%dmid**3)
             aero(ii,jj,dd)%vratiolo = aero(ii,jj,dd)%vlolim/(pi6*aero(ii,jj,dd)%dmid**3)
          END DO

          !-- 3) size regime 3: --------------------------------------
          !  - bin mid *dry* diameter [m]
          ratio = reglim(4)/reglim(3)   ! section spacing

          DO dd = in2a+nbin2,fn2a
             cc = dd - (fn2a-(nbin3-1))

             aero(ii,jj,dd)%vlolim = pi6*(reglim(3)*ratio**(real(cc)/nbin3))**3
             aero(ii,jj,dd)%vhilim = pi6*(reglim(3)*ratio**(real(cc+1)/nbin3))**3
             aero(ii,jj,dd)%dmid = ( (aero(ii,jj,dd)%vhilim + aero(ii,jj,dd)%vlolim) /  &
                                     (2.*pi6) )**(1./3.)
             aero(ii,jj,dd)%vratiohi = aero(ii,jj,dd)%vhilim/(pi6*aero(ii,jj,dd)%dmid**3)
             aero(ii,jj,dd)%vratiolo = aero(ii,jj,dd)%vlolim/(pi6*aero(ii,jj,dd)%dmid**3)
          END DO

          !-- 2.2) same values for subregime 2b
          aero(ii,jj,in2b:fn2b)%vlolim = aero(ii,jj,in2a:fn2a)%vlolim
          aero(ii,jj,in2b:fn2b)%vhilim = aero(ii,jj,in2a:fn2a)%vhilim
          aero(ii,jj,in2b:fn2b)%dmid = aero(ii,jj,in2a:fn2a)%dmid
          aero(ii,jj,in2b:fn2b)%vratiohi = aero(ii,jj,in2a:fn2a)%vratiohi
          aero(ii,jj,in2b:fn2b)%vratiolo = aero(ii,jj,in2a:fn2a)%vratiolo

          ! Initialize the wet diameter with the bin dry diameter to avoid numerical proplems later
          aero(ii,jj,:)%dwet = aero(ii,jj,:)%dmid
          aero(ii,jj,:)%core = pi6*aero(ii,jj,:)%dmid**3

          ! Set volume and number concentrations to zero
          aero(ii,jj,:)%numc = 0.
          DO vv=1,8
             aero(ii,jj,:)%volc(vv) = 0.
          END DO

       END DO !ii
    END DO !!jj

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(aerobins(nbins))
    DO cc = 1,nbins
       aerobins(cc) = (aero(1,1,cc)%vlolim/pi6)**(1./3.)
    END DO
    aerobins = 0.5*aerobins ! to radius

  END SUBROUTINE set_sizebins

  !--------------------------------------------------------------------------
  !
  ! *******************************
  ! SUBROUTINE set_cloudbins
  ! *******************************
  !
  ! Setup of hydrometeor size bins similar to the subroutine *set_sizebins*
  !
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------------

  SUBROUTINE set_cloudbins()
    USE mo_submctl, ONLY : pi6,             &
                               ica,icb,         &
                               fca,fcb,         &
                               ira,fra,         &
                               ncld,nprc,        &
                               in2a,fn2a,       &
                               fn2b,       &
                               cloudbins,       &
                               precpbins
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                                 cloud,precp,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc,nba,nbb

    REAL :: tmplolim(7), tmphilim(7)

    ! Helper arrays to set up precipitation size bins
    tmplolim = (/50.,55.,65.,100.,200.,500.,1000./)*1.e-6
    tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6

    ! Cloud bins are parallel with the 2a and 2b aerosol bins

    ! Number of cloud bins in regime a
    nba = fn2a-in2a+1
    ! Number of cloud bins in regime b
    nbb = nba

    ! Reset cloud bin indices accordingly. The two components give the cloud regime index,
    ! and the aerosol bin index with which they are parallel
    ica%cur = 1;                      ica%par = in2a
    fca%cur = ica%cur + nba-1; fca%par = ica%par + nba-1
    icb%cur = fca%cur + 1;            icb%par = fn2b - nbb + 1
    fcb%cur = icb%cur + nbb-1; fcb%par = icb%par + nbb-1

    ncld = fcb%cur

    ! Rain/drizzle bins
    ira = 1; fra = 7;
    nprc = fra

    ! ----------------------------------------
    ! Allocate cloud and precipitation arrays
    ! ----------------------------------------
    ALLOCATE(cloud(kbdim,klev,ncld), precp(kbdim,klev,nprc))

    DO jj = 1,klev
       DO ii = 1,kbdim

          ! -------------------------------------------------
          ! Set cloud properties (parallel to aerosol bins)
          ! -------------------------------------------------
          cloud(ii,jj,ica%cur:fca%cur)%vhilim = aero(ii,jj,ica%par:fca%par)%vhilim
          cloud(ii,jj,icb%cur:fcb%cur)%vhilim = aero(ii,jj,icb%par:fcb%par)%vhilim

          cloud(ii,jj,ica%cur:fca%cur)%vlolim = aero(ii,jj,ica%par:fca%par)%vlolim
          cloud(ii,jj,icb%cur:fcb%cur)%vlolim = aero(ii,jj,icb%par:fcb%par)%vlolim

          cloud(ii,jj,ica%cur:fca%cur)%vratiohi = aero(ii,jj,ica%par:fca%par)%vratiohi
          cloud(ii,jj,icb%cur:fcb%cur)%vratiohi = aero(ii,jj,icb%par:fcb%par)%vratiohi

          cloud(ii,jj,ica%cur:fca%cur)%vratiolo = aero(ii,jj,ica%par:fca%par)%vratiolo
          cloud(ii,jj,icb%cur:fcb%cur)%vratiolo = aero(ii,jj,icb%par:fcb%par)%vratiolo

          cloud(ii,jj,ica%cur:fca%cur)%dmid = aero(ii,jj,ica%par:fca%par)%dmid
          cloud(ii,jj,icb%cur:fcb%cur)%dmid = aero(ii,jj,icb%par:fcb%par)%dmid

          ! Initialize the droplet diameter ("wet diameter") as the dry
          ! mid diameter of the nucleus to avoid problems later.
          cloud(ii,jj,:)%dwet = cloud(ii,jj,:)%dmid
          cloud(ii,jj,:)%core = pi6*cloud(ii,jj,:)%dmid**3

          ! Initialize the volume and number concentrations for clouds.
          ! First "real" values are only obtained upon the first calculation
          ! of the cloud droplet activation.
          DO cc = 1,8
             cloud(ii,jj,:)%volc(cc) = 0.
          END DO
          cloud(ii,jj,:)%numc = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the precipitation properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the *wet* radius
          ! ---------------------------------------------------------------------------------------
          precp(ii,jj,:)%vhilim = pi6*tmphilim(:)**3
          precp(ii,jj,:)%vlolim = pi6*tmplolim(:)**3
          precp(ii,jj,:)%dmid = ( (precp(ii,jj,:)%vlolim + precp(ii,jj,:)%vhilim) / (2.*pi6) )**(1./3.)
          precp(ii,jj,:)%vratiohi = precp(ii,jj,:)%vhilim / ( pi6*precp(ii,jj,:)%dmid**3 )
          precp(ii,jj,:)%vratiolo = precp(ii,jj,:)%vlolim / ( pi6*precp(ii,jj,:)%dmid**3 )

          ! Initialize the wet diameter as the bin mid diameter
          precp(ii,jj,:)%dwet = precp(ii,jj,:)%dmid
          precp(ii,jj,:)%core = pi6*precp(ii,jj,:)%dmid**3

          DO cc = 1,8
             precp(ii,jj,:)%volc(cc) = 0.
          END DO
          precp(ii,jj,:)%numc = 0.

       END DO
    END DO

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(cloudbins(ncld))
    cloudbins(:) = (cloud(1,1,:)%vlolim/pi6)**(1./3.)
    cloudbins = 0.5*cloudbins ! To radius
    ALLOCATE(precpbins(nprc))
    precpbins(:) = (precp(1,1,:)%vlolim/pi6)**(1./3.)
    precpbins = 0.5*precpbins ! To radius

  END SUBROUTINE set_cloudbins

  !--------------------------------------------------------------------------
  !
  ! *******************************
  ! SUBROUTINE set_icebins
  ! *******************************
  !
  ! Setup of ice size bins similar to the subroutines *set_sizebins* & *set_cloudbins*
  !
  ! Jaakko Ahola (FMI) 2015
  !
  !---------------------------------------------------------------------------
  SUBROUTINE set_icebins()
    USE mo_submctl, ONLY : pi6,             &
                               iia,iib,         &
                               fia,fib,         &
                               isa,fsa,         &
                               nice,nsnw,        &
                               in2a,fn2a,       &
                               fn2b,       &
                               icebins,         &
                               snowbins
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                                 ice,snow,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc,nba,nbb

    REAL :: tmplolim(7), tmphilim(7)

    ! Helper arrays to set up snow size bins
    tmplolim = (/50.,55.,65., 100.,200.,500., 1000./)*1.e-6
    tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6

    ! Number of ice bins in regime a (soluble nuclei)
    nba = fn2a-in2a+1
    ! Number of cloud bins in regime b (insoluble nuclei)
    nbb = nba

    ! Reset ice bin indices accordingly. The two components give the cloud regime index,
    ! and the aerosol bin index with which they are parallel
    iia%cur = 1;                       iia%par = in2a
    fia%cur = iia%cur + nba-1;  fia%par = iia%par + nba-1
    iib%cur = fia%cur + 1;             iib%par = fn2b - nbb + 1
    fib%cur = iib%cur + nbb-1;  fib%par = iib%par + nbb-1

    nice = fib%cur

    ! snow bins
    isa = 1; fsa = 7;
    nsnw = fsa

    ! ----------------------------------------
    ! Allocate ice arrays
    ! ----------------------------------------
    ALLOCATE(ice(kbdim,klev,nice), snow(kbdim,klev,nsnw))

    DO jj = 1,klev
       DO ii = 1,kbdim

          ! -------------------------------------------------
          ! Set iceproperties (parallel to aerosol bins)
          ! -------------------------------------------------
          ice(ii,jj,iia%cur:fia%cur)%vhilim = aero(ii,jj,iia%par:fia%par)%vhilim
          ice(ii,jj,iib%cur:fib%cur)%vhilim = aero(ii,jj,iib%par:fib%par)%vhilim

          ice(ii,jj,iia%cur:fia%cur)%vlolim = aero(ii,jj,iia%par:fia%par)%vlolim
          ice(ii,jj,iib%cur:fib%cur)%vlolim = aero(ii,jj,iib%par:fib%par)%vlolim

          ice(ii,jj,iia%cur:fia%cur)%vratiohi = aero(ii,jj,iia%par:fia%par)%vratiohi
          ice(ii,jj,iib%cur:fib%cur)%vratiohi = aero(ii,jj,iib%par:fib%par)%vratiohi

          ice(ii,jj,iia%cur:fia%cur)%vratiolo = aero(ii,jj,iia%par:fia%par)%vratiolo
          ice(ii,jj,iib%cur:fib%cur)%vratiolo = aero(ii,jj,iib%par:fib%par)%vratiolo

          ice(ii,jj,iia%cur:fia%cur)%dmid = aero(ii,jj,iia%par:fia%par)%dmid
          ice(ii,jj,iib%cur:fib%cur)%dmid = aero(ii,jj,iib%par:fib%par)%dmid

          ! Initialize the "wet" diameter as the dry mid diameter of the nucleus
          ice(ii,jj,:)%dwet = ice(ii,jj,:)%dmid
          ice(ii,jj,:)%core = pi6*ice(ii,jj,:)%dmid**3

          ! Initialize the volume and number concentrations for ice.
          DO cc = 1,8
             ice(ii,jj,:)%volc(cc) = 0.
          END DO
          ice(ii,jj,:)%numc = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the snow properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the "wet" radius
          ! ---------------------------------------------------------------------------------------

          snow(ii,jj,:)%vhilim = pi6*tmphilim(:)**3
          snow(ii,jj,:)%vlolim = pi6*tmplolim(:)**3
          snow(ii,jj,:)%dmid = ( (snow(ii,jj,:)%vlolim + snow(ii,jj,:)%vhilim) / (2.*pi6) )**(1./3.)
          snow(ii,jj,:)%vratiohi = snow(ii,jj,:)%vhilim / ( pi6*snow(ii,jj,:)%dmid**3 )
          snow(ii,jj,:)%vratiolo = snow(ii,jj,:)%vlolim / ( pi6*snow(ii,jj,:)%dmid**3 )

          ! Initialize the wet diameter as the bin mid diameter
          snow(ii,jj,:)%dwet = snow(ii,jj,:)%dmid
          snow(ii,jj,:)%core = pi6*snow(ii,jj,:)%dmid**3

          DO cc = 1,8
             snow(ii,jj,:)%volc(cc) = 0.
          END DO
          snow(ii,jj,:)%numc = 0.

       END DO
    END DO

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(icebins(nice))
    icebins(:) = (ice(1,1,:)%vlolim/pi6)**(1./3.)
    icebins = 0.5*icebins ! To radius
    ALLOCATE(snowbins(nsnw))
    snowbins(:) = (snow(1,1,:)%vlolim/pi6)**(1./3.)
    snowbins = 0.5*snowbins ! To radius

  END SUBROUTINE set_icebins



  !----------------------------------------------------------------------
  !
  ! *************************
  ! SUBROUTINE define_salsa
  ! *************************
  !
  ! Reads logical switches and aerosol/hydrometeor size bin definitions
  ! from a namelist.
  !
  ! Juha Tonttila (FMI) 2014
  !
  !----------------------------------------------------------------------
  SUBROUTINE define_salsa(level)

    USE mo_submctl, ONLY : nlcoag,                &
                               nlcgaa,nlcgcc,nlcgpp,  &
                               nlcgca,nlcgpa,nlcgpc,  &
                               nlcgia,nlcgic,nlcgii,  &
                               nlcgip,nlcgsa,nlcgsc,  &
                               nlcgsi,nlcgsp,nlcgss,  &
                               nlcnd,                 &
                               nlcndgas,              &
                               nlcndh2oae,nlcndh2ocl, &
                               nlcndh2oic,            &
                               nlauto,nlautosnow,     &
                               autoc_rain_zd0, autoc_rain_sigmag, &
                               autoc_snow_zd0, autoc_snow_sigmag, &
                               nlactiv,               &
                               nlactintst,            &
                               nlactbase,            &
                               nlicenucl,             &
                               fixinc,                &
                               ice_hom, ice_imm, ice_dep, &
                               icenucl_tstart,        &
                               nlicmelt,              &
                               stat_b_bins, stat_micro_ts, &
                               nbin,reglim,           &
                               nice,nsnw,             &
                               nspec,listspec,        &
                               volDistA, volDistB,    &
                               nf2a, isdtyp,          &
                               sigmag,dpg,n,          &
                               rhlim

    IMPLICIT NONE

    INTEGER, INTENT(in) :: level

    NAMELIST /salsa/  &
         nlcoag,      & ! Coagulation master switch
         nlcgaa,      & ! Coagulation between aerosols
         nlcgcc,      & ! Collision-coalescence between cloud droplets
         nlcgpp,      & ! Collisions between rain drops
         nlcgca,      & ! Cloud collection of aerosols
         nlcgpa,      & ! Collection of aerosols by precip
         nlcgpc,      & ! Collection of cloud droplets by rain
         nlcgia,      & ! Ice collection of aerosols
         nlcgic,      & ! Collection of cloud droplets by ice particles
         nlcgii,      & ! Collision-coalescence between ice particles
         nlcgip,      & ! Collection of precipitation by ice particles
         nlcgsa,      & ! Collection of aerosols by snow
         nlcgsc,      & ! Collection of cloud droplets by snow
         nlcgsi,      & ! Collection of ice by snow
         nlcgsp,      & ! Collection of precipitation by snow
         nlcgss,      & ! Collision-coalescence between snow particles
         nlcnd,       & ! Switch for condensation subroutine
         nlcndgas,    & ! Condensation of precursor gases
         nlicenucl,   & ! Ice nucleation master switch
         fixinc,      & ! Constant ice number concentration (fixinc > 0 #/kg) is maintained by converting cloud droplets to ice
         ice_hom,     & ! If fixinc is not set or it is not positive, ice nucleation can be modelled based on homogeneous, ...
         ice_imm,     & ! immersion and/or ...
         ice_dep,     & ! deposition freezing mechanisms
         icenucl_tstart, & ! Start time (s) for ice formation
         nlicmelt,    & ! Switch for ice'n'snow melting
         nbin,        & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         nice,        & ! number of ice bins
         nsnw,        & ! number of snow bins
         nlcndh2ocl,    & ! Condensation of water vapour on clouds (drizzle)
         nlcndh2oic,    & ! Condensation of water vapour on ice particles ! ice'n'snow
         nlcndh2oae,    & ! Condensation of water vapour on aerosols (FALSE -> equilibrium calc.)
         nlauto,        & ! Switch for autoconversion of cloud droplets to drizzle and rain
         nlautosnow,    & ! Switch for autoconversion of ice particles to snowing
         autoc_rain_zd0, autoc_rain_sigmag, & ! Cloud to rain autoconversion parameters
         autoc_snow_zd0, autoc_snow_sigmag, & ! Ice to snow autoconversion parameters
         nlactiv,       & ! Master switch for cloud droplet activation
         nlactbase,     & ! Switch for parameterized cloud base activation
         nlactintst,    & ! Switch for interstitial activation based on particle growth and host model S
         stat_b_bins,   & ! Save statistics about SALSA b-bins
         stat_micro_ts, & ! Save statistics (*.ts) about microphysical process rates

         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 - vertical profile, read from file
         reglim,        & ! Low/high diameter limits of the 2 aerosol size regimes (1d table with length 4)
         nbin,          & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         nspec,         & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
                          ! Must be an array of length 7, with empty strings for unused stuff.
         volDistA,      & ! Initial relative contribution [0-1] of each species to particle volume in a-bins. Must be
                          ! an array of length 7, with zero for unused species.
         volDistB,      & ! Same as above but for b-bins
         nf2a,          & ! Number fraction of particles allocated to a-bins in regime 2. b-bins will get 1-nf2a

         rhlim,         & ! Upper limit RH/100 for sals during initialization and spinup 

         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n                ! Number concentration for the 7 initial lognormal modes



    OPEN(11,STATUS='old',FILE='NAMELIST')
    READ(11,NML=salsa)
    CLOSE(11)


    ! if thermodynamical level is less than 5, set all ice process switches to false
    IF(level < 5) THEN
          nlcgia      = .false.
          nlcgic      = .false.
          nlcgii      = .false.
          nlcgip      = .false.
          nlcgsa      = .false.
          nlcgsc      = .false.
          nlcgsi      = .false.
          nlcgsp      = .false.
          nlcgss      = .false.

          nlcndh2oic  = .false.

          nlautosnow  = .false.

          nlicenucl   = .false.
          nlicmelt    = .false.
    END IF !level

  END SUBROUTINE define_salsa



  !-------------------------------------------------------------------------------
  !
  ! *****************************
  ! SUBROUTINE salsa_initialize
  ! *****************************
  !
  ! SALSA initializations. Modified and rewritten for more dynamic control
  ! and LES implementation.
  !
  ! Juha Tonttila (FMI) 2014
  !
  !-------------------------------------------------------------------------------
  SUBROUTINE salsa_initialize()

    !
    !-------------------------------------------------------------------------------

    USE mo_submctl, ONLY : nbin,     &
                               in1a,fn1a,in2a,fn2a,in2b,fn2b,  &
                               nbins, &
                               massacc

    IMPLICIT NONE

    ! Remember to call 'define_salsa' for namelist paramers before calling this subroutine!

    ! --1) Set derived indices

    in1a = 1
    in2a = in1a + nbin(1)
    in2b = in2a + nbin(2)

    fn1a = in2a - 1
    fn2a = fn1a + nbin(2)
    fn2b = fn2a + nbin(2)

    nbins = fn2b

    ! --2) Allocate arrays

    ALLOCATE(massacc(nbins))

    massacc = 1.


    ! -- Aerosol tracers are allocated in *set_sizebins*
    ! -- Hydrometeor tracer in *set_cloudbins*

    ! --3) Call other initialization routines
    CALL set_sizebins()

    CALL set_cloudbins()

    CALL set_icebins()

  END SUBROUTINE salsa_initialize


END MODULE mo_salsa_init
