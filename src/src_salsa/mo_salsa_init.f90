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
         maxnreg,     & ! maximum number of subregimes
         in1a, fn1a,  & ! size regime bin indices: 1a
         in2a, fn2a,  & !     - " -       2a
         in2b, fn2b,  & !     - " -       2b
         nbins,  &
         maxnspec,  &
         aerobins

    USE mo_salsa_driver, ONLY : &
         kbdim,klev,aero

    IMPLICIT NONE

    !-- local variables ----------
    INTEGER :: ii, jj,cc,dd,vv ! loop indices
    REAL ::  ratio ! ratio of regime upper and lower diameter

    ! -------------------------
    ! Allocate aerosol tracers
    !--------------------------
    ALLOCATE(aero(kbdim,klev,nbins))

    ! Juha: Corrected some bugs with the section volume definitions
    !-------------------------------------------------------------------------------

    DO jj = 1,klev
       DO ii = 1,kbdim

          !-- 1) aerosol size regimes 1 and 2
          !  - minimum & maximum *dry* volumes [fxm]
          !  - bin mid *dry* diameter [m]
          dd = 0
          DO vv = 1,maxnreg
              IF (nbin(vv)<=0) EXIT
              ratio = reglim(vv+1)/reglim(vv)

              DO cc = 1,nbin(vv)
                 dd = dd+1
                 aero(ii,jj,dd)%vlolim = pi6*(reglim(vv)*ratio**(real(cc-1)/nbin(vv)))**3
                 aero(ii,jj,dd)%vhilim = pi6*(reglim(vv)*ratio**(real(cc)/nbin(vv)))**3
                 aero(ii,jj,dd)%dmid = ( (aero(ii,jj,dd)%vhilim + aero(ii,jj,dd)%vlolim) /(2.*pi6) )**(1./3.)
                 aero(ii,jj,dd)%vratiohi = aero(ii,jj,dd)%vhilim/(pi6*aero(ii,jj,dd)%dmid**3)
                 aero(ii,jj,dd)%vratiolo = aero(ii,jj,dd)%vlolim/(pi6*aero(ii,jj,dd)%dmid**3)
              END DO
          END DO
          IF (vv<3 .OR. dd/=fn2a) THEN
             WRITE(*,*) cc, nbin(1:cc), dd, fn2a
             WRITE(*,*) 'Error in given aerosol size bin definitions!'
             STOP
          ENDIF

          !-- 2) same values for subregime 2b
          aero(ii,jj,in2b:fn2b)%vlolim = aero(ii,jj,in2a:fn2a)%vlolim
          aero(ii,jj,in2b:fn2b)%vhilim = aero(ii,jj,in2a:fn2a)%vhilim
          aero(ii,jj,in2b:fn2b)%dmid = aero(ii,jj,in2a:fn2a)%dmid
          aero(ii,jj,in2b:fn2b)%vratiohi = aero(ii,jj,in2a:fn2a)%vratiohi
          aero(ii,jj,in2b:fn2b)%vratiolo = aero(ii,jj,in2a:fn2a)%vratiolo

          ! Initialize the wet diameter with the bin dry diameter to avoid numerical proplems later
          aero(ii,jj,:)%dwet = aero(ii,jj,:)%dmid

          ! Set volume and number concentrations to zero
          aero(ii,jj,:)%numc = 0.
          DO vv=1,maxnspec
             aero(ii,jj,:)%volc(vv) = 0.
          END DO

       END DO !ii
    END DO !!jj

    ! Save 1a and 2a bin limits to be delivered e.g. to host model if needed
    ALLOCATE(aerobins(fn2a+1))
    DO cc = 1,fn2a
       aerobins(cc) = (aero(1,1,cc)%vlolim/pi6)**(1./3.)
    END DO
    aerobins(fn2a+1) = (aero(1,1,fn2a)%vhilim/pi6)**(1./3.)
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
    USE mo_submctl, ONLY : pi6,maxnspec,                &
                               in2a,fn2a,in2b,fn2b,     &
                               inp2a,fnp2a,inp2b,fnp2b, &
                               ncld,nprc,               &
                               precpbins,rainbinlim,    &
                               cldbinlim,nout_cld
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                               cloud,precp,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc,nba,nbb

    ! Cloud (and ice) bins are parallel with the 2a and 2b aerosol bins

    ! Number of cloud bins in regime a
    nba = fn2a-in2a+1
    ! Number of cloud bins in regime b
    nbb = fn2b-in2b+1

    ! Indices for the 2a and 2b cloud (and ice) bins
    inp2a = 1;       fnp2a = nba       ! 2a
    inp2b = fnp2a+1; fnp2b = nba + nbb ! 2b

    ncld = nba + nbb

    ! Rain/drizzle bins
    IF (rainbinlim(1)<0.) THEN
        ! Use the default
        rainbinlim(1:8)=(/50.,55.,65.,100.,200.,500.,1000.,2000./)*1.e-6
        nprc=7
    ELSEIF (rainbinlim(2)<0. .OR. rainbinlim(1)>rainbinlim(2) .OR. rainbinlim(1)<0.1 .OR. rainbinlim(1)>1000.) THEN
        ! Bins should start from about 20-80 microns
        WRITE(*,*) 'Bad input rain bin limits!'
        STOP
    ELSE
        ii=2
        DO WHILE (rainbinlim(ii)>0.)
            IF (rainbinlim(ii-1)>rainbinlim(ii)) THEN
                WRITE(*,*) 'Non-monotonic input rain bin limits!'
                STOP
            ENDIF
            ii=ii+1
        END DO
        rainbinlim(1:ii-1)=rainbinlim(1:ii-1)*1.e-6 ! from microns to meters
        nprc=ii-2
    ENDIF

    ! ----------------------------------------
    ! Allocate cloud and precipitation arrays
    ! ----------------------------------------
    ALLOCATE(cloud(kbdim,klev,ncld), precp(kbdim,klev,nprc))

    DO jj = 1,klev
       DO ii = 1,kbdim

          ! -------------------------------------------------
          ! Set cloud properties (parallel to aerosol bins)
          ! -------------------------------------------------
          cloud(ii,jj,1:)%vhilim = aero(ii,jj,in2a:)%vhilim
          cloud(ii,jj,1:)%vlolim = aero(ii,jj,in2a:)%vlolim
          cloud(ii,jj,1:)%vratiohi = aero(ii,jj,in2a:)%vratiohi
          cloud(ii,jj,1:)%vratiolo = aero(ii,jj,in2a:)%vratiolo
          cloud(ii,jj,1:)%dmid = aero(ii,jj,in2a:)%dmid

          ! Initialize the droplet diameter ("wet diameter") as the dry
          ! mid diameter of the nucleus to avoid problems later.
          cloud(ii,jj,:)%dwet = cloud(ii,jj,:)%dmid

          ! Initialize the volume and number concentrations for clouds.
          ! First "real" values are only obtained upon the first calculation
          ! of the cloud droplet activation.
          DO cc = 1,maxnspec
             cloud(ii,jj,:)%volc(cc) = 0.
          END DO
          cloud(ii,jj,:)%numc = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the precipitation properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the *wet* radius
          ! ---------------------------------------------------------------------------------------
          precp(ii,jj,:)%vhilim = pi6*rainbinlim(2:nprc+1)**3
          precp(ii,jj,:)%vlolim = pi6*rainbinlim(1:nprc)**3
          precp(ii,jj,:)%dmid = ( (precp(ii,jj,:)%vlolim + precp(ii,jj,:)%vhilim) / (2.*pi6) )**(1./3.)
          precp(ii,jj,:)%vratiohi = precp(ii,jj,:)%vhilim / ( pi6*precp(ii,jj,:)%dmid**3 )
          precp(ii,jj,:)%vratiolo = precp(ii,jj,:)%vlolim / ( pi6*precp(ii,jj,:)%dmid**3 )

          ! Initialize the wet diameter as the bin mid diameter
          precp(ii,jj,:)%dwet = precp(ii,jj,:)%dmid

          DO cc = 1,maxnspec
             precp(ii,jj,:)%volc(cc) = 0.
          END DO
          precp(ii,jj,:)%numc = 0.

       END DO
    END DO

    ! Save precipitation bin limits to be delivered e.g. to host model if needed
    ALLOCATE(precpbins(nprc+1))
    precpbins(1:nprc) = (precp(1,1,1:nprc)%vlolim/pi6)**(1./3.)
    precpbins(nprc+1) = (precp(1,1,nprc)%vhilim/pi6)**(1./3.)
    precpbins = 0.5*precpbins ! To radius

    ! Cloud bin limits for outputs
    IF (cldbinlim(1)>=0.) THEN
        ii=2
        DO WHILE (cldbinlim(ii)>=0.)
            IF (cldbinlim(ii-1)>cldbinlim(ii)) THEN
                WRITE(*,*) 'Non-monotonic input cloud bin limits!'
                STOP
            ENDIF
            ii=ii+1
        END DO
        cldbinlim(1:ii-1)=cldbinlim(1:ii-1)*0.5e-6 ! from microns to meters and to radius
        nout_cld=ii-2
    ENDIF

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
    USE mo_submctl, ONLY : pi6,maxnspec,             &
                               in2a,ncld,nice,nsnw,  &
                               snowbins, snowbinlim, &
                               icebinlim, nout_ice
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                               ice,snow,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc

    ! Ice bins are identical with the cloud bins
    nice = ncld

    ! snow bins
    IF (snowbinlim(1)<0.) THEN
        ! Use the default
        snowbinlim(1:8)=(/50.,55.,65.,100.,200.,500.,1000.,2000./)*1.e-6
        nsnw=7
    ELSEIF (snowbinlim(2)<0. .OR. snowbinlim(1)>snowbinlim(2) .OR. snowbinlim(1)<0.1 .OR. snowbinlim(1)>1000.) THEN
        ! Bins should start from about 20-80 microns
        WRITE(*,*) 'Bad input snow bin limits!'
        STOP
    ELSE
        ii=2
        DO WHILE (snowbinlim(ii)>0.)
            IF (snowbinlim(ii-1)>snowbinlim(ii)) THEN
                WRITE(*,*) 'Non-monotonic input snow bin limits!'
                STOP
            ENDIF
            ii=ii+1
        END DO
        snowbinlim(1:ii-1)=snowbinlim(1:ii-1)*1.e-6 ! from microns to meters
        nsnw=ii-2
    ENDIF

    ! ----------------------------------------
    ! Allocate ice arrays
    ! ----------------------------------------
    ALLOCATE(ice(kbdim,klev,nice), snow(kbdim,klev,nsnw))

    DO jj = 1,klev
       DO ii = 1,kbdim

          ! -------------------------------------------------
          ! Set ice properties (parallel to aerosol bins)
          ! -------------------------------------------------
          ice(ii,jj,1:)%vhilim = aero(ii,jj,in2a:)%vhilim
          ice(ii,jj,1:)%vlolim = aero(ii,jj,in2a:)%vlolim
          ice(ii,jj,1:)%vratiohi = aero(ii,jj,in2a:)%vratiohi
          ice(ii,jj,1:)%vratiolo = aero(ii,jj,in2a:)%vratiolo
          ice(ii,jj,1:)%dmid = aero(ii,jj,in2a:)%dmid

          ! Initialize the "wet" diameter as the dry mid diameter of the nucleus
          ice(ii,jj,:)%dwet = ice(ii,jj,:)%dmid

          ! Initialize the volume and number concentrations for ice.
          DO cc = 1,maxnspec
             ice(ii,jj,:)%volc(cc) = 0.
          END DO
          ice(ii,jj,:)%numc = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the snow properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the "wet" radius
          ! ---------------------------------------------------------------------------------------

          snow(ii,jj,:)%vhilim = pi6*snowbinlim(2:nsnw+1)**3
          snow(ii,jj,:)%vlolim = pi6*snowbinlim(1:nsnw)**3
          snow(ii,jj,:)%dmid = ( (snow(ii,jj,:)%vlolim + snow(ii,jj,:)%vhilim) / (2.*pi6) )**(1./3.)
          snow(ii,jj,:)%vratiohi = snow(ii,jj,:)%vhilim / ( pi6*snow(ii,jj,:)%dmid**3 )
          snow(ii,jj,:)%vratiolo = snow(ii,jj,:)%vlolim / ( pi6*snow(ii,jj,:)%dmid**3 )

          ! Initialize the wet diameter as the bin mid diameter
          snow(ii,jj,:)%dwet = snow(ii,jj,:)%dmid

          DO cc = 1,maxnspec
             snow(ii,jj,:)%volc(cc) = 0.
          END DO
          snow(ii,jj,:)%numc = 0.

       END DO
    END DO

    ! Save snow bin limits to be delivered e.g. to host model if needed
    ALLOCATE(snowbins(nsnw+1))
    snowbins(1:nsnw) = (snow(1,1,1:nsnw)%vlolim/pi6)**(1./3.)
    snowbins(nsnw+1) = (snow(1,1,nsnw)%vhilim/pi6)**(1./3.)
    snowbins = 0.5*snowbins ! To radius

    ! Ice bin limits for outputs
    IF (icebinlim(1)>=0.) THEN
        ii=2
        DO WHILE (icebinlim(ii)>=0.)
            IF (icebinlim(ii-1)>icebinlim(ii)) THEN
                WRITE(*,*) 'Non-monotonic input ice bin limits!'
                STOP
            ENDIF
            ii=ii+1
        END DO
        icebinlim(1:ii-1)=icebinlim(1:ii-1)*0.5e-6 ! from microns to meters and to radius
        nout_ice=ii-2
    ENDIF

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
                               eddy_dis_rt,           &
                               nlcgrain,              &
                               nlcnd,                 &
                               nlcndgas,              &
                               nlcndh2oae,nlcndh2ocl, &
                               nlcndh2oic,            &
                               rhlim,                 &
                               nlauto,nlautosnow,     &
                               auto_sb,               &
                               autoc_rain_zd0, autoc_rain_sigmag, &
                               autoc_snow_zd0, autoc_snow_sigmag, &
                               nlactiv,               &
                               nlactintst,            &
                               nlactbase,            &
                               nlactprc,             &
                               nlicenucl,             &
                               fixinc, ice_source_opt,&
                               ice_hom, ice_imm, ice_dep, &
                               icenucl_tstart,        &
                               ice_target_opt,        &
                               nlicmelt,              &
                               rainbinlim,            &
                               snowbinlim,            &
                               nbin,reglim,           &
                               nspec,listspec,        &
                               volDistA, volDistB,    &
                               salsa1a_SO4_OC,        &
                               nf2a, isdtyp,          &
                               sigmag,dpg,n,          &
                               msu, disssu, rhosu,    &
                               mno, dissno, rhono,    &
                               mnh, dissnh, rhonh,    &
                               moc, dissoc, rhooc,    &
                               mbc, dissbc, rhobc,    &
                               mss, dissss, rhoss,    &
                               mdu, dissdu, rhodu,    &
                               mwa, disswa, rhowa,    &
                               conc_h2so4, conc_ocnv, &
                               nvbs_setup, laqsoa,    &
                               zvbs_k_OH, zvbs_Eact_p_OH, &
                               ox_prescribed,         &
                               conc_oh, conc_o3, conc_no3, &
                               zdayfac_oh, zdayfac_o3, znightfac_no3, ox_conc_flag, &
                               zphotofac_aqsoa, aqsoa_photo_flag, &
                               conc_voc, conc_vbsg, conc_aqsoag, &
                               cldbinlim, icebinlim

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
         nlcgrain,    & ! Rain formation based on cloud-cloud collisions
         nlcgia,      & ! Ice collection of aerosols
         nlcgic,      & ! Collection of cloud droplets by ice particles
         nlcgii,      & ! Collision-coalescence between ice particles
         nlcgip,      & ! Collection of precipitation by ice particles
         nlcgsa,      & ! Collection of aerosols by snow
         nlcgsc,      & ! Collection of cloud droplets by snow
         nlcgsi,      & ! Collection of ice by snow
         nlcgsp,      & ! Collection of precipitation by snow
         nlcgss,      & ! Collision-coalescence between snow particles
         eddy_dis_rt, & ! Eddy dissipation rate, negative means take from LES

         nlcnd,       & ! Condensation master switch
         nlcndgas,    & ! Condensation of H2SO4 and organic vapors
         nlcndh2ocl,  & ! Condensation of water vapour on clouds and drizzle
         nlcndh2oic,  & ! Condensation of water vapour on ice and snow particles
         nlcndh2oae,  & ! Condensation of water vapour on aerosols (FALSE -> equilibrium calc.)
         rhlim,       & ! Upper limit RH/100 during initialization and spinup

         nlauto,        & ! Switch for autoconversion of cloud droplets to rain
         nlautosnow,    & ! Switch for autoconversion of ice particles to snow
         auto_sb,       & ! Use the Seifert & Beheng (2001) autoconversion method
         autoc_rain_zd0, autoc_rain_sigmag, & ! Cloud to rain autoconversion parameters
         autoc_snow_zd0, autoc_snow_sigmag, & ! Ice to snow autoconversion parameters

         nlactiv,       & ! Master switch for cloud droplet activation
         nlactbase,     & ! Switch for parameterized cloud base activation
         nlactintst,    & ! Switch for interstitial activation based on particle growth and host model S
         nlactprc,      & ! Activated aerosol described with rain bins

         nlicenucl,     & ! Ice nucleation master switch
         fixinc,        & ! Constant ice number concentration (fixinc > 0 #/kg) is maintained by converting cloud droplets to ice
         ice_source_opt,& ! Cloud freezing order: >0: start from the largest bin, 0: all bins evenly, <0: start from the smallest bin
         ice_hom,       & ! If fixinc is not set or it is not positive, ice nucleation can be modelled based on homogeneous, ...
         ice_imm,       & ! immersion and/or ...
         ice_dep,       & ! deposition freezing mechanisms
         icenucl_tstart,& ! Start time (s) for ice formation
         ice_target_opt,& ! Where to put new ice/snow: <0: parallel ice bin, 0: find matching snow bin, >0 snow bin specified by ice_target_opt
         nlicmelt,      & ! Switch for ice and snow melting

         rainbinlim,    & ! Rain bin limits (microns)
         snowbinlim,    & ! Snow bin limits (microns)
         nbin,          & ! Number of bins used for the 1a and 2a aerosol size regimes (1d table with length 2)
         reglim,        & ! Low/high diameter limits for the 1a and 2a aerosol size regimes (1d table with length 4)
         nspec,         & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 - vertical profile, read from file
         volDistA,      & ! Initial relative contribution [0-1] of each species to particle volume in a-bins.
         volDistB,      & ! Same as above but for b-bins
         salsa1a_SO4_OC,& ! Limit 1a composition to OC and/or SO4
         nf2a,          & ! Number fraction of particles allocated to a-bins in regime 2. b-bins will get 1-nf2a
         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n,             & ! Number concentration for the 7 initial lognormal modes

         msu, disssu, rhosu, & ! Physical properties of the species; sulphate
         mno, dissno, rhono, & ! HNO3
         mnh, dissnh, rhonh, & ! NH3
         moc, dissoc, rhooc, & ! organic carbon
         mbc, dissbc, rhobc, & ! black carbon
         mss, dissss, rhoss, & ! sea salt (NaCl)
         mdu, dissdu, rhodu, & ! mineral dust
         mwa, disswa, rhowa, & ! water

         conc_h2so4,    & ! Vapor phase concentration for sulfuric acid (#/kg)
         conc_ocnv,     & ! -||- non-volatile organics

         nvbs_setup,    & ! Detailed secondary organic aerosol formation (VOC+oxidant => BVS(g) <=> VBS(s))
         laqsoa,        & ! Additional aqSOA formation
         zvbs_k_OH, zvbs_Eact_p_OH, & ! Rate coefficient for VBS(g) aging
         ox_prescribed, & ! Oxidant concentrations can be fixed (diagnostic parameter)
         conc_oh, conc_o3, conc_no3, & ! Initial oxidant concentrations (number mixing ratios)
         zdayfac_oh, zdayfac_o3, znightfac_no3, ox_conc_flag, & ! Scaling factors for oxidant concentrations
         zphotofac_aqsoa, aqsoa_photo_flag, & ! Scaling factor for aqSOA photodissociation
         conc_voc, conc_vbsg, conc_aqsoag, & ! Arrays for initial VOC(g), VBS(g) and aqSOA(g) concentrations (mass mixing ratios)

         cldbinlim,     & ! Output cloud bin diameter limits (microns)
         icebinlim        ! Output ice bin diameter limits (microns)


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

    USE mo_submctl, ONLY : nbin,in1a,fn1a,in2a,fn2a,in2b,fn2b,nbins,massacc, &
                           zspec,dens,diss,mws,ih2o,iso,ioc,ibc,idu,iss,inh,ino,nspec,listspec, &
                           rhowa,disswa,mwa,rhosu,disssu,msu,rhooc,dissoc,moc,rhobc,dissbc,mbc, &
                           rhodu,dissdu,mdu,rhoss,dissss,mss,rhono,dissno,mno,rhonh,dissnh,mnh, &
                           nlcndgas,ngases,zgas,mws_gas, &
                           conc_h2so4,conc_ocnv,part_h2so4,part_ocnv,isog,iocg, &
                           nvbs_setup,laqsoa, model_lat, start_doy
    USE step, ONLY : cntlat, strtim
    USE mo_vbs_init, ONLY : init_vbs
    IMPLICIT NONE
    INTEGER :: ss, nvbs

    ! Remember to call 'define_salsa' for namelist paramers before calling this subroutine!

    ! --1) Set derived indices

    in1a = 1
    in2a = in1a + nbin(1)
    in2b = in2a + SUM(nbin(2:))

    fn1a = in2a - 1
    fn2a = fn1a + SUM(nbin(2:))
    fn2b = fn2a + SUM(nbin(2:))

    nbins = fn2b

    ! --2) Allocate arrays

    ALLOCATE(massacc(nbins))

    massacc = 1.

    ! --3) Call other initialization routines
    CALL set_sizebins()

    CALL set_cloudbins()

    CALL set_icebins()

    ! --4) Set the input aerosol species: name, density, dissociation constant, MW and index

    ! Water is always the first species!
    ss=1
    zspec(ss)='H2O'
    dens(ss)=rhowa
    diss(ss)=disswa
    mws(ss)=mwa
    ih2o=ss
    ! .. then normal aerosol species
    nvbs=0
    DO ss=2,nspec+1
         zspec(ss)=listspec(ss-1) ! Inputs do not include water
         SELECT CASE(listspec(ss-1))
            CASE('SO4')
                dens(ss)=rhosu
                diss(ss)=disssu
                mws(ss)=msu
                iso=ss
            CASE('OC')
                dens(ss)=rhooc
                diss(ss)=dissoc
                mws(ss)=moc
                ioc=ss
            CASE('BC')
                dens(ss)=rhobc
                diss(ss)=dissbc
                mws(ss)=mbc
                ibc=ss
            CASE('DU')
                dens(ss)=rhodu
                diss(ss)=dissdu
                mws(ss)=mdu
                idu=ss
            CASE('SS')
                dens(ss)=rhoss
                diss(ss)=dissss
                mws(ss)=mss
                iss=ss
            CASE('NO')
                dens(ss)=rhono
                diss(ss)=dissno
                mws(ss)=mno
                ino=ss
            CASE('NH')
                dens(ss)=rhonh
                diss(ss)=dissnh
                mws(ss)=mnh
                inh=ss
            CASE('VB1','VB2','VB3','VB4','VB5','VB6','VB7','VB8','VB9','AQ1','AQ2')
                ! Volatility Basis Set (VBS) bins and aqSOA species: their properties will be specified
                ! later in the VBS setup. However, these species should be listed here when initial
                ! aerosol-phase volume fractions are specified.
                nvbs=nvbs + 1
            CASE DEFAULT
                WRITE(*,*) 'Unkown species: '//TRIM(zspec(ss))
                STOP
        END SELECT
    END DO
    ! Volatility bin input concentrations are separate from the default aerosol input,
    ! so do not change the counter!
    nspec = nspec - nvbs

    ! --5) Gas phase chemistry
    ngases = 0
    IF (nlcndgas) THEN
        ! Simple SO2 and non-volatile organics
        IF (conc_h2so4>=0.) THEN
            IF (iso<=0) STOP 'Error: sulfate partitioning without aerosol SO4!'
            ngases = ngases + 1
            part_h2so4 = .TRUE.
            zgas(ngases)='SO2'
            isog = ngases
            mws_gas(ngases)=msu
        ENDIF
        IF (conc_ocnv>=0.) THEN
            IF (ioc<=0) STOP 'Error: organic vapor partitioning without aerosol OC!'
            ngases = ngases + 1
            part_ocnv = .TRUE.
            zgas(ngases)='NOA'
            iocg = ngases
            mws_gas(ngases)=moc
        ENDIF

        ! Detailed SOA formation includes VOC(g) -> VBS(g) <=> VBS(s,aq) and optionally also aqSOA
        IF (nvbs_setup>=0) THEN
            ! VBS setup
            !   Export: nvbs_setup, laqsoa
            !   Update (gas): nvocs, nvbs, naqsoa, ngases, ngases_diag, nspec, mws_gas, zgas, id_oh, id_no3, id_o3
            !   Update (aerosol): dens, diss, mws, zspec
            ! Additional VBS parameters
            start_doy=strtim ! Start time as decimal day of year
            model_lat=cntlat ! Center latitude (degrees)

            CALL init_vbs(nvbs_setup, laqsoa)
        ELSE
            laqsoa = .false.
            nvbs_setup = -1
        ENDIF
    ENDIF

  END SUBROUTINE salsa_initialize


END MODULE mo_salsa_init
