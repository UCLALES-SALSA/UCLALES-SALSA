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
         nbin2, nbin3,& !
         nbins,  &
         aerobins

    USE mo_salsa_driver, ONLY : &
         kbdim,klev,  &
         aero

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

          ! Set volume and number concentrations to zero
          aero(ii,jj,:)%numc = 0.
          aero(ii,jj,:)%core = 0.
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
                               ncldbin,         &
                               dmincld,         &
                               in1a,fn2a,       &
                               in2b,fn2b,       &
                               cloudbins,       &
                               precpbins
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                                 cloud,precp,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc,bb

    ! For cloud bins
    REAL :: zvoldiff(fn2a)
    INTEGER :: zindex(fn2a)
    INTEGER :: imin, nba,nbb
    INTEGER :: armin(1)
    LOGICAL :: l_min(fn2a)

    REAL :: tmplolim(7), tmphilim(7)

    ! Helper arrays to set up precipitation size bins
    tmplolim = (/50.,55.,65.,100.,200.,500.,1000./)*1.e-6
    tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6

    ! For setting cloud droplet bins
    DO ii = in1a,fn2a
       zindex(ii) = ii
    END DO

    ! Determine the smallest full aerosol bin with which the smallest cloud droplet bin will coincide
    zvoldiff(in1a:fn2a) = ABS(aero(1,1,in1a:fn2a)%vlolim - pi6*dmincld**3)
    l_min = ( zvoldiff(in1a:fn2a) == MINVAL(zvoldiff(in1a:fn2a)) )
    IF (COUNT(l_min) /= 1) STOP 'Initialization error! Try changing the cloud droplet bin lower bound slightly'
    armin = PACK(zindex,l_min)
    imin = armin(1)
    IF (aero(1,1,imin)%vlolim < pi6*dmincld**3) imin=imin+1
    IF (imin < in1a .OR. imin > fn2a) STOP 'Invalid first cloud bin'

    ! Number of cloud bins in regime a (soluble nuclei)
    nba = fn2a-imin+1
    ! Number of cloud bins in regime b (insoluble nuclei)
    IF (aero(1,1,imin)%dmid > aero(1,1,in2b)%dmid) THEN
       ! All cloud bins are within regime 2
       nbb = nba
    ELSE
       ! The smallest regime a cloud bins go to regime 1
       nbb = fn2b-fn2a
    END IF

    ncldbin(1) = nba
    ncldbin(2) = nbb

    ! Reset cloud bin indices accordingly. The two components give the cloud regime index,
    ! and the aerosol bin index with which they are parallel
    ica%cur = 1;                      ica%par = imin
    fca%cur = ica%cur + ncldbin(1)-1; fca%par = ica%par + ncldbin(1)-1
    icb%cur = fca%cur + 1;            icb%par = fn2b - ncldbin(2) + 1
    fcb%cur = icb%cur + ncldbin(2)-1; fcb%par = icb%par + ncldbin(2)-1

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
          cloud(ii,jj,ica%cur:fcb%cur)%dwet = cloud(ii,jj,ica%cur:fcb%cur)%dmid

          ! Initialize the volume and number concentrations for clouds.
          ! First "real" values are only obtained upon the first calculation
          ! of the cloud droplet activation.
          DO cc = 1,8
             cloud(ii,jj,ica%cur:fcb%cur)%volc(cc) = 0.
          END DO
          cloud(ii,jj,ica%cur:fcb%cur)%numc = 0.
          cloud(ii,jj,ica%cur:fcb%cur)%core = 0.
          cloud(ii,jj,ica%cur:fcb%cur)%veqh2o = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the precipitation properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the *wet* radius
          ! ---------------------------------------------------------------------------------------
          DO bb = ira,fra
             precp(ii,jj,bb)%vhilim = pi6*tmphilim(bb)**3
             precp(ii,jj,bb)%vlolim = pi6*tmplolim(bb)**3
             precp(ii,jj,bb)%dmid = ( (precp(ii,jj,bb)%vlolim + precp(ii,jj,bb)%vhilim) /  &
                                      (2.*pi6) )**(1./3.)
             precp(ii,jj,bb)%vratiohi = precp(ii,jj,bb)%vhilim / ( pi6*precp(ii,jj,bb)%dmid**3 )
             precp(ii,jj,bb)%vratiolo = precp(ii,jj,bb)%vlolim / ( pi6*precp(ii,jj,bb)%dmid**3 )

             ! Initialize the wet diameter as the bin mid diameter
             precp(ii,jj,bb)%dwet = precp(ii,jj,bb)%dmid

             DO cc = 1,8
                precp(ii,jj,bb)%volc(cc) = 0.
             END DO
             precp(ii,jj,bb)%numc = 0.
             precp(ii,jj,bb)%core = 0.
             precp(ii,jj,bb)%veqh2o = 0.

          END DO

       END DO
    END DO

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(cloudbins(ncld))
    DO bb = 1,ncld
       cloudbins(bb) = (cloud(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
    cloudbins = 0.5*cloudbins ! To radius
    ALLOCATE(precpbins(nprc))
    DO bb = 1,nprc
       precpbins(bb) = (precp(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
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
                               nicebin,         &
                               dminice,         &
                               in1a,fn2a,       &
                               in2b,fn2b,       &
                               icebins,         &
                               snowbins
    USE mo_salsa_driver, ONLY : kbdim, klev, &
                                 ice,snow,aero

    IMPLICIT NONE

    INTEGER :: ii,jj,cc,bb

    ! For ice bins
    REAL :: zvoldiff(fn2a)
    INTEGER :: zindex(fn2a)
    INTEGER :: imin, nba,nbb
    INTEGER :: armin(1)
    LOGICAL :: l_min(fn2a)

    REAL :: tmplolim(7), tmphilim(7)

    ! Helper arrays to set up snow size bins
    tmplolim = (/50.,55.,65., 100.,200.,500., 1000./)*1.e-6
    tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6

    ! For setting ice particle bins
    DO ii = in1a,fn2a
       zindex(ii) = ii
    END DO

    ! Determine the smallest full aerosol bin with which the smallest ice particle bin will coincide
    zvoldiff(in1a:fn2a) = ABS(aero(1,1,in1a:fn2a)%vlolim - pi6*dminice**3)
    l_min = ( zvoldiff(in1a:fn2a) == MINVAL(zvoldiff(in1a:fn2a)) )
    IF (COUNT(l_min) /= 1) STOP 'Initialization error! Try changing the ice particle bin lower bound slightly'
    armin = PACK(zindex,l_min)
    imin = armin(1)
    IF (aero(1,1,imin)%vlolim < pi6*dminice**3) imin=imin+1
    IF (imin < in1a .OR. imin > fn2a) STOP 'Invalid first ice bin'

    ! Number of ice bins in regime a (soluble nuclei)
    nba = fn2a-imin+1
    ! Number of cloud bins in regime b (insoluble nuclei) !!huomhuom insoluble
    IF (aero(1,1,imin)%dmid > aero(1,1,in2b)%dmid) THEN
       ! All cloud bins are within regime 2
       nbb = nba
    ELSE
       ! The smallest regime a cloud bins go to regime 1
       nbb = fn2b-fn2a
    END IF

    nicebin(1) = nba
    nicebin(2) = nbb

    ! Reset ice bin indices accordingly. The two components give the cloud regime index,
    ! and the aerosol bin index with which they are parallel
    iia%cur = 1;                       iia%par = imin
    fia%cur = iia%cur + nicebin(1)-1;  fia%par = iia%par + nicebin(1)-1
    iib%cur = fia%cur + 1;             iib%par = fn2b - nicebin(2) + 1
    fib%cur = iib%cur + nicebin(2)-1;  fib%par = iib%par + nicebin(2)-1

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
          ! Set iceproperties (parallel to aerosol bins) !!!huomhuom tehdäänkö jäälle myös näin?
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

          ! Initialize the droplet diameter ("wet diameter") as the dry
          ! mid diameter of the nucleus to avoid problems later. !!huomhuom mitä tehdään jäälle
          ice(ii,jj,iia%cur:fib%cur)%dwet = ice(ii,jj,iia%cur:fib%cur)%dmid

          ! Initialize the volume and number concentrations for ice.
          ! First "real" values are only obtained upon the first calculation
          ! of the cloud droplet activation. !! huomhuom

          DO cc = 1,8
             ice(ii,jj,iia%cur:fib%cur)%volc(cc) = 0.
          END DO
          ice(ii,jj,iia%cur:fib%cur)%numc = 0.
          ice(ii,jj,iia%cur:fib%cur)%core = 0.
          ice(ii,jj,iia%cur:fib%cur)%veqh2o = 0.

          ! ---------------------------------------------------------------------------------------
          ! Set the snow properties; unlike aerosol and cloud bins, the size distribution
          ! goes according to the *wet* radius !!huomhuom minkä säteen mukaan menee?
          ! ---------------------------------------------------------------------------------------

          DO bb = isa,fsa

             snow(ii,jj,bb)%vhilim = pi6*tmphilim(bb)**3
             snow(ii,jj,bb)%vlolim = pi6*tmplolim(bb)**3
             snow(ii,jj,bb)%dmid = ( (snow(ii,jj,bb)%vlolim + snow(ii,jj,bb)%vhilim) /  &
                                      (2.*pi6) )**(1./3.)
             snow(ii,jj,bb)%vratiohi = snow(ii,jj,bb)%vhilim / ( pi6*snow(ii,jj,bb)%dmid**3 )
             snow(ii,jj,bb)%vratiolo = snow(ii,jj,bb)%vlolim / ( pi6*snow(ii,jj,bb)%dmid**3 )

             ! Initialize the wet diameter as the bin mid diameter
             snow(ii,jj,bb)%dwet = snow(ii,jj,bb)%dmid

             DO cc = 1,8
                snow(ii,jj,bb)%volc(cc) = 0.
             END DO
             snow(ii,jj,bb)%numc = 0.
             snow(ii,jj,bb)%core = 0.
             snow(ii,jj,bb)%veqh2o = 0.

          END DO

       END DO
    END DO

    ! Save bin limits to be delivered e.g. to host model if needed
    ALLOCATE(icebins(nice))
    DO bb = 1,nice
       icebins(bb) = (ice(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
    icebins = 0.5*icebins ! To radius
    ALLOCATE(snowbins(nsnw))
    DO bb = 1,nsnw
       snowbins(bb) = (snow(1,1,bb)%vlolim/pi6)**(1./3.)
    END DO
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
  SUBROUTINE define_salsa()

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
                               nlactiv,               &
                               nlactintst,            &
                               nlactbase,            &
                               nlichom,               &
                               nlichet,               &
                               nlicimmers,            &
                               nlicmelt,              &
                               dmincld,nbin,reglim,   &
                               nice,nsnw,             &
                               nspec,listspec,        &
                               volDistA, volDistB,    &
                               initliqice,            &
                               liqFracA, iceFracA,    &
                               liqFracB, iceFracB,    &
                               nf2a, isdtyp,          &
                               sigmag,dpg,n,nldebug,  &
                               rhlim

    IMPLICIT NONE


    NAMELIST /salsa/  &
         nldebug,     & ! switch for debug printing
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
         nlichom,     & ! Switch for homogeneous ice nucleation
         nlichet,     & ! Switch for heterogeneous ice nucleation
         nlicimmers,  & ! Switch for ice nucleation by immersion
         nlicmelt,    & ! Switch for ice'n'snow melting
         nbin,        & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         nice,        & ! number of ice bins
         nsnw,        & ! number of snow bins
         nlcndh2ocl,    & ! Condensation of water vapour on clouds (drizzle)
         nlcndh2oic,    & ! Condensation of water vapour on ice particles ! ice'n'snow
         nlcndh2oae,    & ! Condensation of water vapour on aerosols (FALSE -> equilibrium calc.)
         nlauto,        & ! Switch for autoconversion of cloud droplets to drizzle and rain
         nlautosnow,    & ! Switch for autoconversion of ice particles to snowing
         nlactiv,       & ! Master switch for cloud droplet activation
         nlactbase,     & ! Switch for parameterized cloud base activation
         nlactintst,    & ! Switch for interstitial activation based on particle growth and host model S
         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 - vertical profile, read from file
         reglim,        & ! Low/high diameter limits of the 2 aerosol size regimes (1d table with length 4)
         nbin,          & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         dmincld,       & ! Minimum hydrometeor bin diameter: the lower limit of the smallest hydrometeor bin
                          ! is set at or above. For now this is also the only way to control the number of
                          ! hydrometeor bins as they are assumed fully parallel with the aerosol bins.
                          ! Implementing different bin limits for hydromets than aerosols would require
                          ! somewhat extensive modification of the microphysical SALSA dynamics calculations
         nspec,         & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
                          ! Must be an array of length 7, with empty strings for unused stuff.
         volDistA,      & ! Initial relative contribution [0-1] of each species to particle volume in a-bins. Must be
                          ! an array of length 7, with zero for unused species.
         volDistB,      & ! Same as above but for b-bins
         nf2a,          & ! Number fraction of particles allocated to a-bins in regime 2. b-bins will get 1-nf2a

         ! ------------
         ! -- Juha: These should eventually be replaced with physical processing !!!! Dont use for liquid clouds!
         initliqice,  & ! initialize ice and liquid cloud particles from aerosol bins
         liqFracA,      & ! fraction of aerosols that are activated to liquid cloud droplets in A bins
         iceFracA,      & ! fraction of aerosols that are activated to ice cloud particles in A bins
         liqFracB,      & ! fraction of aerosols that are activated to liquid cloud droplets  in B bins
         iceFracB,      & ! fraction of aerosols that are activated to    ice cloud particles in B bins
         ! ---------------------------------------------------------------------------------------------------
         rhlim,         & ! Upper limit RH/100 for sals during initialization and spinup 

         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n                ! Number concentration for the 7 initial lognormal modes



    OPEN(11,STATUS='old',FILE='NAMELIST')
    READ(11,NML=salsa)
    CLOSE(11)

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

    USE mo_submctl, ONLY : nbin, nbin2, nbin3,     &
                               in1a,fn1a,in2a,fn2a,in2b,fn2b,  &
                               nbins, &
                               massacc

    IMPLICIT NONE

    ! Remember to call 'define_salsa' for namelist paramers before calling this subroutine!

    ! --1) Set derived indices
    nbin2 = 4
    nbin3 = nbin(2) - nbin2

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
