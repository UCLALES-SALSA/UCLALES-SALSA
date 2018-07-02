!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_INIT                                       *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to initialize the particle size grid and aerosol           *
!*   processes                                                  *
!*                                                              *
!****************************************************************
MODULE mo_salsa_init
  USE classSection, ONLY : Section
  USE classSpecies, ONLY : Species

  USE mo_submctl, ONLY : reglim,nbin,in1a,fn1a,in2a,fn2a,in2b,fn2b,           &
                         ica,fca,icb,fcb,ira,fra,iia,fia,iib,fib,isa,fsa,     &  
                         nbins, ncld, nprc, bloPrc, nice, nsnw, ntotal,             &
                         dlaero, dlcloud, dlprecp, dlice, dlsnow,             &
                         aerobins,cloudbins,precpbins,icebins,snowbins,       &
                         spec, nspec_dry, listspec, nlim, prlim, massacc, pi6 

  USE mo_salsa_types, ONLY : aero, cloud, precp, ice, snow, allSALSA, frozen, liquid,        &
                             iaero, faero, icloud, fcloud, iprecp, fprecp,    &
                             iice, fice, isnow, fsnow                      

  USE mo_salsa_driver, ONLY : kbdim,klev
                             
  USE mo_salsa_optical_properties, ONLY : initialize_optical_properties

  IMPLICIT NONE
   
CONTAINS

  !----------------------------------------------------------------------------------
  ! Subroutine set_masterbins: This will allocate the master array "allSALSA"
  !                            that holds size distributions and bin parameters
  !                            for all particle types. Arrays "aero", "cloud" etc.
  !                            are pointers associated with segments of allSALSA
  !                            to allow  easy access to specific particle types.
  !
  ! Juha Tonttila, FMI, 2017
  !-----------------------------

  SUBROUTINE set_masterbins(dumaero, dumcloud, dumprecp, dumice, dumsnow)
    TYPE(Section), INTENT(in) :: dumaero(kbdim,klev,nbins), dumcloud(kbdim,klev,ncld),  &
                                 dumprecp(kbdim,klev,nprc), dumice(kbdim,klev,nice),   &
                                 dumsnow(kbdim,klev,nsnw)

    INTEGER :: lo,hi
    INTEGER :: nspec

    nspec = spec%getNSpec()
    ntotal = nbins+ncld+nprc+nice+nsnw

    ! Allocate the combined particle size distribution array
    ALLOCATE(allSALSA(kbdim,klev,ntotal))
    allSALSA(:,:,:) = Section(0,nlim,dlaero)

    ! Associate pointers for specific particle types
    lo = 1 
    hi = nbins
    aero => allSALSA(:,:,lo:hi)
    iaero = lo; faero = hi

    lo = hi + 1 
    hi = hi + ncld
    cloud => allSALSA(:,:,lo:hi)
    icloud = lo; fcloud = hi

    lo = hi + 1 
    hi = hi + nprc
    precp => allSALSA(:,:,lo:hi)
    iprecp = lo; fprecp = hi

    lo = hi + 1 
    hi = hi + nice
    ice => allSALSA(:,:,lo:hi)
    iice = lo; fice = hi

    lo = hi + 1 
    hi = hi + nsnw
    snow => allSALSA(:,:,lo:hi)
    isnow = lo; fsnow = hi

    ! Associate some (potentially helpful) subcollections 
    ! (note: for this it is necessary to have all the particles containing liquid water consecutively, and then ice containing particles
    !  consecutively so that the indexing works)
    lo = 1
    hi = nbins + ncld + nprc
    liquid => allSALSA(:,:,lo:hi)
    
    lo = nbins + ncld + nprc + 1
    hi = nbins + ncld + nprc + nice + nsnw
    frozen => allSALSA(:,:,lo:hi)

    ! Copy the dummy size distributions to actual size dists.
    aero = dumaero
    cloud = dumcloud
    precp = dumprecp
    ice = dumice
    snow = dumsnow

  END SUBROUTINE set_masterbins



   ! fxm: when dpmid is used for calculating coagulation coefficients
   ! (only?), would it make more sense to use approximate wet radii
   ! e.g. for sea salt particles?
   !********************************************************************
   !
   ! Subroutine SET_SIZEBINS()
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

   SUBROUTINE set_aerobins(dumaero)
     IMPLICIT NONE

     TYPE(Section), INTENT(out), ALLOCATABLE :: dumaero(:,:,:)

     !-- local variables ----------
     INTEGER :: ii, jj,cc,dd ! loop indices
     INTEGER :: nbin2, nbin3
     REAL ::  ratio ! ratio of regime upper and lower diameter
     INTEGER :: nspec

     nspec = spec%getNSpec()

     nbin2 = 4
     nbin3 = nbin(2) - nbin2

     ! -------------------------------------------------
     ! Allocate and initialize the dummy aerosol tracers
     !--------------------------------------------------
     ALLOCATE(dumaero(kbdim,klev,nbins))
     DO jj = 1,klev
        DO ii = 1,kbdim
           DO cc = 1,nbins
              dumaero(ii,jj,cc) = Section(1,nlim,dlaero)
           END DO
        END DO
     END DO

     DO jj = 1, klev
        DO ii = 1, kbdim
           
           !-- 1) size regime 1: --------------------------------------
           !  - minimum & maximum *dry* volumes [fxm]
           !  - bin mid *dry* diameter [m]
           ratio = reglim(2)/reglim(1)   ! section spacing
           
           DO cc = in1a, fn1a
              dumaero(ii,jj,cc)%vlolim = pi6*(reglim(1)*ratio**(REAL(cc-1)/nbin(1)))**3
              dumaero(ii,jj,cc)%vhilim = pi6*(reglim(1)*ratio**(REAL(cc)/nbin(1)))**3
              dumaero(ii,jj,cc)%dmid = ( (dumaero(ii,jj,cc)%vhilim + dumaero(ii,jj,cc)%vlolim) /  &
                   (2.*pi6) )**(1./3.)
              dumaero(ii,jj,cc)%vratiohi = dumaero(ii,jj,cc)%vhilim/(pi6*dumaero(ii,jj,cc)%dmid**3)
              dumaero(ii,jj,cc)%vratiolo = dumaero(ii,jj,cc)%vlolim/(pi6*dumaero(ii,jj,cc)%dmid**3)
           END DO
           
           !-- 2) size regime 2: --------------------------------------
           !  - minimum & maximum *dry* volumes [fxm]
           !  - bin mid *dry* diameter [m]
           
           !-- 2.1) first for subregime 2a
           ratio = reglim(3)/reglim(2)   ! section spacing
           
           DO dd = in2a, in2a+nbin2-1
              cc = dd - in2a
              dumaero(ii,jj,dd)%vlolim = pi6*(reglim(2)*ratio**(REAL(cc)/nbin2))**3
              dumaero(ii,jj,dd)%vhilim = pi6*(reglim(2)*ratio**(REAL(cc+1)/nbin2))**3
              dumaero(ii,jj,dd)%dmid   = ( (dumaero(ii,jj,dd)%vhilim + dumaero(ii,jj,dd)%vlolim) /  &
                   (2.*pi6) )**(1./3.)
              dumaero(ii,jj,dd)%vratiohi = dumaero(ii,jj,dd)%vhilim/(pi6*dumaero(ii,jj,dd)%dmid**3)
              dumaero(ii,jj,dd)%vratiolo = dumaero(ii,jj,dd)%vlolim/(pi6*dumaero(ii,jj,dd)%dmid**3)
           END DO
           
           !-- 3) size regime 3: --------------------------------------
           !  - bin mid *dry* diameter [m]
           ratio = reglim(4)/reglim(3)   ! section spacing
           
           DO dd = in2a+nbin2, fn2a
              cc = dd - (fn2a-(nbin3-1))
               
             dumaero(ii,jj,dd)%vlolim = pi6*(reglim(3)*ratio**(REAL(cc)/nbin3))**3
              dumaero(ii,jj,dd)%vhilim = pi6*(reglim(3)*ratio**(REAL(cc+1)/nbin3))**3
              dumaero(ii,jj,dd)%dmid = ( (dumaero(ii,jj,dd)%vhilim + dumaero(ii,jj,dd)%vlolim) /  &
                   (2.*pi6) )**(1./3.)
              dumaero(ii,jj,dd)%vratiohi = dumaero(ii,jj,dd)%vhilim/(pi6*dumaero(ii,jj,dd)%dmid**3)
              dumaero(ii,jj,dd)%vratiolo = dumaero(ii,jj,dd)%vlolim/(pi6*dumaero(ii,jj,dd)%dmid**3)
           END DO
           
           !-- 2.2) same values for subregime 2b
           dumaero(ii,jj,in2b:fn2b)%vlolim = dumaero(ii,jj,in2a:fn2a)%vlolim
           dumaero(ii,jj,in2b:fn2b)%vhilim = dumaero(ii,jj,in2a:fn2a)%vhilim
           dumaero(ii,jj,in2b:fn2b)%dmid = dumaero(ii,jj,in2a:fn2a)%dmid
           dumaero(ii,jj,in2b:fn2b)%vratiohi = dumaero(ii,jj,in2a:fn2a)%vratiohi
           dumaero(ii,jj,in2b:fn2b)%vratiolo = dumaero(ii,jj,in2a:fn2a)%vratiolo
           
           ! Initialize the wet diameter with the bin dry diameter to avoid numerical proplems later
           dumaero(ii,jj,:)%dwet = dumaero(ii,jj,:)%dmid
           
        END DO !ii
     END DO !!jj
     
     ! Save bin limits to be delivered e.g. to host model if needed
     ALLOCATE(aerobins(nbins))
     DO cc = 1, nbins
        aerobins(cc) = (dumaero(1,1,cc)%vlolim/pi6)**(1./3.)
     END DO
     
   END SUBROUTINE set_aerobins
   
   !--------------------------------------------------------------------------
   !
   ! *******************************
   ! Subroutine set_cloudbins
   ! *******************************
   !
   ! Setup of hydrometeor size bins similar to the subroutine *set_sizebins*
   !
   ! Juha Tonttila (FMI) 2014
   !
   !---------------------------------------------------------------------------

   SUBROUTINE set_cloudbins(dumaero,dumcloud,dumprecp)
     IMPLICIT NONE
     
     TYPE(Section), INTENT(in) :: dumaero(kbdim,klev,nbins)
     TYPE(Section), INTENT(out), ALLOCATABLE :: dumcloud(:,:,:), dumprecp(:,:,:)

     INTEGER :: ii,jj,bb,nba,nbb

     REAL, ALLOCATABLE :: tmplolim(:), tmphilim(:)
     
     INTEGER :: nspec

     nspec = spec%getNSpec()
     nprc = bloPrc%nbins
     
     ! Bin diameter limits for precipitation bins
     CALL buildBinLimits(bloPrc, tmplolim, tmphilim)
     
     ! Cloud bins are parallel with the 2a and 2b aerosol bins

     ! Number of cloud bins in regime a
     nba = fn2a-in2a+1
     ! Number of cloud bins in regime b
     nbb = nba
     
     ! Reset cloud bin indices accordingly. The two components give the cloud regime index,
     ! and the aerosol bin index with which they are parallel
     ica%cur = 1;               ica%par = in2a
     fca%cur = ica%cur + nba-1; fca%par = ica%par + nba-1
     icb%cur = fca%cur + 1;     icb%par = fn2b - nbb + 1
     fcb%cur = icb%cur + nbb-1; fcb%par = icb%par + nbb-1
     
     ncld = fcb%cur
     
     ! Rain/drizzle bins
     ira = 1; fra = nprc;
     
     
     ! ----------------------------------------
     ! Allocate cloud and precipitation arrays
     ! ----------------------------------------
     ALLOCATE(dumcloud(kbdim,klev,ncld), dumprecp(kbdim,klev,nprc))
     DO jj = 1,klev
        DO ii = 1,kbdim
           DO bb = 1,ncld
              dumcloud(ii,jj,bb) = Section(2,nlim,dlcloud)
           END DO
           DO bb = 1,nprc
              dumprecp(ii,jj,bb) = Section(3,prlim,dlprecp)
           END DO
        END DO
     END DO
     
     DO jj = 1, klev
        DO ii = 1, kbdim
           
           ! -------------------------------------------------
           ! Set cloud properties (parallel to aerosol bins)
           ! -------------------------------------------------
           dumcloud(ii,jj,ica%cur:fca%cur)%vhilim = dumaero(ii,jj,ica%par:fca%par)%vhilim
           dumcloud(ii,jj,icb%cur:fcb%cur)%vhilim = dumaero(ii,jj,icb%par:fcb%par)%vhilim
           
           dumcloud(ii,jj,ica%cur:fca%cur)%vlolim = dumaero(ii,jj,ica%par:fca%par)%vlolim
           dumcloud(ii,jj,icb%cur:fcb%cur)%vlolim = dumaero(ii,jj,icb%par:fcb%par)%vlolim
           
           dumcloud(ii,jj,ica%cur:fca%cur)%vratiohi = dumaero(ii,jj,ica%par:fca%par)%vratiohi
           dumcloud(ii,jj,icb%cur:fcb%cur)%vratiohi = dumaero(ii,jj,icb%par:fcb%par)%vratiohi
           
           dumcloud(ii,jj,ica%cur:fca%cur)%vratiolo = dumaero(ii,jj,ica%par:fca%par)%vratiolo
           dumcloud(ii,jj,icb%cur:fcb%cur)%vratiolo = dumaero(ii,jj,icb%par:fcb%par)%vratiolo
           
           dumcloud(ii,jj,ica%cur:fca%cur)%dmid = dumaero(ii,jj,ica%par:fca%par)%dmid
           dumcloud(ii,jj,icb%cur:fcb%cur)%dmid = dumaero(ii,jj,icb%par:fcb%par)%dmid
           
           ! Initialize the droplet diameter ("wet diameter") as the dry
           ! mid diameter of the nucleus to avoid problems later.
           dumcloud(ii,jj,ica%cur:fcb%cur)%dwet = dumcloud(ii,jj,ica%cur:fcb%cur)%dmid
           
           ! ---------------------------------------------------------------------------------------
           ! Set the precipitation properties; unlike aerosol and cloud bins, the size distribution
           ! goes according to the *wet* radius
           ! ---------------------------------------------------------------------------------------
           DO bb = ira, fra
              dumprecp(ii,jj,bb)%vhilim = pi6*tmphilim(bb)**3
              dumprecp(ii,jj,bb)%vlolim = pi6*tmplolim(bb)**3
              dumprecp(ii,jj,bb)%dmid = ( (dumprecp(ii,jj,bb)%vlolim + dumprecp(ii,jj,bb)%vhilim) /  &
                   (2.*pi6) )**(1./3.)
              dumprecp(ii,jj,bb)%vratiohi = dumprecp(ii,jj,bb)%vhilim / ( pi6*dumprecp(ii,jj,bb)%dmid**3 )
              dumprecp(ii,jj,bb)%vratiolo = dumprecp(ii,jj,bb)%vlolim / ( pi6*dumprecp(ii,jj,bb)%dmid**3 )
              
              ! Initialize the wet diameter as the bin mid diameter
              dumprecp(ii,jj,bb)%dwet = dumprecp(ii,jj,bb)%dmid
              
           END DO
           
        END DO
     END DO
     
     ! Save bin limits to be delivered e.g. to host model if needed
     ALLOCATE(cloudbins(ncld))
     DO bb = 1, ncld
        cloudbins(bb) = (dumcloud(1,1,bb)%vlolim/pi6)**(1./3.)
     END DO
     ALLOCATE(precpbins(nprc))
     DO bb = 1, nprc
        precpbins(bb) = (dumprecp(1,1,bb)%vlolim/pi6)**(1./3.)
     END DO
     
   END SUBROUTINE set_cloudbins

   !--------------------------------------------------------------------------
   !
   ! *******************************
   ! Subroutine set_icebins
   ! *******************************
   !
   ! Setup of ice size bins similar to the subroutines *set_sizebins* & *set_cloudbins*
   !
   ! Jaakko Ahola (FMI) 2015
   !
   !---------------------------------------------------------------------------
   SUBROUTINE set_icebins(dumaero,dumice,dumsnow)     
     IMPLICIT NONE
     
     TYPE(Section), INTENT(in) :: dumaero(kbdim,klev,nbins)
     TYPE(Section), INTENT(out), ALLOCATABLE :: dumice(:,:,:), dumsnow(:,:,:)
     
     INTEGER :: ii,jj,bb,nba,nbb
     
     REAL :: tmplolim(7), tmphilim(7)

     INTEGER :: nspec
     
     nspec = spec%getNSpec()

     ! Helper arrays to set up snow size bins
     tmplolim = (/50.,55.,65., 100.,200.,500., 1000./)*1.e-6
     tmphilim = (/55.,65.,100.,200.,500.,1000.,2000./)*1.e-6
     
     ! Number of ice bins in regime a (soluble nuclei)
     nba = fn2a-in2a+1
     ! Number of cloud bins in regime b (insoluble nuclei)
     nbb = nba
     
     ! Reset ice bin indices accordingly. The two components give the cloud regime index,
     ! and the aerosol bin index with which they are parallel
     iia%cur = 1;                iia%par = in2a
     fia%cur = iia%cur + nba-1;  fia%par = iia%par + nba-1
     iib%cur = fia%cur + 1;      iib%par = fn2b - nbb + 1
     fib%cur = iib%cur + nbb-1;  fib%par = iib%par + nbb-1
     
     nice = fib%cur
     
     ! snow bins
     isa = 1; fsa = 7;
     nsnw = fsa
     
     ! ----------------------------------------
     ! Allocate ice arrays
     ! ----------------------------------------
     ALLOCATE(dumice(kbdim,klev,nice), dumsnow(kbdim,klev,nsnw))
     DO jj = 1,klev
        DO ii = 1,kbdim
           DO bb = 1,nice
              dumice(ii,jj,bb) = Section(4,prlim,dlice)
           END DO
           DO bb = 1,nsnw
              dumsnow(ii,jj,bb) = Section(5,prlim,dlsnow)
           END DO
        END DO
     END DO

     DO jj = 1, klev
        DO ii = 1, kbdim
           
           ! -------------------------------------------------
           ! Set iceproperties (parallel to aerosol bins) 
           ! -------------------------------------------------
           dumice(ii,jj,iia%cur:fia%cur)%vhilim = dumaero(ii,jj,iia%par:fia%par)%vhilim
           dumice(ii,jj,iib%cur:fib%cur)%vhilim = dumaero(ii,jj,iib%par:fib%par)%vhilim
           
           dumice(ii,jj,iia%cur:fia%cur)%vlolim = dumaero(ii,jj,iia%par:fia%par)%vlolim
           dumice(ii,jj,iib%cur:fib%cur)%vlolim = dumaero(ii,jj,iib%par:fib%par)%vlolim
           
           dumice(ii,jj,iia%cur:fia%cur)%vratiohi = dumaero(ii,jj,iia%par:fia%par)%vratiohi
           dumice(ii,jj,iib%cur:fib%cur)%vratiohi = dumaero(ii,jj,iib%par:fib%par)%vratiohi
           
           dumice(ii,jj,iia%cur:fia%cur)%vratiolo = dumaero(ii,jj,iia%par:fia%par)%vratiolo
           dumice(ii,jj,iib%cur:fib%cur)%vratiolo = dumaero(ii,jj,iib%par:fib%par)%vratiolo
           
           dumice(ii,jj,iia%cur:fia%cur)%dmid = dumaero(ii,jj,iia%par:fia%par)%dmid
           dumice(ii,jj,iib%cur:fib%cur)%dmid = dumaero(ii,jj,iib%par:fib%par)%dmid
           
           ! Initialize the "wet" diameter as the dry mid diameter of the nucleus
           dumice(ii,jj,iia%cur:fib%cur)%dwet = dumice(ii,jj,iia%cur:fib%cur)%dmid
           
           ! ---------------------------------------------------------------------------------------
           ! Set the snow properties; unlike aerosol and cloud bins, the size distribution
           ! goes according to the *wet* radius !!
           ! ---------------------------------------------------------------------------------------
           
           DO bb = isa, fsa
              
              dumsnow(ii,jj,bb)%vhilim = pi6*tmphilim(bb)**3
              dumsnow(ii,jj,bb)%vlolim = pi6*tmplolim(bb)**3
              dumsnow(ii,jj,bb)%dmid = ( (dumsnow(ii,jj,bb)%vlolim + dumsnow(ii,jj,bb)%vhilim) /  &
                   (2.*pi6) )**(1./3.)
              dumsnow(ii,jj,bb)%vratiohi = dumsnow(ii,jj,bb)%vhilim / ( pi6*dumsnow(ii,jj,bb)%dmid**3 )
              dumsnow(ii,jj,bb)%vratiolo = dumsnow(ii,jj,bb)%vlolim / ( pi6*dumsnow(ii,jj,bb)%dmid**3 )
              
              ! Initialize the wet diameter as the bin mid diameter
              dumsnow(ii,jj,bb)%dwet = dumsnow(ii,jj,bb)%dmid
              
           END DO
           
        END DO
     END DO
     
     ! Save bin limits to be delivered e.g. to host model if needed
     ALLOCATE(icebins(nice))
     DO bb = 1, nice
        icebins(bb) = (dumice(1,1,bb)%vlolim/pi6)**(1./3.)
     END DO
     ALLOCATE(snowbins(nsnw))
     DO bb = 1, nsnw
        snowbins(bb) = (dumsnow(1,1,bb)%vlolim/pi6)**(1./3.)
     END DO
     
   END SUBROUTINE set_icebins
   
   ! ---------------------------------------------------
   
   SUBROUTINE buildBinLimits(blo, dlolim, dhilim)
     USE classBinLayout, ONLY : BinLayout
     TYPE(BinLayout), INTENT(in) :: blo   ! Bin layout instance
     REAL, ALLOCATABLE, INTENT(out) :: dlolim(:), dhilim(:)   ! Vectors for the low and upper limit diameters

     INTEGER :: N,ii
     REAL :: first, ratio

     N = blo%nbins
     ALLOCATE(dlolim(N), dhilim(N))
     dlolim = 0.; dhilim = 0.

     first = blo%dlo
     ratio = blo%vol_ratio
     dlolim(1) = first
     dhilim(1) = first*( ratio**(1./3.) )
     DO ii = 2,N
        dlolim(ii) = dlolim(ii-1)*( ratio**(1./3.) )
        dhilim(ii) = dlolim(ii)*( ratio**(1./3.) )
     END DO
     
   END SUBROUTINE buildBinLimits


   !----------------------------------------------------------------------
   !
   ! *************************
   ! Subroutine define_salsa
   ! *************************
   !
   ! Reads LOGICAL switches and aerosol/hydrometeor size bin definitions
   ! from a namelist.
   !
   ! Juha Tonttila (FMI) 2014
   !
   !----------------------------------------------------------------------
   SUBROUTINE define_salsa(level)

      USE mo_submctl, ONLY : lscoag,                      &
                             lscnd,                       &
                             lsauto,                      &
                             lsautosnow,                  & 
                             lsactiv,                     &
                             lsicenucl,                   &
                             lsicemelt,                   &

                             lscgaa,lscgcc,lscgpp,        &
                             lscgca,lscgpa,lscgpc,        &
                             lscgia,lscgic,lscgii,        &
                             lscgip,lscgsa,lscgsc,        &
                             lscgsi,lscgsp,lscgss,        &

                             lscndgas,                    &
                             lscndh2oae,lscndh2ocl,       &
                             lscndh2oic,                  &

                             lsdistupdate,          &
                             lscheckarrays,         &
                             fixINC,                &
                             ice_hom, ice_imm, ice_dep, &

                             bloPrc,                      &                            
                             nbin,reglim,                 &
                             nice,nsnw,                   &
                             nspec_dry,listspec,          &
                             volDistA, volDistB,          &
                             nf2a, isdtyp,                &
                             sigmag,dpg,n,                &
                             lsfreeRH,rhlim

      IMPLICIT NONE

    INTEGER, INTENT(in) :: level

    NAMELIST /salsa/  &
         ! Master process switches
         lscoag,      &
         lscnd,       &
         lsauto,      &     ! mode 1: parameterized simple autoconversion, mode 2: more elaborate precipitation formation based on coagulation
         lsautosnow,  &
         lsactiv,     &     ! mode 1: Aerosol growth-based activation, mode 2: Parameterized cloud base activation 
         lsicenucl,   &
         lsicemelt,   &

         ! Subprocess switches
         lscgaa,      & ! Coagulation between aerosols
         lscgcc,      & ! Collision-coalescence between cloud droplets
         lscgpp,      & ! Collisions between rain drops
         lscgca,      & ! Cloud collection of aerosols
         lscgpa,      & ! Collection of aerosols by precip
         lscgpc,      & ! Collection of cloud droplets by rain
         lscgia,      & ! Ice collection of aerosols
         lscgic,      & ! Collection of cloud droplets by ice particles
         lscgii,      & ! Collision-coalescence between ice particles
         lscgip,      & ! Collection of precipitation by ice particles
         lscgsa,      & ! Collection of aerosols by snow
         lscgsc,      & ! Collection of cloud droplets by snow
         lscgsi,      & ! Collection of ice by snow
         lscgsp,      & ! Collection of precipitation by snow
         lscgss,      & ! Collision-coalescence between snow particles

         lscndgas,    & ! Condensation of precursor gases
         lscndh2ocl,    & ! Condensation of water vapour on clouds (drizzle)
         lscndh2oic,    & ! Condensation of water vapour on ice particles ! ice'n'snow
         lscndh2oae,    & ! Condensation of water vapour on aerosols (FALSE -> equilibrium calc.)

         lsdistupdate,  & ! Switch for size dsitribution update
         lscheckarrays, & ! Switch for runnin the array check routine in mo_salsa

         fixINC,      & ! fixed ice number concentration #/kg
         ice_hom,     & ! Switch for homogeneous ice nucleation
         ice_imm,     & ! .. for immersio freezing
         ice_dep,     & ! .. for deposition freezing

         nbin,        & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         nice,        & ! number of ice bins  DONT USE YET, WILL NOT WORK
         nsnw,        & ! number of snow bins DONT USE YET, WILL NOT WORK
         
         bloPrc,      & ! Precipitation bin definitions
         
         isdtyp,        & ! Type of initial size distribution: 0 - uniform; 1 - vertical profile, read from file
         reglim,        & ! Low/high diameter limits of the 2 aerosol size regimes (1d table with length 4)
         nbin,          & ! Number of bins used for each of the aerosol size regimes (1d table with length 2)
         nspec_dry,     & ! Number of aerosol species used in the model
         listspec,      & ! List of strings specifying the names of the aerosol species that are active.
                          ! Must be an array of length 7, with empty strings for unused stuff.
         volDistA,      & ! Initial relative contribution [0-1] of each species to particle volume in a-bins. Must be
                          ! an array of length 7, with zero for unused species.
         volDistB,      & ! Same as above but for b-bins
         nf2a,          & ! Number fraction of particles allocated to a-bins in regime 2. b-bins will get 1-nf2a

         lsfreeRH,      & ! Switch for using rhlim
         rhlim,         & ! Upper limit RH/100 for sals during initialization and spinup 

         sigmag,        & ! Stdev for the 7 initial lognormal modes
         dpg,           & ! Mean diameter for the 7 initial lognormal modes
         n                ! Number concentration for the 7 initial lognormal modes


      ! Associate master switch pointers before reading NAMELIST
      CALL associate_master_switches()
      
      ! Set default values to BinLayout objects
      CALL setDefaultBinLayouts()

      ! Read the NAMELIST
      OPEN(11,STATUS='old',FILE='NAMELIST')
      READ(11,NML=salsa)
      CLOSE(11)

      ! if thermodynamical level is 4, set all ice process switches to false
      IF(level == 4) THEN
            lscgia      = .false.
            lscgic      = .false.
            lscgii      = .false.
            lscgip      = .false.
            lscgsa      = .false.
            lscgsc      = .false.
            lscgsi      = .false.
            lscgsp      = .false.
            lscgss      = .false.

            lscndh2oic  = .false.

            lsautosnow%switch = .FALSE.
            lsicenucl%switch = .FALSE.
            lsicemelt%switch = .FALSE.

      END IF !level

   END SUBROUTINE define_salsa

   ! -----------------------------
   ! 
   ! ***********************************
   ! SUBROUTINE associate_master_switches
   ! Associate master switch pointers
   !
   SUBROUTINE associate_master_switches
     USE classProcessSwitch, ONLY : ProcessSwitch
     USE mo_submctl, ONLY : Nmaster, lsmaster, lscoag, lscnd, lsauto,  &
                            lsautosnow, lsactiv, lsicenucl,   &
                            lsicemelt, lsfreeRH
     IMPLICIT NONE
     
     INTEGER :: i

     ! Initialize the values
     DO i = 1,Nmaster
        lsmaster(i) = ProcessSwitch()
     END DO

     ! Associate pointers
     lscoag => lsmaster(1)
     lscnd => lsmaster(2)
     lsauto => lsmaster(3)
     lsautosnow => lsmaster(4)
     lsactiv => lsmaster(5)
     lsicenucl => lsmaster(6)
     lsicemelt => lsmaster(7)

     ! Use this to initialize also some other switches
     lsfreeRH = ProcessSwitch()

   END SUBROUTINE associate_master_switches
   
   SUBROUTINE setDefaultBinLayouts
     USE mo_submctl, ONLY : bloPrc
     USE classBinLayout, ONLY : BinLayout
     IMPLICIT NONE

     bloPrc = BinLayout(20,20.e-6,2.)
    
   END SUBROUTINE setDefaultBinLayouts
   

   !-------------------------------------------------------------------------------
   !
   ! *****************************
   ! Subroutine salsa_initialize
   ! *****************************
   !
   ! SALSA initializations. Modified and rewritten for more dynamic control
   ! and LES implementation.
   !
   ! define_salsa **MUST** be called before salsa_initialize so that the 
   ! NAMELIST-parameters will have an effect.
   !
   ! Juha Tonttila (FMI) 2014
   !
   !-------------------------------------------------------------------------------
   SUBROUTINE salsa_initialize()

      !
      !-------------------------------------------------------------------------------


      IMPLICIT NONE

      ! Dummy size distributions just for setting everything up!!
      ! May not be the smartest or the fastest way, but revise later... 
      TYPE(Section), ALLOCATABLE :: dumaero(:,:,:), dumcloud(:,:,:), dumprecp(:,:,:), &
                                      dumice(:,:,:), dumsnow(:,:,:)
      INTEGER :: nspec

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
      
      ! Initialize and sort pointers to aerosol properties according to the order in which the species are given in the NAMELIST
      spec = Species(nspec_dry,listspec)

      ! -- Aerosol tracers are allocated in *set_aerobins*
      ! -- Hydrometeor tracer in *set_cloudbins*

      ! --3) Call other initialization routines
      CALL set_aerobins(dumaero)

      CALL set_cloudbins(dumaero, dumcloud, dumprecp)

      CALL set_icebins(dumaero, dumice, dumsnow)

      CALL set_masterbins(dumaero, dumcloud, dumprecp, dumice, dumsnow)

      ! Initialize aerosol optical properties - uses settings from "spec"
      CALL initialize_optical_properties()
      
      ! Initialize the coagulation kernel arrays and sink/source term arrays for processes
      nspec = spec%getNSpec()
      !CALL initialize_coagulation_kernels(kbdim,klev)
      !CALL initialize_coagulation_processes(kbdim,klev,nspec)


      IF ( ALLOCATED(dumaero) ) DEALLOCATE(dumaero)
      IF ( ALLOCATED(dumcloud)) DEALLOCATE(dumcloud)
      IF ( ALLOCATED(dumprecp)) DEALLOCATE(dumprecp)
      IF ( ALLOCATED(dumice)) DEALLOCATE(dumice)
      IF ( ALLOCATED(dumsnow)) DEALLOCATE(dumsnow)

   END SUBROUTINE salsa_initialize


END MODULE mo_salsa_init
