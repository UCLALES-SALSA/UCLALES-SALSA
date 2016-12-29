
!****************************************************************
!*                                                              *
!*   module MO_SALSA_DYNAMICS                               *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************

MODULE mo_salsa_dynamics


CONTAINS

  ! this calculated for empty bins too!!!
  ! fxm: test well, esp. self-coagulation (but other bits too!)
  ! AL_note: Diagnostic variables of cond and nucl mass
  !********************************************************************
  !
  ! subroutine COAGULATION(kproma,kbdim,klev, &
  !       pnaero,pvols,pdwet, &
  !       pcore, ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates particle loss and change in size distribution
  !  due to (Brownian) coagulation
  !
  !
  ! Method:
  ! -------
  ! Semi-implicit, non-iterative method:
  !  Volume concentrations of the smaller colliding particles
  !  added to the bin of the larger colliding particles.
  !  Start from first bin and use the updated number and volume
  !  for calculation of following bins. NB! Our bin numbering
  !  does not follow particle size in regime 2.
  !
  !Schematic for bin numbers in different regimes:
  !             1                            2
  !    +-------------------------------------------+
  !  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
  !  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
  !    +-------------------------------------------+
  !
  ! Exact coagulation coefficients for each pressure level
  !  are calculated in subroutine SET_COAGC (in mo_salsa_init)
  !  which is called once at the beginning of the simulation
  !  from model driver. In subroutine COAGULATION, these exact
  !  coefficients are scaled according to current particle wet size
  !  (linear scaling).
  !
  ! Juha: Now also considers coagulation between hydrometeors,
  !       and hydrometeors and aerosols.
  !
  !       Since the bins are organized in terms of the dry size of
  !       of the condensation nucleus, while coagulation kernell is
  !       calculated with the actual hydrometeor size, some assumptions
  !       are laid out:
  !                 1. Cloud droplets from each size bin are lost by
  !                    coagulation with other cloud droplets that have
  !                    larger condensation nucleus.
  !
  !                 2. Cloud droplets from each size bin are lost by
  !                    coagulation with all drizzle bins, regardless of
  !                    the nucleus size in the latter (collection of cloud
  !                    droplets by rain).
  !
  !                 3. Coagulation between drizzle bins acts like 1.
  !
  !       ISSUES:
  !           Process selection should be made smarter - now just lots of IFs
  !           inside loops. Bad.
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005
  ! Harri Kokkola (FMI) 2006
  ! Tommi Bergman (FMI) 2012
  ! Matti Niskanen(FMI) 2012
  ! Anton Laakso  (FMI) 2013
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------------


  SUBROUTINE coagulation(kproma, kbdim,  klev,    &
                         paero,  pcloud, pprecp, pice, psnow,  &
                         ptstep, ptemp,  ppres     )

    USE mo_submctl, ONLY:        &
         t_parallelbin, t_section,   & ! Datatypes for the cloud bin representation
         in1a, fn1a,                 & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 &
         ica,fca,icb,fcb,            &
         ncld, nprc,        &
         iia,fia,iib,fib,            &
         nice, nsnw,         &
         pi6,                        &
         rhosu,rhooc,rhono,rhonh,    &
         rhobc,rhodu,rhoss,rhowa,    &
         rhoic,                      & ! density of ice  !ice'n'snow
         rhosn,                      & ! density of snow
         nlim,prlim,                 &
         lscgaa, lscgcc, lscgca,     &
         lscgpp, lscgpa, lscgpc,     &
         lscgia, lscgic, lscgii, lscgip, &
         lscgsa, lscgsc, lscgsi, lscgsp, lscgss, &
         debug

    IMPLICIT NONE


    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical klev

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &  ! Hydrometeor properties
         paero(kbdim,klev,fn2b),      &  ! Aerosol properties
         pprecp(kbdim,klev,nprc),     &  ! precipitation properties
         pice(kbdim,klev,nice),       & ! ice properties
         psnow(kbdim,klev,nsnw)        ! snow properties

    REAL, INTENT(IN) ::         &
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
    !-- Local variables ------------------------
    INTEGER ::                      &
         ii,jj,kk,ll,mm,nn,cc,      & ! loop indices
         index_2a, index_2b        ! corresponding bin in regime 2a/2b

    REAL ::                     &
         zcc(fn2b,fn2b),            & ! updated coagulation coefficients [m3/s]
         zcccc(ncld,ncld),          & ! - '' - for collision-coalescence [m3/s]
         zccca(fn2b,ncld),          & ! - '' - for cloud collection of aerosols [m3/s]
         zccpc(ncld,nprc),          & ! - '' - for collection of cloud droplets by precip [m3/s]
         zccpa(fn2b,nprc),          & ! - '' - for collection of aerosols by precip
         zccpp(nprc,nprc),          & ! - '' - for collitions between precip particles (neglected?)
         zccia(fn2b,nice),          & ! - '' - for collection of aerosols by ice !!huomhuom
         zccic(ncld,nice),          & ! - '' - for collection of cloud particles droplets by ice !!huomhuom
         zccii(nice,nice),          & ! - '' - for collitions between ice particles !!huomhuom
         zccip(nprc,nice),          & ! - '' - for collection of precip by ice-collision !!huomhuom
         zccsa(fn2b,nsnw),          & ! - '' - for collection of aerosols by snow !!huomhuom
         zccsc(ncld,nsnw),          & ! - '' - for collection of cloud droples by snow !!huomhuom
         zccsi(nice,nsnw),          & ! - '' - for collection of ice by snow !!huomhuom
         zccsp(nprc,nsnw),          & ! - '' - for collection of precip by snow !!huomhuom
         zccss(nsnw,nsnw),          & ! - '' - for collitions between snow particles !!huomhuom
         zminusterm,                & ! coagulation loss in a bin [1/s]
         zplusterm(8)                 ! coagulation gain in a bin [fxm/s]
                                      ! (for each chemical compound)

    REAL :: &
         zmpart(fn2b),   & ! approximate mass of particles [kg]
         zmcloud(ncld),  &    ! approximate mass of cloud droplets [kg]
         zmprecp(nprc),  & ! Approximate mass for rain drops [kg]
         zmice(nice),     & ! approximate mass for ice particles [kg] !huomhuom
         zmsnow(nsnw), &  ! approximate mass for snow particles [kg] !!huomhuom
         zdpart(fn2b),   & ! diameter of particles [kg]
         zdcloud(ncld),  &   ! diameter of cloud droplets [kg]
         zdprecp(nprc),  & ! diameter for rain drops [kg]
         zdice(nice),     & ! diameter for ice particles [kg] !huomhuom
         zdsnow(nsnw)     ! diameter for snow particles [kg] !!huomhuom

    REAL :: &
         temppi,pressi

    LOGICAL :: any_cloud, any_precp, any_ice, any_snow

    !-----------------------------------------------------------------------------
    !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
    !      CoagSink ~ Dp in continuum regime, thus we calculate
    !      'effective' number concentration of coarse particles


    !-- 2) Updating coagulation coefficients -------------------------------------

    IF (debug) WRITE(*,*) 'start coagulation kernels'

     DO jj = 1,klev      ! vertical grid
        DO ii = 1,kbdim ! horizontal kproma in the slab
           ! Which species are included
           any_cloud = ANY(pcloud(ii,jj,:)%numc > nlim)
           any_precp = ANY(pprecp(ii,jj,:)%numc > prlim)
           any_ice = ANY(pice(ii,jj,:)%numc > prlim)
           any_snow = ANY(psnow(ii,jj,:)%numc > prlim)

           !-- Aerosol diameter [m] and mass [kg]; density of 1500 kg/m3 assumed
           zdpart(1:fn2b) = MIN(paero(ii,jj,1:fn2b)%dwet, 30.e-6) ! Limit to 30 um
           zmpart(1:fn2b) = pi6*(zdpart(1:fn2b)**3)*1500. ! Should calculate masses here (%dwet may be inaccurate)

           !-- Cloud droplet diameter and mass; Assume water density
           zdcloud(1:ncld) = pcloud(ii,jj,1:ncld)%dwet ! No size limit?
           zmcloud(1:ncld) = pi6*(zdcloud(1:ncld)**3)*rhowa

           !-- Precipitation droplet diameter and mass
           zdprecp(1:nprc) = MIN(pprecp(ii,jj,1:nprc)%dwet, 2.e-3) ! Limit to 2 mm
           zmprecp(1:nprc) = pi6*(zdprecp(1:nprc)**3)*rhowa

           !-- Ice particle diameter and mass
           zdice(1:nice) = MIN(pice(ii,jj,1:nice)%dwet, 2.e-3) ! Limit to 2 mm
           zmice(1:nice) =   pi6*(zdice(1:nice)**3)*rhoic

           !-- Snow diameter and mass
           zdsnow(1:nsnw) = MIN(psnow(ii,jj,1:nsnw)%dwet, 2.e-3) ! Limit to 2 mm (too low!)
           zmsnow(1:nsnw) =  pi6*(zdsnow(1:nsnw)**3)*rhosn

           temppi=ptemp(ii,jj)
           pressi=ppres(ii,jj)
           zcc = 0.
           zcccc = 0.
           zccca = 0.
           zccpp = 0.
           zccpc = 0.
           zccpa = 0.
           zccia = 0.
           zccic = 0.
           zccii = 0.
           zccip = 0.
           zccsa = 0.
           zccsc = 0.
           zccsi = 0.
           zccsp = 0.
           zccss = 0.

           ! Aero-aero coagulation
           IF (lscgaa) THEN
              DO mm = 1,fn2b         ! smaller colliding particle
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = mm,fn2b            ! larger colliding particle
                    IF (paero(ii,jj,nn)%numc<nlim) cycle
                    zcc(mm,nn) = coagc(zdpart(mm),zdpart(nn),zmpart(mm),zmpart(nn),temppi,pressi,1)
                    zcc(nn,mm) = zcc(mm,nn)
                 END DO
              END DO
           END IF
          ! Collision-coalescence between cloud droplets
           IF (lscgcc .AND. any_cloud) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) cycle
                 DO nn = mm,ncld
                    IF (pcloud(ii,jj,nn)%numc<nlim) cycle
                    zcccc(mm,nn) = coagc(zdcloud(mm),zdcloud(nn),zmcloud(mm),zmcloud(nn),temppi,pressi,2)
                    zcccc(nn,mm) = zcccc(mm,nn)
                 END DO
              END DO
           END IF
           ! Self-collection of rain drops
           IF (lscgpp .AND. any_precp) THEN
              DO mm = 1,nprc
                 IF (pprecp(ii,jj,mm)%numc<prlim) cycle
                 DO nn = mm,nprc
                    IF (pprecp(ii,jj,nn)%numc<prlim) cycle
                    zccpp(mm,nn) =  coagc(zdprecp(mm),zdprecp(nn),zmprecp(mm),zmprecp(nn),temppi,pressi,2)
                    zccpp(nn,mm) = zccpp(mm,nn)
                 END DO
              END DO
           END IF
           ! Cloud collection of aerosols
           IF (lscgca .AND. any_cloud) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,ncld
                    IF (pcloud(ii,jj,nn)%numc<nlim) cycle
                    zccca(mm,nn) = coagc(zdpart(mm),zdcloud(nn),zmpart(mm),zmcloud(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF
           ! Collection of aerosols by rain
           IF (lscgpa .AND. any_precp) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nprc
                    IF (pprecp(ii,jj,nn)%numc<prlim) cycle
                    zccpa(mm,nn) = coagc(zdpart(mm),zdprecp(nn),zmpart(mm),zmprecp(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF
           ! Collection of cloud droplets by rain
           IF (lscgpc .AND. any_cloud .AND. any_precp) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nprc
                    IF (pprecp(ii,jj,nn)%numc<prlim) cycle
                    zccpc(mm,nn) = coagc(zdcloud(mm),zdprecp(nn),zmcloud(mm),zmprecp(nn),temppi,pressi,2)
                  END DO
              END DO
           END IF
           !  collection of aerosols by ice
           IF (lscgia .AND. any_ice) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) cycle
                    zccia(mm,nn) =  coagc(zdpart(mm),zdice(nn),zmpart(mm),zmice(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF
          !  collection of cloud particles droplets by ice
           IF (lscgic .AND. any_ice .AND. any_cloud) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) cycle
                    zccic(mm,nn) = coagc(zdcloud(mm),zdice(nn),zmcloud(mm),zmice(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF
           !  collitions between ice particles
           IF (lscgii .AND. any_ice) THEN
              DO mm = 1,nice
                 IF (pice(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = mm,nice
                    IF (pice(ii,jj,nn)%numc<prlim) CYCLE
                    zccii(mm,nn) = coagc(zdice(mm),zdice(nn),zmice(mm),zmice(nn),temppi,pressi,2)
                    zccii(nn,mm) = zcccc(mm,nn)
                 END DO
              END DO
           END IF
           !  collection of precip by ice-collision
           IF (lscgip .AND. any_precp .AND. any_ice) THEN
              DO mm = 1,nprc
                 IF (pprecp(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) CYCLE
                    zccip(mm,nn) = coagc(zdprecp(mm),zdice(nn),zmprecp(mm),zmice(nn),temppi,pressi,2)
                  END DO
              END DO
           END IF
           ! Self-collection of snow particles
           IF (lscgss .AND. any_snow) THEN
              DO mm = 1,nsnw
                 IF (psnow(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = mm,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccss(mm,nn) =  coagc(zdsnow(mm),zdsnow(nn),zmsnow(mm),zmsnow(nn),temppi,pressi,2)
                    zccss(nn,mm) = zccss(mm,nn)
                 END DO
              END DO
           END IF
           ! Collection of aerosols by snow
           IF (lscgsa .AND. any_snow) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsa(mm,nn) = coagc(zdpart(mm),zdsnow(nn),zmpart(mm),zmsnow(nn),temppi,pressi,2)
                 END DO
              END DO
           END IF
           ! collection of precip by snow
           IF (lscgsp .AND. any_precp .AND. any_snow) THEN
              DO mm = 1,nprc
                 IF (pprecp(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsp(mm,nn) = coagc(zdprecp(mm),zdsnow(nn),zmprecp(mm),zmsnow(nn),temppi,pressi,2)
                  END DO
              END DO
           END IF
           ! collection of cloud droples by snow
           IF (lscgsc .AND. any_cloud .AND. any_snow) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsc(mm,nn) = coagc(zdcloud(mm),zdsnow(nn),zmcloud(mm),zmsnow(nn),temppi,pressi,2)
                  END DO
              END DO
           END IF
           ! collection of ice by snow
           IF (lscgsi .AND. any_ice .AND. any_snow) THEN
              DO mm = 1,nice
                 IF (pice(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsi(mm,nn) = coagc(zdice(mm),zdsnow(nn),zmice(mm),zmsnow(nn),temppi,pressi,2)
                  END DO
              END DO
           END IF

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           !-- 3) New particle and volume concentrations after coagulation -------------

           ! Aerosols in regime 1a
           ! --------------------------------
           DO kk = in1a,fn1a
              IF (paero(ii,jj,kk)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.
              ! Particles lost by coagulation with larger aerosols
              DO ll = kk+1,fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc
              END DO

              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! particles lost by rain collection
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! particles lost by ice collection
              DO ll = 1,nice
                 zminusterm = zminusterm + zccia(kk,ll)*pice(ii,jj,ll)%numc
              END DO

              ! particles lost by snow collection
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsa(kk,ll)*psnow(ii,jj,ll)%numc
              END DO

              DO ll = in1a,kk-1
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:8) = ( paero(ii,jj,kk)%volc(1:8)+ptstep*zplusterm(1:8) * &
                   paero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)

           END DO

           ! Aerosols in regime 2a
           ! ---------------------------------
           DO kk = in2a,fn2a
              IF (paero(ii,jj,kk)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! Find corresponding size bin in subregime 2b
              index_2b = kk - in2a + in2b

              ! Particles lost by larger particles in 2a
              DO ll = kk+1, fn2a
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2a
              END DO

              ! Particles lost by larger particles in 2b
              DO ll = index_2b+1, fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2b
              END DO

              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by rain
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! particles lost by ice collection
              DO ll = 1,nice
                 zminusterm = zminusterm + zccia(kk,ll)*pice(ii,jj,ll)%numc
              END DO

              ! particles lost by snow collection
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsa(kk,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Particle volume gained from smaller particles in regimes 1, 2a and 2b
              DO ll = in1a, kk-1
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              ! Particle volume gained from smaller (and equal) particles in 2b
              DO ll = in2b, index_2b
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:8) = ( paero(ii,jj,kk)%volc(1:8)+ptstep*zplusterm(1:8) *  &
                   paero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)

           END DO

           ! Aerosols in regime 2b
           ! ---------------------------------
           DO kk = in2b,fn2b
              IF (paero(ii,jj,kk)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              !-- Find corresponding size bin in subregime 2a
              index_2a = kk - in2b + in2a

              ! Particles lost to larger particles in regimes 2b
              DO ll = kk+1, fn2b
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc ! 2b
              END DO

              ! Particles lost to larger and equal particles in 2a
              DO ll = index_2a, fn2a
                 zminusterm = zminusterm + zcc(kk,ll)*paero(ii,jj,ll)%numc
              END DO

              ! Particles lost by cloud collection
              DO ll = 1,ncld
                 zminusterm = zminusterm + zccca(kk,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by rain
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpa(kk,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! particles lost by ice collection
              DO ll = 1,nice
                 zminusterm = zminusterm + zccia(kk,ll)*pice(ii,jj,ll)%numc
              END DO

              ! particles lost by snow collection
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsa(kk,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Particle volume gained from smaller particles in 1/2a
              DO ll = in1a, index_2a-1
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              ! Particle volume gained from smaller particles in 2b
              DO ll = in2b, kk-1
                 zplusterm(1:8) = zplusterm(1:8) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:8)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:8) = ( paero(ii,jj,kk)%volc(1:8)+ptstep*zplusterm(1:8) *  &
                   paero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)

           END DO

           ! Cloud droplets, regime a
           ! ------------------------------------------------
           DO cc = ica%cur,fca%cur
              IF (pcloud(ii,jj,cc)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime b cloud droplets
              kk = MAX(cc-fca%cur+ncld,icb%cur) ! Regime a has more bins than b:
                                                     ! Set this at minimum to beginnign of b.

              ! Droplets lost by those with larger nucleus in regime a
              DO ll = cc+1,fca%cur
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by those with larger nucleus in regime b
              DO ll = kk+1,fcb%cur
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by rain drops
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpc(cc,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by ice particles
              DO ll = 1,nice
                 zminusterm = zminusterm + zccic(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by snow particles
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsc(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained from cloud collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller droplets in a
              DO ll = ica%cur,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller or equal droplets in b
              DO ll = icb%cur,kk
                 zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              pcloud(ii,jj,cc)%volc(1:8) = max(0.,( pcloud(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*pcloud(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pcloud(ii,jj,cc)%numc = max(0., pcloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zcccc(cc,cc)*pcloud(ii,jj,cc)%numc ) )

           END DO

           ! Cloud droplets, regime b
           ! -----------------------------------------
           DO cc = icb%cur,fcb%cur
              IF (pcloud(ii,jj,cc)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime a cloud droplets
              kk = cc - ncld + fca%cur

              ! Droplets lost by those with larger nucleus in regime b
              DO ll = cc+1,fcb%cur
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by those with larger nucleus in regime a
              DO ll = kk+1,fca%cur
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by rain drops
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccpc(cc,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by ice
              DO ll = 1,nice
                 zminusterm = zminusterm + zccic(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Droplets lost by collection by snow particles
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsc(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained from cloud collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller droplets in b
              DO ll = icb%cur,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller or equal droplets in a
              DO ll = ica%cur,kk
                 zplusterm(1:8) = zplusterm(1:8) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              pcloud(ii,jj,cc)%volc(1:8) = max(0., ( pcloud(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*pcloud(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
              pcloud(ii,jj,cc)%numc = max(0.,pcloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zcccc(cc,cc)*pcloud(ii,jj,cc)%numc ) )

           END DO

           ! Rain drops
           ! -----------------------------------
           DO cc = 1,nprc
              IF (pprecp(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! Drops lost by coagulation with larger drops
              DO ll = cc+1,nprc
                 zminusterm = zminusterm + zccpp(cc,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Drops lost by collection by snow drops
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsp(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Drops lost by collisions with ice
              DO ll = 1,nice
                 zminusterm = zminusterm + zccip(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Volume gained by collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccpa(ll,cc)*paero(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained by collection of cloud droplets
              DO ll = 1,ncld
                 zplusterm(1:8) = zplusterm(1:8) + zccpc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller drops
              DO ll = 1,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zccpp(ll,cc)*pprecp(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              pprecp(ii,jj,cc)%volc(1:8) = max(0., ( pprecp(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*pprecp(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
              pprecp(ii,jj,cc)%numc = max(0.,pprecp(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccpp(cc,cc)*pprecp(ii,jj,cc)%numc ) )

           END DO

           ! Ice particles, regime a
           ! ------------------------------------------------
           DO cc = iia%cur,fia%cur
              IF (pice(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime b cloud droplets
              kk = MAX(cc-fia%cur+nice,iib%cur) ! Regime a has more bins than b:
                                                     ! Set this at minimum to beginnign of b.

              ! Particles lost by those with larger nucleus in regime a
              DO ll = cc+1,fia%cur
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by those with larger nucleus in regime b
              DO ll = kk+1,fcb%cur
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by snow drops
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsi(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by rain drops !! huomhuom ice'n'precp
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccip(ll,cc)*pice(ii,jj,ll)%numc
              END DO

              ! Volume gained from ice collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccia(ll,cc)*paero(ii,jj,ll)%volc(1:8)*rhowa/rhoic
              END DO

              ! Volume gained from cloud collection
              DO ll = 1,ncld
                 zplusterm(1:8) = zplusterm(1:8) + zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)*rhowa/rhoic
              END DO

              ! Volume gained from smaller droplets in a
              DO ll = iia%cur,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller or equal droplets in b
              DO ll = iib%cur,kk
                 zplusterm(1:8) = zplusterm(1:8) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              pice(ii,jj,cc)%volc(1:8) = max(0., ( pice(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*pice(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pice(ii,jj,cc)%numc = max(0.,pice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccii(cc,cc)*pice(ii,jj,cc)%numc ) )

           END DO

           ! Ice particles, regime b
           ! -----------------------------------------
           DO cc = iib%cur,fib%cur
              IF (pice(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime a cloud droplets
              kk = cc - nice + fia%cur

              ! Particles lost by those with larger nucleus in regime b
              DO ll = cc+1,fib%cur
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by those with larger nucleus in regime a
              DO ll = kk+1,fia%cur
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by snow drops
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsi(cc,ll)*pprecp(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by rain drops !! huomhuom ice'n'precp
              DO ll = 1,nprc
                 zminusterm = zminusterm + zccip(ll,cc)*pice(ii,jj,ll)%numc
              END DO

              ! Volume gained from ice collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccia(ll,cc)*paero(ii,jj,ll)%volc(1:8)*rhowa/rhoic
              END DO

              ! Volume gained from cloud collection
              DO ll = 1,ncld
                 zplusterm(1:8) = zplusterm(1:8) + zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)*rhowa/rhoic
              END DO

              ! Volume gained from smaller droplets in b
              DO ll = iib%cur,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:8)
              END DO

              ! Volume gained from smaller or equal droplets in a
              DO ll = iia%cur,kk
                 zplusterm(1:8) = zplusterm(1:8) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              pice(ii,jj,cc)%volc(1:8) = max(0.,( pice(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*pice(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
              pice(ii,jj,cc)%numc = max(0.,pice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccii(cc,cc)*pice(ii,jj,cc)%numc ) )

           END DO

           ! Snow
           ! -----------------------------------
           DO cc = 1,nsnw
              IF (psnow(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! Drops lost by coagulation with larger snow drops
              DO ll = cc+1,nsnw
                 zminusterm = zminusterm + zccss(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained by collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:8) = zplusterm(1:8) + zccsa(ll,cc)*paero(ii,jj,ll)%volc(1:8)*rhowa/rhosn
              END DO

              ! Volume gained by collection of cloud droplets
              DO ll = 1,ncld
                 zplusterm(1:8) = zplusterm(1:8) + zccsc(ll,cc)*pcloud(ii,jj,ll)%volc(1:8)*rhowa/rhosn
              END DO

              ! Volume gained by collection of ice particles
              DO ll = 1,nice
                 zplusterm(1:8) = zplusterm(1:8) + zccsi(ll,cc)*pice(ii,jj,ll)%volc(1:8)*rhoic/rhosn
              END DO

              ! Volume gained by collection of rain drops
              DO ll = 1,nprc
                 zplusterm(1:8) = zplusterm(1:8) + zccsp(ll,cc)*pprecp(ii,jj,ll)%volc(1:8)*rhowa/rhosn
              END DO

              ! Volume gained by collisions between ice and rain
              nn = min(nice,nprc)
              DO ll = 1,nn
                 zplusterm(1:8) = zplusterm(1:8) + zccip(ll,cc)*pice(ii,jj,ll)%volc(1:8)*rhoic/rhosn* &
                                                                  pprecp(ii,jj,ll)%volc(1:8)*rhowa/rhosn
              END DO

              ! Volume gained from smaller drops
              DO ll = 1,cc-1
                 zplusterm(1:8) = zplusterm(1:8) + zccss(ll,cc)*psnow(ii,jj,ll)%volc(1:8)
              END DO

              ! Update the hydrometeor volume concentrations
              psnow(ii,jj,cc)%volc(1:8) = max(0.,( psnow(ii,jj,cc)%volc(1:8) +  &
                   ptstep*zplusterm(1:8)*psnow(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              psnow(ii,jj,cc)%numc = max(0.,psnow(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccss(cc,cc)*psnow(ii,jj,cc)%numc ) )

           END DO

      END DO ! kbdim
     END DO ! klev

  END SUBROUTINE coagulation


  ! fxm: calculated for empty bins too
  ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
  !      and organic vapours (average values? 'real' values for each?)
  !********************************************************************
  !
  ! subroutine CONDENSATION(kproma, kbdim,  klev,        &
  !                         pnaero, pvols,  pdwet, plwc, &
  !                         pcsa,   pcocnv, pcocsv,      &
  !                         ptemp,  ppres,  ptstep)
  !
  !********************************************************************
  !
  ! Purpose:
  ! --------
  ! Calculates the increase in particle volume and
  !  decrease in gas phase concentrations due to condensation
  !  of sulphuric acid and two organic compounds (non-volatile
  !  and semivolatile)
  !
  !
  ! Method:
  ! -------
  ! Regime 3 particles only act as a sink for condensing vapours
  !  while their size and composition does not change.
  ! Exception: Soluble fraction of regime 3c particles can change
  !  and thus they can be moved to regime 3b
  !
  ! New gas and aerosol phase concentrations calculated according
  !  to Jacobson (1997): Numerical techniques to solve
  !  condensational and dissolutional growth equations
  !  when growth is coupled to reversible reactions,
  !  Aerosol Sci. Tech., 27, pp 491-498.
  !
  ! fxm: one should really couple with vapour production and loss terms as well
  !      should nucleation be coupled here as well????
  !
  ! Juha: Now does the condensation of water vapour on hydrometeors as well,
  !       + the condensation of semivolatile aerosol species on hydromets.
  !       Modified for the new aerosol datatype. LWC is obtained from %volc(8)
  !
  !
  ! Interface:
  ! ----------
  ! Called from main aerosol model
  !
  !
  ! Coded by:
  ! ---------
  ! Hannele Korhonen (FMI) 2005
  ! Harri Kokkola (FMI) 2006
  ! Juha Tonttila (FMI) 2014
  !
  !---------------------------------------------------------------
  !
  ! Following parameterization has been used:
  ! ------------------------------------------
  !
  ! Molecular diffusion coefficient of condensing vapour [m2/s]
  !  (Reid et al. (1987): Properties of gases and liquids,
  !   McGraw-Hill, New York.)
  !
  ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
  !  {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
  !
  ! M_air = 28.965 : molar mass of air [g/mol]
  ! d_air = 19.70  : diffusion volume of air
  ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
  ! d_h2so4 = 51.96  : diffusion volume of h2so4
  !
  !---------------------------------------------------------------

  SUBROUTINE condensation(kproma,  kbdim,  klev,   krow,      &
                          paero,   pcloud, pprecp,            &
                          pice,    psnow,                     &
                          pcsa,                               &
                          pcocnv,  pcocsv, pchno3, pcnh3,     &
                          prv,prs, prsi,ptemp,  ppres,  ptstep,    &
                          ppbl,    prtcl)

    USE mo_salsa_nucleation

    USE mo_submctl,    ONLY :   &
         t_section,                 & ! Data type for the cloud bin representation
         fn2b,                      &
         ncld,nprc,                  &
         nice,nsnw,                 &
         lscndgas,                  & 
         nlcndh2oae, nlcndh2ocl, nlcndh2oic, & ! Condensation to aerosols, clouds and ice particles
         nsnucl                     ! nucleation

    USE class_componentIndex, ONLY : ComponentIndex,IsUsed

    IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

    REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep,                    & ! timestep [s]
         prs(kbdim,klev),           & ! Water vapor saturation mixing ratio
         prsi(kbdim,klev)              ! Saturation mixing ratio    [kg/m3]

    TYPE(ComponentIndex), INTENT(in) :: prtcl  ! Keeps track which substances are used

    INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

    REAL, INTENT(INOUT) ::     &
         prv(kbdim,klev),          & ! Water vapor mixing ratio
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
         pcnh3(kbdim,klev)           ! ammonia concentration [#/m3]

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     & ! Hydrometeor properties
         paero(kbdim,klev,fn2b),      & ! Aerosol properties
         pprecp(kbdim,klev,nprc),     & ! rain properties
         pice(kbdim,klev,nice),       & ! ice properties
         psnow(kbdim,klev,nsnw)        ! snowing properties

    REAL :: zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles
         zxocnv(kbdim,klev), &
         zrh(kbdim,klev)

    zxocnv = 0.
    zxsa = 0.
    zj3n3 = 0.
    zrh(1:kbdim,:) = prv(1:kbdim,:)/prs(1:kbdim,:)


    !------------------------------------------------------------------------------

    ! Nucleation
    IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,   &
                                    paero,  ptemp,  zrh,    ppres,  &
                                    pcsa,   pcocnv, ptstep, zj3n3,  &
                                    zxsa,   zxocnv, ppbl            )

    ! Condensation of H2SO4 and organic vapors
    IF (lscndgas) CALL condgas(kproma,  kbdim,  klev,   krow,      &
                          paero,   pcloud, pprecp,            &
                          pice,    psnow,                     &
                          pcsa, pcocnv, pcocsv, pchno3, pcnh3,     &
                          zxsa, prv,prs, prsi,ptemp,  ppres, ptstep,    &
                          prtcl)

    ! Condensation of water vapour
    IF (nlcndh2ocl .OR. nlcndh2oae .OR. nlcndh2oic) &
        CALL gpparth2o(kproma,kbdim,klev,krow,  &
                   paero, pcloud, pprecp,   &
                   pice, psnow,             & ! ice'n'snow
                   ptemp,ppres,prs,prsi,prv,     &
                   ptstep)

    ! HNO3/NH3 - currently disabled
    !CALL gpparthno3(kproma,kbdim,klev,krow,ppres,ptemp,paero,pcloud,   &
    !                pprecp,pchno3,pcnh3,prv,prs,zbeta,ptstep           )

  END SUBROUTINE condensation

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE condgas(kproma,  kbdim,  klev,   krow,      &
                          paero,   pcloud, pprecp,            &
                          pice,    psnow,                     &
                          pcsa,                               &
                          pcocnv,  pcocsv, pchno3, pcnh3,     &
                          zxsa, prv,prs, prsi,ptemp,  ppres,  ptstep,    &
                          prtcl)

    USE mo_submctl,    ONLY :   &
         pi,                        &
         pi6,                       & ! pi/6
         in1a, in2a,                & ! size bin indices
         fn2b,                &
         t_section,                 & ! Data type for the cloud bin representation
         ncld,                      &
         nprc,                      &
         nice,                      & ! ice'n'snow
         nsnw,                      & ! ice'n'snow
         avog,                      &
         nlim,                      &
         prlim,                     &
         rhowa,                     & ! density of water (kg/m3)
         rhoic,                     & ! density of ice (kg/m3)
         rhosn,                     & ! density of snow (kg/m3)
         rhosu,                     & ! density of sulphate (kg/m3)
         rhooc,                     & ! density of organic carbon (kg/m3)
         rhoss,                     & ! density of sea salt (kg/m3)
         rhono,                     & ! density of nitric acid (kg/m3)
         rhonh,                     & ! density of ammonia (kg/m3)
         rhobc,                     &
         rhodu,                      &
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         msu,                       & ! molar mass of sulphate [kg/mol]
         moc,                       & !       "       organic carbon
         mss,                       & !       "       sea salt
         mno,                       & !       "       nitrate
         mnh,                       & !       "       ammonium
         mbc,                       & !
         mdu,                       &
         mwa,                       & !               water
         mair,                      & !       "       air
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [m3]
         mvnh, mvno, mvwa,          & ! molecular volumes of HNO3 and NH3,H20 [m3]
         d_sa,                      & ! diameter of H2SO4 molecule [m]

         epsoc,                     & ! soluble fraction of organics (scaled to sulphate)
         massacc,                   & ! mass accomodation coefficients in each bin
         n3,                        & ! number of molecules in one 3 nm particle [1]
         surfw0,                    & ! surface tension of water
         surfi0                       ! surface tension of ice

    USE class_componentIndex, ONLY : ComponentIndex,IsUsed

    USE mo_constants,      ONLY: g, avo, alv, rv, als
   IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) ::          &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

    REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep,                    & ! timestep [s]
         prs(kbdim,klev),           & ! Water vapor saturation mixing ratio
         prsi(kbdim,klev)              ! Saturation mixing ratio    [kg/m3]

    TYPE(ComponentIndex), INTENT(in) :: prtcl  ! Keeps track which substances are used

    REAL, INTENT(INOUT) ::     &
         prv(kbdim,klev),          & ! Water vapor mixing ratio
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
         pcnh3(kbdim,klev),        & ! ammonia concentration [#/m3]
         zxsa(kbdim,klev)            ! ratio of sulphuric acid and organic vapor in 3nm particles

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     & ! Hydrometeor properties
         paero(kbdim,klev,fn2b),      & ! Aerosol properties
         pprecp(kbdim,klev,nprc),     & ! rain properties
         pice(kbdim,klev,nice),       & ! ice properties
         psnow(kbdim,klev,nsnw)        ! snowing properties


    !-- Local variables ----------------------
    INTEGER :: ii, jj    ! loop indices

    REAL ::                      &
         zvisc,                      & ! viscosity of air [kg/(m s)]
         zdfvap,                     & ! air diffusion coefficient [m2/s]
         zmfp,                       & ! mean free path of condensing vapour [m]
         zcs_tot,                    & ! total condensation sink [1/s] (gases)
         zcs_ocsv,                   & ! condensation sink for semivolatile organics [1/s]
         zcs_su,                     & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                   & ! condensation sink for nonvolatile organics [1/s]
                                       ! vapour concentration after time step [#/m3]
         zcvap_new1,                 & ! sulphuric acid
         zcvap_new2,                 & ! nonvolatile organics
         zcvap_new3,                 & ! semivolatile organics
                                       ! change in vapour concentration [#/m3]
         zdvap1,                     & ! sulphuric acid
         zdvap2,                     & ! nonvolatile organics
         zdvap3,                     & ! semivolatile organics

         zdfpart(in1a+1),            & ! particle diffusion coefficient

         zknud(fn2b),                & ! particle Knudsen number
         zknca(ncld),                & ! Knudsen number for cloud droplets and aerosol vapours
         zknpa(nprc),                & ! Knudsen number for rain drops and aerosol vapours
         zknia(nice),                & ! Knudsen number for ice particles and aerosol vapours
         zknsa(nsnw),                & ! Knudsen number for snow flakes and aerosol vapours

         zbeta(fn2b),                & ! transitional correction factor for aerosols
         zbetaca(ncld),              & ! - '' - for condensing aerosol vapours on clouds (is this needed?)
         zbetapa(nprc),              & ! - '' - for condensing aerosol vapours on rain drops
         zbetaia(nice),              & ! - '' - for condensing aerosol vapours on ice (is this needed?)
         zbetasa(nsnw),              & ! - '' - for condensing aerosol vapours on snow flakes

         zcolrate(fn2b),             & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),        & ! collision rate of organic molecules to particles [1/s]
         zcolrateca(ncld),           & ! Collision rate of aerosol vapour molecules to cloud drops
         zcolratepa(nprc),           & ! Collision rate of gases to rain drops
         zcolrateia(nice),           & ! Collision rate of aerosol vapour molecules to ice particles
         zcolratesa(nsnw),           & ! Collision rate of gases to rain drops

         zdvolsa(fn2b),              & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),              & !    - " - organics

         zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                    & ! ratio of nucleation of all mass transfer in the smallest bin
         zxocnv(kbdim,klev)


    zj3n3 = 0.
    zxocnv = 0.

    zdvolsa=0.
    zn_vs_c=0.
    DO jj = 1,klev
       DO ii = 1,kbdim

          zdvoloc = 0.

          !-- 1) Properties of air and condensing gases --------------------
          zvisc  = (7.44523e-3*ptemp(ii,jj)**1.5)/(5093.*(ptemp(ii,jj)+110.4))! viscosity of air [kg/(m s)]
          zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
          zmfp   = 3.*zdfvap*sqrt(pi*msu/(8.*rg*ptemp(ii,jj)))                      ! mean free path [m]

          !-- 2) Transition regime correction factor for particles ---------
          !
          !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          !  Topics in current aerosol research, Pergamon.
          !
          !  Size of condensing molecule considered only for
          !  nucleation mode (3 - 20 nm)
          !

          !-- particle Knudsen numbers
          zknud(in1a:in1a+1) = 2.*zmfp/(paero(ii,jj,in1a:in1a+1)%dwet+d_sa)              ! Gases on aerosols
          zknud(in1a+2:fn2b) = 2.*zmfp/paero(ii,jj,in1a+2:fn2b)%dwet

          zknca(1:ncld) = 2.*zmfp/pcloud(ii,jj,1:ncld)%dwet          ! Knudsen number for gases on cloud drplets

          zknpa(1:nprc) = 2.*zmfp/pprecp(ii,jj,1:nprc)%dwet          ! Knudsen number for gases on rain drops

          zknia(1:nice) = 2.*zmfp/pice(ii,jj,1:nice)%dwet          ! Knudsen number for gases on ice particles

          zknsa(1:nsnw) = 2.*zmfp/psnow(ii,jj,1:nsnw)%dwet          ! Knudsen number for gases on snow flakes

          !-- transitional correction factor
          zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &     ! Aerosol + gas
                  (3.*massacc)*(zknud+zknud**2))

          zbetaca = 1. + zknca*( 1.33 + (0.71/zknca) )/( 1. + (1./zknca) ) ! Hydrometeor + gas
          zbetaca = 1./zbetaca

          zbetapa = 1. + zknpa*( 1.33 + (0.71/zknpa) )/( 1. + (1./zknpa) ) ! Rain drop + gas
          zbetapa = 1./zbetapa

          zbetaia = 1. + zknia*( 1.33 + (0.71/zknia) )/( 1. + (1./zknia) ) ! ice + gas
          zbetaia = 1./zbetaia

          zbetasa = 1. + zknsa*( 1.33 + (0.71/zknsa) )/( 1. + (1./zknsa) ) ! Rain drop + gas
          zbetasa = 1./zbetasa
          !-- 3) Collision rate of molecules to particles -------------------
          !
          !  Particle diffusion coefficient considered only for
          !  nucleation mode (3 - 20 nm)
          !

          !-- particle diffusion coefficient [m2/s]
          zdfpart = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &
                    (3.*pi*zvisc*paero(ii,jj,in1a:in1a+1)%dwet)

          !-- collision rate (gases on aerosols) [1/s]
          zcolrate = 0.
          zcolrate(in1a:in1a+1) = MERGE( 2.*pi*(paero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    &
                                              (zdfvap+zdfpart)*zbeta(in1a:in1a+1)*              &
                                              paero(ii,jj,in1a:in1a+1)%numc,                    &
                                              0.,                                            &
                                              paero(ii,jj,in1a:in1a+1)%numc > nlim              )

          zcolrate(in1a+2:fn2b) = MERGE( 2.*pi*paero(ii,jj,in1a+2:fn2b)%dwet*zdfvap*    &
                                              zbeta(in1a+2:fn2b)*paero(ii,jj,in1a+2:fn2b)%numc, &
                                              0.,                                            &
                                              paero(ii,jj,in1a+2:fn2b)%numc > nlim              )

          !-- gases on hydrometeors
          zcolrateca = 0.
          zcolrateca(1:ncld) = MERGE( 2.*pi*pcloud(ii,jj,1:ncld)%dwet*zdfvap*    &
                                           zbetaca(1:ncld)*pcloud(ii,jj,1:ncld)%numc,    &
                                           0.,                                        &
                                           pcloud(ii,jj,1:ncld)%numc > nlim              )

          ! Gases on rain drops
          zcolratepa = 0.
          zcolratepa(1:nprc) = MERGE( 2.*pi*pprecp(ii,jj,1:nprc)%dwet*zdfvap*    &
                                           zbetapa(1:nprc)*pprecp(ii,jj,1:nprc)%numc,    &
                                           0.,                                        &
                                           pprecp(ii,jj,1:nprc)%numc > prlim             )
          !-- gases on ice particles
          zcolrateia = 0.
          zcolrateia(1:nice) = MERGE( 2.*pi*pice(ii,jj,1:nice)%dwet*zdfvap*    &
                                           zbetaia(1:nice)*pice(ii,jj,1:nice)%numc,    &
                                           0.,                                        &
                                           pice(ii,jj,1:ncld)%numc > prlim              )

          ! Gases on snow flakes
          zcolratesa = 0.
          zcolratesa(1:nsnw) = MERGE( 2.*pi*psnow(ii,jj,1:nsnw)%dwet*zdfvap*    &
                                           zbetasa(1:nsnw)*psnow(ii,jj,1:nsnw)%numc,    &
                                           0.,                                        &
                                           psnow(ii,jj,1:nsnw)%numc > prlim             )

          !-- 4) Condensation sink [1/s] -------------------------------------

          zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)+ sum(zcolrateia) + sum(zcolratesa)  ! total sink

          !-- 5) Changes in gas-phase concentrations and particle volume -----
          !
          !--- 5.1) Organic vapours ------------------------

          !---- 5.1.1) Non-volatile organic compound: condenses onto all bins
          IF(pcocnv(ii,jj) > 1.e-10 .and. zcs_tot > 1.e-30 .AND. IsUsed(prtcl,'OC')) THEN

             zn_vs_c = 0.

             IF(zj3n3(ii,jj,2) > 1.) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                                                  pcocnv(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate_ocnv = zcolrate
             zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pcocnv(ii,jj)

             zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pcocnv(ii,jj)                   ! total sink for organic vapor

             zcvap_new2 = pcocnv(ii,jj)/(1.+ptstep*zcs_ocnv)                  ! new gas phase concentration [#/m3]
             zdvap2 = pcocnv(ii,jj) - zcvap_new2                                 ! change in gas concentration [#/m3]
             pcocnv(ii,jj) = zcvap_new2                                          ! updating vapour concentration [#/m3]

             zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*mvoc*zdvap2             ! volume change of particles
                                                                                 !  [m3(OC)/m3(air)]

             paero(ii,jj,in1a:fn2b)%volc(2) = paero(ii,jj,in1a:fn2b)%volc(2) + & !-- change of volume
                                                    zdvoloc                      !   due to condensation in 1a-2b

             ! Condensation on hydromets
             pcloud(ii,jj,1:ncld)%volc(2) = pcloud(ii,jj,1:ncld)%volc(2) +  &
                  zcolrateca(1:ncld)/zcs_ocnv*mvoc*zdvap2

             ! Condensation on rain drops
             pprecp(ii,jj,1:nprc)%volc(2) = pprecp(ii,jj,1:nprc)%volc(2) +  &
                  zcolratepa(1:nprc)/zcs_ocnv*mvoc*zdvap2

             ! Condensation on ice particles
             pice(ii,jj,1:nice)%volc(2) = pice(ii,jj,1:nice)%volc(2) +  &
                  zcolrateia(1:nice)/zcs_ocnv*mvoc*zdvap2

             ! Condensation on snow
             psnow(ii,jj,1:nsnw)%volc(2) = psnow(ii,jj,1:nsnw)%volc(2) +  &
                  zcolratesa(1:nsnw)/zcs_ocnv*mvoc*zdvap2

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
             ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
             IF (zxocnv(ii,jj) > 0.) THEN
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
             END IF

          END IF


          !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
          zcs_ocsv = sum(zcolrate(in2a:fn2b)) +  &       ! sink for semivolatile organics
                     sum(zcolrateca(1:ncld))  +  &       ! ... including condensation on cloud droplets
                     sum(zcolratepa(1:nprc))  +  &       ! and rain drops
                     sum(zcolrateia(1:nice))  +  &       ! and ice particles
                     sum(zcolratesa(1:nsnw))             ! and snow particles

          IF(pcocsv(ii,jj) > 1.e-10 .and. zcs_ocsv > 1.e-30 .AND. IsUsed(prtcl,'OC')) THEN


             zcvap_new3 = pcocsv(ii,jj)/(1.+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
             zdvap3 = pcocsv(ii,jj) - zcvap_new3                  ! change in gas concentration [#/m3]
             pcocsv(ii,jj) = zcvap_new3                           ! updating gas concentration [#/m3]

             zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &          ! volume change of particles
                  zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3        !  [m3(OC)/m3(air)]

             paero(ii,jj,in1a:fn2b)%volc(2) = &                   !-- change of volume due
                  paero(ii,jj,in1a:fn2b)%volc(2) + zdvoloc        !   due to condensation in 1a-2b

             ! Condensation on hydromets
             pcloud(ii,jj,1:ncld)%volc(2) = pcloud(ii,jj,1:ncld)%volc(2)  +  &
                  zcolrateca(1:ncld)/zcs_ocsv*mvoc*zdvap3

             ! Condensation on rain drops
             pprecp(ii,jj,1:nprc)%volc(2) = pprecp(ii,jj,1:nprc)%volc(2)  +  &
                  zcolratepa(1:nprc)/zcs_ocsv*mvoc*zdvap3

             ! Condensation on ice particles
             pice(ii,jj,1:nice)%volc(2) = pice(ii,jj,1:nice)%volc(2)  +  &
                  zcolrateia(1:nice)/zcs_ocsv*mvoc*zdvap3

             ! Condensation on snow particles
             psnow(ii,jj,1:nsnw)%volc(2) = psnow(ii,jj,1:nsnw)%volc(2)  +  &
                  zcolratesa(1:nprc)/zcs_ocsv*mvoc*zdvap3

          END IF


          ! ---- 5.2) Sulphate -------------------------------------------
          IF(pcsa(ii,jj) > 1.e-10 .and. zcs_tot > 1.e-30 .AND. IsUsed(prtcl,'SO4')) THEN

             !-- Ratio of mass transfer between nucleation and condensation

             zn_vs_c = 0.

             IF(zj3n3(ii,jj,1) > 1.) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                              (zj3n3(ii,jj,1) +  &
                                              pcsa(ii,jj) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pcsa(ii,jj)

             zcs_su = zcs_tot + zj3n3(ii,jj,1) / pcsa(ii,jj)      ! total sink for sulfate

             !--- Sulphuric acid -------------------------
             !
             zcvap_new1 = pcsa(ii,jj) /(1.+ptstep*zcs_su)         ! new gas phase concentration [#/m3]
             zdvap1 = pcsa(ii,jj) - zcvap_new1                    ! change in gas concentration [#/m3]
             pcsa(ii,jj) = zcvap_new1                             ! updating vapour concentration [#/m3]

             zdvolsa = zcolrate(in1a:fn2b)/zcs_su*mvsu*zdvap1     ! volume change of particles
             ! [m3(SO4)/m3(air)] by condensation

             !-- Change of volume concentration of sulphate in aerosol [fxm]
             paero(ii,jj,in1a:fn2b)%volc(1) = paero(ii,jj,in1a:fn2b)%volc(1) + zdvolsa

             !-- Clouds
             pcloud(ii,jj,1:ncld)%volc(1) = pcloud(ii,jj,1:ncld)%volc(1)  +  &
                  zcolrateca(1:ncld)/zcs_su*mvsu*zdvap1

             ! Rain drops
             pprecp(ii,jj,1:nprc)%volc(1) = pprecp(ii,jj,1:nprc)%volc(1)  +  &
                  zcolratepa(1:nprc)/zcs_su*mvsu*zdvap1

             !-- Ice clouds
             pice(ii,jj,1:nice)%volc(1) = pice(ii,jj,1:nice)%volc(1)  +  &
                  zcolrateia(1:nice)/zcs_su*mvsu*zdvap1

             ! Snow particles
             psnow(ii,jj,1:nsnw)%volc(1) = psnow(ii,jj,1:nsnw)%volc(1)  +  &
                  zcolratesa(1:nsnw)/zcs_su*mvsu*zdvap1

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxsa(ii,jj) > 0.) THEN
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc +          &
                     zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
             END IF

          END IF

       END DO ! kproma

    END DO ! klev

  END SUBROUTINE condgas

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE gpparth2o(kproma, kbdim,  klev, krow,  &
                       paero,  pcloud, pprecp,      &
                       pice, psnow,                 &
                       ptemp,  ppres,  prs,prsi, prv,    &
                       ptstep)
    
    USE mo_submctl, ONLY : t_section,            &
                               nbins, ncld, nprc,    &
                               nice, nsnw,            &
                               rhowa, rhoic, rhosn,mwa, mair,     &
                               surfw0,surfi0, rg,           &
                               pi, prlim, nlim,      &
                               massacc,avog,pstand,  &
                               in1a,in2a,  &
                               fn2b,            &
                               lscndh2oae, lscndh2ocl, lscndh2oic
    USE mo_constants, ONLY : alv, als  ! ice'n'snow
    USE mo_salsa_properties, ONLY : equilibration
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),  &
                                      pcloud(kbdim,klev,ncld),  &
                                      pprecp(kbdim,klev,nprc),  &
                                      pice(kbdim,klev,nice),    & ! ice'n'snow
                                      psnow(kbdim,klev,nsnw)      ! ice'n'snow

    REAL, INTENT(inout) :: prv(kbdim,klev)

    REAL :: zkelvin(nbins), zkelvincd(ncld), zkelvinpd(nprc), &  ! Kelvin effects
                zkelvinid(nice), zkelvinsd(nsnw)                      ! Kelvin effects ice'n'snow
    REAL :: zka(nbins), zkacd(ncld), zkapd(nprc),        &  ! Activity coefficients
                zkaid(nice), zkasd(nsnw)                        ! Activity coefficients ! ice'n'snow
    REAL :: zcwsurfae(nbins), zcwsurfcd(ncld), zcwsurfpd(nprc), & ! Surface mole concentrations
                zcwsurfid(nice), zcwsurfsd(nsnw)                 ! surface mole concentrations ice'n'snow
    REAL :: zmtae(nbins), zmtcd(ncld), zmtpd(nprc),      & ! Mass transfer coefficients
                zmtid(nice), zmtsd(nsnw)
    REAL :: zwsatae(nbins), zwsatcd(ncld), zwsatpd(nprc), &  ! Water saturation ratios above
                zwsatid(nice), zwsatsd(nsnw)                      ! ice'n'snow
    REAL :: zcwtot                                        ! Total water mole concentration
    REAL :: zcwc, zcwn, zcwint                            ! Current and new water vapour mole concentrations
    REAL :: zcwcae(nbins), zcwnae(nbins), zcwintae(nbins) ! Current and new water mole concentrations in aerosols
    REAL :: zcwccd(ncld), zcwncd(ncld), zcwintcd(ncld)    !     -  ''  -     in cloud droplets
    REAL :: zcwcpd(nprc), zcwnpd(nprc), zcwintpd(nprc)    !     -  ''  -     in rain drops
    REAL :: zcwcid(nice), zcwnid(nice), zcwintid(nice)    !     -  ''  -     in ice particles
    REAL :: zcwcsd(nsnw), zcwnsd(nsnw), zcwintsd(nsnw)    !     -  ''  -     in snow particles
    REAL :: zdfh2o, zthcond,rhoair
    REAL :: zbeta,zknud,zmfph2o
    REAL :: zact, zhlp1,zhlp2,zhlp3
    REAL :: adt,adtc(nbins),ttot
    REAL :: testi(nbins)

    REAL :: zrh(kbdim,klev)

    REAL :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

    INTEGER :: nstr
    INTEGER :: ii,jj,cc


    zhlp1 = 0.
    zrh(:,:) = prv(:,:)/prs(:,:)

    ! Calculate the condensation only for 2a/2b aerosol bins
    nstr = in2a

    ! Save the current aerosol water content
    zaelwc1(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(8),DIM=3)*rhowa

    ! For 1a bins do the equilibrium calculation
    CALL equilibration(kproma,kbdim,klev,      &
                       zrh,ptemp,paero,.FALSE. )

    ! If RH < 98 % OR dynamic condensation for aerosols switched off, do equilibrium for all bins
    IF (zrh(1,1) < 0.98 .OR. .NOT. lscndh2oae)  CALL equilibration(kproma,kbdim,klev,      &
                                                                      zrh,ptemp,paero,.TRUE. )

    ! The new aerosol water content after equilibrium calculation
    zaelwc2(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(8),DIM=3)*rhowa

    prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )*ppres(:,:)*mair/(rg*ptemp(:,:))

    adtc(:) = 0.
    zcwc = 0.; zcwint = 0.; zcwn = 0.
    zcwcae = 0.; zcwccd = 0.; zcwcpd = 0.; zcwcid = 0.; zcwcsd = 0.;
    zcwintae = 0.; zcwintcd = 0.; zcwintpd = 0.; zcwintid = 0.; zcwintsd = 0.
    zcwnae = 0.; zcwncd = 0.; zcwnpd = 0.; zcwnid = 0.; zcwnsd = 0.
    zwsatae = 0.; zwsatcd = 0.; zwsatpd = 0.; zwsatid = 0.; zwsatsd = 0.

    DO jj = 1,klev
       DO ii = 1,kbdim

          rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

          zdfh2o = ( 5./(16.*avog*rhoair*1.e-3*(3.11e-8)**2) ) * &
                   SQRT( rg*1e7*ptemp(ii,jj)*mair*1.e3*(mwa+mair)*1.e3/( 2.*pi*mwa*1.e3 ) )
          zdfh2o = zdfh2o*1.e-4

          zmfph2o = 3.*zdfh2o*sqrt(pi*mwa/(8.*rg*ptemp(ii,jj)))
          zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air

          ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
          zkelvinpd = 1.; zkelvincd = 1.; zkelvin = 1.
          zka = 1.; zkacd = 1.; zkapd = 1.; zkaid = 1.; zkasd = 1. ! Assume activity coefficients as 1 for now.

          ! Kelvin effects
          zkelvin(1:nbins) = exp( 4.*surfw0*mwa /  &
               (rg*ptemp(ii,jj)*rhowa*paero(ii,jj,1:nbins)%dwet) )

          zkelvincd(1:ncld) = exp( 4.*surfw0*mwa /  &
               (rg*ptemp(ii,jj)*rhowa*pcloud(ii,jj,1:ncld)%dwet) )

          zkelvinpd(1:nprc) = exp( 4.*surfw0*mwa /  &
               (rg*ptemp(ii,jj)*rhowa*MIN(pprecp(ii,jj,1:nprc)%dwet,2.e-3)) )

          zkelvinid(1:nice) = exp( 4.*surfi0*mwa /  &          ! ice surface tension
               (rg*ptemp(ii,jj)*rhoic*pice(ii,jj,1:nice)%dwet) )

          zkelvinsd(1:nsnw) = exp( 4.*surfi0*mwa /  &
               (rg*ptemp(ii,jj)*rhosn*MIN(psnow(ii,jj,1:nsnw)%dwet,2.e-3)) ) !! huomhuom onko MIN-lauseke tarpeellinen, plus tarkista tiheys

          ! Cloud droplets --------------------------------------------------------------------------------
          zmtcd(:) = 0.
          zcwsurfcd(:) = 0.
          DO cc = 1,ncld
             IF (pcloud(ii,jj,cc)%numc > nlim .AND. lscndh2ocl) THEN

                ! Activity + Kelvin effect
                zact = acth2o(pcloud(ii,jj,cc))

                ! Saturation mole concentration over flat surface
                zcwsurfcd(cc) = prs(ii,jj)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatcd(cc) = zact*zkelvincd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/pcloud(ii,jj,cc)%dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pcloud(ii,jj,cc)%numc*2.*pi*pcloud(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Rain drops --------------------------------------------------------------------------------
          zmtpd(:) = 0.
          zcwsurfpd(:) = 0.
          DO cc = 1,nprc
             IF (pprecp(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN

                ! Activity + Kelvin effect
                zact = acth2o(pprecp(ii,jj,cc))

                ! Saturation mole concentrations over flat surface
                zcwsurfpd(cc) = prs(ii,jj)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatpd(cc) = zact*zkelvinpd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/pprecp(ii,jj,cc)%dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pprecp(ii,jj,cc)%numc*2.*pi*pprecp(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd(cc)/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Ice particles --------------------------------------------------------------------------------
          zmtid(:) = 0.
          zcwsurfid(:) = 0.
          DO cc = 1,nice
             IF (pice(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN

                ! Activity + Kelvin effect
                zact = acth2o(pice(ii,jj,cc))


                ! Saturation mole concentration over flat surface
                zcwsurfid(cc) = prsi(ii,jj)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatid(cc) = zact*zkelvinid(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/pice(ii,jj,cc)%dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pice(ii,jj,cc)%numc*2.*pi*pice(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*als*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj)) !! huomhuom als
                zhlp3 = ( (als*mwa)/(rg*ptemp(ii,jj)) ) - 1. !! huomhuom als

                zmtid(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Snow particles --------------------------------------------------------------------------------
          zmtsd(:) = 0.
          zcwsurfsd(:) = 0.
          DO cc = 1,nsnw
             IF (psnow(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN

                ! Activity + Kelvin effect
                zact = acth2o(psnow(ii,jj,cc))


                ! Saturation mole concentrations over flat surface
                zcwsurfsd(cc) = prsi(ii,jj)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatsd(cc) = zact*zkelvinsd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/psnow(ii,jj,cc)%dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = psnow(ii,jj,cc)%numc*2.*pi*psnow(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*als*zwsatsd(cc)*zcwsurfsd(cc)/(zthcond*ptemp(ii,jj)) !! huomhuom als
                zhlp3 = ( (als*mwa)/(rg*ptemp(ii,jj)) ) - 1. !! huomhuom als

                zmtsd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! -- Aerosols: ------------------------------------------------------------------------------------
          zmtae(:) = 0.
          zcwsurfae(:) = 0.
          DO cc = 1,nbins
             IF (paero(ii,jj,cc)%numc > nlim .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) THEN

                ! Water activity
                zact = acth2o(paero(ii,jj,cc))
                testi(cc) = zact
                ! Saturation mole concentration over flat surface
                ! Limit the supersaturation to max 1.01 for the mass transfer
                ! EXPERIMENTAL
                zcwsurfae(cc) = MAX(prs(ii,jj),prv(ii,jj)/1.01)*rhoair/mwa

                ! Equilibrium saturation ratio
                zwsatae(cc) = zact*zkelvin(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/paero(ii,jj,cc)%dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.*massacc(cc))*(zknud+zknud**2))

                ! Mass transfer
                zhlp1 = paero(ii,jj,cc)%numc*2.*pi*paero(ii,jj,cc)%dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae(cc)/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Current mole concentrations
          zcwc = prv(ii,jj)*rhoair/mwa
          zcwcae(1:nbins) = paero(ii,jj,1:nbins)%volc(8)*rhowa/mwa
          zcwccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(8)*rhowa/mwa
          zcwcpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(8)*rhowa/mwa
          zcwcid(1:nice) = pice(ii,jj,1:nice)%volc(8)*rhoic/mwa !! ice'n'snow
          zcwcsd(1:nsnw) = psnow(ii,jj,1:nsnw)%volc(8)*rhosn/mwa !! ice'n'snow

          zcwtot = zcwc + SUM(zcwcae) + &
                          SUM(zcwccd) + &
                          SUM(zcwcpd) + &
                          SUM(zcwcid) + &
                          SUM(zcwcsd)
          ttot = 0.
          adtc = 0.

          zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd; zcwintid = zcwcid; zcwintsd = zcwcsd

          ! Substepping loop
          ! ---------------------------------
          zcwint = 0.
          DO WHILE (ttot < ptstep)

             adt=2.e-2
             ! New vapor concentration
             zhlp1 = zcwc + adt * ( SUM(zmtae(nstr:nbins)*zwsatae(nstr:nbins)*zcwsurfae(nstr:nbins))  + &
                                    SUM(zmtcd(1:ncld)*zwsatcd(1:ncld)*zcwsurfcd(1:ncld))              + &
                                    SUM(zmtpd(1:nprc)*zwsatpd(1:nprc)*zcwsurfpd(1:nprc))              + &
                                    SUM(zmtid(1:nice)*zwsatid(1:nice)*zcwsurfid(1:nice))              + &
                                    SUM(zmtsd(1:nsnw)*zwsatsd(1:nsnw)*zcwsurfsd(1:nsnw))              )

             zhlp2 = 1. + adt * ( SUM(zmtae(nstr:nbins)) + SUM(zmtcd(1:ncld)) + SUM(zmtpd(1:nprc)) &
                                   + SUM(zmtid(1:nice)) + SUM(zmtsd(1:nsnw)) )
             zcwint = zhlp1/zhlp2
             zcwint = MIN(zcwint,zcwtot)

             IF ( ANY(paero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 ) THEN
                DO cc = nstr,nbins
                   zcwintae(cc) = zcwcae(cc) + min(max(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                        -0.02*zcwcae(cc)),0.05*zcwcae(cc))
                   zwsatae(cc) = acth2o(paero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                END DO
             END IF
             IF ( ANY(pcloud(ii,jj,:)%numc > nlim) ) THEN
                DO cc = 1,ncld
                   zcwintcd(cc) = zcwccd(cc) + min(max(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd(cc)), &
                        -0.02*zcwccd(cc)),0.05*zcwccd(cc))
                   zwsatcd(cc) = acth2o(pcloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                END DO
             END IF
             IF ( ANY(pprecp(ii,jj,:)%numc > prlim) ) THEN
                DO cc = 1,nprc
                   zcwintpd(cc) = zcwcpd(cc) + min(max(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd(cc)), &
                        -0.02*zcwcpd(cc)),0.05*zcwcpd(cc))
                   zwsatpd(cc) = acth2o(pprecp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                END DO
             END IF
             IF (ANY(pice(ii,jj,:)%numc > prlim) ) THEN
                DO cc = 1,nice
                   zcwintid(cc) = zcwcid(cc) + min(max(adt*zmtid(cc)*(zcwint - zwsatid(cc)*zcwsurfid(cc)), &
                        -0.02*zcwcid(cc)),0.05*zcwcid(cc))
                   zwsatid(cc) = zkelvinid(cc)
                END DO
             END IF
             IF (ANY(psnow(ii,jj,:)%numc > prlim) ) THEN
                DO cc = 1,nsnw
                   zcwintsd(cc) = zcwcsd(cc) + min(max(adt*zmtsd(cc)*(zcwint - zwsatsd(cc)*zcwsurfsd(cc)),&
                        -0.02*zcwcsd(cc)),0.05*zcwcsd(cc))
                   zwsatsd(cc) = zkelvinsd(cc)
                END DO
             END IF
             zcwintae(nstr:nbins) = MAX(zcwintae(nstr:nbins),0.)
             zcwintcd(1:ncld) = MAX(zcwintcd(1:ncld),0.)
             zcwintpd(1:nprc) = MAX(zcwintpd(1:nprc),0.)
             zcwintid(1:nice) = MAX(zcwintid(1:nice),0.)
             zcwintsd(1:nsnw) = MAX(zcwintsd(1:nsnw),0.)

             ! Updae vapor concentration for consistency
             zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                               SUM(zcwintcd(1:ncld))  - &
                               SUM(zcwintpd(1:nprc))  - &
                               SUM(zcwintid(1:nice))  - &
                               SUM(zcwintsd(1:nsnw))

             ! Update "old" values for next cycle
             zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd
             zcwc = zcwint

             ttot = ttot + adt

          END DO ! ADT

          zcwn = zcwint
          zcwnae = zcwintae
          zcwncd = zcwintcd
          zcwnpd = zcwintpd
          zcwnid = zcwintid
          zcwnsd = zcwintsd

          prv(ii,jj) = zcwn*mwa/rhoair

          paero(ii,jj,1:nbins)%volc(8) = max(0.,zcwnae(1:nbins)*mwa/rhowa)
          pcloud(ii,jj,1:ncld)%volc(8) = max(0.,zcwncd(1:ncld)*mwa/rhowa)
          pprecp(ii,jj,1:nprc)%volc(8) = max(0.,zcwnpd(1:nprc)*mwa/rhowa)
          pice(ii,jj,1:nice)%volc(8) = max(0.,zcwnid(1:nice)*mwa/rhowa)
          psnow(ii,jj,1:nsnw)%volc(8) = max(0.,zcwnsd(1:nsnw)*mwa/rhowa)

       END DO !kproma

    END DO ! klev

  END SUBROUTINE gpparth2o
  !-------------------------------------------------------
  REAL FUNCTION acth2o(ppart,pcw)

    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhoss, mss,   &
                               rhowa, mwa,   &
                               rhono, mno,   &
                               rhonh, mnh, eps
    IMPLICIT NONE

    TYPE(t_section), INTENT(in) :: ppart
    REAL, INTENT(in), OPTIONAL :: pcw

    REAL :: zns, znw

    zns =  ( 3.*(ppart%volc(1)*rhosu/msu) +  &
                   (ppart%volc(2)*rhooc/moc) +  &
             2.*(ppart%volc(5)*rhoss/mss) +  &
                   (ppart%volc(6)*rhono/mno) +  &
                   (ppart%volc(7)*rhonh/mnh) )

    IF (PRESENT(pcw)) THEN
       znw = pcw
    ELSE
       znw = ppart%volc(8)*rhowa/mwa
    END IF

    ! Assume activity coefficient of 1 for water...
    acth2o = MAX(0.1,znw/max(eps,(znw+zns)))
  END FUNCTION acth2o

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE gpparthno3(kproma,kbdim,klev,krow,ppres,ptemp,paero,pcloud,    &
                        pprecp,pghno3,pgnh3,prv,prs,pbeta,ptstep)
    
    USE mo_submctl, ONLY : t_section,           &
                               nbins, ncld, nprc,   &
                               surfw0, mvno, mvnh, boltz, &
                               rhono, mno,          &
                               rhonh, mnh,          &
                               rhosu, msu,          &
                               avog, pi,            &
                               pstand,              &
                               nlim, prlim

    IMPLICIT NONE

    INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
    REAL, INTENT(in) :: ptstep
    REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev)
    REAL, INTENT(in) :: prv(kbdim,klev),prs(kbdim,klev)
    REAL, INTENT(in) :: pbeta(kbdim,klev,nbins)

    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,nbins),   &
                                      pcloud(kbdim,klev,ncld),   &
                                      pprecp(kbdim,klev,nprc)
    REAL, INTENT(inout) :: pghno3(kbdim,klev),   &
                               pgnh3(kbdim,klev)

    REAL :: zkelno3ae(nbins), zkelno3cd(ncld), zkelno3pd(nprc)     ! Kelvin effects for HNO3
    REAL :: zkelnh3ae(nbins), zkelnh3cd(ncld), zkelnh3pd(nprc)     ! Kelvin effects for NH3

    REAL :: zcno3cae(nbins), zcno3intae(nbins), zcno3nae(nbins),  & ! Current, intermediate and new HNO3 in aerosols
                zcnh3cae(nbins), zcnh3intae(nbins), zcnh3nae(nbins),  & !  -  ''  - NH3

                zcno3ccd(ncld), zcno3intcd(ncld), zcno3ncd(ncld),  & ! -  ''  - HNO3 in cloud droplets
                zcnh3ccd(ncld), zcnh3intcd(ncld), zcnh3ncd(ncld),  & ! -  ''  - NH3

                zcno3cpd(nprc), zcno3intpd(nprc), zcno3npd(nprc),  & ! -  ''  - HNO3 in precipitation
                zcnh3cpd(nprc), zcnh3intpd(nprc), zcnh3npd(nprc)     ! -  ''  - NH3

    REAL :: zcno3c, zcno3int, zcno3n                             ! Current, intermediate and new HNO3 gas concentration
    REAL :: zcnh3c, zcnh3int, zcnh3n                             ! -  ''  - NH3

    REAL :: zcgnh3eqae(nbins), zcgno3eqae(nbins), &              ! Equilibrium gas concentrations
                zcgnh3eqcd(ncld), zcgno3eqcd(ncld), &
                zcgnh3eqpd(nprc), zcgno3eqpd(nprc)

    REAL :: zacno3ae(nbins), zacno3cd(ncld), zacno3pd(nprc)       ! Activity coefficients for HNO3
    REAL :: zacnh3ae(nbins), zacnh3cd(ncld), zacnh3pd(nprc)       ! Activity coefficients for NH3
    REAL :: zacnh4hso2ae(nbins), zacnh4hso2cd(ncld), zacnh4hso2pd(nprc)
    REAL :: zachhso4ae(nbins), zachhso4cd(ncld), zachhso4pd(nprc)
    REAL :: zmolsae(nbins,7),zmolscd(ncld,7),zmolspd(nprc,7)  ! Ion molalities from pdfite

    REAL :: zcnh3tot, zcno3tot                                                 ! Total mole concentrations

    REAL :: zmtno3ae(nbins), zmtno3cd(ncld), zmtno3pd(nprc) ! Mass transfer coefficients for HNO3
    REAL :: zmtnh3ae(nbins), zmtnh3cd(ncld), zmtnh3pd(nprc) ! Mass transfer coefficients for NH3

    REAL :: zsathno3ae(nbins), zsathno3cd(ncld), zsathno3pd(nprc)
    REAL :: zsatnh3ae(nbins), zsatnh3cd(ncld), zsatnh3pd(nprc)

    REAL :: zbeta
    REAL :: zdfvap ! Diffusion coefficient for vapors

    REAL :: zrh

    REAL :: zhlp1,zhlp2,zhlp3,zhlp4,zhlp5,zhlp6

    REAL :: adt,adtcae(2,nbins),adtccd(2,ncld),adtcpd(2,nprc) ! Adaptive timestep
    REAL :: adtc2(2,nbins)

    INTEGER :: nstr
    INTEGER :: ii,jj,cc

    nstr = 1
    adtcae(:,:) = 0.
    adtccd(:,:) = 0.
    adtcpd(:,:) = 0.
    adtc2(:,:) = 0.

    DO jj = 1,klev
       DO ii = 1,kbdim

          zkelno3pd = 1.; zkelno3cd = 1.; zkelno3ae = 1.
          zkelnh3pd = 1.; zkelnh3cd = 1.; zkelnh3ae = 1.
          zacno3ae = 1.; zacno3cd = 1.; zacno3pd = 1.
          zacnh3ae = 1.; zacnh3cd = 1.; zacnh3pd = 1.
          zsathno3ae = 1.; zsathno3cd = 1.; zsathno3pd = 1.
          zsatnh3ae = 1.; zsatnh3cd = 1.; zsatnh3pd = 1.

          zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]

          ! Kelvin effects
          zkelno3ae(1:nbins) = exp( 4.*surfw0*mvno /  &
                                    (boltz*ptemp(ii,jj)*paero(ii,jj,1:nbins)%dwet) )
          zkelnh3ae(1:nbins) = exp( 4.*surfw0*mvnh /  &
                                    (boltz*ptemp(ii,jj)*paero(ii,jj,1:nbins)%dwet) )

          zkelno3cd(1:ncld) = exp( 4.*surfw0*mvno /  &
                                   (boltz*ptemp(ii,jj)*pcloud(ii,jj,1:ncld)%dwet) )
          zkelnh3cd(1:ncld) = exp( 4.*surfw0*mvnh /  &
                                   (boltz*ptemp(ii,jj)*pcloud(ii,jj,1:ncld)%dwet) )

          zkelno3pd(1:nprc) = exp( 4.*surfw0*mvno /  &
                                   (boltz*ptemp(ii,jj)*pprecp(ii,jj,1:nprc)%dwet) )
          zkelnh3pd(1:nprc) = exp( 4.*surfw0*mvnh /  &
                                   (boltz*ptemp(ii,jj)*pprecp(ii,jj,1:nprc)%dwet) )

          ! Current gas concentrations
          zcno3c = pghno3(ii,jj)/avog
          zcnh3c = pgnh3(ii,jj)/avog

          ! Current particle concentrations
          zcno3cae(1:nbins) = paero(ii,jj,1:nbins)%volc(6)*rhono/mno
          zcnh3cae(1:nbins) = paero(ii,jj,1:nbins)%volc(7)*rhonh/mnh

          zcno3ccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(6)*rhono/mno
          zcnh3ccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(7)*rhonh/mnh

          zcno3cpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(6)*rhono/mno
          zcnh3cpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(7)*rhonh/mnh

          ! Total concentrations
          zcno3tot = zcno3c + SUM(zcno3cae(1:nbins)) +   &
                              SUM(zcno3ccd(1:ncld))       +   &
                              SUM(zcno3cpd(1:nprc))

          zcnh3tot = zcnh3c + SUM(zcnh3cae(1:nbins)) +   &
                              SUM(zcnh3ccd(1:ncld))       +   &
                              SUM(zcnh3cpd(1:nprc))

          ! BETA PUUTTUU PILVILTÄ
          zbeta = 1.

          ! Mass transfer coefficients
          zmtno3ae = 0.; zmtnh3ae = 0.
          zmtno3cd = 0.; zmtnh3cd = 0.
          zmtno3pd = 0.; zmtnh3pd = 0.
          zmtno3ae(1:nbins) = 2.*pi*paero(ii,jj,1:nbins)%dwet *  &
                              zdfvap*paero(ii,jj,1:nbins)%numc*pbeta(ii,jj,1:nbins)
          zmtnh3ae(1:nbins) = 2.*pi*paero(ii,jj,1:nbins)%dwet *  &
                              zdfvap*paero(ii,jj,1:nbins)%numc*pbeta(ii,jj,1:nbins)

          zmtno3cd(1:ncld) = 2.*pi*pcloud(ii,jj,1:ncld)%dwet *  &
                             zdfvap*pcloud(ii,jj,1:ncld)%numc*zbeta
          zmtnh3cd(1:ncld) = 2.*pi*pcloud(ii,jj,1:ncld)%dwet *  &
                             zdfvap*pcloud(ii,jj,1:ncld)%numc*zbeta

          zmtno3pd(1:nprc) = 2.*pi*pprecp(ii,jj,1:nprc)%dwet *  &
                             zdfvap*pprecp(ii,jj,1:nprc)%numc*zbeta
          zmtnh3pd(1:nprc) = 2.*pi*pprecp(ii,jj,1:nprc)%dwet *  &
                             zdfvap*pprecp(ii,jj,1:nprc)%numc*zbeta

          zrh = prv(ii,jj)/prs(ii,jj)

          zmolsae = 0.
          zmolscd = 0.
          zmolspd = 0.

          ! Get the equilibrium concentrations
          ! Aerosols
          CALL NONHEquil(nbins,zrh,ptemp(ii,jj),paero(ii,jj,:),    &
                         zcgno3eqae,zcgnh3eqae,zacno3ae,zacnh3ae,  &
                         zacnh4hso2ae,zachhso4ae,zmolsae           )

          ! Cloud droplets
          CALL NONHEquil(ncld,zrh,ptemp(ii,jj),pcloud(ii,jj,:),    &
                         zcgno3eqcd,zcgnh3eqcd,zacno3cd,zacnh3cd,  &
                         zacnh4hso2cd,zachhso4cd,zmolscd           )

          ! Precipitation
          CALL NONHEquil(nprc,zrh,ptemp(ii,jj),pprecp(ii,jj,:),    &
                         zcgno3eqpd,zcgnh3eqpd,zacno3pd,zacnh3pd,  &
                         zacnh4hso2pd,zachhso4pd,zmolspd           )

          zsatnh3ae = 1.; zsathno3ae = 1.
          zsatnh3cd = 1.; zsathno3cd = 1.
          zsatnh3pd = 1.; zsathno3pd = 1.
          ! NH4/HNO3 saturation ratios for
          ! aerosols
          CALL SVsat(nbins,ptemp(ii,jj),paero(ii,jj,:),zacno3ae,zacnh3ae,zacnh4hso2ae,  &
               zachhso4ae,zcgno3eqae,zcno3cae,zcnh3cae,zkelno3ae,zkelnh3ae,       &
               zsathno3ae,zsatnh3ae,zmolsae,nlim                                          )
          ! clouds
          CALL SVsat(ncld,ptemp(ii,jj),pcloud(ii,jj,:),zacno3cd,zacnh3cd,zacnh4hso2cd,  &
               zachhso4cd,zcgno3eqcd,zcno3ccd,zcnh3ccd,zkelno3cd,zkelnh3cd,       &
               zsathno3cd,zsatnh3cd,zmolscd,nlim                                          )
          ! precipitation
          CALL SVsat(nprc,ptemp(ii,jj),pprecp(ii,jj,:),zacno3pd,zacnh3pd,zacnh4hso2pd,  &
               zachhso4pd,zcgno3eqpd,zcno3cpd,zcnh3cpd,zkelno3pd,zkelnh3pd,       &
               zsathno3pd,zsatnh3pd,zmolspd,prlim                                         )

          adt = ptstep

          zhlp1 = SUM( zcno3cae(nstr:nbins) /  &
               (1. + adt*zmtno3ae(nstr:nbins)*zsathno3ae(nstr:nbins)) )
          zhlp2 = SUM( zcno3ccd(1:ncld)     /  &
               (1. + adt*zmtno3cd(1:ncld)*zsathno3cd(1:ncld)) )
          zhlp3 = SUM( zcno3cpd(1:nprc)     /  &
               (1. + adt*zmtno3pd(1:nprc)*zsathno3pd(1:nprc)) )
          zhlp4 = SUM( zmtno3ae(nstr:nbins)/(1. + adt*zmtno3ae(nstr:nbins)*zsathno3ae(nstr:nbins)) )
          zhlp5 = SUM( zmtno3cd(1:ncld)/(1. + adt*zmtno3cd(1:ncld)*zsathno3cd(1:ncld)) )
          zhlp6 = SUM( zmtno3pd(1:nprc)/(1. + adt*zmtno3pd(1:nprc)*zsathno3pd(1:nprc)) )
          zcno3int = ( zcno3tot - (zhlp1 + zhlp2 + zhlp3) )/( 1. + adt*(zhlp4 + zhlp5 + zhlp6) )

          zhlp1 = SUM( zcnh3cae(nstr:nbins) / &
               (1. + adt*zmtnh3ae(nstr:nbins)*zsatnh3ae(nstr:nbins)) )
          zhlp2 = SUM( zcnh3ccd(1:ncld)     / &
               (1. + adt*zmtnh3cd(1:ncld)*zsatnh3cd(1:ncld)) )
          zhlp3 = SUM( zcnh3cpd(1:nprc)     / &
               (1. + adt*zmtnh3pd(1:nprc)*zsatnh3pd(1:nprc)) )
          zhlp4 = SUM( zmtnh3ae(nstr:nbins)/(1. + adt*zmtnh3ae(nstr:nbins)*zsatnh3ae(nstr:nbins)) )
          zhlp5 = SUM( zmtnh3cd(1:ncld)/(1. + adt*zmtnh3cd(1:ncld)*zsatnh3cd(1:ncld))  )
          zhlp6 = SUM( zmtnh3pd(1:nprc)/(1. + adt*zmtnh3pd(1:nprc)*zsatnh3pd(1:nprc))  )
          zcnh3int = ( zcnh3tot - (zhlp1 + zhlp2 + zhlp3) )/( 1. + adt*(zhlp4 + zhlp5 + zhlp6) )

          zcno3int = MIN(zcno3int, zcno3tot)
          zcnh3int = MIN(zcnh3int, zcnh3tot)

          ! Calculate the new particle concentrations
          zcno3intae(:) = zcno3cae(:)
          zcno3intcd(:) = zcno3ccd(:)
          zcno3intpd(:) = zcno3cpd(:)
          zcnh3intae(:) = zcnh3cae(:)
          zcnh3intcd(:) = zcnh3ccd(:)
          zcnh3intpd(:) = zcnh3cpd(:)
          DO cc = nstr,nbins
             zcno3intae(cc) = ( zcno3cae(cc) + adt*zmtno3ae(cc)*zcno3int ) /  &
                  ( 1. + adt*zmtno3ae(cc)*zsathno3ae(cc) )
             zcnh3intae(cc) = ( zcnh3cae(cc) + adt*zmtnh3ae(cc)*zcnh3int ) /  &
                  ( 1. + adt*zmtnh3ae(cc)*zsatnh3ae(cc) )
          END DO
          DO cc = 1,ncld
             zcno3intcd(cc) = ( zcno3ccd(cc) + adt*zmtno3cd(cc)*zcno3int ) /  &
                  ( 1. + adt*zmtno3cd(cc)*zsathno3cd(cc) )
             zcnh3intcd(cc) = ( zcnh3ccd(cc) + adt*zmtnh3cd(cc)*zcnh3int ) /  &
                  ( 1. + adt*zmtnh3cd(cc)*zsatnh3cd(cc) )
          END DO
          DO cc = 1,nprc
             zcno3intpd(cc) = ( zcno3cpd(cc) + adt*zmtno3pd(cc)*zcno3int ) /  &
                  ( 1. + adt*zmtno3pd(cc)*zsathno3pd(cc) )
             zcnh3intpd(cc) = ( zcnh3cpd(cc) + adt*zmtnh3pd(cc)*zcnh3int ) /  &
                  ( 1. + adt*zmtnh3pd(cc)*zsatnh3pd(cc) )
          END DO
          zcno3intae(1:nbins) = MAX(zcno3intae(1:nbins),0.)
          zcno3intcd(1:ncld) = MAX(zcno3intcd(1:ncld),0.)
          zcno3intpd(1:nprc) = MAX(zcno3intpd(1:nprc),0.)
          zcnh3intae(1:nbins) = MAX(zcnh3intae(1:nbins),0.)
          zcnh3intcd(1:ncld) = MAX(zcnh3intcd(1:ncld),0.)
          zcnh3intpd(1:nprc) = MAX(zcnh3intpd(1:nprc),0.)

          !zcno3int = zcno3tot - SUM(zcno3intae(1:nbins)) - &
          !     SUM(zcno3intcd(1:ncld))  - &
          !     SUM(zcno3intpd(1:nprc))

          !zcnh3int = zcnh3tot - SUM(zcnh3intae(1:nbins)) - &
          !     SUM(zcnh3intcd(1:ncld))  - &
          !     SUM(zcnh3intpd(1:nprc))

          zcno3n = zcno3int
          zcno3nae(:) = zcno3intae(:); zcno3ncd(:) = zcno3intcd(:); zcno3npd(:) = zcno3intpd(:)
          zcnh3n = zcnh3int
          zcnh3nae(:) = zcnh3intae(:); zcnh3ncd(:) = zcnh3intcd(:); zcnh3npd(:) = zcnh3intpd(:)

          ! Model timestep reached - update the new arrays
          pghno3(ii,jj) = zcno3n*avog
          pgnh3(ii,jj) = zcnh3n*avog

          paero(ii,jj,1:nbins)%volc(6) = zcno3nae(1:nbins)*mno/rhono
          pcloud(ii,jj,1:ncld)%volc(6) = zcno3ncd(1:ncld)*mno/rhono
          pprecp(ii,jj,1:nprc)%volc(6) = zcno3npd(1:nprc)*mno/rhono

          paero(ii,jj,1:nbins)%volc(7) = zcnh3nae(1:nbins)*mnh/rhonh
          pcloud(ii,jj,1:ncld)%volc(7) = zcnh3ncd(1:ncld)*mnh/rhonh
          pprecp(ii,jj,1:nprc)%volc(7) = zcnh3npd(1:nprc)*mnh/rhonh

       END DO

    END DO

  END SUBROUTINE gpparthno3

  ! ---------------------------------------------------------------

  REAL FUNCTION acthno3(ppart,pgamma,pchno3p)
    
    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhobc, mbc,   &
                               rhodu, mdu,   &
                               rhoss, mss,   &
                               rhono, mno,   &
                               rhonh, mnh,   &
                               rhowa, mwa

    IMPLICIT NONE

    TYPE(t_section), INTENT(in) :: ppart
    REAL, INTENT(in), OPTIONAL :: pchno3p  ! Current particle HNO3 mole concentration
    REAL, INTENT(in) :: pgamma

    REAL :: zns,znhno3

    ! Solute ~everything else (soluble)?
    zns = ( 3.*(ppart%volc(1)*rhosu/msu)  + &
                  (ppart%volc(2)*rhooc/moc)  + &
            2.*(ppart%volc(5)*rhoss/mss)  + &
                  (ppart%volc(7)*rhonh/mnh)  + &
                  (ppart%volc(8)*rhowa/mwa))

    IF (PRESENT(pchno3p)) THEN
       znhno3 = pchno3p
    ELSE
       znhno3 = ppart%volc(6)*rhono/mno
    END IF

    acthno3 = MAX(znhno3/(zns+znhno3),1.e-3)*pgamma
  END FUNCTION acthno3
  ! -------------------------------------------------------
  REAL FUNCTION actnh3(ppart,pgamma,pcnh3p)
    
    USE mo_submctl, ONLY : t_section,  &
                               rhosu, msu,   &
                               rhooc, moc,   &
                               rhobc, mbc,   &
                               rhodu, mdu,   &
                               rhoss, mss,   &
                               rhono, mno,   &
                               rhonh, mnh,   &
                               rhowa, mwa

    IMPLICIT NONE

    TYPE(t_section), INTENT(in) :: ppart
    REAL, INTENT(in), OPTIONAL :: pcnh3p  ! Current particle HNO3 mole concentration
    REAL, INTENT(in) :: pgamma

    REAL :: zns,znnh3

    ! Solute ~everything else (soluble)?
    zns = ( 3.*(ppart%volc(1)*rhosu/msu)  + &
                  (ppart%volc(2)*rhooc/moc)  + &
            2.*(ppart%volc(5)*rhoss/mss)  + &
                  (ppart%volc(6)*rhono/mno)  + &
                  (ppart%volc(8)*rhowa/mwa))

    IF (PRESENT(pcnh3p)) THEN
       znnh3 = pcnh3p
    ELSE
       znnh3 = ppart%volc(7)*rhonh/mnh
    END IF

    actnh3 = MAX(znnh3/(zns+znnh3),1.e-3)*pgamma
  END FUNCTION actnh3


  ! ------------------------------------------------------------------
  SUBROUTINE NONHEquil(nb,prh,ptemp,ppart,pcgno3eq,pcgnh3eq,      &
                       pgammano,pgammanh,pgammanh4hso2,pgammahhso4,pmols)
    
    USE mo_submctl, ONLY : t_section,    &
                               rhosu,msu,    &
                               rhoss,mss,    &
                               rhono,mno,    &
                               rhonh,mnh,    &
                               rhowa,mwa,    &
                               rg,nlim
    USE aerosol_thermodynamics, ONLY : inorganic_pdfite
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nb
    TYPE(t_section), INTENT(in) :: ppart(nb)
    REAL, INTENT(in) :: prh,ptemp

    REAL, INTENT(out) :: pgammano(nb),pgammanh(nb),pgammanh4hso2(nb),pgammahhso4(nb)
    REAL, INTENT(out) :: pcgno3eq(nb),pcgnh3eq(nb)
    REAL, INTENT(out) :: pmols(nb,7)

    REAL :: zions(7)                ! mol/m3

    REAL :: zwatertotal,        &   ! Total water in particles (mol/m3) ???
                zphno3,zphcl,zpnh3, &   ! Equilibrium vapor pressures (Pa??)
                zgammas(7)              ! Activity coefficients

    REAL :: zhlp

    INTEGER :: cc

    pmols = 0.

    DO cc = 1,nb

       ! 2*H2SO4 + CL + NO3 - Na - NH4

       zhlp = 2.*ppart(cc)%volc(1)*rhosu/msu + ppart(cc)%volc(5)*rhoss/mss + &
              ppart(cc)%volc(6)*rhono/mno - ppart(cc)%volc(5)*rhoss/mss -  &
              ppart(cc)%volc(7)*rhonh/mnh

       zhlp = MAX(zhlp, 1.e-30)

       zions(1) = zhlp !0. ! H+
       zions(2) = ppart(cc)%volc(7)*rhonh/mnh ! NH4
       zions(3) = ppart(cc)%volc(5)*rhoss/mss ! Na
       zions(4) = ppart(cc)%volc(1)*rhosu/msu ! SO4
       zions(5) = 0. ! HSO4
       zions(6) = ppart(cc)%volc(6)*rhono/mno ! NO3
       zions(7) = ppart(cc)%volc(5)*rhoss/mss ! Cl

       zwatertotal = ppart(cc)%volc(8)*rhowa/mwa

       CALL inorganic_pdfite(prh,ptemp,zions,zwatertotal,zphno3,zphcl,zpnh3,zgammas,pmols(cc,:))

       pgammano(cc) = zgammas(1)
       pgammanh(cc) = zgammas(3)
       pgammanh4hso2(cc) = zgammas(6)
       pgammahhso4(cc) = zgammas(7)

       pcgno3eq(cc) = zphno3/(rg*ptemp)
       pcgnh3eq(cc) = zpnh3/(rg*ptemp)

    END DO

  END SUBROUTINE NONHEquil

  ! ------------------------------------------------------------------

  SUBROUTINE SVsat(nb,ptemp,ppart,pachno3,pacnh3,pacnh4hso2,   &
                   pachhso4,pchno3eq,pchno3,pcnh3,pkelhno3,    &
                   pkelnh3,psathno3,psatnh3,pmols,plim         )
    
    USE mo_submctl, ONLY : t_section,   &
                               rhosu,msu,   &
                               rhooc,moc,   &
                               rhoss,mss,   &
                               rhono,mno,   &
                               rhonh,mnh,   &
                               rhowa,mwa,   &
                               rg, nlim
    IMPLICIT NONE

    ! Calculates the saturation ratio for semivolatile species

    INTEGER, INTENT(in) :: nb
    REAL, INTENT(in) :: ptemp
    TYPE(t_section), INTENT(in) :: ppart(nb)
    REAL, INTENT(in) :: pachno3(nb),pacnh3(nb),pacnh4hso2(nb),pachhso4(nb) ! Activity coefficients
    REAL, INTENT(in) :: pchno3eq(nb) ! Equolibrium surface concentration of HNO3
    REAL, INTENT(in) :: pchno3(nb)   ! Current particle concentration of HNO3
    REAL, INTENT(in) :: pcnh3(nb)    ! Current particle concentration of NH3
    REAL, INTENT(in) :: pkelhno3(nb), pkelnh3(nb)  ! Kelvin effects
    REAL, INTENT(in) :: pmols(nb,7)
    REAL, INTENT(in) :: plim
    REAL, INTENT(out) :: psathno3(nb), psatnh3(nb)

    REAL :: KllH2O, KllNH3, KglNH3, KglHNO3

    REAL, PARAMETER :: zt0 = 298.15    ! Reference temp
    REAL, PARAMETER :: zatm = 101325.  ! Unit atm in Pa

    REAL, PARAMETER :: K01 = 1.01e-14,   &
                           K02 = 1.81e-5,    &
                           K03 = 57.64,      &
                           K04 = 2.51e6

    REAL, PARAMETER :: a1 = -22.52,      &
                           a2 = 1.50,        &
                           a3 = 13.79,       &
                           a4 = 29.47,       &
                           b1 = 26.92,       &
                           b2 = 26.92,       &
                           b3 = -5.39,       &
                           b4 = 16.84

    REAL :: zmolno3     ! molality of NO3-
    REAL :: zmolhp      ! molality of H+
    REAL :: zmolso4,  & ! molality of SO4
                zmolcl,   & ! molality of Cl
                zmolnh4,  & ! Molality of NH4
                zmolna      ! Molality of Na
    REAL :: zhlp1,zhlp2,zhlp3, zxi

    INTEGER :: cc

    zhlp1 = 0.; zhlp2 = 0.; zhlp3 = 0.

    zmolno3 = 0.
    zmolhp = 0.

    ! Calculates K^ll_h20, K^ll_NH3, K^gl_NH3, K^gl_HNO3
    zhlp1 = zt0/ptemp
    zhlp2 = zhlp1-1.
    zhlp3 = 1. + LOG(zhlp1) - zhlp1

    KllH2O = K01*EXP( a1*zhlp2 + b1*zhlp3 )
    KllNH3 = K02*EXP( a2*zhlp2 + b2*zhlp3 )
    KglNH3 = k03*EXP( a3*zhlp2 + b3*zhlp3 )
    KglHNO3 = k04*EXP( a4*zhlp2 + b4*zhlp3 )

    ! Get NO3- and H+ molality
    DO cc = 1,nb

       IF ( ppart(cc)%numc > 1.e-40) THEN

          zhlp1 = pcnh3(cc)*mnh + ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(2)*rhooc +  &
               ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
          zmolno3 = pchno3(cc) / zhlp1

          zxi = ( pcnh3(cc) + ppart(cc)%volc(5)*rhoss/mss ) / &
                ( ppart(cc)%volc(1)*rhosu/msu )
          IF ( zxi <= 2. ) THEN
             ! Molality of SO4
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
             zmolso4 = (ppart(cc)%volc(1)*rhosu/msu)/zhlp1
             ! Molality of Cl
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(8)*rhowa
             zmolcl = (ppart(cc)%volc(5)*rhoss/mss)/zhlp1
             ! Molality of NH4
             zhlp1 =  pchno3(cc)*mno + ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(5)*rhoss + ppart(cc)%volc(8)*rhowa
             zmolnh4 = pcnh3(cc)/zhlp1
             ! Molality of Na
             zhlp1 = pcnh3(cc)*mnh + pchno3(cc)*mno + ppart(cc)%volc(2)*rhooc + &
                  ppart(cc)%volc(1)*rhosu + ppart(cc)%volc(8)*rhowa
             zmolna = (ppart(cc)%volc(5)*rhoss/mss)/zhlp1

             zmolhp = 2.*zmolso4 + zmolno3 + zmolcl - (zmolnh4 + zmolna)

          ELSE

             zhlp2 = pkelhno3(cc)*zmolno3*pachno3(cc)**2.
             zmolhp = KglHNO3*pchno3eq(cc)/zhlp2

          END IF

          zhlp1 = ppart(cc)%volc(8)*rhowa*rg*ptemp*KglHNO3
          ! Saturation ratio for NH3
          zhlp2 = pkelnh3(cc)/(zhlp1*zmolhp)
          zhlp3 = KllH2O/(KllNH3+KglNH3)
          psatnh3(cc) = zhlp2 * ( (pacnh4hso2(cc)/pachhso4(cc))**2. ) * zhlp3
          ! Saturation ratio for HNO3
          psathno3(cc) = ( pkelhno3(cc)*zmolhp*pachno3(cc)**2 ) / zhlp1

       ELSE

          psatnh3(cc) = 0.
          psathno3(cc) = 0.
       END IF

    END DO

  END SUBROUTINE SVsat


  ! ------------------------------------------------------------------

  FUNCTION GetTstep(nb,zcg,zcs,zmt,zconst) RESULT(tscale)
    
    IMPLICIT NONE

    INTEGER, INTENT(in) :: nb
    REAL, INTENT(in) :: zcg
    REAL, INTENT(in) :: zcs(nb), zmt(nb)
    REAL, INTENT(in) :: zconst

    REAL :: th(nb)
    REAL :: tscale

    INTEGER :: cc

    DO cc = 1,nb
       th(cc) = (zcg - zcs(cc))/MAX(zcg,zcs(cc))
    END DO

    tscale = zconst/SUM( ABS(zmt(:)*th(:)),MASK=(th(:)*zmt(:) /= 0.) )

  END FUNCTION GetTstep

!
! ----------------------------------------------------------------------------------------------------------
!

  FUNCTION satvaph2o(ptemp) RESULT(psat)
    !-----------------------------------------------------------------
    ! Saturation vapour pressure of water vapour
    ! This is a local function for the subroutine *cloud_condensation*
    !
    ! J. Tonttila, FMI, 03/2014
    !-----------------------------------------------------------------
    
    IMPLICIT NONE

    REAL, INTENT(in) :: ptemp

    REAL, PARAMETER ::        &
         za0 = 6.107799961,    &
         za1 = 4.436518521e-1, &
         za2 = 1.428945805e-2, &
         za3 = 2.650648471e-4, &
         za4 = 3.031240396e-6, &
         za5 = 2.034080948e-8, &
         za6 = 6.136820929e-11

    REAL :: zt

    REAL :: psat

    zt = ptemp - 273.15

    psat = za0 + za1*zt + za2*zt**2 + za3*zt**3 +   &
           za4*zt**4 + za5*zt**5 + za6*zt**6

    ! To Pascals
    psat = psat*100.

  END FUNCTION satvaph2o
!
! --------------------------------------------------------------------------------------------------------------
!

  !------------------------------------------------
  !
  ! ***************
  ! FUNCTION coagc
  ! ***************
  !
  ! Calculation of coagulation coefficients.
  ! Extended version of the function originally
  ! found in mo_salsa_init. This is now placed
  ! here to avoid cyclic dependencies between
  ! modules upon coupling with UCLALES.
  !
  ! J. Tonttila, FMI, 05/2014
  !
  !-------------------------------------------------
  REAL FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres,kernel)

    USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav
    USE mo_constants, ONLY : rd, amd

    IMPLICIT NONE

    !-- Input variables ----------
    REAL, INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres        ! ambient pressure [fxm]

    INTEGER, INTENT(in) :: kernel ! Select the type of kernel: 1 - aerosol-aerosol coagulation (the original version)
                                  !                            2 - hydrometeor-aerosol or hydrometeor-hydrometeor coagulation

    !-- Output variables ---------

    !-- Local variables ----------
    REAL ::  &
         visc,   &   ! viscosity of air [kg/(m s)]
         vkin,   &   ! Kinematic viscosity of air [m2 s-1]
         zrhoa,  &   ! Density of air [kg m-3]
         mfp,    &   ! mean free path of air molecules [m]
         mdiam,  &   ! mean diameter of colliding particles [m]
         fmdist, &   ! distance of flux matching [m]
         zecoll, &   ! Collition efficiency for graviational collection
         zev,    &   !
         zea,    &
         zbrown, &   ! Components for coagulation kernel; Brownian
         zbrconv,&   !                                    Convective diffusion enhancement
         zgrav       !                                    Gravitational collection


    REAL, DIMENSION (2) :: &
         diam,   &   ! diameters of particles [m]
         mpart,  &   ! masses of particles [kg]
         knud,   &   ! particle knudsen number [1]
         beta,   &   ! Cunningham correction factor [1]
         zrhop,  &   ! Particle density [kg m-3]
         dfpart, &   ! particle diffusion coefficient [m2/s]
         mtvel,  &   ! particle mean thermal velocity [m/s]
         termv,  &   ! Particle terminal velocity
         omega,  &   !
         tva,    &   ! temporary variable [m]
         flux        ! flux in continuum and free molec. regime [m/s]


    REAL ::  &
         schm(2), &   ! Schmidt nubmer
         reyn(2), &    ! Reynolds number
         stok             ! Stokes number
    REAL :: rhoref               ! Reference air density in STP conditions
    INTEGER :: lrg,sml,i

    zbrown = 0.
    zbrconv = 0.
    zgrav = 0.
    zev = 0.
    coagc = 0.

    rhoref = 1.01325e5/(287.*273.15)

    !-------------------------------------------------------------------------------

    !-- 0) Initializing particle and ambient air variables --------------------
    diam = (/ diam1, diam2 /)       ! particle diameters [m]
    mpart = (/ mass1, mass2 /)       ! particle masses [kg]

    visc = (7.44523e-3*temp**1.5)/(5093.*(temp+110.4)) ! viscosity of air [kg/(m s)]

    mfp = (1.656e-10*temp+1.828e-8)*pstand/pres ! mean free path of air [m]


    !-- 2) Slip correction factor for small particles -------------------------

    knud = 2.*mfp/diam                                    ! Knudsen number
    beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))! Cunningham correction factor
    ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

    !-- 3) Particle properties ------------------------------------------------

    dfpart = beta*boltz*temp/(3.*pi*visc*diam)  ! diffusion coefficient [m2/s]
    mtvel = sqrt((8.*boltz*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
    omega = 8.*dfpart/(pi*mtvel)

    mdiam = 0.5*(diam(1)+diam(2))               ! mean diameter [m]

    !-- 4) Calculation of fluxes and flux matching ----------------------------

    flux(1) = 4.*pi*mdiam*( dfpart(1)+dfpart(2) )    ! flux in continuum regime [m3/s]
    flux(2) = pi*sqrt((mtvel(1)**2)+(mtvel(2)**2))*(mdiam**2) !  -"- in free molec. regime [m3/s]

    tva(1) = ((mdiam+omega(1))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(1)**2)* &
         sqrt((mdiam**2+omega(1)**2)))/ &
         (3.*mdiam*omega(1)) - mdiam

    tva(2) = ((mdiam+omega(2))**3 - &              ! temporary variable [m]
         (mdiam**2+omega(2)**2)* &
         sqrt((mdiam**2+omega(2)**2)))/ &
         (3.*mdiam*omega(2)) - mdiam

    fmdist = sqrt(tva(1)**2+tva(2)**2)             ! flux matching distance [m]

    SELECT CASE(kernel)
       CASE(1)

          ! Aerosol-Aerosol coagulation - like the f version
          !-- 5) Coagulation coefficient [m3/s] -------------------------------------
          coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))

       CASE(2)

          ! Which particle is larger?
          sml = 1; lrg = 2
          IF (diam(1) >= diam(2)) THEN
             lrg = 1; sml = 2
          END IF

          zrhoa = pres/(rd*temp)   ! Density of air
          zrhop = mpart/(pi6*diam**3)             ! Density of particles
          vkin = visc/zrhoa   ! Kinematic viscosity of air [m2 s-1]
          DO i = 1,2
             IF (diam(i) < 40.e-6) THEN
                termv(i) = ( (diam(i)**2) * (zrhop(i) - zrhoa) * grav * beta(i) )/( 18.*visc  ) ![m s-1]
             ELSE IF (diam(i) >= 40.e-6) THEN
                termv(i) = 2.e3*diam(i)*(rhoref/zrhoa)**2
             END IF
          END DO
          ! Reynolds number
          reyn = diam*termv/vkin
          ! Schmidt number for the smaller particle
          schm = vkin/dfpart
          ! Stokes number
          stok = 2.*termv(sml)*ABS(termv(1) - termv(2))/( diam(lrg)*grav )

          !Brownian component
          zbrown = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))

          ! Convective enhancement
          IF (reyn(lrg) <= 1.) THEN
             zbrconv = 0.45*zbrown*( reyn(lrg)**(1./3.) )*( schm(sml)**(1./3.) )
          ELSE IF (reyn(lrg) > 1.) THEN
             zbrconv = 0.45*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1./3.) )
          END IF

          ! gravitational collection
          zea = stok**2/( stok + 0.5 )**2
          IF (stok > 1.214) THEN
             zev = 0.75*LOG(2.*stok)/(stok - 1.214)
             zev = (1. + zev)**(-2.)
          ELSE IF (stok <= 1.214) THEN
             zev = 0.
          END IF

          zecoll = (60.*zev + zea*reyn(lrg))/(60. + reyn(lrg))
          zgrav = zecoll * pi * mdiam**2
          zgrav = zgrav * ABS(termv(1)-termv(2))

          ! Total coagulation kernel
          coagc = zbrown  + zbrconv + zgrav

    END SELECT

  END FUNCTION coagc



END MODULE mo_salsa_dynamics
