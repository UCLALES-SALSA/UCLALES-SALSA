
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
  ! subroutine COAGULATION(kbdim,klev, &
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


  SUBROUTINE coagulation(kbdim,  klev,    &
                         paero,  pcloud, pprecp, pice, psnow,  &
                         ptstep, ptemp,  ppres, pedr     )

    USE mo_submctl, ONLY:        &
         t_section,   & ! Datatypes for the cloud bin representation
         in1a, fn1a,                 & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 &
         ncld, nprc, nice, nsnw,     &
         inp2a, fnp2a, inp2b, fnp2b, &
         pi6,                        &
         nlim,prlim,                 &
         lscgaa, lscgcc, lscgca,     &
         lscgpp, lscgpa, lscgpc,     &
         lscgia, lscgic, lscgii, lscgip, &
         lscgsa, lscgsc, lscgsi, lscgsp, lscgss, &
         nspec, CalcDimension, CalcMass, lscgrain, &
         nlsip_hm, rime_volc_ice, rime_volc_snw, &
         hm_dmin_drop, hm_dmin_ice, &
         nlsip_iibr, coll_rate_ii, coll_rate_si, coll_rate_ss

    IMPLICIT NONE


    !-- Input and output variables -------------
    INTEGER, INTENT(IN) ::          &
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
         ppres(kbdim,klev),         &
         pedr(kbdim,klev)          ! eddy dissipation rate from LES
    !-- Local variables ------------------------
    INTEGER ::                      &
         ii,jj,kk,ll,mm,nn,cc,      & ! loop indices
         index_2a, index_2b,     & ! corresponding bin in regime 2a/2b
         nt

    REAL ::                     &
         zcc(fn2b,fn2b),            & ! updated coagulation coefficients [m3/s]
         zcccc(ncld,ncld),          & ! - '' - for collision-coalescence [m3/s]
         zccca(fn2b,ncld),          & ! - '' - for cloud collection of aerosols [m3/s]
         zccpc(ncld,nprc),          & ! - '' - for collection of cloud droplets by precip [m3/s]
         zccpa(fn2b,nprc),          & ! - '' - for collection of aerosols by precip
         zccpp(nprc,nprc),          & ! - '' - for collisions between precip particles
         zccia(fn2b,nice),          & ! - '' - for collection of aerosols by ice
         zccic(ncld,nice),          & ! - '' - for collection of cloud droplets by ice
         zccii(nice,nice),          & ! - '' - for collisions between ice
         zccip(nprc,nice),          & ! - '' - for collection of precip by ice
         zccsa(fn2b,nsnw),          & ! - '' - for collection of aerosols by snow
         zccsc(ncld,nsnw),          & ! - '' - for collection of cloud droples by snow
         zccsi(nice,nsnw),          & ! - '' - for collection of ice by snow
         zccsp(nprc,nsnw),          & ! - '' - for collection of precip by snow
         zccss(nsnw,nsnw),          & ! - '' - for collisions between snow
         zminusterm,                & ! coagulation loss in a bin [1/s]
         zplusterm(nspec+1)           ! coagulation gain in a bin [fxm/s] (for each chemical compound)

    REAL :: &
         zmpart(fn2b),   & ! approximate mass of aerosol [kg]
         zmcloud(ncld),  & ! approximate mass of cloud droplets [kg]
         zmprecp(nprc),  & ! approximate mass for rain drops [kg]
         zmice(nice),    & ! approximate mass for ice [kg]
         zmsnow(nsnw),   & ! approximate mass for snow [kg]
         zdpart(fn2b),   & ! diameter of aerosol [m]
         zdcloud(ncld),  & ! diameter of cloud droplets [m]
         zdprecp(nprc),  & ! diameter for rain drops [m]
         zdice(nice),    & ! diameter for ice [m]
         zdsnow(nsnw)      ! diameter for snow [m]

    REAL :: temppi, pressi, edri
    REAL :: num_coag(nprc), vol_cltd(nprc,nspec+1), V_new, zminusterm_c(ncld)

    LOGICAL :: any_cloud, any_precp, any_ice, any_snow

    !-----------------------------------------------------------------------------

    IF (.NOT.ALLOCATED(rime_volc_ice)) ALLOCATE(rime_volc_ice(kbdim,klev,nice), &
        rime_volc_snw(kbdim,klev,nsnw))
    rime_volc_ice(:,:,:) = 0.
    rime_volc_snw(:,:,:) = 0.
    IF (.NOT.ALLOCATED(coll_rate_ii)) ALLOCATE(coll_rate_ii(kbdim,klev,nice,nice), &
        coll_rate_si(kbdim,klev,nsnw,nice), coll_rate_ss(kbdim,klev,nsnw,nsnw))
    coll_rate_ii(:,:,:,:) = 0.
    coll_rate_si(:,:,:,:) = 0.
    coll_rate_ss(:,:,:,:) = 0.

    nt = nspec + 1 ! Total number of spcecies + water

     DO jj = 1,klev      ! vertical grid
        DO ii = 1,kbdim ! dimension for arrays

           !-- 1) Updating coagulation coefficients -------------------------------------

           ! Which species are included
           any_cloud = ANY(pcloud(ii,jj,:)%numc > nlim)
           any_precp = ANY(pprecp(ii,jj,:)%numc > prlim)
           any_ice = ANY(pice(ii,jj,:)%numc > prlim)
           any_snow = ANY(psnow(ii,jj,:)%numc > prlim)

           !-- Aerosol diameter [m] and mass [kg]; density of 1500 kg/m3 assumed
           CALL CalcDimension(fn2b,paero(ii,jj,1:fn2b),nlim,1)
           zdpart(1:fn2b) = paero(ii,jj,1:fn2b)%dwet
           CALL CalcMass(zmpart,fn2b,paero(ii,jj,:),nlim,1)

           !-- Cloud droplet diameter and mass; Assume water density
           CALL CalcDimension(ncld,pcloud(ii,jj,1:ncld),nlim,2)
           zdcloud(1:ncld) = pcloud(ii,jj,1:ncld)%dwet
           CALL CalcMass(zmcloud,ncld,pcloud(ii,jj,:),nlim,2)

           !-- Precipitation droplet diameter and mass
           CALL CalcDimension(nprc,pprecp(ii,jj,1:nprc),prlim,3)
           zdprecp(1:nprc) = pprecp(ii,jj,1:nprc)%dwet
           CALL CalcMass(zmprecp,nprc,pprecp(ii,jj,:),prlim,3)

           !-- Ice particle diameter and mass - may not be spherical
           CALL CalcDimension(nice,pice(ii,jj,1:nice),prlim,4)
           zdice(1:nice) = pice(ii,jj,1:nice)%dwet
           CALL CalcMass(zmice,nice,pice(ii,jj,:),prlim,4)

           !-- Snow diameter and mass - may not be spherical
           CALL CalcDimension(nsnw,psnow(ii,jj,1:nsnw),prlim,5)
           zdsnow(1:nsnw) = psnow(ii,jj,1:nsnw)%dwet
           CALL CalcMass(zmsnow,nsnw,psnow(ii,jj,:),prlim,5)

           temppi=ptemp(ii,jj)
           pressi=ppres(ii,jj)
           edri = pedr(ii,jj)
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
                    zcc(mm,nn) = coagc(zdpart(mm),zdpart(nn),zmpart(mm),zmpart(nn),temppi,pressi,edri,1,1)
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
                    zcccc(mm,nn) = coagc(zdcloud(mm),zdcloud(nn),zmcloud(mm),zmcloud(nn),temppi,pressi,edri,2,2)
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
                    zccpp(mm,nn) =  coagc(zdprecp(mm),zdprecp(nn),zmprecp(mm),zmprecp(nn),temppi,pressi,edri,3,3)
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
                    zccca(mm,nn) = coagc(zdpart(mm),zdcloud(nn),zmpart(mm),zmcloud(nn),temppi,pressi,edri,1,2)
                 END DO
              END DO
           END IF
           ! Collection of aerosols by rain
           IF (lscgpa .AND. any_precp) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nprc
                    IF (pprecp(ii,jj,nn)%numc<prlim) cycle
                    zccpa(mm,nn) = coagc(zdpart(mm),zdprecp(nn),zmpart(mm),zmprecp(nn),temppi,pressi,edri,1,3)
                 END DO
              END DO
           END IF
           ! Collection of cloud droplets by rain
           IF (lscgpc .AND. any_cloud .AND. any_precp) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nprc
                    IF (pprecp(ii,jj,nn)%numc<prlim) cycle
                    zccpc(mm,nn) = coagc(zdcloud(mm),zdprecp(nn),zmcloud(mm),zmprecp(nn),temppi,pressi,edri,2,3)
                  END DO
              END DO
           END IF
           !  collection of aerosols by ice
           IF (lscgia .AND. any_ice) THEN
              DO mm = 1,fn2b
                 IF (paero(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) cycle
                    zccia(mm,nn) =  coagc(zdpart(mm),zdice(nn),zmpart(mm),zmice(nn),temppi,pressi,edri,1,4)
                 END DO
              END DO
           END IF
          !  collection of cloud droplets by ice
           IF (lscgic .AND. any_ice .AND. any_cloud) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) cycle
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) cycle
                    zccic(mm,nn) = coagc(zdcloud(mm),zdice(nn),zmcloud(mm),zmice(nn),temppi,pressi,edri,2,4)
                 END DO
              END DO
           END IF
           !  collisions between ice particles
           IF (lscgii .AND. any_ice) THEN
              DO mm = 1,nice
                 IF (pice(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = mm,nice
                    IF (pice(ii,jj,nn)%numc<prlim) CYCLE
                    zccii(mm,nn) = coagc(zdice(mm),zdice(nn),zmice(mm),zmice(nn),temppi,pressi,edri,4,4)
                    zccii(nn,mm) = zccii(mm,nn)
                 END DO
              END DO
           END IF
           !  collection of precip by ice-collision
           IF (lscgip .AND. any_precp .AND. any_ice) THEN
              DO mm = 1,nprc
                 IF (pprecp(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nice
                    IF (pice(ii,jj,nn)%numc<prlim) CYCLE
                    zccip(mm,nn) = coagc(zdprecp(mm),zdice(nn),zmprecp(mm),zmice(nn),temppi,pressi,edri,3,4)
                  END DO
              END DO
           END IF
           ! Self-collection of snow
           IF (lscgss .AND. any_snow) THEN
              DO mm = 1,nsnw
                 IF (psnow(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = mm,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccss(mm,nn) =  coagc(zdsnow(mm),zdsnow(nn),zmsnow(mm),zmsnow(nn),temppi,pressi,edri,5,5)
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
                    zccsa(mm,nn) = coagc(zdpart(mm),zdsnow(nn),zmpart(mm),zmsnow(nn),temppi,pressi,edri,1,5)
                 END DO
              END DO
           END IF
           ! collection of precip by snow
           IF (lscgsp .AND. any_precp .AND. any_snow) THEN
              DO mm = 1,nprc
                 IF (pprecp(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsp(mm,nn) = coagc(zdprecp(mm),zdsnow(nn),zmprecp(mm),zmsnow(nn),temppi,pressi,edri,3,5)
                  END DO
              END DO
           END IF
           ! collection of cloud droples by snow
           IF (lscgsc .AND. any_cloud .AND. any_snow) THEN
              DO mm = 1,ncld
                 IF (pcloud(ii,jj,mm)%numc<nlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsc(mm,nn) = coagc(zdcloud(mm),zdsnow(nn),zmcloud(mm),zmsnow(nn),temppi,pressi,edri,2,5)
                  END DO
              END DO
           END IF
           ! collection of ice by snow
           IF (lscgsi .AND. any_ice .AND. any_snow) THEN
              DO mm = 1,nice
                 IF (pice(ii,jj,mm)%numc<prlim) CYCLE
                 DO nn = 1,nsnw
                    IF (psnow(ii,jj,nn)%numc<prlim) CYCLE
                    zccsi(mm,nn) = coagc(zdice(mm),zdsnow(nn),zmice(mm),zmsnow(nn),temppi,pressi,edri,4,5)
                 END DO
              END DO
           END IF

           !-- 2) New particle and volume concentrations after coagulation --------------

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

              ! Particle volume gained from smaller particles in regime 1a
              DO ll = in1a,kk-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:nt) = ( paero(ii,jj,kk)%volc(1:nt)+ptstep*zplusterm(1:nt) * &
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

              ! Particle volume gained from smaller particles in regimes 1a and 2a
              DO ll = in1a, kk-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              ! Particle volume gained from smaller (and equal) particles in 2b
              DO ll = in2b, index_2b
                 zplusterm(1:nt) = zplusterm(1:nt) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:nt) = ( paero(ii,jj,kk)%volc(1:nt)+ptstep*zplusterm(1:nt) *  &
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

              ! Particle volume gained from smaller particles in 1a and 2a
              DO ll = in1a, index_2a-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              ! Particle volume gained from smaller particles in 2b
              DO ll = in2b, kk-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcc(ll,kk)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              !-- Volume and number concentrations after coagulation update [fxm]
              paero(ii,jj,kk)%volc(1:nt) = ( paero(ii,jj,kk)%volc(1:nt)+ptstep*zplusterm(1:nt) *  &
                   paero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

              paero(ii,jj,kk)%numc = paero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                   0.5*ptstep*zcc(kk,kk)*paero(ii,jj,kk)%numc)

           END DO

           ! Cloud droplets, regime a
           ! ------------------------------------------------
           DO cc = inp2a,fnp2a
              IF (pcloud(ii,jj,cc)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime b cloud droplets
              kk = fnp2a + cc

              ! Droplets lost by those with larger nucleus in regime a
              DO ll = cc+1,fnp2a
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by those with larger nucleus in regime b
              DO ll = kk+1,fnp2b
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
                 zplusterm(1:nt) = zplusterm(1:nt) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller droplets in a
              DO ll = inp2a,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller or equal droplets in b
              DO ll = inp2b,kk
                 zplusterm(1:nt) = zplusterm(1:nt) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
              END DO

              ! Update the hydrometeor volume concentrations
              pcloud(ii,jj,cc)%volc(1:nt) = max(0.,( pcloud(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*pcloud(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pcloud(ii,jj,cc)%numc = max(0., pcloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zcccc(cc,cc)*pcloud(ii,jj,cc)%numc ) )

           END DO

           ! Cloud droplets, regime b
           ! -----------------------------------------
           DO cc = inp2b,fnp2b
              IF (pcloud(ii,jj,cc)%numc<nlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime a cloud droplets
              kk = cc - fnp2a

              ! Droplets lost by those with larger nucleus in regime b
              DO ll = cc+1,fnp2b
                 zminusterm = zminusterm + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
              END DO

              ! Droplets lost by those with larger nucleus in regime a
              DO ll = kk+1,fnp2a
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
                 zplusterm(1:nt) = zplusterm(1:nt) + zccca(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller droplets in b
              DO ll = inp2b,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller or equal droplets in a
              DO ll = inp2a,kk
                 zplusterm(1:nt) = zplusterm(1:nt) + zcccc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
              END DO

              ! Update the hydrometeor volume concentrations
              pcloud(ii,jj,cc)%volc(1:nt) = max(0., ( pcloud(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*pcloud(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
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

              ! Drops lost by collection by snow (assume freezing of rain drops at any T)
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsp(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Drops lost by collisions with ice (assume freezing of rain drops at any T)
              DO ll = 1,nice
                 zminusterm = zminusterm + zccip(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Volume gained by collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:nt) = zplusterm(1:nt) + zccpa(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained by collection of cloud droplets
              DO ll = 1,ncld
                 zplusterm(1:nt) = zplusterm(1:nt) + zccpc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller drops
              DO ll = 1,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zccpp(ll,cc)*pprecp(ii,jj,ll)%volc(1:nt)
              END DO

              ! Update the hydrometeor volume concentrations
              pprecp(ii,jj,cc)%volc(1:nt) = max(0., ( pprecp(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*pprecp(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pprecp(ii,jj,cc)%numc = max(0.,pprecp(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccpp(cc,cc)*pprecp(ii,jj,cc)%numc ) )

           END DO

           ! Rain formation from cloud-cloud collisions
           ! ------------------------------------------------
           IF (lscgrain .AND. any_cloud) THEN
              zminusterm_c(:)=0.
              vol_cltd(:,:)=0.
              num_coag(:)=0.
              DO cc = 1, ncld
                 IF (pcloud(ii,jj,cc)%numc<nlim) CYCLE
                 DO ll = cc, ncld
                    IF (pcloud(ii,jj,ll)%numc<nlim) CYCLE
                    !
                    V_new=pi6*(pcloud(ii,jj,cc)%dwet**3+pcloud(ii,jj,ll)%dwet**3)
                    IF (V_new>pprecp(ii,jj,1)%vlolim) THEN
                       ! Rain bin
                       mm=COUNT(V_new>pprecp(ii,jj,:)%vlolim)
                       ! New rain mass and number
                       vol_cltd(mm,1:nt) = vol_cltd(mm,1:nt) + (pcloud(ii,jj,ll)%volc(1:nt)*pcloud(ii,jj,cc)%numc + &
                                    pcloud(ii,jj,cc)%volc(1:nt)*pcloud(ii,jj,ll)%numc)*zcccc(ll,cc)
                       num_coag(mm) = num_coag(mm) + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc*pcloud(ii,jj,cc)%numc
                       ! Remove the mass from the original cloud bins
                       zminusterm_c(cc) = zminusterm_c(cc) + zcccc(cc,ll)*pcloud(ii,jj,ll)%numc
                       zminusterm_c(ll) = zminusterm_c(ll) + zcccc(cc,ll)*pcloud(ii,jj,cc)%numc
                    ENDIF
                 ENDDO
                 ! Update cloud bin cc
                 pcloud(ii,jj,cc)%numc = pcloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm_c(cc))
                 pcloud(ii,jj,cc)%volc(1:nt) = pcloud(ii,jj,cc)%volc(1:nt)/(1. + ptstep*zminusterm_c(cc))
              ENDDO
              ! Update rain
              DO ll = 1,nprc
                 pprecp(ii,jj,ll)%volc(1:nt) = pprecp(ii,jj,ll)%volc(1:nt) + ptstep*vol_cltd(ll,1:nt)
                 pprecp(ii,jj,ll)%numc = pprecp(ii,jj,ll)%numc + ptstep*num_coag(ll)
              END DO
           ENDIF


           ! Ice particles, regime a
           ! ------------------------------------------------
           DO cc = inp2a,fnp2a
              IF (pice(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime b ice
              kk = fnp2a + cc

              ! Particles lost by those with larger nucleus in regime a
              DO ll = cc+1,fnp2a
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by those with larger nucleus in regime b
              DO ll = kk+1,fnp2b
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by snow
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsi(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained from aerosol collection
              DO ll = in1a,fn2b
                 zplusterm(1:nt) = zplusterm(1:nt) + zccia(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdpart(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccia(ll,cc)*paero(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from cloud collection
              DO ll = 1,ncld
                 zplusterm(1:nt) = zplusterm(1:nt) + zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdcloud(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from rain drops
              DO ll = 1,nprc
                 zplusterm(1:nt) = zplusterm(1:nt) + zccip(ll,cc)*pprecp(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdprecp(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccip(ll,cc)*pprecp(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from smaller ice particles in regime a
              DO ll = inp2a,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller or equal ice particles in regime b
              DO ll = inp2b,kk
                 zplusterm(1:nt) = zplusterm(1:nt) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:nt)
              END DO

              ! Save ice-ice collisions for collisional breakup
              IF (nlsip_iibr) THEN
                 ! Smaller and equal (self) from regime a
                 DO ll = inp2a,cc
                    coll_rate_ii(ii,jj,cc,ll) = &
                       ptstep*zccii(ll,cc)*pice(ii,jj,ll)%numc*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
                 ! Smaller and equal from regime b
                 DO ll = inp2b,kk
                    coll_rate_ii(ii,jj,cc,ll) = &
                       ptstep*zccii(ll,cc)*pice(ii,jj,ll)%numc*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
              ENDIF

              ! Update the hydrometeor volume concentrations
              pice(ii,jj,cc)%volc(1:nt) = max(0., ( pice(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*pice(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pice(ii,jj,cc)%numc = max(0.,pice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccii(cc,cc)*pice(ii,jj,cc)%numc ) )

           END DO

           ! Ice particles, regime b
           ! -----------------------------------------
           DO cc = inp2b,fnp2b
              IF (pice(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! corresponding index for regime a
              kk = cc - fnp2a

              ! Particles lost by those with larger nucleus in regime b
              DO ll = cc+1,fnp2b
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by those with larger nucleus in regime a
              DO ll = kk+1,fnp2a
                 zminusterm = zminusterm + zccii(cc,ll)*pice(ii,jj,ll)%numc
              END DO

              ! Particles lost by collection by snow
              DO ll = 1,nsnw
                 zminusterm = zminusterm + zccsi(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained from aerosol collection
              DO ll = in1a,fn2b
                 zplusterm(1:nt) = zplusterm(1:nt) + zccia(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdpart(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccia(ll,cc)*paero(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from cloud collection
              DO ll = 1,ncld
                 zplusterm(1:nt) = zplusterm(1:nt) + zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdcloud(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccic(ll,cc)*pcloud(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from rain drops
              DO ll = 1,nprc
                 zplusterm(1:nt) = zplusterm(1:nt) + zccip(ll,cc)*pprecp(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdprecp(ll)>hm_dmin_drop .AND. zdice(cc)>hm_dmin_ice) THEN
                    rime_volc_ice(ii,jj,cc) = rime_volc_ice(ii,jj,cc) + &
                        ptstep*zccip(ll,cc)*pprecp(ii,jj,ll)%volc(1)*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained from smaller ice particles in b
              DO ll = inp2b,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller ice particles in a
              DO ll = inp2a,kk-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zccii(ll,cc)*pice(ii,jj,ll)%volc(1:nt)
              END DO

              ! Save ice-ice collisions for collisional breakup
              IF (nlsip_iibr) THEN
                 ! Smaller from regime a
                 DO ll = inp2b,cc-1
                    coll_rate_ii(ii,jj,cc,ll) = &
                       ptstep*zccii(ll,cc)*pice(ii,jj,ll)%numc*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
                 ! Smaller and equal (self) from regime b
                 DO ll = inp2a,kk
                    coll_rate_ii(ii,jj,cc,ll) = &
                       ptstep*zccii(ll,cc)*pice(ii,jj,ll)%numc*pice(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
              ENDIF

              ! Update the hydrometeor volume concentrations
              pice(ii,jj,cc)%volc(1:nt) = max(0.,( pice(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*pice(ii,jj,cc)%numc ) /         &
                   (1. + ptstep*zminusterm) )

              ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
              pice(ii,jj,cc)%numc = max(0.,pice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                   0.5*ptstep*zccii(cc,cc)*pice(ii,jj,cc)%numc ) )

           END DO

           ! Snow
           ! -----------------------------------
           DO cc = 1,nsnw
              IF (psnow(ii,jj,cc)%numc<prlim) CYCLE

              zminusterm = 0.
              zplusterm(:) = 0.

              ! Drops lost by coagulation with larger snow
              DO ll = cc+1,nsnw
                 zminusterm = zminusterm + zccss(cc,ll)*psnow(ii,jj,ll)%numc
              END DO

              ! Volume gained by collection of aerosols
              DO ll = in1a,fn2b
                 zplusterm(1:nt) = zplusterm(1:nt) + zccsa(ll,cc)*paero(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdpart(ll)>hm_dmin_drop .AND. zdsnow(cc)>hm_dmin_ice) THEN
                    rime_volc_snw(ii,jj,cc) = rime_volc_snw(ii,jj,cc) + &
                        ptstep* zccsa(ll,cc)*paero(ii,jj,ll)%volc(1)*psnow(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained by collection of cloud droplets
              DO ll = 1,ncld
                 zplusterm(1:nt) = zplusterm(1:nt) + zccsc(ll,cc)*pcloud(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdcloud(ll)>hm_dmin_drop .AND. zdsnow(cc)>hm_dmin_ice) THEN
                    rime_volc_snw(ii,jj,cc) = rime_volc_snw(ii,jj,cc) + &
                        ptstep*zccsc(ll,cc)*pcloud(ii,jj,ll)%volc(1)*psnow(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained by collection of rain drops
              DO ll = 1,nprc
                 zplusterm(1:nt) = zplusterm(1:nt) + zccsp(ll,cc)*pprecp(ii,jj,ll)%volc(1:nt)
                 ! Save rime for Hallett-Mossop
                 IF (nlsip_hm .AND. zdprecp(ll)>hm_dmin_drop .AND. zdsnow(cc)>hm_dmin_ice) THEN
                    rime_volc_snw(ii,jj,cc) = rime_volc_snw(ii,jj,cc) + &
                        ptstep*zccsp(ll,cc)*pprecp(ii,jj,ll)%volc(1)*psnow(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 ENDIF
              END DO

              ! Volume gained by collection of ice particles
              DO ll = 1,nice
                 zplusterm(1:nt) = zplusterm(1:nt) + zccsi(ll,cc)*pice(ii,jj,ll)%volc(1:nt)
              END DO

              ! Volume gained from smaller snow
              DO ll = 1,cc-1
                 zplusterm(1:nt) = zplusterm(1:nt) + zccss(ll,cc)*psnow(ii,jj,ll)%volc(1:nt)
              END DO

              ! Save snow-ice and snow-snow collisions for collisional breakup
              IF (nlsip_iibr) THEN
                 ! Ice
                 DO ll = 1,nice
                    coll_rate_si(ii,jj,cc,ll) = &
                       ptstep*zccsi(ll,cc)*pice(ii,jj,ll)%numc*psnow(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
                 ! Smaller and equal (self) snow
                 DO ll = 1,cc
                    coll_rate_ss(ii,jj,cc,ll) = &
                       ptstep*zccss(ll,cc)*psnow(ii,jj,ll)%numc*psnow(ii,jj,cc)%numc/(1.+ptstep*zminusterm)
                 END DO
              END IF

              ! Update the hydrometeor volume concentrations
              psnow(ii,jj,cc)%volc(1:nt) = max(0.,( psnow(ii,jj,cc)%volc(1:nt) +  &
                   ptstep*zplusterm(1:nt)*psnow(ii,jj,cc)%numc ) /         &
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
  ! subroutine CONDENSATION(kbdim,  klev,                &
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
  !       Modified for the new aerosol datatype. LWC is obtained from %volc(1)
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

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE condgas(kbdim,  klev,                          &
                     paero,  pcloud, pprecp, pice,  psnow,  &
                     pc_gas, ngas,   ptemp,  ppres, ptstep  )

    USE mo_submctl,    ONLY :   &
         pi, avog, boltz, rg,       & ! constants
         in1a, fn2b,                & ! size bin indices
         t_section,                 & ! data type for the cloud bin representation
         ncld, nprc, nice, nsnw,    & ! number of bins
         nlim, prlim,               & ! concentration limits
         CalcDimension,             & ! function for updating wet sizes
         pstand,                    & ! standard pressure [Pa]
         msu,moc,rhosu,rhooc,       & ! molar mass [kg/mol] and density [kg/m3] of sulphate and OC
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         massacc,                   & ! mass accomodation coefficients in each bin
         n3,                        & ! number of molecules in one 3 nm particle [1]
         part_h2so4, iso, isog,     & ! secondary sulfate aerosol formation
         part_ocnv,  ioc, iocg        ! secondary organic aerosol formation: simple non-volatile species

   IMPLICIT NONE

    !-- Input and output variables ----------
    INTEGER, INTENT(IN) :: &
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         ngas                         ! number of gases

    REAL, INTENT(IN) :: &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep                       ! timestep [s]

    REAL, INTENT(INOUT) :: &
         pc_gas(kbdim,klev,ngas)      ! gas concentrations [mol/m3]

    TYPE(t_section), INTENT(inout) :: &
         paero(kbdim,klev,fn2b),    & ! aerosol properties
         pcloud(kbdim,klev,ncld),   & ! cloud properties
         pprecp(kbdim,klev,nprc),   & ! rain properties
         pice(kbdim,klev,nice),     & ! ice properties
         psnow(kbdim,klev,nsnw)       ! snow properties


    !-- Local variables ----------------------
    INTEGER :: ii, jj    ! loop indices

    REAL :: &
         zvisc,                      & ! viscosity of air [kg/(m s)]
         zdfvap,                     & ! air diffusion coefficient [m2/s]
         zmfp,                       & ! mean free path of condensing vapour [m]
         zcs_tot,                    & ! total condensation sink [1/s] (gases)
         zcs_nucl,                   & ! that including also nucleation [1/s]
         zcvap_new,                  & ! vapour concentration after time step [mol/m3]
         zdvap,                      & ! change in vapour concentration [mol/m3]
         zdfpart(in1a+1),            & ! particle diffusion coefficient (up to bin in1a+1)
         ! Knudsen number [-] for aerosols, cloud droplets, rain drops, ice and snow
         zknud(fn2b),zknca(ncld),zknpa(nprc),zknia(nice),zknsa(nsnw), &
         ! Transitional correction factor [-] for aerosols, cloud droplets, rain drops, ice and snow
         zbeta(fn2b),zbetaca(ncld),zbetapa(nprc),zbetaia(nice),zbetasa(nsnw), &
         ! Collision rate [1/s] of molecules to aerosols, cloud droplets, rain drops, ice and snow
         zcolrate(fn2b),zcolrateca(ncld),zcolratepa(nprc),zcolrateia(nice),zcolratesa(nsnw), &
         zcolrateaa(fn2b),           & ! collision rate of molecules to aerosol accounting for nucleation
         zj3n3(kbdim,klev,2),        & ! formation rate of molecules in nucleation [molec/m3s] (H2SO4 and organic vapor)
         zn_vs_c,                    & ! ratio of nucleation of all mass transfer in the smallest bin
         zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles
         zxocnv(kbdim,klev)

    ! Nucleation disabled!
    zj3n3 = 0.
    zxsa = 0.
    zxocnv = 0.

    DO jj = 1,klev
       DO ii = 1,kbdim
          ! Calculate wet diameters (%dwet)
          CALL CalcDimension(fn2b,paero(ii,jj,:),nlim,1)
          CALL CalcDimension(ncld,pcloud(ii,jj,:),nlim,2)
          CALL CalcDimension(nprc,pprecp(ii,jj,:),prlim,3)
          CALL CalcDimension(nice,pice(ii,jj,:),prlim,4)
          CALL CalcDimension(nsnw,psnow(ii,jj,:),prlim,5)

          !-- 1) Properties of air and condensing gases --------------------
          zvisc  = (7.44523e-3*ptemp(ii,jj)**1.5)/(5093.*(ptemp(ii,jj)+110.4))! viscosity of air [kg/(m s)]
          zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)          ! diffusion coefficient of H2SO4 in air [m2/s]
          zmfp   = 3.*zdfvap*sqrt(pi*msu/(8.*rg*ptemp(ii,jj)))                ! mean free path of H2SO4 [m]

          !-- 2) Transition regime correction factor for particles ---------
          !
          !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
          !  Topics in current aerosol research, Pergamon.
          !
          !  Size of condensing molecule considered only for
          !  nucleation mode (3 - 20 nm)
          !

          !-- Knudsen numbers
          ! Aerosol
          zknud(:) = 2.*zmfp/paero(ii,jj,:)%dwet
          ! the first two 1a bins
          zknud(in1a:in1a+1) = 2.*zmfp/(paero(ii,jj,in1a:in1a+1)%dwet+d_sa) ! note: this assumes sulfate
          ! Cloud droplets, rain drops, ice and snow
          zknca(:) = 2.*zmfp/pcloud(ii,jj,:)%dwet
          zknpa(:) = 2.*zmfp/pprecp(ii,jj,:)%dwet
          zknia(:) = 2.*zmfp/pice(ii,jj,:)%dwet
          zknsa(:) = 2.*zmfp/psnow(ii,jj,:)%dwet

          !-- transitional correction factors
          ! Aerosol
          zbeta = (zknud + 1.)/(0.377*zknud+1.+4./(3.*massacc)*(zknud+zknud**2))
          ! Cloud droplets, rain drops, ice and snow
          zbetaca = 1./( 1. + zknca*( 1.33 + (0.71/zknca) )/( 1. + (1./zknca) ) )
          zbetapa = 1./( 1. + zknpa*( 1.33 + (0.71/zknpa) )/( 1. + (1./zknpa) ) )
          zbetaia = 1./( 1. + zknia*( 1.33 + (0.71/zknia) )/( 1. + (1./zknia) ) )
          zbetasa = 1./( 1. + zknsa*( 1.33 + (0.71/zknsa) )/( 1. + (1./zknsa) ) )

          !-- 3) Collision rate of molecules to particles -------------------
          !
          !  Particle diffusion coefficient considered only for
          !  nucleation mode (3 - 20 nm)
          !

          !-- particle diffusion coefficient [m2/s] - the first two 1a bins
          zdfpart(in1a:in1a+1) = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/(3.*pi*zvisc*paero(ii,jj,in1a:in1a+1)%dwet)

          !-- collision rates [1/s]
          ! Aerosol
          zcolrate(:)   = MERGE( 2.*pi*paero(ii,jj,:)%dwet*zdfvap*zbeta(:)*paero(ii,jj,:)%numc,      &
                                 0., paero(ii,jj,:)%numc > nlim )
          ! the first two 1a bins
          zcolrate(in1a:in1a+1) = MERGE( 2.*pi*(paero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    &
                                 (zdfvap+zdfpart(in1a:in1a+1))*zbeta(in1a:in1a+1)*paero(ii,jj,in1a:in1a+1)%numc, &
                                 0.,paero(ii,jj,in1a:in1a+1)%numc > nlim )
          ! Cloud droplets, rain drops, ice and snow
          zcolrateca(:) = MERGE( 2.*pi*pcloud(ii,jj,:)%dwet*zdfvap*zbetaca(:)*pcloud(ii,jj,:)%numc,  &
                                 0., pcloud(ii,jj,:)%numc > nlim )
          zcolratepa(:) = MERGE( 2.*pi*pprecp(ii,jj,:)%dwet*zdfvap*zbetapa(:)*pprecp(ii,jj,:)%numc,  &
                                 0., pprecp(ii,jj,:)%numc > prlim )
          zcolrateia(:) = MERGE( 2.*pi*pice(ii,jj,:)%dwet*zdfvap*zbetaia(:)*pice(ii,jj,:)%numc,      &
                                 0., pice(ii,jj,:)%numc > prlim )
          zcolratesa(:) = MERGE( 2.*pi*psnow(ii,jj,:)%dwet*zdfvap*zbetasa(:)*psnow(ii,jj,:)%numc,    &
                                 0., psnow(ii,jj,:)%numc > prlim )

          !-- 4) Total condensation sink [1/s] -------------------------------------

          zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)+ sum(zcolrateia) + sum(zcolratesa)

          !-- 5) Changes in gas-phase concentrations and particle volume -----
          !
          !--- 5.1) Organic vapours ------------------------

          !---- 5.1.1) Non-volatile organic compound: condenses onto all bins
          IF(pc_gas(ii,jj,iocg) > 1.e-20 .and. zcs_tot > 1.e-30 .and. part_ocnv) THEN

             zn_vs_c = 0.
             ! zj3n3 is the formation rate of particles (#/m3/s)
             IF(zj3n3(ii,jj,2) > 1.) zn_vs_c = zj3n3(ii,jj,2)/(zj3n3(ii,jj,2) + &
                                              avog * pc_gas(ii,jj,iocg) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrateaa = zcolrate
             zcolrateaa(in1a) = zcolrateaa(in1a) + zj3n3(ii,jj,2)/pc_gas(ii,jj,iocg)/avog

             zcs_nucl = zcs_tot + zj3n3(ii,jj,2)/pc_gas(ii,jj,iocg)/avog  ! total sink for organic vapor [1/s]

             zcvap_new = pc_gas(ii,jj,iocg)/(1.+ptstep*zcs_nucl)          ! new gas phase concentration [mol/m3]
             zdvap = pc_gas(ii,jj,iocg) - zcvap_new                       ! change in gas concentration [mol/m3]
             pc_gas(ii,jj,iocg) = zcvap_new                               ! updating vapour concentration [mol/m3]

             ! Change in total particle volume concentration (m3/m3) divided by the total sink [1/s]
             zdvap = zdvap*(moc/rhooc)/zcs_nucl

             !-- Change of volume concentration of organics in aerosol, cloud droplets, rain drops, ice and snow
             paero(ii,jj,:)%volc(ioc)  = paero(ii,jj,:)%volc(ioc)  + zcolrateaa(:)*zdvap
             pcloud(ii,jj,:)%volc(ioc) = pcloud(ii,jj,:)%volc(ioc) + zcolrateca(:)*zdvap
             pprecp(ii,jj,:)%volc(ioc) = pprecp(ii,jj,:)%volc(ioc) + zcolratepa(:)*zdvap
             pice(ii,jj,:)%volc(ioc)   = pice(ii,jj,:)%volc(ioc)   + zcolrateia(:)*zdvap
             psnow(ii,jj,:)%volc(ioc)  = psnow(ii,jj,:)%volc(ioc)  + zcolratesa(:)*zdvap

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxocnv(ii,jj) > 0.) THEN
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc + &
                     zn_vs_c * zcolrateaa(in1a)*zdvap/(n3*zxocnv(ii,jj))
             END IF

          END IF


          ! ---- 5.2) Sulphate -------------------------------------------
          IF(pc_gas(ii,jj,isog) > 1.e-20 .and. zcs_tot > 1.e-30 .and. part_h2so4) THEN

             zn_vs_c = 0.
             IF(zj3n3(ii,jj,1) > 1.) zn_vs_c = zj3n3(ii,jj,1)/(zj3n3(ii,jj,1) + &
                                              avog * pc_gas(ii,jj,isog) * zcolrate(in1a))

             !   collision rate in the smallest bin, including nucleation and condensation
             !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
             !   equation (16.73)
             zcolrateaa = zcolrate
             zcolrateaa(in1a) = zcolrateaa(in1a) + zj3n3(ii,jj,1)/pc_gas(ii,jj,isog)/avog

             zcs_nucl = zcs_tot + zj3n3(ii,jj,1)/pc_gas(ii,jj,isog)/avog  ! total sink for sulfate [1/s]

             zcvap_new = pc_gas(ii,jj,isog)/(1.+ptstep*zcs_nucl)          ! new gas phase concentration [mol/m3]
             zdvap = pc_gas(ii,jj,isog) - zcvap_new                       ! change in gas concentration [mol/m3]
             pc_gas(ii,jj,isog) = zcvap_new                               ! updating vapour concentration [mol/m3]

             ! Change in total particle volume concentration (m3/m3) divided by the total sink [1/s]
             zdvap = zdvap*(msu/rhosu)/zcs_nucl

             !-- Change of volume concentration of sulphate in aerosol, cloud droplets, rain drops, ice and snow
             paero(ii,jj,:)%volc(iso)  = paero(ii,jj,:)%volc(iso)  + zcolrateaa(:)*zdvap
             pcloud(ii,jj,:)%volc(iso) = pcloud(ii,jj,:)%volc(iso) + zcolrateca(:)*zdvap
             pprecp(ii,jj,:)%volc(iso) = pprecp(ii,jj,:)%volc(iso) + zcolratepa(:)*zdvap
             pice(ii,jj,:)%volc(iso)   = pice(ii,jj,:)%volc(iso)   + zcolrateia(:)*zdvap
             psnow(ii,jj,:)%volc(iso)  = psnow(ii,jj,:)%volc(iso)  + zcolratesa(:)*zdvap

             !-- Change of number concentration in the smallest bin caused by nucleation
             !   Jacobson (2005), equation (16.75)
             IF (zxsa(ii,jj) > 0.) THEN
                paero(ii,jj,in1a)%numc = paero(ii,jj,in1a)%numc +          &
                     zn_vs_c * zcolrateaa(in1a)*zdvap/(n3*zxsa(ii,jj))
             END IF

          END IF

       END DO ! kbdim

    END DO ! klev

  END SUBROUTINE condgas

!
! ----------------------------------------------------------------------------------------------------------
!

  SUBROUTINE gpparth2o(kbdim,  klev,                &
                       paero,  pcloud, pprecp,      &
                       pice, psnow,                 &
                       ptemp,  ppres,  prs,prsi, prv,    &
                       ptstep)
    
    USE mo_submctl, ONLY : t_section,            &
                               nbins, ncld, nprc,    &
                               nice, nsnw,            &
                               rhowa, mwa, mair,     &
                               surfw0, surfi0, rg,           &
                               pi, pi6, prlim, nlim,      &
                               massacc,avog,  &
                               in1a,in2a,  &
                               fn2b,            &
                               lscndh2oae, lscndh2ocl, lscndh2oic, &
                               alv, als, CalcDimension
    USE mo_salsa_properties, ONLY : equilibration
    IMPLICIT NONE

    INTEGER, INTENT(in) :: kbdim,klev
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
    REAL :: zcwsurfae, zcwsurfcd, zcwsurfpd, zcwsurfid, zcwsurfsd ! Saturation mole concentrations
    REAL :: zmtae(nbins), zmtcd(ncld), zmtpd(nprc),      & ! Mass transfer coefficients
                zmtid(nice), zmtsd(nsnw)
    REAL :: zwsatae(nbins), zwsatcd(ncld), zwsatpd(nprc), &  ! Water saturation ratios above
                zwsatid(nice), zwsatsd(nsnw)                      ! ice'n'snow
    REAL :: zcwtot                                        ! Total water mole concentration
    REAL :: zcwc, zcwint                   ! Current and new water vapour mole concentrations
    REAL :: zcwcae(nbins), zcwintae(nbins) ! Current and new water mole concentrations in aerosols
    REAL :: zcwccd(ncld), zcwintcd(ncld)   !     -  ''  -     in cloud droplets
    REAL :: zcwcpd(nprc), zcwintpd(nprc)   !     -  ''  -     in rain drops
    REAL :: zcwcid(nice), zcwintid(nice)   !     -  ''  -     in ice particles
    REAL :: zcwcsd(nsnw), zcwintsd(nsnw)   !     -  ''  -     in snow particles
    REAL :: zdfh2o, zthcond,rhoair
    REAL :: zbeta,zknud,zmfph2o
    REAL :: zact, zhlp1,zhlp2,zhlp3
    REAL :: adt
    REAL :: dwet, cap
    REAL :: zrh(kbdim,klev)

    REAL :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

    INTEGER :: nstr,nstep,istep
    INTEGER :: ii,jj,cc
    LOGICAL aero_eq, any_aero, any_cloud, any_prec, any_ice, any_snow

    zrh(:,:) = prv(:,:)/prs(:,:)

    ! Calculate the condensation only for 2a/2b aerosol bins
    nstr = in2a

    ! Save the current aerosol water content
    zaelwc1(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(1),DIM=3)*rhowa

    ! If RH < 98 % or dynamic condensation for aerosol is switched off,
    ! do equilibrium for all aerosol bins, but otherwise just 1a.
    aero_eq = zrh(1,1) < 0.98 .OR. .NOT. lscndh2oae ! Equilibrium for 2a and 2b?
    CALL equilibration(kbdim,klev,zrh,ptemp,paero,aero_eq)

    ! The new aerosol water content after equilibrium calculation
    zaelwc2(:,:) = SUM(paero(:,:,in1a:fn2b)%volc(1),DIM=3)*rhowa

    prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )/( ppres(:,:)*mair/(rg*ptemp(:,:)) )

    ! Steps in the substepping loop (default adt=2.e-2)
    nstep=MAX(1,NINT(ptstep/2.e-2))
    adt=ptstep/REAL(nstep)

    DO jj = 1,klev
       DO ii = 1,kbdim

          any_aero = ANY(paero(ii,jj,:)%numc > nlim)
          any_cloud = ANY(pcloud(ii,jj,:)%numc > nlim)
          any_prec = ANY(pprecp(ii,jj,:)%numc > prlim)
          any_ice = ANY(pice(ii,jj,:)%numc > prlim)
          any_snow = ANY(psnow(ii,jj,:)%numc > prlim)

          IF ( .NOT. ( &
                ((any_cloud .OR. any_prec) .AND. lscndh2ocl) .OR. &
                ((any_ice .OR. any_snow) .AND. lscndh2oic) .OR. &
                (any_aero .AND. .NOT.aero_eq) &
                ) ) CYCLE

          rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

          zdfh2o = ( 5./(16.*avog*rhoair*1.e-3*(3.11e-8)**2) ) * &
                   SQRT( rg*1e7*ptemp(ii,jj)*mair*1.e3*(mwa+mair)*1.e3/( 2.*pi*mwa*1.e3 ) )
          zdfh2o = zdfh2o*1.e-4

          zmfph2o = 3.*zdfh2o*sqrt(pi*mwa/(8.*rg*ptemp(ii,jj)))
          zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air

          ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
          zkelvinpd = 1.; zkelvincd = 1.; zkelvin = 1.; zkelvinid = 1.; zkelvinsd = 1.
          zwsatae = 0.; zwsatcd = 0.; zwsatpd = 0.; zwsatid = 0.; zwsatsd = 0.
          zmtae = 0.; zmtcd = 0.; zmtpd = 0.; zmtid = 0.; zmtsd = 0.

          ! Cloud droplets --------------------------------------------------------------------------------
          ! Saturation mole concentration over flat surface
          zcwsurfcd  = prs(ii,jj)*rhoair/mwa
          DO cc = 1,ncld
             IF (pcloud(ii,jj,cc)%numc > nlim .AND. lscndh2ocl) THEN
                ! Wet diameter
                dwet=( SUM(pcloud(ii,jj,cc)%volc(:))/pcloud(ii,jj,cc)%numc/pi6 )**(1./3.)

                ! Activity + Kelvin effect
                zact = acth2o(pcloud(ii,jj,cc))
                zkelvincd(cc) = exp( 4.*surfw0*mwa / (rg*ptemp(ii,jj)*rhowa*dwet) )

                ! Equilibrium saturation ratio
                zwsatcd(cc) = zact*zkelvincd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pcloud(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Rain drops --------------------------------------------------------------------------------
          ! Saturation mole concentration over flat surface
          zcwsurfpd = prs(ii,jj)*rhoair/mwa
          DO cc = 1,nprc
             IF (pprecp(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN
                ! Wet diameter
                dwet=( SUM(pprecp(ii,jj,cc)%volc(:))/pprecp(ii,jj,cc)%numc/pi6 )**(1./3.)

                ! Activity + Kelvin effect
                zact = acth2o(pprecp(ii,jj,cc))
                zkelvinpd(cc) = exp( 4.*surfw0*mwa / (rg*ptemp(ii,jj)*rhowa*dwet) )

                ! Equilibrium saturation ratio
                zwsatpd(cc) = zact*zkelvinpd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pprecp(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Ice particles --------------------------------------------------------------------------------
          ! Saturation mole concentration over flat surface
          zcwsurfid = prsi(ii,jj)*rhoair/mwa
          DO cc = 1,nice
             IF (pice(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                ! Dimension
                CALL CalcDimension(1,pice(ii,jj,cc),prlim,4)
                dwet=pice(ii,jj,cc)%dwet

                ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                cap=0.5*dwet

                ! Activity + Kelvin effect - edit when needed
                !   Can be calculated just like for sperical homogenous particle or just ignored,
                !   because these are not known for solid, irregular and non-homogenous particles.
                !   Ice may not be that far from a sphere, but most particles are large and at least
                !   growing particles are covered by a layer of pure ice.
                zact = 1.0 ! Note: acth2o does not work for ice or snow!
                zkelvinid(cc) = exp( 4.*surfi0*mwa / (rg*ptemp(ii,jj)*rhowa*dwet) )

                ! Equilibrium saturation ratio
                zwsatid(cc) = zact*zkelvinid(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = pice(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*als*zwsatid(cc)*zcwsurfid/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (als*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtid(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Snow particles --------------------------------------------------------------------------------
          ! Saturation mole concentration over flat surface
          zcwsurfsd= prsi(ii,jj)*rhoair/mwa
          DO cc = 1,nsnw
             IF (psnow(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                ! Dimension
                CALL CalcDimension(1,psnow(ii,jj,cc),prlim,5)
                dwet=psnow(ii,jj,cc)%dwet

                ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                cap=0.5*dwet

                ! Activity + Kelvin effect
                !   Can be calculated just like for sperical homogenous particle or just ignored,
                !   because these are not known for solid, irregular and non-homogenous particles.
                !   Especially snow is typically highly irregular (e.g. dendrite).
                zact = 1.0 ! Note: acth2o does not work for ice or snow!
                zkelvinsd(cc) = exp( 4.*surfi0*mwa / (rg*ptemp(ii,jj)*rhowa*dwet) )

                ! Equilibrium saturation ratio
                zwsatsd(cc) = zact*zkelvinsd(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.)*(zknud+zknud**2))

                ! Mass transfer according to Jacobson
                zhlp1 = psnow(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*als*zwsatsd(cc)*zcwsurfsd/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (als*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtsd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! -- Aerosols: ------------------------------------------------------------------------------------
          ! Saturation mole concentration over flat surface
          ! Limit the supersaturation to max 1.01 for the mass transfer EXPERIMENTAL
          zcwsurfae =MAX(prs(ii,jj),prv(ii,jj)/1.01)*rhoair/mwa
          DO cc = nstr,nbins
             IF (paero(ii,jj,cc)%numc > nlim .AND. .NOT.aero_eq) THEN
                ! Wet diameter
                dwet=( SUM(paero(ii,jj,cc)%volc(:))/paero(ii,jj,cc)%numc/pi6 )**(1./3.)

                ! Water activity + Kelvin effect
                zact = acth2o(paero(ii,jj,cc))
                zkelvin(cc) = exp( 4.*surfw0*mwa / (rg*ptemp(ii,jj)*rhowa*dwet) )

                ! Equilibrium saturation ratio
                zwsatae(cc) = zact*zkelvin(cc)

                !-- transitional correction factor
                zknud = 2.*zmfph2o/dwet
                zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                     (3.*massacc(cc))*(zknud+zknud**2))

                ! Mass transfer
                zhlp1 = paero(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                zhlp2 = mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae/(zthcond*ptemp(ii,jj))
                zhlp3 = ( (alv*mwa)/(rg*ptemp(ii,jj)) ) - 1.

                zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

             END IF
          END DO

          ! Current mole concentrations
          zcwc = prv(ii,jj)*rhoair/mwa
          zcwcae(1:nbins) = paero(ii,jj,1:nbins)%volc(1)*rhowa/mwa
          zcwccd(1:ncld) = pcloud(ii,jj,1:ncld)%volc(1)*rhowa/mwa
          zcwcpd(1:nprc) = pprecp(ii,jj,1:nprc)%volc(1)*rhowa/mwa
          zcwcid(1:nice) = pice(ii,jj,1:nice)%volc(1)*rhowa/mwa
          zcwcsd(1:nsnw) = psnow(ii,jj,1:nsnw)%volc(1)*rhowa/mwa

          zcwtot = zcwc + SUM(zcwcae) + &
                          SUM(zcwccd) + &
                          SUM(zcwcpd) + &
                          SUM(zcwcid) + &
                          SUM(zcwcsd)

          zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd; zcwintid = zcwcid; zcwintsd = zcwcsd

          ! Substepping loop
          ! ---------------------------------
          DO istep=1,nstep

             ! New vapor concentration
             zhlp1 = zcwc + adt * ( SUM(zmtae(nstr:nbins)*zwsatae(nstr:nbins)*zcwsurfae) + &
                                    SUM(zmtcd(1:ncld)*zwsatcd(1:ncld)*zcwsurfcd)         + &
                                    SUM(zmtpd(1:nprc)*zwsatpd(1:nprc)*zcwsurfpd)         + &
                                    SUM(zmtid(1:nice)*zwsatid(1:nice)*zcwsurfid)         + &
                                    SUM(zmtsd(1:nsnw)*zwsatsd(1:nsnw)*zcwsurfsd)         )

             zhlp2 = 1. + adt * ( SUM(zmtae(nstr:nbins)) + SUM(zmtcd(1:ncld)) + SUM(zmtpd(1:nprc)) &
                                   + SUM(zmtid(1:nice)) + SUM(zmtsd(1:nsnw)) )
             zcwint = zhlp1/zhlp2
             zcwint = MIN(zcwint,zcwtot)

             IF ( any_aero .AND. .NOT.aero_eq ) THEN
                DO cc = nstr,nbins
                   zcwintae(cc) = zcwcae(cc) + min(max(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae), &
                        -0.02*zcwcae(cc)),0.05*zcwcae(cc))
                   zwsatae(cc) = acth2o(paero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                END DO
                zcwintae(nstr:nbins) = MAX(zcwintae(nstr:nbins),0.)
             END IF
             IF ( any_cloud ) THEN
                DO cc = 1,ncld
                   zcwintcd(cc) = zcwccd(cc) + min(max(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd), &
                        -0.02*zcwccd(cc)),0.05*zcwccd(cc))
                   zwsatcd(cc) = acth2o(pcloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                END DO
                zcwintcd(1:ncld) = MAX(zcwintcd(1:ncld),0.)
             END IF
             IF ( any_prec ) THEN
                DO cc = 1,nprc
                   zcwintpd(cc) = zcwcpd(cc) + min(max(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd), &
                        -0.02*zcwcpd(cc)),0.05*zcwcpd(cc))
                   zwsatpd(cc) = acth2o(pprecp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                END DO
                zcwintpd(1:nprc) = MAX(zcwintpd(1:nprc),0.)
             END IF
             IF ( any_ice ) THEN
                DO cc = 1,nice
                   zcwintid(cc) = zcwcid(cc) + min(max(adt*zmtid(cc)*(zcwint - zwsatid(cc)*zcwsurfid), &
                        -0.02*zcwcid(cc)),0.05*zcwcid(cc))
                   zwsatid(cc) = zkelvinid(cc)
                END DO
                zcwintid(1:nice) = MAX(zcwintid(1:nice),0.)
             END IF
             IF ( any_snow ) THEN
                DO cc = 1,nsnw
                   zcwintsd(cc) = zcwcsd(cc) + min(max(adt*zmtsd(cc)*(zcwint - zwsatsd(cc)*zcwsurfsd),&
                        -0.02*zcwcsd(cc)),0.05*zcwcsd(cc))
                   zwsatsd(cc) = zkelvinsd(cc)
                END DO
                zcwintsd(1:nsnw) = MAX(zcwintsd(1:nsnw),0.)
             END IF

             ! Update vapor concentration for consistency
             zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                               SUM(zcwintcd(1:ncld))  - &
                               SUM(zcwintpd(1:nprc))  - &
                               SUM(zcwintid(1:nice))  - &
                               SUM(zcwintsd(1:nsnw))

             ! Update "old" values for next cycle
             zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd; zcwcid =zcwintid; zcwcsd = zcwintsd;
             zcwc = zcwint

          END DO ! ADT

          prv(ii,jj) = zcwint*mwa/rhoair

          paero(ii,jj,1:nbins)%volc(1) = max(0.,zcwintae(1:nbins)*mwa/rhowa)
          pcloud(ii,jj,1:ncld)%volc(1) = max(0.,zcwintcd(1:ncld)*mwa/rhowa)
          pprecp(ii,jj,1:nprc)%volc(1) = max(0.,zcwintpd(1:nprc)*mwa/rhowa)
          pice(ii,jj,1:nice)%volc(1) = max(0.,zcwintid(1:nice)*mwa/rhowa)
          psnow(ii,jj,1:nsnw)%volc(1) = max(0.,zcwintsd(1:nsnw)*mwa/rhowa)

       END DO !kbdim

    END DO ! klev

  END SUBROUTINE gpparth2o
  !-------------------------------------------------------
  REAL FUNCTION acth2o(ppart,pcw)

    USE mo_submctl, ONLY : t_section, eps, diss, dens, mws
    IMPLICIT NONE

    TYPE(t_section), INTENT(in) :: ppart
    REAL, INTENT(in), OPTIONAL :: pcw

    REAL :: zns, znw

    zns = SUM( ppart%volc(2:)*diss(2:)*dens(2:)/mws(2:) )

    IF (PRESENT(pcw)) THEN
       znw = pcw
    ELSE
       znw = ppart%volc(1)*dens(1)/mws(1) ! This is not valid for ice and snow!
    END IF

    ! Assume activity coefficient of 1 for water and that there is always some soluble material
    acth2o = MAX(0.1,znw/max(eps,(znw+zns)))
  END FUNCTION acth2o

!
! ----------------------------------------------------------------------------------------------------------
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
  REAL FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres,eddy_dis,flag1,flag2)

    USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav, rda, terminal_vel

    IMPLICIT NONE

    !-- Input variables ----------
    REAL, INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres,   &   ! ambient pressure [fxm]
         eddy_dis    ! eddy dissipation rate

    INTEGER, INTENT(in) :: flag1,flag2 ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)

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
         zgrav  ,&   !                                    Gravitational collection
         ztshear, &     ! turbulent shear
         zturbinert  ! turbulent inertia
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
    INTEGER :: lrg,sml

    coagc = 0.

    ! Mass 1e-30 kg means no particle(s)
    IF (mass1<1e-29 .OR. mass2<1e-29) RETURN

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


          ! Which particle is larger?
          sml = 1; lrg = 2
          IF (diam(1) >= diam(2)) THEN
             lrg = 1; sml = 2
          END IF

          zrhoa = pres/(rda*temp)   ! Density of air
          zrhop = mpart/(pi6*diam**3)             ! Density of particles
          vkin = visc/zrhoa   ! Kinematic viscosity of air [m2 s-1]

          termv(1) = terminal_vel(diam(1)/2.,zrhop(1),zrhoa,visc,beta(1),flag1)
          termv(2) = terminal_vel(diam(2)/2.,zrhop(2),zrhoa,visc,beta(2),flag2)

          ! Reynolds number
          reyn = diam*termv/vkin
          ! Schmidt number for the smaller particle
          schm = vkin/dfpart
          ! Stokes number
          stok = 2.*termv(sml)*ABS(termv(1) - termv(2))/( diam(lrg)*grav )

          !Brownian component
          zbrown = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))

          ! Convective enhancement
          zbrconv = 0.
          IF (reyn(lrg) <= 1.) THEN
             zbrconv = 0.45*zbrown*( reyn(lrg)**(1./3.) )*( schm(sml)**(1./3.) )
          ELSE IF (reyn(lrg) > 1.) THEN
             zbrconv = 0.45*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1./3.) )
          END IF

          ! Turbulent Shear
          ztshear=(8.*pi*eddy_dis/(15.*vkin))**(1./2.)*(0.5*(diam(1)+diam(2)))**3.
          ! Turbulent inertial motion
          zturbinert = pi*eddy_dis**(3./4.)/(grav*vkin**(1./4.)) &
               *(0.5*(diam(1)+diam(2)))**2.* ABS(termv(1)-termv(2))

          ! gravitational collection
          zea = stok**2/( stok + 0.5 )**2
          zev = 0.
          IF (stok > 1.214) THEN
             zev = 0.75*LOG(2.*stok)/(stok - 1.214)
             zev = (1. + zev)**(-2.)
          END IF

          zecoll = (60.*zev + zea*reyn(lrg))/(60. + reyn(lrg))
          zgrav = zecoll * pi * mdiam**2
          zgrav = zgrav * ABS(termv(1)-termv(2))

          ! Total coagulation kernel
          coagc = zbrown  + zbrconv + (zgrav**2+ ztshear**2+ zturbinert**2)**(1./2.)

  END FUNCTION coagc



END MODULE mo_salsa_dynamics
