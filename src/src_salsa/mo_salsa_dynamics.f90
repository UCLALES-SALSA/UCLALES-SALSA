
!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_DYNAMICS                               *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************

MODULE mo_salsa_dynamics
   IMPLICIT NONE


CONTAINS

   ! this calculated for empty bins too!!!
   ! fxm: test well, esp. self-coagulation (but other bits too!)
   ! AL_note: Diagnostic variables of cond and nucl mass
   !********************************************************************
   !
   ! Subroutine COAGULATION(kproma,kbdim,klev, &
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
   !           Process selection should be made smarter - now just lots of ifs
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


   SUBROUTINE coagulation(kproma,kbdim,  klev,    &
                          allSALSA,  &
                          ptstep, ptemp,  ppres    )

      USE classSection
      USE mo_submctl, ONLY:        &
         t_parallelbin,   & ! Datatypes for the cloud bin representation
         in1a, fn1a,                 & ! size bin indices
         in2a, fn2a,                 &
         in2b, fn2b,                 &
         ica,fca,icb,fcb,            &
         ncld, nprc,                 &
         iia,fia,iib,fib,            &
         nice, nsnw,                 &
         ntotal,                     &
         spec,                       &
         pi6,                        &
         nlim,prlim,                 &
         lscgaa, lscgcc, lscgca,     &
         lscgpp, lscgpa, lscgpc,     &
         lscgia, lscgic, lscgii, lscgip, &
         lscgsa, lscgsc, lscgsi, lscgsp, lscgss, &
         CalcDimension,              &
         aero, cloud, precp, ice, snow

      IMPLICIT NONE


      !-- Input and output variables -------------
      INTEGER, INTENT(IN) ::        &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical klev

      TYPE(Section), INTENT(inout) :: &
         allSALSA(kbdim,klev,ntotal)  ! Size distributions for all particle types

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
         zccpp(nprc,nprc),          & ! - '' - for collisions between precip particles 
         zccia(fn2b,nice),          & ! - '' - for collection of aerosols by ice 
         zccic(ncld,nice),          & ! - '' - for collection of cloud particles droplets by ice 
         zccii(nice,nice),          & ! - '' - for collisions between ice 
         zccip(nprc,nice),          & ! - '' - for collection of precip by ice
         zccsa(fn2b,nsnw),          & ! - '' - for collection of aerosols by snow 
         zccsc(ncld,nsnw),          & ! - '' - for collection of cloud droples by snow 
         zccsi(nice,nsnw),          & ! - '' - for collection of ice by snow 
         zccsp(nprc,nsnw),          & ! - '' - for collection of precip by snow 
         zccss(nsnw,nsnw),          & ! - '' - for collisions between snow
         zminusterm                   ! coagulation loss in a bin [1/s]
      REAL, ALLOCATABLE ::          &
         zplusterm(:)                 ! coagulation gain in a bin [fxm/s]
                                      ! (for each chemical compound)

      REAL :: &
         zmpart(fn2b),   & ! approximate mass of particles [kg]
         zmcloud(ncld),  &    ! approximate mass of cloud droplets [kg]
         zmprecp(nprc),  & ! Approximate mass for rain drops [kg]
         zmice(nice),     & ! approximate mass for ice particles [kg] 
         zmsnow(nsnw), &  ! approximate mass for snow particles [kg] 
         zdpart(fn2b),   & ! diameter of particles [m]
         zdcloud(ncld),  &   ! diameter of cloud droplets [m]
         zdprecp(nprc),  & ! diameter for rain drops [m]
         zdice(nice),    & ! diameter for ice [m]
         zdsnow(nsnw)      ! diameter for snow [m]

      REAL :: temppi,pressi

      LOGICAL :: any_cloud, any_precp, any_ice, any_snow

      INTEGER :: nspec, ndry, iwa  ! Total number of compounds, number of dry compounds (shoudl be total-1) and index of water

      nspec = spec%getNSpec(type="wet")
      ndry = spec%getNSpec(type="dry")
      iwa = spec%getIndex("H2O")

      ALLOCATE(zplusterm(nspec))
      zplusterm = 0.

      !-----------------------------------------------------------------------------
      !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
      !      CoagSink ~ Dp in continuum regime, thus we calculate
      !      'effective' number concentration of coarse particles


      !-- 2) Updating coagulation coefficients -------------------------------------

      DO jj = 1, klev      ! vertical grid
         DO ii = 1, kbdim ! horizontal kproma in the slab

            ! Which species are included
            any_cloud = ANY(cloud(ii,jj,:)%numc > nlim)
            any_precp = ANY(precp(ii,jj,:)%numc > prlim)
            any_ice = ANY(ice(ii,jj,:)%numc > prlim)
            any_snow = ANY(snow(ii,jj,:)%numc > prlim)

            !-- Aerosol diameter [m] and mass [kg]; density of            1500 kg/m3 assumed
            CALL CalcDimension(fn2b,aero(ii,jj,1:fn2b),nlim,zdpart(1:fn2b),1)
            zdpart(1:fn2b) = MIN(zdpart(1:fn2b), 30.e-6) ! Limit to 30 um
            zmpart(1:fn2b) = pi6*(zdpart(1:fn2b)**3)*1500.

              !-- Cloud droplet diameter and mass; Assume water density
            CALL CalcDimension(ncld,cloud(ii,jj,1:ncld),nlim,zdcloud(1:ncld),2)
            ! No size limit?
            zmcloud(1:ncld) = pi6*(zdcloud(1:ncld)**3)*spec%rhowa

             !-- Precipitation droplet diameter and mass
            CALL CalcDimension(nprc,precp(ii,jj,1:nprc),prlim,zdprecp(1:nprc),3)
            zdprecp(1:nprc) = MIN(zdprecp(1:nprc), 2.e-3) ! Limit to 2 mm
            zmprecp(1:nprc) = pi6*(zdprecp(1:nprc)**3)*spec%rhowa

             !-- Ice particle diameter and mass
            CALL CalcDimension(nice,ice(ii,jj,1:nice),prlim,zdice(1:nice),4)
            zdice(1:nice) = MIN(zdice(1:nice), 2.e-3) ! Limit to 2 mm
            zmice(1:nice) = pi6*(zdice(1:nice)**3)*spec%rhoic

             !-- Snow diameter and mass
            CALL CalcDimension(nsnw,snow(ii,jj,1:nsnw),prlim,zdsnow(1:nsnw),5)
            zdsnow(1:nsnw) = MIN(zdsnow(1:nsnw), 10.e-3) ! Limit to 10 mm 
            zmsnow(1:nsnw) = pi6*(zdsnow(1:nsnw)**3)*spec%rhosn

            temppi = ptemp(ii,jj)
            pressi = ppres(ii,jj)
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
               DO mm = 1, fn2b         ! smaller colliding particle
                  IF (aero(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = mm, fn2b            ! larger colliding particle
                     IF (aero(ii,jj,nn)%numc < nlim) CYCLE
                     zcc(mm,nn) = coagc(zdpart(mm),zdpart(nn),zmpart(mm),zmpart(nn),temppi,pressi,1,1,1)
                     zcc(nn,mm) = zcc(mm,nn)
                  END DO
               END DO
            END IF

            ! Collision-coalescence between cloud droplets
            IF (lscgcc .AND. any_cloud) THEN
               DO mm = 1, ncld
                  IF (cloud(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = mm, ncld
                     IF (cloud(ii,jj,nn)%numc < nlim) CYCLE
                     zcccc(mm,nn) = coagc(zdcloud(mm),zdcloud(nn),zmcloud(mm),zmcloud(nn),temppi,pressi,2,2,2)
                     zcccc(nn,mm) = zcccc(mm,nn)                  
                  END DO
               END DO
            END IF

            ! Self-collection of rain drops
            IF (lscgpp .AND. any_precp) THEN
               DO mm = 1, nprc
                  IF (precp(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = mm, nprc
                     IF (precp(ii,jj,nn)%numc < prlim) CYCLE
                     zccpp(mm,nn) = coagc(zdprecp(mm),zdprecp(nn),zmprecp(mm),zmprecp(nn),temppi,pressi,2,3,3)
                     zccpp(nn,mm) = zccpp(mm,nn)
                  END DO
               END DO
            END IF

            ! Cloud collection of aerosols
            IF (lscgca .AND. any_cloud) THEN
               DO mm = 1, fn2b
                  IF (aero(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, ncld
                     IF (cloud(ii,jj,nn)%numc < nlim) CYCLE
                     zccca(mm,nn) = coagc(zdpart(mm),zdcloud(nn),zmpart(mm),zmcloud(nn),temppi,pressi,2,1,2)
                  END DO
               END DO
            END IF

            ! Collection of aerosols by rain
            IF (lscgpa .AND. any_precp) THEN
               DO mm = 1, fn2b
                  IF (aero(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nprc
                     IF (precp(ii,jj,nn)%numc < prlim) CYCLE
                     zccpa(mm,nn) = coagc(zdpart(mm),zdprecp(nn),zmpart(mm),zmprecp(nn),temppi,pressi,2,1,3)
                  END DO
               END DO
            END IF

            ! Collection of cloud droplets by rain
            IF (lscgpc .AND. any_cloud .AND. any_precp) THEN
               DO mm = 1, ncld
                  IF (cloud(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nprc
                     IF (precp(ii,jj,nn)%numc < prlim) CYCLE
                     zccpc(mm,nn) = coagc(zdcloud(mm),zdprecp(nn),zmcloud(mm),zmprecp(nn),temppi,pressi,2,2,3)
                  END DO
               END DO
            END IF

            !  collection of aerosols by ice
            IF (lscgia .AND. any_ice) THEN
               DO mm = 1, fn2b
                  IF (aero(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nice
                     IF (ice(ii,jj,nn)%numc < prlim) CYCLE
                     zccia(mm,nn) = coagc(zdpart(mm),zdice(nn),zmpart(mm),zmice(nn),temppi,pressi,2,1,4)
                  END DO
               END DO
            END IF

            !  collection of cloud particles droplets by ice
            IF (lscgic .AND. any_ice .AND. any_cloud) THEN
               DO mm = 1, ncld
                  IF (cloud(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nice
                     IF (ice(ii,jj,nn)%numc < prlim) CYCLE
                     zccic(mm,nn) = coagc(zdcloud(mm),zdice(nn),zmcloud(mm),zmice(nn),temppi,pressi,2,2,4)
                  END DO
               END DO
            END IF

            !  collisions between ice particles
            IF (lscgii .AND. any_ice) THEN
               DO mm = 1, nice
                  IF (ice(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = mm, nice
                     IF (ice(ii,jj,nn)%numc < prlim) CYCLE
                     zccii(mm,nn) = coagc(zdice(mm),zdice(nn),zmice(mm),zmice(nn),temppi,pressi,2,4,4)
                     zccii(nn,mm) = zccii(mm,nn)
                  END DO
               END DO
            END IF

            !  collection of precip by ice-collision
            IF (lscgip .AND. any_precp .AND. any_ice) THEN
               DO mm = 1, nprc
                  IF (precp(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = 1, nice
                     IF (ice(ii,jj,nn)%numc < prlim) CYCLE
                     zccip(mm,nn) = coagc(zdprecp(mm),zdice(nn),zmprecp(mm),zmice(nn),temppi,pressi,2,3,4)
                  END DO
               END DO
            END IF

            ! Self-collection of snow particles
            IF (lscgss .AND. any_snow) THEN
               DO mm = 1, nsnw
                  IF (snow(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = mm, nsnw
                     IF (snow(ii,jj,nn)%numc < prlim) CYCLE
                     zccss(mm,nn) = coagc(zdsnow(mm),zdsnow(nn),zmsnow(mm),zmsnow(nn),temppi,pressi,2,5,5)
                     zccss(nn,mm) = zccss(mm,nn)
                  END DO
               END DO
            END IF

            ! Collection of aerosols by snow
            IF (lscgsa .AND. any_snow) THEN
               DO mm = 1, fn2b
                  IF (aero(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nsnw
                     IF (snow(ii,jj,nn)%numc < prlim) CYCLE
                     zccsa(mm,nn) = coagc(zdpart(mm),zdsnow(nn),zmpart(mm),zmsnow(nn),temppi,pressi,2,1,5)
                  END DO
               END DO
            END IF

            ! collection of precip by snow
            IF (lscgsp .AND. any_precp .AND. any_snow) THEN
               DO mm = 1, nprc
                  IF (precp(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = 1, nsnw
                     IF (snow(ii,jj,nn)%numc < prlim) CYCLE
                     zccsp(mm,nn) = coagc(zdprecp(mm),zdsnow(nn),zmprecp(mm),zmsnow(nn),temppi,pressi,2,3,5)
                  END DO
               END DO
            END IF

            ! collection of cloud droples by snow
            IF (lscgsc .AND. any_cloud .AND. any_snow) THEN
               DO mm = 1, ncld
                  IF (cloud(ii,jj,mm)%numc < nlim) CYCLE
                  DO nn = 1, nsnw
                     IF (snow(ii,jj,nn)%numc < prlim) CYCLE
                     zccsc(mm,nn) = coagc(zdcloud(mm),zdsnow(nn),zmcloud(mm),zmsnow(nn),temppi,pressi,2,2,5)
                  END DO
               END DO
            END IF

            ! collection of ice by snow
            IF (lscgsi .AND. any_ice .AND. any_snow) THEN
               DO mm = 1, nice
                  IF (ice(ii,jj,mm)%numc < prlim) CYCLE
                  DO nn = 1, nsnw
                     IF (snow(ii,jj,nn)%numc < prlim) CYCLE
                     zccsi(mm,nn) = coagc(zdice(mm),zdsnow(nn),zmice(mm),zmsnow(nn),temppi,pressi,2,4,5)
                  END DO
               END DO
            END IF

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !-- 3) New particle and volume concentrations after coagulation -------------


            ! Rain drops
            ! -----------------------------------
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc < prlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! Drops lost by coagulation with larger drops
               DO ll = cc+1, nprc
                  zminusterm = zminusterm + zccpp(cc,ll)*precp(ii,jj,ll)%numc
               END DO

               ! Drops lost by collection by snow drops
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsp(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Drops lost by collisions with ice
               DO ll = 1, nice
                  zminusterm = zminusterm + zccip(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Volume gained by collection of aerosols
               DO ll = in1a, fn2b
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccpa(ll,cc)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained by collection of cloud droplets
               DO ll = 1, ncld
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccpc(ll,cc)*cloud(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller drops
               DO ll = 1, cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccpp(ll,cc)*precp(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               precp(ii,jj,cc)%volc(1:nspec) = max(0., ( precp(ii,jj,cc)%volc(1:nspec) +  &
                                                      ptstep*zplusterm(1:nspec)*precp(ii,jj,cc)%numc ) / &
                                                      (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
               precp(ii,jj,cc)%numc = max(0.,precp(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                           0.5*ptstep*zccpp(cc,cc)*precp(ii,jj,cc)%numc ) )

            END DO


            ! Aerosols in regime 1a
            ! --------------------------------
            DO kk = in1a, fn1a
               IF (aero(ii,jj,kk)%numc < nlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.
               ! Particles lost by coagulation with larger aerosols
               DO ll = kk+1, fn2b
                  zminusterm = zminusterm + zcc(kk,ll)*aero(ii,jj,ll)%numc
               END DO

               ! Particles lost by cloud collection
               DO ll = 1, ncld
                  zminusterm = zminusterm + zccca(kk,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! particles lost by rain collection
               DO ll = 1, nprc
                  zminusterm = zminusterm + zccpa(kk,ll)*precp(ii,jj,ll)%numc
               END DO

               ! particles lost by ice collection
               DO ll = 1, nice
                  zminusterm = zminusterm + zccia(kk,ll)*ice(ii,jj,ll)%numc
               END DO

               ! particles lost by snow collection
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsa(kk,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Particle volume gained from smaller particles in regime 1a
               DO ll = in1a, kk-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcc(ll,kk)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               !-- Volume and number concentrations after coagulation update [fxm]
               aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+ptstep*zplusterm(1:nspec) * &
                                             aero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

               aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                                                            0.5*ptstep*zcc(kk,kk)*aero(ii,jj,kk)%numc)

            END DO

            ! Aerosols in regime 2a
            ! ---------------------------------
            DO kk = in2a, fn2a
               IF (aero(ii,jj,kk)%numc < nlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! Find corresponding size bin in subregime 2b
               index_2b = kk - in2a + in2b

               ! Particles lost by larger particles in 2a
               DO ll = kk+1, fn2a
                  zminusterm = zminusterm + zcc(kk,ll)*aero(ii,jj,ll)%numc ! 2a
               END DO

               ! Particles lost by larger particles in 2b
               DO ll = index_2b+1, fn2b
                  zminusterm = zminusterm + zcc(kk,ll)*aero(ii,jj,ll)%numc ! 2b
               END DO

               ! Particles lost by cloud collection
               DO ll = 1, ncld
                  zminusterm = zminusterm + zccca(kk,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Particles lost by collection by rain
               DO ll = 1, nprc
                  zminusterm = zminusterm + zccpa(kk,ll)*precp(ii,jj,ll)%numc
               END DO

               ! particles lost by ice collection
               DO ll = 1, nice
                  zminusterm = zminusterm + zccia(kk,ll)*ice(ii,jj,ll)%numc
               END DO

               ! particles lost by snow collection
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsa(kk,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Particle volume gained from smaller particles in regimes 1, 2a
               DO ll = in1a, kk-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcc(ll,kk)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Particle volume gained from smaller (and equal) particles in 2b
               DO ll = in2b, index_2b
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcc(ll,kk)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               !-- Volume and number concentrations after coagulation update [fxm]
               aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+ptstep*zplusterm(1:nspec) *  &
                                             aero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

               aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                                                            0.5*ptstep*zcc(kk,kk)*aero(ii,jj,kk)%numc)

            END DO

            ! Aerosols in regime 2b
            ! ---------------------------------
            DO kk = in2b, fn2b
               IF (aero(ii,jj,kk)%numc < nlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               !-- Find corresponding size bin in subregime 2a
               index_2a = kk - in2b + in2a

               ! Particles lost to larger particles in regimes 2b
               DO ll = kk+1, fn2b
                  zminusterm = zminusterm + zcc(kk,ll)*aero(ii,jj,ll)%numc ! 2b
               END DO

               ! Particles lost to larger and equal particles in 2a
               DO ll = index_2a, fn2a
                  zminusterm = zminusterm + zcc(kk,ll)*aero(ii,jj,ll)%numc
               END DO

               ! Particles lost by cloud collection
               DO ll = 1, ncld
                  zminusterm = zminusterm + zccca(kk,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Particles lost by collection by rain
               DO ll = 1, nprc
                  zminusterm = zminusterm + zccpa(kk,ll)*precp(ii,jj,ll)%numc
               END DO

                  ! particles lost by ice collection
               DO ll = 1, nice
                  zminusterm = zminusterm + zccia(kk,ll)*ice(ii,jj,ll)%numc
               END DO

               ! particles lost by snow collection
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsa(kk,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Particle volume gained from smaller particles in 1/2a
               DO ll = in1a, index_2a-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcc(ll,kk)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Particle volume gained from smaller particles in 2b
               DO ll = in2b, kk-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcc(ll,kk)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               !-- Volume and number concentrations after coagulation update [fxm]
               aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+ptstep*zplusterm(1:nspec) *  &
                                             aero(ii,jj,kk)%numc ) / (1. + ptstep*zminusterm)

               aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + ptstep*zminusterm  + &
                                                            0.5*ptstep*zcc(kk,kk)*aero(ii,jj,kk)%numc)

            END DO


            ! Cloud droplets, regime a
            ! ------------------------------------------------
            DO cc = ica%cur, fca%cur
               IF (cloud(ii,jj,cc)%numc < nlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! corresponding index for regime b cloud droplets
               kk = MAX(cc-fca%cur+ncld,icb%cur) ! Regime a has more bins than b:
                                                      ! Set this at minimum to beginnign of b.

               ! Droplets lost by those with larger nucleus in regime a
               DO ll = cc+1, fca%cur
                  zminusterm = zminusterm + zcccc(cc,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Droplets lost by those with larger nucleus in regime b
               DO ll = kk+1, fcb%cur
                  zminusterm = zminusterm + zcccc(cc,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by rain drops
               DO ll = 1, nprc
                  zminusterm = zminusterm + zccpc(cc,ll)*precp(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by ice particles
               DO ll = 1, nice
                  zminusterm = zminusterm + zccic(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by snow particles
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsc(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Volume gained from cloud collection of aerosols
               DO ll = in1a, fn2b
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccca(ll,cc)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller droplets in a
               DO ll = ica%cur, cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcccc(ll,cc)*cloud(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller or equal droplets in b
               DO ll = icb%cur, kk
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcccc(ll,cc)*cloud(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               cloud(ii,jj,cc)%volc(1:nspec) = max(0.,( cloud(ii,jj,cc)%volc(1:nspec) +  &
                                                ptstep*zplusterm(1:nspec)*cloud(ii,jj,cc)%numc ) /  &
                                                (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
               cloud(ii,jj,cc)%numc = max(0., cloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                           0.5*ptstep*zcccc(cc,cc)*cloud(ii,jj,cc)%numc ) )

            END DO

            ! Cloud droplets, regime b
            ! -----------------------------------------
            DO cc = icb%cur, fcb%cur
               IF (cloud(ii,jj,cc)%numc < nlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! corresponding index for regime a cloud droplets
               kk = cc - ncld + fca%cur

               ! Droplets lost by those with larger nucleus in regime b
               DO ll = cc+1, fcb%cur
                  zminusterm = zminusterm + zcccc(cc,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Droplets lost by those with larger nucleus in regime a
               DO ll = kk+1, fca%cur
                  zminusterm = zminusterm + zcccc(cc,ll)*cloud(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by rain drops
               DO ll = 1, nprc
                  zminusterm = zminusterm + zccpc(cc,ll)*precp(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by ice
               DO ll = 1, nice
                  zminusterm = zminusterm + zccic(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Droplets lost by collection by snow particles
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsc(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Volume gained from cloud collection of aerosols
               DO ll = in1a, fn2b
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccca(ll,cc)*aero(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller droplets in b
               DO ll = icb%cur, cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcccc(ll,cc)*cloud(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller or equal droplets in a
               DO ll = ica%cur, kk
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zcccc(ll,cc)*cloud(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               cloud(ii,jj,cc)%volc(1:nspec) = max(0., ( cloud(ii,jj,cc)%volc(1:nspec) +  &
                                                      ptstep*zplusterm(1:nspec)*cloud(ii,jj,cc)%numc ) /     &
                                                      (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
               cloud(ii,jj,cc)%numc = max(0.,cloud(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                           0.5*ptstep*zcccc(cc,cc)*cloud(ii,jj,cc)%numc ) )

            END DO

            ! Ice particles, regime a
            ! ------------------------------------------------
            DO cc = iia%cur, fia%cur
               IF (ice(ii,jj,cc)%numc < prlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! corresponding index for regime b ice
               kk = MAX(cc-fia%cur+nice, iib%cur) ! Regime a has more bins than b:
                                                      ! Set this at minimum to beginning of b.

               ! Particles lost by those with larger nucleus in regime a
               DO ll = cc+1, fia%cur
                  zminusterm = zminusterm + zccii(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Particles lost by those with larger nucleus in regime b
               DO ll = kk+1, fib%cur
                  zminusterm = zminusterm + zccii(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Particles lost by collection by snow
               DO ll = 1, nsnw
                  zminusterm = zminusterm + zccsi(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Volume gained from aerosol collection
               DO ll = in1a,fn2b
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccia(ll,cc)*aero(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccia(ll,cc)*aero(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from cloud collection
               DO ll = 1,ncld
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccic(ll,cc)*cloud(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccic(ll,cc)*cloud(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from rain drops
               DO ll = 1,nprc
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccip(ll,cc)*precp(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccip(ll,cc)*precp(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from smaller ice particles in regime a
               DO ll = iia%cur,cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccii(ll,cc)*ice(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller or equal ice particles in regime b
               DO ll = iib%cur,kk
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccii(ll,cc)*ice(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               ice(ii,jj,cc)%volc(1:nspec) = max(0., ( ice(ii,jj,cc)%volc(1:nspec) +  &
                                                    ptstep*zplusterm(1:nspec)*ice(ii,jj,cc)%numc ) / &
                                                    (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
               ice(ii,jj,cc)%numc = max(0.,ice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                                                  0.5*ptstep*zccii(cc,cc)*ice(ii,jj,cc)%numc ) )

            END DO

            ! Ice particles, regime b
            ! -----------------------------------------
            DO cc = iib%cur, fib%cur
               IF (ice(ii,jj,cc)%numc < prlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! corresponding index for regime a
               kk = cc - nice + fia%cur

               ! Particles lost by those with larger nucleus in regime b
               DO ll = cc+1, fib%cur
                  zminusterm = zminusterm + zccii(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Particles lost by those with larger nucleus in regime a
               DO ll = kk+1, fia%cur
                  zminusterm = zminusterm + zccii(cc,ll)*ice(ii,jj,ll)%numc
               END DO

               ! Particles lost by collection by snow
               DO ll = 1,nsnw
                  zminusterm = zminusterm + zccsi(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Volume gained from aerosol collection
               DO ll = in1a,fn2b
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccia(ll,cc)*aero(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccia(ll,cc)*aero(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from cloud collection
               DO ll = 1,ncld
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccic(ll,cc)*cloud(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccic(ll,cc)*cloud(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from rain drops
               DO ll = 1,nprc
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccip(ll,cc)*precp(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccip(ll,cc)*precp(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhoic
               END DO

               ! Volume gained from smaller ice particles in b
               DO ll = iib%cur,cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccii(ll,cc)*ice(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Volume gained from smaller ice particles in a
               DO ll = iia%cur,kk-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccii(ll,cc)*ice(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               ice(ii,jj,cc)%volc(1:nspec) = max(0.,( ice(ii,jj,cc)%volc(1:nspec) +  &
                                                   ptstep*zplusterm(1:nspec)*ice(ii,jj,cc)%numc ) / &
                                                   (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
               ice(ii,jj,cc)%numc = max(0.,ice(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                         0.5*ptstep*zccii(cc,cc)*ice(ii,jj,cc)%numc ) )

            END DO

            ! Snow
            ! -----------------------------------
            DO cc = 1, nsnw
               IF (snow(ii,jj,cc)%numc < prlim) CYCLE

               zminusterm = 0.
               zplusterm(:) = 0.

               ! Drops lost by coagulation with larger snow
               DO ll = cc+1, nsnw
                  zminusterm = zminusterm + zccss(cc,ll)*snow(ii,jj,ll)%numc
               END DO

               ! Volume gained by collection of aerosols
               DO ll = in1a,fn2b
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccsa(ll,cc)*aero(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccsa(ll,cc)*aero(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhosn
               END DO

               ! Volume gained by collection of cloud droplets
               DO ll = 1,ncld
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccsc(ll,cc)*cloud(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccsc(ll,cc)*cloud(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhosn
               END DO

               ! Volume gained by collection of rain drops
               DO ll = 1,nprc
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccsp(ll,cc)*precp(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccsp(ll,cc)*precp(ii,jj,ll)%volc(iwa)*spec%rhowa/spec%rhosn
               END DO

               ! Volume gained by collection of ice particles
               DO ll = 1,nice
                  zplusterm(1:ndry) = zplusterm(1:ndry) + zccsi(ll,cc)*ice(ii,jj,ll)%volc(1:ndry)
                  zplusterm(iwa) = zplusterm(iwa) + zccsi(ll,cc)*ice(ii,jj,ll)%volc(iwa)*spec%rhoic/spec%rhosn
               END DO

               ! Volume gained from smaller snow
               DO ll = 1,cc-1
                  zplusterm(1:nspec) = zplusterm(1:nspec) + zccss(ll,cc)*snow(ii,jj,ll)%volc(1:nspec)
               END DO

               ! Update the hydrometeor volume concentrations
               snow(ii,jj,cc)%volc(1:nspec) = max(0.,( snow(ii,jj,cc)%volc(1:nspec) +  &
                                                    ptstep*zplusterm(1:nspec)*snow(ii,jj,cc)%numc ) / &
                                                    (1. + ptstep*zminusterm) )

               ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
               snow(ii,jj,cc)%numc = max(0.,snow(ii,jj,cc)%numc/( 1. + ptstep*zminusterm +  &
                                          0.5*ptstep*zccss(cc,cc)*snow(ii,jj,cc)%numc ) )

            END DO

         END DO ! kbdim
      END DO ! klev

      DEALLOCATE(zplusterm)

   END SUBROUTINE coagulation


   ! fxm: calculated for empty bins too
   ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
   !      and organic vapours (average values? 'real' values for each?)
   !********************************************************************
   !
   ! Subroutine CONDENSATION(kproma, kbdim,  klev,        &
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
                           level,allSALSA,                           &
                           pcsa,                               &
                           pcocnv,  pcocsv, pchno3, pcnh3,     &
                           prv,prs, prsi,ptemp,  ppres,  ptstep,    &
                           ppbl)

      USE mo_salsa_nucleation
      USE classSection
      USE mo_submctl, ONLY :      &
         fn2b,                      &
         nbins,ncld,nprc,                 &
         nice,nsnw,ntotal,           &
         lscndgas,                  & 
         lscndh2oae, lscndh2ocl, lscndh2oic, & ! Condensation to aerosols, clouds and ice particles
         aero,cloud,precp,ice,snow,    &
         spec,                         &
         nsnucl                     ! nucleation

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow
      INTEGER, INTENT(in) :: level

      REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep,                    & ! timestep [s]
         prs(kbdim,klev),           & ! Water vapor saturation mixing ratio
         prsi(kbdim,klev)              ! Saturation mixing ratio    [kg/m3]

      INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

      REAL, INTENT(INOUT) ::     &
         prv(kbdim,klev),          & ! Water vapor mixing ratio [kg/kg]
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
         pcnh3(kbdim,klev)           ! ammonia concentration [#/m3]

      TYPE(Section), INTENT(inout) :: &
           allSALSA(kbdim,klev,ntotal)

      REAL :: zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
              zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles
              zxocnv(kbdim,klev),         &
              zrh(kbdim,klev)

      zxocnv = 0.
      zxsa = 0.
      zj3n3 = 0.
      zrh(1:kbdim,:) = prv(1:kbdim,:)/prs(1:kbdim,:)

      !------------------------------------------------------------------------------

      ! Nucleation
      IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,  &
                                      ptemp,  zrh,    ppres,  &
                                      pcsa,   pcocnv, ptstep, zj3n3,  &
                                      zxsa,   zxocnv, ppbl            )

      ! Condensation of H2SO4 and organic vapors
      IF (lscndgas) CALL condgas(kproma,  kbdim,  klev,    krow,      &
                                 allSALSA,                            &
                                 pcsa, pcocnv, pcocsv,     &
                                 zxsa, ptemp,  ppres, ptstep )

      ! Condensation of water vapour
      IF (lscndh2ocl .OR. lscndh2oae .OR. lscndh2oic) &
         CALL gpparth2o(kproma, kbdim, klev, krow,  &
                        level,allSALSA,                       &
                        ptemp, ppres, prs, prsi, prv,     &
                        ptstep)

   END SUBROUTINE condensation

   !
   ! ----------------------------------------------------------------------------------------------------------
   !

   SUBROUTINE condgas(kproma,  kbdim,  klev,    krow,      &
                      allSALSA,                           &
                      pcsa,                               &
                      pcocnv,  pcocsv,      &
                      zxsa,ptemp,  ppres,  ptstep )

      USE classSection
      USE mo_submctl, ONLY :      &
         pi,                        &
         in1a, in2a,                & ! size bin indices
         fn2b,                      &
         ncld,                      &
         nprc,                      &
         nice,                      &
         nsnw,                      &
         ntotal,                    &
         nlim,                      &
         prlim,                     &
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [m3]
         spec,                      &
         aero,cloud,precp,ice,snow, &
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         massacc,                   & ! mass accomodation coefficients in each bin
         n3                           ! number of molecules in one 3 nm particle [1]

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

      REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep                       ! timestep [s]

      REAL, INTENT(INOUT) ::     &
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         zxsa(kbdim,klev)            ! ratio of sulphuric acid and organic vapor in 3nm particles

      TYPE(Section), INTENT(inout) :: &
           allSALSA(kbdim,klev,ntotal)


      !-- Local variables ----------------------
      INTEGER :: ii, jj    ! loop indices

      REAL ::                        &
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

      INTEGER :: ioc, iso4


      ioc = spec%getIndex("OC",notFoundValue = 0)
      iso4 = spec%getIndex("SO4",notFoundValue = 0)

      zj3n3 = 0.
      zxocnv = 0.

      zdvolsa = 0.
      zn_vs_c = 0.
      DO jj = 1, klev
         DO ii = 1, kbdim

            zdvoloc = 0.

            !-- 1) Properties of air and condensing gases --------------------
            zvisc  = (7.44523e-3*SQRT(ptemp(ii,jj)**3))/(5093.*(ptemp(ii,jj)+110.4))      ! viscosity of air [kg/(m s)]
            zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
            zmfp   = 3.*zdfvap*sqrt(pi*spec%msu/(8.*rg*ptemp(ii,jj)))                      ! mean free path [m]

            !-- 2) Transition regime correction factor for particles ---------
            !
            !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
            !  Topics in current aerosol research, Pergamon.
            !
            !  Size of condensing molecule considered only for
            !  nucleation mode (3 - 20 nm)
            !

            !-- particle Knudsen numbers
            zknud(in1a:in1a+1) = 2.*zmfp/(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)              ! Gases on aerosols
            zknud(in1a+2:fn2b) = 2.*zmfp/aero(ii,jj,in1a+2:fn2b)%dwet

            zknca(1:ncld) = 2.*zmfp/cloud(ii,jj,1:ncld)%dwet          ! Knudsen number for gases on cloud drplets

            zknpa(1:nprc) = 2.*zmfp/precp(ii,jj,1:nprc)%dwet          ! Knudsen number for gases on rain drops

            zknia(1:nice) = 2.*zmfp/ice(ii,jj,1:nice)%dwet          ! Knudsen number for gases on ice particles

            zknsa(1:nsnw) = 2.*zmfp/snow(ii,jj,1:nsnw)%dwet          ! Knudsen number for gases on snow flakes

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
                      (3.*pi*zvisc*aero(ii,jj,in1a:in1a+1)%dwet)

            !-- collision rate (gases on aerosols) [1/s]
            zcolrate = 0.
            zcolrate(in1a:in1a+1) = MERGE( 2.*pi*(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    &
                                          (zdfvap+zdfpart)*zbeta(in1a:in1a+1)*           &
                                          aero(ii,jj,in1a:in1a+1)%numc,                 &
                                          0.,                                            &
                                          aero(ii,jj,in1a:in1a+1)%numc > nlim        )

            zcolrate(in1a+2:fn2b) = MERGE( 2.*pi*aero(ii,jj,in1a+2:fn2b)%dwet*zdfvap*       &
                                           zbeta(in1a+2:fn2b)*aero(ii,jj,in1a+2:fn2b)%numc, &
                                           0.,                                               &
                                           aero(ii,jj,in1a+2:fn2b)%numc > nlim        )

            !-- gases on hydrometeors
            zcolrateca = 0.
            zcolrateca(1:ncld) = MERGE( 2.*pi*cloud(ii,jj,1:ncld)%dwet*zdfvap*         &
                                        zbetaca(1:ncld)*cloud(ii,jj,1:ncld)%numc,      &
                                        0.,                                             &
                                        cloud(ii,jj,1:ncld)%numc > nlim           )

            ! Gases on rain drops
            zcolratepa = 0.
            zcolratepa(1:nprc) = MERGE( 2.*pi*precp(ii,jj,1:nprc)%dwet*zdfvap*    &
                                        zbetapa(1:nprc)*precp(ii,jj,1:nprc)%numc, &
                                        0.,                                        &
                                        precp(ii,jj,1:nprc)%numc > prlim       )
            !-- gases on ice particles
            zcolrateia = 0.
            zcolrateia(1:nice) = MERGE( 2.*pi*ice(ii,jj,1:nice)%dwet*zdfvap*      &
                                        zbetaia(1:nice)*ice(ii,jj,1:nice)%numc,   &
                                        0.,                                        &
                                        ice(ii,jj,1:ncld)%numc > prlim         )

            ! Gases on snow flakes
            zcolratesa = 0.
            zcolratesa(1:nsnw) = MERGE( 2.*pi*snow(ii,jj,1:nsnw)%dwet*zdfvap*     &
                                        zbetasa(1:nsnw)*snow(ii,jj,1:nsnw)%numc,  &
                                        0.,                                        &
                                        snow(ii,jj,1:nsnw)%numc > prlim        )

            !-- 4) Condensation sink [1/s] -------------------------------------

            zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)+ sum(zcolrateia) + sum(zcolratesa)  ! total sink

            !-- 5) Changes in gas-phase concentrations and particle volume -----
            !
            !--- 5.1) Organic vapours ------------------------

            !---- 5.1.1) Non-volatile organic compound: condenses onto all bins
            IF(pcocnv(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. ioc > 0) THEN

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

               aero(ii,jj,in1a:fn2b)%volc(ioc) = aero(ii,jj,in1a:fn2b)%volc(ioc) + & !-- change of volume
                                                 zdvoloc                             !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc) +  &
                                              zcolrateca(1:ncld)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc) +  &
                                              zcolratepa(1:nprc)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc) +  &
                                            zcolrateia(1:nice)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on snow
               snow(ii,jj,1:nsnw)%volc(ioc) = snow(ii,jj,1:nsnw)%volc(ioc) +  &
                                             zcolratesa(1:nsnw)/zcs_ocnv*mvoc*zdvap2

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
               ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
               IF (zxocnv(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
               END IF

            END IF


            !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
            zcs_ocsv = sum(zcolrate(in2a:fn2b)) +  &       ! sink for semivolatile organics
                       sum(zcolrateca(1:ncld))  +  &       ! ... including condensation on cloud droplets
                       sum(zcolratepa(1:nprc))  +  &       ! and rain drops
                       sum(zcolrateia(1:nice))  +  &       ! and ice particles
                       sum(zcolratesa(1:nsnw))             ! and snow particles

            IF(pcocsv(ii,jj) > 1.e-10 .AND. zcs_ocsv > 1.e-30 .AND. ioc > 0) THEN


               zcvap_new3 = pcocsv(ii,jj)/(1.+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
               zdvap3 = pcocsv(ii,jj) - zcvap_new3                  ! change in gas concentration [#/m3]
               pcocsv(ii,jj) = zcvap_new3                           ! updating gas concentration [#/m3]

               zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &                     ! volume change of particles
                                    zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3   !  [m3(OC)/m3(air)]

               aero(ii,jj,in1a:fn2b)%volc(ioc) = &                   !-- change of volume due
                  aero(ii,jj,in1a:fn2b)%volc(ioc) + zdvoloc        !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc)  +  &
                                              zcolrateca(1:ncld)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc)  +  &
                                              zcolratepa(1:nprc)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc)  +  &
                                            zcolrateia(1:nice)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on snow particles
               snow(ii,jj,1:nsnw)%volc(ioc) = snow(ii,jj,1:nsnw)%volc(ioc)  +  &
                                             zcolratesa(1:nprc)/zcs_ocsv*mvoc*zdvap3

            END IF


            ! ---- 5.2) Sulphate -------------------------------------------
            IF(pcsa(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. iso4 > 0) THEN

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
               aero(ii,jj,in1a:fn2b)%volc(iso4) = aero(ii,jj,in1a:fn2b)%volc(iso4) + zdvolsa

               !-- Clouds
               cloud(ii,jj,1:ncld)%volc(iso4) = cloud(ii,jj,1:ncld)%volc(iso4)  +  &
                                              zcolrateca(1:ncld)/zcs_su*mvsu*zdvap1

               ! Rain drops
               precp(ii,jj,1:nprc)%volc(iso4) = precp(ii,jj,1:nprc)%volc(iso4)  +  &
                                              zcolratepa(1:nprc)/zcs_su*mvsu*zdvap1

               !-- Ice clouds
               ice(ii,jj,1:nice)%volc(iso4) = ice(ii,jj,1:nice)%volc(iso4)  +  &
                                            zcolrateia(1:nice)/zcs_su*mvsu*zdvap1

               ! Snow particles
               snow(ii,jj,1:nsnw)%volc(iso4) = snow(ii,jj,1:nsnw)%volc(iso4)  +  &
                                             zcolratesa(1:nsnw)/zcs_su*mvsu*zdvap1

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               IF (zxsa(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc +          &
                                           zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
               END IF

            END IF

         END DO ! kbdim

      END DO ! klev

   END SUBROUTINE condgas

   !
   ! ----------------------------------------------------------------------------------------------------------
   !

   SUBROUTINE gpparth2o(kproma, kbdim,  klev, krow,  &
                        level,allSALSA,               &
                        ptemp,  ppres,  prs,prsi, prv,    &
                        ptstep)
    
      USE classSection
      USE mo_submctl, ONLY : nbins, ncld, nprc,    &
                             nice, nsnw, ntotal,            &
                             spec,                           &
                             mair,                         &
                             aero,cloud,precp,ice,snow,   &
                             surfw0, surfi0, rg,           &
                             pi, pi6, prlim, nlim,      &
                             massacc, avog,  &
                             in1a, in2a,  &
                             fn2b,            &
                             lscndh2oae, lscndh2ocl, lscndh2oic, &
                             alv, als, CalcDimension
      USE mo_salsa_properties, ONLY : equilibration
      IMPLICIT NONE

      INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
      INTEGER, INTENT(in) :: level
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)
      TYPE(Section), INTENT(inout) :: allSALSA(kbdim,klev,ntotal)

      REAL, INTENT(inout) :: prv(kbdim,klev)

      REAL :: zkelvin(nbins), zkelvincd(ncld), zkelvinpd(nprc), &  ! Kelvin effects
              zkelvinid(nice), zkelvinsd(nsnw)                      ! Kelvin effects ice'n'snow
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
      REAL :: adt,ttot
      REAL :: dwet, cap
      REAL :: zrh(kbdim,klev)
      REAL :: dwice(nice), dwsnow(nsnw)

      REAL :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

      INTEGER :: nstr
      INTEGER :: ii,jj,cc

      INTEGER :: iwa,nspec

      zrh(:,:) = prv(:,:)/prs(:,:)
      
      iwa = spec%getIndex("H2O")
      nspec = spec%getNSpec()

      ! Calculate the condensation only for 2a/2b aerosol bins
      nstr = in2a

      ! Save the current aerosol water content
      zaelwc1(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      ! For 1a bins do the equilibrium calculation
      CALL equilibration(kproma,kbdim,klev,      &
                         zrh,ptemp,.FALSE. )

      ! If RH < 98 % OR dynamic condensation for aerosols switched off, do equilibrium for all bins
      IF (zrh(1,1) < 0.98 .OR. .NOT. lscndh2oae)  CALL equilibration(kproma,kbdim,klev,      &
                                                                     zrh,ptemp,.TRUE. )

      ! The new aerosol water content after equilibrium calculation
      zaelwc2(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )/(ppres(:,:)*mair/(rg*ptemp(:,:)))

      DO jj = 1, klev
         DO ii = 1, kbdim
			! Necessary?
            IF ( .NOT. ( &
                (ANY(cloud(ii,jj,:)%numc > nlim) .OR. ANY(precp(ii,jj,:)%numc > prlim) .AND. lscndh2ocl) .OR. &
                (ANY(ice(ii,jj,:)%numc > prlim) .OR. ANY(snow(ii,jj,:)%numc > prlim) .AND. lscndh2oic) .OR. &
                (ANY(aero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) &
                ) ) CYCLE

            rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

            zdfh2o = ( 5./(16.*avog*rhoair*1.e-3*(3.11e-8)**2) ) * &
               SQRT( rg*1e7*ptemp(ii,jj)*mair*1.e3*(spec%mwa+mair)*1.e3/( 2.*pi*spec%mwa*1.e3 ) )
            zdfh2o = zdfh2o*1.e-4

            zmfph2o = 3.*zdfh2o*sqrt(pi*spec%mwa/(8.*rg*ptemp(ii,jj)))
            zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air

            ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
            zkelvinpd = 1.; zkelvincd = 1.; zkelvin = 1.; zkelvinid = 1.; zkelvinsd = 1.

            zcwc = 0.; zcwint = 0.; zcwn = 0.
            zcwcae = 0.; zcwccd = 0.; zcwcpd = 0.; zcwcid = 0.; zcwcsd = 0.;
            zcwintae = 0.; zcwintcd = 0.; zcwintpd = 0.; zcwintid = 0.; zcwintsd = 0.
            zcwnae = 0.; zcwncd = 0.; zcwnpd = 0.; zcwnid = 0.; zcwnsd = 0.
            zwsatae = 0.; zwsatcd = 0.; zwsatpd = 0.; zwsatid = 0.; zwsatsd = 0.
            
            zmtpd(:) = 0.
            zcwsurfpd(:) = 0.
            zmtcd(:) = 0.
            zcwsurfcd(:) = 0.
            zmtid(:) = 0.
            zcwsurfid(:) = 0.
            zmtsd(:) = 0.
            zcwsurfsd(:) = 0.
            zmtae(:) = 0.
            zcwsurfae(:) = 0.

            ! Cloud droplets --------------------------------------------------------------------------------
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > nlim .AND. lscndh2ocl) THEN

                  ! Wet diameter
                  dwet = ( SUM(cloud(ii,jj,cc)%volc(1:nspec))/cloud(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Activity + Kelvin effect
                  zact = acth2o(cloud(ii,jj,cc))
                  zkelvincd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  zcwsurfcd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatcd(cc) = zact*zkelvincd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = cloud(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Rain drops --------------------------------------------------------------------------------
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > prlim .AND. lscndh2ocl) THEN
                  ! Wet diameter
                  dwet = ( SUM(precp(ii,jj,cc)%volc(1:nspec))/precp(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Activity + Kelvin effect
                  zact = acth2o(precp(ii,jj,cc))
                  zkelvinpd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentrations over flat surface
                  zcwsurfpd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatpd(cc) = zact*zkelvinpd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = precp(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Ice particles --------------------------------------------------------------------------------
            ! Dimension
            CALL CalcDimension(nice,ice(ii,jj,:),prlim,dwice,4)
            DO cc = 1, nice
               IF (ice(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                  dwet=dwice(cc)
                     
                  ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                  cap=0.5*dwet
                     
                  ! Activity + Kelvin effect - edit when needed
                  !   Can be calculated just like for sperical homogenous particle or just ignored,
                  !   because these are not known for solid, irregular and non-homogenous particles.
                  !   Ice may not be that far from a sphere, but most particles are large and at least
                  !   growing particles are covered by a layer of pure ice.
                  zact = 1.0 
                  zkelvinid(cc) = exp( 4.*surfi0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )
                  
                  ! Saturation mole concentration over flat surface
                  zcwsurfid(cc) = prsi(ii,jj)*rhoair/spec%mwa
                  
                  ! Equilibrium saturation ratio
                  zwsatid(cc) = zact*zkelvinid(cc)
                  
                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                       (3.)*(zknud+zknud**2))
                  
                  ! Mass transfer according to Jacobson
                  zhlp1 = ice(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*als*zwsatid(cc)*zcwsurfid(cc)/(zthcond*ptemp(ii,jj)) 
                  zhlp3 = ( (als*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.
                  
                  zmtid(cc) = zhlp1/( zhlp2*zhlp3 + 1. )
                  
               END IF
            END DO
            
            ! Snow particles --------------------------------------------------------------------------------
            ! Dimension
            CALL CalcDimension(nsnw,snow(ii,jj,:),prlim,dwsnow,5)
            DO cc = 1, nsnw
               IF (snow(ii,jj,cc)%numc > prlim .AND. lscndh2oic) THEN
                  dwet=dwsnow(cc)
                  
                  ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                  cap=0.5*dwet
                  
                  ! Activity + Kelvin effect
                  !   Can be calculated just like for sperical homogenous particle or just ignored,
                  !   because these are not known for solid, irregular and non-homogenous particles.
                  !   Especially snow is typically highly irregular (e.g. dendrite).
                  zact = 1.0 
                  zkelvinsd(cc) = exp( 4.*surfi0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )
                  
                  ! Saturation mole concentrations over flat surface
                  zcwsurfsd(cc) = prsi(ii,jj)*rhoair/spec%mwa
                  
                  ! Equilibrium saturation ratio
                  zwsatsd(cc) = zact*zkelvinsd(cc)
                  
                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                       (3.)*(zknud+zknud**2))
                  
                  ! Mass transfer according to Jacobson
                  zhlp1 = snow(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*als*zwsatsd(cc)*zcwsurfsd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (als*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.
                  
                  zmtsd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )
                  
               END IF
            END DO
            
            ! -- Aerosols: ------------------------------------------------------------------------------------
            DO cc = 1, nbins
               IF (aero(ii,jj,cc)%numc > nlim .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) THEN
                  ! Wet diameter
                  dwet = ( SUM(aero(ii,jj,cc)%volc(1:nspec))/aero(ii,jj,cc)%numc/pi6 )**(1./3.)

                  ! Water activity + Kelvin effect
                  zact = acth2o(aero(ii,jj,cc))
                  zkelvin(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  ! Limit the supersaturation to max 1.01 for the mass transfer
                  ! EXPERIMENTAL
                  zcwsurfae(cc) = MAX(prs(ii,jj),prv(ii,jj)/1.01)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatae(cc) = zact*zkelvin(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.*massacc(cc))*(zknud+zknud**2))

                  ! Mass transfer
                  zhlp1 = aero(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Current mole concentrations
            zcwc = prv(ii,jj)*rhoair/spec%mwa
            zcwcae(1:nbins) = aero(ii,jj,1:nbins)%volc(iwa)*spec%rhowa/spec%mwa
            zcwccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(iwa)*spec%rhowa/spec%mwa
            zcwcpd(1:nprc) = precp(ii,jj,1:nprc)%volc(iwa)*spec%rhowa/spec%mwa
            zcwcid(1:nice) = ice(ii,jj,1:nice)%volc(iwa)*spec%rhoic/spec%mwa
            zcwcsd(1:nsnw) = snow(ii,jj,1:nsnw)%volc(iwa)*spec%rhosn/spec%mwa


            zcwtot = zcwc + SUM(zcwcae) + &
                     SUM(zcwccd) + &
                     SUM(zcwcpd) + &
                     SUM(zcwcid) + &
                     SUM(zcwcsd)
            ttot = 0.

            zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd; zcwintid = zcwcid; zcwintsd = zcwcsd

            ! Substepping loop
            ! ---------------------------------
            zcwint = 0.
            DO WHILE (ttot < ptstep)

               adt = 2.e-2
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

               IF ( ANY(aero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 ) THEN
                  DO cc = nstr, nbins
                     zcwintae(cc) = zcwcae(cc) + min(max(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                                                     -0.02*zcwcae(cc)),0.05*zcwcae(cc))
                     zwsatae(cc) = acth2o(aero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                  END DO
               END IF
               IF ( ANY(cloud(ii,jj,:)%numc > nlim) ) THEN
                  DO cc = 1, ncld
                     zcwintcd(cc) = zcwccd(cc) + min(max(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd(cc)), &
                                                     -0.02*zcwccd(cc)),0.05*zcwccd(cc))
                     zwsatcd(cc) = acth2o(cloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                  END DO
               END IF
               IF ( ANY(precp(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nprc
                     zcwintpd(cc) = zcwcpd(cc) + min(max(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd(cc)), &
                                                     -0.02*zcwcpd(cc)),0.05*zcwcpd(cc))
                     zwsatpd(cc) = acth2o(precp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                  END DO
               END IF
               IF (ANY(ice(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nice
                     zcwintid(cc) = zcwcid(cc) + min(max(adt*zmtid(cc)*(zcwint - zwsatid(cc)*zcwsurfid(cc)), &
                          -0.02*zcwcid(cc)),0.05*zcwcid(cc))
                     zwsatid(cc) = zkelvinid(cc)
                  END DO
               END IF
               IF (ANY(snow(ii,jj,:)%numc > prlim) ) THEN
                  DO cc = 1, nsnw
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

               ! Update vapor concentration for consistency
               zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                        SUM(zcwintcd(1:ncld))  - &
                        SUM(zcwintpd(1:nprc))  - &
                        SUM(zcwintid(1:nice))  - &
                        SUM(zcwintsd(1:nsnw))

               ! Update "old" values for next cycle
               zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd; zcwcid =zcwintid; zcwcsd = zcwintsd;
               zcwc = zcwint

               ttot = ttot + adt

            END DO ! ADT

            zcwn = zcwint
            zcwnae = zcwintae
            zcwncd = zcwintcd
            zcwnpd = zcwintpd
            zcwnid = zcwintid
            zcwnsd = zcwintsd

            prv(ii,jj) = zcwn*spec%mwa/rhoair

            aero(ii,jj,1:nbins)%volc(iwa) = max(0.,zcwnae(1:nbins)*spec%mwa/spec%rhowa)
            cloud(ii,jj,1:ncld)%volc(iwa) = max(0.,zcwncd(1:ncld)*spec%mwa/spec%rhowa)
            precp(ii,jj,1:nprc)%volc(iwa) = max(0.,zcwnpd(1:nprc)*spec%mwa/spec%rhowa)
            ice(ii,jj,1:nice)%volc(iwa) = max(0.,zcwnid(1:nice)*spec%mwa/spec%rhoic)
            snow(ii,jj,1:nsnw)%volc(iwa) = max(0.,zcwnsd(1:nsnw)*spec%mwa/spec%rhosn)

         END DO !kproma

      END DO ! klev

   END SUBROUTINE gpparth2o
   !-------------------------------------------------------
   REAL FUNCTION acth2o(ppart,pcw)

      USE classSection
      USE mo_submctl, ONLY : eps, spec
      IMPLICIT NONE

      TYPE(Section), INTENT(in) :: ppart
      REAL, INTENT(in), OPTIONAL  :: pcw

      REAL :: zns, znw
      INTEGER :: ndry, iwa, nn, ss
      CHARACTER(len=3) :: snam
      
      ndry = spec%getNSpec(type="dry")
      iwa = spec%getIndex("H2O")
      
      ! This is only relevant for solution particles so use rholiq
      zns = 0.
      DO nn = 1,ndry  ! Leaves out water and non-soluble species (zero dissociation factor)
         zns = zns + spec%diss(nn)*ppart%volc(nn)*spec%rholiq(nn)/spec%MM(nn)
      END DO 

      IF (present(pcw)) THEN
         znw = pcw
      ELSE
         znw = ppart%volc(iwa)*spec%rholiq(iwa)/spec%MM(iwa)
      END IF

      ! Assume activity coefficient of 1 for water...
      acth2o = MAX(0.1,znw/max(eps,(znw+zns)))

   END FUNCTION acth2o

   ! ------------------------------------------------------------------

   FUNCTION satvaph2o(ptemp) RESULT(psat)
      !-----------------------------------------------------------------
      ! Saturation vapour pressure of water vapour
      ! This is a local function for the subroutine *cloud_condensation*
      !
      ! J. Tonttila, FMI, 03/2014
      !-----------------------------------------------------------------
    
      IMPLICIT NONE

      REAL, INTENT(in) :: ptemp

      REAL, PARAMETER ::       &
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
   ! Function coagc
   ! ***************
   !
   ! Calculation of coagulation coefficients.
   ! Extended version of the function originally
   ! found in mo_salsa_init. This is now placed
   ! here to avoid cyclic dependencies between
   ! MODULEs upon coupling with UCLALES.
   !
   ! J. Tonttila, FMI, 05/2014
   !
   !-------------------------------------------------
   REAL FUNCTION coagc(diam1,diam2,mass1,mass2,temp,pres,kernel,flag1,flag2)

      USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav, rd, terminal_vel

      IMPLICIT NONE

      !-- Input variables ----------
      REAL, INTENT(IN) :: &
         diam1,  &   ! diameters of colliding particles [m]
         diam2,  &   !
         mass1,  &   ! masses -"- [kg]
         mass2,  &
         temp,   &   ! ambient temperature [K]
         pres        ! ambient pressure [fxm]

      INTEGER, INTENT(in) :: kernel ! select the type of kernel: 1 - aerosol-aerosol coagulation (the original version)
                                        !                            2 - hydrometeor-aerosol or hydrometeor-hydrometeor coagulation
      INTEGER, INTENT(in) :: flag1,flag2 ! Parameter for identifying aerosol (1), cloud (2), precipitation (3), ice (4) and snow (5)

      !-- Output variables ---------

      !-- Local variables ----------
      REAL ::  &
         visc,     &   ! viscosity of air [kg/(m s)]
         vkin,     &   ! Kinematic viscosity of air [m2 s-1]
         zrhoa,    &   ! Density of air [kg m-3]
         mfp,      &   ! mean free path of air molecules [m]
         mdiam,    &   ! mean diameter of colliding particles [m]
         fmdist,   &   ! distance of flux matching [m]
         eddy_dis, &   ! Eddy dissipation time
         zecoll,   &   ! Collision efficiency for graviational collection
         zev,      &   !
         zea,      &
         zbrown,   &   ! Components for coagulation kernel; Brownian
         zbrconv,  &   !                                    Convective diffusion enhancement
         zgrav,    &   !                                    Gravitational collection
         ztshear,  &   ! turbulent shear
         zturbinert    ! turbulent inertia

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

      zbrown = 0.
      zbrconv = 0.
      zgrav = 0.
      zev = 0.
      coagc = 0.

      !-------------------------------------------------------------------------------

      !-- 0) Initializing particle and ambient air variables --------------------
      diam  = (/ diam1, diam2 /)       ! particle diameters [m]
      mpart = (/ mass1, mass2 /)       ! particle masses [kg]

      visc = (7.44523e-3*SQRT(temp**3))/(5093.*(temp+110.4)) ! viscosity of air [kg/(m s)]

      mfp = (1.656e-10*temp+1.828e-8)*pstand/pres ! mean free path of air [m]

      !-- 2) Slip correction factor for small particles -------------------------

      knud = 2.*mfp/diam                                    ! Knudsen number
      beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))! Cunningham correction factor
      ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)

      !-- 3) Particle properties ------------------------------------------------

      dfpart = beta*boltz*temp/(3.*pi*visc*diam)  ! diffusion coefficient [m2/s]
      mtvel  = sqrt((8.*boltz*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
      omega  = 8.*dfpart/(pi*mtvel)

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
            vkin  = visc/zrhoa   ! Kinematic viscosity of air [m2 s-1]

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
            IF (reyn(lrg) <= 1.) THEN
               zbrconv = 0.45*zbrown*( (reyn(lrg)*schm(sml))**(1./3.) )
            ELSE IF (reyn(lrg) > 1.) THEN
               zbrconv = 0.45*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1./3.) )
            END IF

            ! Turbulent Shear
            eddy_dis = 10.e-4 ! Values suggested by Sami - could be taken from the LES model?
            ztshear = SQRT(8.*pi*eddy_dis/(15.*vkin))*(0.5*(diam(1)+diam(2)))**3
            zturbinert = pi*eddy_dis**(0.75) /(grav* SQRT(SQRT( vkin )))  &
                         *(0.5*(diam(1)+diam(2)))**2* ABS(termv(1)-termv(2))

            ! gravitational collection
            zea = stok**2/( stok + 0.5 )**2
            IF (stok > 1.214) THEN
               zev = 0.75*LOG(2.*stok)/(stok - 1.214)
               zev = (1. + zev)**(-2)
            ELSE IF (stok <= 1.214) THEN
               zev = 0.
            END IF

            zecoll = (60.*zev + zea*reyn(lrg))/(60. + reyn(lrg))
            zgrav = zecoll * pi * mdiam**2
            zgrav = zgrav * ABS(termv(1)-termv(2))

            ! Total coagulation kernel
            coagc = zbrown  + zbrconv + SQRT(zgrav**2+ ztshear**2+ zturbinert**2)

      END SELECT

   END FUNCTION coagc


END MODULE mo_salsa_dynamics
