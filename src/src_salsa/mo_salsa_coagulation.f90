MODULE mo_salsa_coagulation
  USE classSection, ONLY : Section
  USE mo_submctl, ONLY : spec,nbins,ncld,nprc,nsnw,nice
  IMPLICIT NONE
  
  ! ---------------------------------------------------
  ! Contains subroutines and functions used for 
  ! coagulation calculations in mo_salsa_dynamics.f90
  !
  
  ! Arrays for coagulation coefficients. These are imported to mo_salsa_dynamics, and then given as 
  ! actual arguments to the build_coagulation_kernels subroutine (this avoids cyclic dependencies while 
  ! keeping the code relatively clean).
  REAL, ALLOCATABLE ::          &
       zccaa(:,:,:,:),          & ! updated coagulation coefficients [m3/s]
       zcccc(:,:,:,:),          & ! - '' - for collision-coalescence between cloud droplets [m3/s]
       zccca(:,:,:,:),          & ! - '' - for cloud collection of aerosols [m3/s]
       zccpc(:,:,:,:),          & ! - '' - for collection of cloud droplets by precip [m3/s]
       zccpa(:,:,:,:),          & ! - '' - for collection of aerosols by precip
       zccpp(:,:,:,:),          & ! - '' - for collision-coalescence between precip particles 
       zccia(:,:,:,:),          & ! - '' - for collection of aerosols by ice 
       zccic(:,:,:,:),          & ! - '' - for collection of cloud particles droplets by ice 
       zccii(:,:,:,:),          & ! - '' - for aggregation between ice 
       zccip(:,:,:,:),          & ! - '' - for collection of precip by ice
       zccsa(:,:,:,:),          & ! - '' - for collection of aerosols by snow 
       zccsc(:,:,:,:),          & ! - '' - for collection of cloud droples by snow 
       zccsi(:,:,:,:),          & ! - '' - for collection of ice by snow 
       zccsp(:,:,:,:),          & ! - '' - for collection of precip by snow 
       zccss(:,:,:,:),          & ! - '' - for aggregation between snow
       
    CONTAINS


  ! 
  SUBROUTINE coag_aero

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
    

  END SUBROUTINE coag_aero

  SUBROUTINE coag_cloud
    
  END SUBROUTINE coag_cloud

  SUBROUTINE coag_precp

  END SUBROUTINE coag_precp

  SUBROUTINE coag_ice

  END SUBROUTINE coag_ice

  SUBROUTINE coag_snow

  END SUBROUTINE coag_snow


  !
  SUBROUTINE accumulateSink(kbdim,klev,ck,itrgt,istr,iend,part,sink)
    INTEGER, INTENT(in) :: kbdim,klev,istr,iend,itrgt
    REAL, INTENT(in) :: ck(:,:,:,:) ! Check that the assumed size array works!
    TYPE(Section), INTENT(in) :: part(:,:,:)
    REAL, INTENT(inout) :: sink(:,:)
    
    INTEGER :: ll,ii,jj

    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             sink(ii,jj) = sink(ii,jj) + ck(ii,jj,itrgt,ll)*part(ii,jj,ll)%numc
          END DO
       END DO
    END DO
    
  END SUBROUTINE accumulateSink

  SUBROUTINE acccumulateSource(kbdim,klev,nspec,ck,itrgt,istr,iend,part,source)
    INTEGER, INTENT(in) :: kbdim,klev,nspec,istr,iend,itrgt
    REAL, INTENT(in) :: ck(:,:,:,:)
    TYPE(Section), INTENT(in) :: part(:,:,:)
    REAL, INTENT(inout) :: source(:,:,:)

    INTEGER :: ll,ii,jj

    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + ck(ii,jj,itrgt,ll)*part(ii,jj,ll)%volc(1:nspec)
          END DO
       END DO
    END DO

  END SUBROUTINE accumulateSource



  !
  ! ------------------------------------------------------------------
  ! Allocate the coagulation kernel arrays and initialize as zeros
  !
  SUBROUTINE initialize_coagulation_kernels(kbdim,klev)
    INTEGER, INTENT(in) :: kbdim,klev


    ! For bin dimensions, the first one is always the collected (or "smaller") and the second is the collector (or the "larger")
    ALLOCATE (                                                                                  &
         zcc(kbdim,klev,nbins,nbins),zcccc(kbdim,klev,ncld,ncld),zccpp(kbdim,klev,nprc,nprc),   &
         zccii(kbdim,klev,nice,nice),zcss(kbdim,klev,nsnw,nsnw),                                &     
         zccca(kbdim,klev,nbins,ncld),                                                          &
         zccpa(kbdim,klev,nbins,nprc), zccpc(kbdim,klev,ncld,nprc),                             &
         zccia(kbdim,klev,nbins,nice), zccic(kbdim,klev,ncld,nice), zccip(kbdim,klev,nprc,nice),&
         zccsa(kbdim,klev,nbins,nsnw), zccsc(kbdim,klev,ncld,nsnw), zccsp(kbdim,klev,nprc,nsnw),&
         zccsi(kbdim,klev,nice,nsnw)                                                            &
             )
          
    
  END SUBROUTINE initialize_coagulation_kernels

  !
  ! ----------------------------------------------------------
  ! Update the coagulation kernels for all active particles
  !
  SUBROUTINE update_coagulation_kernels(kbdim,klev)
    USE mo_submctl, ONLY : nlim, prlim, aero, cloud, precp, ice, snow

    LOGICAL :: any_aero, any_cloud, any_precp, any_ice, any_snow

    any_aero = ANY(aero(:,:,:)%numc > nlim)
    any_cloud = ANY(cloud(:,:,:)%numc > nlim)
    any_precp = ANY(precp(:,:,:)%numc > prlim)
    any_ice = ANY(ice(:,:,:)%numc > prlim)
    any_snow = ANY(snow(:,:,:)%numc > prlim)
    
    IF (lscgaa .AND. any_aero) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nbins,nlim,aero,            &
                                        nbins,nlim,aero,            &
                                        zccaa                       )
    
    ! Collision-coalescence between cloud droplets
    IF (lscgcc .AND. any_cloud) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        ncld,nlim,cloud,            &
                                        ncld,nlim,cloud,            &
                                        zcccc                       ) 


    ! Self-collection of rain drops
    IF (lscgpp .AND. any_precp) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nprc,prlim,precp,           &
                                        nprc,prlim,precp,           &
                                        zccpp                       ) 
         
    ! Self-collection of snow particles
    IF (lscgss .AND. any_snow) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nsnw,prlim,snow,            &
                                        nsnw,prlim,snow,            &
                                        zccss                       ) 
    
    !  collisions between ice particles
    IF (lscgii .AND. any_ice) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nice,prlim,ice,             &
                                        nice,prlim,ice,             &
                                        zccii                       ) 
    
    ! Cloud collection of aerosols
    IF (lscgca .AND. any_cloud .AND. any_aero) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nbins,nlim,aero,            &
                                        ncld,nlim,cloud,            &
                                        zccca                       ) 

    ! Collection of aerosols by rain
    IF (lscgpa .AND. any_aero .AND. any_precp) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nbins,nlim,aero,            &
                                        nprc,prlim,precp,           &
                                        zccpa                       ) 
    
    ! Collection of cloud droplets by rain
    IF (lscgpc .AND. any_cloud .AND. any_precp) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        ncld,nlim,cloud,            &
                                        nprc,prlim,precp,           &
                                        zccpc                       ) 
        
    !  collection of aerosols by ice
    IF (lscgia .AND. any_aero .AND. any_ice) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nbins,nlim,aero,            &
                                        nice,prlim,ice,             &
                                        zccia                       ) 
    
    !  collection of cloud particles droplets by ice
    IF (lscgic .AND. any_cloud .AND. any_ice) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        ncld,nlim,cloud,            &
                                        nice,prlim,ice,             &
                                        zccic                       ) 
    
    !  collection of precip by ice-collision
    IF (lscgip .AND. any_precp .AND. any_ice) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nprc,prlim,precp,           &
                                        nice,prlim,ice,             &
                                        zccip                       ) 
    
    ! Collection of aerosols by snow
    IF (lscgsa .AND. any_aero .AND. any_snow) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nbins,nlim,aero,            &
                                        nsnw,prlim,snow,            &
                                        zccsa                       ) 
    
    ! collection of precip by snow
    IF (lscgsp .AND. any_precp .AND. any_snow) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nprc,prlim,precp,           &
                                        nsnw,prlim,snow,            &
                                        zccsp                       )     

    ! collection of cloud droples by snow
    IF (lscgsc .AND. any_cloud .AND. any_snow) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        ncld,nlim,cloud,            &
                                        nsnw,prlim,snow,            &
                                        zccsc                       ) 
    
    ! collection of ice by snow
    IF (lscgsi .AND. any_ice .AND. any_snow ) &
         CALL build_coagulation_kernels(kbdim,klev,ptemp,ppres,     &
                                        nice,prlim,ice,             &
                                        nsnw,prlim,snow,            &
                                        zccsi                       ) 
       
  END SUBROUTINE update_coagulation_kernels
       

  !
  ! ------------------------------------------------------------------
  ! Calculate coagulation kernels for given particles
  !
  SUBROUTINE build_coagulation_kernels(kbdim,klev,temp, press,    &    
                                       nb1,nlim1,part1,           &
                                       nb2,nlim2,part2,           &
                                       ck                        )
    
    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: temp(kbdim,klev), press(kbdim,klev)
    INTEGER, INTENT(in) :: nb1,nb2
    REAL, INTENT(in) :: nlim1,nlim2
    TYPE(Section), INTENT(in) :: part1(kbdim,klev,nb1), part2(kbdim,klev,nb2)
    REAL, INTENT(out) :: ck(kbdim,klev,nb1,nb2)
    
    REAL :: mass1(nb1),mass2(nb2),diam1(nb1),diam2(nb2)
    INTEGER :: phase1,phase2
    
    INTEGER :: mm,nn,ii,jj,nspec
    INTEGER :: bb
    
    nspec = spec%getNSpec()
    ck = 0.
    
    DO jj = 1,klev
       DO ii = 1,kbdim
          
          !-- diameter [m] and mass [kg]
          phase1 = part1(ii,jj,1)%phase
          CALL CalcDiamSALSA(nb1,part1(ii,jj,:),nlim1,diam1,phase1)
          mass1 = 0.
          DO bb = 1,nb1
             IF (part1(ii,jj,bb)%numc < nlim1) CYCLE
             mass1(bb) = SUM((part1(ii,jj,bb)%volc(1:nspec)/part1(ii,jj,bb)%numc)*spec%rhos(phase1,1:nspec))
          END DO
          
          phase2 = part2(ii,jj,1)%phase
          CALL CalcDiamSALSA(nb2,part2(ii,jj,:),nlim2,diam2,phase2)
          mass2 = 0.
          DO bb = 1,nb2
             IF (part2(ii,jj,bb)%numc < nlim2) CYCLE
             mass2(bb) = SUM((part2(ii,jj,bb)%volc(1:nspec)/part2(ii,jj,bb)%numc)*spec%rhos(phase2,1:nspec))
          END DO
          
          DO mm = 1, nb1         ! smaller colliding particle
             IF (part1(ii,jj,mm)%numc < nlim) CYCLE
             DO nn = 1, nb2            ! larger colliding particle
                IF (part2(ii,jj,nn)%numc < nlim) CYCLE
                ck(ii,jj,mm,nn) = coagc(diam1(mm),diam2(nn),mass1(mm),mass2(nn),temp(ii,jj),pres(ii,jj),2,phase1,phase2)
             END DO
          END DO
          
       END DO ! ii
       
    END DO ! jj
    
    
  END SUBROUTINE build_coagulation_kernels
  
  
  
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
    
    USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav, rd 
    USE mo_particle_external_properties, ONLY : terminal_vel
    
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
       
       termv(1) = terminal_vel(diam(1),zrhop(1),zrhoa,visc,beta(1),flag1)
       termv(2) = terminal_vel(diam(2),zrhop(2),zrhoa,visc,beta(2),flag2)
       
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
  
 
END MODULE mo_salsa_coagulation
