MODULE mo_salsa_coagulation_processes

  USE mo_submctl, ONLY : in1a,fn1a,in2a,fn2a,in2b,fn2b,   &
                         ica,fca,icb,fcb,iia,fia,iib,fib, &
                         nbins,ncld,nprc,nice,nsnw,       &
                         nlim,prlim,                      &
                         aero,cloud,precp,ice,snow,       &
                         spec,                            &
                         lscgaa, lscgcc, lscgpp, lscgii, lscgss,        &
                         lscgca, lscgpa, lscgpc, lscgia, lscgic,        &
                         lscgip, lscgsa, lscgsc, lscgsp, lscgsi
  USE classSection, ONLY : Section
  IMPLICIT NONE
 
  REAL, ALLOCATABLE :: zminusterm(:,:)
  REAL, ALLOCATABLE :: zplusterm(:,:,:)            ! shape == number of active aerosol species, kbdim,klev
 
  CONTAINS
    
  SUBROUTINE initialize_coagulation_processes(kbdim,klev,nn)
    INTEGER, INTENT(in) :: nn ! number of active aerosol species
    INTEGER, INTENT(in) :: kbdim,klev
    ALLOCATE(zplusterm(nn,kbdim,klev),   &
             zminusterm(kbdim,klev))

  END SUBROUTINE initialize_coagulation_processes

  !
  ! Aerosol coagulation
  ! ---------------------
  !
  SUBROUTINE coag_aero(kbdim,klev,tstep,    &
                       zccaa,zccca,zccpa,zccia,zccsa)

    ! Here, the aerosol particles are also allowed to collect cloud droplets, if they
    ! have a smaller dry ccn diameter. While this technically implies an artificial 
    ! deactivation (of a very small number) of cloud droplets, it preserves the mass also in the
    ! presense of giant ccn particles, which have been seen to cause problems with the more
    ! intuitive approach, where cloud droplets collect aerosol, but not the other way around.

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: zccaa(kbdim,klev,nbins,nbins),zccca(kbdim,klev,nbins,ncld),  &
                        zccpa(kbdim,klev,nbins,nprc),zccia(kbdim,klev,nbins,nice),   &
                        zccsa(kbdim,klev,nbins,nsnw)

    INTEGER :: index_2a,index_2b,index_cld
    INTEGER :: ii,jj,kk

    INTEGER :: nspec
    nspec = spec%getNSpec()
 
    ! Aerosols in regime 1a
    ! --------------------------------
    DO kk = in1a, fn1a
       IF (ALL(aero(:,:,kk)%numc < nlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.

       ! Particles lost by coagulation with larger aerosols
       IF (lscgaa) &
            CALL accumulateSink(kbdim,klev,zccaa,kk,kk+1,fn2b,aero,nbins,nbins,zminusterm)

       ! Particles lost by cloud collection (cloud bins a) by larger cloud bins
       IF (lscgca) THEN
          ! Cloud bins a
          index_cld = MAX(kk - ica%par + ica%cur, ica%cur)
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fca%cur,cloud,nbins,ncld,zminusterm)
          ! Cloud bins b
          index_cld = MAX(fca%cur + kk - icb%par + icb%cur, icb%cur)
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fcb%cur,cloud,nbins,ncld,zminusterm)
       END IF

       ! particles lost by rain collection
       IF (lscgpa) &
            CALL accumulateSink(kbdim,klev,zccpa,kk,1,nprc,precp,nbins,nprc,zminusterm)

       ! particles lost by ice collection
       IF (lscgia) &
            CALL accumulateSink(kbdim,klev,zccia,kk,1,nice,ice,nbins,nice,zminusterm)

       ! particles lost by snow collection
       IF (lscgsa) &
            CALL accumulateSink(kbdim,klev,zccsa,kk,1,nsnw,snow,nbins,nsnw,zminusterm)

       ! Particle volume gained from smaller particles in regime 1a
       IF (kk > in1a .AND. lscgaa) &
            CALL accumulateSource(kbdim,klev,nspec,zccaa,kk,in1a,kk-1,aero,nbins,nbins,zplusterm)

       ! Particles gained from cloud droplets with smaller dry ccn diameter (usually not relevant for the smallest particles but
       ! the possibility included for generity)
       IF (lscgca) THEN
          IF (kk > ica%par) THEN
             ! Get the corresponding cloud index minus 1
             index_cld = kk-ica%par
             CALL accumulateSourceReverse(kbdim,klev,nspec,zccca,kk,ica%cur,index_cld,cloud,nbins,ncld,zplusterm)
          END IF 
       END IF
      
       DO jj = 1,klev
          DO ii = 1,kbdim
             !-- Volume and number concentrations after coagulation update [fxm]
             aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+tstep*zplusterm(1:nspec,ii,jj) * &
                  aero(ii,jj,kk)%numc ) / (1. + tstep*zminusterm(ii,jj))
             
             aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + tstep*zminusterm(ii,jj)  + &
                  0.5*tstep*zccaa(ii,jj,kk,kk)*aero(ii,jj,kk)%numc)
          END DO
       END DO

    END DO
    
    IF (.FALSE.) THEN
    ! Aerosols in regime 2a
    ! ---------------------------------
    DO kk = in2a, fn2a
       IF (ALL(aero(:,:,kk)%numc < nlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! Find corresponding size bin in subregime 2b
       index_2b = kk - in2a + in2b
       
       ! Particles lost by larger particles in 2a
       IF (lscgaa) &
            CALL accumulateSink(kbdim,klev,zccaa,kk,kk+1,fn2a,aero,nbins,nbins,zminusterm)
       
       ! Particles lost by larger particles in 2b
       IF (lscgaa) &
            CALL accumulateSink(kbdim,klev,zccaa,kk,index_2b+1,fn2b,aero,nbins,nbins,zminusterm)
       
       ! Particles lost by cloud collection by larger or equal cloud bins
       IF (lscgca) THEN
          ! cloud bins a
          index_cld = MAX(kk - ica%par + ica%cur, ica%cur)
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fca%cur,cloud,nbins,ncld,zminusterm)
          ! cloud bins b
          index_cld = MAX(index_2b - icb%par + icb%cur, icb%cur)
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fcb%cur,cloud,nbins,ncld,zminusterm)
       END IF

       ! Particles lost by collection by rain
       IF (lscgpa) &
            CALL accumulateSink(kbdim,klev,zccpa,kk,1,nprc,precp,nbins,nprc,zminusterm)
       
       ! particles lost by ice collection
       IF (lscgia) &
            CALL accumulateSink(kbdim,klev,zccia,kk,1,nice,ice,nbins,nice,zminusterm)

       ! particles lost by snow collection
       IF (lscgsa) &
            CALL accumulateSink(kbdim,klev,zccsa,kk,1,nsnw,snow,nbins,nsnw,zminusterm)

       ! Particle volume gained from smaller particles in regimes 1, 2a
       IF (lscgaa) &
            CALL accumulateSource(kbdim,klev,nspec,zccaa,kk,in1a,kk-1,aero,nbins,nbins,zplusterm)
       
       ! Particle volume gained from smaller (and equal) particles in 2b
       IF (lscgaa) &
            CALL accumulateSource(kbdim,klev,nspec,zccaa,kk,in2b,index_2b,aero,nbins,nbins,zplusterm)

       ! Particles gained from cloud droplets with smaller dry ccn diameter
       IF (lscgca) THEN
          IF ( kk > ica%par ) THEN
             ! Get the corresponding cloud index minus 1
             index_cld = MAX(kk - ica%par + ica%cur - 1, ica%cur)
             CALL accumulateSourceReverse(kbdim,klev,nspec,zccca,kk,ica%cur,index_cld,cloud,nbins,ncld,zplusterm)
          END IF
          IF ( index_2b > icb%par ) THEN
             ! Get the corresponding cloud index minus 1
             index_cld = MAX(index_2b - icb%par + icb%cur - 1, icb%cur)
             CALL accumulateSourceReverse(kbdim,klev,nspec,zccca,kk,icb%cur,index_cld,cloud,nbins,ncld,zplusterm)
          END IF
       END IF
       
       DO jj = 1,klev
          DO ii = 1,kbdim
             !-- Volume and number concentrations after coagulation update [fxm]
             aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+tstep*zplusterm(1:nspec,ii,jj) *  &
                  aero(ii,jj,kk)%numc ) / (1. + tstep*zminusterm(ii,jj))
             
             aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + tstep*zminusterm(ii,jj)  + &
                  0.5*tstep*zccaa(ii,jj,kk,kk)*aero(ii,jj,kk)%numc)
          END DO
       END DO

    END DO
    
    ! Aerosols in regime 2b
    ! ---------------------------------
    DO kk = in2b, fn2b
       IF (ALL(aero(:,:,kk)%numc < nlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       !-- Find corresponding size bin in subregime 2a
       index_2a = kk - in2b + in2a
              
       ! Particles lost to larger particles in regimes 2b
       IF (lscgaa) &
            CALL accumulateSink(kbdim,klev,zccaa,kk,kk+1,fn2b,aero,nbins,nbins,zminusterm)

       ! Particles lost to larger and equal particles in 2a
       IF (lscgaa) &
            CALL accumulateSink(kbdim,klev,zccaa,kk,index_2a,fn2a,aero,nbins,nbins,zminusterm)
       
       ! Particles lost by cloud collection by larger cloud bins
       IF (lscgca) THEN
          ! Cloud bins a
          index_cld = MAX( index_2a - ica%par + ica%cur, ica%cur )
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fca%cur,cloud,nbins,ncld,zminusterm)
          ! Cloud bins b
          index_cld = MAX( kk - icb%par + icb%cur, icb%cur )
          CALL accumulateSink(kbdim,klev,zccca,kk,index_cld,fcb%cur,cloud,nbins,ncld,zminusterm)
       END IF

       ! Particles lost by collection by rain
       IF (lscgpa) &
            CALL accumulateSink(kbdim,klev,zccpa,kk,1,nprc,precp,nbins,nprc,zminusterm)
       
       ! particles lost by ice collection
       IF (lscgia) &
            CALL accumulateSink(kbdim,klev,zccia,kk,1,nice,ice,nbins,nice,zminusterm)
       
       ! particles lost by snow collection
       IF (lscgsa) &
            CALL accumulateSink(kbdim,klev,zccsa,kk,1,nsnw,snow,nbins,nsnw,zminusterm)
       
       ! Particle volume gained from smaller particles in 1/2a
       IF (lscgaa) &
            CALL accumulateSource(kbdim,klev,nspec,zccaa,kk,in1a,index_2a-1,aero,nbins,nbins,zplusterm)
       
       ! Particle volume gained from smaller particles in 2b
       IF (lscgaa) &
            CALL accumulateSource(kbdim,klev,nspec,zccaa,kk,in2b,kk-1,aero,nbins,nbins,zplusterm)

       ! Particle volume gained from smaller cloud bins
       IF (lscgca) THEN
          IF ( index_2a > ica%par ) THEN
             ! Get the corresponding cloud index minus 1
             index_cld = MAX( index_2a - ica%par + ica%cur - 1, ica%cur )
             CALL accumulateSourceReverse(kbdim,klev,nspec,zccca,kk,ica%cur,index_cld,cloud,nbins,ncld,zplusterm)
          END IF
          IF ( kk > icb%par ) THEN
             ! Get the corresponding cloud index minus 1
             index_cld = MAX( kk - icb%par + icb%cur - 1, icb%cur )
             CALL accumulateSourceReverse(kbdim,klev,nspec,zccca,kk,icb%cur,index_cld,cloud,nbins,ncld,zplusterm)
          END IF

       END IF

       DO jj = 1,klev
          DO ii = 1,kbdim
             !-- Volume and number concentrations after coagulation update [fxm]
             aero(ii,jj,kk)%volc(1:nspec) = ( aero(ii,jj,kk)%volc(1:nspec)+tstep*zplusterm(1:nspec,ii,jj) *  &
                  aero(ii,jj,kk)%numc ) / (1. + tstep*zminusterm(ii,jj))
       
             aero(ii,jj,kk)%numc = aero(ii,jj,kk)%numc/(1. + tstep*zminusterm(ii,jj)  + &
                  0.5*tstep*zccaa(ii,jj,kk,kk)*aero(ii,jj,kk)%numc)
          END DO
       END DO
    END DO
    END IF

  END SUBROUTINE coag_aero

  !
  ! Cloud droplet coagulation
  ! --------------------------
  !
  SUBROUTINE coag_cloud(kbdim,klev,tstep,    &
                        zcccc,zccca,zccpc,zccic,zccsc)

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: zcccc(kbdim,klev,ncld,ncld),zccca(kbdim,klev,nbins,ncld),  &
                        zccpc(kbdim,klev,ncld,nprc),zccic(kbdim,klev,ncld,nice),   &
                        zccsc(kbdim,klev,ncld,nsnw)

    INTEGER :: ii,jj,kk,index_a, index_b,index_aero
    INTEGER :: nspec
    nspec = spec%getNSpec()
    
    ! Cloud bins a

    DO kk = ica%cur, fca%cur
       IF (ALL(cloud(:,:,kk)%numc < nlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! corresponding index for regime b cloud droplets
       index_b = MAX(kk-fca%cur+ncld,icb%cur) ! If regime a has more bins than b:
                                              ! Set this at minimum to beginnign of b.

       ! Droplets lost to larger aerosol bins
       IF (lscgca) THEN
          ! Aerosol bins a, parallel aerosol bin plus 1
          index_aero = ica%par + kk
          CALL accumulateSinkReverse(kbdim,klev,zccca,kk,index_aero,fn2a,aero,ncld,nbins,zminusterm)
          ! Aerosol bins b, aerosol bin plus 1
          index_aero = icb%par + (index_b-fca%cur)
          CALL accumulateSinkReverse(kbdim,klev,zccca,kk,index_aero,fn2b,aero,ncld,nbins,zminusterm)
       END IF

       ! Droplets lost by those with larger nucleus in regime a
       IF (lscgcc) &
            CALL accumulateSink(kbdim,klev,zcccc,kk,kk+1,fca%cur,cloud,ncld,ncld,zminusterm)
       
       ! Droplets lost by those with larger nucleus in regime b
       IF (lscgcc) &
            CALL accumulateSink(kbdim,klev,zcccc,kk,index_b+1,fcb%cur,cloud,ncld,ncld,zminusterm)
       
       ! Droplets lost by collection by rain drops
       IF (lscgpc) &
            CALL accumulateSink(kbdim,klev,zccpc,kk,1,nprc,precp,ncld,nprc,zminusterm)
       
       ! Droplets lost by collection by ice particles
       IF (lscgic) &
            CALL accumulateSink(kbdim,klev,zccic,kk,1,nice,ice,ncld,nice,zminusterm)
       
       ! Droplets lost by collection by snow particles
       IF (lscgsc) &
            CALL accumulateSink(kbdim,klev,zccsc,kk,1,nsnw,snow,ncld,nsnw,zminusterm)
       
       ! Volume gained from cloud collection of smaller and equal aerosol bins
       IF (lscgca) THEN
          ! corresponding index aerosol a
          index_aero = ica%par + (kk - ica%cur)
          CALL accumulateSource(kbdim,klev,nspec,zccca,kk,in1a,index_aero,aero,ncld,nbins,zplusterm)
          ! corresponding index aerosol b
          index_aero = icb%par + (index_b - icb%cur)
          CALL accumulateSource(kbdim,klev,nspec,zccca,kk,in2b,index_aero,aero,ncld,nbins,zplusterm)
       END IF 
          
       ! Volume gained from smaller droplets in a
       IF (lscgcc) &
            CALL accumulateSource(kbdim,klev,nspec,zcccc,kk,ica%cur,kk-1,cloud,ncld,ncld,zplusterm)

       ! Volume gained from smaller or equal droplets in b
       IF (lscgcc) &
            CALL accumulateSource(kbdim,klev,nspec,zcccc,kk,icb%cur,index_b,cloud,ncld,ncld,zplusterm)

       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             cloud(ii,jj,kk)%volc(1:nspec) = max(0.,( cloud(ii,jj,kk)%volc(1:nspec) +  &
                  tstep*zplusterm(1:nspec,ii,jj)*cloud(ii,jj,kk)%numc ) /  &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
             cloud(ii,jj,kk)%numc = max(0., cloud(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zcccc(ii,jj,kk,kk)*cloud(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO

    ! Cloud bins b
    DO kk = icb%cur, fcb%cur
       IF (ALL(cloud(:,:,kk)%numc < nlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
               
       ! corresponding index for regime a cloud droplets
       index_a = kk - ncld + fca%cur

       ! Droplets lost to larger aerosol bins
       IF (lscgca) THEN
          ! Aerosol bins a, parallel aerosol bin plus 1
          index_aero = ica%par + index_a
          CALL accumulateSinkReverse(kbdim,klev,zccca,kk,index_aero,fn2a,aero,ncld,nbins,zminusterm)
          ! Aerosol bins b, aerosol bin plus 1
          index_aero = icb%par + (kk-fca%cur)
          CALL accumulateSinkReverse(kbdim,klev,zccca,kk,index_aero,fn2b,aero,ncld,nbins,zminusterm)
       END IF
       
       ! Droplets lost by those with larger nucleus in regime b
       IF (lscgcc) &
            CALL accumulateSink(kbdim,klev,zcccc,kk,kk+1,fcb%cur,cloud,ncld,ncld,zminusterm)       
       
       ! Droplets lost by those with larger nucleus in regime a
       IF (lscgcc) &
            CALL accumulateSink(kbdim,klev,zcccc,kk,index_a+1,fca%cur,cloud,ncld,ncld,zminusterm)
       
       ! Droplets lost by collection by rain drops
       IF (lscgpc) &
            CALL accumulateSink(kbdim,klev,zccpc,kk,1,nprc,precp,ncld,nprc,zminusterm)
       
       ! Droplets lost by collection by ice
       IF (lscgic) &
            CALL accumulateSink(kbdim,klev,zccic,kk,1,nice,ice,ncld,nice,zminusterm)
       
       ! Droplets lost by collection by snow particles
       IF (lscgsc) &
            CALL accumulateSink(kbdim,klev,zccsc,kk,1,nsnw,snow,ncld,nsnw,zminusterm)
       
       ! Volume gained from cloud collection of smaller and equal aerosol bins
       IF (lscgca) THEN
          ! corresponding index aerosol a
          index_aero = ica%par + (index_a - ica%cur)
          CALL accumulateSource(kbdim,klev,nspec,zccca,kk,in1a,index_aero,aero,ncld,nbins,zplusterm)
          ! corresponding index aerosol b
          index_aero = icb%par + (kk - icb%cur)
          CALL accumulateSource(kbdim,klev,nspec,zccca,kk,in2b,index_aero,aero,ncld,nbins,zplusterm)
       END IF 
       
       ! Volume gained from smaller droplets in b
       IF (lscgcc) &
            CALL accumulateSource(kbdim,klev,nspec,zcccc,kk,icb%cur,kk-1,cloud,ncld,ncld,zplusterm)
       
       ! Volume gained from smaller or equal droplets in a
       IF (lscgcc) &
            CALL accumulateSource(kbdim,klev,nspec,zcccc,kk,ica%cur,index_a,cloud,ncld,ncld,zplusterm)
       
       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             cloud(ii,jj,kk)%volc(1:nspec) = max(0., ( cloud(ii,jj,kk)%volc(1:nspec) +  &
                  tstep*zplusterm(1:nspec,ii,jj)*cloud(ii,jj,kk)%numc ) /     &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
             cloud(ii,jj,kk)%numc = max(0.,cloud(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zcccc(ii,jj,kk,kk)*cloud(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO
    
  END SUBROUTINE coag_cloud


  !
  ! Precipitation coagulation
  ! ----------------------------
  !
  SUBROUTINE coag_precp(kbdim,klev,tstep,        &
                        zccpp,zccpa,zccpc,zccip,zccsp)

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: zccpp(kbdim,klev,nprc,nprc),zccpa(kbdim,klev,nbins,nprc), &
                        zccpc(kbdim,klev,ncld,nprc),zccip(kbdim,klev,nprc,nice),  &
                        zccsp(kbdim,klev,nprc,nsnw)

    INTEGER :: ii,jj,kk
    INTEGER :: nspec
    nspec = spec%getNSpec()

    DO kk = 1, nprc
       IF (ALL(precp(:,:,kk)%numc < prlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! Drops lost by coagulation with larger drops
       IF (lscgpp) &
            CALL accumulateSink(kbdim,klev,zccpp,kk,kk+1,nprc,precp,nprc,nprc,zminusterm)       
       
       ! Drops lost by collection by snow drops
       IF (lscgsp) &
            CALL accumulateSink(kbdim,klev,zccsp,kk,1,nsnw,snow,nprc,nsnw,zminusterm)       
       
       ! Drops lost by collisions with ice
       IF (lscgip) &
            CALL accumulateSink(kbdim,klev,zccip,kk,1,nice,ice,nprc,nice,zminusterm)       
       
       ! Volume gained by collection of aerosols
       IF (lscgpa) &
            CALL accumulateSource(kbdim,klev,nspec,zccpa,kk,in1a,fn2b,aero,nprc,nbins,zplusterm)
       
       ! Volume gained by collection of cloud droplets
       IF (lscgpc) &
            CALL accumulateSource(kbdim,klev,nspec,zccpc,kk,1,ncld,cloud,nprc,ncld,zplusterm)
       
       ! Volume gained from smaller drops
       IF (lscgpp .AND. precp(1,1,kk)%numc > prlim) THEN
          CALL accumulateSource(kbdim,klev,nspec,zccpp,kk,1,kk-1,precp,nprc,nprc,zplusterm)
       END IF

       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             precp(ii,jj,kk)%volc(1:nspec) = max(0., ( precp(ii,jj,kk)%volc(1:nspec) +  &
                  tstep*zplusterm(1:nspec,ii,jj)*precp(ii,jj,kk)%numc ) / &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
             precp(ii,jj,kk)%numc = max(0.,precp(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zccpp(ii,jj,kk,kk)*precp(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO
    
  END SUBROUTINE coag_precp

  !
  ! Ice coagulation
  ! -----------------
  !
  SUBROUTINE coag_ice(kbdim,klev,tstep,     &
                      zccii,zccia,zccic,zccip,zccsi)

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: zccii(kbdim,klev,nice,nice),zccia(kbdim,klev,nbins,nice),  &
                        zccic(kbdim,klev,ncld,nice),zccip(kbdim,klev,nprc,nice),   &
                        zccsi(kbdim,klev,nice,nsnw)

    INTEGER :: ii,jj,kk,cc
    INTEGER :: nwet,ndry,iwa
    REAL :: rhowa,rhoic

    nwet = spec%getNSpec(type="wet")
    ndry = spec%getNSpec(type="dry")
    iwa = spec%getIndex("H2O")
    rhowa = spec%rhowa
    rhoic = spec%rhoic

    ! ice bins a
    DO kk = iia%cur, fia%cur
       IF (ALL(ice(:,:,kk)%numc < prlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! corresponding index for regime b ice
       cc = MAX(kk-fia%cur+nice, iib%cur) ! Regime a has more bins than b:
                                          ! Set this at minimum to beginning of b.
       
       ! Particles lost by those with larger nucleus in regime a
       IF (lscgii) &
            CALL accumulateSink(kbdim,klev,zccii,kk,kk+1,fia%cur,ice,nice,nice,zminusterm)       

       ! Particles lost by those with larger nucleus in regime b
       IF (lscgii) &
            CALL accumulateSink(kbdim,klev,zccii,kk,kk+1,fib%cur,ice,nice,nice,zminusterm)       
       
       ! Particles lost by collection by snow
       IF (lscgsi) &
            CALL accumulateSink(kbdim,klev,zccsi,kk,1,nsnw,snow,nice,nsnw,zminusterm)       
       
       ! Volume gained from aerosol collection
       IF (lscgia) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccia,kk,in1a,fn2b, &
                                             aero,rhoic,rhowa,nice,nbins,zplusterm)
       
       ! Volume gained from cloud collection
       IF (lscgic) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccic,kk,1,ncld,    &
                                             cloud,rhoic,rhowa,nice,ncld,zplusterm)
       
       ! Volume gained from rain drops
       IF (lscgip) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccip,kk,1,nprc,    & 
                                             precp,rhoic,rhowa,nice,nprc,zplusterm)
       
       ! Volume gained from smaller ice particles in regime a
       IF (lscgii) &
            CALL accumulateSource(kbdim,klev,nwet,zccii,kk,iia%cur,kk-1,ice,nice,nice,zplusterm)
       
       ! Volume gained from smaller or equal ice particles in regime b
       IF (lscgii) &
            CALL accumulateSource(kbdim,klev,nwet,zccii,kk,iib%cur,cc,ice,nice,nice,zplusterm)
       
       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             ice(ii,jj,kk)%volc(1:nwet) = max(0., ( ice(ii,jj,kk)%volc(1:nwet) +  &
                  tstep*zplusterm(1:nwet,ii,jj)*ice(ii,jj,kk)%numc ) / &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
             ice(ii,jj,kk)%numc = max(0.,ice(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zccii(ii,jj,kk,kk)*ice(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO

    DO kk = iib%cur, fib%cur
       IF (ALL(ice(:,:,kk)%numc < prlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! corresponding index for regime a
       cc = kk - nice + fia%cur
       
       ! Particles lost by those with larger nucleus in regime b
       IF (lscgii) &
            CALL accumulateSink(kbdim,klev,zccii,kk,kk+1,fib%cur,ice,nice,nice,zminusterm)       
       
       ! Particles lost by those with larger nucleus in regime a
       IF (lscgii) &
            CALL accumulateSink(kbdim,klev,zccii,kk,cc+1,fia%cur,ice,nice,nice,zminusterm)       
       
       ! Particles lost by collection by snow
       IF (lscgsi) &
            CALL accumulateSink(kbdim,klev,zccsi,kk,1,nsnw,snow,nice,nsnw,zminusterm)       
       
       ! Volume gained from aerosol collection
       IF (lscgia) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccia,kk,in1a,fn2b, & 
                                             aero,rhoic,rhowa,nice,nbins,zplusterm)
       
       ! Volume gained from cloud collection
       IF (lscgic) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccic,kk,1,ncld,    &
                                             cloud,rhoic,rhowa,nice,ncld,zplusterm)
       
       ! Volume gained from rain drops
       IF (lscgip) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccip,kk,1,nprc,    &
                                             precp,rhoic,rhowa,nice,nprc,zplusterm)
       
       ! Volume gained from smaller ice particles in b
       IF (lscgii) &
            CALL accumulateSource(kbdim,klev,nwet,zccii,kk,iib%cur,kk-1,ice,nice,nice,zplusterm)
       
       ! Volume gained from smaller ice particles in a
       IF (lscgii) &
            CALL accumulateSource(kbdim,klev,nwet,zccii,kk,iia%cur,cc-1,ice,nice,nice,zplusterm)

       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             ice(ii,jj,kk)%volc(1:nwet) = max(0.,( ice(ii,jj,kk)%volc(1:nwet) +  &
                  tstep*zplusterm(1:nwet,ii,jj)*ice(ii,jj,kk)%numc ) / &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with lrger bins and self)
             ice(ii,jj,kk)%numc = max(0.,ice(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zccii(ii,jj,kk,kk)*ice(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO
    
  END SUBROUTINE coag_ice

  SUBROUTINE coag_snow(kbdim,klev,tstep,          &
                       zccss,zccsa,zccsc,zccsp,zccsi)

    INTEGER, INTENT(in) :: kbdim,klev
    REAL, INTENT(in) :: tstep
    REAL, INTENT(in) :: zccss(kbdim,klev,nsnw,nsnw),zccsa(kbdim,klev,nbins,nsnw),  &
                        zccsc(kbdim,klev,ncld,nsnw),zccsp(kbdim,klev,nprc,nsnw),   &
                        zccsi(kbdim,klev,nice,nsnw)

    INTEGER :: ii,jj,kk,cc
    INTEGER :: nwet,ndry,iwa
    REAL :: rhosn,rhoic,rhowa

    nwet = spec%getNSpec(type="wet")
    ndry = spec%getNSpec(type="dry")
    iwa = spec%getIndex("H2O")
    rhosn = spec%rhosn
    rhoic = spec%rhoic
    rhowa = spec%rhowa

    DO kk = 1, nsnw
       IF (ALL(snow(:,:,kk)%numc < prlim)) CYCLE
       
       zminusterm(:,:) = 0.
       zplusterm(:,:,:) = 0.
       
       ! Drops lost by coagulation with larger snow
       IF (lscgss) &
            CALL accumulateSink(kbdim,klev,zccss,kk,kk+1,nsnw,snow,nsnw,nsnw,zminusterm)       
       
       ! Volume gained by collection of aerosols  
       IF (lscgsa) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccsa,kk,in1a,fn2b,  &
                                             aero,rhosn,rhowa,nsnw,nbins,zplusterm)
       
       ! Volume gained by collection of cloud droplets
       IF (lscgsc) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccsc,kk,1,ncld,     &
                                             cloud,rhosn,rhowa,nsnw,ncld,zplusterm)
       
       ! Volume gained by collection of rain drops
       IF (lscgsp) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccsp,kk,1,nprc,     &
                                             precp,rhosn,rhowa,nsnw,nprc,zplusterm)
       
       ! Volume gained by collection of ice particles
       IF (lscgsi) &
            CALL accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,zccsi,kk,1,nice,     &
                                             ice,rhosn,rhoic,nsnw,nice,zplusterm)
       
       ! Volume gained from smaller snow
       IF (lscgss) &
            CALL accumulateSource(kbdim,klev,nwet,zccss,kk,1,kk-1,snow,nsnw,nsnw,zplusterm)

       DO jj = 1,klev
          DO ii = 1,kbdim
             ! Update the hydrometeor volume concentrations
             snow(ii,jj,kk)%volc(1:nwet) = max(0.,( snow(ii,jj,kk)%volc(1:nwet) +  &
                  tstep*zplusterm(1:nwet,ii,jj)*snow(ii,jj,kk)%numc ) / &
                  (1. + tstep*zminusterm(ii,jj)) )
             
             ! Update the hydrometeor number concentration (Removal by coagulation with larger bins and self)
             snow(ii,jj,kk)%numc = max(0.,snow(ii,jj,kk)%numc/( 1. + tstep*zminusterm(ii,jj) +  &
                  0.5*tstep*zccss(ii,jj,kk,kk)*snow(ii,jj,kk)%numc ) )
          END DO
       END DO

    END DO
    

  END SUBROUTINE coag_snow

  SUBROUTINE accumulateSink(kbdim,klev,ck,itrgt,istr,iend,part,nbtrgt,nbpart,sink)

    INTEGER, INTENT(in) :: kbdim,klev,istr,iend,itrgt
    INTEGER, INTENT(in) :: nbtrgt,nbpart
    REAL, INTENT(in) :: ck(kbdim,klev,nbtrgt,nbpart) 
    TYPE(Section), INTENT(in) :: part(kbdim,klev,nbpart)
    REAL, INTENT(inout) :: sink(kbdim,klev)
    
    INTEGER :: ll,ii,jj
    
    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             sink(ii,jj) = sink(ii,jj) + ck(ii,jj,itrgt,ll)*part(ii,jj,ll)%numc
          END DO
       END DO
    END DO
    
  END SUBROUTINE accumulateSink

  ! ------------------------------------------------
  SUBROUTINE accumulateSinkReverse(kbdim,klev,ck,itrgt,istr,iend,part,nbtrgt,nbpart,sink)

    INTEGER, INTENT(in) :: kbdim,klev,istr,iend,itrgt
    INTEGER, INTENT(in) :: nbtrgt,nbpart
    REAL, INTENT(in) :: ck(kbdim,klev,nbpart,nbtrgt)
    TYPE(Section), INTENT(in) :: part(kbdim,klev,nbpart)
    REAL, INTENT(inout) :: sink(kbdim,klev)

    INTEGER :: ll,ii,jj
    
    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             sink(ii,jj) = sink(ii,jj) + ck(ii,jj,ll,itrgt)*part(ii,jj,ll)%numc
          END DO
       END DO
    END DO
      
  END SUBROUTINE accumulateSinkReverse

  ! ------------------------------------------------

  SUBROUTINE accumulateSource(kbdim,klev,nspec,ck,itrgt,istr,iend,part,nbtrgt,nbpart,source)

    INTEGER, INTENT(in) :: kbdim,klev,nspec,istr,iend,itrgt
    INTEGER, INTENT(in) :: nbtrgt,nbpart
    REAL, INTENT(in) :: ck(kbdim,klev,nbpart,nbtrgt)
    TYPE(Section), INTENT(in) :: part(kbdim,klev,nbpart)
    REAL, INTENT(inout) :: source(nspec,kbdim,klev)

    INTEGER :: ll,ii,jj

    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + ck(ii,jj,ll,itrgt)*part(ii,jj,ll)%volc(1:nspec)
          END DO
       END DO
    END DO

  END SUBROUTINE accumulateSource

  ! -------------------------------------------------


  !
  ! -----------------------------------------
  ! SUBROUTINE accumulateSourceGiantCCN
  !
  ! This is a special case subroutine for simulations with significant contribution from very large
  ! aerosol particles. It will substitute the regular mechanism for cloud collection of aerosol. Since
  ! the aerosol and cloud bins are set in identical size bins, for each collision, the larger of the
  ! two bins will determine in which cloud size bin the contribution of the collection mechanism will
  ! be placed (in the regular approach it is always the same cloud bin even if the collected aerosol would
  ! be significantly larger). 
  !
  ! The selection between this and the basic approach is given by the switch lscollectGCCN
  !
  SUBROUTINE accumulateSourceReverse(kbdim,klev,nspec,ck,itrgt,istr,iend,part,nbtrgt,nbpart,source)
    ! 
    ! ---------------------------------------------------------------------------------------------------
    ! This is meant for collection processes operating in a "reverse" direction in terms of the general
    ! collision process hierarchy. This is specifically designed for treating the collection processes
    ! between aerosol and cloud droplets, which use identical size bins based on dry particle size.
    !
    INTEGER, INTENT(in) :: kbdim,klev,nspec
    INTEGER, INTENT(in) :: istr,iend,itrgt
    INTEGER, INTENT(in) :: nbtrgt,nbpart  ! Number of bins in the target category and the coagulating category
    REAL, INTENT(in) :: ck(kbdim,klev,nbtrgt,nbpart)
    TYPE(Section), INTENT(in) :: part(kbdim,klev,nbpart)
    REAL, INTENT(inout) :: source(nspec,kbdim,klev)

    INTEGER :: iae_cld, ll, jj, ii, cc
    REAL :: contr(nspec)

    source(:,:,:) = 0.
    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + ck(ii,jj,itrgt,ll)*part(ii,jj,ll)%volc(1:nspec)
          END DO
       END DO
    END DO

  END SUBROUTINE accumulateSourceReverse

  ! -------------------------------------------------

  SUBROUTINE accumulateSourcePhaseChange(kbdim,klev,ndry,iwa,ck,itrgt,istr,iend,  &
                                         part,rhotrgt,rhosource,nbtrgt,nbpart,source)

    INTEGER,INTENT(in) :: kbdim,klev,ndry,iwa,istr,iend,itrgt
    INTEGER, INTENT(in) :: nbtrgt,nbpart
    REAL, INTENT(in) :: ck(kbdim,klev,nbpart,nbtrgt)
    TYPE(Section), INTENT(in) :: part(kbdim,klev,nbpart)
    REAL, INTENT(in) :: rhotrgt,rhosource ! densities of the target and source materials (for phase changes)
    REAL, INTENT(inout) :: source(ndry+1,kbdim,klev)

    INTEGER :: ll,ii,jj

    DO ll = istr,iend
       DO jj = 1,klev
          DO ii = 1,kbdim
             source(1:ndry,ii,jj) = source(1:ndry,ii,jj) + ck(ii,jj,ll,itrgt)*part(ii,jj,ll)%volc(1:ndry)
             source(iwa,ii,jj) = source(iwa,ii,jj) + ck(ii,jj,ll,itrgt)*part(ii,jj,ll)%volc(iwa)*rhosource/rhotrgt
          END DO
       END DO
    END DO

  END SUBROUTINE accumulateSourcePhaseChange

END MODULE mo_salsa_coagulation_processes
