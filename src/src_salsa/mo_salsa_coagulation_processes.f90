MODULE mo_salsa_coagulation_processes
  USE mo_salsa_types, ONLY : aero, cloud, precp, ice, snow
  USE mo_submctl, ONLY : in1a, fn1a, in2a, fn2a, in2b, fn2b, nbins,  &
                         ica, fca, icb, fcb, ncld,                  &
                         nprc,                                      &
                         iia, fia, iib, fib, nice,                    &
                         nsnw,                              &
                         precpbins,                                 &
                         lscgaa, lscgcc, lscgpp, lscgii, lscgss,  &
                         lscgca, lscgpa, lscgia, lscgsa,          &
                         lscgpc, lscgic, lscgsc,                  &
                         lscgip, lscgsp,                          &
                         lscgsi, &
                         lsauto, &
                         spec
  USE classSection, ONLY : Section
  USE mo_particle_external_properties, ONLY : calcDiamSALSA
  IMPLICIT NONE

  CONTAINS

    !
    ! Aerosol coagulation
    ! -----------------------
    SUBROUTINE coag_aero(kbdim,klev,nspec,ptstep,zccaa,zccca,zccpa,zccia,zccsa)
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccaa(kbdim,klev,nbins,nbins), zccca(kbdim,klev,nbins,ncld),    &
                          zccpa(kbdim,klev,nbins,nprc), zccia(kbdim,klev,nbins,nice),     &
                          zccsa(kbdim,klev,nbins,nsnw)
      INTEGER :: kk
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      INTEGER :: index_cld_a, index_cld_b
      INTEGER :: index_a, index_b

      !
      ! Bin regime a 
      ! -----------------
      !
      DO kk = in1a, fn1a
         IF ( ALL(aero(:,:,kk)%numc < aero(:,:,kk)%nlim) ) CYCLE

         index_cld_a = ica%cur + MAX(kk - ica%par, 0) ! Corresponding index in the cloud bins (or the first cloud bin) in regime a
         index_cld_b = icb%cur + (index_cld_a - ica%cur) ! ... and the same for regime b

         zplusterm(:,:,:) = 0.
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.

         ! Collisions with larger aerosols and self collection
         IF (lscgaa) THEN
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2b,zccaa,aero,zminusterm)
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,multp=0.5)
         end if

         ! Collection by larger and equal cloud droplet bins (regime a and b)
         IF (lscgca) THEN
            CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm)
            CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm)
         end if

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia%cur,fib%cur,zccia,ice,zminusterm)
              
         ! Collection by snow
         IF (lscgsa) &
              CALL accumulateSink(kbdim,klev,nbins,nsnw,kk,1,nsnw,zccsa,snow,zminusterm)

         ! Particle volume gained from smaller particles in regime 1a
         IF (lscgaa .AND. kk > in1a) &
              CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,kk-1,zccaa,aero,zplusterm)

         ! Particle volume gained from smaller cloud droplet bins (does not happen with the current configuration for 1a,
         ! but the possiblity is included for generity)
         IF (lscgca) THEN
            index_cld_a = ica%cur + (kk-ica%par) - 1
            index_cld_b = icb%cur + (index_cld_a - ica%cur)
            IF ( index_cld_a >= ica%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,ica%cur,index_cld_a,zccca,cloud,zplusterm)
            IF ( index_cld_b >= icb%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,icb%cur,index_cld_b,zccca,cloud,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update [fxm]
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,aero,ptstep,zplusterm,zminusterm,zminus_self)

      END DO !kk

      !
      ! Bin regime 2a
      ! ---------------------
      !
      DO kk = in2a,fn2a
         IF ( ALL(aero(:,:,kk)%numc < aero(:,:,kk)%nlim) ) CYCLE

         ! Corresponding index in bin regime 2b
         index_b = in2b + (kk - in2a)
         ! Corresponding cloud bins
         index_cld_a = ica%cur + MAX(kk - ica%par, 0)
         index_cld_b = icb%cur + MAX(kk - ica%par, 0)

         zplusterm(:,:,:) = 0.
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.

         IF (lscgaa) THEN
            ! Collision with larger aerosol in regime 2a and larger or equal aerosol in 2b
            IF ( kk < fn2a) THEN
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2a,zccaa,aero,zminusterm)
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,index_b,fn2b,zccaa,aero,zminusterm)
            END IF

            ! Self collection
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,multp=0.5)
         END IF

         ! Collection by larger or equal cloud droplet bins
         IF (lscgca) THEN
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm)
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm)
         END IF

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia%cur,fib%cur,zccia,ice,zminusterm)
              
         ! Collection by snow
         IF (lscgsa) &
              CALL accumulateSink(kbdim,klev,nbins,nsnw,kk,1,nsnw,zccsa,snow,zminusterm)

         ! Particle volume gained from smaller particles in regime a and b
         IF (lscgaa) THEN
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,kk-1,zccaa,aero,zplusterm)
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in2b,index_b-1,zccaa,aero,zplusterm)
         END IF

         ! Volume gained from smaller cloud droplet bins
         IF (lscgca) THEN
            index_cld_a = ica%cur + (kk-ica%par) - 1
            index_cld_b = icb%cur + (index_cld_a - ica%cur)
            IF ( index_cld_a >= ica%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,ica%cur,index_cld_a,zccca,cloud,zplusterm)
            IF ( index_cld_b >= icb%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,icb%cur,index_cld_b,zccca,cloud,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update [fxm]
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,aero,ptstep,zplusterm,zminusterm,zminus_self)

      END DO

      !
      ! Bin regime 2b
      ! ---------------------
      !
      DO kk = in2b,fn2b
         IF ( ALL(aero(:,:,kk)%numc < aero(:,:,kk)%nlim) ) CYCLE

         ! Corresponding index in bin regime 2b
         index_a = in2a + (kk - in2b)
         ! Corresponding cloud bins
         index_cld_a = ica%cur + MAX(kk - icb%par, 0)
         index_cld_b = icb%cur + MAX(kk - icb%par, 0)

         zplusterm(:,:,:) = 0.
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.

         IF (lscgaa) THEN
            ! Collision with larger aerosol in regime 2b and larger or equal in 2a
            IF ( kk < fn2b) THEN
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2b,zccaa,aero,zminusterm)
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,index_a,fn2a,zccaa,aero,zminusterm)
            END IF

            ! Self collection
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,multp=0.5)
         END IF

         ! Collection by larger or equal cloud droplet bins
         IF (lscgca) THEN
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm)
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm)
         END IF

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia%cur,fib%cur,zccia,ice,zminusterm)
              
         ! Collection by snow
         IF (lscgsa) &
              CALL accumulateSink(kbdim,klev,nbins,nsnw,kk,1,nsnw,zccsa,snow,zminusterm)

         ! Particle volume gained from smaller particles in regime a and b
         IF (lscgaa) THEN
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in2b,kk-1,zccaa,aero,zplusterm)
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,index_a-1,zccaa,aero,zplusterm)
         END IF

         ! Volume gained from smaller cloud droplet bins (reverse collection)
         IF (lscgca) THEN
            index_cld_b = icb%cur + (kk-icb%par) - 1
            index_cld_a = ica%cur + (index_cld_b-icb%cur)
            IF ( index_cld_a >= ica%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,ica%cur,index_cld_a,zccca,cloud,zplusterm)
            IF ( index_cld_b >= icb%cur) &
                 CALL accumulateSourceReverse(kbdim,klev,nbins,ncld,nspec,kk,icb%cur,index_cld_b,zccca,cloud,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update [fxm]
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,aero,ptstep,zplusterm,zminusterm,zminus_self)

      END DO

    END SUBROUTINE coag_aero

    !
    ! Cloud droplet coagulation
    ! -----------------------------
    !
    SUBROUTINE coag_cloud(kbdim,klev,nspec,ptstep,zcccc,zccca,zccpc,zccic,zccsc)
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zcccc(kbdim,klev,ncld,ncld), zccca(kbdim,klev,nbins,ncld),    &
                          zccpc(kbdim,klev,ncld,nprc), zccic(kbdim,klev,ncld,nice),     &
                          zccsc(kbdim,klev,ncld,nsnw)

      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      ! These are necessary only for coagulation based precip formation
      REAL :: zvolsink_slf(nspec,kbdim,klev)
      REAL :: zvol_prc(nspec,kbdim,klev,nprc)
      REAL :: znum_prc(kbdim,klev,nprc)
      REAL :: zdia(ncld)  ! Gonna be deprecated
      ! ---

      INTEGER :: index_aero_a, index_aero_b
      INTEGER :: index_a, index_b
      INTEGER :: ii,jj,kk

      ! Update the cloud droplet diameters as they are needed later; THIS CAN BE DONE IN A CLEANER WAY IN SUBSEQUENT VERSIONS
      IF (lsauto%state .AND. lsauto%mode == 1) THEN
         DO jj = 1,klev
            DO ii = 1,kbdim
               CALL calcDiamSALSA(ncld,cloud(ii,jj,:),zdia)
               DO kk = 1,ncld
                  cloud(ii,jj,kk)%dwet = MIN(zdia(kk),cloud(ii,jj,kk)%dlim)
               END DO
            END DO
         END DO
      END IF

      DO kk = ica%cur,fca%cur
         IF ( ALL(cloud(:,:,kk)%numc < cloud(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         zvolsink_slf(:,:,:) = 0.
         zvol_prc(:,:,:,:) = 0.
         znum_prc(:,:,:) = 0.

         ! Corresponding index in the regime b droplets
         index_b = icb%cur + (kk-ica%cur)

         ! Collection by larger aerosol bins (reverse collection)
         IF (lscgca .AND. kk < fca%cur) THEN
            ! Corresponding aerosol bin + 1 in a and b
            index_aero_a = ica%par + (kk - ica%cur) + 1
            index_aero_b = icb%par + (kk - ica%cur) + 1
            CALL accumulateSinkReverse(kbdim,klev,ncld,nbins,kk,index_aero_a,fn2a,zccca,aero,zminusterm)
            CALL accumulateSinkReverse(kbdim,klev,ncld,nbins,kk,index_aero_b,fn2b,zccca,aero,zminusterm)
         END IF
         
         ! Collection by larger droplet bins in a and b, and self collection
         IF (lscgcc) THEN
            IF (kk < fca%cur) THEN
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk+1,fca%cur,zcccc,cloud,zminusterm)
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_b+1,fcb%cur,zcccc,cloud,zminusterm)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,zcccc,              &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc )
            ELSE
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk,kk,zcccc,cloud,zminus_self,multp=0.5)
            END IF
         END IF
         
         ! Collection by precipitation
         IF (lscgpc) &
              CALL accumulateSink(kbdim,klev,ncld,nprc,kk,1,nprc,zccpc,precp,zminusterm)
         
         ! Collection by ice
         IF (lscgic) &
              CALL accumulateSink(kbdim,klev,ncld,nice,kk,1,nice,zccic,ice,zminusterm)

         ! Collection by snow
         IF (lscgsc) &
              CALL accumulateSink(kbdim,klev,ncld,nsnw,kk,1,nsnw,zccsc,snow,zminusterm)

         ! Volume gained from smaller and equal aerosol bins
         IF (lscgca) THEN
            ! Corresponding aerosol bins in a and b
            index_aero_a = ica%par + (kk - ica%cur)
            index_aero_b = icb%par + (kk - ica%cur)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in1a,index_aero_a,zccca,aero,zplusterm)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in2b,index_aero_b,zccca,aero,zplusterm)
         END IF

         ! Volume gained from smaller cloud droplets (smaller or equal for b)
         IF (lscgcc .AND. kk > ica%cur) THEN

            IF ( lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,zcccc,cloud,zplusterm)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,zcccc,cloud,zplusterm)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,cloud,ptstep,zplusterm,zminusterm,zminus_self)

         IF (lsauto%state .AND. lsauto%mode == 1) &
              CALL applyCoagPrecipFormation(kbdim,klev,nspec,kk,ptstep,zvolsink_slf,zvol_prc,znum_prc)

      END DO
      
      !
      ! B- bins
      !
      DO kk = icb%cur,fcb%cur
         IF ( ALL(cloud(:,:,kk)%numc < cloud(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         zvolsink_slf(:,:,:) = 0.
         zvol_prc(:,:,:,:) = 0.
         znum_prc(:,:,:) = 0.

         ! Corresponding index in the regime b droplets
         index_a = ica%cur + (kk-icb%cur)

         ! Collection by larger aerosol bins (reverse collection)
         IF (lscgca .AND. kk < fcb%cur) THEN
            ! Corresponding aerosol bin + 1 in a and b
            index_aero_a = ica%par + (kk - icb%cur) + 1
            index_aero_b = icb%par + (kk - icb%cur) + 1
            CALL accumulateSinkReverse(kbdim,klev,ncld,nbins,kk,index_aero_a,fn2a,zccca,aero,zminusterm)
            CALL accumulateSinkReverse(kbdim,klev,ncld,nbins,kk,index_aero_b,fn2b,zccca,aero,zminusterm)
         END IF
         
         ! Collection by larger droplet bins in a and b, and self collection
         IF (lscgcc) THEN
            IF (kk < fcb%cur) THEN
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk+1,fcb%cur,zcccc,cloud,zminusterm)
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_a+1,fca%cur,zcccc,cloud,zminusterm)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,zcccc,              &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc )
            ELSE
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk,kk,zcccc,cloud,zminus_self,multp=0.5)
            END IF
         END IF
         
         ! Collection by precipitation
         IF (lscgpc) &
              CALL accumulateSink(kbdim,klev,ncld,nprc,kk,1,nprc,zccpc,precp,zminusterm)
         
         ! Collection by ice
         IF (lscgic) &
              CALL accumulateSink(kbdim,klev,ncld,nice,kk,1,nice,zccic,ice,zminusterm)

         ! Collection by snow
         IF (lscgsc) &
              CALL accumulateSink(kbdim,klev,ncld,nsnw,kk,1,nsnw,zccsc,snow,zminusterm)

         ! Volume gained from smaller and equal aerosol bins
         IF (lscgca) THEN
            ! Corresponding aerosol bins in a and b
            index_aero_a = ica%par + (kk - icb%cur)
            index_aero_b = icb%par + (kk - icb%cur)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in1a,index_aero_a,zccca,aero,zplusterm)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in2b,index_aero_b,zccca,aero,zplusterm)
         END IF

         ! Volume gained from smaller cloud droplets (smaller or equal for b)
         IF (lscgcc .AND. kk > icb%cur) THEN

            IF ( lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,zcccc,cloud,zplusterm)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,zcccc,cloud,zplusterm)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,cloud,ptstep,zplusterm,zminusterm,zminus_self)

         IF (lsauto%state .AND. lsauto%mode == 1) &
              CALL applyCoagPrecipFormation(kbdim,klev,nspec,kk,ptstep,zvolsink_slf,zvol_prc,znum_prc)


      END DO

    END SUBROUTINE coag_cloud

    !
    ! Precipitation coagulation
    ! --------------------------
    !
    SUBROUTINE coag_precp(kbdim,klev,nspec,ptstep,zccpp,zccpa,zccpc,zccip,zccsp)
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccpp(kbdim,klev,nprc,nprc), zccpa(kbdim,klev,nbins,nprc),    &
                          zccpc(kbdim,klev,ncld,nprc), zccip(kbdim,klev,nprc,nice),     &
                          zccsp(kbdim,klev,nprc,nsnw)
      INTEGER :: kk
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      DO kk = 1,nprc
         IF ( ALL(precp(:,:,kk)%numc < precp(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.

         ! Collection by larger precip and self collection
         IF (lscgpp) THEN
            IF ( kk < nprc) &
                 CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk+1,nprc,zccpp,precp,zminusterm)

            CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk,kk,zccpp,precp,zminus_self,multp=0.5)
         END IF

         ! Collection by ice
         IF (lscgip) &
              CALL accumulateSink(kbdim,klev,nprc,nice,kk,1,nice,zccip,ice,zminusterm)

         ! Collection by snow
         IF (lscgsp) &
              CALL accumulateSink(kbdim,klev,nprc,nsnw,kk,1,nsnw,zccsp,snow,zminusterm)

         ! Volume gained from collection of aerosol
         IF (lscgpa) &
              CALL accumulateSource(kbdim,klev,nprc,nbins,nspec,kk,in1a,fn2b,zccpa,aero,zplusterm)
         
         ! Volume gained from collection of cloud droplets
         IF (lscgpc) &
              CALL accumulateSource(kbdim,klev,nprc,ncld,nspec,kk,1,ncld,zccpc,cloud,zplusterm)

         ! Volume gained from smaller precp
         IF (lscgpp .AND. kk > 1) &
              CALL accumulateSource(kbdim,klev,nprc,nprc,nspec,kk,1,kk-1,zccpp,precp,zplusterm)

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,nprc,nspec,kk,precp,ptstep,zplusterm,zminusterm,zminus_self)

      END DO

    END SUBROUTINE coag_precp

    !
    ! Ice coagulation
    ! ------------------
    ! 
    SUBROUTINE coag_ice(kbdim,klev,nspec,ptstep,zccii,zccia,zccic,zccip,zccsi)

      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccii(kbdim,klev,nprc,nprc), zccia(kbdim,klev,nbins,nprc),    &
                          zccic(kbdim,klev,ncld,nprc), zccip(kbdim,klev,nprc,nice),     &
                          zccsi(kbdim,klev,nprc,nsnw)

      INTEGER :: kk, index_b, index_a
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)
      INTEGER :: nwet,ndry,iwa
      REAL :: rhowa,rhoic

      nwet = spec%getNSpec(type="wet")
      ndry = spec%getNSpec(type="dry")
      iwa = spec%getIndex("H2O")
      rhowa = spec%rhowa
      rhoic = spec%rhoic

      DO kk = iia%cur, fia%cur
         IF (ALL(ice(:,:,kk)%numc < ice(:,:,kk)%nlim)) CYCLE
         
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.

         ! Corresponding index for ice in regime b
         index_b = iib%cur + (kk-iia%cur)
         
         ! Collection by larger ice and self coagulation
         IF (lscgii) THEN
            IF (kk < fia%cur) THEN
               CALL accumulateSink(kbdim,klev,nice,nice,kk,kk+1,fia%cur,zccii,ice,zminusterm)
               CALL accumulateSink(kbdim,klev,nice,nice,kk,index_b+1,fib%cur,zccii,ice,zminusterm)
            END IF
            CALL accumulateSink(kbdim,klev,nice,nice,kk,kk,kk,zccii,ice,zminus_self,multp=0.5)
         END IF

         ! Collection by snow
         IF (lscgsi) &
              CALL accumulateSink(kbdim,klev,nice,nsnw,kk,1,nsnw,zccsi,snow,zminusterm)
         
         ! Volume gained from aerosol collection
         IF (lscgia) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nbins,ndry,iwa,kk,in1a,fn2b,  &
                                               rhoic,rhowa,zccia,aero,zplusterm)

         ! Volume gained from cloud collection
         IF (lscgic) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,ncld,ndry,iwa,kk,ica%cur,fcb%cur, &
                                               rhoic,rhowa,zccic,cloud,zplusterm)

         ! Volume gained from precip collection
         IF (lscgip) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nprc,ndry,iwa,kk,1,nprc,  &
                                               rhoic,rhowa,zccip,precp,zplusterm)

         ! Volume gained from smaller ice particles
         IF (lscgii .AND. kk > iia%cur) THEN
            CALL accumulateSource(kbdim,klev,nice,nice,nspec,kk,iia%cur,kk-1,zccii,ice,zplusterm)
            CALL accumulateSource(kbdim,klev,nice,nice,nspec,kk,iib%cur,index_b-1,zccii,ice,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,nice,nspec,kk,ice,ptstep,zplusterm,zminusterm,zminus_self)

      END DO

      DO kk = iib%cur, fib%cur
         IF (ALL(ice(:,:,kk)%numc < ice(:,:,kk)%nlim)) CYCLE
         
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.

         ! Corresponding index for ice in regime b
         index_a = iia%cur + (kk-iib%cur)
         
         ! Collection by larger ice and self coagulation
         IF (lscgii) THEN
            IF (kk < fib%cur) THEN
               CALL accumulateSink(kbdim,klev,nice,nice,kk,kk+1,fib%cur,zccii,ice,zminusterm)
               CALL accumulateSink(kbdim,klev,nice,nice,kk,index_a+1,fia%cur,zccii,ice,zminusterm)
            END IF
            CALL accumulateSink(kbdim,klev,nice,nice,kk,kk,kk,zccii,ice,zminus_self,multp=0.5)
         END IF

         ! Collection by snow
         IF (lscgsi) &
              CALL accumulateSink(kbdim,klev,nice,nsnw,kk,1,nsnw,zccsi,snow,zminusterm)
         
         ! Volume gained from aerosol collection
         IF (lscgia) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nbins,ndry,iwa,kk,in1a,fn2b,  &
                                               rhoic,rhowa,zccia,aero,zplusterm)

         ! Volume gained from cloud collection
         IF (lscgic) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,ncld,ndry,iwa,kk,ica%cur,fcb%cur, &
                                               rhoic,rhowa,zccic,cloud,zplusterm)

         ! Volume gained from precip collection
         IF (lscgip) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nprc,ndry,iwa,kk,1,nprc,  &
                                               rhoic,rhowa,zccip,precp,zplusterm)

         ! Volume gained from smaller ice particles
         IF (lscgii .AND. kk > iib%cur) THEN
            CALL accumulateSource(kbdim,klev,nice,nice,nspec,kk,iib%cur,kk-1,zccii,ice,zplusterm)
            CALL accumulateSource(kbdim,klev,nice,nice,nspec,kk,iia%cur,index_a-1,zccii,ice,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,nice,nspec,kk,ice,ptstep,zplusterm,zminusterm,zminus_self)

      END DO

    END SUBROUTINE coag_ice
    
    ! -----------------------------------------------------------------

    SUBROUTINE coag_snow(kbdim,klev,nspec,ptstep,zccss,zccsa,zccsc,zccsp,zccsi)
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccss(kbdim,klev,nsnw,nsnw), zccsa(kbdim,klev,nbins,nsnw),    &
                          zccsc(kbdim,klev,ncld,nsnw), zccsp(kbdim,klev,nprc,nsnw),     &
                          zccsi(kbdim,klev,nice,nsnw)
      INTEGER :: kk
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      INTEGER :: nwet,ndry,iwa
      REAL :: rhosn,rhoic,rhowa

      nwet = spec%getNSpec(type='wet')
      ndry = spec%getNSpec(type='dry')
      iwa = spec%getIndex('H2O')
      
      rhosn = spec%rhosn
      rhoic = spec%rhoic
      rhowa = spec%rhowa

      DO kk = 1,nsnw
         IF (ALL(snow(:,:,kk)%numc < snow(:,:,kk)%nlim)) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.

         ! Collection by larger snow and self coagulation
         IF (lscgss) THEN
            IF (kk < nsnw) &
                 CALL accumulateSink(kbdim,klev,nsnw,nsnw,kk,kk+1,nsnw,zccss,snow,zminusterm)
            CALL accumulateSink(kbdim,klev,nsnw,nsnw,kk,kk,kk,zccss,snow,zminus_self,multp=0.5)
         END IF

         ! Volume gained from aerosol
         IF (lscgsa) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nsnw,nbins,ndry,iwa,kk,in1a,fn2b,  &
                                               rhosn,rhowa,zccsa,aero,zplusterm)
         
         ! Volume gained from cloud droplets
         IF (lscgsc) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nsnw,ncld,ndry,iwa,kk,ica%cur,fcb%cur, &
                                               rhosn,rhowa,zccsc,cloud,zplusterm)

         ! Volume gained from precip
         IF (lscgsp) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nsnw,nprc,ndry,iwa,kk,1,nprc, &
                                               rhosn,rhowa,zccsp,precp,zplusterm)              
         
         ! Volume gained from ice (Btw the density change does not really work this way for snow...)
         IF (lscgsi) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nsnw,nice,ndry,iwa,kk,iia%cur,fib%cur, &
                                               rhosn,rhoic,zccsi,ice,zplusterm)

         ! Volume gained from smaller snow
         IF (lscgss .AND. kk > 1) &
              CALL accumulateSource(kbdim,klev,nsnw,nsnw,nspec,kk,1,kk-1,zccss,snow,zplusterm)

         CALL applyCoag(kbdim,klev,nsnw,nspec,kk,snow,ptstep,zplusterm,zminusterm,zminus_self)

      END DO
     
    END SUBROUTINE coag_snow

    ! -----------------------------------------------------------------

    SUBROUTINE applyCoag(kbdim,klev,nb,nspec,itrgt,part,ptstep,source,sink,sink_self)
      INTEGER,INTENT(in) :: kbdim,klev,nb,nspec,itrgt
      REAL,INTENT(in)    :: ptstep
      TYPE(Section),INTENT(inout) :: part(kbdim,klev,nb)
      REAL,INTENT(in) :: source(nspec,kbdim,klev), sink(kbdim,klev), sink_self(kbdim,klev)
 
      INTEGER :: ii,jj

      DO jj = 1,klev
         DO ii = 1,kbdim
            part(ii,jj,itrgt)%volc(1:nspec) =          &
                 ( part(ii,jj,itrgt)%volc(1:nspec) +   &
                   ptstep*source(1:nspec,ii,jj)*part(ii,jj,itrgt)%numc ) / &
                 ( 1. + ptstep*sink(ii,jj) )

            part(ii,jj,itrgt)%numc = part(ii,jj,itrgt)%numc / ( 1. + ptstep*(sink(ii,jj) + sink_self(ii,jj)) )

         END DO
      END DO

    END SUBROUTINE applyCoag

    ! -----------------------------------------------------------------

    SUBROUTINE applyCoagPrecipFormation(kbdim,klev,nspec,itrgt,ptstep,volsink_slf,vol_prc,num_prc)
      ! 
      ! Make the necessary contributions from the coagulation based precip formation method
      !
      INTEGER, INTENT(in) :: kbdim,klev
      INTEGER, INTENT(in) :: nspec  ! Number of compounds
      INTEGER, INTENT(in) :: itrgt
      REAL, INTENT(in)    :: ptstep
      REAL, INTENT(in)    :: volsink_slf(nspec,kbdim,klev)
      REAL, INTENT(in)    :: vol_prc(nspec,kbdim,klev,nprc)
      REAL, INTENT(in)    :: num_prc(kbdim,klev,nprc)
       
      INTEGER :: ii,jj,cc

      DO jj = 1,klev
         DO ii = 1,kbdim
            ! Sink in cloud droplet volume due to precipitation formation via self coagulation
            cloud(ii,jj,itrgt)%volc(1:nspec) = cloud(ii,jj,itrgt)%volc(1:nspec) -   &
                 ptstep*volsink_slf(1:nspec,ii,jj)
         END DO
      END DO
      
      DO cc = 1,nprc
         DO jj = 1,klev
            DO ii = 1,kbdim
               ! Do the contribution to the precipitation bins
               precp(ii,jj,cc)%volc(1:nspec) = precp(ii,jj,cc)%volc(1:nspec) +     &
                    ptstep*vol_prc(1:nspec,ii,jj,cc)

               precp(ii,jj,cc)%numc = precp(ii,jj,cc)%numc +    &
                    ptstep*num_prc(ii,jj,cc)
            END DO
         END DO
      END DO

    END SUBROUTINE applyCoagPrecipFormation


    ! -----------------------------------------------------------------

    SUBROUTINE accumulateSink(kbdim,klev,nbtrgt,nbcoll,itrgt,istr,iend,zcc,coll,sink,multp)
      ! 
      ! The "direct" method, i.e. "larger" particle category collects "smaller" category.
      ! Here, the target refers always to the collected "smaller" particle category.
      !
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt,nbcoll    ! Number of bins in the target and the collector categories
      INTEGER,INTENT(in) :: itrgt,istr,iend  ! Index of the target bin, start and end indices of the collector bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbtrgt,nbcoll)
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll) ! Collector particles properties
      REAL,INTENT(inout) :: sink(kbdim,klev)
      REAL,INTENT(in), OPTIONAL :: multp

      INTEGER :: ll,ii,jj
      REAL :: xx

      ! For self collection
      IF ( PRESENT(multp) ) THEN
         xx = multp
      ELSE
         xx = 1.0
      END IF

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               sink(ii,jj) = sink(ii,jj) + xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSink
    ! --------------------------------------------
    SUBROUTINE accumulateSinkReverse(kbdim,klev,nbtrgt,nbcoll,itrgt,istr,iend,zcc,coll,sink)
      !
      ! The reverse method, where the "smaller" particle category collects the larger category.
      ! The target refers always to the collected "larger" category
      !
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt,nbcoll   ! Number of bins in the target and the collector categories
      INTEGER,INTENT(in) :: itrgt,istr,iend ! Index of the target bin, start and end indices of the collector bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbcoll,nbtrgt)
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll)
      REAL,INTENT(inout) :: sink(kbdim,klev)
      
      INTEGER :: ll,ii,jj

      DO ll = istr,iend
         DO jj =  1,klev
            DO ii = 1,kbdim
               sink(ii,jj) = sink(ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%numc
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSinkReverse
    ! --------------------------------------------
    SUBROUTINE accumulateSource(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,coll,source)
      !  
      ! The direct method, where the "larger" particle category collects the "smaller" category.
      ! The target refers always to the "larger", collector category.
      ! 
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt, nbcoll   ! Number of bins in the target and collectee categories
      INTEGER,INTENT(in) :: nspec           ! Number of chemical compounds
      INTEGER,INTENT(in) :: itrgt,istr,iend  ! Index of the target bin, start and end indices of the collectee bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbcoll,nbtrgt)
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll) ! Collected particle properties
      REAL, INTENT(inout) :: source(nspec,kbdim,klev)

      INTEGER :: ll,ii,jj

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)               
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSource
    ! ---------
    SUBROUTINE accumulateSourceReverse(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,coll,source)
      !
      ! The reverse method, where the "smaller" particle category collects the larger category.
      ! The target refers always to the "smaller", collector category
      !
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt, nbcoll
      INTEGER,INTENT(in) :: nspec
      INTEGER,INTENT(in) :: itrgt,istr,iend
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbtrgt,nbcoll)
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll)
      REAL, INTENT(inout) :: source(nspec,kbdim,klev)

      INTEGER :: ll,ii,jj

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%volc(1:nspec)
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSourceReverse
    !----------
    SUBROUTINE accumulateSourcePhaseChange(kbdim,klev,nbtrgt,nbcoll,ndry,iwa,itrgt,istr,iend,  &
                                           rhotrgt,rhocoll,zcc,coll,source)
      !
      ! The direct method, where the "larger" particle category collects the "smaller" category.
      ! The target refers always to the "larger", collector category.
      ! 
      ! This subroutine takes into account the change in water density due to freezing upon collection
      ! 
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt, nbcoll   ! Number of bins in the target and collectee categories
      INTEGER,INTENT(in) :: ndry,iwa           ! Number of dry checmical compounds, index for water 
      INTEGER,INTENT(in) :: itrgt,istr,iend  ! Index of the target bin, start and end indices of the collectee bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbcoll,nbtrgt)
      REAL, INTENT(in)   :: rhotrgt,rhocoll      ! Water densities for the target and collected categories 
                                                 ! (typically liquid and frozen, respectively).
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll) ! Collected particle properties
      REAL, INTENT(inout) :: source(iwa,kbdim,klev)

      INTEGER :: ll,ii,jj

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:ndry,ii,jj) = source(1:ndry,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:ndry)
               source(iwa,ii,jj) = source(iwa,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*rhocoll/rhotrgt 
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSourcePhaseChange

    ! ---------------------------------
    
    ! Category specific processes
    ! ---------------------------------
    SUBROUTINE accumulatePrecipFormation(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,     &
                                         source, sink, volsink_slf, vol_prc, num_prc             )
      !
      ! This method for collision-coalescence transfers the resulting
      ! droplets directly to precipitation bins, if the resulting diameter is 
      ! larger than the first precip bin diameter. This is done via separate source terms
      ! that are output from this subroutine.
      !
      INTEGER, INTENT(in) :: kbdim,klev          
      INTEGER, INTENT(in) :: nbtrgt, nbcoll   ! Number of bins in the target and collected categories
      INTEGER, INTENT(in) :: nspec            ! Number of compounds
      INTEGER, INTENT(in) :: itrgt,istr,iend  ! Index of the target (cloud droplet) bin, start and end indices for the collected bins
      REAL, INTENT(in)      :: zcc(kbdim,klev,nbcoll,nbtrgt)    ! Collision kernels
      REAL, INTENT(inout)   :: source(nspec,kbdim,klev)         ! Regular coagulation source term
      REAL, INTENT(inout)   :: sink(kbdim,klev)           ! Regular coagulation sink term
      REAL, INTENT(inout)   :: volsink_slf(nspec,kbdim,klev)    ! Volume sink for cloud droplets upon self collection resulting in precipitation formation
      REAL, INTENT(inout)   :: vol_prc(nspec,kbdim,klev,nprc)   ! Volume source for precipitation as the result of collision coalescence
      REAL, INTENT(inout)   :: num_prc(kbdim,klev,nprc)         ! Number source for precipitation as the result of collision coalescence
      
      REAL :: D_new
      INTEGER :: trgt_prc
      INTEGER :: ii,jj,ll
      REAL :: xx
      LOGICAL :: selfcoll

      selfcoll = .FALSE.
      IF (itrgt == istr .AND. itrgt == iend) THEN
         ! Self collection
         selfcoll = .TRUE.
      END IF
 
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               ! The estimated diameter of the droplets after collision
                D_new = (cloud(ii,jj,itrgt)%dwet**3 + cloud(ii,jj,ll)%dwet**3)**(1./3.)

                ! Check out to which precip bin this belogs (if any)
                trgt_prc = 0
                trgt_prc = COUNT( D_new > precpbins(:) )
                
                IF ( trgt_prc > 0) THEN

                   ! The resulting droplets are large -> put the collision coalescence contribution to precipitation
                   IF ( selfcoll ) THEN
                      ! Self collection
                      vol_prc(1:nspec,ii,jj,trgt_prc) = vol_prc(1:nspec,ii,jj,trgt_prc) +   &
                           cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc*zcc(ii,jj,ll,itrgt) 

                      num_prc(ii,jj,trgt_prc) = num_prc(ii,jj,trgt_prc) +        &
                           (0.5*zcc(ii,jj,ll,itrgt)*cloud(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc) 
                   
                   ELSE
                      vol_prc(1:nspec,ii,jj,trgt_prc) = vol_prc(1:nspec,ii,jj,trgt_prc) +    &
                           ( cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc +         &
                             cloud(ii,jj,itrgt)%volc(1:nspec)*cloud(ii,jj,ll)%numc ) * zcc(ii,jj,ll,itrgt)

                      num_prc(ii,jj,trgt_prc) = num_prc(ii,jj,trgt_prc) +        &
                           (zcc(ii,jj,ll,itrgt)*cloud(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc) 

                   END IF
                   
                   IF ( selfcoll ) THEN
                      ! Change in cloud droplet volume due to precip formation in self collection
                      volsink_slf(1:nspec,ii,jj) = volsink_slf(1:nspec,ii,jj) +   &
                           cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc*zcc(ii,jj,ll,itrgt)
                      
                      ! Contribution of precip formation due to self collection to the regular sink term
                      sink(ii,jj) = sink(ii,jj) + zcc(ii,jj,ll,itrgt)*cloud(ii,jj,itrgt)%numc

                   ELSE
                      ! Precip formation consumes particles from the target cloud droplet bin also when not self collection, which is not accounted 
                      ! for by the regular sink term accumulation
                      sink(ii,jj) = sink(ii,jj) + zcc(ii,jj,itrgt,ll)*cloud(ii,jj,ll)%numc

                   END IF

                ELSE
                   
                   ! No precipitation formation -> regular approach
                   IF ( selfcoll ) THEN
                      ! Self collection
                      sink(ii,jj) = sink(ii,jj) + 0.5*zcc(ii,jj,itrgt,ll)*cloud(ii,jj,ll)%numc
                   ELSE
                      ! Not self collection:                
                      source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + zcc(ii,jj,ll,itrgt)*cloud(ii,jj,ll)%volc(1:nspec)

                   END IF

                END IF

            END DO
         END DO
      END DO

    END SUBROUTINE accumulatePrecipFormation
    
END MODULE mo_salsa_coagulation_processes
