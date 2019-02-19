MODULE mo_salsa_coagulation_processes
  USE mo_salsa_types, ONLY : aero, cloud, precp, ice, rateDiag
  USE mo_submctl, ONLY : in1a, fn1a, in2a, fn2a, in2b, fn2b, nbins,  &
                         ica, fca, icb, fcb, ncld,                  &
                         nprc,                                      &
                         iia, fia, nice, & 
                         precpbins,                                 &
                         lscgaa, lscgcc, lscgpp, lscgii, & 
                         lscgca, lscgpa, lscgia, & 
                         lscgpc, lscgic, & 
                         lscgip, & 
                         lsauto, &
                         spec
  USE classSection, ONLY : Section
  IMPLICIT NONE

  CONTAINS

    !
    ! Aerosol coagulation
    ! -----------------------
    SUBROUTINE coag_aero(kbdim,klev,nspec,ptstep,zccaa,zccca,zccpa,zccia) 
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccaa(kbdim,klev,nbins,nbins), zccca(kbdim,klev,nbins,ncld),    &
                          zccpa(kbdim,klev,nbins,nprc), zccia(kbdim,klev,nbins,nice)

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
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2b,zccaa,aero,zminusterm,1)
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,1,multp=0.5)
         end if

         ! Collection by larger and equal cloud droplet bins (regime a and b)
         IF (lscgca) THEN
            CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm,1)
            CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm,1)
         end if

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm,1)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia,fia,zccia,ice,zminusterm,1)
              
         ! Particle volume gained from smaller particles in regime 1a
         IF (lscgaa .AND. kk > in1a) &
              CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,kk-1,zccaa,aero,zplusterm,1)

         ! Particle volume gained from smaller cloud droplet bins (does not happen with the current configuration,
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
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,ptstep,aero,zplusterm,zminusterm,zminus_self)

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
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2a,zccaa,aero,zminusterm,1)
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,index_b,fn2b,zccaa,aero,zminusterm,1)
            END IF

            ! Self collection
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,1,multp=0.5)
         END IF

         ! Collection by larger or equal cloud droplet bins
         IF (lscgca) THEN
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm,1)
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm,1)
         END IF

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm,1)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia,fia,zccia,ice,zminusterm,1)

         ! Particle volume gained from smaller particles in regime a and b
         IF (lscgaa) THEN
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,kk-1,zccaa,aero,zplusterm,1)
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in2b,index_b-1,zccaa,aero,zplusterm,1)
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
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,ptstep,aero,zplusterm,zminusterm,zminus_self)

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
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk+1,fn2b,zccaa,aero,zminusterm,1)
               CALL accumulateSink(kbdim,klev,nbins,nbins,kk,index_a,fn2a,zccaa,aero,zminusterm,1)
            END IF

            ! Self collection
            CALL accumulateSink(kbdim,klev,nbins,nbins,kk,kk,kk,zccaa,aero,zminus_self,1,multp=0.5)
         END IF

         ! Collection by larger or equal cloud droplet bins
         IF (lscgca) THEN
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_a,fca%cur,zccca,cloud,zminusterm,1)
              CALL accumulateSink(kbdim,klev,nbins,ncld,kk,index_cld_b,fcb%cur,zccca,cloud,zminusterm,1)
         END IF

         ! Collection by precipitation
         IF (lscgpa) &
              CALL accumulateSink(kbdim,klev,nbins,nprc,kk,1,nprc,zccpa,precp,zminusterm,1)

         ! Collection by ice
         IF (lscgia) &
              CALL accumulateSink(kbdim,klev,nbins,nice,kk,iia,fia,zccia,ice,zminusterm,1)
              
         ! Particle volume gained from smaller particles in regime a and b
         IF (lscgaa) THEN
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in2b,kk-1,zccaa,aero,zplusterm,1)
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,index_a-1,zccaa,aero,zplusterm,1)
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
         CALL applyCoag(kbdim,klev,nbins,nspec,kk,ptstep,aero,zplusterm,zminusterm,zminus_self)

      END DO

    END SUBROUTINE coag_aero

    !
    ! Cloud droplet coagulation
    ! -----------------------------
    !
    SUBROUTINE coag_cloud(kbdim,klev,nspec,ptstep,zcccc,zccca,zccpc,zccic)
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zcccc(kbdim,klev,ncld,ncld), zccca(kbdim,klev,nbins,ncld),    &
                          zccpc(kbdim,klev,ncld,nprc), zccic(kbdim,klev,ncld,nice)

      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      ! These are necessary only for coagulation based precip formation
      REAL :: zvolsink_slf(nspec,kbdim,klev)
      REAL :: zvol_prc(nspec,kbdim,klev,nprc)
      REAL :: znum_prc(kbdim,klev,nprc)
      ! ---

      INTEGER :: index_aero_a, index_aero_b
      INTEGER :: index_a, index_b
      INTEGER :: ii,jj,kk

      ! Update the cloud droplet diameters as they are needed later; THIS CAN BE DONE IN A CLEANER WAY IN SUBSEQUENT VERSIONS
      IF (lsauto%state .AND. lsauto%mode == 1) THEN
         DO kk = 1,ncld
            DO jj = 1,klev
               DO ii = 1,kbdim
                  CALL cloud(ii,jj,kk)%updateDiameter(.TRUE.)
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
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk+1,fca%cur,zcccc,cloud,zminusterm,2)
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_b+1,fcb%cur,zcccc,cloud,zminusterm,2)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,zcccc,              &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc )
            ELSE
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk,kk,zcccc,cloud,zminus_self,2,multp=0.5)
            END IF
         END IF
         
         ! Collection by precipitation
         IF (lscgpc) &
              CALL accumulateSink(kbdim,klev,ncld,nprc,kk,1,nprc,zccpc,precp,zminusterm,2)
         
         ! Collection by ice
         IF (lscgic) &
              CALL accumulateSink(kbdim,klev,ncld,nice,kk,iia,fia,zccic,ice,zminusterm,2)

         ! Volume gained from smaller and equal aerosol bins
         IF (lscgca) THEN
            ! Corresponding aerosol bins in a and b
            index_aero_a = ica%par + (kk - ica%cur)
            index_aero_b = icb%par + (kk - ica%cur)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in1a,index_aero_a,zccca,aero,zplusterm,2)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in2b,index_aero_b,zccca,aero,zplusterm,2)
         END IF

         ! Volume gained from smaller cloud droplets (smaller or equal for b)
         IF (lscgcc .AND. kk > ica%cur) THEN

            IF ( lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,zcccc,cloud,zplusterm,2)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,zcccc,cloud,zplusterm,2)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,ptstep,cloud,zplusterm,zminusterm,zminus_self)

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
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk+1,fcb%cur,zcccc,cloud,zminusterm,2)
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_a+1,fca%cur,zcccc,cloud,zminusterm,2)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,zcccc,              &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc )
            ELSE
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,kk,kk,zcccc,cloud,zminus_self,2,multp=0.5)
            END IF
         END IF
         
         ! Collection by precipitation
         IF (lscgpc) &
              CALL accumulateSink(kbdim,klev,ncld,nprc,kk,1,nprc,zccpc,precp,zminusterm,2)
         
         ! Collection by ice
         IF (lscgic) &
              CALL accumulateSink(kbdim,klev,ncld,nice,kk,iia,fia,zccic,ice,zminusterm,2)

         ! Volume gained from smaller and equal aerosol bins
         IF (lscgca) THEN
            ! Corresponding aerosol bins in a and b
            index_aero_a = ica%par + (kk - icb%cur)
            index_aero_b = icb%par + (kk - icb%cur)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in1a,index_aero_a,zccca,aero,zplusterm,2)
            CALL accumulateSource(kbdim,klev,ncld,nbins,nspec,kk,in2b,index_aero_b,zccca,aero,zplusterm,2)
         END IF

         ! Volume gained from smaller cloud droplets (smaller or equal for b)
         IF (lscgcc .AND. kk > icb%cur) THEN

            IF ( lsauto%state .AND. lsauto%mode == 1) THEN
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,zcccc,cloud,zplusterm,2)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,zcccc,cloud,zplusterm,2)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,ptstep,cloud,zplusterm,zminusterm,zminus_self)

         IF (lsauto%state .AND. lsauto%mode == 1) &
              CALL applyCoagPrecipFormation(kbdim,klev,nspec,kk,ptstep,zvolsink_slf,zvol_prc,znum_prc)

      END DO

    END SUBROUTINE coag_cloud

    !
    ! Precipitation coagulation
    ! --------------------------
    !
    SUBROUTINE coag_precp(kbdim,klev,nspec,ptstep,zccpp,zccpa,zccpc,zccip) 
      
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccpp(kbdim,klev,nprc,nprc), zccpa(kbdim,klev,nbins,nprc),    &
                          zccpc(kbdim,klev,ncld,nprc), zccip(kbdim,klev,nprc,nice)

      INTEGER :: kk,ii,jj
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)
      REAL :: dzplus(nspec,kbdim,klev)
      
      DO kk = 1,nprc
         IF ( ALL(precp(:,:,kk)%numc < precp(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         dzplus(:,:,:) = 0.
         
         ! Collection by larger precip and self collection
         IF (lscgpp) THEN
            IF ( kk < nprc ) &
                 CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk+1,nprc,zccpp,precp,zminusterm,3)

            CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk,kk,zccpp,precp,zminus_self,3,multp=0.5)
         END IF

         ! Collection by ice
         IF (lscgip) &
              CALL accumulateSink(kbdim,klev,nprc,nice,kk,iia,fia,zccip,ice,zminusterm,3)

         ! Volume gained from collection of aerosol
         IF (lscgpa) &
              CALL accumulateSource(kbdim,klev,nprc,nbins,nspec,kk,in1a,fn2b,zccpa,aero,zplusterm,3)
         dzplus = zplusterm
         
         ! Volume gained from collection of cloud droplets
         IF (lscgpc) &
              CALL accumulateSource(kbdim,klev,nprc,ncld,nspec,kk,1,ncld,zccpc,cloud,zplusterm,3)
         
         ! Volume gained from smaller precp
         IF (lscgpp .AND. kk > 1) &
              CALL accumulateSource(kbdim,klev,nprc,nprc,nspec,kk,1,kk-1,zccpp,precp,zplusterm,3)

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,nprc,nspec,kk,ptstep,precp,zplusterm,zminusterm,zminus_self)

      END DO

    END SUBROUTINE coag_precp

    !
    ! Ice coagulation
    ! ------------------
    ! 
    SUBROUTINE coag_ice(kbdim,klev,nspec,ptstep,zccii,zccia,zccic,zccip)

      INTEGER, INTENT(in) :: kbdim,klev,nspec   ! nspec should contain all compounds including unrimed and rimed ice!
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: zccii(kbdim,klev,nice,nice), zccia(kbdim,klev,nbins,nice),    &
                          zccic(kbdim,klev,ncld,nice), zccip(kbdim,klev,nprc,nice)

      INTEGER :: kk, index_b, index_a
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev) 
      INTEGER :: iwa,irim
      REAL :: rhowa,rhoic,rhorime

      iwa = spec%getIndex("H2O")
      irim = spec%getIndex("rime")
      rhowa = spec%rhowa
      rhoic = spec%rhoic
      rhorime = spec%rhori

      DO kk = iia, fia
         IF (ALL(ice(:,:,kk)%numc < ice(:,:,kk)%nlim)) CYCLE
         
         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         
         ! Collection by larger ice and self coagulation
         IF (lscgii) THEN
            IF (kk < fia) THEN
               CALL accumulateSink(kbdim,klev,nice,nice,kk,kk+1,fia,zccii,ice,zminusterm,4)
            END IF
            CALL accumulateSink(kbdim,klev,nice,nice,kk,kk,kk,zccii,ice,zminus_self,4,multp=0.5)
         END IF

         ! Volume gained from aerosol collection. Assume this to produce pristine ice
         IF (lscgia) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nbins,nspec,iwa,iwa,kk,in1a,fn2b,  &
                                               rhoic,rhowa,zccia,aero,zplusterm)

         ! Volume gained from cloud collection. Produces rimed ice
         IF (lscgic) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,ncld,nspec,irim,iwa,kk,ica%cur,fcb%cur, &
                                               rhorime,rhowa,zccic,cloud,zplusterm)

         ! Volume gained from precip collection. Produces rimed ice
         IF (lscgip) &
              CALL accumulateSourcePhaseChange(kbdim,klev,nice,nprc,nspec,irim,iwa,kk,1,nprc,     &
                                               rhorime,rhowa,zccip,precp,zplusterm)

         ! Volume gained from smaller ice particles.
         IF (lscgii .AND. kk > iia) THEN
            CALL accumulateSourceIce(kbdim,klev,nice,nice,nspec,kk,iia,kk-1, &
                                     zccii,ice,zplusterm)
         END IF

         !-- Volume and number concentrations after coagulation update. Use the Ice variant, which 
         !   calculates the volume mean densities etc correctly for the pristine and rimed ice contributions.
         CALL applyCoagIce(kbdim,klev,nice,nspec,iwa,irim,kk,ptstep,ice,      &
                           zplusterm,zminusterm,zminus_self  )

      END DO

    END SUBROUTINE coag_ice
    
    ! -----------------------------------------------------------------

    SUBROUTINE applyCoag(kbdim,klev,nb,nspec,itrgt,ptstep,part,source,sink,sink_self)
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
 ! ----------------------------------------
    SUBROUTINE applyCoagIce(kbdim,klev,nb,nspec,iwa,irim,itrgt,ptstep,part,  &
                            source,sink,sink_self)
      ! --------------------------------------------------------------------------------------
      ! Note: The mass conversions and effective densities only account for the ice part.
      ! Not the dry aerosol!!! Is this a good enough approximation?
      ! --------------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: kbdim,klev,nb,nspec,iwa,irim,itrgt  ! nspec should contain all compounds including both unrimed and rimed ice. Necessary for memory allocation of input arrays
      REAL, INTENT(in)    :: ptstep
      TYPE(Section), INTENT(inout) :: part(kbdim,klev,nb)
      REAL, INTENT(in) :: source(nspec,kbdim,klev), sink(kbdim,klev), sink_self(kbdim,klev)

      REAL :: mtrgt_t, mtrgt_r, mtot, mrime ! The ice mass in target particle, the mass source term for total ice
                                            ! and rimed ice
      INTEGER :: ii,jj
      INTEGER :: ndry

      ndry = spec%getNSpec(type="dry")  ! Number of dry species
      
      DO jj = 1,klev
         DO ii = 1,kbdim

            ! Apply the change due to coagulation to the dry components. Standard procedure
            part(ii,jj,itrgt)%volc(1:ndry) =                &
                 ( part(ii,jj,itrgt)%volc(1:ndry) +          &
                   ptstep*source(1:ndry,ii,jj)*part(ii,jj,itrgt)%numc ) /  &
                 ( 1. + ptstep*sink(ii,jj) )

            ! Apply the coagulation sources and sinks of particle volume concentrations
            part(ii,jj,itrgt)%volc(irim) = ( part(ii,jj,itrgt)%volc(irim) +     &
                                        ptstep*source(irim,ii,jj)*part(ii,jj,itrgt)%numc ) / &
                                      ( 1. + ptstep*sink(ii,jj) )
            part(ii,jj,itrgt)%volc(iwa) = ( part(ii,jj,itrgt)%volc(iwa) +  &
                                            ptstep*source(iwa,ii,jj)*part(ii,jj,itrgt)%numc ) / &
                                          ( 1. + ptstep*sink(ii,jj) )

            ! Apply the sink term for number concentration
            part(ii,jj,itrgt)%numc = part(ii,jj,itrgt)%numc / ( 1. + ptstep*(sink(ii,jj) + sink_self(ii,jj)) )


            ! Update the mean particle density
            CALL part(ii,jj,itrgt)%updateRhomean()

         END DO
      END DO

      
    END SUBROUTINE applyCoagIce
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
      REAL :: vrate(nspec), nrate

      vrate = 0.
      nrate = 0.
      
      DO jj = 1,klev
         DO ii = 1,kbdim
            ! Sink in cloud droplet volume due to precipitation formation via self coagulation
            vrate(1:nspec) = volsink_slf(1:nspec,ii,jj)
            cloud(ii,jj,itrgt)%volc(1:nspec) = cloud(ii,jj,itrgt)%volc(1:nspec) -   &
                 ptstep*vrate(1:nspec)
         END DO
      END DO
      
      DO cc = 1,nprc
         DO jj = 1,klev
            DO ii = 1,kbdim
               ! Do the contribution to the precipitation bins
               vrate(1:nspec) = vol_prc(1:nspec,ii,jj,cc)
               precp(ii,jj,cc)%volc(1:nspec) = precp(ii,jj,cc)%volc(1:nspec) +     &
                    ptstep*vrate(1:nspec)

               nrate = num_prc(ii,jj,cc)
               precp(ii,jj,cc)%numc = precp(ii,jj,cc)%numc +    &
                    ptstep*nrate
            END DO
         END DO
      END DO

    END SUBROUTINE applyCoagPrecipFormation


    ! -----------------------------------------------------------------

    SUBROUTINE accumulateSink(kbdim,klev,nbtrgt,nbcoll,itrgt,istr,iend,zcc,coll,sink,trgtphase,multp)
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
      INTEGER, INTENT(in) :: trgtphase    ! phase indentifier for the target particles
      REAL,INTENT(in), OPTIONAL :: multp

      INTEGER :: ll,ii,jj
      REAL :: xx
      REAL :: nrate(kbdim,klev)
      REAL :: dnum
      
      nrate = 0.
      dnum = 0.
      
      ! For self collection
      IF ( PRESENT(multp) ) THEN
         xx = multp
      ELSE
         xx = 1.0
      END IF

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               dnum = xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
               sink(ii,jj) = sink(ii,jj) + dnum
               nrate(ii,jj) = nrate(ii,jj) + dnum
            END DO
         END DO
      END DO

      ! Diagnostics
      IF (coll(1,1,1)%phase == 3 .AND. trgtphase == 2) THEN
         ! Accretion -- number sink term
         nrate(1,1) = nrate(1,1) * cloud(1,1,itrgt)%numc
         CALL rateDiag%Accretion%accumulate(n=nrate(1,1))
      ELSE IF (coll(1,1,1)%phase == 2 .AND. trgtphase == 1) THEN
         ! Cloud collection of aerosol -- number sink
         nrate(1,1) = nrate(1,1) * aero(1,1,itrgt)%numc
         CALL rateDiag%ACcoll%accumulate(n=nrate(1,1))
      ELSE IF (coll(1,1,1)%phase == 3 .AND. trgtphase == 1) THEN
         ! Precipitation collection of aerosol -- number sink
         nrate(1,1) = nrate(1,1) * aero(1,1,itrgt)%numc
         CALL rateDiag%APcoll%accumulate(n=nrate(1,1))
      ELSE IF (coll(1,1,1)%phase == 4 .AND. trgtphase == 1) THEN
         ! Ice collection of aerosol -- number sink
         nrate(1,1) = nrate(1,1) * aero(1,1,itrgt)%numc
         CALL rateDiag%AIcoll%accumulate(n=nrate(1,1))
      END IF
      
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

    SUBROUTINE accumulateSource(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,coll,source,trgtphase)
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
      INTEGER, INTENT(in) :: trgtphase   ! Phase identifier of the target bin, 1 aerosol, 2 clouds, 3 precip, 4 ice

      INTEGER :: ll,ii,jj
      REAL :: vrate(nspec,kbdim,klev)
      REAL :: dvol(nspec)
      vrate = 0.
      dvol = 0.
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)    
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + dvol(1:nspec)               
               vrate(1:nspec,ii,jj) = vrate(1:nspec,ii,jj) + dvol(1:nspec) ! For diagnostics
            END DO
         END DO
      END DO

      !Diagnostics
      IF (trgtphase == 3 .AND. coll(1,1,1)%phase == 2) THEN
         ! Accretion - volume source term
         vrate(1:nspec,1,1) = vrate(1:nspec,1,1) * precp(1,1,itrgt)%numc
         CALL rateDiag%Accretion%Accumulate(v=vrate(1:nspec,1,1))
      ELSE IF (trgtphase == 2 .AND. coll(1,1,1)%phase == 1) THEN
         ! Cloud collection of aerosol - volume source term
         vrate(1:nspec,1,1) = vrate(1:nspec,1,1) * cloud(1,1,itrgt)%numc
         CALL rateDiag%ACcoll%Accumulate(v=vrate(1:nspec,1,1))
      ELSE IF (trgtphase == 3 .AND. coll(1,1,1)%phase == 1) THEN
         ! Precipitation collection of aerosol -- volume source term
         vrate(1:nspec,1,1) = vrate(1:nspec,1,1) * precp(1,1,itrgt)%numc
         CALL rateDiag%APcoll%Accumulate(v=vrate(1:nspec,1,1))
      ELSE IF (trgtphase == 4 .AND. coll(1,1,1)%phase == 1) THEN
         ! Ice collection of aerosol -- volume source term
         vrate(1:nspec,1,1) = vrate(1:nspec,1,1) * ice(1,1,itrgt)%numc
         CALL rateDiag%AIcoll%Accumulate(v=vrate(1:nspec,1,1))
      END IF
      
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
    SUBROUTINE accumulateSourcePhaseChange(kbdim,klev,nbtrgt,nbcoll,nspec,iice,iwa,itrgt,istr,iend,  &
                                           rhotrgt,rhocoll,zcc,coll,source)
      !
      ! The direct method, where the "larger" particle category collects the "smaller" category.
      ! The target refers always to the "larger", collector category.
      ! 
      ! This subroutine takes into account the change in water density due to freezing upon collection
      ! 
      INTEGER,INTENT(in) :: kbdim,klev
      INTEGER,INTENT(in) :: nbtrgt, nbcoll   ! Number of bins in the target and collectee categories
      INTEGER,INTENT(in) :: nspec           ! Number of dry checmical compounds, including both unrimed and rimed ice 
      INTEGER,INTENT(in) :: iice,iwa            ! Index of the formed ice type, index of liquid water
      INTEGER,INTENT(in) :: itrgt,istr,iend  ! Index of the target bin, start and end indices of the collectee bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbcoll,nbtrgt)
      REAL, INTENT(in)   :: rhotrgt,rhocoll      ! Water densities for the target and collected categories 
                                                 ! (typically liquid and frozen, respectively).
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll) ! Collected particle properties
      REAL, INTENT(inout) :: source(nspec,kbdim,klev)

      INTEGER :: ll,ii,jj
      INTEGER :: ndry

      ndry = spec%getNSpec(type="dry")

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:ndry,ii,jj) = source(1:ndry,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:ndry)
               source(iice,ii,jj) = source(iice,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*rhocoll/rhotrgt 
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSourcePhaseChange

    ! ------------------------------------------

    SUBROUTINE accumulateSourceIce(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,coll,source)
      !
      ! This subroutine is for collection between different sized ice bins. It's similar to accumulateSourceRime
      ! in every other aspect, except the rime contribution has to be taken from its dedicated field. All the source
      ! terms regardles of the ice type are scaled to bulk pristine ice density for an easy conversion to mass in
      ! applyCoagIce, so that the mass weighted densities can be obtained correctly.
      !
      INTEGER, INTENT(in) :: kbdim,klev
      INTEGER, INTENT(in) :: nbtrgt,nbcoll
      INTEGER, INTENT(in) :: nspec    ! should contain all compounds including unrimed and rimed ice
      INTEGER, INTENT(in) :: itrgt,istr,iend
      REAL, INTENT(in)    :: zcc(kbdim,klev,nbcoll,nbtrgt)   !  
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll)
      REAL, INTENT(inout) :: source(nspec,kbdim,klev) ! Source term for all compounds

      INTEGER :: ll,ii,jj

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSourceIce


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

                  ! Diagnostics
                  CALL rateDiag%Autoconversion%Accumulate(n=num_prc(ii,jj,trgt_prc),    &
                                                          v=vol_prc(1:nspec,ii,jj,trgt_prc))
                  
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
