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
                         lssecice, &
                         spec
  USE classSection, ONLY : Section
  IMPLICIT NONE

  CONTAINS

    !
    ! Aerosol coagulation
    ! -----------------------
    SUBROUTINE coag_aero(kbdim,klev,nspec,ptstep) 
      USE mo_salsa_types, ONLY : zccaa, zccca, zccpa, zccia
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep

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
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in2b,index_b,zccaa,aero,zplusterm,1)
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
            CALL accumulateSource(kbdim,klev,nbins,nbins,nspec,kk,in1a,index_a,zccaa,aero,zplusterm,1)
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
    SUBROUTINE coag_cloud(kbdim,klev,nspec,ptstep)
      USE mo_salsa_types, ONLY : zcccc, zccca, zccpc, zccic
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep

      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)

      ! These are necessary only for coagulation based precip formation
      REAL :: zvolsink_slf(nspec,kbdim,klev)
      REAL :: zvol_prc(nspec,kbdim,klev,nprc)
      REAL :: znum_prc(kbdim,klev,nprc)
      REAL :: zINF_prc(kbdim,klev,nprc)
      ! ---

      INTEGER :: index_aero_a, index_aero_b
      INTEGER :: index_a, index_b
      INTEGER :: ii,jj,kk
      REAL :: fix_coag

      DO kk = ica%cur,fca%cur
         IF ( ALL(cloud(:,:,kk)%numc < cloud(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         zvolsink_slf(:,:,:) = 0.
         zvol_prc(:,:,:,:) = 0.
         znum_prc(:,:,:) = 0.
         zINF_prc(:,:,:) = 0.

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
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_b,fcb%cur,zcccc,cloud,zminusterm,2)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               fix_coag = (max(1. - ptstep*(0.5*zcccc(1,1,kk,kk)*cloud(1,1,kk)%numc**2 / &
                    (cloud(1,1,kk)%numc)),0.1))
               ! This handles also self collection at the same time
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,fix_coag,zcccc,              &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc,zINF_prc    )
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
               fix_coag =    &
                    max(1. - ptstep*sum(zcccc(1,1,kk,ica%cur:fca%cur)*cloud(1,1,ica%cur:fca%cur)%numc), 0.1)
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,fix_coag,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc,zINF_prc  )
               fix_coag =   &
                    max(1. - ptstep*sum(zcccc(1,1,kk,icb%cur:fcb%cur)*cloud(1,1,icb%cur:fcb%cur)%numc), 0.1)
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,fix_coag,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc,zINF_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,kk-1,zcccc,cloud,zplusterm,2)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,index_b,zcccc,cloud,zplusterm,2)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,ptstep,cloud,zplusterm,zminusterm,zminus_self)

         IF (lsauto%state .AND. lsauto%mode == 1) &
              CALL applyCoagPrecipFormation(kbdim,klev,nspec,kk,ptstep,zvolsink_slf,zvol_prc,znum_prc,zINF_prc)

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
         zINF_prc(:,:,:) = 0.
         
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
               CALL accumulateSink(kbdim,klev,ncld,ncld,kk,index_a,fca%cur,zcccc,cloud,zminusterm,2)
            END IF

            IF (lsauto%state .AND. lsauto%mode == 1) THEN
               fix_coag = (max(1. - ptstep*(0.5*zcccc(1,1,kk,kk)*cloud(1,1,kk)%numc**2 / &
                    (cloud(1,1,kk)%numc)),0.1))
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,kk,kk,fix_coag,zcccc,  &
                                              zplusterm,zminus_self,zvolsink_slf,zvol_prc,znum_prc,zINF_prc )
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
               fix_coag =    &
                    max(1. - ptstep*sum(zcccc(1,1,kk,icb%cur:fcb%cur)*cloud(1,1,icb%cur:fcb%cur)%numc), 0.1)
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,fix_coag,zcccc,    &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc,zINF_prc  )
               fix_coag =    &
                    max(1. - ptstep*sum(zcccc(1,1,kk,ica%cur:fca%cur)*cloud(1,1,ica%cur:fca%cur)%numc), 0.1)
               CALL accumulatePrecipFormation(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,fix_coag,zcccc, &
                                              zplusterm,zminusterm,zvolsink_slf,zvol_prc,znum_prc,zINF_prc  )
            ELSE
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,icb%cur,kk-1,zcccc,cloud,zplusterm,2)
               CALL accumulateSource(kbdim,klev,ncld,ncld,nspec,kk,ica%cur,index_a,zcccc,cloud,zplusterm,2)
            END IF

         END IF

         !-- Volume and number concentrations after coagulation update 
         CALL applyCoag(kbdim,klev,ncld,nspec,kk,ptstep,cloud,zplusterm,zminusterm,zminus_self)

         IF (lsauto%state .AND. lsauto%mode == 1) &
              CALL applyCoagPrecipFormation(kbdim,klev,nspec,kk,ptstep,zvolsink_slf,zvol_prc,znum_prc,zINF_prc)

      END DO

    END SUBROUTINE coag_cloud

    !
    ! Precipitation coagulation
    ! --------------------------
    !
    SUBROUTINE coag_precp(kbdim,klev,nspec,ptstep) 
      USE mo_salsa_types, ONLY: zccpp, zccpa, zccpc, zccip
      INTEGER, INTENT(in) :: kbdim,klev,nspec
      REAL, INTENT(in) :: ptstep
      
      INTEGER :: kk,ii,jj
      REAL :: zplusterm(nspec,kbdim,klev), zminusterm(kbdim,klev), zminus_self(kbdim,klev)
      REAL :: znum_slf(kbdim,klev), zvol_slf(nspec,kbdim,klev)

      DO kk = 1,nprc
         IF ( ALL(precp(:,:,kk)%numc < precp(:,:,kk)%nlim) ) CYCLE

         zminusterm(:,:) = 0.
         zminus_self(:,:) = 0.
         zplusterm(:,:,:) = 0.
         znum_slf(:,:) = 0.
         zvol_slf(:,:,:) = 0.
         
         ! Collection by larger precip and self collection
         IF (lscgpp) THEN
            IF ( kk < nprc ) THEN
               ! Collection by larger
               CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk+1,nprc,zccpp,precp,zminusterm,3)
               ! Self collection for all but the largest bin
               CALL precipSelfCoag(kbdim,klev,nprc,nspec,kk,zccpp,znum_slf,zvol_slf)
            ELSE
               ! Self collection for the largest bin -- standard approach
               CALL accumulateSink(kbdim,klev,nprc,nprc,kk,kk,kk,zccpp,precp,zminus_self,3,multp=0.5)
            END IF
         END IF

         ! Collection by ice
         IF (lscgip) &
              CALL accumulateSink(kbdim,klev,nprc,nice,kk,iia,fia,zccip,ice,zminusterm,3)

         ! Volume gained from collection of aerosol
         IF (lscgpa) &
              CALL accumulateSource(kbdim,klev,nprc,nbins,nspec,kk,in1a,fn2b,zccpa,aero,zplusterm,3)
         
         ! Volume gained from collection of cloud droplets
         IF (lscgpc) &
              CALL accumulateSource(kbdim,klev,nprc,ncld,nspec,kk,1,ncld,zccpc,cloud,zplusterm,3)
         
         ! Volume gained from smaller precp
         IF (lscgpp .AND. kk > 1) &
              CALL accumulateSource(kbdim,klev,nprc,nprc,nspec,kk,1,kk-1,zccpp,precp,zplusterm,3)
         
         !-- Volume and number concentrations after coagulation update 
         CALL applyCoagPrecp(kbdim,klev,nprc,nspec,kk,ptstep,precp,zplusterm,zminusterm,zminus_self,znum_slf,zvol_slf)

      END DO
      
    END SUBROUTINE coag_precp

    !
    ! Ice coagulation
    ! ------------------
    ! 
    SUBROUTINE coag_ice(kbdim,klev,nspec,ptstep)
      USE mo_salsa_types, ONLY : zccii, zccia, zccic, zccip
      INTEGER, INTENT(in) :: kbdim,klev,nspec   ! nspec should contain all compounds including unrimed and rimed ice!
      REAL, INTENT(in) :: ptstep
      
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

    SUBROUTINE applyCoagPrecp(kbdim,klev,nb,nspec,itrgt,ptstep,part,     &
                              source,sink,sink_self,num_slf,vol_slf      )
      INTEGER, INTENT(in) :: kbdim,klev,nb,nspec,itrgt
      REAL, INTENT(in)    :: ptstep
      TYPE(Section), INTENT(inout) :: part(kbdim,klev,nb)
      REAL, INTENT(in) :: source(nspec,kbdim,klev), sink(kbdim,klev),     &
                          sink_self(kbdim,klev), num_slf(kbdim,klev),     &
                          vol_slf(nspec,kbdim,klev)
      INTEGER :: ii,jj
      
      DO jj = 1,klev
         DO ii = 1,kbdim
            part(ii,jj,itrgt)%volc(1:nspec) =              &
                 ( part(ii,jj,itrgt)%volc(1:nspec) +       &
                   ptstep*source(1:nspec,ii,jj)*part(ii,jj,itrgt)%numc )  / &
                 ( 1. + ptstep*sink(ii,jj) )

            part(ii,jj,itrgt)%numc = part(ii,jj,itrgt)%numc / ( 1. + ptstep*(sink(ii,jj) + sink_self(ii,jj)) )            
         END DO
      END DO
      
      ! For all but the largest bin, force the products from self coagulation to the next bin for more realistic growth rate
      ! to rain drops
      IF ( itrgt < nb ) THEN
         DO jj = 1,klev
            DO ii = 1,kbdim
               part(ii,jj,itrgt+1)%volc(1:nspec) =          &
                    part(ii,jj,itrgt+1)%volc(1:nspec) + ptstep*vol_slf(1:nspec,ii,jj)
               part(ii,jj,itrgt+1)%numc = part(ii,jj,itrgt+1)%numc + ptstep*num_slf(ii,jj)

               part(ii,jj,itrgt)%volc(1:nspec) =            &
                    part(ii,jj,itrgt)%volc(1:nspec) - ptstep*vol_slf(1:nspec,ii,jj)
               part(ii,jj,itrgt)%numc = part(ii,jj,itrgt)%numc - ptstep*num_slf(ii,jj)                              
            END DO
         END DO
         
      END IF
           
    END SUBROUTINE applyCoagPrecp

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

  SUBROUTINE applyCoagPrecipFormation(kbdim,klev,nspec,itrgt,ptstep,volsink_slf,vol_prc,num_prc,INF_prc)
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
      REAL, INTENT(in)    :: INF_prc(kbdim,klev,nprc)
      
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

               IF ( precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim )  &
                    precp(ii,jj,cc)%INdef = ( precp(ii,jj,cc)%INdef * (precp(ii,jj,cc)%numc - ptstep*nrate) +   &
                                              INF_prc(ii,jj,cc) * ptstep*nrate ) / precp(ii,jj,cc)%numc

            END DO
         END DO
      END DO

    END SUBROUTINE applyCoagPrecipFormation


    ! -----------------------------------------------------------------

    SUBROUTINE accumulateSink(kbdim,klev,nbtrgt,nbcoll,itrgt,istr,iend,zcc,coll,sink,trgtphase,multp)
      USE mo_salsa_secondary_ice, ONLY : nfrzn_rs, nfrzn_df, dlliq_df, dlice_rs, dlliq_rs ! REMOVE dlice_df
      USE mo_submctl, ONLY : Eiagg_max, Eiagg_min, lssipdropfrac
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
      
      INTEGER :: ll,ii,jj,ix,ex
      REAL :: xx
      REAL :: nrate(kbdim,klev), nrate80(kbdim,klev), nrate_au80(kbdim,klev),  & ! For diagnostics
              nrate50(kbdim,klev), nrate_au50(kbdim,klev)
      REAL :: dnum, D
      REAL :: Eagg(kbdim,klev,nbcoll)   !! Aggregation efficiency, for now only used for ice-ice
      REAL :: rimfr

      REAL :: fix_coag(kbdim,klev)  ! correction factor for the direct forward diagnostic calculations, see if works
      
      INTEGER :: iri,iwa
      
      iwa = spec%getIndex('H2O')
      iri = spec%getIndex('rime')
      
      nrate = 0.
      nrate50 = 0.
      nrate80 = 0.
      nrate_au50 = 0.
      nrate_au80 = 0.
      dnum = 0.
      
      ! For self collection
      IF ( PRESENT(multp) ) THEN
         xx = multp
      ELSE
         xx = 1.0
      END IF

      ! If ice-ice collision, determine aggregation efficiency.
      ! Min and max values given from namelist. Assume the actual
      ! value to be inversely proportional to rime fraction in this
      ! range.



      Eagg = 1.
      IF ( coll(1,1,1)%phase == 4 .AND. trgtphase == 4 ) THEN
         DO ll = istr,iend
            DO jj = 1,klev
               DO ii = 1,kbdim
                  rimfr = MAX( coll(ii,jj,ll)%getRimeFraction(),    &
                               coll(ii,jj,itrgt)%getRimeFraction()  )
                  Eagg(ii,jj,ll) = Eiagg_max - rimfr * (Eiagg_max - Eiagg_min)
               END DO               
            END DO
         END DO
      END IF

      
      ! Collection sink term and rate diagnostics according to bin regime limits
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               dnum = Eagg(ii,jj,ll)*xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
               sink(ii,jj) = sink(ii,jj) + dnum
               nrate(ii,jj) = nrate(ii,jj) + dnum
            END DO
         END DO
      END DO
      
      ! Additional diagnostics for Accretion by drizzle D>80um and autoconversion of droplets past 80um
      IF (coll(1,1,1)%phase == 3) THEN    ! For autoconversion do only the case for growth of drizzle embryos,
                                          !cloud droplets already taken care of in sourcePrecipFormation
         ! collect cloud droplets or drizzle embryos
         IF ( trgtphase == 2 .OR.            &
              (trgtphase == 3 .AND.           &
               precpbins(MIN(itrgt,nprc)) < 80.e-6) ) THEN 
            ! Accretion loop
            ix = COUNT(precpbins < 80.e-6)+1 ! Collector bin loop to start from 80um
            DO ll = ix,iend
               DO jj = 1,klev
                  DO ii = 1,kbdim
                     dnum = Eagg(ii,jj,ll)*xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
                     nrate80(ii,jj) = nrate80(ii,jj) + dnum
                  END DO
               END DO
            END DO
            ! Autoconversion loop
            ex = MAX(istr,COUNT(precpbins < 80.e-6)) ! Collector bin loop to end to 80um (low limit for drizzle/rain )
                                                     ! istr should be ok, since itrgt is limited to < 80 um
            IF (trgtphase == 2) THEN
               D = cloud(1,1,itrgt)%dwet
            ELSE IF (trgtphase == 3) THEN
               D = precp(1,1,itrgt)%dwet
            END IF
            DO ll = istr,ex
               DO jj = 1,klev
                  DO ii = 1,kbdim
                     ! The resulting drop should have D >= 80um
                     IF ( D**3 + coll(ii,jj,ll)%dwet**3 >= (80.e-6)**3 ) THEN
                        dnum = Eagg(ii,jj,ll)*xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
                        nrate_au80(ii,jj) = nrate_au80(ii,jj) + dnum
                     END IF
                  END DO
               END DO                  
            END DO
            
         END IF
      END IF

      ! Additional diagnostics for Accretion by drizzle D>50um and autoconversion of droplets past 50um
      IF (coll(1,1,1)%phase == 3) THEN    ! For autoconversion do only the case for growth of drizzle embryos,
                                          !cloud droplets already taken care of in sourcePrecipFormation
         ! collect cloud droplets or drizzle embryos
         IF ( trgtphase == 2 .OR.            &
              (trgtphase == 3 .AND.           &
               precpbins(MIN(itrgt,nprc)) < 50.e-6) ) THEN 
            ! Accretion loop
            ix = COUNT(precpbins < 50.e-6)+1 ! Collector bin loop to start from 50um
            DO ll = ix,iend
               DO jj = 1,klev
                  DO ii = 1,kbdim
                     dnum = Eagg(ii,jj,ll)*xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
                     nrate50(ii,jj) = nrate50(ii,jj) + dnum
                  END DO
               END DO
            END DO
            ! Autoconversion loop
            ex = MAX(istr,COUNT(precpbins < 50.e-6)) ! Collector bin loop to end to 50um (low limit for drizzle/rain )
                                                     ! istr should be ok, since itrgt is limited to < 50 um
            ! NOTE THESE CONDITIONS (AND MANY OTHER...) HAVE TO BE REVISED IF THIS VERSION OF SALSA
            ! IS TO BE USED IN ANY OTHER THAN BOX MODEL CONFIG
            IF (trgtphase == 2) THEN
               D = cloud(1,1,itrgt)%dwet
            ELSE IF (trgtphase == 3) THEN
               D = precp(1,1,itrgt)%dwet
            END IF
            DO ll = istr,ex
               DO jj = 1,klev
                  DO ii = 1,kbdim
                     ! The resulting drop should have D >= 50um
                     IF ( D**3 + coll(ii,jj,ll)%dwet**3 >= (50.e-6)**3 ) THEN
                        dnum = Eagg(ii,jj,ll)*xx*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc
                        nrate_au50(ii,jj) = nrate_au50(ii,jj) + dnum
                     END IF
                  END DO
               END DO                  
            END DO
            
         END IF

      END IF

      ! Diagnostics for secondary ice parameterizations
      IF (coll(1,1,1)%phase == 4 .AND. lssecice%state) THEN

         ! CHECK IF THIS WORKS...
         fix_coag = 1.
         DO jj = 1,klev
            DO ii = 1,kbdim
               fix_coag(ii,jj) = MAX( 1. - SUM( zcc(ii,jj,itrgt,1:nice)*coll(ii,jj,1:nice)%numc ), 0.1 ) 
            END DO
         END DO
            
         DO ll = istr,iend
            DO jj = 1,klev
               DO ii = 1,kbdim
                  IF (trgtphase == 3 .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim ) THEN
                     
                     ! Drop fracturing: large drops collected by small ice; Possible for all drop fragment parameterizations
                     IF ( coll(ii,jj,ll)%dwet < precp(ii,jj,itrgt)%dwet .AND. precp(ii,jj,itrgt)%dwet > dlliq_df .AND.  &    ! SWITCH dlice_df -> precp%dwet
                          precp(ii,jj,itrgt)%numc > precp(ii,jj,itrgt)%nlim .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim) THEN                  
                        nfrzn_df(ii,jj,itrgt,ll) = nfrzn_df(ii,jj,itrgt,ll) +      &
                             Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*precp(ii,jj,itrgt)%numc*   &
                             fix_coag(ii,jj)
                     END IF
                     
                     ! Drop fracturing: drop collected by more massive ice; Possible for the full 2-mode Phillips et al
                     IF (lssipdropfrac%mode==3 .AND. &
                         coll(ii,jj,ll)%dwet >= precp(ii,jj,itrgt)%dwet .AND. precp(ii,jj,itrgt)%dwet > dlliq_df .AND.  &    ! SWITCH dlice_df -> precp%dwet
                         precp(ii,jj,itrgt)%numc > precp(ii,jj,itrgt)%nlim .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim) THEN
                        nfrzn_df(ii,jj,itrgt,ll) = nfrzn_df(ii,jj,itrgt,ll) +      &
                           Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*precp(ii,jj,itrgt)%numc*   &
                           fix_coag(ii,jj)
                     END IF

                     ! Hallet-Mossop with precp
                     IF ( coll(ii,jj,ll)%dwet > precp(ii,jj,itrgt)%dwet  .AND. coll(ii,jj,ll)%dwet > dlice_rs .AND.   &
                          precp(ii,jj,itrgt)%numc > precp(ii,jj,itrgt)%nlim ) THEN
                        nfrzn_rs(ii,jj,itrgt,ll) = nfrzn_rs(ii,jj,itrgt,ll) +      &
                             Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*precp(ii,jj,itrgt)%numc*   &
                             fix_coag(ii,jj)
                     END IF
                     
                  ELSE IF (trgtphase == 2 .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim ) THEN

                     ! Drop fracturing: large drops collected by small ice; Possible for all parameterizations
                     IF (coll(ii,jj,ll)%dwet < cloud(ii,jj,itrgt)%dwet .AND. cloud(ii,jj,itrgt)%dwet > dlliq_df .AND.  &  ! SWITCH dlice_df -> cloud%dwet
                         cloud(ii,jj,itrgt)%numc > cloud(ii,jj,itrgt)%nlim .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim) THEN
                        nfrzn_df(ii,jj,1,ll) = nfrzn_df(ii,jj,1,ll) +      &
                             Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc*   &
                             fix_coag(ii,jj)
                     END IF
                     
                     ! Drop fracturing: drop collected by more massive ice; Possible for the full 2-mode Phillips et al
                     IF (lssipdropfrac%mode==3 .AND. &
                         coll(ii,jj,ll)%dwet >= cloud(ii,jj,itrgt)%dwet .AND. cloud(ii,jj,itrgt)%dwet > dlliq_df .AND.  &    ! SWITCH dlice_df -> precp%dwet
                         cloud(ii,jj,itrgt)%numc > cloud(ii,jj,itrgt)%nlim .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim) THEN
                        nfrzn_df(ii,jj,1,ll) = nfrzn_df(ii,jj,1,ll) +      &
                           Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc*   &
                           fix_coag(ii,jj)
                     END IF

                     ! Hallet-Mossop with cloud droplets
                     IF ( coll(ii,jj,ll)%dwet > cloud(ii,jj,itrgt)%dwet  .AND. coll(ii,jj,ll)%dwet > dlice_rs .AND.  &
                         cloud(ii,jj,itrgt)%numc > cloud(ii,jj,itrgt)%nlim ) THEN
                        nfrzn_rs(ii,jj,1,ll) = nfrzn_rs(ii,jj,1,ll) +      &   !!! CHECK BINNING, THIS WILL BE WRONG!
                             Eagg(ii,jj,ll)*zcc(ii,jj,itrgt,ll)*coll(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc*   &
                             fix_coag(ii,jj)
                     END IF
                     
                  END IF                  
               END DO
            END DO
         END DO
      END IF

      
      ! Accumulate autoconversion and accretion diagnostics according to custom size limits
      IF (trgtphase == 2) THEN
         nrate80(1,1) = nrate80(1,1) * cloud(1,1,itrgt)%numc
         nrate50(1,1) = nrate50(1,1) * cloud(1,1,itrgt)%numc
         nrate_au80(1,1) = nrate_au80(1,1) * cloud(1,1,itrgt)%numc
         nrate_au50(1,1) = nrate_au50(1,1) * cloud(1,1,itrgt)%numc
         CALL rateDiag%Accretion80%Accumulate(n=nrate80(1,1))
         CALL rateDiag%Accretion50%Accumulate(n=nrate50(1,1))
         CALL rateDiag%Autoconversion80%Accumulate(n=nrate_au80(1,1))
         CALL rateDiag%Autoconversion50%Accumulate(n=nrate_au50(1,1))
      ELSE IF (trgtphase == 3) THEN
         nrate80(1,1) = nrate80(1,1) * precp(1,1,itrgt)%numc
         nrate50(1,1) = nrate50(1,1) * precp(1,1,itrgt)%numc
         nrate_au80(1,1) = nrate_au80(1,1) * precp(1,1,itrgt)%numc
         nrate_au50(1,1) = nrate_au50(1,1) * precp(1,1,itrgt)%numc
         CALL rateDiag%Accretion80%Accumulate(n=nrate80(1,1))
         CALL rateDiag%Accretion50%Accumulate(n=nrate50(1,1))
         CALL rateDiag%Autoconversion80%Accumulate(n=nrate_au80(1,1))
         CALL rateDiag%Autoconversion50%Accumulate(n=nrate_au50(1,1))
      END IF
            
      !  Accumulate diagnostics depending on bin regime limits  
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

      INTEGER :: ll,ii,jj,ex
      REAL :: vrate(nspec,kbdim,klev), vrate80(nspec,kbdim,klev), vrate_au80(nspec,kbdim,klev),  &
              vrate50(nspec,kbdim,klev), vrate_au50(nspec,kbdim,klev)
      REAL :: dvol(nspec)
      REAL :: D
      vrate = 0.
      vrate50 = 0
      vrate80 = 0.
      vrate_au50 = 0.
      vrate_au80 = 0.
      dvol = 0.

      ! Source term and rate diagnostics according to bin regime limits
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)
               vrate(1:nspec,ii,jj) = vrate(1:nspec,ii,jj) + dvol(1:nspec)
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) + dvol(1:nspec)               
            END DO
         END DO
      END DO

      ! Diagnostics
      ! Autoconversion for drops past 80um
      D = 0.
      ex = iend
      !IF (trgtphase == 2) THEN      ! This case is already accounted for in the sourcePrecipFormation
      !   D = cloud(1,1,itrgt)%dwet
      IF (trgtphase == 3) THEN
         D = precp(1,1,itrgt)%dwet
      END IF      
      IF (coll(1,1,1)%phase == 2) THEN
         ex = iend
      ELSE IF (coll(1,1,1)%phase == 3) THEN
         ex = MIN(iend,COUNT(precpbins < 80.e-6)) ! This should never actually trigger because of the
                                                  ! call structure
      END IF      
      IF ( D > 0. .AND. D < 80.e-6 .AND.               &
          (coll(1,1,1)%phase==2 .OR. coll(1,1,1)%phase==3) ) THEN
         DO ll = istr,ex
            DO jj = 1,klev
               DO ii = 1,kbdim
                  IF (D**3 + coll(ii,jj,ll)%dwet**3 > (80.e-6)**3) THEN
                     dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)    
                     vrate_au80(1:nspec,ii,jj) = vrate_au80(1:nspec,ii,jj) + dvol(1:nspec) ! For Collection and accretion terms
                  END IF
               END DO
            END DO
         END DO
      END IF

      ! Autoconversion for drops past 50 um
      D = 0.
      ex = iend
      !IF (trgtphase == 2) THEN      ! This case is already accounted for in the sourcePrecipFormation
      !   D = cloud(1,1,itrgt)%dwet
      IF (trgtphase == 3) THEN
         D = precp(1,1,itrgt)%dwet
      END IF      
      IF (coll(1,1,1)%phase == 2) THEN
         ex = iend
      ELSE IF (coll(1,1,1)%phase == 3) THEN
         ex = MIN(iend,COUNT(precpbins < 50.e-6)) ! This should never actually trigger because of the
                                                  ! call structure
      END IF      
      IF ( D > 0. .AND. D < 50.e-6 .AND.               &
          (coll(1,1,1)%phase==2 .OR. coll(1,1,1)%phase==3) ) THEN
         DO ll = istr,ex
            DO jj = 1,klev
               DO ii = 1,kbdim
                  IF (D**3 + coll(ii,jj,ll)%dwet**3 > (50.e-6)**3) THEN
                     dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)    
                     vrate_au50(1:nspec,ii,jj) = vrate_au50(1:nspec,ii,jj) + dvol(1:nspec) ! For Collection and accretion terms
                  END IF
               END DO
            END DO
         END DO
      END IF

      
      ! Additional diagnostics for Accretion by drizzle D > 80um
      IF (trgtphase == 3) THEN
         IF (precpbins(itrgt) > 80.e-6) THEN
            ex = iend
            IF (coll(1,1,1)%phase == 3) ex = MIN(iend,COUNT(precpbins < 80.e-6)) 
            IF ( coll(1,1,1)%phase == 2 .OR.   &
                 coll(1,1,1)%phase == 3        ) THEN

               DO ll = istr,ex
                  DO jj = 1,klev
                     DO ii = 1,kbdim
                        dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec) 
                        vrate80(1:nspec,ii,jj) = vrate80(1:nspec,ii,jj) + dvol(1:nspec)
                     END DO
                  END DO
               END DO
            
            END IF
         END IF
      END IF

      ! Additional diagnostics for Accretion by drizzle D > 50um
      IF (trgtphase == 3) THEN
         IF (precpbins(itrgt) > 50.e-6) THEN
            ex = iend
            IF (coll(1,1,1)%phase == 3) ex = MIN(iend,COUNT(precpbins < 50.e-6)) 
            IF ( coll(1,1,1)%phase == 2 .OR.   &
                 coll(1,1,1)%phase == 3        ) THEN

               DO ll = istr,ex
                  DO jj = 1,klev
                     DO ii = 1,kbdim
                        dvol(1:nspec) = zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec) 
                        vrate50(1:nspec,ii,jj) = vrate50(1:nspec,ii,jj) + dvol(1:nspec)
                     END DO
                  END DO
               END DO
            
            END IF
         END IF
      END IF
      
      ! Accumulate autoconversion and accretion diagnostics according to custom size limits
      IF (trgtphase == 3) THEN
         vrate80(1:nspec,1,1) = vrate80(1:nspec,1,1) * precp(1,1,itrgt)%numc
         vrate50(1:nspec,1,1) = vrate50(1:nspec,1,1) * precp(1,1,itrgt)%numc
         vrate_au80(1:nspec,1,1) = vrate_au80(1:nspec,1,1) * precp(1,1,itrgt)%numc
         vrate_au50(1:nspec,1,1) = vrate_au50(1:nspec,1,1) * precp(1,1,itrgt)%numc         
         CALL rateDiag%Accretion80%Accumulate(v=vrate80(1:nspec,1,1))
         CALL rateDiag%Accretion50%Accumulate(v=vrate50(1:nspec,1,1))
         CALL rateDiag%Autoconversion80%Accumulate(v=vrate_au80(1:nspec,1,1))
         CALL rateDiag%Autoconversion50%Accumulate(v=vrate_au50(1:nspec,1,1))
      END IF
         
      !Accumulate diagnostics depending on bin regime limits
      IF ( trgtphase == 3 .AND. coll(1,1,1)%phase == 2 ) THEN
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
      USE mo_salsa_secondary_ice, ONLY : mfrzn_df, mfrzn_rs, dlliq_df, dlice_rs, dlliq_rs ! REMOVE dlice_df
      USE mo_submctl, ONLY : lssipdropfrac
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
      INTEGER,INTENT(in) :: itrgt,istr,iend  ! Index of the target bin, start and end indices of the collected bins
      REAL, INTENT(in)   :: zcc(kbdim,klev,nbcoll,nbtrgt)
      REAL, INTENT(in)   :: rhotrgt,rhocoll      ! Water densities for the target and collected categories 
                                                 ! (typically frozen and liquid, respectively).
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll) ! Collected particle properties
      REAL, INTENT(inout) :: source(nspec,kbdim,klev)
      
      REAL :: fix_coag(kbdim,klev)    ! Correction term to limit the direct forward diagnostic towards the semi-implicit solution (?)

      INTEGER :: ll,ii,jj, ix
      INTEGER :: ndry, nb

      ndry = spec%getNSpec(type="dry")

      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:ndry,ii,jj) = source(1:ndry,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:ndry)
               source(iice,ii,jj) = source(iice,ii,jj) + zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*rhocoll/rhotrgt 
            END DO
         END DO
      END DO

      ! Diagnostics for secondary ice production
      IF (lssecice%state) THEN

         !! CHECK IF THIS ACTUALLY WORKS...
         fix_coag = 1.
         IF (coll(1,1,1)%phase == 2) THEN
            DO jj = 1,klev
               DO ii = 1,kbdim
                  fix_coag(ii,jj) = MAX( 1. - SUM( zcc(ii,jj,1:ncld,itrgt)*coll(ii,jj,1:ncld)%numc ), 0.1 ) 
               END DO
            END DO
         ELSE IF (coll(1,1,1)%phase == 3) THEN
            DO jj = 1,klev
               DO ii = 1,kbdim
                  fix_coag(ii,jj) = MAX( 1. - SUM( zcc(ii,jj,1:nprc,itrgt)*coll(ii,jj,1:nprc)%numc ), 0.1 ) 
               END DO
            END DO            
         END IF
         
         DO ll = istr,iend
            DO jj = 1,klev
               DO ii = 1,kbdim
                  IF ( ANY(coll(ii,jj,ll)%phase ==  [2,3])) THEN  
                     ! Drop fracturing: large drops collected by small ice; Possible for all parameterizations
                     IF ( ice(ii,jj,itrgt)%dwet < coll(ii,jj,ll)%dwet .AND. coll(ii,jj,ll)%dwet > dlliq_df .AND. &
                          ice(ii,jj,itrgt)%numc > ice(ii,jj,itrgt)%nlim .AND. coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim) THEN  ! SWITCH dlice_df -> coll%dwet
                        IF (coll(ii,jj,ll)%phase == 3) THEN
                           ix = ll
                        ELSE IF (coll(ii,jj,ll)%phase == 2) THEN
                           ix = 1            ! CHECK THIS BINNING, GET FROM ACTUAL DWET IN CASE OF CLOUD DROPS??
                        END IF
                        mfrzn_df(ii,jj,ix,itrgt) = mfrzn_df(ii,jj,ix,itrgt) +     &
                             zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*ice(ii,jj,itrgt)%numc*rhocoll * fix_coag(ii,jj)
                     END IF

                     ! Drop fracturing: drop collected by more massive ice; Possible for the full 2-mode Phillips et al
                     IF ( lssipdropfrac%mode==3 .AND. &
                          ice(ii,jj,itrgt)%dwet >= coll(ii,jj,ll)%dwet .AND. coll(ii,jj,ll)%dwet > dlliq_df .AND.  &   
                          coll(ii,jj,ll)%numc > coll(ii,jj,ll)%nlim .AND. ice(ii,jj,itrgt)%numc > ice(ii,jj,itrgt)%nlim ) THEN
                        IF (coll(ii,jj,ll)%phase == 3) THEN
                           ix = ll
                        ELSE IF (coll(ii,jj,ll)%phase == 2) THEN
                           ix = 1            ! CHECK THIS BINNING, GET FROM ACTUAL DWET IN CASE OF CLOUD DROPS??
                        END IF
                        mfrzn_df(ii,jj,ix,itrgt) = mfrzn_df(ii,jj,ix,itrgt) +     &
                              zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*ice(ii,jj,itrgt)%numc*rhocoll * fix_coag(ii,jj)
                     END IF

                     ! Hallet-Mossop: small drops collected by large ice
                     IF ( ice(ii,jj,itrgt)%dwet > dlice_rs .AND. coll(ii,jj,ll)%dwet < ice(ii,jj,itrgt)%dwet ) THEN
                        IF (coll(ii,jj,ll)%phase == 3) THEN
                           ix = ll
                        ELSE IF (coll(ii,jj,ll)%phase ==2 ) THEN
                           ix = 1
                        END IF
                        mfrzn_rs(ii,jj,ix,itrgt) = mfrzn_rs(ii,jj,ix,itrgt) +     &
                             zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(iwa)*ice(ii,jj,itrgt)%numc*rhocoll * fix_coag(ii,jj)
                     END IF
                                          
                  END IF

               END DO
            END DO
         END DO                  
      END IF
         
    END SUBROUTINE accumulateSourcePhaseChange

    ! ------------------------------------------

    SUBROUTINE accumulateSourceIce(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,zcc,coll,source)
      USE mo_submctl, ONLY : Eiagg_max, Eiagg_min
      !
      ! This subroutine is for collection between different sized ice bins.
      !
      INTEGER, INTENT(in) :: kbdim,klev
      INTEGER, INTENT(in) :: nbtrgt,nbcoll
      INTEGER, INTENT(in) :: nspec    ! should contain all compounds including unrimed and rimed ice
      INTEGER, INTENT(in) :: itrgt,istr,iend
      REAL, INTENT(in)    :: zcc(kbdim,klev,nbcoll,nbtrgt)   !  
      TYPE(Section), INTENT(in) :: coll(kbdim,klev,nbcoll)
      REAL, INTENT(inout) :: source(nspec,kbdim,klev) ! Source term for all compounds

      INTEGER :: ll,ii,jj

      REAL :: Eagg(kbdim,klev,nbcoll)
      REAL :: rimfr
      
     ! If ice-ice collision, determine aggregation efficiency.
      ! Min and max values given from namelist. Assume the actual
      ! value to be inversely proportional to rime fraction in this
      ! range.
      Eagg = 1.
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               rimfr = MAX( coll(ii,jj,ll)%getRimeFraction(),    &
                    coll(ii,jj,itrgt)%getRimeFraction()  )
               Eagg(ii,jj,ll) = Eiagg_max - rimfr * (Eiagg_max - Eiagg_min)
            END DO
         END DO
      END DO
      
      DO ll = istr,iend
         DO jj = 1,klev
            DO ii = 1,kbdim
               source(1:nspec,ii,jj) = source(1:nspec,ii,jj) +   &
                    Eagg(ii,jj,ll)*zcc(ii,jj,ll,itrgt)*coll(ii,jj,ll)%volc(1:nspec)
            END DO
         END DO
      END DO

    END SUBROUTINE accumulateSourceIce


    ! Category specific processes
    ! ---------------------------------
    SUBROUTINE accumulatePrecipFormation(kbdim,klev,nbtrgt,nbcoll,nspec,itrgt,istr,iend,fix_coag,zcc,     &
                                         source, sink, volsink_slf, vol_prc, num_prc, INF_prc             )
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
      REAL, INTENT(in)    :: fix_coag
      REAL, INTENT(in)      :: zcc(kbdim,klev,nbcoll,nbtrgt)    ! Collision kernels
      REAL, INTENT(inout)   :: source(nspec,kbdim,klev)         ! Regular coagulation source term
      REAL, INTENT(inout)   :: sink(kbdim,klev)           ! Regular coagulation sink term
      REAL, INTENT(inout)   :: volsink_slf(nspec,kbdim,klev)    ! Volume sink for cloud droplets upon self collection resulting in precipitation formation
      REAL, INTENT(inout)   :: vol_prc(nspec,kbdim,klev,nprc)   ! Volume source for precipitation as the result of collision coalescence
      REAL, INTENT(inout)   :: num_prc(kbdim,klev,nprc)         ! Number source for precipitation as the result of collision coalescence
      REAL, INTENT(inout)   :: INF_prc(kbdim,klev,nprc)         ! Contribution to the IN nucleated fraction from the source bins
                                                                ! (needed for ice nucleation if active)      
      REAL :: D_new,dnum,dvol(nspec)
      INTEGER :: trgt_prc
      INTEGER :: ii,jj,ll
      LOGICAL :: selfcoll

      !fix_coag = 1.
      
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
               dvol = 0.
               dnum = 0.
               
               IF ( D_new > precpbins(1)) THEN
                  
                  ! The resulting droplets are large -> put the collision coalescence contribution to precipitation
                  IF ( selfcoll ) THEN
                     ! Self collection
                     dvol(1:nspec) = cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc*zcc(ii,jj,ll,itrgt) * fix_coag
                     vol_prc(1:nspec,ii,jj,trgt_prc) = vol_prc(1:nspec,ii,jj,trgt_prc) + dvol(1:nspec)

                     dnum = (0.5*zcc(ii,jj,ll,itrgt)*cloud(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc) * fix_coag
                     num_prc(ii,jj,trgt_prc) = num_prc(ii,jj,trgt_prc) + dnum

                     ! Contribution to the IN nucleated fraction. Take the value from the self-coagulating cloud droplets
                     IF (num_prc(ii,jj,trgt_prc) > 0.) &  
                          INF_prc(ii,jj,trgt_prc) = ( INF_prc(ii,jj,trgt_prc) * (num_prc(ii,jj,trgt_prc)-dnum) +   &
                                                      cloud(ii,jj,ll)%INdef * dnum ) / num_prc(ii,jj,trgt_prc)
                     
                  ELSE
                     dvol(1:nspec) = ( cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc +         &
                          cloud(ii,jj,itrgt)%volc(1:nspec)*cloud(ii,jj,ll)%numc ) * zcc(ii,jj,ll,itrgt) * fix_coag
                     vol_prc(1:nspec,ii,jj,trgt_prc) = vol_prc(1:nspec,ii,jj,trgt_prc) + dvol(1:nspec)

                     dnum = (zcc(ii,jj,ll,itrgt)*cloud(ii,jj,ll)%numc*cloud(ii,jj,itrgt)%numc) * fix_coag
                     num_prc(ii,jj,trgt_prc) = num_prc(ii,jj,trgt_prc) + dnum

                     ! Contribution to the IN nucleated fraction. Take the value as the mean from the colliding droplet bins
                     IF (num_prc(ii,jj,trgt_prc) > 0.) &
                          INF_prc(ii,jj,trgt_prc) = ( INF_prc(ii,jj,trgt_prc) * (num_prc(ii,jj,trgt_prc)-dnum) +      &
                                                      0.5*(cloud(ii,jj,ll)%INdef+cloud(ii,jj,itrgt)%INdef) * dnum ) / &
                                                      num_prc(ii,jj,trgt_prc)


                     
                  END IF

                  ! Diagnostics
                  CALL rateDiag%Autoconversion%Accumulate(n=dnum,    &
                                                          v=dvol(1:nspec))

                  IF (D_new > 50.e-6) &
                       CALL rateDiag%Autoconversion50%Accumulate(n=dnum,    &
                                                                 v=dvol(1:nspec))
                  
                  IF (D_new > 80.e-6) &
                       CALL rateDiag%Autoconversion80%Accumulate(n=dnum,    &
                                                                 v=dvol(1:nspec))

                  
                  IF ( selfcoll ) THEN
                     ! Change in cloud droplet volume due to precip formation in self collection
                     volsink_slf(1:nspec,ii,jj) = volsink_slf(1:nspec,ii,jj) +   &
                          cloud(ii,jj,ll)%volc(1:nspec)*cloud(ii,jj,itrgt)%numc*zcc(ii,jj,ll,itrgt) * fix_coag 
                      ! Taahan on kaiken lisaksi sama kuin ylla laskettu???
                     
                     ! Contribution of precip formation due to self collection to the regular sink term
                     sink(ii,jj) = sink(ii,jj) + 0.5*zcc(ii,jj,ll,itrgt)*cloud(ii,jj,itrgt)%numc
                     
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

    ! -----------------------------------------------

    SUBROUTINE precipSelfCoag(kbdim,klev,nb,nspec,ibin,zcc,num_slf,vol_slf)
      INTEGER, INTENT(in) :: kbdim, klev
      INTEGER, INTENT(in) :: nb, nspec    ! Number of bins, number of species
      INTEGER, INTENT(in) :: ibin         ! the current bin for self-coagulation. The resulting particle should be moved to next bin
      REAL, INTENT(in)    :: zcc(kbdim,klev,nprc,nprc)
      REAL, INTENT(inout) :: vol_slf(nspec,kbdim,klev)
      REAL, INTENT(inout) :: num_slf(kbdim,klev)
      INTEGER :: ii,jj


      DO jj = 1,klev
         DO ii = 1, kbdim
            ! Self collection
            vol_slf(1:nspec,ii,jj) = vol_slf(1:nspec,ii,jj) +   &
                 precp(ii,jj,ibin)%volc(1:nspec)*precp(ii,jj,ibin)%numc*zcc(ii,jj,ibin,ibin) 
            
            num_slf(ii,jj) = num_slf(ii,jj) +                   &
                 (0.5*zcc(ii,jj,ibin,ibin)*precp(ii,jj,ibin)%numc*precp(ii,jj,ibin)%numc) 

         END DO
      END DO
      
      
    END SUBROUTINE precipSelfCoag
    
END MODULE mo_salsa_coagulation_processes
