MODULE mo_salsa_secondary_ice
  USE mo_salsa_types, ONLY : ice,cloud,precp, rateDiag
  USE mo_submctl, ONLY : nice, pi6, spec, icebins, lsicedropfrac
  IMPLICIT NONE

  SAVE

  PUBLIC  :: rimesplintering, dropfracturing
  PRIVATE :: df_lawson, df_phillips
  
  
  ! Arrays to track the number and mass of frozen drops due to ice collection
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init.
  REAL, ALLOCATABLE :: nfrzn_rs(:,:,:), mfrzn_rs(:,:,:)
  REAL, ALLOCATABLE :: nfrzn_df(:,:,:), mfrzn_df(:,:,:)

  ! Ice and liquid drop diameter limits for drop fracturing
  REAL :: dlice_df = 1.e-3,    &  ! Max diameter for ice in drop fracturing
          dlliq_df = 80.e-6       ! Min diameter for liquid in drop fracturing

  REAL :: dlice_rs = 50.e-6,   &  ! Min diameter for ice in hallet-mossop
          dlliq_rs = 500.e-6      ! Max diameter for liquid in Hallet mossop; check appropriate value!!

  
  ! Diameter limit for ice and liquid bins. For Hallet-Mossop, require ice bin diameter > dlimit
  ! and frozen drop diameter < dlimit. For drop fracturing, ice bin diameter < dlimit and frozen
  ! drop diameter > dlimit. This is similar to Qu et al. (2022) and avoids overlap between the
  ! two processes
  !REAL, PARAMETER :: dlimit = 100.e-6   
  
  CONTAINS

    SUBROUTINE rimesplintering(kbdim,kproma,klev,nspec,ptemp,ptstep)
      ! The Hallet-Mossop secondary ice production by rime splintering
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec  ! nspec should contain active compounds + rime
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      
      REAL, PARAMETER :: c1 = 3.6e8   ! Splinters generated per milligram of new rime
      REAL, PARAMETER :: tmax = 270.16, tmid = 268.16, tmin = 265.16
      REAL, PARAMETER :: Dsplint = 10.e-6  ! Assumed splinter diameter

      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice)
      
      INTEGER :: ii,jj,bb
      REAL :: dN  ! Number of splintered rime
      REAL :: dV  ! Mass of splintered rime
      INTEGER :: iwa, iri,ndry
      LOGICAL :: lt13umgt25um(kbdim,klev)
      INTEGER :: splbin  ! Target bin for splinters
      REAL :: frac
      
      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")
      ndry = spec%getNSpec(type="dry")
      
      ! Convert freezing rates to changes over timestep
      mfrzn_rs = mfrzn_rs * ptstep
      nfrzn_rs = nfrzn_rs * ptstep
      
      ! Mask for where there are suitable size droplets present, i.e. smaller than 13um and larger than 25um
      ! The coagulation should have already been applied to the liquid bins, but with short timestep this is
      ! unlikely to fully empty the bin, so this is probably ok...
      lt13umgt25um = .FALSE.
      DO jj = 1,klev
         DO ii = 1,kproma            
            lt13umgt25um(ii,jj) = &
                 ( ANY( cloud(ii,jj,:)%numc > cloud(ii,jj,:)%nlim .AND. cloud(ii,jj,:)%dwet < 13.e-6 ) .OR.    &
                   ANY( precp(ii,jj,:)%numc > precp(ii,jj,:)%nlim .AND. precp(ii,jj,:)%dwet < 13.e-6 ) ) .AND. &
                 ( ANY( cloud(ii,jj,:)%numc > cloud(ii,jj,:)%nlim .AND. cloud(ii,jj,:)%dwet > 25.e-6 ) .OR.    &
                   ANY( precp(ii,jj,:)%numc > precp(ii,jj,:)%nlim .AND. precp(ii,jj,:)%dwet > 25.e-6 ) )
         END DO
      END DO

      ! Assuming splinter diameter as Dsplint, find the corresponding ice bin
      splbin = MAX( COUNT(icebins < Dsplint), 1 )      

      ! Initialize arrays
      fragvolc = 0.; sinkvolc = 0.
      fragnumc = 0.
      
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma

               IF ( nfrzn_rs(ii,jj,bb) > 0. .AND. ice(ii,jj,bb)%volc(iri) > 0. .AND. lt13umgt25um(ii,jj) ) THEN
                  dN = 0.
                  ! Number of generated splinters. The maximum rate is at 268.16 K
                  IF ( ptemp(ii,jj) < tmax .AND. ptemp(ii,jj) >= tmid ) THEN                     
                     dN = c1 * mfrzn_rs(ii,jj,bb) * (tmax - ptemp(ii,jj))/(tmax-tmid) ! linear slope across the temp range                     
                  ELSE IF ( ptemp(ii,jj) < tmid .AND. ptemp(ii,jj) >= tmin ) THEN                     
                     dN = c1 * mfrzn_rs(ii,jj,bb) * (ptemp(ii,jj) - tmin)/(tmid-tmin)                     
                  END IF

                  ! This will assume that the splinters consist of frozen spheres 10 um in diameter.
                  dV = dN*pi6*Dsplint**3

                  IF (ice(ii,jj,bb)%volc(iri) < dV) WRITE(*,*) "SECICE HM FAIL " 

                  ! Allocate the fragments to temporary ice bins 
                  fragnumc(ii,jj,splbin) = fragnumc(ii,jj,splbin) + dN
                     
                  fragvolc(ii,jj,splbin,1:nspec) = fragvolc(ii,jj,splbin,1:nspec) +     &
                       ice(ii,jj,bb)%volc(1:nspec)*MIN( dV/SUM(ice(ii,jj,bb)%volc(1:nspec)), 1. )                  

                  sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) +   &
                       ice(ii,jj,bb)%volc(1:nspec)* MIN( dV/SUM(ice(ii,jj,bb)%volc(1:nspec)), 1. )  
                
                  ! Secondary ice diagnostics
                  ice(ii,jj,1)%SIP_rmspl = ice(ii,jj,1)%SIP_rmspl + dN                  
                  CALL rateDiag%rmsplrate%Accumulate(n=dN/ptstep)

               END IF
                  
            END DO
         END DO
      END DO

      ! Apply changes to bins
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               IF (fragnumc(ii,jj,bb) < 0.) WRITE(*,*) 'fragnumc < 0'
               IF (ANY(fragvolc(ii,jj,bb,:) < 0.) ) WRITE(*,*) 'fragvolc < 0'
               
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc + fragnumc(ii,jj,bb)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) + fragvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) - sinkvolc(ii,jj,bb,1:nspec)
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'RIME SPLNT NEGA END', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
            END DO
         END DO
      END DO
      
      ! Reset the tracking arrays
      mfrzn_rs = 0.
      nfrzn_rs = 0.
      
    END SUBROUTINE rimesplintering

    ! -------

    SUBROUTINE dropfracturing(kbdim,kproma,klev,nspec,ptemp,ptstep)
      !
      ! -------------------------------------------------------
      ! The drop fracturing SIP from Lawson et al. 2015
      ! the inputs include the mass and number of frozen drizzle in ice particle bins.
      ! It is assumed this mass accumulation to ice is present identically as it comes out of the
      ! coagulation routines, i.e. the bin redistribution should NOT be calculated between sec ice
      ! and coagulation. Doing so would result in the loss of the required information on drop freezing.
      ! 
      ! -----------------------------------------------
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec   ! nspec should contain active compounds + rime
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      REAL, PARAMETER :: tmin = 248.15, tmax = 271.15  ! Check these so the conform with the different formulations!!!

      REAL :: Nnorm          ! Normalization factor for distributing fragments
      REAL :: ddmean         ! Mean diameter of frozen drops per ice bin
      REAL :: dN,dm          ! Total number and mass of fragments generated per ice bin
      REAL :: dNb(nice)      ! Number of fragments distributed to ice bins
      REAL :: dVb(nice)      ! Volume of fragments distributed to ice bins
      INTEGER :: bb,bb1,ii,jj,iri,iwa, nimax, npmax
      REAL :: icediams(nice), icebw(nice)
      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice)
      REAL :: frac
      REAL, PARAMETER :: inf = HUGE(1.)
      
      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")

      ! Convert freezing rates to changes over timestep
      mfrzn_df = mfrzn_df * ptstep
      nfrzn_df = nfrzn_df * ptstep
      
      icediams = 0.  ! Need the ice bin center diameters and bin widths, is there a better way for this?
      icebw = 0.
      DO bb = 1,nice
         icediams(bb) = ice(1,1,bb)%dmid
         icebw(bb) = ( (ice(1,1,bb)%vhilim/pi6)**(1./3) - (ice(1,1,bb)%vlolim/pi6)**(1./3))
      END DO
      
      ! Index of the largest ice bin which is assumed to trigger drop fracturing
      nimax = COUNT(icediams <= dlice_df)

      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'DROP FRAC NEGA BEG', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
            END DO
         END DO
      END DO      
      
      ! Initialize arrays
      fragvolc = 0.; sinkvolc = 0.
      fragnumc = 0.
      
      DO bb = 1,nimax       ! Assume drop fracturing to take place from ice bins < dlice_df
         DO jj = 1,klev
            DO ii = 1,kproma
               
               IF ( ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR.  &    ! Outside temperature range, see Keinert et al 2020
                    nfrzn_df(ii,jj,bb) < 1.e-6 .OR. SUM(ice(ii,jj,bb)%volc(:)) < 1.e-23 .OR. &
                    ice(ii,jj,bb)%numc < ice(ii,jj,bb)%nlim ) CYCLE ! no collection/empty bin
               
               ! Mean diameter of the frozen drops on current ice bin
               ddmean = (mfrzn_df(ii,jj,bb)/nfrzn_df(ii,jj,bb)/spec%rhowa/pi6)**(1./3.)

               ! Limit to 5 mm and require the freezing drop diameter to be larger tha dlliq_df
               ddmean = MIN(ddmean,5.e-3)                 
               IF ( ddmean < dlliq_df ) CYCLE  
               
               ! Ice bin index corresponding to the mean frozen drop diameter minus one; Fragments are distributed to ice bins 1:npmax
               npmax = MAX(COUNT(icediams <= ddmean) - 1, 1) 
              
               ! Calculate the number of fragments generated per freezing droplet for current bin
               IF (lsicedropfrac%mode == 1) THEN
                  dN = df_lawson(nfrzn_df(ii,jj,bb),ddmean)
               ELSE IF (lsicedropfrac%mode == 2) THEN
                  dN = df_phillips(ptemp(ii,jj),nfrzn_df(ii,jj,bb),ddmean)
               END IF

               ! Impose maximum value for the number of fragments. This is somewhat arbitrary, but perhaps the easiest way to constrain the
               ! total mass of fragments, reducing the risk of overshooting the source bin concentration.
               dN = MIN( dN, 1.e4 )
               
               ! Assume the mass of fragments distributed evenly to ice bins 1:npmax (Lawson et al 2015).
               ! For this, first distribute dN as d**-3.
               dNb = 0.
               dVb = 0.
               dNb(1:npmax) = 1./(icediams(1:npmax)**3)              ! density function               
               Nnorm = SUM(dNb(1:npmax)*icebw(1:npmax))              ! Normalization factor

               dNb(1:npmax) = dN * dNb(1:npmax)*icebw(1:npmax)/Nnorm ! Distributed bin concentrations of fragments
               dNb = MERGE(0., dNb, dNb < ice(ii,jj,bb)%nlim)        ! Cutoff for spuriously small fragment numbers
               dVb(1:npmax) = dNb(1:npmax) * pi6*icediams(1:npmax)**3  ! Determine the fragment mass based on the ice bin diameters
                                             
               ! Allocate the fragments to temporary ice bins 
               DO bb1 = 1,npmax
                  fragnumc(ii,jj,bb1) = fragnumc(ii,jj,bb1) + dNb(bb1)

                  fragvolc(ii,jj,bb1,1:nspec) = fragvolc(ii,jj,bb1,1:nspec) +     &
                       ice(ii,jj,bb)%volc(1:nspec)*( dVb(bb1)/SUM(ice(ii,jj,bb)%volc(1:nspec)) )                  
               END DO

               sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) +   &
                    ice(ii,jj,bb)%volc(1:nspec)* MIN((SUM(dVb(1:npmax))/SUM(ice(ii,jj,bb)%volc(1:nspec))), 1.)  

               IF ( SUM(dVb(1:npmax))/SUM(ice(ii,jj,bb)%volc(1:nspec)) > 1.)  &
                    WRITE(*,*) 'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS', &
                    SUM(dVb(1:npmax)), SUM(ice(ii,jj,bb)%volc(1:nspec))

               IF ( SUM(dVb(1:npmax)) > SUM(ice(ii,jj,bb)%volc(1:nspec)) )     &
                    WRITE(*,*)  'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS 2', &
                    SUM(dVb(1:npmax)), SUM(ice(ii,jj,bb)%volc(1:nspec))
               
               ! for diagnostics
               ice(ii,jj,1:npmax)%SIP_drfr = ice(ii,jj,1:npmax)%SIP_drfr + dNb(1:npmax)
               CALL rateDiag%drfrrate%Accumulate(n=dN/ptstep)  
            END DO
         END DO
      END DO
          
      ! Apply changes to bins
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               IF (fragnumc(ii,jj,bb) < 0.) WRITE(*,*) 'fragnumc < 0'
               IF (ANY(fragvolc(ii,jj,bb,:) < 0.) ) WRITE(*,*) 'fragvolc < 0'
               
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc + fragnumc(ii,jj,bb)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) + fragvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) - sinkvolc(ii,jj,bb,1:nspec)
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'DROP FRAC NEGA END', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
            END DO
         END DO
      END DO
      
      ! Reset the tracking arrays
      mfrzn_df = 0.
      nfrzn_df = 0.
      
    END SUBROUTINE dropfracturing
    ! -----
    PURE REAL FUNCTION df_lawson(nfrzn,ddmean)
      REAL, INTENT(in) :: nfrzn, ddmean
      REAL, PARAMETER :: c1 = 2.5e-11, cexp = 4.
      df_lawson = nfrzn * (c1 * (ddmean*1.e6)**cexp)
    END FUNCTION df_lawson
    ! -----
    PURE REAL FUNCTION df_phillips(ptemp,nfrzn,ddmean)
      REAL, INTENT(in) :: ptemp, nfrzn, ddmean
      REAL :: hT
      REAL, PARAMETER :: tlims(5) = [-24., -20., -16., -10., -6.]+273.15 
      REAL, PARAMETER :: hTv(5) = [0.6, 0.24, 2.6, 0.43, 0.35]
      REAL, PARAMETER :: c2 = (4./9.)*1.e4
      
      hT = 0.
      IF ( ptemp >= tlims(1) .AND. ptemp < tlims(2) ) THEN
         hT = (htv(2)-htv(1)) * ((ptemp-tlims(1))/4.) + htv(1)
      ELSE IF ( ptemp >= tlims(2) .AND. ptemp < tlims(3) ) THEN
         hT = (htv(3)-htv(2)) * ((ptemp-tlims(3))/4.) + htv(3)
      ELSE IF ( ptemp >= tlims(3) .AND. ptemp < tlims(4) ) THEN
         hT = (htv(4)-htv(3)) * ((ptemp-tlims(4))/6.) + htv(4)
      ELSE IF ( ptemp >= tlims(4) .AND. ptemp < tlims(5) ) THEN
         hT = (htv(5)-htv(4)) * ((ptemp-tlims(5))/4.) + htv(5)
      END IF

      df_phillips = nfrzn * (c2 * hT * ddmean )           
    END FUNCTION df_phillips




    
END MODULE mo_salsa_secondary_ice
