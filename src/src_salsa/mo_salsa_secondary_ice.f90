MODULE mo_salsa_secondary_ice


  IMPLICIT NONE

  SAVE

  PRIVATE
  PUBLIC  :: rimesplintering, dropfracturing, nfrzn_df, mfrzn_df, nfrzn_rs, mfrzn_rs, &
             dlliq_df, dlice_rs, dlliq_rs
  
  ! Arrays to track the number and mass of frozen drops due to ice collection
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init. Bin dimensions will be (nprc,nice). Possible contribution by 
  ! cloud droplets will be put to the first bin or smth?? Should binnihs for rime splintering be
  ! something else?? Check later!
  REAL, ALLOCATABLE :: nfrzn_rs(:,:,:,:), mfrzn_rs(:,:,:,:)
  REAL, ALLOCATABLE :: nfrzn_df(:,:,:,:), mfrzn_df(:,:,:,:)

  ! Ice and liquid drop diameter limits for drop fracturing
  !REAL :: dlice_df = 1.e-3,    &  ! Max diameter for ice in drop fracturing
  !        dlliq_df = 80.e-6       ! Min diameter for liquid in drop fracturing
  REAL :: dlliq_df = 100.e-6          ! Min droplet diameter for drop fracturing. Ice is expected to be more massive
                                   ! than the freezing drop.

  REAL :: dlice_rs = 100.e-6,   &  ! Min diameter for ice in hallet-mossop
          dlliq_rs = 2000.e-6      ! Max diameter for liquid in Hallet mossop; check appropriate value!!

  
  
  CONTAINS

    SUBROUTINE rimesplintering(kbdim,kproma,klev,nspec,ptemp,ptstep)
      USE mo_salsa_types, ONLY : ice,cloud,precp, rateDiag
      USE mo_submctl, ONLY : nprc,nice, pi6, spec, icebins
      ! The Hallet-Mossop secondary ice production by rime splintering
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec  ! nspec should contain active compounds + rime
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      
      REAL, PARAMETER :: c1 = 3.6e8   ! Splinters generated per milligram of new rime
      REAL, PARAMETER :: tmax = 270.16, tmid = 268.16, tmin = 265.16
      REAL, PARAMETER :: Dsplint = 10.e-6  ! Assumed splinter diameter

      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice)
      
      INTEGER :: ii,jj,bb,cc
      REAL :: dN  ! Number of splintered rime
      REAL :: dV  ! Mass of splintered rime
      INTEGER :: iwa, iri,ndry
      LOGICAL :: lt13umgt25um(kbdim,klev)
      INTEGER :: splbin  ! Target bin for splinters
      REAL :: ddrop ! Freezing drop diameter
      
      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")
      ndry = spec%getNSpec(type="dry")
      
      ! Convert freezing rates to changes over timestep
      mfrzn_rs = mfrzn_rs * ptstep
      nfrzn_rs = nfrzn_rs * ptstep
      
      ! Mask for where there are suitable size droplets present, i.e. smaller than 13um and larger than 25um
      ! The coagulation should have already been applied to the liquid bins, but with short timestep this is
      ! unlikely to fully empty the bin, so this is probably ok...
!!      lt13umgt25um = .FALSE.
!!      DO jj = 1,klev
!!         DO ii = 1,kproma            
!!            lt13umgt25um(ii,jj) = &
!!                 ( ANY( cloud(ii,jj,:)%numc > cloud(ii,jj,:)%nlim .AND. cloud(ii,jj,:)%dwet < 13.e-6 ) .OR.    &
!!                   ANY( precp(ii,jj,:)%numc > precp(ii,jj,:)%nlim .AND. precp(ii,jj,:)%dwet < 13.e-6 ) ) .AND. &
!!                 ( ANY( cloud(ii,jj,:)%numc > cloud(ii,jj,:)%nlim .AND. cloud(ii,jj,:)%dwet > 25.e-6 ) .OR.    &
!!                   ANY( precp(ii,jj,:)%numc > precp(ii,jj,:)%nlim .AND. precp(ii,jj,:)%dwet > 25.e-6 ) )
!!         END DO
!!      END DO

      ! Assuming splinter diameter as Dsplint, find the corresponding ice bin
      splbin = MAX( COUNT(icebins < Dsplint), 1 )      

      ! Initialize arrays
      fragvolc = 0.; sinkvolc = 0.
      fragnumc = 0.
      
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               DO cc = 1,nprc

                  IF ( ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR.  &    ! Outside temperature range, see Keinert et al 2020
                       nfrzn_rs(ii,jj,cc,bb) < 1.e-6 .OR. SUM(ice(ii,jj,bb)%volc(:)) < 1.e-23 .OR. &
                       ice(ii,jj,bb)%numc < ice(ii,jj,bb)%nlim ) CYCLE ! no collection/empty bin


                  ! Diameter of the frozen drops on current ice bin
                  ddrop = (mfrzn_rs(ii,jj,cc,bb)/nfrzn_rs(ii,jj,cc,bb)/spec%rhowa/pi6)**(1./3.)

                  ! Require freezing drop diameter to be larger than 25um
                  IF (ddrop < 25.e-6) CYCLE

                  dN = 0.
                  ! Number of generated splinters. The maximum rate is at 268.16 K
                  IF ( ptemp(ii,jj) < tmax .AND. ptemp(ii,jj) >= tmid ) THEN                     
                     dN = c1 * mfrzn_rs(ii,jj,cc,bb) * (tmax - ptemp(ii,jj))/(tmax-tmid) ! linear slope across the temp range                     
                  ELSE IF ( ptemp(ii,jj) < tmid .AND. ptemp(ii,jj) >= tmin ) THEN                     
                     dN = c1 * mfrzn_rs(ii,jj,cc,bb) * (ptemp(ii,jj) - tmin)/(tmid-tmin)                     
                  END IF

                  ! This will assume that the splinters consist of frozen spheres 10 um in diameter.
                  dV = dN*pi6*Dsplint**3

                  ! Allocate the fragments to temporary ice bins 
                  fragnumc(ii,jj,splbin) = fragnumc(ii,jj,splbin) + dN
                     
                  fragvolc(ii,jj,splbin,1:nspec) = fragvolc(ii,jj,splbin,1:nspec) +     &
                       ice(ii,jj,bb)%volc(1:nspec)*MIN( dV/SUM(ice(ii,jj,bb)%volc(1:nspec)), 1. )                  

                  sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) +   &
                       ice(ii,jj,bb)%volc(1:nspec)*MIN( dV/SUM(ice(ii,jj,bb)%volc(1:nspec)), 1. )  
                
                  ! Secondary ice diagnostics
                  ice(ii,jj,1)%SIP_rmspl = ice(ii,jj,1)%SIP_rmspl + dN                  
                  CALL rateDiag%rmsplrate%Accumulate(n=dN/ptstep)

               END DO
            END DO
         END DO
      END DO
           
      ! Apply changes to bins
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
!!               IF (fragnumc(ii,jj,bb) < 0.) WRITE(*,*) 'fragnumc < 0'
!!               IF (ANY(fragvolc(ii,jj,bb,:) < 0.) ) WRITE(*,*) 'fragvolc < 0'
               
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc + fragnumc(ii,jj,bb)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) + fragvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) - sinkvolc(ii,jj,bb,1:nspec)
!!               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
!!                    WRITE(*,*) 'RIME SPLNT NEGA END', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
            END DO
         END DO
      END DO
      
      ! IMPORTANT: Reset the collection tracking arrays
      mfrzn_rs = 0.
      nfrzn_rs = 0.
      
   END SUBROUTINE rimesplintering

    ! -------

   SUBROUTINE dropfracturing(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
      USE mo_salsa_types, ONLY : ice, rateDiag
      USE mo_submctl, ONLY : nprc, nice, pi6, spec, icebins, lssipdropfrac
      !
      ! -------------------------------------------------------
      ! The drop fracturing SIP from Lawson et al. 2015
      ! the inputs include the mass and number of frozen drizzle in ice particle bins.
      ! It is assumed this mass accumulation to ice is present identically as it comes out of the
      ! coagulation routines, i.e. the bin redistribution should NOT be calculated between sec ice
      ! and coagulation. Doing so would result in the loss of the required information on drop freezing.
      ! 
      ! -----------------------------------------------
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec   ! nspec should contain active compounds + rime, i.e. "total"
      REAL, INTENT(in) :: ppres(kbdim,klev),ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      REAL, PARAMETER :: tmin = 248.15, tmax = 271.15  ! Check these so the conform with the different formulations!!!

      REAL :: Nnorm          ! Normalization factor for distributing fragments
      REAL :: ddmean         ! Mean diameter of frozen drops per ice bin
      REAL :: dN,dm          ! Total number and mass of fragments generated per ice bin
      REAL :: dNb(nice)      ! Number of fragments distributed to ice bins
      REAL :: dVb(nice)      ! Volume of fragments distributed to ice bins
      INTEGER :: cc,bb,bb1,ii,jj,iri,iwa, nimax, npmax
      REAL :: icediams(nice), icebw(nice)
      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice), sinknumc(kbdim,klev,nice)  ! Number to be added and removed
      REAL :: fragv_loc(nice,nspec)  !! Local fragment vol contributions per ice bin
      REAL :: fragn_loc(nice)       !! Local fragment num contributions per ice bin
      REAL :: v_i                    ! Volume of single ice particle in a bin 
      REAL :: sinkv(nspec)           ! sink volume for single collision
      REAL :: frconst                ! constraining fraction for limiting the mass sink to ragments
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

      ! REMOVE IF NO LONGER USING FIXED ICE LIMITS
      ! Index of the largest ice bin which is assumed to trigger drop fracturing
      !nimax = COUNT(icediams <= dlice_df)

      ! POISTA
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'DROP FRAC NEGA BEG', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
            END DO
         END DO
      END DO      
      ! --------------------

      
      ! Initialize arrays
      fragvolc = 0.; sinkvolc = 0.
      fragnumc = 0.; sinknumc = 0.
      sinkv = 0.

      ! NO MORE FIXED ICE LIMITS -> CHANGE NIMAX TO NICE
      DO bb = 1,nice      
         DO jj = 1,klev
            DO ii = 1,kproma
               fragv_loc = 0.
               fragn_loc = 0.
               DO cc = 1,nprc
                  IF ( ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR.  &    ! Outside temperature range, see Keinert et al 2020
                       nfrzn_df(ii,jj,cc,bb) < 1.e-6 .OR. SUM(ice(ii,jj,bb)%volc(:)) < 1.e-23 .OR. &
                       ice(ii,jj,bb)%numc < ice(ii,jj,bb)%nlim ) CYCLE ! no collection/empty bin
               
                  ! Diameter of the frozen drops on current ice bin
                  ddmean = (mfrzn_df(ii,jj,cc,bb)/nfrzn_df(ii,jj,cc,bb)/spec%rhowa/pi6)**(1./3.)

                  !Require the freezing drop diameter to be larger tha dlliq_df
                  IF ( ddmean < dlliq_df ) CYCLE  

                  ! POISTA
                  IF (ddmean < 20.e-6 .OR. ddmean > 1.e3) WRITE(*,*) 'ddmean error ',ddmean,bb,nimax,icediams(bb),icebw(bb)
                  IF (ddmean /= ddmean) WRITE(*,*) 'ddmean nan ',ddmean,bb,nimax,icediams(bb),icebw(bb)
                  ! -----------------
               
                  ! Ice bin index corresponding to the mean frozen drop diameter minus one; Fragments are distributed to ice bins 1:npmax
                  npmax = MAX(COUNT(icediams <= ddmean) - 1, 1) 
              
                  ! Calculate the number of fragments generated per freezing droplet for current bin
                  IF (lssipdropfrac%mode == 1) THEN
                     dN = df_lawson(ptemp(ii,jj),nfrzn_df(ii,jj,cc,bb),ddmean)
                  ELSE IF (lssipdropfrac%mode == 2) THEN
                     dN = df_phillips_simple(ptemp(ii,jj),nfrzn_df(ii,jj,cc,bb),ddmean)
                  ELSE IF (lssipdropfrac%mode == 3) THEN
                     ! Updated diameter needed here
                     CALL ice(ii,jj,bb)%updateDiameter(.TRUE.,type="all")
                     dN = df_phillips_full(nspec,ppres(ii,jj),ptemp(ii,jj),ddmean,ice(ii,jj,bb))
                  END IF
     
                  ! Assume the mass of fragments distributed evenly to ice bins 1:npmax (Lawson et al 2015).
                  ! For this, first distribute dN as d**-3.
                  dNb = 0.
                  dVb = 0.
                  dNb(1:npmax) = 1./(icediams(1:npmax)**3)              ! density function               
                  Nnorm = SUM(dNb(1:npmax)*icebw(1:npmax))              ! Normalization factor

                  dNb(1:npmax) = dN * dNb(1:npmax)*icebw(1:npmax)/Nnorm ! Distributed bin concentrations of fragments
                  dVb(1:npmax) = dNb(1:npmax) * pi6*icediams(1:npmax)**3  ! Determine the fragment mass based on the ice bin diameters                  
                  
                  ! Allocate the fragments to temporary ice bins 
                  DO bb1 = 1,npmax
                     fragn_loc(bb1) = fragn_loc(bb1) + dNb(bb1)
                     fragv_loc(bb1,1:nspec) = fragv_loc(bb1,1:nspec) +    &
                          ice(ii,jj,bb)%volc(1:nspec)*( dVb(bb1)/SUM(ice(ii,jj,bb)%volc(1:nspec)) )                                      
                  END DO

                  ! Sink of volume from current bin
                  sinkv(1:nspec) = ice(ii,jj,bb)%volc(1:nspec)* SUM(dVb(1:npmax))/SUM(ice(ii,jj,bb)%volc(1:nspec))                  
                  sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) + sinkv(1:nspec)

                  ! Volume of a single ice particle in current bin for calculating the number concentration sink. This does not necessarily 
                  ! provide an exact representation for the fracturing particle size, but works as a first approximation.
                  v_i  = mfrzn_df(ii,jj,cc,bb)/nfrzn_df(ii,jj,cc,bb)/spec%rhori  !SUM(ice(ii,jj,bb)%volc(1:nspec))/ice(ii,jj,bb)%numc
               
                  ! Sink of number concentration from current bin - assume that the volume of single ice crystal stays constant through the process
                  sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) + SUM( sinkv(1:nspec) ) / v_i
               
                  ! for diagnostics
                  ice(ii,jj,1:npmax)%SIP_drfr = ice(ii,jj,1:npmax)%SIP_drfr + dNb(1:npmax)
                  !if (dN > 1.) WRITE(*,*) 'hephep ', dN,ptstep,SUM(dNb)
                  CALL rateDiag%drfrrate%Accumulate(n=SUM(dNb)/ptstep)    ! miks tanne tulee 0??? NOTE: syotin vakioarvoa subroutinen alussa, se kylla toimi.
               END DO
               
               !! Safeguard: Allow the fragments to take up to 90% of the source ice bin mass
               IF ( SUM(sinkvolc(ii,jj,bb,1:nspec)) > 0.9 * SUM(ice(ii,jj,bb)%volc(1:nspec)) ) THEN
                  frconst = 0.9 * SUM(ice(ii,jj,bb)%volc(1:nspec)) / SUM(sinkvolc(ii,jj,bb,1:nspec))
                  fragv_loc = fragv_loc * frconst
                  fragn_loc = fragn_loc * frconst
                  sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) * frconst
                  sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) * frconst
               END IF

               sinknumc(ii,jj,bb) = MIN(sinknumc(ii,jj,bb), 0.9*ice(ii,jj,bb)%numc) !! Additional constrain because for some reason
                                                                                    !! this still failed in the last bin...
               
               fragnumc(ii,jj,:) = fragnumc(ii,jj,:) + fragn_loc(:)
               fragvolc(ii,jj,:,:) = fragvolc(ii,jj,:,:) + fragv_loc(:,:)

               ! POISTA
               IF ( SUM(sinkvolc(ii,jj,bb,:))/SUM(ice(ii,jj,bb)%volc(1:nspec)) > 1.)  &
                    WRITE(*,*) 'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS', &
                    SUM(sinkvolc(ii,jj,bb,:)), SUM(fragvolc(ii,jj,:,:)), SUM(ice(ii,jj,bb)%volc(1:nspec))
               
               IF ( SUM(sinkvolc(ii,jj,bb,:)) > 0.95*SUM(ice(ii,jj,bb)%volc(1:nspec)) )     &
                    WRITE(*,*)  'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS 2', & 
                    SUM(sinkvolc(ii,jj,bb,:)), SUM(fragvolc(ii,jj,:,:)), SUM(ice(ii,jj,bb)%volc(1:nspec))

               IF (0.95*ice(ii,jj,bb)%numc < sinknumc(ii,jj,bb)) &
                    WRITE(*,*) 'SEC ICE ERROR: NUMBER SINK EXCEEED BIN NUMBER',  &
                    ice(ii,jj,bb)%numc, sinknumc(ii,jj,bb), bb, SUM(fragnumc(ii,jj,:)) 
               ! ---------------------------------------


               
            END DO
         END DO
      END DO
          
      ! Apply changes to bins
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               ! POISTA
               IF (fragnumc(ii,jj,bb) < 0.) WRITE(*,*) 'fragnumc < 0'
               IF ( ANY(fragvolc(ii,jj,bb,:) < 0.) ) WRITE(*,*) 'fragvolc < 0'
               IF (fragnumc(ii,jj,bb) /= fragnumc(ii,jj,bb)) &
                    WRITE(*,*) 'fragnumc nan',bb,dlliq_df
               IF ( ANY(fragvolc(ii,jj,bb,:) /= fragvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'fragvolc nan ',bb,dlliq_df,fragvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) < 0. ) ) &
                    WRITE(*,*) 'sinkvolc nega ',bb,dlliq_df,sinkvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) /= sinkvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'sinkvolc nan ',  bb,dlliq_df,sinkvolc(ii,jj,bb,:)
               IF (fragnumc(ii,jj,bb) > 1.e5) WRITE(*,*) 'fragnumc > 1e5 ',bb,dlliq_df,fragnumc(ii,jj,bb),    &
                    (SUM(mfrzn_df(ii,jj,:,bb))/SUM(nfrzn_df(ii,jj,:,bb))/spec%rhowa/pi6)**(1./3.), &
                    SUM(nfrzn_df(ii,jj,:,bb)), ice(ii,jj,bb)%numc
               ! ---------------------
               
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc + fragnumc(ii,jj,bb)
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc - sinknumc(ii,jj,bb)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) + fragvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) - sinkvolc(ii,jj,bb,1:nspec)
               ! POISTA
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'DROP FRAC NEGA END', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
               ! ---------------------------
            END DO
         END DO
      END DO
      
      ! IMPORTANT: Reset the collection tracking arrays
      mfrzn_df = 0.
      nfrzn_df = 0.
      
    END SUBROUTINE dropfracturing

    ! -----

    REAL FUNCTION df_lawson(ptemp,nfrzn,ddmean)
      USE math_functions, ONLY : f_gauss
      ! ---------------------------------------------
      ! Lawson et al. 2015 drop fracturing rate
      !
      REAL, INTENT(in) :: ptemp
      REAL, INTENT(in) :: nfrzn, ddmean
      REAL, PARAMETER :: c1 = 2.5e-11, c2 = 0.2, cexp = 4., T0 = 258., Tsig = 10.
      REAL :: hT
      hT = f_gauss(ptemp,Tsig,T0)/f_gauss(T0,Tsig,T0)
      IF (hT > 1.0 .OR. hT < 1.e-8) WRITE(*,*) 'HT VAARIN ',hT 
      df_lawson = nfrzn * c2*hT * c1*(MIN(ddmean,3.e-3)*1.e6)**cexp ! c2*hT according to Sullivan et al. 2018
    END FUNCTION df_lawson

    ! -----

    REAL FUNCTION df_phillips_simple(ptemp,nfrzn,ddmean)
      ! ------------------------------------------------------------
      ! Simplified drop fracturing rate from Phillips et al 2018
      !
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

      df_phillips_simple = nfrzn * (c2 * hT * ddmean )           
    END FUNCTION df_phillips_simple

    ! ------------------------------------

   REAL FUNCTION df_phillips_full(nspec,ppres,ptemp,ddmean,pice)
      USE classSection, ONLY : Section
      USE mo_submctl, ONLY : spec,pi6

      INTEGER, INTENT(in) :: nspec    ! Should contain nwet + rime, i.e. "total"
      REAL, INTENT(in) :: ppres,ptemp,ddmean
      TYPE(Section), INTENT(in) :: pice

      REAL, PARAMETER :: dmin1=50.e-6, dmin2=150.e-6, Tmin = 267.15

      REAL :: mrim,mpri,ncice  ! rimed and unrimed bin ice mix rats, ice number concentration
      REAL :: mip,mdp          ! Masses of single ice crystal, single freezing drop
      REAL :: rhoip            ! Bin mean ice density
      REAL :: ddmeanx

      df_phillips_full = 0.

      mrim = pice%volc(nspec) * spec%rhori
      mpri = SUM(pice%volc(1:nspec-1)) * spec%rhoic ! Cutting a little corners here with the volc...
      ncice = pice%numc
      
      ! Single particle and drop masses
      mip = (mrim+mpri)/ncice
      mdp = spec%rhowa * pi6 * ddmean**3

      IF (mdp > mip) THEN
         !! Mode 1 drop fragmentation

         !! This will take care of the "step functions" in Eq1 @ Phillips et al 2018
         IF ( ddmean < dmin1 .AND. ptemp > Tmin ) RETURN 

         ddmeanx = MIN(ddmean,1.6)
      
         df_phillips_full = df_phillips_mode1(ddmeanx,ptemp)

      ELSE IF (mdp <= mip) THEN
         !! Mode 2 

         IF (ddmean < dmin2) RETURN
         ddmeanx = ddmean
         
         rhoip = ( mrim*spec%rhori + mpri*spec%rhoic ) / ( mrim + mpri )

         df_phillips_full = df_phillips_mode2(ppres,ptemp,ddmeanx,pice%dwet,pice%dnsp,mrim,mpri,ncice,rhoip,spec%rhowa)

      END IF

   END FUNCTION df_phillips_full


   REAL FUNCTION df_phillips_mode1(ddmean,ptemp)     
      ! ---------------------------------------------------------------------------
      ! Mode 1 (small ice, big drop) drop fracturing rate from Phillips et al 2018
      !
      REAL, INTENT(in) :: ddmean, ptemp  !! ddmean in m, ptemp in K      
      REAL :: T0, zeta, eta, beta
      REAL :: tc
      tc = ptemp-273.15      
      df_phillips_mode1 = 0.      

      T0 = ph_T0(ddmean)
      zeta = ph_zeta(ddmean)
      eta = ph_eta(ddmean)
      beta = ph_beta(ddmean)      
      df_phillips_mode1 = beta*tc + (zeta * eta**2) / &
                          ( (tc-T0)**2 + eta**2 )               

   END FUNCTION df_phillips_mode1

   REAL FUNCTION df_phillips_mode2(ppres,ptemp,ddmean,disph,dinsph,mrim,mpri,ncice,rhoip,rhowa)
      USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
      USE mo_particle_external_properties, ONLY : terminal_vel
      USE mo_submctl, ONLY : rd, pstand,pi6,pi,surfw0,cwa,alf
      ! ---------------------------------------------------------------------------
      ! Mode 2 (big ice, small drop) drop fracturing rate from Phillips et al 2018
      !
      REAL, INTENT(in) :: ptemp,ppres    ! ptemp in K
      REAL, INTENT(in) :: ddmean         ! Freezing drop diameter in m
      REAL, INTENT(in) :: disph, dinsph  ! Spherical quivalent and non-spherical (max) diameters of ice particles
      REAL, INTENT(in) :: mrim, mpri      ! rimed and unrimed ice bin mixing ratios
      REAL, INTENT(in) :: ncice          ! Ice bin number concentration
      REAL, INTENT(in) :: rhoip          ! bin mean ice density
      REAL, INTENT(in) :: rhowa          ! Water density

      TYPE(t_shape_coeffs) :: ishape     ! Ice shape coefficients
      REAL :: mip,mdp  ! Masses of single ice particle and the freezing drop
      REAL :: vti,vtd  ! Terminal velocities of ice and drop
      REAL :: rhoa     ! air density

      ! This is repeating a LOT of the stuff already done once in coagulation kernels,
      ! which is BS and sad... But can't do much about it currently.
      REAL :: visc             ! Viscosity of air
      REAL :: mfp, knud, beta  ! Mean free path, knudsen number and cunningham correction
      REAL :: K0, DE, fT, tc

      df_phillips_mode2 = 0.

      rhoa = ppres/(rd*ptemp)
      visc = (7.44523e-3*SQRT(ptemp**3))/(5093.*(ptemp+110.4)) ! viscosity of air [kg/(m s)]
      mfp = (1.656e-10*ptemp+1.828e-8)*pstand/ppres

      ! Get the terminal velocities
      ! Ice
      knud = 2.*mfp/dinsph
      beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))
      CALL getShapeCoefficients(ishape,mpri,mrim,ncice)
      vti = terminal_vel(disph,rhoip,rhoa,visc,beta,4,ishape,dinsph)

      ! The droplet
      knud = 2.*mfp/ddmean
      beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))      
      vtd = terminal_vel(ddmean,rhowa,rhoa,visc,beta,3)

      mip = (mrim+mpri)/ncice
      mdp = rhowa*pi6*ddmean**3
      K0 = 0.5 * (mip*mdp/(mdp + mip)) * (vtd - vti)**2
      DE = K0 / (surfw0*pi*ddmean**2)
      tc = ptemp-273.15
      fT = -cwa*tc/alf

      df_phillips_mode2 = 3.*MIN(4.*fT,1.) * (1.-fT) * MAX(DE-0.2,0.)

   END FUNCTION df_phillips_mode2

    REAL FUNCTION ph_beta(ddmean)
      ! Polynomial for beta in Phillips et al 2018
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = -0.1839, c2 = -0.2017, c3 = -0.0512
      REAL :: dx
      ph_beta = 0.
      IF (ddmean >= 0.4e-3 ) THEN
         dx = LOG( MIN(ddmean*1.e3, 1.6) )
         ph_beta = (c1*dx**2) + (c2*dx) + c3
      END IF
    END FUNCTION ph_beta
    
    REAL FUNCTION ph_zeta(ddmean)
      ! Polynomial for zeta in Phillips et al. 2018
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 2.4268, c2 = 3.3274, c3 = 2.0783, c4 = 1.2927
      REAL :: dx, logzeta
      dx = LOG( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      logzeta = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
      ph_zeta = 10.**logzeta
    END FUNCTION ph_zeta
    
    REAL FUNCTION ph_eta(ddmean)
      ! polynomial for eta in Phillips et al. 2018
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 0.1242, c2 = -0.2316, c3 = -0.9874, c4 = -0.0827
      REAL :: dx, logeta
      dx = LOG( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      logeta = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
      ph_eta = 10.**logeta
    END FUNCTION ph_eta

    REAL FUNCTION ph_T0(ddmean)
      ! Polynomia for T0 in Phillips et al. 2018
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = -1.3999, c2 = -5.3285, c3 = -3.9847, c4 = -15.0332
      REAL :: dx
      dx = LOG( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      ph_T0 = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4      
    END FUNCTION ph_T0
    

    
END MODULE mo_salsa_secondary_ice
