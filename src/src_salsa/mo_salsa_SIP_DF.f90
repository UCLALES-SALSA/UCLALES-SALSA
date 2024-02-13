MODULE mo_salsa_SIP_DF


  IMPLICIT NONE

  SAVE

  PRIVATE
  PUBLIC  :: dropfracturing,nfrzn_df, mfrzn_df, dlliq_df
  
  ! Arrays to track the number and mass of frozen drops due to ice collection
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init. Bin dimensions will be (nprc,nice). Possible contribution by 
  ! cloud droplets will be put to the first bin or smth?? 
  
  REAL, ALLOCATABLE :: nfrzn_df(:,:,:,:), mfrzn_df(:,:,:,:)

  ! liquid drop diameter limits for drop fracturing
  REAL ::  dlliq_df = 100.e-6      ! Min droplet diameter for drop fracturing.
  ! In droplet-ice collisions freezing drops are always assumed to be
  ! smaller than ice particles except in the case of 
  ! SIP parameterization of Phillips-full in function df_phillips_full_total
  ! If changes from published values are needed go to line 323
    
  CONTAINS

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
      REAL :: dN             ! Total number of fragments generated per ice bin per drop bin
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
      REAL :: dNbig
      
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

      !
       DO bb = 1,nice      
         DO jj = 1,klev
            DO ii = 1,kproma
               fragv_loc = 0.
               fragn_loc = 0.
               DO cc = 1,nprc
                  IF ( ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR.  &    ! Outside temperature range, see Keinert et al 2020
                       nfrzn_df(ii,jj,cc,bb) < 1.e-12 .OR. SUM(ice(ii,jj,bb)%volc(:)) < 1.e-23 .OR. &
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
                     dN = df_sullivan(ptemp(ii,jj),nfrzn_df(ii,jj,cc,bb),ddmean)
                  ELSE IF (lssipdropfrac%mode == 3) THEN
                     dN = df_phillips_simple(ptemp(ii,jj),nfrzn_df(ii,jj,cc,bb),ddmean)
                  ELSE IF (lssipdropfrac%mode == 4) THEN
                     ! Updated diameter needed here
                     CALL ice(ii,jj,bb)%updateDiameter(.TRUE.,type="all")
                     dN = df_phillips_full_total(nspec,ppres(ii,jj),ptemp(ii,jj),ddmean,ice(ii,jj,bb)) &
                          * nfrzn_df(ii,jj,cc,bb)
                     ! Functions to calculate the number of big fragments are also coded here. We do not use them.
                     !dNbig = df_phillips_full_big(nspec,ppres(ii,jj),ptemp(ii,jj),ddmean,ice(ii,jj,bb)) * nfrzn_df(ii,jj,cc,bb)
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
               !IF ( SUM(sinkvolc(ii,jj,bb,:))/MAX(SUM(ice(ii,jj,bb)%volc(1:nspec)),1.e-23) > 1.)  &
               !     WRITE(*,*) 'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS', &
               !     SUM(sinkvolc(ii,jj,bb,:)), SUM(fragvolc(ii,jj,:,:)), SUM(ice(ii,jj,bb)%volc(1:nspec))
               
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
      ! Lawson, R. P., Woods, S., & Morrison, H. (2015).
      ! The Microphysics of Ice and Precipitation Development in Tropical Cumulus Clouds.
      ! Journal of the Atmospheric Sciences, 72(6), 2429–2445. https://doi.org/10.1175/JAS-D-14-0274.1
      !
      REAL, INTENT(in) :: ptemp
      REAL, INTENT(in) :: nfrzn, ddmean
      REAL, PARAMETER :: c1 = 2.5e-11, cexp = 4.
      df_lawson = nfrzn * c1*(MIN(ddmean,3.e-3)*1.e6)**cexp
    END FUNCTION df_lawson

    ! -----
    
    REAL FUNCTION df_sullivan(ptemp,nfrzn,ddmean)
      USE math_functions, ONLY : f_gauss
      ! ---------------------------------------------
      ! Sullivan, S. C., Hoose, C., Kiselev, A., Leisner, T., & Nenes, A. (2018).
      ! Initiation of secondary ice production in clouds.
      ! Atmospheric Chemistry and Physics, 18(3), 1593–1610. https://doi.org/10.5194/acp-18-1593-2018
      !
      REAL, INTENT(in) :: ptemp
      REAL, INTENT(in) :: nfrzn, ddmean
      REAL, PARAMETER :: c1 = 2.5e-11, c2 = 0.2, cexp = 4., T0 = 258., Tsig = 10.
      REAL :: hT
      hT = f_gauss(ptemp,Tsig,T0)/f_gauss(T0,Tsig,T0)
      IF (hT > 1.0 .OR. hT < 1.e-8) WRITE(*,*) 'HT VAARIN ',hT 
      df_sullivan = nfrzn * c2*hT * c1*(MIN(ddmean,3.e-3)*1.e6)**cexp ! c2*hT according to Sullivan et al. 2018
    END FUNCTION df_sullivan

    ! -----
    
    REAL FUNCTION df_phillips_simple(ptemp,nfrzn,ddmean)
      ! ------------------------------------------------------------
      ! Simplified drop fracturing rate from Phillips et al 2018
      ! Equation 16
      ! Phillips, V. T. J., Patade, S., Gutierrez, J., & Bansemer, A. (2018).
      ! Secondary Ice Production by Fragmentation of Freezing Drops: Formulation and Theory.
      ! Journal of the Atmospheric Sciences, 75(9), 3031–3070. https://doi.org/10.1175/JAS-D-17-0190.1
      !
      REAL, INTENT(in) :: ptemp, nfrzn, ddmean
      REAL :: hT
      REAL, PARAMETER :: tlims(5) = [-24., -20., -16., -10., -6.]+273.15 
      REAL, PARAMETER :: hTv(5) = [0.6, 0.24, 2.6, 0.43, 0.35]
      REAL, PARAMETER :: c2 = (4./9.)*1.e4 ! tensile strength of pure ice at 2GPa(-15degC) 
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

    REAL FUNCTION df_phillips_full_total(nspec,ppres,ptemp,ddmean,pice)
      ! --------------------------------------------------------------------
      ! Phillips, V. T. J., Patade, S., Gutierrez, J., & Bansemer, A. (2018).
      ! Secondary Ice Production by Fragmentation of Freezing Drops: Formulation and Theory.
      ! Journal of the Atmospheric Sciences, 75(9), 3031–3070. https://doi.org/10.1175/JAS-D-17-0190.1
      ! Mode 1: Equation 1
      ! Mode 2: Equation 7  
      !
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

      df_phillips_full_total = 0.

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

         ! Total number of fragments or secondary ice particles
         df_phillips_full_total = df_phillips_mode1_total(ddmeanx,ptemp)

       ELSE IF (mdp <= mip) THEN
         !! Mode 2 

         IF (ddmean < dmin2) RETURN
         ddmeanx = ddmean
         
         rhoip = ( mrim*spec%rhori + mpri*spec%rhoic ) / ( mrim + mpri )

         df_phillips_full_total = df_phillips_mode2(ppres,ptemp,ddmeanx,pice%dwet,pice%dnsp,mrim,mpri,ncice,rhoip,spec%rhowa)

      END IF

    END FUNCTION df_phillips_full_total


       REAL FUNCTION df_phillips_full_big(nspec,ppres,ptemp,ddmean,pice)
      ! --------------------------------------------------------------------
      ! Phillips, V. T. J., Patade, S., Gutierrez, J., & Bansemer, A. (2018).
      ! Secondary Ice Production by Fragmentation of Freezing Drops: Formulation and Theory.
      ! Journal of the Atmospheric Sciences, 75(9), 3031–3070. https://doi.org/10.1175/JAS-D-17-0190.1
      ! Mode 1: Equation 1
      ! Mode 2: Equation 7  
      !
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

      df_phillips_full_big = 0.

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

         ! Total number of fragments or secondary ice particles
         df_phillips_full_big = df_phillips_mode1_big(ddmeanx,ptemp)

       ELSE IF (mdp <= mip) THEN
         !! Mode 2 

         IF (ddmean < dmin2) RETURN

         df_phillips_full_big = 0.

      END IF

   END FUNCTION df_phillips_full_big


   REAL FUNCTION df_phillips_mode1_total(ddmean,ptemp)     
      ! ---------------------------------------------------------------------------
      ! Mode 1 (small ice, big drop) drop fracturing rate from Phillips et al 2018
      ! Equation 1 in Phillips et al. (2018)
      ! It gives the total number of fragments or secondary ice particles
     
      REAL, INTENT(in) :: ddmean, ptemp  !! ddmean in m, ptemp in K      
      REAL :: T0, zeta, eta, beta
      REAL :: tc
      tc = ptemp-273.15      
      df_phillips_mode1_total = 0.      

      T0 = ph_T0(ddmean)
      zeta = ph_zeta(ddmean)
      eta = ph_eta(ddmean)
      beta = ph_beta(ddmean)
     
      df_phillips_mode1_total = beta*tc + (zeta * eta**2) / &
                          ( (tc-T0)**2 + eta**2 )               

    END FUNCTION df_phillips_mode1_total

   REAL FUNCTION df_phillips_mode1_big(ddmean,ptemp)     
      ! ---------------------------------------------------------------------------
      ! Mode 1 (small ice, big drop) drop fracturing rate from Phillips et al 2018
      ! Equation 3 in Phillips et al. (2018)
      ! It gives the number of big fragments among secondary ice particles
     
      REAL, INTENT(in) :: ddmean, ptemp  !! ddmean in m, ptemp in K      
      REAL :: T0, zeta, eta
      REAL :: tc
      REAL :: df_phillips_full_total
      
      tc = ptemp-273.15      
      df_phillips_mode1_big = 0.      

      T0 = ph_T0B(ddmean)
      zeta = ph_zetaB(ddmean)
      eta = ph_etaB(ddmean)

      df_phillips_full_total = df_phillips_mode1_total(ddmean,ptemp)
      
      df_phillips_mode1_big = MIN((zeta * eta**2) / &
                          ( (tc-T0)**2 + eta**2 ), df_phillips_full_total)               

   END FUNCTION df_phillips_mode1_big
    

   REAL FUNCTION df_phillips_mode2(ppres,ptemp,ddmean,disph,dinsph,mrim,mpri,ncice,rhoip,rhowa)
      USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
      USE mo_particle_external_properties, ONLY : terminal_vel
      USE mo_submctl, ONLY : rd, pstand,pi6,pi,surfw0,cwa,alf
      ! ---------------------------------------------------------------------------
      ! Mode 2 (big ice, small drop) drop fracturing rate from Phillips et al 2018
      ! Equation 7 in Phillips et al. (2018) 
      !
      REAL, INTENT(in) :: ptemp,ppres    ! ptemp in K and ppress in Pa
      REAL, INTENT(in) :: ddmean         ! Freezing drop diameter in m
      REAL, INTENT(in) :: disph, dinsph  ! Spherical equivalent and non-spherical (max) diameters of ice particles
      REAL, INTENT(in) :: mrim, mpri     ! rimed and unrimed ice bin mixing ratios
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

      ! Droplet
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
      ! Polynomial expression for beta in Phillips et al 2018
      ! See Table 3
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = -0.1839, c2 = -0.2017, c3 = -0.0512
      REAL :: dx
      ph_beta = 0.
      IF (ddmean >= 0.4e-3 ) THEN
         dx = LOG10( MIN(ddmean*1.e3, 1.6) )
         ph_beta = (c1*dx**2) + (c2*dx) + c3
      END IF
    END FUNCTION ph_beta
    
    REAL FUNCTION ph_zeta(ddmean)
      ! Polynomial expression for zeta in Phillips et al. 2018
      ! See Table 3
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 2.4268, c2 = 3.3274, c3 = 2.0783, c4 = 1.2927
      REAL :: dx, logzeta
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      logzeta = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
      ph_zeta = 10.**logzeta
    END FUNCTION ph_zeta
    
    REAL FUNCTION ph_eta(ddmean)
      ! Polynomial expression for eta in Phillips et al. 2018
      ! See Table 3 
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 0.1242, c2 = -0.2316, c3 = -0.9874, c4 = -0.0827
      REAL :: dx, logeta
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      logeta = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
      ph_eta = 10.**logeta
    END FUNCTION ph_eta

    REAL FUNCTION ph_T0(ddmean)
      ! Polynomial expression for T0 in Phillips et al. 2018
      ! See Table 3
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = -1.3999, c2 = -5.3285, c3 = -3.9847, c4 = -15.0332
      REAL :: dx
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      ph_T0 = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4      
    END FUNCTION ph_T0

       
    REAL FUNCTION ph_zetaB(ddmean)
      ! Polynomial expression for zetaB in Phillips et al. 2018
      ! See Table 4
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = -0.4651, c2 = -1.1072, c3 = -0.4539, c4 = 0.5137
      REAL :: dx
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      ph_zetaB = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
    END FUNCTION ph_zetaB
    
    REAL FUNCTION ph_etaB(ddmean)
      ! Polynomial expression for etaB in Phillips et al. 2018
      ! See Table 4
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 28.5888, c2 = 49.8504, c3 = 22.4873, c4 = 8.0481
      REAL :: dx
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      ph_etaB= (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4
    END FUNCTION ph_etaB

    REAL FUNCTION ph_T0B(ddmean)
      ! Polynomial expression for T0B in Phillips et al. 2018
      ! See Table 4
      REAL, INTENT(in) :: ddmean
      REAL, PARAMETER :: c1 = 13.3588, c2 = 15.7432, c3 = -2.6543, c4 = -18.4875
      REAL :: dx
      dx = LOG10( MAX( MIN(ddmean*1.e3, 1.6), 0.06 ) )
      ph_T0B = (c1*dx**3) + (c2*dx**2) + (c3*dx) + c4      
    END FUNCTION ph_T0B
    

    
END MODULE mo_salsa_SIP_DF
