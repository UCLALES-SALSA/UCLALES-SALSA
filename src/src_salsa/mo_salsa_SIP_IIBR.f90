MODULE mo_salsa_SIP_IIBR


  IMPLICIT NONE

  SAVE

  PRIVATE
  PUBLIC  :: iceicecollbreak, nii_ibr, mii_ibr
  
  ! Arrays to track the number and mass of ice particles during ice collection (i.e. ice-ice collisions)
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init. Bin dimensions will be (nice,nice).
  
  REAL, ALLOCATABLE :: nii_ibr(:,:,:,:), mii_ibr(:,:,:,:)

  ! Ice and liquid drop diameter limits for ice-ice collisional breakup
  ! There are none
  
 
  CONTAINS
    
   SUBROUTINE iceicecollbreak(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
      USE mo_salsa_types, ONLY : ice, rateDiag
      USE mo_submctl, ONLY : nice, pi6, spec, icebins, lssecice, lssipicecollbreak
      USE classSection, ONLY : Section
      ! ---------------------------------------------------------------------------------------------
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec   ! nspec should contain active compounds + rime
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ppres(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      REAL, PARAMETER :: tmin = 248.15, tmax = 270.15  ! following Takahashi et al. (1995) observations 
      ! Check temperature limits  so the conform with the different formulations!!!
      ! Takahashi, T., Nagao, Y., & Kushiyama, Y. (1995).
      ! Possible High Ice Particle Production during Graupel–Graupel Collisions.
      ! Journal of Atmospheric Sciences, 52(24), 4523–4527.
      ! https://doi.org/10.1175/1520-0469(1995)052<4523:PHIPPD>2.0.CO;2

      REAL :: IMF            ! Ice multiplication factor or number of secondary ice particles per SIP event
      REAL :: dN,dm          ! Total number and mass of fragments generated per ice bin
      REAL :: dNb(nice)      ! Number of fragments distributed to ice bins
      REAL :: dVb(nice)      ! Volume of fragments distributed to ice bins
      REAL :: Nnorm          ! Normalization factor used during the mass distribution of secondary ice particles
      INTEGER :: cc,bb,bb1,ii,jj,npmax
      REAL :: icediams(nice), icebw(nice)
      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice), sinknumc(kbdim,klev,nice)  ! Number to be added and removed
      REAL :: fragv_loc(nice,nspec)  !! Local fragment vol contributions per ice bin
      REAL :: fragn_loc(nice)       !! Local fragment num contributions per ice bin
      REAL :: v_i                    ! Volume of single ice particle in a bin 
      REAL :: sinkv(nspec)           ! sink volume for single collision
      REAL :: frconst                ! constraining fraction for limiting the mass sink to fragments
      REAL, PARAMETER :: inf = HUGE(1.)
      
      REAL :: dinsphmin  ! Non-spherical diameter (maximum length) of the smaller ice particle in the colliding pair or most fragile
      REAL :: disphmin   ! Spherical equivalent diameter of the smaller ice particle in the colliding pair or most fragile

      ! These are for checking purposes
      ! REAL :: dMean(nice) ! Mean size of the fragment in each bin
      ! REAL :: icediamslo(nice), icediamshi(nice) ! Low and high limits in bin diameter
      
      ! Convert collision rates to changes over timestep  (ice-ice collision rates)
      mii_ibr = mii_ibr * ptstep ! mass gained from smaller ice particles
      nii_ibr = nii_ibr * ptstep ! collisions from the accumulateSink as collection by larger ice and self-coagulation

      
      icediams = 0.  ! Need the ice bin center diameters and bin widths, is there a better way for this?
      icebw = 0.
      DO bb = 1,nice
         icediams(bb) = ice(1,1,bb)%dmid
         icebw(bb) = ( (ice(1,1,bb)%vhilim/pi6)**(1./3) - (ice(1,1,bb)%vlolim/pi6)**(1./3))
         !icediamslo(bb) = (ice(1,1,bb)%vlolim/pi6)**(1./3) 
         !icediamshi(bb) = (ice(1,1,bb)%vhilim/pi6)**(1./3)
      END DO
 
      ! Initialize arrays
      fragvolc = 0.; sinkvolc = 0.
      fragnumc = 0.; sinknumc = 0.
      sinkv = 0.

      ! Cycling through all possible collisions --> smaller particle kk with bb , larger particles kk+1 with cc
      DO bb = 1,nice      ! smaller
         DO jj = 1,klev
            DO ii = 1,kproma
               fragv_loc = 0.
               fragn_loc = 0.
              DO cc = bb, nice ! larger particles
               
               IF((ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax) .OR. & ! Outside temperature range, see Takahashi et al. 1995 
                       (nii_ibr(ii,jj,cc,bb) <1.e-12) .OR. &
                       (mii_ibr(ii,jj,cc,bb) <1.e-15) .OR. &
                       (SUM(ice(ii,jj,bb)%volc(:)) < 1.e-15).OR.(ice(ii,jj,bb)%numc  < ice(ii,jj,bb)%nlim).OR. &
                       (SUM(ice(ii,jj,cc)%volc(:)) < 1.e-15).OR.(ice(ii,jj,cc)%numc <  ice(ii,jj,cc)%nlim)) CYCLE ! no collection/empty bin
                   
               ! If colliding particles have the same size, SIP can still occur
               ! Phillips el. 2017 was corrected by Sotiropoulou et al. 2021
               ! using an expression to account for underestimates of the collision energy in these cases
               ! Ice-ice collisions in the mo_salsa_coagulation_processes.f90 are also corrected for these cases

               ! Ice bin index corresponding to the most fragile particle (smaller) minus one  because fragments 
               ! are smaller and must be distributed to ice bins 1:npmax
               ! Ice multiplication factors are written in terms of the maximum length of the smaller particle
               ! We will assume that the maximum length is the spherical equivalent diameter obtained from the mass 
               ! estimated with the effective density 
               ! Effective ice density. This should take into account non-spherical shape as well as the bulk ice composition
               !dinsphmin  = ((mii_ibr(ii,jj,bb,cc) / nii_ibr(ii,jj,bb,cc) / ice(ii,jj,cc)%rhoeff) / pi6)**(1./3.)
               !write( *, *) 'Effective density ', ice(ii,jj,bb)%rhoeff ! It is zero at this point
               !write(*,*) 'Mean density', ice(ii,jj,bb)%rhomean
               !write(*,*) 'nii_ibr', nii_ibr(ii,jj,cc,bb)
               !write(*,*) 'mii_ibr', mii_ibr(ii,jj,cc,bb)

               ! Volume of a single ice particle in current bin for calculating the number concentration sink. This does not necessarily 
               ! provide an exact representation for the fracturing particle size, but works as a first approximation
               !  
               CALL ice(ii,jj,bb)%updateDiameter(.TRUE.,type="all") 
                       
               v_i  = mii_ibr(ii,jj,cc,bb) / nii_ibr(ii,jj,cc,bb) / ice(ii,jj,bb)%rhomean !to be consistent with dinsphmin calculation

               dinsphmin  = (v_i / pi6)**(1./3.)             
               
               ! dinsphmin is the spherical equivalent diameter for the smallest particle with rhomean
               !write(*,*) 'Vi ', v_i
               !write(*,*) 'Nonspherical diameter ', dinsphmin

               npmax = COUNT(icediams <= dinsphmin)-1              

               IF (npmax <= 2) CYCLE ! ice particles smaller than 4um cannot experience fragmentation to even smaller particles
               !write(*,*) 'Size bin below dinsphmin ', npmax
                  
               dN = 0.
            
               ! Calculate the ice multiplication factor or number of fragments generated per ice-ice collision event for current bin
               ! Imposing an upper limit for IMF equal to 100. as Sotiropoulou et al. 2021

               IF (lssipicecollbreak%mode == 1) THEN ! Sullivan et al 2018 based on Takahashi et al. 1995
                     IMF = imf_sullivan(ptemp(ii,jj))
		     IMF = min(IMF, 100.0)
                     dN  = IMF *nii_ibr(ii,jj,cc,bb)
               
               ELSE IF (lssipicecollbreak%mode == 2) THEN  ! Sotiropoulou et al 2021 based on Sullivan et al 2018                   
                     IMF = imf_sotiropoulou(ptemp(ii,jj),dinsphmin)
                     IMF = min(IMF, 100.0)
                     dN  = IMF *nii_ibr(ii,jj,cc,bb) 
               
               ELSE IF (lssipicecollbreak%mode == 3) THEN  ! Phillips et al 2017 corrected by Sotiropoulou et al 2020
                     ! disphmin is the spherical equivalent diameter for the smallest particle using rhomean
                     ! rhomean is the mean ice density for frozen particles. Takes into account only the bulk ice composition                     
                     disphmin  =  ((mii_ibr(ii,jj,cc,bb) / nii_ibr(ii,jj,cc,bb) / ice(ii,jj,bb)%rhomean) / pi6)**(1./3.)
                     IMF = imf_phillips(ppres(ii,jj), ptemp(ii,jj),ice(ii,jj,bb), ice(ii,jj,cc), dinsphmin,disphmin) 
                     IMF = min(IMF, 100.0)
                     dN  = IMF *nii_ibr(ii,jj,cc,bb)

               END IF

               !write(*,*) 'IMF ', IMF
               !write(*,*) 'dN ',  dN
               
               ! Assume the mass of fragments distributed evenly to ice bins 1:npmax (Lawson et al 2015).
               ! For this, first distribute dN as d**-3.
               dNb = 0.
               dVb = 0.
               dNb(1:npmax) = 1./(icediams(1:npmax)**3)              ! density function               
               Nnorm = SUM(dNb(1:npmax)*icebw(1:npmax))              ! Normalization factor

               dNb(1:npmax) = dN * dNb(1:npmax)*icebw(1:npmax)/Nnorm ! Distributed bin concentrations of fragments
               dVb(1:npmax) = dNb(1:npmax) * pi6*icediams(1:npmax)**3! Determine the fragment mass based on the ice bin diameters                          
               
               ! Allocate the fragments to temporary ice bins 
               DO bb1 = 1,npmax
                     fragn_loc(bb1) = fragn_loc(bb1) + dNb(bb1)
                     fragv_loc(bb1,1:nspec) = fragv_loc(bb1,1:nspec) +    &
                  ice(ii,jj,bb)%volc(1:nspec)*( dVb(bb1)/SUM(ice(ii,jj,bb)%volc(1:nspec)) )                                      
               END DO

              ! Sink of volume from current bin
               sinkv(1:nspec) = ice(ii,jj,bb)%volc(1:nspec)* SUM(dVb(1:npmax))/SUM(ice(ii,jj,bb)%volc(1:nspec))                  
               sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) + sinkv(1:nspec)
 
              ! Sink of number concentration from current bin - assume that the volume of single ice crystal stays constant through the process
              sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) + SUM( sinkv(1:nspec) ) / v_i

              ! for diagnostics
              ice(ii,jj,1:npmax)%SIP_iibr = ice(ii,jj,1:npmax)%SIP_iibr + dNb(1:npmax)
            
              CALL rateDiag%iibrrate%Accumulate(n=SUM(dNb)/ptstep)    ! miks tanne tulee 0??? NOTE: syotin vakioarvoa subroutinen alussa, se kylla toimi.
              
            END DO

	    
            fragnumc(ii,jj,:) = fragnumc(ii,jj,:) + fragn_loc(:)
            fragvolc(ii,jj,:,:) = fragvolc(ii,jj,:,:) + fragv_loc(:,:)

            ! POISTA           
            IF ( SUM(sinkvolc(ii,jj,bb,:)) > SUM(ice(ii,jj,bb)%volc(1:nspec)) )     &
                  WRITE(*,*)  'SIP-IIBR ERROR: FRAGMENT MASS EXCEEDS BIN MASS 2', & 
                  SUM(sinkvolc(ii,jj,bb,:)), SUM(fragvolc(ii,jj,:,:)), SUM(ice(ii,jj,bb)%volc(1:nspec))

            IF (ice(ii,jj,bb)%numc < sinknumc(ii,jj,bb)) THEN
                  WRITE(*,*) 'SIP-IIBR ERROR: NUMBER SINK EXCEEED BIN NUMBER',  &
                  ice(ii,jj,bb)%numc, sinknumc(ii,jj,bb), bb, SUM(fragnumc(ii,jj,:)) 
                  sinknumc(ii,jj,bb) =  ice(ii,jj,bb)%numc
                  sinkvolc(ii,jj,bb,1:nspec) = ice(ii,jj,bb)%volc(1:nspec)
            END IF
            ! ---------------------------------------       
             !! Safeguard: Allow the fragments to take up to 99% of the source ice bin mass
            IF ( SUM(sinkvolc(ii,jj,bb,1:nspec)) > SUM(ice(ii,jj,bb)%volc(1:nspec)) ) THEN
               frconst = 0.99* SUM(ice(ii,jj,bb)%volc(1:nspec)) / SUM(sinkvolc(ii,jj,bb,1:nspec))
               fragv_loc = fragv_loc * frconst
               fragn_loc = fragn_loc * frconst
               sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) * frconst
               sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) * frconst
               sinknumc(ii,jj,bb) = MIN(sinknumc(ii,jj,bb), 0.99*ice(ii,jj,bb)%numc)
            END IF

            ! Confirming that the possible issue was solved
            IF ( SUM(sinkvolc(ii,jj,bb,1:nspec)) > SUM(ice(ii,jj,bb)%volc(1:nspec)) ) THEN
               WRITE(*,*) 'SIP-IIBR: SUM(sinkvolc(ii,jj,bb,1:nspec)) > SUM(ice(ii,jj,bb)%volc(1:nspec))'
               WRITE(*,*) SUM(sinkvolc(ii,jj,bb,1:nspec)), SUM(ice(ii,jj,bb)%volc(1:nspec))
            END IF	
               
            ! ---------------------------------------
           END DO
      END DO
   END DO
          
      ! Apply changes to bins
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma
               ! POISTA
               IF (fragnumc(ii,jj,bb) < 0.) WRITE(*,*) 'SIP-IIBR fragnumc < 0'
               IF ( ANY(fragvolc(ii,jj,bb,:) < 0.) ) WRITE(*,*) 'SIP-IIBR fragvolc < 0'
               IF (fragnumc(ii,jj,bb) /= fragnumc(ii,jj,bb)) &
                    WRITE(*,*) 'SIP-IIBR fragnumc nan',bb
               IF ( ANY(fragvolc(ii,jj,bb,:) /= fragvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'SIP-IIBR fragvolc nan ',bb,fragvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) < 0. ) ) &
                    WRITE(*,*) 'SIP-IIBR sinkvolc nega ',bb,sinkvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) /= sinkvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'SIP-IIBR sinkvolc nan ',  bb,sinkvolc(ii,jj,bb,:)
               !IF (fragnumc(ii,jj,bb) > 1.e5) WRITE(*,*) 'SIP-IIBR fragnumc > 1e5 ',bb,cc, fragnumc(ii,jj,bb), &
               !     nii_ibr(ii,jj,cc,bb), mii_ibr(ii,jj,cc,bb), ice(ii,jj,bb)%numc, &
               !     ice(ii,jj,bb)%dwet , ice(ii,jj,bb)%dnsp
               ! ---------------------
               
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc + fragnumc(ii,jj,bb)
               ice(ii,jj,bb)%numc = ice(ii,jj,bb)%numc - sinknumc(ii,jj,bb)
               ice(ii,jj,bb)%numc = MAX(0.,ice(ii,jj,bb)%numc)
               
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) + fragvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = ice(ii,jj,bb)%volc(1:nspec) - sinkvolc(ii,jj,bb,1:nspec)
               ice(ii,jj,bb)%volc(1:nspec) = MAX(0., ice(ii,jj,bb)%volc(1:nspec))
               
               ! POISTA
               IF ( ANY(ice(ii,jj,bb)%volc(1:nspec) < 0.) )  &
                    WRITE(*,*) 'DROP FRAC NEGA END', SUM(ice(ii,jj,bb)%volc(1:nspec)), ice(ii,jj,bb)%numc, bb
               ! ---------------------------
            END DO
         END DO
      END DO
      
      ! IMPORTANT: Reset the collection tracking arrays
      mii_ibr = 0.
      nii_ibr = 0.
      
    END SUBROUTINE iceicecollbreak
    
    ! ----------------------------------------------------------------------------------------------------------------------
    PURE REAL FUNCTION imf_sullivan(ptemp)
      
      ! Sullivan, S. C., Hoose, C., Kiselev, A., Leisner, T., & Nenes, A. (2018).
      ! Initiation of secondary ice production in clouds.
      ! Atmospheric Chemistry and Physics, 18(3), 1593–1610. https://doi.org/10.5194/acp-18-1593-2018
      
      REAL, INTENT(in) :: ptemp
      REAL, PARAMETER :: fbr = 280., Tmin = 252.
      
      ! fbr : leading coefficient of the fragment number generated per collision based upon data from Takahashi et al. (1995)
      ! Tmin: minimal temperature for ice-ice collision breakup to occur
      ! imf_sullivan: ice multiplication factor or number of secondary ice particles generated by ice-ice collision event

      ! IMF becomes negative below Tmin = 252K in this parameterization, then is set to zero
      
      imf_sullivan = max(0.,fbr*(ptemp-Tmin)**1.2*exp(-(ptemp-Tmin)/5.)) !Table 1 in Sullivan et al. (2018)
      
    END FUNCTION imf_sullivan
    ! -----

    ! ------------------------------------------------------------------------------------------------------------------------
    PURE REAL FUNCTION imf_sotiropoulou(ptemp,dinsphmin)
      ! Sotiropoulou, G., Ickes, L., Nenes, A., & Ekman, A. M. L. (2021).
      ! Ice multiplication from ice–ice collisions in the high Arctic: sensitivity to ice habit, rimed fraction,
      ! ice type and uncertainties in the numerical description of the process.
      ! Atmospheric Chemistry and Physics, 21(12), 9741–9760. https://doi.org/10.5194/acp-21-9741-2021
      
      REAL, INTENT(in) :: ptemp
      REAL, INTENT(in) :: dinsphmin  ! Non-spherical diameter (maximum length) of the smaller ice particle in the colliding pair or most fragile
      
      REAL, PARAMETER :: fbr = 280., Tmin = 252., D0= 0.02
      
      ! fbr : leading coefficient of the fragment number generated per collision based upon data from Takahashi et al. (1995)
      ! Tmin: minimal temperature for ice-ice collision breakup to occur
      ! imf_sotiropoulou: ice multiplication factor or number of secondary ice particles generated by ice-ice collision
      ! D: D (in meters) is the size of the ice particle that undergoes fracturing, or smaller ice particle in the ice-ice colliding pair
      ! D0:  0.02 m is the size of hail balls used by Takahashi et al. (1995)
      ! D is assumed to be the nonspherical diameter or maximum length of the smaller ice particle in the ice-ice colliding pair

      ! question: Could we avoid the temperature comparison using the max. IMF becomes negative below Tmin = 252K
      
      imf_sotiropoulou = MAX(0., fbr*(ptemp-Tmin)**1.2*exp(-(ptemp-Tmin)/5.)*dinsphmin/D0)
      
    END FUNCTION imf_sotiropoulou
    
    ! ------------------------------------------------------------------------------------------------------------------------

    REAL FUNCTION imf_phillips(ppres, ptemp,icelarge, icesmall, dinsphmin, disphmin) 
      ! Phillips, V. T. J., Yano, J.-I., & Khain, A. (2017).
      ! Ice Multiplication by Breakup in Ice–Ice Collisions. Part I: Theoretical Formulation.
      ! Journal of the Atmospheric Sciences, 74(6), 1705–1719. https://doi.org/10.1175/JAS-D-
      !  All equations reported in Table 1
            
      USE mo_submctl, ONLY : pi, spec
      USE classSection, ONLY : Section
      
      REAL, INTENT(in) :: ppres, ptemp
      
      TYPE(Section), INTENT(in) :: icelarge
      TYPE(Section), INTENT(in) :: icesmall !ice particle that undergoes fracturing, or smaller ice particle in the ice-ice colliding pair

      REAL :: dinsphmin ! Maximum length of the ice particle that undergoes fracturing, or smaller ice particle in the ice-ice colliding pair
      REAL :: disphmin  ! Spherical-equivalent diameter of the smaller ice particle (most fragile)  
 
      REAL :: K0        ! Kinetic collision energy of the ice-ice pair
      REAL :: mrim,mpri ! rimed and unrimed bin ice mass mix rats
      REAL :: rimfrac   ! fraction of rimed ice in the size bin of the smaller ice particle in the colliding pair      
      REAL :: alpha     ! Equivalent-spherical surface area of the colliding particle with the smaller maximum dimension in 1/m2 

       
     
      REAL, PARAMETER :: psi= 3.5e-3 ! fraction correcting the field observations by Vardiman for sublimation
      REAL, PARAMETER :: dmin=5.e-4, dmax=5.e-3   ! Limiting diameters for collisions type I graupel-hail
      REAL, PARAMETER :: Tmin=256.15, Tmax=261.15 ! Limiting temperatures for collisions type II hail-hail
      REAL, PARAMETER :: T0=258.15   ! T0 = -15 celsius Minimal temperature for ice-ice collision breakup to occur
      REAL :: C=0.                ! parameters(1) C, Asperity-fragility coefficient in 1/J
      REAL :: g=0.                ! parameters(2) gamma, Exponent in scheme for breakup (Eq.13), dimensionless
      REAL :: Nmax=0.             ! parameters(3) Nmax, Maximum number of fragments per ice-ice collision
      REAL :: Am=0.               ! parameters(4) Am, Measure of number density of breakable asperities in region of contact in 1/m^2
      REAL :: a0=0.               ! Maximum of Am or number density of breakable asperities in region of contact          
      INTEGER :: iwa,iri

      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")
      
      ! Getting the rimed fraction (by mass) of the more fragile particle in the colliding pair
      mrim = icesmall%volc(iri) * spec%rhori
      mpri = SUM(icesmall%volc(1:iwa)) * spec%rhoic ! Cutting a little corners here with the volc...
      rimfrac = mrim / (mrim + mpri)
      
      !  Getting the model parameters needed to calculate Am or number density of breakable asperities in region of contact
      !  params = getPhillipsparameters(ptemp, dinsphmin,rimfrac)
      ! From Phillips et al. (2017): when values are outside the valid ranges for D and rimfrac
      ! the inputs to the scheme are set to the nearest limit of the range
      
      IF (dinsphmin < dmin) THEN
         dinsphmin = dmin
      ELSE IF (dinsphmin > dmax) THEN
         dinsphmin = dmax
      END IF
      
      ! Calculating Am or number density of breakable asperities in region of contact
      IF (rimfrac>0.5) THEN ! Type I : graupel-hail or hail-hail 
         IF (dinsphmin >= dmin .AND. dinsphmin <= dmax ) THEN ! the smallest ice particle is graupel
            a0   = 3.78e-4*(1+0.0079/dinsphmin**1.5)
            C    = 6.3e6*psi ! C
            g    = 0.3  ! g
            Nmax = 100. ! Nmax
            Am   = a0/3.+MAX(2./3.*a0-1./9.*a0*ABS(ptemp-T0),0.) ! Am
         ELSE ! the smallest particle is hail 
            a0   = 4.35e5
            C    = 3.31e5 ! C
            g    = 0.54   ! g
            Nmax = 1000.  ! Nmax
            Am   = a0/3.+MAX(2./3.*a0-1./9.*a0*ABS(ptemp-T0),0.) ! Am
         END IF
      ELSE ! rimfrac<0.5
           ! "If the cloud model in which the scheme is implemented cannot resolve
           ! habits of ice, then all snow/crystals may be treated as dendrites
           ! between -12oC and -17oC and as spatial planar particles at other
           ! subzero temperatures" Section 5 in Phillips et al. 2017
          IF (ptemp>Tmin .AND. ptemp<Tmax) THEN !treating ice particles as dendrites
             C    = 3.09E6*psi ! C
             g    = 0.5-0.25*rimfrac ! g
             Nmax = 100. ! Nmax
             Am   = 1.41E6*(1+100*rimfrac**2.)*(1+3.98e-5/dinsphmin**1.5) ! Am
          ELSE ! treating ice as spatial planar particles
             C    = 7.08E6*psi ! C
             g    = 0.5-0.25*rimfrac ! g
             Nmax = 100. ! Nmax
             Am   = 1.58E7*(1+100*rimfrac**2.)*(1+1.33e-4/dinsphmin**1.5) ! Am
          END IF
      END IF
      
      ! Equivalent-spherical surface area of the colliding particle with the smaller maximum dimension in 1/m2 
      alpha = pi * disphmin**2.
      ! Get K0
      K0 = kinetic_collision_energy(ppres,ptemp,icelarge,icesmall)
      
      ! Ice multiplication factor or number of secondary ice particles produced per ice-ice collision    
      imf_phillips =  MIN(alpha*Am*(1-exp(-(C*K0/alpha/Am)**g)),Nmax)
        
     END FUNCTION imf_phillips

   ! ------------------------------------------------------------------------------------------------------------
   
   
   ! ------------------------------------------------------------------------------------------------------------
     REAL FUNCTION kinetic_collision_energy(ppres,ptemp,pice1,pice2) result(K0) 
     
        USE mo_submctl, ONLY : rd, pstand, pi, spec
        USE classSection, ONLY : Section

        REAL, INTENT(in)  :: ptemp,ppres    ! ptemp in K and ppres in Pa

        TYPE(Section), INTENT(in) :: pice1
        TYPE(Section), INTENT(in) :: pice2
      
        REAL :: m1, m2  ! Masses of  ice particles
        REAL :: v1,v2   ! Terminal velocities of ice particles
      
        ! This is repeating a LOT of the stuff already done once in coagulation kernels,
        ! which is BS and sad... But can't do much about it currently.
        REAL :: visc             ! Viscosity of air [kg/(m s)]
        REAL :: rhoa             ! air density      [kg/(m3)]      
        REAL :: mfp              ! air mean free path [m]
            
        K0 = 0.
        rhoa = ppres/(rd*ptemp)
        visc = (7.44523e-3*SQRT(ptemp**3))/(5093.*(ptemp+110.4)) 
        mfp = (1.656e-10*ptemp+1.828e-8)*pstand/ppres
	
        ! Hydrometeor 1 in the colliding-pair
        ! Get the ice particle terminal velocity 
        v1 = getvelocity(pice1,visc,rhoa,mfp)
        ! Get the ice particle mass
        m1 = mip(pice1) 
      
        ! Hydrometeor 2 in the colliding-pair
        ! Get the ice particle terminal velocity 
        v2 = getvelocity(pice2, visc,rhoa, mfp)
        ! Get the ice particle mass
        m2 = mip(pice2) 
      
        ! Get the collision kinetic energy 
        ! K0 = 0.5 * (m1*m2/(m1 + m2)) * (v1 - v2)**2
        ! Sotipoulou et al 2020 includes a correction factor in these expressions
        ! to account for underestimates when the terminal velocities are
        ! too close or equal
        K0 = 0.5 * (m1*m2/(m1+m2)) * (ABS((1.7*v1 - v2)**2.0 - 0.3*v1*v2))**0.5
      
     END FUNCTION kinetic_collision_energy
     
     ! ----------------------------------------------------------------------------------------------------------
     
     REAL FUNCTION getvelocity(pice, visc, rhoa, mfp) result(vt)
                
        USE mo_submctl, ONLY : spec
        USE mo_particle_external_properties, ONLY : terminal_vel
        USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
        USE classSection, ONLY : Section
        
        TYPE(Section), INTENT(in) :: pice
        
        REAL :: visc             ! Viscosity of air [kg/(m s)]
        REAL :: rhoa             ! air density      [kg/(m3)]      
        REAL :: mfp              ! air mean free path [m]
        
        TYPE(t_shape_coeffs) :: ishape ! Ice shape coefficients
        REAL :: knud, beta       ! Particle knudsen number and Cunningham correction
        REAL :: rhoip            ! rimed fraction weighted average density of ice particle   
        REAL :: mrim,mpri,ncice  ! rimed and unrimed bin ice mass mix rats, ice number concentration

        INTEGER :: iwa, iri      ! compound index for water and rimed ice in the spec derived data 
         
        iwa = spec%getIndex('H2O')
        iri = spec%getIndex('rime')

        knud  = 2.*mfp/pice%dnsp
        beta  = 1.+knud*(1.142+0.558*exp(-0.999/knud)) 
        mrim  = pice%volc(iri) * spec%rhori
        mpri  = SUM(pice%volc(1:iwa)) * spec%rhoic ! Cutting a little corners here with the volc...
        ncice = pice%numc
        rhoip = (mrim*spec%rhori + mpri*spec%rhoic ) / ( mrim + mpri )
        CALL getShapeCoefficients(ishape,mpri,mrim,ncice)
        
        vt    = terminal_vel(pice%dwet,rhoip,rhoa,visc,beta,4,ishape,pice%dnsp)
     
     END FUNCTION getvelocity
    ! ----------------------------------------------------------------------------------------------------------
    
    REAL FUNCTION mip(pice)
         
       USE mo_submctl, ONLY : spec
       USE classSection, ONLY : Section
        
       TYPE(Section), INTENT(in) :: pice 
       REAL :: mrim,mpri,ncice  ! rimed and unrimed bin ice mass mix rats, ice number concentration in [1/kg]
       INTEGER :: iwa, iri      ! compound index for water and rimed ice in the spec derived data 
         
       iwa = spec%getIndex('H2O')
       iri = spec%getIndex('rime')
         
       mrim = pice%volc(iri) * spec%rhori
       mpri = SUM(pice%volc(1:iwa)) * spec%rhoic ! Cutting a little corners here with the volc...
       ncice = pice%numc
         
       ! Single particle mass
       mip = (mrim+mpri)/ ncice ! Count mean mass for ice particles in the size bin [kg]
             
    END FUNCTION mip
    ! -----------------------------------------------------------------------------------------------------------
    
END MODULE mo_salsa_SIP_IIBR
