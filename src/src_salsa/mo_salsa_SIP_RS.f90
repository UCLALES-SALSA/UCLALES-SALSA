MODULE mo_salsa_SIP_RS


  IMPLICIT NONE

  SAVE

  PRIVATE
  PUBLIC  :: rimesplintering,nfrzn_rs, mfrzn_rs, dlliq_rs, dltemp_rs
  
  ! Arrays to track the number and mass of frozen drops due to ice collection
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init. Bin dimensions will be (nprc,nice). Possible contribution by 
  ! cloud droplets will be put to the first bin or smth?? 
  REAL, ALLOCATABLE :: nfrzn_rs(:,:,:,:), mfrzn_rs(:,:,:,:)

  REAL :: dlliq_rs = 24.0e-6   ! Min droplet diameter for rime splintering. 

  INTEGER :: dltemp_rs = 5   ! Temperature-dependence for SIP-rime splintering efficiency
  
  REAL, PARAMETER :: Dsplint = 10.e-6  ! Assumed splinter diameter

    
  CONTAINS

  ! -------

   SUBROUTINE rimesplintering(kbdim,kproma,klev,nspec,ppres,ptemp,ptstep)
      USE mo_salsa_types, ONLY : ice, rateDiag
      USE mo_submctl, ONLY : nprc, nice, pi6, spec, icebins, lssiprimespln
      USE classSection, ONLY : Section
      !
      ! -------------------------------------------------------
      ! Secondary ice production following Hallet-Mossop (1976) and Mossop(1974) experiments
      ! the inputs include the mass and number of frozen drizzle in ice particle bins.
      ! It is assumed this mass accumulation to ice is present identically as it comes out of the
      ! coagulation routines, i.e. the bin redistribution should NOT be calculated between sec ice
      ! and coagulation. Doing so would result in the loss of the required information on drop freezing.
      ! 
      ! -----------------------------------------------
      INTEGER, INTENT(in) :: kbdim,kproma,klev,nspec   ! nspec should contain active compounds + rime, i.e. "total"
      REAL, INTENT(in) :: ppres(kbdim,klev),ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      REAL, PARAMETER :: tmin = 238.15, tmax = 273.15  ! From homogenous nucleation limit to 0degC
                                                       ! to widen temperature range and follow implementations
                                                       ! by Sullivan et al. (2018), Sotiropoulou et al. (2020)

      REAL :: ddmean      ! Mean diameter of frozen drops per ice bin
      REAL :: dN,dV       ! Total number and  volume of fragments generated per ice bin
   
      INTEGER :: cc,bb,bb1,ii,jj,iri,iwa, nimax
      REAL :: icediams(nice), icebw(nice)
      REAL :: fragvolc(kbdim,klev,nice,nspec), sinkvolc(kbdim,klev,nice,nspec) ! Volume to be added and removed
      REAL :: fragnumc(kbdim,klev,nice), sinknumc(kbdim,klev,nice)  ! Number to be added and removed
      REAL :: fragv_loc(nice,nspec)  !! Local fragment vol contributions per ice bin
      REAL :: fragn_loc(nice)       !! Local fragment num contributions per ice bin
      REAL :: v_i                    ! Volume of single ice particle in a bin 
      REAL :: sinkv(nspec)           ! sink volume for single collision
      REAL :: frconst                ! constraining fraction for limiting the mass sink to fragments
      REAL, PARAMETER :: inf = HUGE(1.)
      INTEGER :: splbin  ! Target bin for splinters
 
      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")

      ! Assuming splinter diameter as Dsplint, find the corresponding ice bin
      splbin = MAX( COUNT(icebins < Dsplint), 1 )     

      ! Convert freezing rates to changes over timestep
      mfrzn_rs = mfrzn_rs * ptstep
      nfrzn_rs = nfrzn_rs * ptstep
      
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

      ! NO MORE FIXED ICE LIMITS -> CHANGE NIMAX TO NICE
      DO bb = 1,nice      
         DO jj = 1,klev
            DO ii = 1,kproma
               fragv_loc = 0.
               fragn_loc = 0.
               DO cc = 1,nprc
                 IF ( ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR.  & 
                     nfrzn_rs(ii,jj,cc,bb) < 1.e-12 .OR. SUM(ice(ii,jj,bb)%volc(:)) < 1.e-15 .OR. &
                     ice(ii,jj,bb)%numc < ice(ii,jj,bb)%nlim ) CYCLE ! no collection/empty bin
               
                 ! Spherical equivalent diameter of the frozen drops on current ice bin
                 ddmean = (mfrzn_rs(ii,jj,cc,bb)/nfrzn_rs(ii,jj,cc,bb)/spec%rhowa/pi6)**(1./3.)
                 
                 !Require the freezing drop diameter to be larger than dlliq_rs
                 IF ( ddmean < dlliq_rs ) CYCLE                 
              
                 ! Calculate the number of splinters generated per freezing droplet for current bin
                 IF (lssiprimespln%mode == 1) THEN
                    dN = rs_halletmossop(ptemp(ii,jj),mfrzn_rs(ii,jj,cc,bb))
                 ELSE IF (lssiprimespln%mode == 2) THEN
                    dN = rs_mossop(ptemp(ii,jj),nfrzn_rs(ii,jj,cc,bb))
                 END IF

                 ! Volume of a single ice particle in current bin for calculating the number concentration sink.
                 v_i  = pi6*Dsplint**3
           
                 ! This will assume that the splinters consist of frozen spheres 10 um in diameter.
                 dV = dN * v_i

                 ! Allocate the fragments to temporary ice bins 
                 fragnumc(ii,jj,splbin) = fragnumc(ii,jj,splbin) + dN
                     
                 fragvolc(ii,jj,splbin,1:nspec) = fragvolc(ii,jj,splbin,1:nspec) +     &
                       ice(ii,jj,bb)%volc(1:nspec)*MIN( dV/SUM(ice(ii,jj,bb)%volc(1:nspec)), 1. )

                 ! Sink of volume from current bin
                 sinkv(1:nspec) = ice(ii,jj,bb)%volc(1:nspec)* MIN(dV/SUM(ice(ii,jj,bb)%volc(1:nspec)),1.0)
                  
                 sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) + sinkv(1:nspec)
                
                 ! Sink of number concentration from current bin - assume that the splinter volume is equivalent to that of a sphere of 10 um
                 sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) + SUM( sinkv(1:nspec) ) / v_i               
                  
                 ! Secondary ice diagnostics
                 ice(ii,jj,splbin)%SIP_rmspl = ice(ii,jj,splbin)%SIP_rmspl + dN
                  
                 CALL rateDiag%rmsplrate%Accumulate(n=dN/ptstep)
                  
               END DO

           
               !! Safeguard: Allow the fragments to take up to 99% of the source ice bin mass
               IF ( SUM(sinkvolc(ii,jj,bb,1:nspec)) > 0.99 * SUM(ice(ii,jj,bb)%volc(1:nspec)) ) THEN
                  frconst = 0.99 * SUM(ice(ii,jj,bb)%volc(1:nspec)) / SUM(sinkvolc(ii,jj,bb,1:nspec))
                  fragv_loc = fragv_loc * frconst
                  fragn_loc = fragn_loc * frconst
                  sinkvolc(ii,jj,bb,1:nspec) = sinkvolc(ii,jj,bb,1:nspec) * frconst
                  sinknumc(ii,jj,bb) = sinknumc(ii,jj,bb) * frconst
               END IF

               sinknumc(ii,jj,bb) = MIN(sinknumc(ii,jj,bb), 0.99*ice(ii,jj,bb)%numc) !! Additional constrain because for some reason
                                                                                    !! this still failed in the last bin...
               
               fragnumc(ii,jj,:) = fragnumc(ii,jj,:) + fragn_loc(:)
               fragvolc(ii,jj,:,:) = fragvolc(ii,jj,:,:) + fragv_loc(:,:)

               ! POISTA           
               IF ( SUM(sinkvolc(ii,jj,bb,:)) > 0.99*SUM(ice(ii,jj,bb)%volc(1:nspec)) )     &
                    WRITE(*,*)  'SEC ICE ERROR: FRAGMENT MASS EXCEEDS BIN MASS 2', & 
                    SUM(sinkvolc(ii,jj,bb,:)), SUM(fragvolc(ii,jj,:,:)), SUM(ice(ii,jj,bb)%volc(1:nspec))

               IF (0.99*ice(ii,jj,bb)%numc < sinknumc(ii,jj,bb)) &
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
                    WRITE(*,*) 'fragnumc nan',bb,dlliq_rs
               IF ( ANY(fragvolc(ii,jj,bb,:) /= fragvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'fragvolc nan ',bb,dlliq_rs,fragvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) < 0. ) ) &
                    WRITE(*,*) 'sinkvolc nega ',bb,dlliq_rs,sinkvolc(ii,jj,bb,:)
               IF ( ANY(sinkvolc(ii,jj,bb,:) /= sinkvolc(ii,jj,bb,:)) ) &
                    WRITE(*,*) 'sinkvolc nan ',  bb,dlliq_rs,sinkvolc(ii,jj,bb,:)
               IF (fragnumc(ii,jj,bb) > 1.e5) WRITE(*,*) 'fragnumc > 1e5 ',bb,dlliq_rs,fragnumc(ii,jj,bb),    &
                    (SUM(mfrzn_rs(ii,jj,:,bb))/SUM(nfrzn_rs(ii,jj,:,bb))/spec%rhowa/pi6)**(1./3.), &
                    SUM(nfrzn_rs(ii,jj,:,bb)), ice(ii,jj,bb)%numc
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
      mfrzn_rs = 0.
      nfrzn_rs = 0.
            
    END SUBROUTINE rimesplintering

    !----------------------------------------------------------------------------------------------
    REAL FUNCTION rs_halletmossop(ptemp,mfrzn_rs)
    ! Secondary ice production parameterization based on the experiments by
    ! Hallet, J., & Mossop, S. C. (1974).
    ! Production of secondary ice particles during the riming process.
    ! Nature, 249(5452), 26–28.
    ! https://doi.org/10.1038/249026a0
    !
    REAL :: dN             ! Total number of fragments (splinters) generated per ice bin
    REAL, INTENT(in) :: ptemp, mfrzn_rs  ! Temperature, mass of freezing droplets collected by ice
    REAL, PARAMETER :: IMF_RS = 3.5E8 ! Number of splinters per kg of rime accreted by ice larger than 24 um
    !REAL :: hT             ! Temperature-dependent efficiency for rime splintering

    IF (dltemp_rs == 1) THEN          ! function shape tower
       rs_halletmossop = IMF_RS*mfrzn_rs*hT_ferrier(ptemp)
    ELSE IF (dltemp_rs == 2) THEN      ! function shape parabole
       rs_halletmossop = IMF_RS*mfrzn_rs*hT_ziegler(ptemp)    
    ELSE IF (dltemp_rs == 3) THEN      ! function shape triangle
       rs_halletmossop = IMF_RS*mfrzn_rs*ht_cotton(ptemp)
    ELSE IF (dltemp_rs == 4) THEN      ! 100 % efficiency, 5% outside temperature range
       rs_halletmossop = IMF_RS*mfrzn_rs*ht_sullivan(ptemp)
    ELSE IF (dltemp_rs == 5) THEN      ! 100 % efficiency no temperature dependence
       rs_halletmossop = IMF_RS*mfrzn_rs 
    END IF
    
    END FUNCTION rs_halletmossop


    
    !---------------------------------------------------------------------------------------------
    REAL FUNCTION rs_mossop(ptemp,nfrzn_rs)
      ! Secondary ice production parameterization based on the experiments by
      ! Mossop, S. C. (1976). Production of secondary ice particles
      ! during the growth of graupel by riming.
      ! Quarterly Journal of the Royal Meteorological Society, 102(431), 45–57.
      ! https://doi.org/https://doi.org/10.1002/qj.49710243104

    REAL :: dN             ! Total number of fragments (splinters) generated per ice bin
    REAL, INTENT(in) :: ptemp, nfrzn_rs ! Temperature and number of freezing droplets collected by ice
    REAL, PARAMETER :: IMF_RS = 4.0E-3 !Number of splinters per 250 droplets with diameter larger than 24 um
    !REAL :: hT             ! Temperature-dependent efficiency for rime splintering

    IF (dltemp_rs == 1) THEN          ! function shape tower
       rs_mossop = IMF_RS*nfrzn_rs*hT_ferrier(ptemp)
    ELSE IF (dltemp_rs == 2) THEN      ! function shape parabole
       rs_mossop = IMF_RS*nfrzn_rs*hT_ziegler(ptemp)    
    ELSE IF (dltemp_rs == 3) THEN      ! function shape triangle
       rs_mossop = IMF_RS*nfrzn_rs*ht_cotton(ptemp)
    ELSE IF (dltemp_rs == 4) THEN      ! 100 % efficiency, 5% outside temperature range
       rs_mossop = IMF_RS*nfrzn_rs*ht_sullivan(ptemp)
    ELSE IF (dltemp_rs == 5) THEN      ! 100 % efficiency no temperature dependence
       rs_mossop = IMF_RS*nfrzn_rs
    END IF
    
    END FUNCTION rs_mossop

    

    !----------------------------------------------------------------------------------------------
    ! Temperature-dependence of the rime splintering mechanism
    ! All functions in the next section comes from the experimental data
    ! presented by
    ! Hallet, J., & Mossop, S. C. (1974).
    ! Production of secondary ice particles during the riming process.
    ! Nature, 249(5452), 26–28.
    ! https://doi.org/10.1038/249026a0

    REAL FUNCTION hT_ferrier(ptemp)
      ! hT includes the temperature-dependence on the rime-splintering mechanism
      ! It comes from Hallet's experiments (1974) and was presented by 
      ! Ferrier, B.D. A Double-Moment Multiple-Phase Four-Class Bulk Ice Scheme. Part I: Description
      ! doi: 10.1175/1520-0469(1994)051<0249:ADMMPF>2.0.CO;2
      ! It was implemented by
      ! Sotiropoulou, G., Sullivan, S., Savre, J., Lloyd, G., Lachlan-Cope, T., Ekman, A. M. L., & Nenes, A. (2020).
      ! The impact of secondary ice production on Arctic stratocumulus.
      ! Atmospheric Chemistry and Physics, 20(3), 1301–1316.
      ! https://doi.org/10.5194/acp-20-1301-2020
      ! Function shape is 
      REAL, INTENT(in) :: ptemp
      REAL, PARAMETER :: tlims(4) = [265.15, 267.15, 269.15, 271.15] ! [-8.0, -6.0, -4.0, -2.0]+273.15 
     
      hT_ferrier = 0.
      IF ( ptemp > tlims(4) .OR. ptemp < tlims(1) ) THEN
         hT_ferrier = 0.05
      ELSE IF ( ptemp >= tlims(1) .AND. ptemp < tlims(2) ) THEN
         hT_ferrier = 0.50
      ELSE IF ( ptemp >= tlims(2) .AND. ptemp <= tlims(3) ) THEN
         hT_ferrier = 1.00
      ELSE IF ( ptemp > tlims(3) .AND. ptemp <= tlims(4) ) THEN
         hT_ferrier = 0.50
      END IF      

    END FUNCTION  hT_ferrier


   REAL FUNCTION hT_ziegler(ptemp)
     ! hT includes the temperature dependence on the rime-splintering mechanism
     ! It comes from Hallet's experiments (1974) and was presented by
     ! Ziegler, C. L., Ray, P. S., & MacGorman, D. R. (1986).
     ! Relations of Kinematics, Microphysics and Electrification in an Isolated Mountain Thunderstorm.
     ! Journal of Atmospheric Sciences, 43(19), 2098–2115.
     ! https://doi.org/https://doi.org/10.1175/1520-0469(1986)043<2098:ROKMAE>2.0.CO;2
     ! Function shape is parabolic upward or open down with a maximum at 268.16 K = 5 C
 
     REAL, INTENT(in) :: ptemp
     REAL, PARAMETER :: tlims(2) = [265.15, 271.15] ![-8.0, -2.0 ] ! Celsius degrees
     REAL :: t_degC

     hT_ziegler = 0.0
     t_degC = ptemp - 273.15
     
      IF ( ptemp > tlims(2) .OR. ptemp < tlims(1) ) THEN
         hT_ziegler = 0.0
      ELSE IF (ptemp >= tlims(1) .AND. ptemp <= tlims(2) ) THEN
         hT_ziegler = -(t_degC+2.0)*(t_degC+8.0) / 9.0
      END IF

    END FUNCTION ht_ziegler


    REAL FUNCTION hT_sullivan(ptemp)
     ! hT includes the temperature dependence on the rime-splintering mechanism
     ! It comes from Hallet's experiments (1974) and was presented by
     ! Sullivan, S. C., Hoose, C., Kiselev, A., Leisner, T., & Nenes, A. (2018).
     ! Initiation of secondary ice production in clouds.
     ! Atmospheric Chemistry and Physics, 18(3), 1593–1610.
     ! https://doi.org/10.5194/acp-18-1593-2018
     ! Function shape is rectangular with maximum values between -3C and -8C

      REAL, INTENT(in) :: ptemp
      REAL, PARAMETER :: tlims(2) = [265.15, 270.15] ![-8.0, -3.0]+273.15

      hT_sullivan = 0.0
      IF ( ptemp > tlims(2) .OR. ptemp < tlims(1) ) THEN
         hT_sullivan = 0.05
      ELSE IF ( ptemp >= tlims(1) .AND. ptemp <= tlims(2) ) THEN
         hT_sullivan = 1.0
      END IF
      
    END FUNCTION hT_sullivan

    REAL FUNCTION hT_cotton(ptemp)
     ! hT includes the temperature dependence on the rime-splintering mechanism
     ! It comes from Hallet's experiments (1974) and was presented by
     ! Cotton, W. R., Tripoli, G. J., Rauber, R. M., & Mulvihill, E. A. (1986).
     ! Numerical Simulation of the Effects of Varying Ice Crystal Nucleation Rates
     ! and Aggregation Processes on Orographic Snowfall.
     ! Journal of Applied Meteorology and Climatology, 25(11), 1658–1680.
     ! https://doi.org/10.1175/1520-0450(1986)025<1658:NSOTEO>2.0.CO;2
     ! Function shape is triangular with a maximum at 268.16 K = 5 C

      REAL, INTENT(in) :: ptemp
      REAL, PARAMETER :: tlims(3) = [265.15, 268.15, 270.15]

      hT_cotton = 0.0

      IF ( ptemp > tlims(3) .OR. ptemp < tlims(1) ) THEN
         hT_cotton = 0.0
      ELSE IF ( ptemp >= tlims(2) .AND. ptemp <= tlims(3) ) THEN
         hT_cotton = (tlims(3) - ptemp) / 2.0
      ELSE IF ( ptemp >= tlims(1) .AND. ptemp < tlims(2) ) THEN
         hT_cotton = (ptemp - tlims(1)) / 3.0
      END IF   

    END FUNCTION hT_cotton

    
   END MODULE mo_salsa_SIP_RS
