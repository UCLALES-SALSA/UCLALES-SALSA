MODULE mo_salsa_secondary_ice
  USE mo_salsa_types, ONLY : ice,cloud,precp, rateDiag
  USE mo_submctl, ONLY : nice, pi6, spec, icebins
  IMPLICIT NONE


  ! Arrays to track the number and mass of frozen drops due to ice collection
  ! diagnosed from the coagulation routines for each timestep. Initialized in
  ! mo_salsa_init.
  REAL, ALLOCATABLE :: nfrzn_hm(:,:,:), mfrzn_hm(:,:,:)
  REAL, ALLOCATABLE :: nfrzn_df(:,:,:), mfrzn_df(:,:,:)

  ! Diameter limit for ice and liquid bins. For Hallet-Mossop, require ice bin diameter > dlimit
  ! and frozen drop diameter < dlimit. For drop fracturing, ice bin diameter < dlimit and frozen
  ! drop diameter > dlimit. This is similar to Qu et al. (2022) and avoids overlap between the
  ! two processes
  REAL, PARAMETER :: dlimit = 100.e-6   
  
  CONTAINS

    SUBROUTINE rimesplintering(kbdim,kproma,klev,ptemp,ptstep)
      ! The Hallet-Mossop secondary ice production by rime splintering
      INTEGER, INTENT(in) :: kbdim,kproma,klev
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      
      REAL, PARAMETER :: c1 = 3.6e8   ! Splinters generated per milligram of new rime
      REAL, PARAMETER :: tmax = 270.16, tmid = 268.16, tmin = 265.16
      REAL, PARAMETER :: Dsplint = 10.e-6  ! Assumed splinter diameter
      
      INTEGER :: ii,jj,bb
      REAL :: dN  ! Number of splintered rime
      REAL :: dm  ! Mass of splintered rime
      INTEGER :: iwa, iri
      LOGICAL :: lt13umgt25um(kbdim,klev)
      INTEGER :: splbin  ! Target bin for splinters
      
      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")

      ! Convert freezing rates to changes over timestep
      mfrzn_hm = mfrzn_hm * ptstep
      nfrzn_hm = nfrzn_hm * ptstep
      
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
      splbin = COUNT(icebins < Dsplint)      

      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma

               IF ( nfrzn_hm(ii,jj,bb) > 0. .AND. lt13umgt25um(ii,jj) ) THEN
                  dN = 0.
                  ! Number of generated splinters. The maximum rate is at 268.16 K
                  IF ( ptemp(ii,jj) < tmax .AND. ptemp(ii,jj) >= tmid ) THEN                     
                     dN = c1 * mfrzn_hm(ii,jj,bb) * (tmax - ptemp(ii,jj))/(tmax-tmid) ! linear slope across the temp range                     
                  ELSE IF ( ptemp(ii,jj) <= tmid .AND. ptemp(ii,jj) > tmin ) THEN                     
                     dN = c1 * mfrzn_hm(ii,jj,bb) * (ptemp(ii,jj) - tmin)/(tmid-tmin)                     
                  END IF

                  ! This will assume that the splinters consist of frozen spheres 10 um in diameter.
                  dm = spec%rhori*dN*pi6*Dsplint**3

                  IF (ice(ii,jj,bb)%volc(iri) < dm/spec%rhori) WRITE(*,*) "SECICE HM FAIL " 
                  
                  ice(ii,jj,splbin)%numc = ice(ii,jj,splbin)%numc + dN
                  ice(ii,jj,splbin)%volc(iwa) = ice(ii,jj,splbin)%volc(iwa) + dm/spec%rhoic
                  ice(ii,jj,bb)%volc(iri) = ice(ii,jj,bb)%volc(iri) - dm/spec%rhori

                  ! Secondary ice diagnostics
                  ice(ii,jj,1)%SIP_rmspl = ice(ii,jj,1)%SIP_rmspl + dN                  
                  CALL rateDiag%rmsplrate%Accumulate(n=dN/ptstep)

               END IF
                  
            END DO
         END DO
      END DO

      ! Reset the tracking arrays
      mfrzn_hm = 0.
      nfrzn_hm = 0.
      
    END SUBROUTINE rimesplintering

    ! -------

    SUBROUTINE dropfracturing(kbdim,kproma,klev,ptemp,ptstep)
      !
      ! -------------------------------------------------------
      ! The drop fracturing SIP from Lawson et al. 2015
      ! the inputs include the mass and number of frozen drizzle in ice particle bins.
      ! It is assumed this mass accumulation to ice is present identically as it comes out of the
      ! coagulation routines, i.e. the bin redistribution should NOT be calculated between sec ice
      ! and coagulation. Doing so would result in the loss of the required information on drop freezing.
      ! 
      ! -----------------------------------------------
      INTEGER, INTENT(in) :: kbdim,kproma,klev
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: ptstep
      REAL, PARAMETER :: tmin = 248.15, tmax = 271.15
      REAL, PARAMETER :: c1 = 2.5e-11

      REAL :: Nnorm          ! Normalization factor for distributing fragments
      REAL :: ddmean         ! Mean diameter of frozen drops per ice bin
      REAL :: dN,dm          ! Total number and mass of fragments generated per ice bin
      REAL :: dNb(nice)      ! Number of fragments distributed to ice bins
      REAL :: dVb(nice)      ! Volume of fragments distributed to ice bins
      INTEGER :: bb,bb1,ii,jj,iri,iwa, nimax, npmax
      REAL :: icediams(nice), icebw(nice)

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
      
      ! index of the largest ice bin where fragments are distributed... could use this or icediams?
      nimax = COUNT(icebins <= dlimit)

      DO bb = 1,nimax       ! Assume drop fracturing to take place from ice bins < 100 um (Qu et al 2022)  
         DO jj = 1,klev
            DO ii = 1,kproma
               
               IF (ptemp(ii,jj) < tmin .OR. ptemp(ii,jj) > tmax .OR. nfrzn_df(ii,jj,bb) < 1.e-10) CYCLE ! Outside temperature range, see Keinert et al 2020
               
               ! Mean diameter of the frozen drops on current ice bin
               ddmean = (mfrzn_df(ii,jj,bb)/nfrzn_df(ii,jj,bb)/spec%rhowa/pi6)**(1./3.)
               ! THIS IS NOW CHECKED ALREADY IN COAGULATION ROUTINES
               IF ( ddmean < dlimit ) CYCLE  ! Require the freezing drop diameter to be > 100 um (as in Qu et al. 2022)

               ! Ice bin index corresponding to the mean frozen drop diameter minus one; for fragment distribution
               npmax = MAX(COUNT(icebins <= ddmean) - 1, 1) ! Max just in case, though it should never be anywhere close to 1 at this point
                              
               ! Total number of fragments generated from current ice bin (Lawson et al 2015)
               dN = nfrzn_df(ii,jj,bb) * c1 * (ddmean*1.e6)**4

               ! Assume the (mass of!) fragments is distributed evenly to ice bins smaller than the frozen drop (Lawson et al 2015).
               ! For this, first distribute dN as d**-3 into small ice bins.
               dNb = 0.
               dVb = 0.
               dNb(1:npmax) = 1./icediams(1:npmax)**3              ! density function               
               Nnorm = SUM(dNb(1:npmax)*icebw(1:npmax))              ! Normalization factor
               dNb(1:npmax) = dN * dNb(1:npmax)*icebw(1:npmax)/Nnorm ! Distributed bin concentrations of fragments
               dVb(1:npmax) = dNb(1:npmax) * pi6*icediams(1:npmax)**3        ! Distributed volume of fragments: should be approx. indetical
               dm = SUM(dVb(1:npmax)*spec%rhori)                    ! Total mass of fragments assumimng they consist of rime

               ! Remove the fragment mass
               ! Total mass in fragments shouldn't exceed rime mass in current bin (something wrong if it does)
               IF (ice(ii,jj,bb)%volc(iri) < dm/spec%rhori) WRITE(*,*) "SECICE DF FAIL " 
               ice(ii,jj,bb)%volc(iri) = ice(ii,jj,bb)%volc(iri) - dm/spec%rhori
               
               ! Allocate the fragments to ice bins
               ice(ii,jj,1:npmax)%numc = ice(ii,jj,1:npmax)%numc + dNb(1:npmax)  ! Should we also reduce the number in the ice source bins?
                                                                                 ! Probably this doesn't matter very much and we have NO information
                                                                                 ! about how many drops specifically fragment in the first place.
               ice(ii,jj,1:npmax)%volc(iri) = ice(ii,jj,1:npmax)%volc(iri) + &   ! Again, fragments assumed to consist of rime.
                    dVb(1:npmax)  

               ! for diagnostics
               ice(ii,jj,1:npmax)%SIP_drfr = ice(ii,jj,1:npmax)%SIP_drfr + dNb(1:npmax)
               CALL rateDiag%drfrrate%Accumulate(n=dN/ptstep)
            END DO
         END DO
      END DO

      ! Reset the tracking arrays
      mfrzn_df = 0.
      nfrzn_df = 0.
      
    END SUBROUTINE dropfracturing
    
END MODULE mo_salsa_secondary_ice
