
!****************************************************************
!*                                                              *
!*   MODULE MO_SALSA_DYNAMICS                               *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate aerosol dynamics                              *
!*                                                              *
!****************************************************************

MODULE mo_salsa_dynamics
   IMPLICIT NONE


CONTAINS

   ! this calculated for empty bins too!!!
   ! fxm: test well, esp. self-coagulation (but other bits too!)
   ! AL_note: Diagnostic variables of cond and nucl mass
   !********************************************************************
   !
   ! Subroutine COAGULATION(kproma,kbdim,klev, &
   !       pnaero,pvols,pdwet, &
   !       pcore, ptstep)
   !
   !********************************************************************
   !
   ! Purpose:
   ! --------
   ! Calculates particle loss and change in size distribution
   !  due to (Brownian) coagulation
   !
   !
   ! Method:
   ! -------
   ! Semi-implicit, non-iterative method:
   !  Volume concentrations of the smaller colliding particles
   !  added to the bin of the larger colliding particles.
   !  Start from first bin and use the updated number and volume
   !  for calculation of following bins. NB! Our bin numbering
   !  does not follow particle size in regime 2.
   !
   !Schematic for bin numbers in different regimes:
   !             1                            2
   !    +-------------------------------------------+
   !  a | 1 | 2 | 3 || 4 | 5 | 6 | 7 |  8 |  9 | 10||
   !  b |           ||11 |12 |13 |14 | 15 | 16 | 17||
   !    +-------------------------------------------+
   !
   ! Exact coagulation coefficients for each pressure level
   !  are calculated in subroutine SET_COAGC (in mo_salsa_init)
   !  which is called once at the beginning of the simulation
   !  from model driver. In subroutine COAGULATION, these exact
   !  coefficients are scaled according to current particle wet size
   !  (linear scaling).
   !
   ! Juha: Now also considers coagulation between hydrometeors,
   !       and hydrometeors and aerosols.
   !
   !       Since the bins are organized in terms of the dry size of
   !       of the condensation nucleus, while coagulation kernell is
   !       calculated with the actual hydrometeor size, some assumptions
   !       are laid out:
   !                 1. Cloud droplets from each size bin are lost by
   !                    coagulation with other cloud droplets that have
   !                    larger condensation nucleus.
   !
   !                 2. Cloud droplets from each size bin are lost by
   !                    coagulation with all drizzle bins, regardless of
   !                    the nucleus size in the latter (collection of cloud
   !                    droplets by rain).
   !
   !                 3. Coagulation between drizzle bins acts like 1.
   !
   !
   ! Interface:
   ! ----------
   ! Called from main aerosol model
   !
   !
   ! Coded by:
   ! ---------
   ! Hannele Korhonen (FMI) 2005
   ! Harri Kokkola (FMI) 2006
   ! Tommi Bergman (FMI) 2012
   ! Matti Niskanen(FMI) 2012
   ! Anton Laakso  (FMI) 2013
   ! Juha Tonttila (FMI) 2014
   !
   !---------------------------------------------------------------------

   SUBROUTINE coagulation(kproma,kbdim,klev,   &
                          ptstep,ptemp,ppres   )

     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, allSALSA,      &
                                zccaa, zcccc, zccpp, zccii,             &
                                zccca, zccpa, zccia,                    &
                                zccpc, zccic,                           &
                                zccip

     USE mo_submctl, ONLY: ntotal,nbins,ncld,nprc,nice, &
                           spec,   &
                           lscgaa, lscgcc, lscgpp, lscgii, & 
                           lscgca, lscgpa, lscgia, & 
                           lscgpc, lscgic, &
                           lscgip,                         &
                           lcgupdt, lscoag

      USE mo_salsa_coagulation_kernels

      USE mo_salsa_coagulation_processes
      
      IMPLICIT NONE

      !-- Input and output variables -------------
      INTEGER, INTENT(IN) ::        &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev                         ! number of vertical klev

      REAL, INTENT(IN) ::         &
         ptstep,                    & ! time step [s]
         ptemp(kbdim,klev),         &
         ppres(kbdim,klev)
      !-- Local variables ------------------------

      INTEGER :: nspec, iri
      INTEGER :: ii,jj,bb

      LOGICAL :: any_aero, any_cloud, any_precp, any_ice
      


      !-----------------------------------------------------------------------------
      !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
      !      CoagSink ~ Dp in continuum regime, thus we calculate
      !      'effective' number concentration of coarse particles
 
      !-- 2) Updating coagulation coefficients -------------------------------------
            
      ! Since this is done here, it won't really be necessary in the subsequent coagulation routines
      ! (its called at least in coagulation_kernels)
      DO bb = 1,ntotal
         DO jj = 1,klev
            DO ii = 1,kproma
               CALL allSALSA(ii,jj,bb)%updateDiameter(limit=.TRUE.,type="all")
            END DO
         END DO
      END DO

      ! Calculate new kernels every timestep if low freq updating is NOT used,
      ! or when low freq IS used AND it is the update timestep.
      IF (lscoag%mode == 1 .OR. lcgupdt ) &
           CALL update_coagulation_kernels(kbdim,klev,ppres,ptemp)
      
      any_aero = ANY( [lscgaa,lscgca,lscgpa,lscgia] )
      any_cloud = ANY( [lscgcc,lscgca,lscgpc,lscgic] ) 
      any_precp = ANY( [lscgpp,lscgpa,lscgpc,lscgip] )
      any_ice = ANY( [lscgii,lscgia,lscgic,lscgip] )
            
      !-- 3) New particle and volume concentrations after coagulation -------------
      !                 GENERALIZE THE PTEMP STATEMENT
      IF (any_ice .AND. ALL(ptemp < 273.15)) THEN

         nspec = spec%getNSpec(type="total")  ! Includes rime         
         CALL coag_ice(kbdim,klev,nspec,ptstep) 
         
      END IF

      ! Get new nspec that omits the rime. 
      nspec = spec%getNSpec(type="wet")  
     
      IF (any_precp) &
           CALL coag_precp(kbdim,klev,nspec,ptstep)
         
      IF (any_aero) &
           CALL coag_aero(kbdim,klev,nspec,ptstep)

      IF (any_cloud) &
           CALL coag_cloud(kbdim,klev,nspec,ptstep)


      ! THIS SHOULD NOT OCCUR, PROBABLY A SETUP PROBLEM? SEE IF THIS COULD BE REMOVED.
      ! Sometimes small negative concentrations due to numerical inaccuracy occur;
      ! Put them to zero here and print a warning if larger negative values are present
      DO bb = 1,ntotal
         DO jj = 1,klev
            DO ii = 1,kproma
               IF ( allSALSA(ii,jj,bb)%numc < -1.e-6)  &
                    WRITE(*,*) 'WARNING SALSA_DYNAMICS; numc < -1e-6 ',ii,jj,bb,allSALSA(ii,jj,bb)%numc  
               allSALSA(ii,jj,bb)%numc = MAX(0.,allSALSA(ii,jj,bb)%numc)
               allSALSA(ii,jj,bb)%volc(:) = MAX(0.,allSALSA(ii,jj,bb)%volc(:))
            END DO
         END DO
      END DO              
      

      
   END SUBROUTINE coagulation


   ! fxm: calculated for empty bins too
   ! fxm: same diffusion coefficients and mean free paths used for sulphuric acid
   !      and organic vapours (average values? 'real' values for each?)
   !********************************************************************
   !
   ! Subroutine CONDENSATION(kproma, kbdim,  klev,        &
   !                         pnaero, pvols,  pdwet, plwc, &
   !                         pcsa,   pcocnv, pcocsv,      &
   !                         ptemp,  ppres,  ptstep)
   !
   !********************************************************************
   !
   ! Purpose:
   ! --------
   ! Calculates the increase in particle volume and
   !  decrease in gas phase concentrations due to condensation
   !  of sulphuric acid and two organic compounds (non-volatile
   !  and semivolatile)
   !
   !
   ! Method:
   ! -------
   ! Regime 3 particles only act as a sink for condensing vapours
   !  while their size and composition does not change.
   ! Exception: Soluble fraction of regime 3c particles can change
   !  and thus they can be moved to regime 3b
   !
   ! New gas and aerosol phase concentrations calculated according
   !  to Jacobson (1997): Numerical techniques to solve
   !  condensational and dissolutional growth equations
   !  when growth is coupled to reversible reactions,
   !  Aerosol Sci. Tech., 27, pp 491-498.
   !
   ! fxm: one should really couple with vapour production and loss terms as well
   !      should nucleation be coupled here as well????
   !
   ! Juha: Now does the condensation of water vapour on hydrometeors as well,
   !       + the condensation of semivolatile aerosol species on hydromets.
   !       Modified for the new aerosol datatype. LWC is obtained from %volc(8)
   !
   !
   ! Interface:
   ! ----------
   ! Called from main aerosol model
   !
   !
   ! Coded by:
   ! ---------
   ! Hannele Korhonen (FMI) 2005
   ! Harri Kokkola (FMI) 2006
   ! Juha Tonttila (FMI) 2014
   !
   !---------------------------------------------------------------
   !
   ! Following parameterization has been used:
   ! ------------------------------------------
   !
   ! Molecular diffusion coefficient of condensing vapour [m2/s]
   !  (Reid et al. (1987): Properties of gases and liquids,
   !   McGraw-Hill, New York.)
   !
   ! D = {1.d-7*sqrt(1/M_air + 1/M_gas)*T^1.75} / &
   !  {p_atm/p_stand * (d_air^(1/3) + d_gas^(1/3))^2 }
   !
   ! M_air = 28.965 : molar mass of air [g/mol]
   ! d_air = 19.70  : diffusion volume of air
   ! M_h2so4 = 98.08  : molar mass of h2so4 [g/mol]
   ! d_h2so4 = 51.96  : diffusion volume of h2so4
   !
   !---------------------------------------------------------------

   SUBROUTINE condensation(kproma,  kbdim,  klev,    krow,      &
                           pcsa,   pcocnv,  pcocsv,    &
                           pchno3,  pcnh3,  prv,prs, prsi,      &
                           ptemp,   ppres,  ptstep,  ppbl       )

      USE mo_salsa_nucleation
      USE mo_submctl, ONLY :      &
         lscndgas,                  & 
         lscndh2oae, lscndh2ocl, lscndh2oic, & ! Condensation to aerosols, clouds and ice particles
         nsnucl                     ! nucleation

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

      REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep,                    & ! timestep [s]
         prs(kbdim,klev),           & ! Water vapor saturation mixing ratio
         prsi(kbdim,klev)              ! Saturation mixing ratio    [kg/m3]

      INTEGER :: ppbl(kbdim)           ! Planetary boundary layer top level

      REAL, INTENT(INOUT) ::     &
         prv(kbdim,klev),          & ! Water vapor mixing ratio [kg/kg]
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         pchno3(kbdim,klev),       & ! nitric acid concentration [#/m3]
         pcnh3(kbdim,klev)           ! ammonia concentration [#/m3]

      REAL :: zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
              zxsa(kbdim,klev),           & ! ratio of sulphuric acid and organic vapor in 3nm particles
              zxocnv(kbdim,klev),         &
              zrh(kbdim,klev)

      zxocnv = 0.
      zxsa = 0.
      zj3n3 = 0.
      zrh(1:kbdim,:) = prv(1:kbdim,:)/prs(1:kbdim,:)

      !------------------------------------------------------------------------------

      ! Nucleation
      IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,  &
                                      ptemp,  zrh,    ppres,  &
                                      pcsa,   pcocnv, ptstep, zj3n3,  &
                                      zxsa,   zxocnv, ppbl            )

      ! Condensation of H2SO4 and organic vapors
      IF (lscndgas) CALL condgas(kproma,  kbdim,  klev,    krow,      &
                                 pcsa, pcocnv, pcocsv,     &
                                 zxsa, ptemp,  ppres, ptstep )

      ! Condensation of water vapour
      IF (lscndh2ocl .OR. lscndh2oae .OR. lscndh2oic) &
         CALL gpparth2o(kproma, kbdim, klev,  krow,  &
                        ptemp,  ppres, prs,   prsi,  &
                        prv,   ptstep        )

   END SUBROUTINE condensation

   !
   ! ----------------------------------------------------------------------------------------------------------
   !

   SUBROUTINE condgas(kproma,  kbdim,  klev,    krow,      &
                      pcsa,    pcocnv, pcocsv,  zxsa,      &
                      ptemp,   ppres,  ptstep              )

     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, allSALSA
      USE mo_submctl, ONLY :      &
         pi,                        &
         in1a, in2a,                & ! size bin indices
         fn2b,                      &
         ncld,                      &
         nprc,                      &
         nice,                      &
         ntotal,                    &
         nlim,                      &
         prlim,                     &
         boltz,                     & ! Boltzmann constant [J/K]
         rg,                        & ! molar gas constant [J/(mol K)]
         pstand,                    & ! standard pressure [Pa]
         mvsu, mvoc,                & ! molecular volumes of sulphate and OC [m3]
         spec,                      &
         d_sa,                      & ! diameter of H2SO4 molecule [m]
         massacc,                   & ! mass accomodation coefficients in each bin
         n3                           ! number of molecules in one 3 nm particle [1]

      IMPLICIT NONE

      !-- Input and output variables ----------
      INTEGER, INTENT(IN) ::      &
         kproma,                    & ! number of horiz. grid kproma
         kbdim,                     & ! dimension for arrays
         klev,                      & ! number of vertical klev
         krow

      REAL, INTENT(IN) ::         &
         ptemp(kbdim,klev),         & ! ambient temperature [K]
         ppres(kbdim,klev),         & ! ambient pressure [Pa]
         ptstep                       ! timestep [s]

      REAL, INTENT(INOUT) ::     &
         pcsa(kbdim,klev),         & ! sulphuric acid concentration [#/m3]
         pcocnv(kbdim,klev),       & ! non-volatile organic concentration [#/m3]
         pcocsv(kbdim,klev),       & ! semivolatile organic concentration [#/m3]
         zxsa(kbdim,klev)            ! ratio of sulphuric acid and organic vapor in 3nm particles

      !-- Local variables ----------------------
      INTEGER :: ii, jj    ! loop indices

      REAL ::                        &
         zvisc,                      & ! viscosity of air [kg/(m s)]
         zdfvap,                     & ! air diffusion coefficient [m2/s]
         zmfp,                       & ! mean free path of condensing vapour [m]
         zcs_tot,                    & ! total condensation sink [1/s] (gases)
         zcs_ocsv,                   & ! condensation sink for semivolatile organics [1/s]
         zcs_su,                     & ! condensation sink for sulfate [1/s]
         zcs_ocnv,                   & ! condensation sink for nonvolatile organics [1/s]
                                       ! vapour concentration after time step [#/m3]
         zcvap_new1,                 & ! sulphuric acid
         zcvap_new2,                 & ! nonvolatile organics
         zcvap_new3,                 & ! semivolatile organics
                                       ! change in vapour concentration [#/m3]
         zdvap1,                     & ! sulphuric acid
         zdvap2,                     & ! nonvolatile organics
         zdvap3,                     & ! semivolatile organics

         zdfpart(in1a+1),            & ! particle diffusion coefficient

         zknud(fn2b),                & ! particle Knudsen number
         zknca(ncld),                & ! Knudsen number for cloud droplets and aerosol vapours
         zknpa(nprc),                & ! Knudsen number for rain drops and aerosol vapours
         zknia(nice),                & ! Knudsen number for ice particles and aerosol vapours

         zbeta(fn2b),                & ! transitional correction factor for aerosols
         zbetaca(ncld),              & ! - '' - for condensing aerosol vapours on clouds (is this needed?)
         zbetapa(nprc),              & ! - '' - for condensing aerosol vapours on rain drops
         zbetaia(nice),              & ! - '' - for condensing aerosol vapours on ice (is this needed?)

         zcolrate(fn2b),             & ! collision rate of molecules to particles [1/s]
         zcolrate_ocnv(fn2b),        & ! collision rate of organic molecules to particles [1/s]
         zcolrateca(ncld),           & ! Collision rate of aerosol vapour molecules to cloud drops
         zcolratepa(nprc),           & ! Collision rate of gases to rain drops
         zcolrateia(nice),           & ! Collision rate of aerosol vapour molecules to ice particles

         zdvolsa(fn2b),              & ! change of sulphate volume in each bin [fxm]
         zdvoloc(fn2b),              & !    - " - organics

         zj3n3(kbdim,klev,2),        & ! Formation massrate of molecules in nucleation, [molec/m3s].  (kbdim,klev,1) for H2SO4 and (kbdim,klev,2) for Organic vapor
         zn_vs_c,                    & ! ratio of nucleation of all mass transfer in the smallest bin
         zxocnv(kbdim,klev)

      INTEGER :: ioc, iso4, bb


      ioc = spec%getIndex("OC",notFoundValue = 0)
      iso4 = spec%getIndex("SO4",notFoundValue = 0)

      zj3n3 = 0.
      zxocnv = 0.

      zdvolsa = 0.
      zn_vs_c = 0.
      DO jj = 1, klev
         DO ii = 1, kbdim

            DO bb = 1,ntotal
               CALL allSALSA(ii,jj,bb)%updateDiameter(.TRUE.,type="wet")
            END DO
            
            zdvoloc = 0.

            !-- 1) Properties of air and condensing gases --------------------
            zvisc  = (7.44523e-3*SQRT(ptemp(ii,jj)**3))/(5093.*(ptemp(ii,jj)+110.4))      ! viscosity of air [kg/(m s)]
            zdfvap = 5.1111e-10*ptemp(ii,jj)**1.75*pstand/ppres(ii,jj)                ! diffusion coefficient [m2/s]
            zmfp   = 3.*zdfvap*sqrt(pi*spec%msu/(8.*rg*ptemp(ii,jj)))                      ! mean free path [m]

            !-- 2) Transition regime correction factor for particles ---------
            !
            !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
            !  Topics in current aerosol research, Pergamon.
            !
            !  Size of condensing molecule considered only for
            !  nucleation mode (3 - 20 nm)
            !

            !-- particle Knudsen numbers
            zknud(in1a:in1a+1) = 2.*zmfp/(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)              ! Gases on aerosols
            zknud(in1a+2:fn2b) = 2.*zmfp/aero(ii,jj,in1a+2:fn2b)%dwet

            zknca(1:ncld) = 2.*zmfp/cloud(ii,jj,1:ncld)%dwet          ! Knudsen number for gases on cloud drplets

            zknpa(1:nprc) = 2.*zmfp/precp(ii,jj,1:nprc)%dwet          ! Knudsen number for gases on rain drops

            zknia(1:nice) = 2.*zmfp/ice(ii,jj,1:nice)%dwet          ! Knudsen number for gases on ice particles

            !-- transitional correction factor
            zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &     ! Aerosol + gas
               (3.*massacc)*(zknud+zknud**2))

            zbetaca = 1. + zknca*( 1.33 + (0.71/zknca) )/( 1. + (1./zknca) ) ! Hydrometeor + gas
            zbetaca = 1./zbetaca

            zbetapa = 1. + zknpa*( 1.33 + (0.71/zknpa) )/( 1. + (1./zknpa) ) ! Rain drop + gas
            zbetapa = 1./zbetapa

            zbetaia = 1. + zknia*( 1.33 + (0.71/zknia) )/( 1. + (1./zknia) ) ! ice + gas
            zbetaia = 1./zbetaia

            !-- 3) Collision rate of molecules to particles -------------------
            !
            !  Particle diffusion coefficient considered only for
            !  nucleation mode (3 - 20 nm)
            !

            !-- particle diffusion coefficient [m2/s]
            zdfpart = boltz*ptemp(ii,jj)*zbeta(in1a:in1a+1)/ &
                      (3.*pi*zvisc*aero(ii,jj,in1a:in1a+1)%dwet)

            !-- collision rate (gases on aerosols) [1/s]
            zcolrate = 0.
            zcolrate(in1a:in1a+1) = MERGE( 2.*pi*(aero(ii,jj,in1a:in1a+1)%dwet+d_sa)*    &
                                          (zdfvap+zdfpart)*zbeta(in1a:in1a+1)*           &
                                          aero(ii,jj,in1a:in1a+1)%numc,                 &
                                          0.,                                            &
                                          aero(ii,jj,in1a:in1a+1)%numc > nlim        )

            zcolrate(in1a+2:fn2b) = MERGE( 2.*pi*aero(ii,jj,in1a+2:fn2b)%dwet*zdfvap*       &
                                           zbeta(in1a+2:fn2b)*aero(ii,jj,in1a+2:fn2b)%numc, &
                                           0.,                                               &
                                           aero(ii,jj,in1a+2:fn2b)%numc > nlim        )

            !-- gases on hydrometeors
            zcolrateca = 0.
            zcolrateca(1:ncld) = MERGE( 2.*pi*cloud(ii,jj,1:ncld)%dwet*zdfvap*         &
                                        zbetaca(1:ncld)*cloud(ii,jj,1:ncld)%numc,      &
                                        0.,                                             &
                                        cloud(ii,jj,1:ncld)%numc > nlim           )

            ! Gases on rain drops
            zcolratepa = 0.
            zcolratepa(1:nprc) = MERGE( 2.*pi*precp(ii,jj,1:nprc)%dwet*zdfvap*    &
                                        zbetapa(1:nprc)*precp(ii,jj,1:nprc)%numc, &
                                        0.,                                        &
                                        precp(ii,jj,1:nprc)%numc > prlim       )
            !-- gases on ice particles
            zcolrateia = 0.
            zcolrateia(1:nice) = MERGE( 2.*pi*ice(ii,jj,1:nice)%dwet*zdfvap*      &
                                        zbetaia(1:nice)*ice(ii,jj,1:nice)%numc,   &
                                        0.,                                        &
                                        ice(ii,jj,1:nice)%numc > prlim         )


            !-- 4) Condensation sink [1/s] -------------------------------------

            zcs_tot = sum(zcolrate) + sum(zcolrateca) + sum(zcolratepa)+ sum(zcolrateia)   ! total sink

            !-- 5) Changes in gas-phase concentrations and particle volume -----
            !
            !--- 5.1) Organic vapours ------------------------

            !---- 5.1.1) Non-volatile organic compound: condenses onto all bins
            IF(pcocnv(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. ioc > 0) THEN

               zn_vs_c = 0.

               IF(zj3n3(ii,jj,2) > 1.) zn_vs_c = (zj3n3(ii,jj,2))/(zj3n3(ii,jj,2) + &
                                                 pcocnv(ii,jj) * zcolrate(in1a))

               !   collision rate in the smallest bin, including nucleation and condensation
               !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
               !   equation (16.73)
               zcolrate_ocnv = zcolrate
               zcolrate_ocnv(in1a) = zcolrate_ocnv(in1a) + zj3n3(ii,jj,2)/pcocnv(ii,jj)

               zcs_ocnv = zcs_tot + zj3n3(ii,jj,2)/pcocnv(ii,jj)                   ! total sink for organic vapor

               zcvap_new2 = pcocnv(ii,jj)/(1.+ptstep*zcs_ocnv)                  ! new gas phase concentration [#/m3]
               zdvap2 = pcocnv(ii,jj) - zcvap_new2                                 ! change in gas concentration [#/m3]
               pcocnv(ii,jj) = zcvap_new2                                          ! updating vapour concentration [#/m3]

               zdvoloc = zcolrate_ocnv(in1a:fn2b)/zcs_ocnv*mvoc*zdvap2             ! volume change of particles
                                                                                   !  [m3(OC)/m3(air)]

               aero(ii,jj,in1a:fn2b)%volc(ioc) = aero(ii,jj,in1a:fn2b)%volc(ioc) + & !-- change of volume
                                                 zdvoloc                             !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc) +  &
                                              zcolrateca(1:ncld)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc) +  &
                                              zcolratepa(1:nprc)/zcs_ocnv*mvoc*zdvap2

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc) +  &
                                            zcolrateia(1:nice)/zcs_ocnv*mvoc*zdvap2

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               ! If zxocnv = 0, then the chosen nucleation mechanism does not take into account
               ! the nonvolatile organic vapors and thus the pnaero does not have to be updated.
               IF (zxocnv(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc + &
                     zn_vs_c * zdvoloc(in1a)/mvoc/(n3*zxocnv(ii,jj))
               END IF

            END IF


            !---- 5.1.2) Semivolatile organic compound: regimes 1, 2 and 3
            zcs_ocsv = sum(zcolrate(in2a:fn2b)) +  &       ! sink for semivolatile organics
                       sum(zcolrateca(1:ncld))  +  &       ! ... including condensation on cloud droplets
                       sum(zcolratepa(1:nprc))  +  &       ! and rain drops
                       sum(zcolrateia(1:nice))  !+  &       ! and ice particles

            IF(pcocsv(ii,jj) > 1.e-10 .AND. zcs_ocsv > 1.e-30 .AND. ioc > 0) THEN


               zcvap_new3 = pcocsv(ii,jj)/(1.+ptstep*zcs_ocsv)   ! new gas phase concentration [#/m3]
               zdvap3 = pcocsv(ii,jj) - zcvap_new3                  ! change in gas concentration [#/m3]
               pcocsv(ii,jj) = zcvap_new3                           ! updating gas concentration [#/m3]

               zdvoloc(in2a:fn2b) = zdvoloc(in2a:fn2b) + &                     ! volume change of particles
                                    zcolrate(in2a:fn2b)/zcs_ocsv*mvoc*zdvap3   !  [m3(OC)/m3(air)]

               aero(ii,jj,in1a:fn2b)%volc(ioc) = &                   !-- change of volume due
                  aero(ii,jj,in1a:fn2b)%volc(ioc) + zdvoloc        !   due to condensation in 1a-2b

               ! Condensation on hydromets
               cloud(ii,jj,1:ncld)%volc(ioc) = cloud(ii,jj,1:ncld)%volc(ioc)  +  &
                                              zcolrateca(1:ncld)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on rain drops
               precp(ii,jj,1:nprc)%volc(ioc) = precp(ii,jj,1:nprc)%volc(ioc)  +  &
                                              zcolratepa(1:nprc)/zcs_ocsv*mvoc*zdvap3

               ! Condensation on ice particles
               ice(ii,jj,1:nice)%volc(ioc) = ice(ii,jj,1:nice)%volc(ioc)  +  &
                                            zcolrateia(1:nice)/zcs_ocsv*mvoc*zdvap3

            END IF


            ! ---- 5.2) Sulphate -------------------------------------------
            IF(pcsa(ii,jj) > 1.e-10 .AND. zcs_tot > 1.e-30 .AND. iso4 > 0) THEN

               !-- Ratio of mass transfer between nucleation and condensation

               zn_vs_c = 0.

               IF(zj3n3(ii,jj,1) > 1.) zn_vs_c = (zj3n3(ii,jj,1)) / &
                                                 (zj3n3(ii,jj,1) +  &
                                                 pcsa(ii,jj) * zcolrate(in1a))

               !   collision rate in the smallest bin, including nucleation and condensation
               !   see Mark Z. Jacobson, Fundamentals of Atmospheric Modeling, Second Edition (2005)
               !   equation (16.73)
               zcolrate(in1a) = zcolrate(in1a) + zj3n3(ii,jj,1) / pcsa(ii,jj)

               zcs_su = zcs_tot + zj3n3(ii,jj,1) / pcsa(ii,jj)      ! total sink for sulfate

               !--- Sulphuric acid -------------------------
               !
               zcvap_new1 = pcsa(ii,jj) /(1.+ptstep*zcs_su)         ! new gas phase concentration [#/m3]
               zdvap1 = pcsa(ii,jj) - zcvap_new1                    ! change in gas concentration [#/m3]
               pcsa(ii,jj) = zcvap_new1                             ! updating vapour concentration [#/m3]

               zdvolsa = zcolrate(in1a:fn2b)/zcs_su*mvsu*zdvap1     ! volume change of particles
               ! [m3(SO4)/m3(air)] by condensation

               !-- Change of volume concentration of sulphate in aerosol [fxm]
               aero(ii,jj,in1a:fn2b)%volc(iso4) = aero(ii,jj,in1a:fn2b)%volc(iso4) + zdvolsa

               !-- Clouds
               cloud(ii,jj,1:ncld)%volc(iso4) = cloud(ii,jj,1:ncld)%volc(iso4)  +  &
                                              zcolrateca(1:ncld)/zcs_su*mvsu*zdvap1

               ! Rain drops
               precp(ii,jj,1:nprc)%volc(iso4) = precp(ii,jj,1:nprc)%volc(iso4)  +  &
                                              zcolratepa(1:nprc)/zcs_su*mvsu*zdvap1

               !-- Ice clouds
               ice(ii,jj,1:nice)%volc(iso4) = ice(ii,jj,1:nice)%volc(iso4)  +  &
                                            zcolrateia(1:nice)/zcs_su*mvsu*zdvap1

               !-- Change of number concentration in the smallest bin caused by nucleation
               !   Jacobson (2005), equation (16.75)
               IF (zxsa(ii,jj) > 0.) THEN
                  aero(ii,jj,in1a)%numc = aero(ii,jj,in1a)%numc +          &
                                           zn_vs_c * zdvolsa(in1a)/mvsu/(n3*zxsa(ii,jj))
               END IF

            END IF

         END DO ! kbdim

      END DO ! klev

   END SUBROUTINE condgas

   !
   ! ----------------------------------------------------------------------------------------------------------
   !

   SUBROUTINE gpparth2o(kproma, kbdim,  klev,  krow,     &
                        ptemp,  ppres,  prs,   prsi,     &
                        prv,   ptstep    )
    
     USE mo_salsa_types, ONLY : aero, cloud, precp, ice, rateDiag, allSALSA
     USE mo_submctl, ONLY : nbins, ncld, nprc,    &
          nice, ntotal, &
          spec,                           &
          mair,                         &
          surfw0, surfi0, rg,           &
          pi, pi6, prlim, nlim,      &
          massacc, avog,  &
          in1a, in2a,  &
          fn2b,            &
          lscndh2oae, lscndh2ocl, lscndh2oic, &
          alv, als 
     USE mo_salsa_properties, ONLY : equilibration
     IMPLICIT NONE

      INTEGER, INTENT(in) :: kproma,kbdim,klev,krow
      REAL, INTENT(in) :: ptstep
      REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev), prs(kbdim,klev), prsi(kbdim,klev)

      REAL, INTENT(inout) :: prv(kbdim,klev)

      REAL :: zkelvin(nbins), zkelvincd(ncld), zkelvinpd(nprc), &  ! Kelvin effects
              zkelvinic(nice)
      REAL :: zcwsurfae(nbins), zcwsurfcd(ncld), zcwsurfpd(nprc), & ! Surface mole concentrations
              zcwsurfic(nice)
      REAL :: zmtae(nbins), zmtcd(ncld), zmtpd(nprc),      & ! Mass transfer coefficients
              zmtic(nice)
      REAL :: zwsatae(nbins), zwsatcd(ncld), zwsatpd(nprc), &  ! Water saturation ratios above
              zwsatic(nice)
      REAL :: zcwtot                                        ! Total water mole concentration
      REAL :: zcwc, zcwn, zcwint                            ! Current and new water vapour mole concentrations
      REAL :: zcwcae(nbins), zcwnae(nbins), zcwintae(nbins) ! Current and new water mole concentrations in aerosols
      REAL :: zcwccd(ncld), zcwncd(ncld), zcwintcd(ncld)    !     -  ''  -     in cloud droplets
      REAL :: zcwcpd(nprc), zcwnpd(nprc), zcwintpd(nprc)    !     -  ''  -     in rain drops
      REAL :: zcwcit(nice), zcwnit(nice), zcwintit(nice)     !     -  ''  -     in total (pristine+rimed) ice
      
      REAL :: zorgic(nice)                                  ! Original total ice mole concentration (not sure if this is really necessary,
                                                            ! could just use the "current" value if it wasn't updated in the substepping loop)
      REAL :: zorgri(nice)                                  ! The same for rime
      
      REAL :: zdfh2o, zthcond,rhoair
      REAL :: zbeta,zknud,zmfph2o
      REAL :: zact, zhlp1,zhlp2,zhlp3
      REAL :: adt,ttot
      REAL :: dwet, cap
      REAL :: zrh(kbdim,klev)

      REAL :: zaelwc1(kbdim,klev), zaelwc2(kbdim,klev)

      REAL :: dvice, dvrime, dvitot ! Volume change for pristine and rimed ice

      INTEGER :: nstr
      INTEGER :: ii,jj,cc
      INTEGER :: counter      
      INTEGER :: iwa,irim,nspec

      REAL, ALLOCATABLE :: vrate(:)
      
      zrh(:,:) = prv(:,:)/prs(:,:)
      
      iwa = spec%getIndex("H2O")
      irim = spec%getIndex("rime")
      nspec = spec%getNSpec(type="total")

      ! For diagnostics
      ALLOCATE(vrate(nspec))
      
      ! Calculate the condensation only for 2a/2b aerosol bins
      nstr = in2a

      ! Save the current aerosol water content
      zaelwc1(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      ! For 1a bins do the equilibrium calculation
      CALL equilibration(kproma,kbdim,klev,      &
                         zrh,ptemp,.FALSE. )

      ! If RH < 98 % OR dynamic condensation for aerosols switched off, do equilibrium for all bins
      IF (zrh(1,1) < 0.98 .OR. .NOT. lscndh2oae)  CALL equilibration(kproma,kbdim,klev,      &
                                                                     zrh,ptemp,.TRUE. )

      ! The new aerosol water content after equilibrium calculation
      zaelwc2(:,:) = SUM(aero(:,:,in1a:fn2b)%volc(iwa),DIM=3)*spec%rhowa

      prv(:,:) = prv(:,:) - ( zaelwc2(:,:) - zaelwc1(:,:) )/(ppres(:,:)*mair/(rg*ptemp(:,:)))

      DO jj = 1, klev
         DO ii = 1, kbdim
            ! Necessary?
            IF ( .NOT. ( &
                (ANY(cloud(ii,jj,:)%numc > nlim) .OR. ANY(precp(ii,jj,:)%numc > prlim) .AND. lscndh2ocl) .OR. &
                (ANY(ice(ii,jj,:)%numc > prlim) .AND. lscndh2oic) .OR.  &
                (ANY(aero(ii,jj,:)%numc > nlim) .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) &
                ) ) CYCLE

            rhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))

            ! Diffusion coef
            zdfh2o = ( 5./(16.*avog*rhoair*1.e-3*(3.11e-8)**2) ) * &
               SQRT( rg*1.e7*ptemp(ii,jj)*mair*1.e3*(spec%mwa+mair)*1.e3/( 2.*pi*spec%mwa*1.e3 ) )
            zdfh2o = zdfh2o*1.e-4

            zmfph2o = 3.*zdfh2o*sqrt(pi*spec%mwa/(8.*rg*ptemp(ii,jj))) ! mean free path
            zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air

            ! -- Water vapour (Follows the analytical predictor method by Jacobson 2005)
            zkelvinpd = 1.; zkelvincd = 1.; zkelvin = 1.; zkelvinic = 1.

            zcwc = 0.; zcwint = 0.; zcwn = 0.
            zcwcae = 0.; zcwccd = 0.; zcwcpd = 0.; zcwcit = 0.
            zcwintae = 0.; zcwintcd = 0.; zcwintpd = 0.; zcwintit = 0.
            zcwnae = 0.; zcwncd = 0.; zcwnpd = 0.; zcwnit = 0.
            zwsatae = 0.; zwsatcd = 0.; zwsatpd = 0.; zwsatic = 0.
            
            zmtpd(:) = 0.
            zcwsurfpd(:) = 0.
            zmtcd(:) = 0.
            zcwsurfcd(:) = 0.
            zmtic(:) = 0.
            zcwsurfic(:) = 0.
            zmtae(:) = 0.
            zcwsurfae(:) = 0.

            ! Update particle diameters
            DO cc = 1,ntotal
               CALL allSALSA(ii,jj,cc)%updateDiameter(limit=.TRUE.,type="all")
               CALL allSALSA(ii,jj,cc)%updateRhomean()
            END DO


            ! Cloud droplets --------------------------------------------------------------------------------
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim .AND. lscndh2ocl) THEN
                  ! Wet diameter
                  dwet = cloud(ii,jj,cc)%dwet

                  ! Activity + Kelvin effect
                  zact = acth2o(cloud(ii,jj,cc))
                  zkelvincd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  zcwsurfcd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatcd(cc) = zact*zkelvincd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = cloud(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatcd(cc)*zcwsurfcd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtcd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Rain drops --------------------------------------------------------------------------------
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim .AND. lscndh2ocl) THEN
                  ! Wet diameter
                  dwet = precp(ii,jj,cc)%dwet

                  ! Activity + Kelvin effect
                  zact = acth2o(precp(ii,jj,cc))
                  zkelvinpd(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentrations over flat surface
                  zcwsurfpd(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatpd(cc) = zact*zkelvinpd(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.)*(zknud+zknud**2))

                  ! Mass transfer according to Jacobson
                  zhlp1 = precp(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatpd(cc)*zcwsurfpd(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtpd(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Ice particles --------------------------------------------------------------------------------
            DO cc = 1, nice
               IF (ice(ii,jj,cc)%numc > ice(ii,jj,cc)%nlim .AND. lscndh2oic .AND. ptemp(ii,jj) < 273.15) THEN
                  ! Wet diameter
                  dwet = ice(ii,jj,cc)%dnsp
                  
                  ! Capacitance (analogous to the liquid radius for spherical particles) - edit when needed
                  cap=0.5*dwet
                     
                  ! Activity + Kelvin effect - edit when needed
                  !   Can be calculated just like for sperical homogenous particle or just ignored,
                  !   because these are not known for solid, irregular and non-homogenous particles.
                  !   Ice may not be that far from a sphere, but most particles are large and at least
                  !   growing particles are covered by a layer of pure ice.
                  zact = 1.0
                  ! 
                  zkelvinic(cc) = exp( 4.*surfi0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )
                  
                  ! Saturation mole concentration over flat surface
                  zcwsurfic(cc) = prsi(ii,jj)*rhoair/spec%mwa
                  
                  ! Equilibrium saturation ratio
                  zwsatic(cc) = zact*zkelvinic(cc)
                  
                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                       (3.)*(zknud+zknud**2))
                  
                  ! Mass transfer according to Jacobson
                  zhlp1 = ice(ii,jj,cc)%numc*4.*pi*cap*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*als*zwsatic(cc)*zcwsurfic(cc)/(zthcond*ptemp(ii,jj)) 
                  zhlp3 = ( (als*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.
                  
                  zmtic(cc) = zhlp1/( zhlp2*zhlp3 + 1. )
                  
               END IF
            END DO
            
            ! -- Aerosols: ------------------------------------------------------------------------------------
            DO cc = 1, nbins
               IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim .AND. zrh(ii,jj) > 0.98 .AND. lscndh2oae) THEN
                  ! Wet diameter
                  dwet = aero(ii,jj,cc)%dwet

                  ! Water activity + Kelvin effect
                  zact = acth2o(aero(ii,jj,cc))
                  zkelvin(cc) = exp( 4.*surfw0*spec%mwa / (rg*ptemp(ii,jj)*spec%rhowa*dwet) )

                  ! Saturation mole concentration over flat surface
                  zcwsurfae(cc) = prs(ii,jj)*rhoair/spec%mwa

                  ! Equilibrium saturation ratio
                  zwsatae(cc) = zact*zkelvin(cc)

                  !-- transitional correction factor
                  zknud = 2.*zmfph2o/dwet
                  zbeta = (zknud + 1.)/(0.377*zknud+1.+4./ &
                          (3.*massacc(cc))*(zknud+zknud**2))

                  ! Mass transfer
                  zhlp1 = aero(ii,jj,cc)%numc*2.*pi*dwet*zdfh2o*zbeta
                  zhlp2 = spec%mwa*zdfh2o*alv*zwsatae(cc)*zcwsurfae(cc)/(zthcond*ptemp(ii,jj))
                  zhlp3 = ( (alv*spec%mwa)/(rg*ptemp(ii,jj)) ) - 1.

                  zmtae(cc) = zhlp1/( zhlp2*zhlp3 + 1. )

               END IF
            END DO

            ! Current mole concentrations
            zcwc = prv(ii,jj)*rhoair/spec%mwa
            zcwcae(1:nbins) = aero(ii,jj,1:nbins)%volc(iwa)*spec%rhowa/spec%mwa
            zcwccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(iwa)*spec%rhowa/spec%mwa
            zcwcpd(1:nprc) = precp(ii,jj,1:nprc)%volc(iwa)*spec%rhowa/spec%mwa

            ! Treat the ice types as one mass during the condensation process.
            ! Further assumptions about compositional changes take place after substepping loop
            zcwcit(1:nice) = ice(ii,jj,1:nice)%volc(iwa)*spec%rhoic/spec%mwa
            IF (spec%isUsed("rime")) &
                 zcwcit(1:nice) = zcwcit(1:nice) + ice(ii,jj,1:nice)%volc(irim)*spec%rhori/spec%mwa

            ! Store original values of pristine and rimed ice to preserve info about composition
            zorgic = 0.; zorgri = 0.
            zorgic = ice(ii,jj,1:nice)%volc(iwa)*spec%rhoic/spec%mwa
            IF (spec%isUsed("rime")) &
                 zorgri = ice(ii,jj,1:nice)%volc(irim)*spec%rhori/spec%mwa
            
            zcwtot = zcwc + SUM(zcwcae) + &
                            SUM(zcwccd) + &
                            SUM(zcwcpd) + &
                            SUM(zcwcit)
            ttot = 0.

            zcwintae = zcwcae; zcwintcd = zcwccd; zcwintpd = zcwcpd
            zcwintit = zcwcit
            
            ! Substepping loop
            ! ---------------------------------
            zcwint = 0.
            counter = 0
            DO WHILE (ttot < ptstep)

               adt = 2.e-2
               ! New vapor concentration
               zhlp1 = zcwc + adt * ( SUM(zmtae(nstr:nbins)*zwsatae(nstr:nbins)*zcwsurfae(nstr:nbins))  + &
                                      SUM(zmtcd(1:ncld)*zwsatcd(1:ncld)*zcwsurfcd(1:ncld))              + &
                                      SUM(zmtpd(1:nprc)*zwsatpd(1:nprc)*zcwsurfpd(1:nprc))              + &
                                      SUM(zmtic(1:nice)*zwsatic(1:nice)*zcwsurfic(1:nice)) )         

               zhlp2 = 1. + adt * ( SUM(zmtae(nstr:nbins)) + SUM(zmtcd(1:ncld)) + SUM(zmtpd(1:nprc)) &
                                  + SUM(zmtic(1:nice)) ) 
               zcwint = zhlp1/zhlp2
               zcwint = MIN(zcwint,zcwtot)

               IF ( ANY(aero(ii,jj,:)%numc > aero(ii,jj,:)%nlim) .AND. zrh(ii,jj) > 0.98 ) THEN
                  DO cc = nstr, nbins
                     !zcwintae(cc) = zcwcae(cc) + MIN(MAX(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                     !                                -1.e-2*zcwtot), 1.e-2*zcwtot)
                     zcwintae(cc) = zcwcae(cc) + MIN(MAX(adt*zmtae(cc)*(zcwint - zwsatae(cc)*zcwsurfae(cc)), &
                                                      -1.e-1*zcwcae(cc)), 1.e-1*zcwcae(cc))  
                     zwsatae(cc) = acth2o(aero(ii,jj,cc),zcwintae(cc))*zkelvin(cc)
                  END DO
               END IF
               IF ( ANY(cloud(ii,jj,:)%numc > cloud(ii,jj,:)%nlim) ) THEN
                  DO cc = 1, ncld
                     zcwintcd(cc) = zcwccd(cc) + MIN(MAX(adt*zmtcd(cc)*(zcwint - zwsatcd(cc)*zcwsurfcd(cc)), &
                                                        -1.e-1*zcwccd(cc)), 1.e-1*zcwccd(cc)) !-1.e-2*zcwtot), 1.e-2*zcwtot) 
                     zwsatcd(cc) = acth2o(cloud(ii,jj,cc),zcwintcd(cc))*zkelvincd(cc)
                  END DO
               END IF
               IF ( ANY(precp(ii,jj,:)%numc > precp(ii,jj,:)%nlim) ) THEN
                  DO cc = 1, nprc
                     zcwintpd(cc) = zcwcpd(cc) + MIN(MAX(adt*zmtpd(cc)*(zcwint - zwsatpd(cc)*zcwsurfpd(cc)), &
                                                        -1.e-1*zcwcpd(cc)), 1.e-1*zcwcpd(cc)) !-1.e-2*zcwtot), 1.e-2*zcwtot) 
                     zwsatpd(cc) = acth2o(precp(ii,jj,cc),zcwintpd(cc))*zkelvinpd(cc)
                  END DO
               END IF
               IF ( ANY(ice(ii,jj,:)%numc > ice(ii,jj,:)%nlim) ) THEN
                  DO cc = 1, nice
                     zcwintit(cc) = zcwcit(cc) + MIN(MAX(adt*zmtic(cc)*(zcwint - zwsatic(cc)*zcwsurfic(cc)), &
                                                         -1.e-1*zcwcit(cc)), 1.e-1*zcwcit(cc)) !-1.e-2*zcwtot), 1.e-2*zcwtot) 
                     zwsatic(cc) = zkelvinic(cc)
                  END DO
               END IF

               zcwintae(nstr:nbins) = MAX(zcwintae(nstr:nbins),0.)
               zcwintcd(1:ncld) = MAX(zcwintcd(1:ncld),0.)
               zcwintpd(1:nprc) = MAX(zcwintpd(1:nprc),0.)
               zcwintit(1:nice) = MAX(zcwintit(1:nice),0.)

               ! Update vapor concentration for consistency
               zcwint = zcwtot - SUM(zcwintae(1:nbins)) - &
                                 SUM(zcwintcd(1:ncld))  - &
                                 SUM(zcwintpd(1:nprc))  - &
                                 SUM(zcwintit(1:nice)) 

               ! Update "old" values for next cycle
               zcwcae = zcwintae; zcwccd = zcwintcd; zcwcpd = zcwintpd; zcwcit = zcwintit
               zcwc = zcwint

               ttot = ttot + adt

               counter = counter + 1
               IF (counter > 100) WRITE(*,*) "CONDENSATION SUBSTEPING: SOMETHING'S WRONG"

            END DO ! ADT

            zcwn = zcwint
            zcwnae = zcwintae
            zcwncd = zcwintcd
            zcwnpd = zcwintpd
            zcwnit = zcwintit
            
            prv(ii,jj) = zcwn*spec%mwa/rhoair

            ! Update particle concentrations and diagnostics
            vrate = 0.
            DO cc = 1,nbins
               vrate(iwa) = max(0.,zcwnae(cc)*spec%mwa/spec%rhowa) - aero(ii,jj,cc)%volc(iwa)
               CALL rateDiag%Cond_a%Accumulate(v=vrate(1:nspec))
               aero(ii,jj,cc)%volc(iwa) = max(0.,zcwnae(cc)*spec%mwa/spec%rhowa)
            END DO
            vrate = 0.
            DO cc = 1,ncld
               vrate(iwa) = max(0.,zcwncd(cc)*spec%mwa/spec%rhowa) - cloud(ii,jj,cc)%volc(iwa)
               CALL rateDiag%Cond_c%Accumulate(v=vrate(1:nspec))
               cloud(ii,jj,cc)%volc(iwa) = max(0.,zcwncd(cc)*spec%mwa/spec%rhowa)
            END DO
            vrate = 0.
            DO cc = 1,nprc
               vrate(iwa) = max(0.,zcwnpd(cc)*spec%mwa/spec%rhowa) - precp(ii,jj,cc)%volc(iwa)
               CALL rateDiag%Cond_p%Accumulate(v=vrate(1:nspec))
               precp(ii,jj,cc)%volc(iwa) = max(0.,zcwnpd(cc)*spec%mwa/spec%rhowa)
            END DO
            
            ! Ice particles: Assume deposition to only contribute to pristine ice.
            ! Sublimation will remove identical fractions from both pristine and
            ! rimed ice contents
            vrate = 0.
            DO cc = 1,nice

               ! Take the total change in ice concentration:
               dvitot = zcwnit(cc) - (zorgic(cc) + zorgri(cc))

               IF ( dvitot > 0. ) THEN
                  ! Deposition - assume to increase only the pristine ice
                  dvice = dvitot
                  ice(ii,jj,cc)%volc(iwa) = ice(ii,jj,cc)%volc(iwa) + dvice*spec%mwa/spec%rhoic
                  vrate(iwa) = dvice*spec%mwa/spec%rhoic
                  CALL rateDiag%Cond_i%accumulate(v=vrate(1:nspec))
               ELSE IF ( dvitot < 0. ) THEN
                  ! Sublimation takes an equal fraction out of both pristine and rimed ice
                  dvice = ( zcwnit(cc)/(zorgic(cc)+zorgri(cc)) ) - 1.
                  dvice = dvice * zorgic(cc)
                  dvrime = ( zcwnit(cc)/(zorgic(cc)+zorgri(cc)) ) - 1.
                  dvrime = dvrime * zorgri(cc)
                  vrate(iwa) = dvice*spec%mwa/spec%rhoic
                  vrate(irim) = dvrime*spec%mwa/spec%rhori
                  CALL rateDiag%Cond_i%accumulate(v=vrate(1:nspec))
                  ice(ii,jj,cc)%volc(iwa) = ice(ii,jj,cc)%volc(iwa) + dvice*spec%mwa/spec%rhoic
                  IF(spec%isUsed("rime")) &
                       ice(ii,jj,cc)%volc(irim) = ice(ii,jj,cc)%volc(irim) + dvrime*spec%mwa/spec%rhori                  
               END IF
                  
            END DO

            DO cc = 1,ntotal
               IF (allSALSA(ii,jj,cc)%numc > allSALSA(ii,jj,cc)%nlim) &
                    CALL allSALSA(ii,jj,cc)%updateRhomean()
            END DO
            
         END DO !kproma

      END DO ! klev

      DEALLOCATE(vrate)
      
   END SUBROUTINE gpparth2o
   !-------------------------------------------------------
   REAL FUNCTION acth2o(ppart,pcw)

      USE classSection, ONLY : Section
      USE mo_submctl, ONLY : eps, spec
      IMPLICIT NONE

      TYPE(Section), INTENT(in) :: ppart
      REAL, INTENT(in), OPTIONAL  :: pcw

      REAL :: zns, znw
      INTEGER :: ndry, iwa, nn
      
      ndry = spec%getNSpec(type="dry")
      iwa = spec%getIndex("H2O")
      
      ! This is only relevant for solution particles so use rholiq
      zns = 0.
      DO nn = 1,ndry  ! Leaves out water and non-soluble species (zero dissociation factor)
         zns = zns + spec%diss(nn)*ppart%volc(nn)*spec%rholiq(nn)/spec%MM(nn)
      END DO 

      IF (present(pcw)) THEN
         znw = pcw
      ELSE
         znw = ppart%volc(iwa)*spec%rholiq(iwa)/spec%MM(iwa)
      END IF

      ! Assume activity coefficient of 1 for water...
      acth2o = MAX(0.1,znw/max(eps,(znw+zns)))

   END FUNCTION acth2o

END MODULE mo_salsa_dynamics
