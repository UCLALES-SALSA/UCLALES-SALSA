
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
   !       ISSUES:
   !           Process selection should be made smarter - now just lots of ifs
   !           inside loops. Bad.
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
                           lssecice, ice_halmos,           &
                           lcgupdt, lscoag


      USE mo_salsa_coagulation_kernels

      USE mo_salsa_coagulation_processes

      USE mo_salsa_secondary_ice, ONLY : halletmossop
      
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

      LOGICAL :: any_lt13(kbdim,klev), any_gt25(kbdim,klev)
      
      ! For Hallet-Mossop
      REAL :: drimdt(kbdim,klev,nice)  ! Volume change in rime due to liquid collection in the presense of liquid
                                       ! hydrometeors with diameters both < 13um and >25um
      REAL :: pre1, pre2, post1, post2

      
      !-----------------------------------------------------------------------------
      !-- 1) Coagulation to coarse mode calculated in a simplified way: ------------
      !      CoagSink ~ Dp in continuum regime, thus we calculate
      !      'effective' number concentration of coarse particles
 
      !-- 2) Updating coagulation coefficients -------------------------------------
      
      nspec = spec%getNSpec(type="total")  ! Note this includes the rime index even if it is not used/present.
                                           ! In classSection the volume concentration array is hardcoded to allocate
                                           ! everything, so this doesn't matter, but is ofcourse extra work. So we
                                           ! need extra work for this...
      iri = spec%getIndex("rime")
      
      ! Since this is done here, it won't really be necessary in the subsequent coagulation routines
      ! (its called at least in coagulation_kernels)
      DO bb = 1,ntotal
         DO jj = 1,klev
            DO ii = 1,kproma
               CALL allSALSA(ii,jj,bb)%updateDiameter(limit=.TRUE.,type="all")
               CALL allSALSA(ii,jj,bb)%updateRhomean()
            END DO
         END DO
      END DO

      ! Calculate new kernels every timestep if low freq updating is NOT used,
      ! or when low freq IS used AND it is the update timestep.
      IF (lscoag%mode == 1 .OR. lcgupdt ) &
           CALL update_coagulation_kernels(kbdim,klev,ppres,ptemp)
      
      any_aero = ANY( aero(:,:,:)%numc > aero(:,:,:)%nlim ) .AND. &
                 ANY( [lscgaa,lscgca,lscgpa,lscgia] )
      any_cloud = ANY( cloud(:,:,:)%numc > cloud(:,:,:)%nlim ) .AND. &
                  ANY( [lscgcc,lscgca,lscgpc,lscgic] ) 
      any_precp = ANY( precp(:,:,:)%numc > precp(:,:,:)%nlim ) .AND. &
                  ANY( [lscgpp,lscgpa,lscgpc,lscgip])
      any_ice = ANY( ice(:,:,:)%numc > ice(:,:,:)%nlim ) .AND. &
                ANY( [lscgii,lscgia,lscgic,lscgip] )

      any_lt13 = ANY( cloud(:,:,:)%numc > cloud(:,:,:)%nlim .AND. cloud(:,:,:)%dwet < 13.e-6, DIM=3 )
      any_gt25 = ANY( cloud(:,:,:)%numc > cloud(:,:,:)%nlim .AND. cloud(:,:,:)%dwet > 25.e-6, DIM = 3 )   &
            .OR. ANY( precp(:,:,:)%numc > precp(:,:,:)%nlim .AND. precp(:,:,:)%dwet > 25.e-6, DIM = 3 )
            
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !-- 3) New particle and volume concentrations after coagulation -------------
      !                 GENERALIZE THE PTEMP STATEMENT
      IF (any_ice .AND. ALL(ptemp < 273.15)) THEN

         ! For H-M: Store the "old" rime volumes
         drimdt(:,:,:) = ice(:,:,:)%volc(iri)
         
         CALL coag_ice(kbdim,klev,nspec,ptstep) 

         ! For H-M: Take the change in rime after collection processes
         drimdt(:,:,:) = ice(:,:,:)%volc(iri) - drimdt(:,:,:)

         ! H-M rime splintering
         IF (lssecice%state .AND. ice_halmos) &
              CALL halletmossop(kbdim,kproma,klev,(any_lt13 .AND. any_gt25),ptemp,drimdt)
         
      END IF

      ! Get new nspec that omits the rime. See comments in the beginning..
      nspec = spec%getNSpec(type="wet")

      IF (any_precp) THEN 
         CALL coag_precp(kbdim,klev,nspec,ptstep)

      END IF
     
      IF (any_aero) &
           CALL coag_aero(kbdim,klev,nspec,ptstep)

      IF (any_cloud) &
           CALL coag_cloud(kbdim,klev,nspec,ptstep)

      
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
                           pcsa,    pcocnv, pcocsv,    &
                           pchno3,  pcnh3,  prv,prs, prsi,      &
                           ptemp,   ppres,  ptstep,  ppbl       )

      USE mo_salsa_nucleation
      USE mo_submctl, ONLY :      &
         ntotal, lscndgas,                  & 
         lscndh2oae, lscndh2ocl, lscndh2oic, & ! Condensation to aerosols, clouds and ice particles
         nsnucl                     ! nucleation
      USE mo_salsa_types, ONLY : allSALSA
      
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

      ! For semivolatile/water substepping...
      REAL :: ztstep2, zsubtime
      INTEGER :: ii,jj,bb
      
      zxocnv = 0.
      zxsa = 0.
      zj3n3 = 0.
      zrh(1:kbdim,:) = prv(1:kbdim,:)/prs(1:kbdim,:)

      ! Do this once here, no need to redo in subsequent condensation subroutines?
      DO bb = 1,ntotal
         DO jj = 1,klev
            DO ii = 1,kproma
               CALL allSALSA(ii,jj,bb)%updateDiameter(limit=.TRUE.,type="all")
               CALL allSALSA(ii,jj,bb)%updateRhomean()
            END DO
         END DO
      END DO

      
      !------------------------------------------------------------------------------

      ! Nucleation
      IF (nsnucl > 0) CALL nucleation(kproma, kbdim,  klev,   krow,   &
                                      ptemp,  zrh,    ppres,          &
                                      pcsa,   pcocnv, ptstep, zj3n3,  &
                                      zxsa,   zxocnv, ppbl            )

      ! Condensation of H2SO4 and organic vapors
      IF (lscndgas) CALL condgas(kproma,  kbdim,  klev,    krow,      &
                                 pcsa, pcocnv, pcocsv,     &
                                 zxsa, ptemp,  ppres, ptstep )


      IF (lscndh2ocl .OR. lscndh2oae .OR. lscndh2oic) &
           CALL partitioning(kproma,kbdim,klev,krow,ppres,ptemp,   &
                             pchno3,pcnh3,prv,prs,ptstep          )
      
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

            ! Done already in the main program
            !DO bb = 1,ntotal
            !   CALL allSALSA(ii,jj,bb)%updateDiameter(.TRUE.,type="wet")
            !END DO
            
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
   SUBROUTINE partitioning(kproma,kbdim,klev,krow,ppres,ptemp,   &
                           pghno3,pgnh3,prv,prs,ptstep)
     
     ! Notes: - Many of the particle loops could be combined as one
     !        - Substepping should be placed inside this subroutine
     !        - Added checks whether semivolatiles in use; the code should be grouped
     !          better to make this less verbose and more efficient
     !        - The process switches of condensing different compounds to different particles
     !          are at the moment essentially broken...
     
     USE mo_salsa_types, ONLY : allSALSA, aero, cloud, precp, ice
     USE classSection
     USE mo_submctl, ONLY : spec,           &
          nbins, ncld, nprc,   &
          in1a, fn2b,          &
          rhowa, mwa, mair,    &   ! Check if rhowa and mwa etc from submctl become deprecated
          surfw0, mvno, mvnh,  &
          mvwa, boltz, rg,     &
          massacc,      &
          avog, pi,pi6,            &
          pstand,              &
          nlim, prlim, d_sa,ntotal, &
          alv
            
     USE aerosol_thermodynamics, ONLY : inorganic_pdfite
     USE mo_salsa_properties, ONLY : equilibration
      
      IMPLICIT NONE
      
      ! Equation numbers refer to those in Jacobson: Fundamentals of Atmospheric Modelling, Second Edition (2005)
      
      INTEGER, INTENT(in)  :: kproma,kbdim,klev,krow
      REAL, INTENT(in) :: ptstep                                 !! Timestep
      REAL, INTENT(in) :: ptemp(kbdim,klev), ppres(kbdim,klev)
      REAL, INTENT(in) :: prs(kbdim,klev)
      
      REAL, INTENT(inout) :: prv(kbdim,klev),      &
                             pghno3(kbdim,klev),   &
                             pgnh3(kbdim,klev)

      REAL :: zph(kbdim,klev,nbins+ncld+nprc)
      
      REAL :: zkelwaae(nbins),  zkelwacd(ncld),  zkelwapd(nprc)      ! Kelvin effects for water
      REAL :: zkelno3ae(nbins), zkelno3cd(ncld), zkelno3pd(nprc)     ! Kelvin effects for HNO3
      REAL :: zkelnh3ae(nbins), zkelnh3cd(ncld), zkelnh3pd(nprc)     ! Kelvin effects for NH3
      
      REAL :: zcwacae(nbins),   zcwanae(nbins),   & ! Current, intermediate and new water in aerosols
              zcno3cae(nbins),  zcno3nae(nbins),  & !  -  ''  - HNO3
              zcnh3cae(nbins),  zcnh3nae(nbins),  & !  -  ''  - NH3
           
              zcwaccd(ncld),    zcwancd(ncld),   & ! -  ''  - water in cloud droplets
              zcno3ccd(ncld),   zcno3ncd(ncld),  & ! -  ''  - HNO3 in cloud droplets
              zcnh3ccd(ncld),   zcnh3ncd(ncld),  & ! -  ''  - NH3 
           
              zcwacpd(nprc),    zcwanpd(nprc),   & ! -  ''  - water in precipitation
              zcno3cpd(nprc),   zcno3npd(nprc),  & ! -  ''  - HNO3 in precipitation
              zcnh3cpd(nprc),   zcnh3npd(nprc)     ! -  ''  - NH3
      
      REAL :: zcwatot,zcnh3tot,zcno3tot
      
      REAL :: zcwac, zcwan                                         ! Current and new water gas concentration
      REAL :: zcno3c, zcno3n                                       ! -  ''  - HNO3
      REAL :: zcnh3c, zcnh3n                                       ! -  ''  - NH3
      
      REAL :: zcgnh3eqae(nbins), zcgno3eqae(nbins), zcgwaeqae(nbins), & ! Equilibrium gas concentrations 
              zcgnh3eqcd(ncld), zcgno3eqcd(ncld),   zcgwaeqcd(ncld), &
              zcgnh3eqpd(nprc), zcgno3eqpd(nprc),   zcgwaeqpd(nprc)
      
      REAL :: zmtwaae(nbins),  zmtwacd(ncld),  zmtwapd(nprc)  ! Mass transfer coefficients for H2O
      REAL :: zmtno3ae(nbins), zmtno3cd(ncld), zmtno3pd(nprc) ! Mass transfer coefficients for HNO3
      REAL :: zmtnh3ae(nbins), zmtnh3cd(ncld), zmtnh3pd(nprc) ! Mass transfer coefficients for NH3
            
      REAL :: zbeta  ! transition correction factor
      REAL :: zdfvap ! Diffusion coefficient for vapors
      
      REAL :: zaw(nbins) ! water activity
      
      REAL :: zHpno3ae(nbins),zHpnh3ae(nbins),zHph2oae(nbins),    & ! H' (Eq (17.99)) for aerosol, 1 = NO3-, 2=NH4+, 3=H2O
              zHpno3cd(ncld),zHpnh3cd(ncld),zHph2ocd(ncld),       & ! H' (Eq (17.99)) for clouds
              zHpno3pd(nprc),zHpnh3pd(nprc),zHph2opd(nprc),       & ! H' (Eq (17.99)) for precipitation
              zexpterm_ae(nbins), & ! exponent term in (17.104) for aerosol bins
              zexpterm_cd(ncld), &  ! exponent term in (17.104) for cloud bins
              zexpterm_pd(nprc)     ! exponent term in (17.104) for precipitation bins
            
      REAL :: zrhoair,            & ! air density [kg m-3]
              zthcond,            & ! thermal conductivity of air
              zdfh2o,zdfno3,zdfnh3                ! diffusion coefficient of water
      
      REAL :: adtc2(2,nbins)
      REAL :: telp,ttot ! Elapsed time
      REAL :: zsum1, zsum2, zhlp1, zhlp2, zhlp3 ! temporary variables 
      REAL :: zDeff_ae, zDeff_cd, zDeff_pd, & ! effective diffusion coefficient
              zDp_ae,   zDp_cd,   zDp_pd
      REAL :: dwet_ae(nbins),   dwet_cd(ncld),   dwet_pd(nprc)
      INTEGER :: nstr,ph_switch= 0
      INTEGER :: ii,jj,kk,cc,tt,new,water, nspec

      ! Juha: zbeta* and zkn* don't need to be arrays - removed allocation
      REAL :: zbetaae,zbetacd,zbetapd,fmax,fmin
      REAL :: zknae,zkncd,zknpd
      REAL :: h2oactae,h2oactcd,h2oactpd
      REAL :: zmfp_h2o, zmfp_no3, zmfp_nh3
      REAL :: rhlimx, alpha_h2o, alpha_no3, alpha_nh3
      INTEGER :: iso4,ioc,ibc,idu,iss,iwa,ino,inh

      !!! Local variables for timestepping
      REAL :: adt ! timestep
      REAL :: subtime

      
      ! Most of these are not needed?
      iso4 = spec%getIndex("SO4",notFoundValue = 0)
      ino = spec%getIndex("NO",notFoundValue = 0)
      inh = spec%getIndex("NH",notFoundValue = 0)
      iwa = spec%getIndex("H2O")
      nspec = spec%getNSpec(type="wet")  ! Includes water, but not rime (ice bins)
      
      nstr = 1
      ! initialize
      !adt = ptstep  ! Cleanup different variables for the same timestep...  !! REORGANIZED
    
      DO jj = 1,klev

         DO ii = 1,kproma
            
            ! density of air
            zrhoair = mair*ppres(ii,jj)/(rg*ptemp(ii,jj))
            zkelno3pd = 1.
            zkelno3cd = 1.
            zkelno3ae = 1.
            zkelnh3pd = 1.
            zkelnh3cd = 1.
            zkelnh3ae = 1.
            zkelwaae = 1.
            zkelwacd = 1.
            zkelwapd = 1.
            
            dwet_ae = 0.
            dwet_cd = 0.
            dwet_pd = 0.

            ! Done already in the main program
            !DO cc = 1,ntotal
            !   CALL allSALSA(ii,jj,cc)%updateDiameter(limit=.TRUE.,type="all")
            !END DO
            
            ! Kelvin effects; These loops could be combined to one
            ! Diameters can be obtained directly from the Section instance. Perhaps should also implement
            ! the Kelvin parameters there?
            DO cc = 1, nbins
               IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN
                  dwet_ae(cc) = aero(ii,jj,cc)%dwet
                  !( SUM(aero(ii,jj,cc)%volc(1:nspec))/aero(ii,jj,cc)%numc/pi6 )**(1./3.)
                  
                  zkelwaae(cc) = MIN(2., EXP( 4.*surfw0*spec%mwa /     &
                       (rg*ptemp(ii,jj)*spec%rhowa*dwet_ae(cc)) ) )
                  IF (ino > 0) &
                       zkelno3ae(cc) = MIN(2., EXP( 4.*surfw0*spec%mno /  &
                       (rg*ptemp(ii,jj)*spec%rhono*dwet_ae(cc)) ) )
                  IF (inh > 0) &
                       zkelnh3ae(cc) = MIN(2., EXP( 4.*surfw0*spec%mnh /  &
                       (rg*ptemp(ii,jj)*spec%rhonh*dwet_ae(cc)) ) )
               END IF
            END DO
            
            DO cc = 1, ncld
               IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN 
                  dwet_cd(cc) = cloud(ii,jj,cc)%dwet
                  !( SUM(cloud(ii,jj,cc)%volc(1:nspec))/cloud(ii,jj,cc)%numc/pi6 )**(1./3.)
                  
                  zkelwacd(cc) = MIN(2., EXP( 4.*surfw0*spec%mwa /        &
                       (rg*ptemp(ii,jj)*spec%rhowa*dwet_cd(cc)) ) )
                  IF (ino > 0)  &
                       zkelno3cd(cc) = MIN(2., EXP( 4.*surfw0*spec%mno /  & 
                       (rg*ptemp(ii,jj)*spec%rhono*dwet_cd(cc)) ) )
                  IF (inh > 0)  &
                       zkelnh3cd(cc) = MIN(2., EXP( 4.*surfw0*spec%mnh /  &
                       (rg*ptemp(ii,jj)*spec%rhonh*dwet_cd(cc)) ) )
               END IF
            END DO
            
            DO cc = 1, nprc
               IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN 
                  dwet_pd(cc) = precp(ii,jj,cc)%dwet
                  !( SUM(precp(ii,jj,cc)%volc(1:nspec))/precp(ii,jj,cc)%numc/pi6 )**(1./3.)
                  
                  zkelwapd(cc) = MIN(2., EXP( 4.*surfw0*spec%mwa /        &
                       (rg*ptemp(ii,jj)*spec%rhowa*dwet_pd(cc)) ) )
                  IF (ino > 0)  &
                       zkelno3pd(cc) = MIN(2., EXP( 4.*surfw0*spec%mno /  &
                       (rg*ptemp(ii,jj)*spec%rhono*dwet_pd(cc)) ) )
                  IF (inh > 0)  &
                       zkelnh3pd(cc) = MIN(2., EXP( 4.*surfw0*spec%mnh /  &
                       (rg*ptemp(ii,jj)*spec%rhonh*dwet_pd(cc)) ) )
               END IF
            END DO
            
            ! Current gas concentrations
            IF (ino > 0) &
                 zcno3c = pghno3(ii,jj)/avog
            IF (inh > 0) &
                 zcnh3c = pgnh3(ii,jj)/avog
            zcwac = (prv(ii,jj)/spec%mwa)*zrhoair
            
            ! Current particle concentrations
            IF (ino>0) THEN
               zcno3cae(1:nbins) = aero(ii,jj,1:nbins)%volc(ino)*spec%rhono/spec%mno
               zcno3ccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(ino)*spec%rhono/spec%mno
               zcno3cpd(1:nprc) = precp(ii,jj,1:nprc)%volc(ino)*spec%rhono/spec%mno
               zcno3tot = zcno3c + SUM(zcno3cae(1:nbins)) + SUM(zcno3ccd(1:ncld)) + SUM(zcno3cpd(1:nprc))
            END IF
            IF (inh>0) THEN
               zcnh3cae(1:nbins) = aero(ii,jj,1:nbins)%volc(inh)*spec%rhonh/spec%mnh
               zcnh3ccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(inh)*spec%rhonh/spec%mnh
               zcnh3cpd(1:nprc) = precp(ii,jj,1:nprc)%volc(inh)*spec%rhonh/spec%mnh
               zcnh3tot = zcnh3c + SUM(zcnh3cae(1:nbins)) + SUM(zcnh3ccd(1:ncld)) + SUM(zcnh3cpd(1:nprc))
            END IF
            zcwacae(1:nbins) = aero(ii,jj,1:nbins)%volc(iwa)*spec%rhowa/spec%mwa
            zcwaccd(1:ncld) = cloud(ii,jj,1:ncld)%volc(iwa)*spec%rhowa/spec%mwa
            zcwacpd(1:nprc) = precp(ii,jj,1:nprc)%volc(iwa)*spec%rhowa/spec%mwa
            zcwatot = zcwac + SUM(zcwacae(1:nbins)) + SUM(zcwaccd(1:ncld)) + SUM(zcwacpd(1:nprc)) ! total amount of water
               
            ! Thermodynamics calls needed for water activity and eq mole concentrations for semivolatiles
            ! aerosol 
            !CALL thermoequil(aero(ii,jj,:),nbins,aero(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),          &
            !     zcgno3eqae,zcgnh3eqae,zcgwaeqae,h2oactae,zpH(ii,jj,1:nbins), &
            !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep              )
            ! cloud droplets            
            !CALL thermoequil(cloud(ii,jj,:),ncld,cloud(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),                     &
            !     zcgno3eqcd,zcgnh3eqcd,zcgwaeqcd,h2oactcd,zpH(ii,jj,nbins+1:nbins+ncld), &
            !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep                         )
            ! precipitation
            !CALL thermoequil(precp(ii,jj,:),nprc,precp(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),                              &
            !     zcgno3eqpd,zcgnh3eqpd,zcgwaeqpd,h2oactpd,zpH(ii,jj,nbins+ncld+1:nbins+ncld+nprc), &
            !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep                                   )
            
            ! 1) Condensation / evaporation of water
            ! Diffusion coefficients
            zdfh2o = DIFFC(ptemp(ii,jj),ppres(ii,jj),1)   ! VAIHDETTU ALK PER !!!!!!!!!!!!!!!!!!!!!!
            !zdfh2o = ( 5./(16.*avog*zrhoair*1.e-3*(3.11e-8)**2) ) * &
            !   SQRT( rg*1.e7*ptemp(ii,jj)*mair*1.e3*(spec%mwa+mair)*1.e3/( 2.*pi*spec%mwa*1.e3 ) )
            !zdfh2o = zdfh2o*1.e-4
            
            IF (ino > 0)  &
                 zdfno3 = DIFFC(ptemp(ii,jj),ppres(ii,jj),2) 
            IF (inh > 0)  &
                 zdfnh3 = DIFFC(ptemp(ii,jj),ppres(ii,jj),4) 
            
            zthcond = 0.023807 + 7.1128e-5*(ptemp(ii,jj) - 273.16) ! Thermal conductivity of air 
            
            !-- 2) Transition regime correction factor for particles ---------
            !  
            !  Fuchs and Sutugin (1971), In: Hidy et al. (ed.)
            !  Topics in current aerosol research, Pergamon.  
            !
            !  Size of condensing molecule considered only for 
            !  nucleation mode (3 - 20 nm) 
            !
            ! mean free path [m]
            zmfp_h2o = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(3.11e-8)**2))*   &
                 SQRT(mair/(mair+spec%mwa))
            IF (ino > 0)  &
                 zmfp_no3 = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(3.07e-8)**2))*   &
                 SQRT(mair/(mair+spec%mno))
            IF (inh > 0)  &
                 zmfp_nh3 = 1.e-2*(mair*1.e3/(pi*avog*zrhoair*1.e-3*(4.32e-8)**2))*   &
                 SQRT(mair/(mair+spec%mnh))
            alpha_h2o = 1. !mass accomodation coeficient Davidovits, P., et al. 2004
            alpha_no3 = 1. !jacobson
            alpha_nh3 = 1. !jacobson
            
            ! Initialization of variables  
            zcwan = 0.
            zcwanae = 0.
            zcwancd = 0.
            zcwanpd = 0.
            rhlimx = 0.95
            fmax = 0.3
            fmin = -0.2
               
            !! Subtimestep loop
            adt = 1.e-4
            subtime = 0.
            DO WHILE (subtime < ptstep)

               ! Mass transfer coefficients 
               zmtwaae = 0.; zmtwacd = 0.; zmtwapd = 0.
               zmtno3ae = 0.; zmtnh3ae = 0.
               zmtno3cd = 0.; zmtnh3cd = 0.
               zmtno3pd = 0.; zmtnh3pd = 0.
            
               ! Get the equilibrium concentrations before condensation of water;;; nh3,no3 dont matter here, place these in the beginning of each loop?
               !zcgno3eqae=0.; zcgnh3eqae=0.; zcgwaeqae=0.
               !zcgno3eqcd=0.; zcgnh3eqcd=0.; zcgwaeqcd=0.
               !zcgno3eqpd=0.; zcgnh3eqpd=0.; zcgwaeqpd=0.
               
               zHph2oae = 1.
               zmtwaae = 0.
               zcgwaeqae = 0.
               ! 1) Condensation / evaporation of water; THIS IS NOW DONE IN ANY CASE; EQ SOLUTION NOT NEEDED ANYMORE?
               ! aerosol bins
               DO cc = 1, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN 
                     h2oactae = acth2o(aero(ii,jj,cc),zcwacae(cc))    !! these should be taken from the thermodynamics model!!
                     zcgwaeqae(cc)  = h2oactae*prs(ii,jj)*zrhoair/spec%mwa
                     zknae = zmfp_h2o/(0.5*dwet_ae(cc))                                            !16.20
                     zbetaae = 1. + zknae*(( 1.33 + (0.71/zknae))/( 1. + (1./zknae)) +   &
                          4.*(1.-alpha_h2o)/(3.*alpha_h2o) )
                     zbetaae = 1./zbetaae  
                     
                     zDp_ae = zdfh2o*zbetaae
                     zDeff_ae = zDp_ae/           &                                            ! (16.55)
                          ( (spec%mwa*zDp_ae*alv*zkelwaae(cc)*zcgwaeqae(cc)/(zthcond*ptemp(ii,jj))) *  &
                          ((alv*spec%mwa/(rg*ptemp(ii,jj))) - 1.) + 1. )
                     zmtwaae(cc) = aero(ii,jj,cc)%numc*2.*pi*dwet_ae(cc)*zDeff_ae                  ! (16.64)
                     
                     IF(zcgwaeqae(cc) > 0. .AND. zcwacae(cc) > 0.)  &
                          zHph2oae(cc) = zcwacae(cc)/zcgwaeqae(cc)                                     ! (17.99)
                     
                  END IF
               END DO

               zHph2ocd = 1.
               zmtwacd = 0.
               zcgwaeqcd = 0.
               ! cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     h2oactcd = acth2o(cloud(ii,jj,cc),zcwaccd(cc))  
                     zcgwaeqcd(cc)  = h2oactcd*prs(ii,jj)*zrhoair/spec%mwa
                     zkncd = zmfp_h2o/(0.5*dwet_cd(cc))  
                     zbetacd = 1. + zkncd*(( 1.33 + (0.71/zkncd))/( 1. + (1./zkncd)) + & 
                          4.*(1.-alpha_h2o)/(3.*alpha_h2o) )                                    
                     zbetacd = 1./zbetacd
                     
                     zDp_cd = zdfh2o*zbetacd
                     zDeff_cd = zDp_cd /          &                                            ! (16.55)
                          ( (spec%mwa*zDp_cd*alv*zkelwacd(cc)*zcgwaeqcd(cc)/(zthcond*ptemp(ii,jj))) *   &
                          ((alv*spec%mwa/(rg*ptemp(ii,jj))) - 1.) + 1. )
                     zmtwacd(cc) = cloud(ii,jj,cc)%numc*2.*pi*dwet_cd(cc)*zDeff_cd                ! (16.64)
                     
                     IF(zcgwaeqcd(cc) > 0. .AND. zcwaccd(cc) > 0.)    &
                          zHph2ocd(cc) = zcwaccd(cc)/zcgwaeqcd(cc)                                     ! (17.99)
                     
                  END IF
               END DO
               
               zHph2opd = 1.
               zmtwapd = 0.
               zcgwaeqpd = 0.
               ! precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     h2oactpd = acth2o(precp(ii,jj,cc),zcwacpd(cc))   
                     zcgwaeqpd(cc) = h2oactpd*prs(ii,jj)*zrhoair/spec%mwa
                     zknpd = zmfp_h2o/(0.5*dwet_pd(cc))                                            !16.20
                     zbetapd = 1. + zknpd*(( 1.33 + (0.71/zknpd))/( 1. + (1./zknpd)) + &
                          4.*(1.-alpha_h2o)/(3.*alpha_h2o) )
                     zbetapd = 1./zbetapd
                     
                     zDp_pd = zdfh2o*zbetapd
                     zDeff_pd = zDp_pd /           &                                                   ! (16.55)
                          ( (spec%mwa*zDp_pd*alv*zkelwapd(cc)*zcgwaeqpd(cc)/(zthcond*ptemp(ii,jj))) *   &
                          ((alv*spec%mwa/(rg*ptemp(ii,jj))) - 1.) + 1. )
                     zmtwapd(cc) = precp(ii,jj,cc)%numc*2.*pi*dwet_pd(cc)*zDeff_pd                    ! (16.64)
                   
                     
                     IF(zcgwaeqpd(cc) > 0. .AND. zcwacpd(cc) > 0.)      &
                          zHph2opd(cc) = zcwacpd(cc)/zcgwaeqpd(cc)                                     ! (17.99)
                     
                  END IF
               END DO
               
               ! update the gas phase concentration [mol/m3] of water
               IF( prv(ii,jj)/prs(ii,jj) < rhlimx ) THEN 
                  zexpterm_ae(:) = EXP(-adt*zkelwaae(:)*zmtwaae(:)/zHph2oae(:))
                  zexpterm_cd(:) = EXP(-adt*zkelwacd(:)*zmtwacd(:)/zHph2ocd(:))
                  zexpterm_pd(:) = EXP(-adt*zkelwapd(:)*zmtwapd(:)/zHph2opd(:))
                  
                  zhlp1 = SUM( zcwacae(:)*(1.-zexpterm_ae(:)) ) +    &
                       SUM( zcwaccd(:)*(1.-zexpterm_cd(:)) ) +    &
                       SUM( zcwacpd(:)*(1.-zexpterm_pd(:)) ) 
                  zhlp2 = SUM( (zHph2oae(:)/zkelwaae(:))*(1.-zexpterm_ae(:)) ) +   &
                       SUM( (zHph2ocd(:)/zkelwacd(:))*(1.-zexpterm_cd(:)) ) +   &
                       SUM( (zHph2opd(:)/zkelwapd(:))*(1.-zexpterm_pd(:)) )                        
               ELSE
                  zhlp1 = SUM( zmtwaae(:)*zkelwaae(:)*zcgwaeqae(:) )  +     &
                       SUM( zmtwacd(:)*zkelwacd(:)*zcgwaeqcd(:) )  +     &
                       SUM( zmtwapd(:)*zkelwapd(:)*zcgwaeqpd(:) )  
                  zhlp1 = zhlp1 * adt
                  
                  zhlp2 = SUM( zmtwaae(:) ) + SUM( zmtwacd(:) ) + SUM( zmtwapd(:) )
                  zhlp2 = zhlp2 * adt               
               END IF
               
               zcwan = MIN( zcwatot, (zcwac + zhlp1)/(1. + zhlp2) )
               
               ! update the particle phase concentration of water in each bin
               
               !aerosol bins
               DO cc = 1, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN
                     IF(prv(ii,jj)/prs(ii,jj) < rhlimx ) THEN                                                 ! APD
                        zcwanae(cc) = zHph2oae(cc)*zcwan/zkelwaae(cc) +    &                    ! (17.102)
                             (zcwacae(cc) - zHph2oae(cc)*zcwan/zkelwaae(cc))*zexpterm_ae(cc)
                     ELSE                                                                                     ! APC
                        zcwanae(cc) = zcwacae(cc) + MIN(MAX(adt*zmtwaae(cc)*(zcwan-zkelwaae(cc)*zcgwaeqae(cc)), &
                             fmin*zcwacae(cc)),fmax*zcwacae(cc))                               !16.69
                     END IF
                  END IF
               END DO
               
               !cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     IF(prv(ii,jj)/prs(ii,jj) < rhlimx ) THEN                                                  ! APD
                        zcwancd(cc) = zHph2ocd(cc)*zcwan/zkelwacd(cc) +                     &  ! (17.102)
                             (zcwaccd(cc) - zHph2ocd(cc)*zcwan/zkelwacd(cc))*zexpterm_cd(cc)
                     ELSE                                                                                      ! APC
                        zcwancd(cc) = zcwaccd(cc) + MIN(MAX(adt*zmtwacd(cc)*(zcwan-zkelwacd(cc)*zcgwaeqcd(cc)), &
                             fmin*zcwaccd(cc)),fmax*zcwaccd(cc))                               !16.69
                     END IF
                  END IF
               END DO
               
               !precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     IF(prv(ii,jj)/prs(ii,jj) < rhlimx ) THEN                                                  ! APD
                        zcwanpd(cc) = zHph2opd(cc)*zcwan/zkelwapd(cc) +                     &  ! (17.102)
                             (zcwacpd(cc) - zHph2opd(cc)*zcwan/zkelwapd(cc))*zexpterm_pd(cc)
                     ELSE                                                                                      ! APC
                        zcwanpd(cc) = zcwacpd(cc) + MIN(MAX(adt*zmtwapd(cc)*(zcwan-zkelwapd(cc)*zcgwaeqpd(cc)), &
                             fmin*zcwacpd(cc)),fmax*zcwacpd(cc))                               !16.69
                     END IF
                  END IF
               END DO
               
               ! Update the "current" concentrations
               zcwac = zcwan
               zcwacae(:) = zcwanae(:)
               zcwaccd(:) = zcwancd(:)
               zcwacpd(:) = zcwanpd(:)

               subtime  = subtime + adt
               adt = MIN( 0.2, adt*1.5, ptstep-subtime)             
               
            END DO ! Subtimesteploop
            

            ! END WATER CONDENSATION
            ! --------------------------------------------------------------------------------------

            IF (.FALSE.) THEN  !!! POISTA !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! 1) Condensation / evaporation of HNO3 
            ! Initialization of variables
            IF ( ALL([ino,inh,iso4] > 0) ) THEN   ! In practice this should probably be IF ino > 0 AND inh > 0 AND iso4?? CHECK
               zsum1 = 0.
               zsum2 = 0.
               zcno3nae = 0.
               zcno3ncd = 0.
               zcno3npd = 0.
               ! aerosol bins
               DO cc = 1, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN
                     zknae = zmfp_no3/(0.5*dwet_ae(cc))                                 !16.20
                     zbetaae = 1. + zknae*(( 1.33 + (0.71/zknae))/( 1. + (1./zknae)) +  &
                                   4.*(1-alpha_no3)/(3.*alpha_no3) )
                     zbetaae = 1./zbetaae
                     zHpno3ae(cc) = 1.e0                                                    ! initialization
                     IF(zcgno3eqae(cc) > 0. .and. zcno3cae(cc) > 0.) &
                          zHpno3ae(cc) = zcno3cae(cc)/zcgno3eqae(cc)                        ! (17.99)

                     zmtno3ae(cc) = 2.*pi*dwet_ae(cc) *  &                                  ! (16.64)
                                    zdfno3*aero(ii,jj,cc)%numc*zbetaae
                     zhlp3 = max(-200.,-adt*zkelno3ae(cc)*zmtno3ae(cc)/zHpno3ae(cc))
                     zexpterm_ae(cc) = exp(zhlp3)                                           ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcno3cae(cc)*(1.-zexpterm_ae(cc))                      ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpno3ae(cc)/zkelno3ae(cc))*(1.-zexpterm_ae(cc))        ! sum term in Eq (17.104) denominator
                  END IF
               END DO
            
               ! cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     zkncd = zmfp_no3/(0.5*dwet_cd(cc))                                 !16.20
                     zbetacd = 1. + zkncd*(( 1.33 + (0.71/zkncd))/( 1. + (1./zkncd)) + &
                                   4.*(1-alpha_no3)/(3.*alpha_no3) )
                     zbetacd = 1./zbetacd
                     zHpno3cd(cc) = 1.e0                                                    ! initialization
                     IF(zcgno3eqcd(cc) > 0. .and. zcno3ccd(cc) > 0.)   &
                          zHpno3cd(cc) = zcno3ccd(cc)/zcgno3eqcd(cc)                        ! (17.99)

                     zmtno3cd(cc) = 2.*pi*dwet_cd(cc) *  &
                                    zdfno3*cloud(ii,jj,cc)%numc*zbetacd
                     zhlp3 = max(-200.,-adt*zkelno3cd(cc)*zmtno3cd(cc)/zHpno3cd(cc))
                     zexpterm_cd(cc) = exp(zhlp3)                                           ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcno3ccd(cc)*(1.-zexpterm_cd(cc))                      ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpno3cd(cc)/zkelno3cd(cc))*(1.-zexpterm_cd(cc))        ! sum term in Eq (17.104) denominator
                  END IF
               END DO
          
               ! precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     zknpd = zmfp_no3/(0.5*dwet_pd(cc))                                 !16.20
                     zbetapd = 1. + zknpd*(( 1.33 + (0.71/zknpd))/( 1. + (1./zknpd)) + &
                                   4.*(1-alpha_no3)/(3.*alpha_no3) )
                     zbetapd = 1./zbetapd
                     zHpno3pd(cc) = 1.e0                                                    ! initialization
                     IF(zcgno3eqpd(cc) > 0. .and. zcno3cpd(cc) > 0.)   &
                          zHpno3pd(cc) = zcno3cpd(cc)/zcgno3eqpd(cc)                        ! (17.99)

                     zmtno3pd(cc) = 2.*pi*dwet_pd(cc) *  &
                                    zdfno3*precp(ii,jj,cc)%numc*zbetapd
                     zhlp3 = max(-200.,-adt*zkelno3pd(cc)*zmtno3pd(cc)/zHpno3pd(cc)) 
                     zexpterm_pd(cc) = exp(zhlp3)                                           ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcno3cpd(cc)*(1.-zexpterm_pd(cc))                      ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpno3pd(cc)/zkelno3pd(cc))*(1.-zexpterm_pd(cc))        ! sum term in Eq (17.104) denominator
                  END IF
               END DO
     
               ! update the gas phase concentration [mol/m3] of HNO3
               zcno3n = (zcno3c + zsum1)/(1.  + zsum2) ! Eq (17.104)
               
               ! update the particle phase concentration of NO3- in each bin
               !aerosol bins
               DO cc = 1, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN
                     zcno3nae(cc) = max(zcnh3cae(cc)-aero(ii,jj,cc)%volc(iso4)*spec%rhosu/spec%msu*1.99, &  !! WAS volc(1), should refer to SO4?
                                        zHpno3ae(cc)*zcno3n/zkelno3ae(cc) +  &              ! (17.102)
                                        (zcno3cae(cc) - zHpno3ae(cc)*zcno3n/zkelno3ae(cc))*zexpterm_ae(cc))
                  END IF
               END DO
            
               !cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     zcno3ncd(cc) = zHpno3cd(cc)*zcno3n/zkelno3cd(cc) +   &                 ! (17.102)
                                    (zcno3ccd(cc) - zHpno3cd(cc)*zcno3n/zkelno3cd(cc))*zexpterm_cd(cc)
                  END IF
               END DO
               
               !precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     zcno3npd(cc) = zHpno3pd(cc)*zcno3n/zkelno3pd(cc) +   &                 ! (17.102)
                                    (zcno3cpd(cc) - zHpno3pd(cc)*zcno3n/zkelno3pd(cc))*zexpterm_pd(cc)
                  END IF
               END DO
               !CALL THERMODYNAMICS AGAIN
               zcgno3eqae=0.
               zcgnh3eqae=0.
               zcgwaeqae=0.
               zcgno3eqcd=0.
               zcgnh3eqcd=0.
               zcgwaeqcd=0.
               zcgno3eqpd=0.
               zcgnh3eqpd=0.
               zcgwaeqpd=0.

               !!!! LAITETTAVA TAKASIN
               !IF (ALL([ino,inh] > 0)) THEN
               ! aerosol droplets
               !CALL thermoequil(aero(ii,jj,:),nbins,aero(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),           &
               !     zcgno3eqae,zcgnh3eqae,zcgwaeqae,h2oactae,zpH(ii,jj,1:nbins),  &
               !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep               )
               ! cloud droplets
               !CALL thermoequil(cloud(ii,jj,:),ncld,cloud(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),                      &
               !     zcgno3eqcd,zcgnh3eqcd,zcgwaeqcd,h2oactcd,zpH(ii,jj,nbins+1:nbins+ncld),  &
               !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep                          )
               ! precipitation
               !CALL thermoequil(precp(ii,jj,:),nprc,precp(ii,jj,1)%nlim,ptemp(ii,jj),ppres(ii,jj),                              &
               !     zcgno3eqpd,zcgnh3eqpd,zcgwaeqpd,h2oactpd,zpH(ii,jj,nbins+ncld+1:nbins+ncld+nprc), &
               !     ph_switch,prv(ii,jj)/prs(ii,jj),ptstep                                   )
               !END IF
            END IF

            IF ( ALL([ino,inh,iso4] > 0)) THEN   ! Is this ok? 
               ! 2) Condensation / evaporation of NH3
               ! Initialization of variables
               zsum1 = 0.
               zsum2 = 0.
               zcnh3nae = 0.
               zcnh3ncd = 0.
               zcnh3npd = 0.
               ! aerosol bins
               DO cc = nstr, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN
                     zknae = zmfp_nh3/(0.5*dwet_ae(cc))                                  !16.20
                     zbetaae = 1. + zknae*(( 1.33 + (0.71/zknae))/( 1. + (1./zknae)) + &
                                   4.*(1-alpha_nh3)/(3.*alpha_nh3) )
                     zbetaae = 1./zbetaae
                     zHpnh3ae(cc) = 1.e0                                                     ! initialize
                     IF(zcgnh3eqae(cc) > 0. .and. zcnh3cae(cc) > 0.) &
                          zHpnh3ae(cc) = zcnh3cae(cc)/zcgnh3eqae(cc)                         ! (17.99)
                     zmtnh3ae(cc) = 2.*pi*dwet_ae(cc) *  &
                                    zdfnh3*aero(ii,jj,cc)%numc*zbetaae
                     zhlp3 = max(-200.,-adt*zkelnh3ae(cc)*zmtnh3ae(cc)/zHpnh3ae(cc))
                     zexpterm_ae(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcnh3cae(cc)*(1.-zexpterm_ae(cc))                       ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpnh3ae(cc)/zkelnh3ae(cc))*(1.-zexpterm_ae(cc))         ! sum term in Eq (17.104) denominator
                  END IF
               END DO
          
               !cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     zkncd = zmfp_nh3/(0.5*dwet_cd(cc))                                  !16.20
                     zbetacd = 1. + zkncd*(( 1.33 + (0.71/zkncd))/( 1. + (1./zkncd)) + &
                                   4.*(1-alpha_nh3)/(3.*alpha_nh3) )
                     zbetacd = 1./zbetacd
                     zHpnh3cd(cc) = 1.e0                                                     ! initialize
                     IF(zcgnh3eqcd(cc) > 0. .and. zcnh3ccd(cc) > 0.)    &
                          zHpnh3cd(cc) = zcnh3ccd(cc)/zcgnh3eqcd(cc)                         ! (17.99)
                     zmtnh3cd(cc) = 2.*pi*dwet_cd(cc) *  &
                                    zdfnh3*cloud(ii,jj,cc)%numc*zbetacd
                     zhlp3 = max(-200.,-adt*zkelnh3cd(cc)*zmtnh3cd(cc)/zHpnh3cd(cc))
                     zexpterm_cd(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcnh3ccd(cc)*(1.-zexpterm_cd(cc))                       ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpnh3cd(cc)/zkelnh3cd(cc))*(1.-zexpterm_cd(cc))         ! sum term in Eq (17.104) denominator
                  END IF
               END DO
          
               ! precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     zknpd = zmfp_nh3/(0.5*dwet_pd(cc))                                  !16.20
                     zbetapd = 1. + zknpd*(( 1.33 + (0.71/zknpd))/( 1. + (1./zknpd)) + &
                                   4.*(1-alpha_nh3)/(3.*alpha_nh3) )
                     zbetapd = 1./zbetapd
                     zHpnh3pd(cc) = 1.e0                                                     ! initialize
                     IF(zcgnh3eqpd(cc) > 0. .and. zcnh3cpd(cc) > 0.)    &
                          zHpnh3pd(cc) = zcnh3cpd(cc)/zcgnh3eqpd(cc)                         ! (17.99)
                     zmtnh3pd(cc) = 2.*pi*dwet_pd(cc) *  &
                                    zdfnh3*precp(ii,jj,cc)%numc*zbetapd
                     zhlp3 = max(-200.,-adt*zkelnh3pd(cc)*zmtnh3pd(cc)/zHpnh3pd(cc))
                     zexpterm_pd(cc)=exp(zhlp3)                                              ! exponent term in Eq (17.104)
                     zsum1 = zsum1 + zcnh3cpd(cc)*(1.-zexpterm_pd(cc))                       ! sum term in Eq (17.104) numerator
                     zsum2 = zsum2 + (zHpnh3pd(cc)/zkelnh3pd(cc))*(1.-zexpterm_pd(cc))         ! sum term in Eq (17.104) denominator
                  END IF
               END DO
          
               ! update the gas phase concentration [mol/m3] of NH3
               zcnh3n = (zcnh3c + zsum1)/(1. + zsum2) ! (17.104)
               ! update the particle phase concentration of NH3 in each bin
          
               !aerosol bins
               DO cc = 1, nbins
                  IF (aero(ii,jj,cc)%numc > aero(ii,jj,cc)%nlim) THEN   
                     zcnh3nae(cc) = max(min(zHpnh3ae(cc)/zkelnh3ae(cc)*zcnh3n +       &      ! (17.102)
                                            (zcnh3cae(cc) - zHpnh3ae(cc)/zkelnh3ae(cc)*zcnh3n)*zexpterm_ae(cc),     &
                                            zcno3nae(cc)+(2.-1.e-5)*aero(ii,jj,cc)%volc(iso4)*spec%rhosu/spec%msu), &  ! WAS volc(1), now volc(iso4)?
                                        0.95*zcnh3cae(cc))               
                  END IF
               END DO
          
               ! cloud bins
               DO cc = 1, ncld
                  IF (cloud(ii,jj,cc)%numc > cloud(ii,jj,cc)%nlim) THEN
                     zcnh3ncd(cc) = max(min(zHpnh3cd(cc)/zkelnh3cd(cc)*zcnh3n +       &      ! (17.102)
                                            (zcnh3ccd(cc) - zHpnh3cd(cc)/zkelnh3cd(cc)*zcnh3n)*zexpterm_cd(cc),      &
                                            zcno3ncd(cc)+(2.-1.e-5)*cloud(ii,jj,cc)%volc(iso4)*spec%rhosu/spec%msu), & ! WAS volc(1), now volc(iso4)?
                                        0.95*zcnh3ccd(cc))
                  END IF
               END DO
          
               ! precipitation bins
               DO cc = 1, nprc
                  IF (precp(ii,jj,cc)%numc > precp(ii,jj,cc)%nlim) THEN
                     zcnh3npd(cc) = max(min(zHpnh3pd(cc)/zkelnh3pd(cc)*zcnh3n +       &      ! (17.102)
                                            (zcnh3cpd(cc) - zHpnh3pd(cc)/zkelnh3pd(cc)*zcnh3n)*zexpterm_pd(cc),      &
                                            zcno3npd(cc)+(2.-1.e-5)*precp(ii,jj,cc)%volc(iso4)*spec%rhosu/spec%msu), & ! WAS volc(1), now volc(iso4)?
                                        0.95*zcnh3cpd(cc))
                  END IF
               END DO
               
            END IF
            END IF !! IF FALSE POISTA

            
            ! Update the final volume concentrations
            !Make sure that mass balance holds in condensation
            IF (inh > 0) THEN
               zcnh3n = zcnh3tot - ( SUM(zcnh3nae) + SUM(zcnh3ncd) + SUM(zcnh3npd) )
               pgnh3(ii,jj) = zcnh3n*avog        ! convert gas concentration to #/m3
               DO cc = 1,nbins
                  aero(ii,jj,cc)%volc(inh) = zcnh3nae(cc)*spec%mnh/spec%rhonh              ! convert to volume concentration
               END DO
               DO cc = 1,ncld
                  cloud(ii,jj,cc)%volc(inh) = zcnh3ncd(cc)*spec%mnh/spec%rhonh             ! convert to volume concentration
               END DO
               DO cc = 1,nprc
                  precp(ii,jj,cc)%volc(inh) = zcnh3npd(cc)*spec%mnh/spec%rhonh             ! convert to volume concentration
               END DO
            END IF
            IF (ino > 0) THEN
               zcno3n = zcno3tot - ( SUM(zcno3nae) + SUM(zcno3ncd) + SUM(zcno3npd) )
               pghno3(ii,jj) = zcno3n*avog       ! convert gas phase concentration to #/m3
               DO cc = 1,nbins
                  aero(ii,jj,cc)%volc(ino) = zcno3nae(cc)*spec%mno/spec%rhono              ! convert to volume concentration
               END DO
               DO cc = 1,ncld
                  cloud(ii,jj,cc)%volc(ino) = zcno3ncd(cc)*spec%mno/spec%rhono             ! convert to volume concentration
               END DO
               DO cc = 1,nprc
                  precp(ii,jj,cc)%volc(ino) = zcno3npd(cc)*spec%mno/spec%rhono             ! convert to volume concentration
               END DO
            END IF
            zcwan = zcwatot - ( SUM(zcwanae) + SUM(zcwancd) + SUM(zcwanpd) )
            prv(ii,jj) = zcwan*spec%mwa/zrhoair ! convert water vapor concentration to kg/kg
            DO cc = 1,nbins
               aero(ii,jj,cc)%volc(iwa) = zcwanae(cc)*spec%mwa/spec%rhowa                  ! convert to volume concentration
               !CALL aero(ii,jj,cc)%updateRhomean()
            END DO
            DO cc = 1,ncld
               cloud(ii,jj,cc)%volc(iwa) = zcwancd(cc)*spec%mwa/spec%rhowa                 ! convert to volume concentration
               !CALL cloud(ii,jj,cc)%updateRhomean()
            END DO
            DO cc = 1,nprc
               precp(ii,jj,cc)%volc(iwa) = zcwanpd(cc)*spec%mwa/spec%rhowa                 ! convert to volume concentration
               !CALL precp(ii,jj,cc)%updateRhomean()
            END DO

         END DO
      END DO
  
    END SUBROUTINE partitioning

    !
    ! ---------------------------------------------------------------------------------------------------------
    !

    SUBROUTINE thermoequil(ppart,nb,nlim,ptemp,ppress, chno3g, cnh3g, ch2og, h2oact,acidity,ph_switch,prh,ptstep)
    USE mo_salsa_types, ONLY : section
    USE mo_submctl, ONLY : & !rhosu,msu,    &
                           !rhoss,mss,    &
                           !rhono,mno,    &
                           !rhonh,mnh,    &
                           !rhowa,mwa,    &
                           rg, avog, spec

    USE aerosol_thermodynamics, ONLY : inorganic_pdfite

    IMPLICIT NONE

    INTEGER, INTENT(in) :: nb
    TYPE(section), INTENT(in) :: ppart(nb)
    REAL,INTENT(in) :: nlim

    REAL, INTENT(in)  :: ptemp,ppress ,ptstep
    INTEGER, INTENT(in) :: ph_switch
    REAL, INTENT(in) :: prh
    ! equilibrium gas phase concentrations over a flat surface
    REAL, INTENT(out) :: chno3g(nb), &
                         cnh3g(nb),  &
                         ch2og(nb) , h2oact(nb) 
    REAL, INTENT(out) :: acidity(nb)
    !<testing
    INTEGER, PARAMETER :: NCmax=3, NAmax=5, NNmax=1
    REAL :: MOLAL(-NAmax:NCmax) 
    !>
    LOGICAL  :: detailed_thermo
    REAL :: zions(7)                ! mol/m3

    REAL :: zwatertotal,        &   ! Total water in particles (mol/m3) ???
            chcl,               &   ! dummy variable for HCl concentration (not in use)
            zgammas(7)              ! Activity coefficients

    INTEGER :: cc

    REAL :: pmols(nb,7), paw(nb) 
    REAL :: c_ions(7),SO4_rat(nb)

    
    REAL :: zKr,               & ! Equilibrium constants (see 
            zKeq,              & ! Jacobson (1999),  Atmos Environ 33, 3635 - 3649), Table 3
            dx,                & ! Change in ion concentration
            zcwl,              & ! Liquid water mol/m3-air
            zhlp,esl,rsl,H_eff_no3          !

    REAL, PARAMETER ::ztemp0   = 298.15      ! Reference temperature (K)    

    ! FROM/TO ISOROPIA
    REAL :: zwi_S(5)          ! Total species concentrations in moles/m**3 air
    REAL :: zcntrl_s(2)       ! nug for different types of problems solved
                              ! different state of aerosols (deliquescent or
                              ! metastable)
    REAL :: zwt_s(5)          ! ?
    REAL :: zaerliq_S(12)     ! Aqueous-phase concentration array in moles/m**3air
    REAL :: zaersld_S(9)      ! Solid-phase concentration array in moles/m**3 air  
    REAL :: zother_s(6)       ! Solution information array
    REAL :: zgas_s(3)         ! Gas-phase concentration array in moles/m**3 air
    REAL :: Kwa               
    CHARACTER(len=15) ::  scasi_s          ! Returns the subcase which the input corresponds to
    REAL :: GAMA(6), gamma_no3, gamma_nh3, atm = 101325., H_eff_nh3

    ! Switch for thermodynamic model. Only PDFite is available with this distribution.
    INTEGER :: thermo_model=2 !1 = ISOROPIA, 2 = PDFiTE, !3 = AIM !0 = No thermodynamic model

    INTEGER :: iso4,iss,iwa,ino,inh

    ! The calls to thermodynamics models include some hardcoded inputs,
    ! fo must keep track which compounds are configured active in SALSA.
    iso4 = spec%getIndex("SO4",notFoundValue = 0)
    iss = spec%getIndex("SS",notFoundValue = 0)
    ino = spec%getIndex("NO",notFoundValue = 0)
    inh = spec%getIndex("NH",notFoundValue = 0)
    iwa = spec%getIndex("H2O",notFoundValue = 0)

    ! choose how thermodynamical equilibrium is calculated
    ! .true. for detailed PD-Fite calculations of activity coefficients
    ! .false. for fast calculation, assuming ideality for all species
    detailed_thermo = .TRUE.
                        
    ! If there is no HNO3, PD-Fite is not required
    IF (ino == 0) &
         detailed_thermo = .FALSE.        
    IF (ino /= 0) THEN 
       IF(SUM(ppart(:)%volc(ino)*spec%rhono/spec%mno) == 0. .AND. thermo_model==2)  &
            detailed_thermo = .FALSE.
    END IF
    
    ! initialize
    pmols = 0.   ! Only needed for PD-Fite interface but not otherwise used?
    h2oact  = 0.
    ch2og = 0.     ! THIS IS NOT TREATED CORRECTLY AND ENDS UP 0. HOWEVER, IT'S CURRENTLY NOT USED ANYWHERE...
    chno3g = 0.
    cnh3g = 0.
    acidity = 0.   
    DO cc = 1, nb
       IF(ppart(cc)%numc > nlim .and. ppart(cc)%volc(iwa)>0.) THEN
          
          !Calculate initial concentrations of individual ions
          ! Accumulate H+ concentration
          zhlp = 0.
          IF (iso4 > 0) zhlp = zhlp + 2.*ppart(cc)%volc(iso4)*spec%rhosu/spec%msu
          IF (ino > 0) zhlp = zhlp + ppart(cc)%volc(ino)*spec%rhono/spec%mno
          IF (inh > 0) zhlp = zhlp - ppart(cc)%volc(inh)*spec%rhonh/spec%mnh
          zhlp = MAX(zhlp, 1.e-30)
          
          c_ions= 0.
          c_ions(1) = zhlp                                         ! H+
          IF (iso4 > 0) c_ions(2) = ppart(cc)%volc(iso4)*spec%rhosu/spec%msu     ! SO42-
          c_ions(3) = 0.                                           ! HSO4
          IF (ino > 0) c_ions(4) = ppart(cc)%volc(ino)*spec%rhono/spec%mno                                  ! NO3
          IF (inh > 0) c_ions(5) = MIN(ppart(cc)%volc(inh)*spec%rhonh/spec%mnh,c_ions(2)*2.+c_ions(4)-zhlp) ! NH4
          if (iss>0) c_ions(6) = ppart(cc)%volc(iss)*spec%rhoss/spec%mss                                    ! Na
          if (iss>0) c_ions(7) = ppart(cc)%volc(iss)*spec%rhoss/spec%mss                                    ! Cl
 	
	  !getting the activity of water
          MOLAL     = 0.
          MOLAL(-2) = c_ions(2)                                      ! SO42-
          MOLAL(-3) = c_ions(4)                                      ! NO3-
          MOLAL(2)  = MIN((2.-1.e-10)*MOLAL(-2)+MOLAL(-3),c_ions(5)) ! NH4+ 
          MOLAL(1)  = MOLAL(-3)+2.*MOLAL(-2)-MOLAL(2)
          MOLAL     = MOLAL/(ppart(cc)%volc(iwa)*spec%rhowa) 

          zcwl = ppart(cc)%volc(iwa)*spec%rhowa/spec%mwa
          
          IF (detailed_thermo) THEN

             !zcwl = ppart(cc)%volc(iwa)*spec%rhowa/spec%mwa
             paw(cc) = 0.
             
             IF (zcwl/(zcwl+c_ions(2)) > 0.00001) THEN
                SELECT CASE(thermo_model)
                CASE(1)
		   ! -----------------------------------
		   ! THERMODYNAMICS WITH ISORROPIA	
		   !
                   ! We are not licensed to distribute ISORROPIA, hence we leaving out all the components 
 		   ! of it from this distribution. If interested in using ISORROPIA with this code, 
		   ! download the ISORROPIA source codes and place them in src/src_salsa/ directory, 
		   ! include the dependencies to src/Makefile and uncomment the lines below. 

                   ! WI moles/m3, aerosol phase, NA, SO4, HN4, NO3, CL
                   !!zwi_s=0
                   !!IF (iss>0) zwi_s(1)= ppart(cc)%volc(iss)*spec%rhoss/spec%mss       !Na
                   !!IF (iso4>0) zwi_s(2)= ppart(cc)%volc(iso4)*spec%rhosu/spec%msu     !SO4
                   !!IF (inh>0) zwi_s(3)= ppart(cc)%volc(inh)*spec%rhonh/spec%mnh       !NH4
                   !!IF (ino>0) zwi_s(4)= ppart(cc)%volc(ino)*spec%rhono/spec%mno       !NO3
                   !!IF (iss>0) zwi_s(5)= ppart(cc)%volc(iss)*spec%rhoss/spec%mss       !Cl

                   !!zcntrl_s(1)=1 ! reverse problem solved with WI only aerosol phase 
                   !!zcntrl_s(2)=1  ! The aerosol can have both solid+liquid phases (deliquescent) 

                   !!CALL ISOROPIA(zwi_s, min(prh,1.), ptemp,  zcntrl_s,   &
                   !!zwt_s, zgas_s, zaerliq_s, zaersld_s, scasi_s, zother_s)                   

                   !!chno3g(cc) = zgas_s(2)
                   !!cnh3g(cc)  = zgas_s(1)

                   !!IF (zcwl/(zcwl+SUM(c_ions)) > 0.995) THEN
                   !!   H_eff_no3 = 3.2e6*exp(8700*(1./ptemp-1./298.))/(c_ions(1)/ppart(cc)%volc(iwa)*1.e-3)
                   !!   chno3g(cc) = (zwi_s(4)/ppart(cc)%volc(iwa)*1.e-3)/H_eff_no3*atm/(rg*ptemp)
                   !!   Kwa=1.0e-14

                   !!   H_eff_nh3 = 62.0*1.7e-5*exp(450*(1./ptemp-1./298.))*(c_ions(1)/ppart(cc)%volc(iwa)*1.e-3)/Kwa
                   !!   cnh3g(cc) = (zwi_s(3)/ppart(cc)%volc(iwa)*1.e-3)/H_eff_nh3*atm/(rg*ptemp)
                   !!END IF

                CASE(2)
                   zions = 0.
                   zions(1) = zhlp                        ! H+
                   IF (inh>0) zions(2) = ppart(cc)%volc(inh)*spec%rhonh/spec%mnh   ! NH4
                   IF (iss>0) zions(3) = ppart(cc)%volc(iss)*spec%rhoss/spec%mss   ! Na
                   IF (iso4>0) zions(4) = ppart(cc)%volc(iso4)*spec%rhosu/spec%msu ! SO4
                   zions(5) = 0. ! HSO4
                   IF (ino>0) zions(6) = ppart(cc)%volc(ino)*spec%rhono/spec%mno ! NO3; WAS zions(ino), now zions(6)
                   if (iss>0) zions(7) = ppart(cc)%volc(iss)*spec%rhoss/spec%mss ! Cl; WAS zions(inh), now zions(7)

                   ! calculate thermodynamical equilibrium in the particle phase
                   CALL inorganic_pdfite(.9,ptemp,zions,zcwl,chno3g(cc),chcl,cnh3g(cc),zgammas,pmols(cc,:))  ! FIXED RH??

                   !!! WHERE DOES ch2og COME FROM????
                   ch2og(cc) =ch2og(cc)/(rg*ptemp)
                   chno3g(cc)=chno3g(cc)/(rg*ptemp)
                   cnh3g(cc) =cnh3g(cc)/(rg*ptemp)

                CASE(3)
		   ! --------------------------------------------------------------
                   ! THERMODYNAMIC WITH AIM
		   !
                   ! We are not licensed to distribute AIM, hence we leaving out all the components of it from this distribution.
                   ! If interested in using AIM with this code, download the AIM source codes and place them in src/src_salsa/ 
		   ! directory, include dependencies to src/Makefile and uncomment the lines below.

                   !!CALL THERMO(ptemp,MOLAL,ch2og(cc),chno3g(cc),cnh3g(cc),paw(cc))
                   !!ch2og(cc) =ch2og(cc)/(rg*ptemp)
                   !!chno3g(cc)=chno3g(cc)/(rg*ptemp)
                   !!cnh3g(cc) =cnh3g(cc)/(rg*ptemp)
                CASE DEFAULT
                   !do nothing
                END SELECT
                
             END IF
             
             h2oact(cc) = zcwl/(zcwl+sum(c_ions))

          ELSE
             
             ! Equilibrium coefficient (mol/kg)
             zKeq = 1.02e-2*EXP(8.84*(ztemp0/ptemp-1.0)+25.15*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp))) ! Table B7 
          
             ! Convert zKeq to zKr (molm-3) , Equation from Table 3 in Jacobson (1999) 
             !                                Atmos Environ 33, 3635 - 3649                                      
             zKr   = zKeq*(zcwl*spec%mwa)
             ! Eq (17.13)
             dx = (-c_ions(2)-c_ions(1)-zKr &                                   ! Eq (17), Jacobson (1999)
                   +SQRT((c_ions(2)+c_ions(1)+zKr)**2. & 
                   -4.0*(c_ions(1)*c_ions(2)-c_ions(3)*zKr)))/2.0
             cnh3g(cc) = 0.
             chno3g(cc) = 0.

             !Change concentrations
             c_ions(1) = c_ions(1) + dx
             c_ions(2) = c_ions(2) + dx
             c_ions(3) = c_ions(3) - dx
             ch2og(cc) = 0. !zcwl/(zcwl+sum(c_ions))*satvaph2o(ptemp)/(rg*ptemp)  ! Currently saturation mix.rat. is already brought here so there's some redundancy
                                                                              ! THIS CAN LEAD TO AN INCONSISTENCY W.R.P TO HOST;
             
             ! vapor pressure of HNO3 at the droplet surface:
             zKeq=2.5e6*EXP(29.17*(ztemp0/ptemp-1.0)+16.83*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp)))/101325. ! Table B.7 
             zKr = (zKeq*(zcwl*spec%mwa)**2*rg*ptemp)                                                ! Table 3 in Jacobson (1999) 
             chno3g(cc) = 0.
             IF(c_ions(4) > 0.) chno3g(cc) = c_ions(4)*c_ions(1)/zKr                                   !    "
             ! vapor pressure of NH3 at the droplet surface:
             zKeq=2.58e17*EXP(64.02*(ztemp0/ptemp-1.0)+11.44*(1.-ztemp0/ptemp+LOG(ztemp0/ptemp)))/101325.**2 ! Table B.7 
             zKr = (zKeq*(zcwl*spec%mwa*rg*ptemp)**2)                              ! Table 3 in Jacobson (1999) 
             cnh3g(cc) = 0.
             IF(chno3g(cc) > 0.) cnh3g(cc) = c_ions(4)*c_ions(1)/(chno3g(cc)*zKr)                      ! 
             h2oact(cc) = zcwl/(zcwl+sum(c_ions))

             
             
             IF (ph_switch == 1) THEN
		!CALL THERMO(ptemp,MOLAL,ch2og(cc),chno3g(cc),cnh3g(cc),paw(cc))
                acidity(cc) = -log10(MOLAL(1))
             END IF
          END IF
       END IF

    END DO

  END SUBROUTINE thermoequil
    
   
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
                  ! Limit the supersaturation to max 1.01 for the mass transfer
                  ! EXPERIMENTAL
                  zcwsurfae(cc) = MAX(prs(ii,jj),prv(ii,jj)/1.01)*rhoair/spec%mwa

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

    ! ------------------------------------------------------

    FUNCTION DIFFC(Temp,Press,I)
      !
      ! This function is used to calculate the diffusion coefficient from the
      ! vapor in the gas as a function of temperature and pressure after expans.
      ! I is species number.
      !
      Real :: Temp, Press
      !     
      Real :: MYY(6), DN(6), DIFFC
      INTEGER :: I,J
      DATA ( DN(J),MYY(J), J=1,6 ) /  &
           22.0748e-6, 1.6658e+0, 12.98e-6, 1.75e+0, 9.380307e-6,    &
           1.75e+0, 20.0e-6, 1.75e+0, 14.76518e-6, 1.75e+0, 0.0, 0.0 /
      !   
      DIFFC = DN(I) * ( (TEMP/273.15)**MYY(I) ) * 101325. / Press
      !     
      RETURN
    END FUNCTION DIFFC
    
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


   
   !
   ! -------------------------
   !

   FUNCTION satvaph2o(ptemp) RESULT(psat)
      !-----------------------------------------------------------------
      ! Saturation vapour pressure of water vapour
      ! This is a local function for the subroutine *cloud_condensation*
      !
      ! J. Tonttila, FMI, 03/2014
      !-----------------------------------------------------------------
    
      IMPLICIT NONE

      REAL, INTENT(in) :: ptemp

      REAL, PARAMETER ::       &
         za0 = 6.107799961,    &
         za1 = 4.436518521e-1, &
         za2 = 1.428945805e-2, &
         za3 = 2.650648471e-4, &
         za4 = 3.031240396e-6, &
         za5 = 2.034080948e-8, &
         za6 = 6.136820929e-11

      REAL :: zt

      REAL :: psat

      zt = ptemp - 273.15

      psat = za0 + za1*zt + za2*zt**2 + za3*zt**3 +   &
             za4*zt**4 + za5*zt**5 + za6*zt**6

      ! To Pascals
      psat = psat*100.

   END FUNCTION satvaph2o
   
END MODULE mo_salsa_dynamics
