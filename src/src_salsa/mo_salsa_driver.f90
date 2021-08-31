  MODULE mo_salsa_driver
  USE classSection, ONLY : Section, CoagCoe
  USE util, ONLY : getMassIndex !!! IS it good to import this here??? The function is anyway handy here too.
  USE mo_submctl
  USE classFieldArray, ONLY : FieldArray
  USE mo_structured_datatypes
  USE omp_lib
  IMPLICIT NONE

   !---------------------------------------------------------------
   ! 
   ! MO_SALSA_DRIVER:
   ! Contains the primary SALSA input/output variables as well as
   ! Subroutines used to call the main SALSA routine.
   ! 
   ! Juha Tonttila, FMI, 2014
   ! 
   !---------------------------------------------------------------

  SAVE
  
  ! Derived type for taking copies of the "old" values
  TYPE old
     REAL, ALLOCATABLE :: d(:)
  END TYPE old
  INTERFACE old
     PROCEDURE old_cnstr
  END INTERFACE old
  
   ! JT: Variables from SALSA
   ! --------------------------------------------
   ! grid points for SALSA; Note that SALSA is called here as a box model, so these are 1.
   ! If this is ever to change, the driver and possibly parts of the code inside SALSA will need to be revised.
   INTEGER, PARAMETER :: kproma = 1
   INTEGER, PARAMETER :: kbdim = 1
   INTEGER, PARAMETER :: klev = 1
   INTEGER, PARAMETER :: krow = 1
   REAL, PARAMETER    :: init_rh(kbdim,klev) = 0.3
   
   ! --------------------------------------------

   ! Variables for coagulation calls in varying temporal intervals
   ! --------------------------------------------------------------
   ! Storage of the coagulation kernels:
   REAL, ALLOCATABLE :: sto_aa(:,:,:,:,:), sto_cc(:,:,:,:,:), sto_pp(:,:,:,:,:), sto_ii(:,:,:,:,:),  &
                        sto_ca(:,:,:,:,:), sto_pa(:,:,:,:,:), sto_ia(:,:,:,:,:),   &
                        sto_pc(:,:,:,:,:), sto_ic(:,:,:,:,:), sto_ip(:,:,:,:,:)
                        
   LOGICAL :: cgsto_initialized = .FALSE.
   
CONTAINS

   ! 
   !----------------------------------------------------
   ! RUN_SALSA
   ! Performs necessary unit and dimension conversion between
   ! the host model and SALSA module, and calls the main SALSA
   ! routine
   ! 
   ! Partially adobted form the original SALSA boxmodel version.
   ! 
   ! Now takes masses in as kg/kg from LES!! Converted to m3/m3 for SALSA

   ! Juha Tonttila, FMI, 2014
   ! Jaakko Ahola, FMI, 2016
   ! 

  SUBROUTINE run_SALSA(Diag, Prog, nzp, nxp, nyp, ns, wp,     &
                       pa_nactd, pa_vactd, tstep, time,istp,  &
                       level, initialize                      )

      USE mo_salsa_types, ONLY : iaero, faero, icloud, fcloud, iprecp, fprecp, iice, fice, &
                                 rateDiag, allCOAGcoe, aero, cloud, precp, ice, liquid, frozen, allSALSA
      USE mo_salsa, ONLY : salsa
      USE mo_salsa_properties, ONLY  : equilibration
      IMPLICIT NONE
      
      TYPE(FieldArray), INTENT(inout) :: Diag, Prog
      
      INTEGER, INTENT(in) :: nzp,nxp,nyp,ns                       ! dimensions: x,y,z,number of chemical species
      REAL, INTENT(in)    :: tstep                                ! Model timestep length
      REAL, INTENT(in)    :: time                                 ! current time (seconds)
      INTEGER, INTENT(in) :: istp                                 ! Number of current timestep
      LOGICAL, INTENT(in) :: initialize                      
      REAL, INTENT(in)    :: wp(nzp,nxp,nyp)
      REAL, INTENT(out)   :: pa_vactd(nzp,nxp,nyp,ns*ncld) ! mass concentrations of newly activated droplets for calculating the
                                                           ! actual tendency due to new droplet formation.
      REAL, INTENT(out)   :: pa_nactd(nzp,nxp,nyp,ncld)   ! Same for number concentration
      INTEGER, INTENT(in) :: level                         ! thermodynamical level
      
      ! To help coupling and backing up for tendencies
      TYPE(old)          :: npart(4), mpart(4)  ! To store the old values
      TYPE(FloatArray1d) :: ntend(4), mtend(4)  ! Arranged pointers to tendencies
      
      ! Named pointers to get SALSA variables
      TYPE(FloatArray3d), POINTER :: press => NULL(), tk => NULL(), rv => NULL(), &
                                     rt => NULL(), rs => NULL(), rsi => NULL(),   &
                                     pdn => NULL()
      TYPE(FloatArray4d), POINTER :: naerop => NULL(), ncloudp => NULL(),   &
                                     nprecpp => NULL(), nicep => NULL(),    &
                                     naerot => NULL(), ncloudt => NULL(),   &
                                     nprecpt => NULL(), nicet => NULL(),    &

                                     maerop => NULL(), mcloudp => NULL(),   &
                                     mprecpp => NULL(), micep => NULL(),    &
                                     maerot => NULL(), mcloudt => NULL(),   &
                                     mprecpt => NULL(), micet => NULL(),    &

                                     indefp => NULL(),                      &
                                     
                                     gaerop => NULL(), gaerot => NULL()
           
      REAL :: in_p(kbdim,klev), in_t(kbdim,klev), in_rv(kbdim,klev), in_rs(kbdim,klev),&
              in_w(kbdim,klev), in_rsi(kbdim,klev)
      REAL :: rv_old(kbdim,klev)

      ! -- Local gas compound tracers [# m-3] WHY ARE THESE GLOBAL?
      REAL :: zgso4(kbdim,klev),   &
              zghno3(kbdim,klev),  &
              zgnh3(kbdim,klev),   &
              zgocnv(kbdim,klev),  &
              zgocsv(kbdim,klev)
      
      ! For process rate diagnostics
      TYPE(FloatArray3d), POINTER :: autoconversion => NULL(), autoconversion50 => NULL(),      &
                                     autoconversion80 => NULL(),                                &
                                     accretion => NULL(), accretion50 => NULL(),                &
                                     accretion80 => NULL(),                                     &
                                     ACcoll => NULL(), APcoll => NULL(), AIcoll => NULL(),      &
                                     activation => NULL(), Icehom => NULL(), Icedep => NULL(),  &
                                     Iceimm => NULL(), Conda => NULL(), Condc => NULL(),        &
                                     Condp => NULL(), Condi => NULL()

      ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! As a local variable, it would be better to declare this info as REAL. Or then need to take care to initialize this,
      ! which is missing
      TYPE(Section) :: actd(kbdim,klev,ncld) ! Activated droplets - for interfacing with SALSA

      INTEGER :: jj,ii,kk,str,endd,nc,nb,nbloc,ndry,nwet,iwa,irim,icat

      INTEGER :: bin_numbers(4), bin_starts(4) ! Number of bins and the starting indices in the
                                               ! allSALSA array according to the phase indexing in classSection.f90
      INTEGER :: catnbins, lo, hi
      
      aero => NULL()
      cloud => NULL()
      precp => NULL()
      ice => NULL()
      liquid => NULL()
      frozen => NULL()
      
      ndry = spec%getNSpec(type="dry")
      nwet = spec%getNSpec(type="wet")  ! excludes rimed ice
      iwa = spec%getIndex("H2O")    ! water/unrimed ice
      irim = spec%getIndex("rime")  ! rimed ice; returns 0 if level<5 and rime not used

      ! The order here should follow the phase indentifier indexing!!
      bin_numbers(1:4) = [nbins,ncld,nprc,nice]
      bin_starts(1:4) = [iaero,icloud,iprecp,iice]

      ! Initialize the arrays for managing "old" values and tendencies
      CALL initialize_arrays(nwet,bin_numbers,npart,mpart,ntend,mtend)
      
      IF (.NOT. cgsto_initialized) &
           CALL initialize_coagstorage(nzp,nxp,nyp)

      ! Get the pointers to variable arrays
      CALL Diag%getData(1,press,name="press")
      CALL Diag%getData(1,tk,name="temp")
      CALL Diag%getData(1,rs,name="rsl")
      IF (level == 5) &
           CALL Diag%getData(1,rsi,name="rsi")
      CALL Diag%getData(1,pdn,name="dn")
      CALL Prog%getData(1,rv,name="rp")
      CALL Prog%getData(2,rt,name="rp")
      ! Note that maybe there could be a "getSeveral" method or something to make this less verbose?
      CALL Prog%getData(1,naerop,name="naero")
      CALL Prog%getData(2,naerot,name="naero")
      CALL Prog%getData(1,ncloudp,name="ncloud")
      CALL Prog%getData(2,ncloudt,name="ncloud")
      CALL Prog%getData(1,nprecpp,name="nprecp")
      CALL Prog%getData(2,nprecpt,name="nprecp")
      IF (level == 5) THEN
           CALL Prog%getData(1,nicep,name="nice")      
           CALL Prog%getData(2,nicet,name="nice")
      END IF
           
      CALL Prog%getData(1,maerop,name="maero")
      CALL Prog%getData(2,maerot,name="maero")
      CALL Prog%getData(1,mcloudp,name="mcloud")
      CALL Prog%getData(2,mcloudt,name="mcloud")
      CALL Prog%getData(1,mprecpp,name="mprecp")
      CALL Prog%getData(2,mprecpt,name="mprecp")
      IF (level == 5) THEN
         CALL Prog%getData(1,micep,name="mice")
         CALL Prog%getData(2,micet,name="mice")
         IF (ice_theta_dist) CALL Prog%getData(1,indefp,name="indef")
      END IF
         
      CALL Prog%getData(1,gaerop,name="gaero")
      CALL Prog%getData(2,gaerot,name="gaero")

      ! Get pointers to process rate output streams
      CALL Diag%getData(1,autoconversion,name="s_m_autoc")
      CALL Diag%getData(1,autoconversion50,name="s_m_autoc50")
      CALL Diag%getData(1,autoconversion80,name="s_m_autoc80")
      CALL Diag%getData(1,accretion,name="s_m_accr")
      CALL Diag%getData(1,accretion50,name="s_m_accr50")
      CALL Diag%getData(1,accretion80,name="s_m_accr80")
      CALL Diag%getData(1,ACcoll,name="s_m_ACcoll_dry")
      CALL Diag%getData(1,APcoll,name="s_m_APcoll_dry")
      CALL Diag%getData(1,activation,name="s_n_activ")
      IF (level == 5) THEN
         CALL Diag%getData(1,AIcoll,name="s_m_AIcoll_dry")
         CALL Diag%getData(1,Icehom,name="s_n_icehom")
         CALL Diag%getData(1,Icedep,name="s_n_icedep")
         CALL Diag%getData(1,Iceimm,name="s_n_iceimm")
         CALL Diag%getData(1,condi,name="s_m_condi")
      END IF
      CALL Diag%getData(1,conda,name="s_m_conda")
      CALL Diag%getData(1,condc,name="s_m_condc")
      CALL Diag%getData(1,condp,name="s_m_condp")
      
      in_p(:,:) = 0.; in_t(:,:) = 0.; in_rs(:,:) = 0.; in_rsi(:,:) = 0.; in_w(:,:) = 0.
      in_rv(:,:) = 0.; rv_old(:,:) = 0.
      
      ! Set the SALSA runtime config 
      CALL set_salsa_runtime(time,tstep,istp)
      
      ! Sets up the parallel region and privates nessecary variables. These are directly carried into subroutines as arguements for that routine
      !$OMP PARALLEL PRIVATE(aero,cloud,precp,ice,liquid,frozen,ntend,mtend) &
      !$OMP PRIVATE(nb,nc,icat,str,catnbins,nbloc,endd,zghno3,zgnh3,zgocnv,zgocsv,zgso4,lo,hi) &
      !$OMP FIRSTPRIVATE(in_p,in_t,in_rs,in_rsi,in_w,in_rv,rv_old,allCOAGcoe,allSALSA,npart,mpart)
      
      ! Sets privated pointers to point at the private copies of allSALSA
      lo = 1
      hi = nbins
      aero => allSALSA(:,:,lo:hi)
      
      lo = hi + 1
      hi = hi + ncld
      cloud => allSALSA(:,:,lo:hi)

      lo = hi + 1
      hi = hi + nprc
      precp => allSALSA(:,:,lo:hi)

      lo = hi + 1
      hi = hi + nice
      ice => allSALSA(:,:,lo:hi)

      lo = 1
      hi = nbins + ncld + nprc
      liquid => allSALSA(:,:,lo:hi)

      lo = nbins + ncld + nprc + 1
      hi = nbins + ncld + nprc + nice
      frozen => allSALSA(:,:,lo:hi)
      
      ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by substance density) 
      !$OMP DO SCHEDULE(GUIDED)
      DO jj = 3, nyp-2
         DO ii = 3, nxp-2
            DO kk = nzp-1, 2, -1
            
               ! Set inputs
               in_p(1,1) = press%d(kk,ii,jj)
               in_t(1,1) = tk%d(kk,ii,jj)
               in_rs(1,1) = rs%d(kk,ii,jj)
               IF (level == 5) &
                    in_rsi(1,1) = rsi%d(kk,ii,jj)
               in_w(1,1) = wp(kk,ii,jj)

               ! For initialization and spinup, limit the RH with the parameter rhlim (assign in namelist.salsa)
               IF ( lsfreeRH%state ) THEN
                  in_rv(1,1) = rv%d(kk,ii,jj)
               ELSE
                  in_rv(1,1) = MIN(rv%d(kk,ii,jj), rs%d(kk,ii,jj)*rhlim)
               END IF
               rv_old(1,1) = in_rv(1,1)
               
               ! Update volume concentrations
               ! ---------------------------------------------------------------------------------------------------
               ! Take a copies for the old values. Note that the order of the categories in these arrays is CRUCIAL.
               ! It follows the phase indexing in classSection. Use the derived type "old"
               ! Also make use of these when preparing SALSA inputs below. this is a bit verbose, but what are you gonna do...
               npart(1)%d = naerop%d(kk,ii,jj,:)
               npart(2)%d = ncloudp%d(kk,ii,jj,:)
               npart(3)%d = nprecpp%d(kk,ii,jj,:)
               IF (level == 5) &
                    npart(4)%d = nicep%d(kk,ii,jj,:)
               mpart(1)%d = maerop%d(kk,ii,jj,:)
               mpart(2)%d = mcloudp%d(kk,ii,jj,:)
               mpart(3)%d = mprecpp%d(kk,ii,jj,:)
               IF (level == 5) &
                    mpart(4)%d = micep%d(kk,ii,jj,:)
               
               ! Make arranged pointers (not copies!) to tendency arrays in the same manner
               ntend(1)%d => naerot%d(kk,ii,jj,:)
               ntend(2)%d => ncloudt%d(kk,ii,jj,:)
               ntend(3)%d => nprecpt%d(kk,ii,jj,:)
               IF (level == 5) &
                    ntend(4)%d => nicet%d(kk,ii,jj,:)
               mtend(1)%d => maerot%d(kk,ii,jj,:)
               mtend(2)%d => mcloudt%d(kk,ii,jj,:)
               mtend(3)%d => mprecpt%d(kk,ii,jj,:)
               IF (level == 5) &
                    mtend(4)%d => micet%d(kk,ii,jj,:)               
               
               ! Update SALSA input arrays
               DO nb = 1,ntotal
                  icat = allSALSA(1,1,nb)%phase  ! Phase indentifier, this should correspond to index in bin_starts and bin_numbers, mpart and npart, ntend and mtend
                  IF (icat == 4 .AND. level < 5) EXIT  ! skip ice if level < 5; since ice gets the last index in classSection, we can use EXIT
                  
                  nbloc = nb - bin_starts(icat) + 1
                  catnbins = bin_numbers(icat)
                  
                  ! Update volume concentrations. Note, rime handled separately below
                  IF (icat < 4) THEN
                     DO nc = 1,nwet
                        str = getMassIndex(catnbins,nbloc,nc) 
                        allSALSA(1,1,nb)%volc(nc) = mpart(icat)%d(str)*pdn%d(kk,ii,jj)/spec%rholiq(nc)                     
                     END DO
                     
                     IF (level == 5 .AND. ice_theta_dist) THEN
                        allSALSA(1,1,nb)%INdef = indefp%d(kk,ii,jj,nb)
                     END IF
                     
                  ELSE IF (icat == 4) THEN
                     DO nc = 1,nwet
                        str = getMassIndex(catnbins,nbloc,nc)
                        allSALSA(1,1,nb)%volc(nc) = mpart(icat)%d(str)*pdn%d(kk,ii,jj)/spec%rhoice(nc)
                     END DO
                     str = getMassIndex(catnbins,nbloc,irim)
                     allSALSA(1,1,nb)%volc(irim) = mpart(icat)%d(str)*pdn%d(kk,ii,jj)/spec%rhori

                  END IF
                  
                  ! Update number concentrations
                  allSALSA(1,1,nb)%numc = npart(icat)%d(nbloc)*pdn%d(kk,ii,jj)
                  
                  ! Update particle size and densities
                  !CALL allSALSA(1,1,nb)%updateDiameter(.TRUE.,type="all")
                  CALL allSALSA(1,1,nb)%updateRhomean()
                  
                  IF (allSALSA(1,1,nb)%numc > allSALSA(1,1,nb)%nlim) THEN
                     allSALSA(1,1,nb)%core = SUM(allSALSA(1,1,nb)%volc(1:ndry))/allSALSA(1,1,nb)%numc
                  ELSE
                     allSALSA(1,1,nb)%core = pi6*(allSALSA(1,1,nb)%dmid)**3
                  END IF
                  
               END DO
               
               ! If this is an initialization call, calculate the equilibrium particle
               If (initialize) CALL equilibration(kproma,kbdim,klev,init_rh, &
                                                  in_t, aero, .TRUE.         )
               
               ! Convert to #/m3
               zgso4(1,1) = gaerop%d(kk,ii,jj,1)*pdn%d(kk,ii,jj)
               zghno3(1,1) = gaerop%d(kk,ii,jj,2)*pdn%d(kk,ii,jj)
               zgnh3(1,1) = gaerop%d(kk,ii,jj,3)*pdn%d(kk,ii,jj)
               zgocnv(1,1) = gaerop%d(kk,ii,jj,4)*pdn%d(kk,ii,jj)
               zgocsv(1,1) = gaerop%d(kk,ii,jj,5)*pdn%d(kk,ii,jj)

               ! Reset process rate diagnostics
               CALL rateDiag%Reset()

               ! If reduced coagulation kernel update freq, and NOT update timestep,
               ! copy kernels from memory
               
               IF ( lscoag%mode == 2 .AND. .NOT. lcgupdt ) &
                     CALL fetch_coag(kk,ii,jj,allCOAGcoe ) 
               
               ! ***************************************!
               !                Run SALSA               !
               ! ***************************************!
               CALL salsa(kproma, kbdim,  klev,   krow,   in_p,   in_rv,  in_rs,  in_rsi,   &
                          in_t,   tstep,  zgso4,  zgocnv, zgocsv, zghno3, zgnh3,  actd,     &
                          in_w,   level,  allSALSA, aero, cloud,  precp,  ice, liquid,  allCOAGcoe )
               
               ! Update tendency arrays
               DO nb = 1,ntotal
                  icat = allSALSA(1,1,nb)%phase   ! Phase indentifier, this should correspond to index in npart and mpart
                  
                  IF (icat == 4 .AND. level < 5) EXIT  ! skip ice if level < 5; since ice gets the last index in classSection, we can use EXIT
                  
                  nbloc = nb - bin_starts(icat) + 1
                  catnbins = bin_numbers(icat)

                  ! Update mass tendencies
                  IF (icat < 4) THEN
                     DO nc = 1,nwet
                        str = getMassIndex(catnbins,nbloc,nc)
                        mtend(icat)%d(str) = mtend(icat)%d(str) + &
                             ( (allSALSA(1,1,nb)%volc(nc)*spec%rholiq(nc)/pdn%d(kk,ii,jj)) -   &
                               mpart(icat)%d(str) ) / tstep
                     END DO

                     IF (level == 5 .AND. ice_theta_dist) THEN
                        indefp%d(kk,ii,jj,nb) = allSALSA(1,1,nb)%INdef 
                     END IF
                      
                  ELSE IF (icat == 4) THEN
                     DO nc = 1,nwet
                        str = getMassIndex(catnbins,nbloc,nc)
                        mtend(icat)%d(str) = mtend(icat)%d(str) + &
                             ( (allSALSA(1,1,nb)%volc(nc)*spec%rhoice(nc)/pdn%d(kk,ii,jj)) -   &
                               mpart(icat)%d(str) ) / tstep
                     END DO
                     str = getMassIndex(catnbins,nbloc,irim)
                     mtend(icat)%d(str) = mtend(icat)%d(str) + &
                          ( (allSALSA(1,1,nb)%volc(irim)*spec%rhori/pdn%d(kk,ii,jj)) -     &
                            mpart(icat)%d(str) ) / tstep
                  END IF
                  
                  ! Update number tendencies
                  ntend(icat)%d(nbloc) = ntend(icat)%d(nbloc) + &
                       ( (allSALSA(1,1,nb)%numc/pdn%d(kk,ii,jj)) - npart(icat)%d(nbloc) )/tstep
                  
               END DO
               
               ! Activated droplets
               pa_nactd(kk,ii,jj,1:ncld) = actd(1,1,1:ncld)%numc/pdn%d(kk,ii,jj)
               
               DO nc = 1,iwa
                  ! Activated droplets (in case of cloud base activation)
                  str = getMassIndex(ncld,1,nc)
                  endd = getMassIndex(ncld,ncld,nc)
                  pa_vactd(kk,ii,jj,str:endd) = actd(1,1,1:ncld)%volc(nc)*spec%rholiq(nc)/pdn%d(kk,ii,jj)                  
               END DO

               IF (lscndgas) THEN
                  gaerot%d(kk,ii,jj,1) = gaerot%d(kk,ii,jj,1) + &
                                          ( (zgso4(1,1)/pdn%d(kk,ii,jj)) - gaerop%d(kk,ii,jj,1) )/tstep

                  gaerot%d(kk,ii,jj,2) = gaerot%d(kk,ii,jj,2) + &
                                          ( (zghno3(1,1)/pdn%d(kk,ii,jj)) - gaerop%d(kk,ii,jj,2) )/tstep

                  gaerot%d(kk,ii,jj,3) = gaerot%d(kk,ii,jj,3) + &
                                          ( (zgnh3(1,1)/pdn%d(kk,ii,jj)) - gaerop%d(kk,ii,jj,3) )/tstep

                  gaerot%d(kk,ii,jj,4) = gaerot%d(kk,ii,jj,4) + &
                                          ( (zgocnv(1,1)/pdn%d(kk,ii,jj)) - gaerop%d(kk,ii,jj,4) )/tstep

                  gaerot%d(kk,ii,jj,5) = gaerot%d(kk,ii,jj,5) + &
                                          ( (zgocsv(1,1)/pdn%d(kk,ii,jj)) - gaerop%d(kk,ii,jj,5) )/tstep
               END IF
               
               ! Tendency of water vapour mixing ratio 
               rt%d(kk,ii,jj) = rt%d(kk,ii,jj) + &
                  ( in_rv(1,1) - rv_old(1,1) )/tstep

               ! Store the process rate diagnostics for output
               autoconversion%d(kk,ii,jj) = rateDiag%Autoconversion%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               autoconversion80%d(kk,ii,jj) = rateDiag%Autoconversion80%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               autoconversion50%d(kk,ii,jj) = rateDiag%Autoconversion50%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               accretion%d(kk,ii,jj) = rateDiag%Accretion%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               accretion80%d(kk,ii,jj) = rateDiag%Accretion80%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               accretion50%d(kk,ii,jj) = rateDiag%Accretion50%volc(iwa)*spec%rhowa/pdn%d(kk,ii,jj)
               ACcoll%d(kk,ii,jj) = SUM(rateDiag%ACcoll%volc(1:ndry)*spec%rholiq(1:ndry))/pdn%d(kk,ii,jj)
               APcoll%d(kk,ii,jj) = SUM(rateDiag%APcoll%volc(1:ndry)*spec%rholiq(1:ndry))/pdn%d(kk,ii,jj)
               activation%d(kk,ii,jj) = rateDiag%Activation%numc/pdn%d(kk,ii,jj)
               conda%d(kk,ii,jj) = rateDiag%Cond_a%volc(iwa)/pdn%d(kk,ii,jj)
               condc%d(kk,ii,jj) = rateDiag%Cond_c%volc(iwa)/pdn%d(kk,ii,jj)
               condp%d(kk,ii,jj) = rateDiag%Cond_p%volc(iwa)/pdn%d(kk,ii,jj)
               IF (level == 5) THEN
                  AIcoll%d(kk,ii,jj) = SUM(rateDiag%AIcoll%volc(1:ndry)*spec%rhoice(1:ndry))/pdn%d(kk,ii,jj)
                  Icehom%d(kk,ii,jj) = rateDiag%Ice_hom%numc/pdn%d(kk,ii,jj)
                  Icedep%d(kk,ii,jj) = rateDiag%Ice_dep%numc/pdn%d(kk,ii,jj)
                  Iceimm%d(kk,ii,jj) = rateDiag%Ice_imm%numc/pdn%d(kk,ii,jj)
                  Condi%d(kk,ii,jj) = ( rateDiag%Cond_i%volc(iwa)*spec%rhoic    +  &
                                        rateDiag%Cond_i%volc(irim)*spec%rhori ) /  &
                                        pdn%d(kk,ii,jj)
               END IF
               
               ! If reduced coagulation kernel update frequency and is update timestep, store the updated kernels
               IF (lscoag%mode == 2 .AND. lcgupdt) &
                    CALL store_coag(kk,ii,jj,allCOAGcoe)
               
            END DO !kk
         END DO ! ii
      END DO ! jj
      !$OMP END DO
      !$OMP END PARALLEL
      
      ! Nullify pointers and clean up
      DO icat = 1,4
         ntend(icat)%d => NULL()
         mtend(icat)%d => NULL()
      END DO
      
      press => NULL(); tk => NULL(); rv => NULL()
      rt => NULL(); rs => NULL(); rsi => NULL()
      pdn => NULL()
      naerop => NULL(); ncloudp => NULL()
      nprecpp => NULL(); nicep => NULL()
      naerot => NULL(); ncloudt => NULL()
      nprecpt => NULL(); nicet => NULL()
      
      maerop => NULL(); mcloudp => NULL()
      mprecpp => NULL(); micep => NULL()
      maerot => NULL(); mcloudt => NULL()
      mprecpt => NULL(); micet => NULL()
      
      gaerop => NULL(); gaerot => NULL()
      indefp => NULL()
      
   END SUBROUTINE run_SALSA

   ! 
   !---------------------------------------------------------------
   ! SET_SALSA_RUNTIME
   ! Set the master process %state:s based on the values of %switch and %delay.
   ! Added some new time-dependent switches
   ! 
   ! Juha Tonttila, FMI, 2014;2020
   ! 
   SUBROUTINE set_SALSA_runtime(time, tstep, istp)
     IMPLICIT NONE
     REAL, INTENT(in) :: time, tstep
     INTEGER, INTENT(in) :: istp
     INTEGER :: i

     ! All the non-interface parameters used below are found in mo_submctl
     
     DO i = 1,Nmaster
        IF( lsmaster(i)%switch .AND. time >= lsmaster(i)%delay ) lsmaster(i)%state = .TRUE.
     END DO

     ! Some other switches
     ! Constraining RH in condensation
     IF ( lsfreeRH%switch .AND. time > lsfreeRH%delay ) lsfreeRH%state = .TRUE.
     ! Constraining contact angle in ice nucleation
     IF ( lsFreeTheta%switch .AND. time > lsFreeTheta%delay) lsFreeTheta%state = .TRUE.
     ! Determine coagulation kernel update status if low freq updating is active
     lcgupdt = ( lscoag%state .AND. lscoag%mode == 2 .AND.    &
                 (MOD(time,cgintvl) < tstep .OR. istp <= 1)   ) ! Making sure the kernels
                                                                ! are calculated on the first
                                                                ! timestep of HISTORY runs
     
   END SUBROUTINE set_SALSA_runtime

   ! 
   ! -----------------------------------------------------------
   ! INITIALIZE_ARRAYS
   ! Initializes the storage for arranging SALSA interface data
   ! 
   SUBROUTINE initialize_arrays(nwet,bins,npart,mpart,ntend,mtend)
     TYPE(old), INTENT(inout) :: npart(4),mpart(4)
     TYPE(FloatArray1d), INTENT(inout) :: ntend(4),mtend(4)
     INTEGER, INTENT(in) :: nwet
     INTEGER, INTENT(in) :: bins(4)     
     INTEGER :: icat
     DO icat = 1,4
        ntend(icat) = FloatArray1d("ntend")   ! Dunno if this mess is really usefull...-Juha
        mtend(icat) = FloatArray1d("mtend")
        npart(icat) = old(bins(icat))
        IF (icat < 4) THEN
           mpart(icat) = old(bins(icat)*nwet)
        ELSE
           mpart(icat) = old(bins(icat)*(nwet+1))
        END IF
     END DO
   END SUBROUTINE initialize_arrays

   ! 
   ! ----------------------------------------------------
   ! Constructor for the datatype storing "old" values
   ! used to back up the tendencies.
   ! 
   FUNCTION old_cnstr(n)
     TYPE(old) :: old_cnstr
     INTEGER, INTENT(in) :: n
     ALLOCATE(old_cnstr%d(n))
     old_cnstr%d(:) = 0.
   END FUNCTION old_cnstr

   ! 
   ! ----------------------------------------------------
   ! Initializes the coagulation kernel storage 
   ! 
   SUBROUTINE initialize_coagstorage(nzp,nxp,nyp)
     INTEGER, INTENT(in) :: nxp,nyp,nzp
     
     IF (lscgaa) THEN
        ALLOCATE(sto_aa(nzp,nxp,nyp,nbins,nbins))
        sto_aa = 0.
     END IF
     IF (lscgcc) THEN
        ALLOCATE(sto_cc(nzp,nxp,nyp,ncld,ncld))
        sto_cc = 0.
     END IF
     IF (lscgpp) THEN
        ALLOCATE(sto_pp(nzp,nxp,nyp,nprc,nprc))
        sto_pp = 0.
     END IF
     IF (lscgii) THEN
        ALLOCATE(sto_ii(nzp,nxp,nyp,nice,nice))
        sto_ii = 0.
     END IF
        
     IF (lscgca) THEN
        ALLOCATE(sto_ca(nzp,nxp,nyp,nbins,ncld))
        sto_ca = 0.
     END IF
     IF (lscgpa) THEN
        ALLOCATE(sto_pa(nzp,nxp,nyp,nbins,nprc))
        sto_pa = 0.
     END IF
     IF (lscgia) THEN
        ALLOCATE(sto_ia(nzp,nxp,nyp,nbins,nice))
        sto_ia = 0.
     END IF

     IF (lscgpc) THEN
        ALLOCATE(sto_pc(nzp,nxp,nyp,ncld,nprc))
        sto_pc = 0.
     END IF
     IF (lscgic) THEN
        ALLOCATE(sto_ic(nzp,nxp,nyp,ncld,nice))
        sto_ic = 0.
     END IF
        
     IF (lscgip) THEN
        ALLOCATE(sto_ip(nzp,nxp,nyp,nice,nprc))
        sto_ip = 0.     
     END IF

     cgsto_initialized = .TRUE.
     
   END SUBROUTINE initialize_coagstorage
   
   ! 
   ! ------------------------------------------------------
   ! Store updated coagulation kernels
   ! 
   SUBROUTINE store_coag(kk,ii,jj,allCOAGcoe)
                         
     TYPE(CoagCoe), INTENT(in) :: allCOAGcoe(:)
     INTEGER, INTENT(in) :: ii,jj,kk
     
     IF (lscgaa) sto_aa(kk,ii,jj,:,:) = allCOAGcoe(1)%zccaa(1,1,:,:)
     
     IF (lscgcc) sto_cc(kk,ii,jj,:,:) = allCOAGcoe(2)%zcccc(1,1,:,:)
     IF (lscgpp) sto_pp(kk,ii,jj,:,:) = allCOAGcoe(3)%zccpp(1,1,:,:)
     IF (lscgii) sto_ii(kk,ii,jj,:,:) = allCOAGcoe(4)%zccii(1,1,:,:)

     IF (lscgca) sto_ca(kk,ii,jj,:,:) = allCOAGcoe(5)%zccca(1,1,:,:)
     IF (lscgpa) sto_pa(kk,ii,jj,:,:) = allCOAGcoe(6)%zccpa(1,1,:,:)
     IF (lscgia) sto_ia(kk,ii,jj,:,:) = allCOAGcoe(7)%zccia(1,1,:,:)

     IF (lscgpc) sto_pc(kk,ii,jj,:,:) = allCOAGcoe(8)%zccpc(1,1,:,:)
     IF (lscgic) sto_ic(kk,ii,jj,:,:) = allCOAGcoe(9)%zccic(1,1,:,:)

     IF (lscgip) sto_ip(kk,ii,jj,:,:) = allCOAGcoe(10)%zccip(1,1,:,:)
     
   END SUBROUTINE store_coag

   ! 
   ! -----------------------------------------------------
   ! Copy coagulation kernels from storage for current timestep,
   ! i .e. avoid calculating new ones
   ! 
   SUBROUTINE fetch_coag(kk,ii,jj,allCOAGcoe )
   
     TYPE(CoagCoe), INTENT(inout) :: allCOAGcoe(:)
     INTEGER, INTENT(in) :: ii,jj,kk
     
     IF (lscgaa) allCOAGcoe(1)%zccaa(1,1,:,:) = sto_aa(kk,ii,jj,:,:)
     IF (lscgcc) allCOAGcoe(2)%zcccc(1,1,:,:) = sto_cc(kk,ii,jj,:,:) 
     IF (lscgpp) allCOAGcoe(3)%zccpp(1,1,:,:) = sto_pp(kk,ii,jj,:,:)
     IF (lscgii) allCOAGcoe(4)%zccii(1,1,:,:) = sto_ii(kk,ii,jj,:,:)

     IF (lscgca) allCOAGcoe(5)%zccca(1,1,:,:) = sto_ca(kk,ii,jj,:,:)
     IF (lscgpa) allCOAGcoe(6)%zccpa(1,1,:,:) = sto_pa(kk,ii,jj,:,:) 
     IF (lscgia) allCOAGcoe(7)%zccia(1,1,:,:) = sto_ia(kk,ii,jj,:,:) 

     IF (lscgpc) allCOAGcoe(8)%zccpc(1,1,:,:) = sto_pc(kk,ii,jj,:,:) 
     IF (lscgic) allCOAGcoe(9)%zccic(1,1,:,:) = sto_ic(kk,ii,jj,:,:) 

     IF (lscgip) allCOAGcoe(10)%zccip(1,1,:,:) = sto_ip(kk,ii,jj,:,:)
     
   END SUBROUTINE fetch_coag
     
END MODULE mo_salsa_driver
