MODULE mo_salsa_driver
  USE classSection, ONLY : Section
  USE util, ONLY : getMassIndex !!! IS it good to import this here??? The function is anyway handy here too.
  USE mo_salsa_types, ONLY : aero, cloud, precp, ice, allSALSA, iaero, icloud, iprecp, iice
  USE mo_submctl
  USE classFieldArray, ONLY : FieldArray
  USE mo_structured_datatypes
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

   ! JT: Variables from SALSA
   ! --------------------------------------------
   ! grid points for SALSA
   INTEGER, PARAMETER :: kproma = 1
   INTEGER, PARAMETER :: kbdim = 1
   INTEGER, PARAMETER :: klev = 1
   INTEGER, PARAMETER :: krow = 1

   REAL, PARAMETER    :: init_rh(kbdim,klev) = 0.3

   ! -- Local gas compound tracers [# m-3]
   REAL :: zgso4(kbdim,klev),   &
           zghno3(kbdim,klev),  &
           zgnh3(kbdim,klev),   &
           zgocnv(kbdim,klev),  &
           zgocsv(kbdim,klev)

 ! --------------------------------------------

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
   !
   ! 05/2016 Juha: This routine is still pretty much in its original shape.
   !               It's dumb as a mule and twice as ugly, so implementation of
   !               an improved solution is necessary sooner or later.
   !
   ! Juha Tonttila, FMI, 2014
   ! Jaakko Ahola, FMI, 2016
   !
  SUBROUTINE run_SALSA(Diag, Prog, nzp, nxp, nyp, ns, wp,    &
                       pa_nactd, pa_vactd, tstep, time,      &
                       level, initialize                     )

      USE mo_salsa, ONLY : salsa
      USE mo_salsa_properties, ONLY  : equilibration
      IMPLICIT NONE
      
      TYPE(FieldArray), INTENT(inout) :: Diag, Prog
      
      INTEGER, INTENT(in) :: nzp,nxp,nyp,ns                       ! dimensions: x,y,z,number of chemical species
      REAL, INTENT(in)    :: tstep                                ! Model timestep length
      REAL, INTENT(in)    :: time
      LOGICAL, INTENT(in) :: initialize                      
      REAL, INTENT(in)    :: wp(nzp,nxp,nyp)
      REAL, INTENT(out)   :: pa_vactd(nzp,nxp,nyp,ns*ncld) ! mass concentrations of newly activated droplets for calculating the
                                                           ! actual tendency due to new droplet formation.
      REAL, INTENT(out)   :: pa_nactd(nzp,nxp,nyp,ncld)   ! Same for number concentration
      INTEGER, INTENT(in) :: level                         ! thermodynamical level

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

                                     gaerop => NULL(), gaerot => NULL()
      TYPE(FloatArray1d) :: npart(4), mpart(4)  ! To store the old values
      TYPE(FloatArray1d) :: ntend(4), mtend(4)  ! Arranged pointers to tendencies      
      
      REAL :: in_p(kbdim,klev), in_t(kbdim,klev), in_rv(kbdim,klev), in_rs(kbdim,klev),&
              in_w(kbdim,klev), in_rsi(kbdim,klev)
      REAL :: rv_old(kbdim,klev)

      REAL :: mrim(nice), mprist(nice), mrim_old(nice), mprist_old(nice)

      TYPE(Section) :: actd(kbdim,klev,ncld) ! Activated droplets - for interfacing with SALSA

      INTEGER :: jj,ii,kk,str,end,nc,nb,nbloc,ndry,nwet,iwa,irim,icat

      INTEGER :: bin_numbers(4), bin_starts(4) ! Number of bins and the starting indices in the
                                               ! allSALSA array according to the phase indexing in classSection.f90
      INTEGER :: catnbins
      
      ndry = spec%getNSpec(type="dry")
      nwet = spec%getNSpec(type="wet")  ! excludes rimed ice
      iwa = spec%getIndex("H2O")    ! water/unrimed ice
      irim = spec%getIndex("rime")  ! rimed ice; returns 0 if level<5 and rime not used

      ! The order here should follow the phase indentifier indexing!!
      bin_numbers(1:4) = [nbins,ncld,nprc,nice]
      bin_starts(1:4) = [iaero,icloud,iprecp,iice]
      
      ! Get the pointers to variable arrays
      CALL Diag%getData(1,press,name="press")
      CALL Diag%getData(1,tk,name="temp")
      CALL Diag%getData(1,rs,name="rsl")
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
      CALL Prog%getData(1,nicep,name="nice")
      CALL Prog%getData(2,nicet,name="nice")

      CALL Prog%getData(1,maerop,name="maero")
      CALL Prog%getData(2,maerot,name="maero")
      CALL Prog%getData(1,mcloudp,name="mcloud")
      CALL Prog%getData(2,mcloudt,name="mcloud")
      CALL Prog%getData(1,mprecpp,name="mprecp")
      CALL Prog%getData(2,mprecpt,name="mprecp") 
      CALL Prog%getData(1,micep,name="mice")
      CALL Prog%getData(2,micet,name="mice")

      CALL Prog%getData(1,gaerop,name="gaero")
      CALL Prog%getData(2,gaerot,name="gaero")

      in_p(:,:) = 0.; in_t(:,:) = 0.; in_rs(:,:) = 0.; in_rsi(:,:) = 0.; in_w(:,:) = 0.
      in_rv(:,:) = 0.; rv_old(:,:) = 0.
      
      ! Set the SALSA runtime config 
      CALL set_salsa_runtime(time)

      ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by substance density)
      DO jj = 3, nyp-2
         DO ii = 3, nxp-2
            DO kk = nzp-1, 2, -1

               ! Set inputs
               in_p(1,1) = press%d(kk,ii,jj)
               in_t(1,1) = tk%d(kk,ii,jj)
               in_rs(1,1) = rs%d(kk,ii,jj)
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
               ! It follows the phase indexing in classSection.
               ! Also make use of these when preparing SALSA inputs below.
               npart(1) = FloatArray1d(naerop%d(kk,ii,jj,:),store=.TRUE.)
               npart(2) = FloatArray1d(ncloudp%d(kk,ii,jj,:),store=.TRUE.)
               npart(3) = FloatArray1d(nprecpp%d(kk,ii,jj,:),store=.TRUE.)
               npart(4) = FloatArray1d(nicep%d(kk,ii,jj,:),store=.TRUE.)
               mpart(1) = FloatArray1d(maerop%d(kk,ii,jj,:),store=.TRUE.)
               mpart(2) = FloatArray1d(mcloudp%d(kk,ii,jj,:),store=.TRUE.)
               mpart(3) = FloatArray1d(mprecpp%d(kk,ii,jj,:),store=.TRUE.)
               mpart(4) = FloatArray1d(micep%d(kk,ii,jj,:),store=.TRUE.)

               ! Make arranged pointers (not copies!) to tendency arrays in the same manner
               ntend(1) = FloatArray1d(naerot%d(kk,ii,jj,:),store=.FALSE.)
               ntend(2) = FloatArray1d(ncloudt%d(kk,ii,jj,:),store=.FALSE.)
               ntend(3) = FloatArray1d(nprecpt%d(kk,ii,jj,:),store=.FALSE.)
               ntend(4) = FloatArray1d(nicet%d(kk,ii,jj,:),store=.FALSE.)
               mtend(1) = FloatArray1d(maerot%d(kk,ii,jj,:),store=.FALSE.)
               mtend(2) = FloatArray1d(mcloudt%d(kk,ii,jj,:),store=.FALSE.)
               mtend(3) = FloatArray1d(mprecpt%d(kk,ii,jj,:),store=.FALSE.)
               mtend(4) = FloatArray1d(micet%d(kk,ii,jj,:),store=.FALSE.)               
               
               ! Update SALSA input arrays
               DO nb = 1,ntotal
                  icat = allSALSA(1,1,nb)%phase   ! Phase indentifier, this should correspond to index in bin_starts and bin_numbers, mpart and npart, ntend and mtend
                  IF (icat == 4 .AND. level < 5) EXIT  ! skip ice if level < 5; since ice gets the last index in classSection, we can use EXIT
                  
                  nbloc = nb - bin_starts(icat) + 1
                  catnbins = bin_numbers(icat)
                  
                  ! Update volume concentrations. Note, rime handled separately below
                  IF (icat < 4) THEN
                     DO nc = 1,nwet
                        str = getMassIndex(catnbins,nbloc,nc) 
                        allSALSA(1,1,nb)%volc(nc) = mpart(icat)%d(str)*pdn%d(kk,ii,jj)/spec%rholiq(nc)                     
                     END DO
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
                  CALL allSALSA(1,1,nb)%updateDiameter(.TRUE.,type="all")
                  CALL allSALSA(1,1,nb)%updateRhomean()

                  IF (allSALSA(1,1,nb)%numc > allSALSA(1,1,nb)%nlim) THEN
                     allSALSA(1,1,nb)%core = SUM(allSALSA(1,1,nb)%volc(1:ndry))/allSALSA(1,1,nb)%numc
                  ELSE
                     allSALSA(1,1,nb)%core = pi6*(allSALSA(1,1,nb)%dmid)**3
                  END IF 
                  
               END DO
                       
               ! If this is an initialization call, calculate the equilibrium particle
               If (initialize) CALL equilibration(kproma,kbdim,klev,   &
                                                  init_rh,in_t,.TRUE.)

               ! Convert to #/m3
               zgso4(1,1) = gaerop%d(kk,ii,jj,1)*pdn%d(kk,ii,jj)
               zghno3(1,1) = gaerop%d(kk,ii,jj,2)*pdn%d(kk,ii,jj)
               zgnh3(1,1) = gaerop%d(kk,ii,jj,3)*pdn%d(kk,ii,jj)
               zgocnv(1,1) = gaerop%d(kk,ii,jj,4)*pdn%d(kk,ii,jj)
               zgocsv(1,1) = gaerop%d(kk,ii,jj,5)*pdn%d(kk,ii,jj)
               
               ! ***************************************!
               !                Run SALSA               !
               ! ***************************************!
               CALL salsa(kproma, kbdim,  klev,   krow,     &
                          in_p,   in_rv,  in_rs,  in_rsi,   &
                          in_t,   tstep,  zgso4,  zgocnv,   &
                          zgocsv, zghno3, zgnh3,  actd,     &
                          in_w,   level                     )
               
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
                  end = getMassIndex(ncld,ncld,nc)
                  pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(nc)*spec%rholiq(nc)/pdn%d(kk,ii,jj)                  
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

               !IF (kk == 40 .AND. ii == 3 .AND. jj == 23) THEN
               !   WRITE(*,*) 'DRIVER, arvo vanha', rv_old(1,1)
               !   WRITE(*,*) 'DRIVER, arvo uusi', in_rv(1,1)
               !END IF
                  
               ! Tendency of water vapour mixing ratio 
               rt%d(kk,ii,jj) = rt%d(kk,ii,jj) + &
                  ( in_rv(1,1) - rv_old(1,1) )/tstep

               !IF (kk == 40 .AND. ii == 3 .AND. jj == 23) THEN
               !   WRITE(*,*) 'DRIVER, tend ', rt%d(kk,ii,jj)
               !END IF
               
            END DO !kk
         END DO ! ii
      END DO ! jj

      ! Nullify pointers
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

      
   END SUBROUTINE run_SALSA

   !
   !---------------------------------------------------------------
   ! SET_SALSA_RUNTIME
   ! Set the master process %state:s based on the values of %switch and %delay
   !
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE set_SALSA_runtime(time)
     USE mo_submctl, ONLY : Nmaster, lsmaster, lsfreeRH
     IMPLICIT NONE
     REAL, INTENT(in) :: time
     INTEGER :: i

     DO i = 1,Nmaster
        IF( lsmaster(i)%switch .AND. time > lsmaster(i)%delay ) lsmaster(i)%state = .TRUE.
     END DO

     ! Some other switches
     IF ( lsfreeRH%switch .AND. time > lsfreeRH%delay ) lsfreeRH%state = .TRUE.
     
   END SUBROUTINE set_SALSA_runtime


END MODULE mo_salsa_driver
