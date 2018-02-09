MODULE mo_salsa_driver
  USE classSection
  USE util, ONLY : getMassIndex !!! IS it good to import this here??? The function is anyway handy here too.
  USE mo_submctl
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

   ! -- Local hydrometeor properties. All the setup and pointer associations will be done in mo_aero_init
   TYPE(Section), ALLOCATABLE, TARGET :: allSALSA(:,:,:)   ! Parent array holding all particle and hydrometeor types consecutively 

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
   SUBROUTINE run_SALSA(pnx, pny, pnz, n4, press, tk, rv, rt, rs, rsi, wp, pdn,   &
                        pa_naerop,  pa_naerot,  pa_maerop,  pa_maerot,   &
                        pa_ncloudp, pa_ncloudt, pa_mcloudp, pa_mcloudt,  &
                        pa_nprecpp, pa_nprecpt, pa_mprecpp, pa_mprecpt,  &
                        pa_nicep,   pa_nicet,   pa_micep,   pa_micet,    &
                        pa_nsnowp,  pa_nsnowt,  pa_msnowp,  pa_msnowt,   &
                        pa_nactd,   pa_vactd,   pa_gaerop,  pa_gaerot,   &
                        prunmode, tstep, level)

      USE mo_submctl, ONLY : nbins,ncld,nprc,pi6,          &
                             nice,nsnw,ntotal,                    &
                             !rhoic,rhosn,                  &
                             !rhowa, rhosu, rhobc, rhooc,   &
                             !rhono, rhonh, rhoss, rhodu,   &
                             !rhlim, lscndgas, nlim, prlim, &
                             aero, cloud, precp, ice, snow  ! The specific particle types are here for convenience, where needed.
      USE mo_salsa, ONLY : salsa
      USE mo_salsa_properties, ONLY  : equilibration
      IMPLICIT NONE

      INTEGER, INTENT(in) :: pnx,pny,pnz,n4                       ! dimensions: x,y,z,number of chemical species
      REAL, INTENT(in)    :: tstep                                ! Model timestep length

      REAL, INTENT(in)    :: press(pnz,pnx,pny), &            ! Pressure (Pa)
                             tk(pnz,pnx,pny),    &            ! Temperature (K)
                             rv(pnz,pnx,pny),    &            ! Water vapor mixing ratio
                             rs(pnz,pnx,pny),    &            ! Water vapour saturation mixing ratio
                             rsi(pnz,pnx,pny),   &            ! water vapour sat mix rat over ice
                             wp(pnz,pnx,pny)                  ! Vertical velocity (m s-1)

      REAL, INTENT(in)    :: pdn(pnz,pnx,pny)             ! Air density (for normalizing concentrations)

      REAL, INTENT(in)    :: pa_naerop(pnz,pnx,pny,nbins),        & ! aerosol number concentration (# kg-1)
                             pa_maerop(pnz,pnx,pny,n4*nbins),     & ! aerosol volume concentration (m3 kg-1)
                             pa_ncloudp(pnz,pnx,pny,ncld),        & ! Cloud droplet number concentration (# kg-1)
                             pa_mcloudp(pnz,pnx,pny,n4*ncld),     & ! Cloud droplet volume concentration (m3 kg-1)
                             pa_nprecpp(pnz,pnx,pny,nprc),        & ! Rain drop number concentration (# kg-1)
                             pa_mprecpp(pnz,pnx,pny,n4*nprc),     & ! Rain drop volume concentration (m3 kg-1)
                             pa_nicep(pnz,pnx,pny,nice),          & ! ice number concentration (# kg-1)
                             pa_micep(pnz,pnx,pny,n4*nice),       & ! ice volume concentration (m3 kg-1)
                             pa_nsnowp(pnz,pnx,pny,nsnw),         & ! snow precipitation number concentration (# kg-1)
                             pa_msnowp(pnz,pnx,pny,n4*nsnw)           ! snow precipitation volume concentration (m3 kg-1)

      REAL, INTENT(in)    :: pa_gaerop(pnz,pnx,pny,5)         ! Gaseous tracers [# kg]

      INTEGER, INTENT(in) :: prunmode                      ! 1: Initialization call
                                                           ! 2: Spinup period call
                                                           ! 3: Regular runtime call'
      INTEGER, INTENT(in) :: level                         ! thermodynamical level

      REAL, INTENT(inout) :: pa_naerot(pnz,pnx,pny,nbins),      & ! Aerosol number tendency
                             pa_maerot(pnz,pnx,pny,n4*nbins),   & ! Aerosol volume tendency
                             pa_ncloudt(pnz,pnx,pny,ncld),      & ! Cloud droplet number tendency
                             pa_mcloudt(pnz,pnx,pny,n4*ncld),   & ! Cloud droplet volume tendency
                             pa_nprecpt(pnz,pnx,pny,nprc),      & ! Rain drop number tendency
                             pa_mprecpt(pnz,pnx,pny,n4*nprc),   &  ! Rain drop volume tendency
                             pa_nicet(pnz,pnx,pny,nice),        & ! Ice particle number tendency
                             pa_micet(pnz,pnx,pny,n4*nice),     & ! Ice particle volume tendency
                             pa_nsnowt(pnz,pnx,pny,nsnw),       & ! snow flake number tendency
                             pa_msnowt(pnz,pnx,pny,n4*nsnw)         ! snow flake volume tendecy

      REAL, INTENT(inout) :: pa_gaerot(pnz,pnx,pny,5)         ! Gaseous tracer tendency
      REAL, INTENT(inout) :: rt(pnz,pnx,pny)                  ! Water vapour tendency

      REAL, INTENT(out) :: pa_vactd(pnz,pnx,pny,n4*ncld) ! Volume concentrations of newly activated droplets for calculating the
                                                           ! actual tendency due to new droplet formation.
      REAL, INTENT(out) :: pa_nactd(pnz,pnx,pny,ncld)   ! Same for number concentration

      TYPE(Section) :: actd(kbdim,klev,ncld) ! Activated droplets - for interfacing with SALSA

      ! Helper arrays for calculating the rates of change
      TYPE(Section) :: aero_old(1,1,nbins), cloud_old(1,1,ncld), precp_old(1,1,nprc), ice_old(1,1,nice), snow_old(1,1,nsnw)

      INTEGER :: jj,ii,kk,ss,str,end, nc,vc, ndry
      REAL :: in_p(kbdim,klev), in_t(kbdim,klev), in_rv(kbdim,klev), in_rs(kbdim,klev),&
              in_w(kbdim,klev), in_rsi(kbdim,klev), in_pdn(kbdim,klev)
      REAL :: rv_old(kbdim,klev)

      actd(:,:,:) = Section()
      aero_old(:,:,:) = Section()
      cloud_old(:,:,:) = Section()
      precp_old(:,:,:) = Section()
      ice_old(:,:,:) = Section()
      snow_old(:,:,:) = Section()

      ! Not needed?
      ! Number is always set, but mass can be uninitialized
      !DO ss = 1, spec%getNSpec()
      !   actd(:,:,:)%volc(ss) = 0.
      !   aero(:,:,:)%volc(ss) = 0.
      !   cloud(:,:,:)%volc(ss) = 0.
      !   precp(:,:,:)%volc(ss) = 0.
      !   ice(:,:,:)%volc(ss) = 0.
      !   snow(:,:,:)%volc(ss) = 0.
      !   aero_old(:,:,:)%volc(ss) = 0.
      !   cloud_old(:,:,:)%volc(ss) = 0.
      !   precp_old(:,:,:)%volc(ss) = 0.
      !   ice_old(:,:,:)%volc(ss) = 0.
      !   snow_old(:,:,:)%volc(ss) = 0.
      !END DO

      ! Set the SALSA runtime config 
      CALL set_salsa_runtime(prunmode)

      ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by substance density)
      DO jj = 3, pny-2
         DO ii = 3, pnx-2
            DO kk = pnz-1, 2, -1

               ! Set inputs
               in_p(1,1) = press(kk,ii,jj)
               in_t(1,1) = tk(kk,ii,jj)
               in_rs(1,1) = rs(kk,ii,jj)
               in_rsi(1,1) = rsi(kk,ii,jj)
               in_w(1,1) = wp(kk,ii,jj)

               ! For initialization and spinup, limit the RH with the parameter rhlim (assign in namelist.salsa)
               IF (prunmode < 3) THEN
                  in_rv(1,1) = MIN(rv(kk,ii,jj), rs(kk,ii,jj)*rhlim)
               ELSE
                  in_rv(1,1) = rv(kk,ii,jj)
               END IF
               rv_old(1,1) = in_rv(1,1)
       
               ! Update volume concentrations
               ! ---------------------------------------------------------------------------------------------------
               DO nc = 1,spec%getNSpec()
                  str = getMassIndex(nbins,1,nc)
                  end = getMassIndex(nbins,nbins,nc)
                  aero(1,1,1:nbins)%volc(nc) =  pa_maerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/spec%rholiq(nc)

                  str = getMassIndex(ncld,1,nc)
                  end = getMassIndex(ncld,ncld,nc)
                  cloud(1,1,1:ncld)%volc(nc) = pa_mcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/spec%rholiq(nc)

                  str = getMassIndex(nprc,1,nc)
                  end = getMassIndex(nprc,nprc,nc)
                  precp(1,1,1:nprc)%volc(nc) = pa_mprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/spec%rholiq(nc)
                  
                  str = getMassIndex(nice,1,nc)
                  end = getMassIndex(nice,nice,nc)
                  ice(1,1,1:nice)%volc(nc) = pa_micep(kk,ii,jj,str:end)*pdn(kk,ii,jj)/spec%rhoice(nc)

                  str = getMassIndex(nsnw,1,nc)
                  end = getMassIndex(nsnw,nsnw,nc)
                  snow(1,1,1:nsnw)%volc(nc) = pa_msnowp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/spec%rhosnow(nc)
               END DO
               ! -------------------------------
               
               ! Update number concentrations and particle sizes
               ! ----------------------------------------------------------------------------
               aero(1,1,1:nbins)%numc = pa_naerop(kk,ii,jj,1:nbins)*pdn(kk,ii,jj)
               cloud(1,1,1:ncld)%numc = pa_ncloudp(kk,ii,jj,1:ncld)*pdn(kk,ii,jj)
               precp(1,1,1:nprc)%numc = pa_nprecpp(kk,ii,jj,1:nprc)*pdn(kk,ii,jj)
               IF (level > 4) THEN
                  ice(1,1,1:nice)%numc = pa_nicep(kk,ii,jj,1:nice)*pdn(kk,ii,jj)
                  snow(1,1,1:nsnw)%numc = pa_nsnowp(kk,ii,jj,1:nsnw)*pdn(kk,ii,jj)
               END IF

               ! Need the number of dry species below
               ndry = spec%getNSpec(type="dry")
               DO ss = 1,ntotal
                  ! For simplicity use prlim here? Doesn't matter so much in this part whether it's prlim or nlim?
                  IF (allSALSA(1,1,ss)%numc > prlim) THEN
                     allSALSA(1,1,ss)%core = SUM(allSALSA(1,1,ss)%volc(1:ndry))/allSALSA(1,1,ss)%numc
                     allSALSA(1,1,ss)%dwet = (SUM(allSALSA(1,1,ss)%volc(:))/allSALSA(1,1,ss)%numc/pi6 )**(1./3.)
                  ELSE
                     allSALSA(1,1,ss)%dwet = allSALSA(1,1,ss)%dmid
                     allSALSA(1,1,ss)%core = pi6*(allSALSA(1,1,ss)%dmid)**3
                  END IF 
               END DO
               ! --------------------------------
               
               ! Take a copy of current concentrations to convert to tendencies after SALSA call
               aero_old = aero; cloud_old = cloud; precp_old = precp
               IF (level > 4) THEN
                  ice_old = ice; snow_old = snow
               END IF
        
               ! If this is an initialization call, calculate the equilibrium particle
               If (prunmode == 1) CALL equilibration(kproma,kbdim,klev,   &
                                                     init_rh,in_t,.TRUE.)

               ! Convert to #/m3
               zgso4(1,1) = pa_gaerop(kk,ii,jj,1)*pdn(kk,ii,jj)
               zghno3(1,1) = pa_gaerop(kk,ii,jj,2)*pdn(kk,ii,jj)
               zgnh3(1,1) = pa_gaerop(kk,ii,jj,3)*pdn(kk,ii,jj)
               zgocnv(1,1) = pa_gaerop(kk,ii,jj,4)*pdn(kk,ii,jj)
               zgocsv(1,1) = pa_gaerop(kk,ii,jj,5)*pdn(kk,ii,jj)

               ! ***************************************!
               !                Run SALSA               !
               ! ***************************************!
               CALL salsa(kproma, kbdim,  klev,   krow,     &
                          in_p,   in_rv,  in_rs,  in_rsi,   &
                          in_t,   tstep,                    &
                          zgso4,  zgocnv, zgocsv, zghno3,   &
                          zgnh3,  allSALSA,                 &
                          actd,   in_w, level               )



               ! Calculate tendencies (convert back to #/kg or kg/kg)
               pa_naerot(kk,ii,jj,1:nbins) = pa_naerot(kk,ii,jj,1:nbins) + &
                    ( aero(1,1,1:nbins)%numc - aero_old(1,1,1:nbins)%numc )/pdn(kk,ii,jj)/tstep
               pa_ncloudt(kk,ii,jj,1:ncld) = pa_ncloudt(kk,ii,jj,1:ncld) + &
                    ( cloud(1,1,1:ncld)%numc - cloud_old(1,1,1:ncld)%numc )/pdn(kk,ii,jj)/tstep
               pa_nprecpt(kk,ii,jj,1:nprc) = pa_nprecpt(kk,ii,jj,1:nprc) + &
                    ( precp(1,1,1:nprc)%numc - precp_old(1,1,1:nprc)%numc )/pdn(kk,ii,jj)/tstep
               pa_nicet(kk,ii,jj,1:nice) = pa_nicet(kk,ii,jj,1:nice) + &
                    ( ice(1,1,1:nice)%numc - ice_old(1,1,1:nice)%numc )/pdn(kk,ii,jj)/tstep
               pa_nsnowt(kk,ii,jj,1:nsnw) = pa_nsnowt(kk,ii,jj,1:nsnw) + &
                    ( snow(1,1,1:nsnw)%numc - snow_old(1,1,1:nsnw)%numc )/pdn(kk,ii,jj)/tstep
               ! Activated droplets
               pa_nactd(kk,ii,jj,1:ncld) = actd(1,1,1:ncld)%numc/pdn(kk,ii,jj)


               ! Get mass tendencies
               DO nc = 1,spec%getNSpec()
                  
                  str = getMassIndex(nbins,1,nc)
                  end = getMassIndex(nbins,nbins,nc)
                  pa_maerot(kk,ii,jj,str:end) = pa_maerot(kk,ii,jj,str:end) + &
                       ( aero(1,1,1:nbins)%volc(nc) - aero_old(1,1,1:nbins)%volc(nc) )*spec%rholiq(nc)/pdn(kk,ii,jj)/tstep

                  str = getMassIndex(ncld,1,nc)
                  end = getMassIndex(ncld,ncld,nc)
                  pa_mcloudt(kk,ii,jj,str:end) = pa_mcloudt(kk,ii,jj,str:end) + &
                       ( cloud(1,1,1:ncld)%volc(nc) - cloud_old(1,1,1:ncld)%volc(nc) )*spec%rholiq(nc)/pdn(kk,ii,jj)/tstep

                  ! Activated droplets (in case of cloud base activation)
                  pa_vactd(kk,ii,jj,str:end) = actd(1,1,1:ncld)%volc(nc)*spec%rholiq(nc)/pdn(kk,ii,jj)
                  
                  str = getMassIndex(nprc,1,nc)
                  end = getMassIndex(nprc,nprc,nc)
                  pa_mprecpt(kk,ii,jj,str:end) = pa_mprecpt(kk,ii,jj,str:end) + &
                       ( precp(1,1,1:nprc)%volc(nc) - precp_old(1,1,1:nprc)%volc(nc) )*spec%rholiq(nc)/pdn(kk,ii,jj)/tstep

                  str = getMassIndex(nice,1,nc)
                  end = getMassIndex(nice,nice,nc)
                  pa_micet(kk,ii,jj,str:end) = pa_micet(kk,ii,jj,str:end) + &
                       ( ice(1,1,1:nice)%volc(nc) - ice_old(1,1,1:nice)%volc(nc) )*spec%rhoice(nc)/pdn(kk,ii,jj)/tstep

                  str = getMassIndex(nsnw,1,nc)
                  end = getMassIndex(nsnw,nsnw,nc)
                  pa_msnowt(kk,ii,jj,str:end) = pa_msnowt(kk,ii,jj,str:end) + &
                       ( snow(1,1,1:nsnw)%volc(nc) - snow_old(1,1,1:nsnw)%volc(nc) )*spec%rhosnow(nc)/pdn(kk,ii,jj)/tstep

               END DO

               IF (lscndgas) THEN
                  pa_gaerot(kk,ii,jj,1) = pa_gaerot(kk,ii,jj,1) + &
                                          ( zgso4(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,1) )/tstep

                  pa_gaerot(kk,ii,jj,2) = pa_gaerot(kk,ii,jj,2) + &
                                          ( zghno3(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,2) )/tstep

                  pa_gaerot(kk,ii,jj,3) = pa_gaerot(kk,ii,jj,3) + &
                                          ( zgnh3(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,3) )/tstep

                  pa_gaerot(kk,ii,jj,4) = pa_gaerot(kk,ii,jj,4) + &
                                          ( zgocnv(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,4) )/tstep

                  pa_gaerot(kk,ii,jj,5) = pa_gaerot(kk,ii,jj,5) + &
                                          ( zgocsv(1,1)/pdn(kk,ii,jj) - pa_gaerop(kk,ii,jj,5) )/tstep
               END IF

               ! Tendency of water vapour mixing ratio 
               rt(kk,ii,jj) = rt(kk,ii,jj) + &
                  ( in_rv(1,1) - rv_old(1,1) )/tstep

            END DO ! kk
         END DO ! ii
      END DO ! jj

   END SUBROUTINE run_SALSA

   !
   !---------------------------------------------------------------
   ! SET_SALSA_RUNTIME
   ! Set LOGICAL switches according to the host model state and
   ! user-specified NAMELIST options.
   !
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE set_SALSA_runtime(prunmode)
      USE mo_submctl, ONLY : nlcoag,                 &
                             nlcgaa,nlcgcc,nlcgpp,   &
                             nlcgca,nlcgpa,nlcgpc,   &
                             nlcnd,                  &
                             nlcndgas,               &
                             nlcndh2oae, nlcndh2ocl, &
                             nlcndh2oic,             &
                             nlauto,nlautosnow,      &
                             nlactiv,                &
                             nlactbase,nlactintst,   &

                             lscoag,                 &
                             lscgaa,lscgcc,lscgpp,   &
                             lscgca,lscgpa,lscgpc,   &
                             lscnd,                  &
                             lscndgas,               &
                             lscndh2oae, lscndh2ocl, &
                             lscndh2oic,             &
                             lsauto,lsautosnow,      &
                             lsactiv,                &
                             lsactbase,lsactintst,   &

                             nlcgia,nlcgic,nlcgii,   &
                             nlcgip,nlcgsa,nlcgsc,   &
                             nlcgsi,nlcgsp,nlcgss,   &
                             nlcnd,                  &
                             nlicenucl,               &
                             nlicmelt,               &
                             nlfixinc,               &

                             lscgia,lscgic,lscgii,   &
                             lscgip,lscgsa,lscgsc,   &
                             lscgsi,lscgsp,lscgss,   &
                             lsicenucl,                &
                             lsicmelt,               &
                             lsfixinc

      IMPLICIT NONE

      INTEGER, INTENT(in) :: prunmode

      ! Apply runtime settings

      lscoag      = nlcoag
      lscgaa      = nlcgaa
      lscgcc      = nlcgcc
      lscgpp      = nlcgpp
      lscgca      = nlcgca
      lscgpa      = nlcgpa
      lscgpc      = nlcgpc
      lscgia      = nlcgia
      lscgic      = nlcgic
      lscgii      = nlcgii
      lscgip      = nlcgip
      lscgsa      = nlcgsa
      lscgsc      = nlcgsc
      lscgsi      = nlcgsi
      lscgsp      = nlcgsp
      lscgss      = nlcgss

      lscnd       = nlcnd
      lscndgas    = nlcndgas
      lscndh2oae  = nlcndh2oae
      lscndh2ocl  = nlcndh2ocl
      lscndh2oic  = nlcndh2oic

      lsauto      = nlauto
      lsautosnow  = nlautosnow

      lsactiv     = nlactiv
      lsactbase   = nlactbase
      lsactintst  = nlactintst

      lsicenucl  = ( nlicenucl .AND. ( .NOT. nlfixinc ) )
      lsfixinc    = nlfixinc
      lsicmelt    = nlicmelt

      ! Adjustments for initialization and spinup

      SELECT CASE(prunmode)

         CASE(1) ! Initialization

            lscoag      = .FALSE.
            lsauto      = .FALSE.
            lsautosnow  = .FALSE.
            lsactbase   = .FALSE.
            lsactintst  = nlactintst
            lsicenucl  = .FALSE.
            lsfixinc    = .FALSE.
            lsicmelt    = .FALSE.

         CASE(2)  ! Spinup period

            lscoag      = .FALSE.
            lsauto      = .FALSE.
            lsautosnow  = .FALSE.

      END SELECT

   END SUBROUTINE set_SALSA_runtime


END MODULE mo_salsa_driver
