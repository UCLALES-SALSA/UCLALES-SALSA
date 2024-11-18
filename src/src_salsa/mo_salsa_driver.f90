MODULE mo_salsa_driver

USE mo_submctl, ONLY : t_section
IMPLICIT NONE

!---------------------------------------------------------------
!
! MO_SALSA_DRIVER:
! Contains the primary SALSA input/output variables as well as
! subroutines used to call the main SALSA routine.
!
! Juha Tonttila, FMI, 2014
!
!---------------------------------------------------------------


  ! JT: Variables from SALSA
  ! --------------------------------------------
  ! grid points for SALSA
  INTEGER, PARAMETER :: kbdim = 1
  INTEGER, PARAMETER :: klev = 1

  ! -- Local hydrometeor properties (set up in aero initialize)
  TYPE(t_section), ALLOCATABLE :: cloud(:,:,:) ! cloud properties
  TYPE(t_section), ALLOCATABLE :: aero(:,:,:)  ! Aerosol properties
  TYPE(t_section), ALLOCATABLE :: precp(:,:,:) ! Precipitation properties
  TYPE(t_section), ALLOCATABLE :: ice(:,:,:) ! ice properties
  TYPE(t_section), ALLOCATABLE :: snow(:,:,:) ! snow aka. ice precip. properties

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
  SUBROUTINE run_SALSA(pnx, pny, pnz, n4, nbins, ncld, nprc, nice, nsnw, &
                       press, tk, rv, rt, rs, rsi, pdn, edr, &
                       pa_naerop,  pa_naerot,  pa_maerop,  pa_maerot,   &
                       pa_ncloudp, pa_ncloudt, pa_mcloudp, pa_mcloudt,  &
                       pa_nprecpp, pa_nprecpt, pa_mprecpp, pa_mprecpt,  &
                       pa_nicep,   pa_nicet,   pa_micep,   pa_micet,    &
                       pa_nsnowp,  pa_nsnowt,  pa_msnowp,  pa_msnowt,   &
                       pa_gasp,  pa_gast, prunmode, tstep, time, level, &
                       sflg, nstat, slist, sdata)

    USE mo_submctl, ONLY : dens, rhlim, lscndgas, ngases, mws_gas, nspec, &
                               ngases_diag, zgas_diag, set_vbs_diag, eddy_dis_rt
    USE mo_salsa, ONLY : salsa
    IMPLICIT NONE

    INTEGER, INTENT(in) :: pnx,pny,pnz,n4, &                ! Dimensions: x,y,z,number of chemical species
                               nbins,ncld,nprc,nice,nsnw    ! ... and number of *prognostic* bins
    REAL, INTENT(in)    :: tstep, time                      ! Model timestep length and time

    REAL, INTENT(in)    :: press(pnz,pnx,pny), &            ! Pressure (Pa)
                               tk(pnz,pnx,pny),    &            ! Temperature (K)
                               rv(pnz,pnx,pny),    &            ! Water vapor mixing ratio
                               rs(pnz,pnx,pny),    &            ! Water vapour saturation mixing ratio
                               rsi(pnz,pnx,pny)                 ! Water vapour saturation mixing ratio over ice

    REAL, INTENT(in)    :: pdn(pnz,pnx,pny)             ! Air density (for normalizing concentrations)
    REAL, INTENT(inout) :: edr(pnz,pnx,pny)             ! Eddy dissipation rate for coagulation

    REAL, INTENT(in)    :: pa_naerop(pnz,pnx,pny,nbins),        & ! Aerosol number concentration (# kg-1)
                               pa_maerop(pnz,pnx,pny,n4*nbins),     & ! Aerosol mass concentration (kg kg-1)
                               pa_ncloudp(pnz,pnx,pny,ncld),        & ! Cloud droplet number concentration (# kg-1)
                               pa_mcloudp(pnz,pnx,pny,n4*ncld),     & ! Cloud droplet mass concentration (kg kg-1)
                               pa_nprecpp(pnz,pnx,pny,nprc),        & ! Rain drop number concentration (# kg-1)
                               pa_mprecpp(pnz,pnx,pny,n4*nprc),     & ! Rain drop mass concentration (kg kg-1)
                               pa_nicep(pnz,pnx,pny,nice),          & ! Ice number concentration (# kg-1)
                               pa_micep(pnz,pnx,pny,n4*nice),       & ! Ice mass concentration (kg kg-1)
                               pa_nsnowp(pnz,pnx,pny,nsnw),         & ! Snow number concentration (# kg-1)
                               pa_msnowp(pnz,pnx,pny,n4*nsnw)         ! Snow mass concentration (kg kg-1)

    REAL, INTENT(in)    :: pa_gasp(pnz,pnx,pny,ngases)   ! Gaseous tracer concentration [kg kg-1)

    LOGICAL, INTENT(in) :: prunmode                      ! Spinup (T) or regular (F) call
    INTEGER, INTENT(in) :: level                         ! thermodynamical level

    LOGICAL, INTENT(IN) :: sflg                          ! statistics sampling flag
    INTEGER, INTENT(IN) :: nstat                         ! the number of requested outputs
    CHARACTER(LEN=7), DIMENSION(:), INTENT(IN) :: slist  ! names of the output variables
    REAL, OPTIONAL, INTENT(OUT) :: sdata(pnz,pnx,pny,nstat) ! output data (needed when saving data)

    REAL, INTENT(inout)   :: pa_naerot(pnz,pnx,pny,nbins),      & ! Aerosol number tendency
                                 pa_maerot(pnz,pnx,pny,n4*nbins),   & ! Aerosol mass tendency
                                 pa_ncloudt(pnz,pnx,pny,ncld),      & ! Cloud droplet number tendency
                                 pa_mcloudt(pnz,pnx,pny,n4*ncld),   & ! Cloud droplet mass tendency
                                 pa_nprecpt(pnz,pnx,pny,nprc),      & ! Rain drop number tendency
                                 pa_mprecpt(pnz,pnx,pny,n4*nprc),   & ! Rain drop mass tendency
                                 pa_nicet(pnz,pnx,pny,nice),        & ! Ice number tendency
                                 pa_micet(pnz,pnx,pny,n4*nice),     & ! Ice mass tendency
                                 pa_nsnowt(pnz,pnx,pny,nsnw),       & ! Snow number tendency
                                 pa_msnowt(pnz,pnx,pny,n4*nsnw)       ! Snow mass tendency

    REAL, INTENT(inout)   :: pa_gast(pnz,pnx,pny,ngases)      ! Gaseous tracer tendency
    REAL, INTENT(inout)   :: rt(pnz,pnx,pny)                  ! Water vapour tendency

    ! -- Local gas concentrations [mol m-3]
    REAL :: zgas(kbdim,klev,ngases+ngases_diag)

    INTEGER :: jj,ii,kk,str,end,nc
    REAL, DIMENSION(kbdim,klev) :: in_p, in_t, in_rv, in_rs, in_rsi, in_edr
    REAL :: rv_old(kbdim,klev), rho
    REAL :: out_sdata(kbdim,klev,nstat)

    IF (sflg .AND. nstat>0) sdata(:,:,:,:) = 0.

    ! Set the SALSA runtime config
    CALL set_salsa_runtime(prunmode,time)

    ! VBS settings for diagnostic oxidants and aqSOA photodissociation
    CALL set_vbs_diag(time)

    ! If eddy_dis_rt<0, then use that from LES
    IF (eddy_dis_rt >= 0.) edr(:,:,:) = eddy_dis_rt

    ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by substance density)
    DO jj = 3,pny-2
       DO ii = 3,pnx-2
          DO kk = pnz-1,2,-1

             ! Set inputs
             in_p(1,1) = press(kk,ii,jj)
             in_t(1,1) = tk(kk,ii,jj)
             in_rs(1,1) = rs(kk,ii,jj)
             in_rsi(1,1) = rsi(kk,ii,jj)
             in_edr(1,1) = edr(kk,ii,jj)

             ! For initialization and spinup, limit the RH with the parameter rhlim (assign in namelist.salsa)
             IF (prunmode) THEN
                IF (rhlim>2.0) THEN
                    ! Time-dependent water vapor mixing ratio limit: initially limited to saturation and exponential
                    ! relaxation towards current mixing ratio. Here rhlim is the relaxation time [s].
                    in_rv(1,1) = MIN( rv(kk,ii,jj), rv(kk,ii,jj)+(rs(kk,ii,jj)-rv(kk,ii,jj))*exp(-time/rhlim) )
                ELSE
                    ! Constant water vapor mixing ratio limit. Here rhlim is the maximum saturation ratio.
                    in_rv(1,1) = MIN(rv(kk,ii,jj), rs(kk,ii,jj)*rhlim)
                ENDIF
             ELSE
                in_rv(1,1) = rv(kk,ii,jj)
             END IF
             rv_old(1,1) = in_rv(1,1)

             ! Set volume concentrations
             DO nc=1,nspec+1
                rho=pdn(kk,ii,jj)/dens(nc) ! rho_air/rho
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(nc) = pa_maerop(kk,ii,jj,str:end)*rho

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(nc) = pa_mcloudp(kk,ii,jj,str:end)*rho

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(nc) = pa_mprecpp(kk,ii,jj,str:end)*rho

                str = (nc-1)*nice+1
                end = nc*nice
                ice(1,1,1:nice)%volc(nc) = pa_micep(kk,ii,jj,str:end)*rho

                str = (nc-1)*nsnw+1
                end = nc*nsnw
                snow(1,1,1:nsnw)%volc(nc) = pa_msnowp(kk,ii,jj,str:end)*rho
             END DO

             ! Number concentrations and particle sizes
             aero(1,1,1:nbins)%numc = pa_naerop(kk,ii,jj,1:nbins)*pdn(kk,ii,jj)
             cloud(1,1,1:ncld)%numc = pa_ncloudp(kk,ii,jj,1:ncld)*pdn(kk,ii,jj)
             precp(1,1,1:nprc)%numc = pa_nprecpp(kk,ii,jj,1:nprc)*pdn(kk,ii,jj)
             ice(1,1,1:nice)%numc = pa_nicep(kk,ii,jj,1:nice)*pdn(kk,ii,jj)
             snow(1,1,1:nsnw)%numc = pa_nsnowp(kk,ii,jj,1:nsnw)*pdn(kk,ii,jj)


             ! Condensable gases (sulfate and organics)
             IF (lscndgas .AND. ngases+ngases_diag>0) THEN
                ! Convert from kg/kg to mol/m^3
                zgas(1,1,1:ngases) = pa_gasp(kk,ii,jj,1:ngases)*pdn(kk,ii,jj)/mws_gas(1:ngases)
                IF (ngases_diag>0) THEN
                    ! Diagnostic gases - units are mol/kg, but convert to mol/m^3
                    zgas(1,1,ngases+1:ngases+ngases_diag) = zgas_diag(1:ngases_diag)*pdn(kk,ii,jj)
                ENDIF
             ENDIF


             ! ***************************************!
             !                Run SALSA               !
             ! ***************************************!
             CALL salsa(kbdim,  klev,                          &
                        in_p,   in_rv,  in_rs,  in_rsi,        &
                        in_t,   in_edr, tstep,                 &
                        zgas,   ngases+ngases_diag,            &
                        aero,   cloud,  precp,                 &
                        ice,    snow,                          &
                        level, (sflg .AND. nstat>0), nstat, out_sdata, slist)


             ! Output statistics (mixing ratios from m^3/m^3 to kg/kg and concentrations from 1/m^3 to 1/kg;
             ! also converted to rates by dividing by the time step)
             IF (sflg .AND. nstat>0) sdata(kk,ii,jj,:) = out_sdata(1,1,:)/pdn(kk,ii,jj)/tstep

             ! Calculate tendencies (convert back to #/kg or kg/kg)
             pa_naerot(kk,ii,jj,1:nbins) = pa_naerot(kk,ii,jj,1:nbins) + &
                  ( aero(1,1,1:nbins)%numc/pdn(kk,ii,jj) - pa_naerop(kk,ii,jj,1:nbins) )/tstep
             pa_ncloudt(kk,ii,jj,1:ncld) = pa_ncloudt(kk,ii,jj,1:ncld) + &
                  ( cloud(1,1,1:ncld)%numc/pdn(kk,ii,jj) - pa_ncloudp(kk,ii,jj,1:ncld) )/tstep
             pa_nprecpt(kk,ii,jj,1:nprc) = pa_nprecpt(kk,ii,jj,1:nprc) + &
                  ( precp(1,1,1:nprc)%numc/pdn(kk,ii,jj) - pa_nprecpp(kk,ii,jj,1:nprc) )/tstep
             pa_nicet(kk,ii,jj,1:nice) = pa_nicet(kk,ii,jj,1:nice) + &
                  ( ice(1,1,1:nice)%numc/pdn(kk,ii,jj) - pa_nicep(kk,ii,jj,1:nice) )/tstep
             pa_nsnowt(kk,ii,jj,1:nsnw) = pa_nsnowt(kk,ii,jj,1:nsnw) + &
                  ( snow(1,1,1:nsnw)%numc/pdn(kk,ii,jj) - pa_nsnowp(kk,ii,jj,1:nsnw) )/tstep

             DO nc=1,nspec+1
                rho=dens(nc)/pdn(kk,ii,jj) ! rho/rho_air
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_maerot(kk,ii,jj,str:end) = pa_maerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(nc)*rho - pa_maerop(kk,ii,jj,str:end) )/tstep

                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_mcloudt(kk,ii,jj,str:end) = pa_mcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(nc)*rho - pa_mcloudp(kk,ii,jj,str:end) )/tstep

                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_mprecpt(kk,ii,jj,str:end) = pa_mprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(nc)*rho - pa_mprecpp(kk,ii,jj,str:end) )/tstep

                str = (nc-1)*nice+1
                end = nc*nice
                pa_micet(kk,ii,jj,str:end) = pa_micet(kk,ii,jj,str:end) + &
                     ( ice(1,1,1:nice)%volc(nc)*rho - pa_micep(kk,ii,jj,str:end) )/tstep

                str = (nc-1)*nsnw+1
                end = nc*nsnw
                pa_msnowt(kk,ii,jj,str:end) = pa_msnowt(kk,ii,jj,str:end) + &
                     ( snow(1,1,1:nsnw)%volc(nc)*rho - pa_msnowp(kk,ii,jj,str:end) )/tstep
             END DO

             IF (lscndgas .AND. ngases>0) THEN
                ! Prognostic gases only
                pa_gast(kk,ii,jj,1:ngases) = pa_gast(kk,ii,jj,1:ngases) + &
                        ( zgas(1,1,1:ngases)/pdn(kk,ii,jj)*mws_gas(1:ngases) - pa_gasp(kk,ii,jj,1:ngases) )/tstep
             ENDIF


             ! Tendency of water vapour mixing ratio is obtained from the change in RH during SALSA run.
             ! Assumes no temperature change during SALSA run.
             rt(kk,ii,jj) = rt(kk,ii,jj) + &
                  ( in_rv(1,1) - rv_old(1,1) )/tstep

          END DO ! kk
       END DO ! ii
    END DO ! jj

  END SUBROUTINE run_SALSA

  !
  !---------------------------------------------------------------
  ! SET_SALSA_RUNTIME
  ! Set logical switches according to the host model state and
  ! user-specified NAMELIST options.
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE set_SALSA_runtime(prunmode,time)
    USE mo_submctl, ONLY : nlcoag,                 &
                               nlcgaa,nlcgcc,nlcgpp,   &
                               nlcgca,nlcgpa,nlcgpc,   &
                               nlcgrain,               &
                               nlcnd,                  &
                               nlcndgas,               &
                               soa_tstart,             &
                               nlauto,nlautosnow,      &
                               nlactiv,                &

                               lscoag,                 &
                               lscgaa,lscgcc,lscgpp,   &
                               lscgca,lscgpa,lscgpc,   &
                               lscgrain,               &
                               lscnd,                  &
                               lscndgas,               &
                               lsauto,lsautosnow,      &
                               lsactiv,                &

                               nlcgia,nlcgic,nlcgii,   &
                               nlcgip,nlcgsa,nlcgsc,   &
                               nlcgsi,nlcgsp,nlcgss,   &
                               nlcnd,                  &
                               nlicenucl,              &
                               nlicmelt,               &
                               icenucl_tstart,         &

                               lscgia,lscgic,lscgii,   &
                               lscgip,lscgsa,lscgsc,   &
                               lscgsi,lscgsp,lscgss,   &
                               lsicenucl,              &
                               lsicmelt

    IMPLICIT NONE

    LOGICAL, INTENT(in) :: prunmode
    REAL, INTENT(in) :: time

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

    lscgrain    = nlcgrain
    lsauto      = nlauto
    lsautosnow  = nlautosnow

    lsactiv     = nlactiv

    lsicenucl  = nlicenucl
    lsicmelt    = nlicmelt


    ! Adjustments for spinup
    IF (prunmode) THEN
          lscoag      = .FALSE.
          lsauto      = .FALSE.
          lsautosnow  = .FALSE.
    ENDIF

    ! Ice formation has an additional spinup time (no ice formation => no autoconversion or melting)
    IF (time<icenucl_tstart) THEN
          lsicenucl  = .FALSE.
          lsautosnow = .FALSE.
          lsicmelt   = .FALSE.
    ENDIF

    ! Spinup for SOA
    IF (time<soa_tstart) lscndgas = .FALSE.

  END SUBROUTINE set_SALSA_runtime


END MODULE mo_salsa_driver
