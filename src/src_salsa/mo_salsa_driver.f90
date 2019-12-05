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

  REAL, PARAMETER :: init_rh(kbdim,klev) = 0.3

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
  SUBROUTINE run_SALSA(pnx, pny, pnz, n4, press, tk, rv, rt, rs, rsi, pdn,   &
                       pa_naerop,  pa_naerot,  pa_maerop,  pa_maerot,   &
                       pa_ncloudp, pa_ncloudt, pa_mcloudp, pa_mcloudt,  &
                       pa_nprecpp, pa_nprecpt, pa_mprecpp, pa_mprecpt,  &
                       pa_nicep,   pa_nicet,   pa_micep,   pa_micet,    &
                       pa_nsnowp,  pa_nsnowt,  pa_msnowp,  pa_msnowt,   &
                       pa_gasp,  pa_gast, prunmode, tstep, time, level, &
                       coag_ra, coag_na, coag_rc, coag_nc, coag_rr,     &
                       coag_nr, coag_ri, coag_ni, coag_rs, coag_ns,     &
                       cond_ra, cond_rc, cond_rr, cond_ri, cond_rs,     &
                       auto_rr, auto_nr, auto_rs, auto_ns,              &
                       cact_rc, cact_nc, nucl_ri, nucl_ni,              &
                       melt_ri, melt_ni, melt_rs, melt_ns)

    USE mo_submctl, ONLY : nbins,ncld,nprc,nice,nsnw,pi6,          &
                               rhoic,rhosn, rhowa, dens, &
                               rhlim, lscndgas, ngases, mws_gas, nlim, prlim, nspec, maxnspec
    USE mo_salsa, ONLY : salsa
    USE mo_salsa_properties, ONLY : equilibration
    IMPLICIT NONE

    INTEGER, INTENT(in) :: pnx,pny,pnz,n4                       ! Dimensions: x,y,z,number of chemical species  
    REAL, INTENT(in)    :: tstep, time                      ! Model timestep length and time

    REAL, INTENT(in)    :: press(pnz,pnx,pny), &            ! Pressure (Pa)
                               tk(pnz,pnx,pny),    &            ! Temperature (K)
                               rv(pnz,pnx,pny),    &            ! Water vapor mixing ratio
                               rs(pnz,pnx,pny),    &            ! Water vapour saturation mixing ratio
                               rsi(pnz,pnx,pny)                 ! Water vapour saturation mixing ratio over ice

    REAL, INTENT(in)    :: pdn(pnz,pnx,pny)             ! Air density (for normalizing concentrations)

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

    INTEGER, INTENT(in) :: prunmode                      ! 1: Initialization call
                                                         ! 2: Spinup period call
                                                         ! 3: Regular runtime call
    INTEGER, INTENT(in) :: level                         ! thermodynamical level

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

    REAL, DIMENSION(pnz,pnx,pny), INTENT(OUT) :: & ! Statistics
                       coag_ra, coag_na, coag_rc, coag_nc, coag_rr, &
                       coag_nr, coag_ri, coag_ni, coag_rs, coag_ns, &
                       cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, &
                       auto_rr, auto_nr, auto_rs, auto_ns, &
                       cact_rc, cact_nc, nucl_ri, nucl_ni, &
                       melt_ri, melt_ni, melt_rs, melt_ns

    ! -- Local gas concentrations [mol m-3]
    REAL :: zgas(kbdim,klev,ngases)

    ! Helper arrays for calculating the rates of change
    TYPE(t_section) :: aero_old(1,1,nbins), cloud_old(1,1,ncld), precp_old(1,1,nprc), &
       ice_old(1,1,nice), snow_old(1,1,nsnw)

    INTEGER :: jj,ii,kk,ss,str,end,nc
    REAL, DIMENSION(kbdim,klev) :: in_p, in_t, in_rv, in_rs, in_rsi, &
                out_coag_va, out_coag_na, out_coag_vc, out_coag_nc, out_coag_vr, &
                out_coag_nr, out_coag_vi, out_coag_ni, out_coag_vs, out_coag_ns, &
                out_cond_va, out_cond_vc, out_cond_vr, out_cond_vi, out_cond_vs, &
                out_auto_vr, out_auto_nr, out_auto_vs, out_auto_ns, &
                out_cact_vc, out_cact_nc, out_nucl_vi, out_nucl_ni, &
                out_melt_vi, out_melt_ni, out_melt_vs, out_melt_ns
    REAL :: rv_old(kbdim,klev), rho

    ! Number is always set, but mass can be uninitialized
    DO ss = 1,maxnspec
       aero(:,:,:)%volc(ss) = 0.
       cloud(:,:,:)%volc(ss) = 0.
       precp(:,:,:)%volc(ss) = 0.
       ice(:,:,:)%volc(ss) = 0.
       snow(:,:,:)%volc(ss) = 0.
       aero_old(:,:,:)%volc(ss) = 0.
       cloud_old(:,:,:)%volc(ss) = 0.
       precp_old(:,:,:)%volc(ss) = 0.
       ice_old(:,:,:)%volc(ss) = 0.
       snow_old(:,:,:)%volc(ss) = 0.
    END DO

    ! Set the SALSA runtime config
    CALL set_salsa_runtime(prunmode,time)

    ! Convert input concentrations for SALSA into #/m3 or m3/m3 instead of kg/kg (multiplied by pdn/divided by substance density)
    DO jj = 3,pny-2
       DO ii = 3,pnx-2
          DO kk = pnz-1,2,-1

             ! Set inputs
             in_p(1,1) = press(kk,ii,jj)
             in_t(1,1) = tk(kk,ii,jj)
             in_rs(1,1) = rs(kk,ii,jj)
             in_rsi(1,1) = rsi(kk,ii,jj)

             ! For initialization and spinup, limit the RH with the parameter rhlim (assign in namelist.salsa)
             IF (prunmode < 3) THEN
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
                rho=dens(nc) ! Density - not valid for water in ice and snow
                str = (nc-1)*nbins+1
                end = nc*nbins
                aero(1,1,1:nbins)%volc(nc) = pa_maerop(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rho
                aero_old(1,1,1:nbins)%volc(nc) = aero(1,1,1:nbins)%volc(nc)

                str = (nc-1)*ncld+1
                end = nc*ncld
                cloud(1,1,1:ncld)%volc(nc) = pa_mcloudp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rho
                cloud_old(1,1,1:ncld)%volc(nc) = cloud(1,1,1:ncld)%volc(nc)

                str = (nc-1)*nprc+1
                end = nc*nprc
                precp(1,1,1:nprc)%volc(nc) = pa_mprecpp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rho
                precp_old(1,1,1:nprc)%volc(nc) = precp(1,1,1:nprc)%volc(nc)

                IF (nc==1) rho=rhoic ! Ice water
                str = (nc-1)*nice+1
                end = nc*nice
                ice(1,1,1:nice)%volc(nc) = pa_micep(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rho
                ice_old(1,1,1:nice)%volc(nc) = ice(1,1,1:nice)%volc(nc)

                IF (nc==1) rho=rhosn ! Snow water
                str = (nc-1)*nsnw+1
                end = nc*nsnw
                snow(1,1,1:nsnw)%volc(nc) = pa_msnowp(kk,ii,jj,str:end)*pdn(kk,ii,jj)/rho
                snow_old(1,1,1:nsnw)%volc(nc) = snow(1,1,1:nsnw)%volc(nc)

             END DO

             ! Number concentrations and particle sizes
             aero(1,1,1:nbins)%numc = pa_naerop(kk,ii,jj,1:nbins)*pdn(kk,ii,jj)
             aero_old(1,1,1:nbins)%numc = aero(1,1,1:nbins)%numc
             DO ss=1,nbins
                IF (aero(1,1,ss)%numc>nlim) THEN
                    aero(1,1,ss)%dwet = ( SUM(aero(1,1,ss)%volc(:))/aero(1,1,ss)%numc/pi6 )**(1./3.)
                ELSE
                    aero(1,1,ss)%dwet = aero(1,1,ss)%dmid
                ENDIF
             ENDDO

             cloud(1,1,1:ncld)%numc = pa_ncloudp(kk,ii,jj,1:ncld)*pdn(kk,ii,jj)
             cloud_old(1,1,1:ncld)%numc = cloud(1,1,1:ncld)%numc
             DO ss=1,ncld
                IF (cloud(1,1,ss)%numc>nlim) THEN
                    cloud(1,1,ss)%dwet = ( SUM(cloud(1,1,ss)%volc(:))/cloud(1,1,ss)%numc/pi6 )**(1./3.)
                ELSE
                    cloud(1,1,ss)%dwet = cloud(1,1,ss)%dmid
                ENDIF
             ENDDO

             precp(1,1,1:nprc)%numc = pa_nprecpp(kk,ii,jj,1:nprc)*pdn(kk,ii,jj)
             precp_old(1,1,1:nprc)%numc = precp(1,1,1:nprc)%numc
             DO ss=1,nprc
                IF (precp(1,1,ss)%numc>prlim) THEN
                    precp(1,1,ss)%dwet = ( SUM(precp(1,1,ss)%volc(:))/precp(1,1,ss)%numc/pi6 )**(1./3.)
                ELSE
                    precp(1,1,ss)%dwet = precp(1,1,ss)%dmid
                ENDIF
             ENDDO

             ice(1,1,1:nice)%numc = pa_nicep(kk,ii,jj,1:nice)*pdn(kk,ii,jj)
             ice_old(1,1,1:nice)%numc = ice(1,1,1:nice)%numc
             DO ss=1,nice
                IF (ice(1,1,ss)%numc>prlim) THEN
                    ice(1,1,ss)%dwet = ( SUM(ice(1,1,ss)%volc(:))/ice(1,1,ss)%numc/pi6 )**(1./3.)
                ELSE
                    ice(1,1,ss)%dwet = ice(1,1,ss)%dmid
                ENDIF
             ENDDO

             snow(1,1,1:nsnw)%numc = pa_nsnowp(kk,ii,jj,1:nsnw)*pdn(kk,ii,jj)
             snow_old(1,1,1:nsnw)%numc = snow(1,1,1:nsnw)%numc
             DO ss=1,nsnw
                IF (snow(1,1,ss)%numc>prlim) THEN
                    snow(1,1,ss)%dwet = ( SUM(snow(1,1,ss)%volc(:))/snow(1,1,ss)%numc/pi6 )**(1./3.)
                ELSE
                    snow(1,1,ss)%dwet = snow(1,1,ss)%dmid
                ENDIF
             ENDDO


             ! Condensable gases (sulfate and organics)
             IF (lscndgas .AND. ngases>0) THEN
                ! Convert from kg/kg to mol/m^3
                zgas(1,1,1:ngases) = pa_gasp(kk,ii,jj,1:ngases)*pdn(kk,ii,jj)/mws_gas(1:ngases)
             ENDIF


             ! If this is an initialization call, calculate the equilibrium particle
             If (prunmode == 1) CALL equilibration(kbdim,klev,   &
                                                    init_rh,in_t,aero,.TRUE.)


             ! ***************************************!
             !                Run SALSA               !
             ! ***************************************!
             CALL salsa(kbdim,  klev,                          &
                        in_p,   in_rv,  in_rs,  in_rsi,        &
                        in_t,   tstep,                         &
                        zgas,   ngases,                        &
                        aero,   cloud,  precp,                 &
                        ice,    snow,                          &
                        level,                                 &
                        out_coag_va, out_coag_na, out_coag_vc, out_coag_nc, out_coag_vr, &
                        out_coag_nr, out_coag_vi, out_coag_ni, out_coag_vs, out_coag_ns, &
                        out_cond_va, out_cond_vc, out_cond_vr, out_cond_vi, out_cond_vs, &
                        out_auto_vr, out_auto_nr, out_auto_vs, out_auto_ns, &
                        out_cact_vc, out_cact_nc, out_nucl_vi, out_nucl_ni, &
                        out_melt_vi, out_melt_ni, out_melt_vs, out_melt_ns)


             ! Output statistics (mixing ratios from m^3/m^3 to kg/kg and concentrations from 1/m^3 to 1/kg;
             ! also converted to rates by dividing by the time step)
             coag_ra(kk,ii,jj)=out_coag_va(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             coag_na(kk,ii,jj)=out_coag_na(1,1)/pdn(kk,ii,jj)/tstep
             coag_rc(kk,ii,jj)=out_coag_vc(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             coag_nc(kk,ii,jj)=out_coag_nc(1,1)/pdn(kk,ii,jj)/tstep
             coag_rr(kk,ii,jj)=out_coag_vr(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             coag_nr(kk,ii,jj)=out_coag_nr(1,1)/pdn(kk,ii,jj)/tstep
             coag_ri(kk,ii,jj)=out_coag_vi(1,1)*rhoic/pdn(kk,ii,jj)/tstep
             coag_ni(kk,ii,jj)=out_coag_ni(1,1)/pdn(kk,ii,jj)/tstep
             coag_rs(kk,ii,jj)=out_coag_vs(1,1)*rhosn/pdn(kk,ii,jj)/tstep
             coag_ns(kk,ii,jj)=out_coag_ns(1,1)/pdn(kk,ii,jj)/tstep
             cond_ra(kk,ii,jj)=out_cond_va(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             cond_rc(kk,ii,jj)=out_cond_vc(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             cond_rr(kk,ii,jj)=out_cond_vr(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             cond_ri(kk,ii,jj)=out_cond_vi(1,1)*rhoic/pdn(kk,ii,jj)/tstep
             cond_rs(kk,ii,jj)=out_cond_vs(1,1)*rhosn/pdn(kk,ii,jj)/tstep
             auto_rr(kk,ii,jj)=out_auto_vr(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             auto_nr(kk,ii,jj)=out_auto_nr(1,1)/pdn(kk,ii,jj)/tstep
             auto_rs(kk,ii,jj)=out_auto_vs(1,1)*rhosn/pdn(kk,ii,jj)/tstep
             auto_ns(kk,ii,jj)=out_auto_ns(1,1)/pdn(kk,ii,jj)/tstep
             cact_rc(kk,ii,jj)=out_cact_vc(1,1)*rhowa/pdn(kk,ii,jj)/tstep
             cact_nc(kk,ii,jj)=out_cact_nc(1,1)/pdn(kk,ii,jj)/tstep
             nucl_ri(kk,ii,jj)=out_nucl_vi(1,1)*rhoic/pdn(kk,ii,jj)/tstep
             nucl_ni(kk,ii,jj)=out_nucl_ni(1,1)/pdn(kk,ii,jj)/tstep
             melt_ri(kk,ii,jj)=out_melt_vi(1,1)*rhoic/pdn(kk,ii,jj)/tstep
             melt_ni(kk,ii,jj)=out_melt_ni(1,1)/pdn(kk,ii,jj)/tstep
             melt_rs(kk,ii,jj)=out_melt_vs(1,1)*rhosn/pdn(kk,ii,jj)/tstep
             melt_ns(kk,ii,jj)=out_melt_ns(1,1)/pdn(kk,ii,jj)/tstep

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

             DO nc=1,nspec+1
                rho=dens(nc)
                str = (nc-1)*nbins+1
                end = nc*nbins
                pa_maerot(kk,ii,jj,str:end) = pa_maerot(kk,ii,jj,str:end) + &
                     ( aero(1,1,1:nbins)%volc(nc) - aero_old(1,1,1:nbins)%volc(nc) )*rho/pdn(kk,ii,jj)/tstep

                str = (nc-1)*ncld+1
                end = nc*ncld
                pa_mcloudt(kk,ii,jj,str:end) = pa_mcloudt(kk,ii,jj,str:end) + &
                     ( cloud(1,1,1:ncld)%volc(nc) - cloud_old(1,1,1:ncld)%volc(nc) )*rho/pdn(kk,ii,jj)/tstep

                str = (nc-1)*nprc+1
                end = nc*nprc
                pa_mprecpt(kk,ii,jj,str:end) = pa_mprecpt(kk,ii,jj,str:end) + &
                     ( precp(1,1,1:nprc)%volc(nc) - precp_old(1,1,1:nprc)%volc(nc) )*rho/pdn(kk,ii,jj)/tstep

                IF (nc==1) rho=rhoic
                str = (nc-1)*nice+1
                end = nc*nice
                pa_micet(kk,ii,jj,str:end) = pa_micet(kk,ii,jj,str:end) + &
                     ( ice(1,1,1:nice)%volc(nc) - ice_old(1,1,1:nice)%volc(nc) )*rho/pdn(kk,ii,jj)/tstep

                IF (nc==1) rho=rhosn
                str = (nc-1)*nsnw+1
                end = nc*nsnw
                pa_msnowt(kk,ii,jj,str:end) = pa_msnowt(kk,ii,jj,str:end) + &
                     ( snow(1,1,1:nsnw)%volc(nc) - snow_old(1,1,1:nsnw)%volc(nc) )*rho/pdn(kk,ii,jj)/tstep
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
                               nlicenucl,              &
                               nlicmelt,               &
                               icenucl_tstart,         &

                               lscgia,lscgic,lscgii,   &
                               lscgip,lscgsa,lscgsc,   &
                               lscgsi,lscgsp,lscgss,   &
                               lsicenucl,              &
                               lsicmelt

    IMPLICIT NONE

    INTEGER, INTENT(in) :: prunmode
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
    lscndh2oae  = nlcndh2oae
    lscndh2ocl  = nlcndh2ocl
    lscndh2oic  = nlcndh2oic

    lsauto      = nlauto
    lsautosnow  = nlautosnow

    lsactiv     = nlactiv
    lsactbase   = nlactbase
    lsactintst  = nlactintst

    lsicenucl  = nlicenucl
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
          lsicmelt    = .FALSE.

       CASE(2)  ! Spinup period

          lscoag      = .FALSE.
          lsauto      = .FALSE.
          lsautosnow  = .FALSE.

    END SELECT

    ! Ice formation has an additional spinup time (no ice formation => no autoconversion or melting)
    IF (time<icenucl_tstart) THEN
          lsicenucl  = .FALSE.
          lsautosnow = .FALSE.
          lsicmelt   = .FALSE.
    ENDIF

  END SUBROUTINE set_SALSA_runtime


END MODULE mo_salsa_driver
