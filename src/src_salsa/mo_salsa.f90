MODULE mo_salsa

  ! --------------------------------------------------------------------------
  ! The SALSA subroutine
  !
  ! Modified for the new aerosol datatype,
  ! Juha Tonttila, FMI, 2014.
  ! --------------------------------------------------------------------------


  PRIVATE

  ! -- subroutines
  PUBLIC :: salsa

CONTAINS

  SUBROUTINE salsa(kbdim,  klev,                        &
                   ppres,  prv,    prs,    prsi,        &
                   ptemp,  ptstep,                      &
                   pc_gas, ngas,                        &
                   paero,  pcloud, pprecp, pice, psnow, &
                   level,                               &
                   coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, &
                   coag_nprecp, coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
                   cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
                   autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
                   act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
                   melt_vice, melt_nice, melt_vsnow, melt_nsnow)

    USE mo_salsa_dynamics, only : coagulation, condensation
    USE mo_salsa_update, ONLY : distr_update
    USE mo_salsa_cloud, only : cloud_activation, autoconv2, autoconv_sb, &
            autosnow, fixed_ice_driver, ice_nucl_driver, ice_melt

    USE mo_submctl, ONLY :      &
         fn2b,ncld,nprc,nice,nsnw,         &
         t_section,                        &
         lscoag,lscnd,                     &
         lsauto,lsautosnow,lsactiv,        &
         lsicenucl,lsicmelt,lsdistupdate,  &
         fixinc, ice_hom, ice_imm, ice_dep

    IMPLICIT NONE

    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         ngas                         ! number of gases


    REAL, INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep                          ! time step [s]

    !-- Input variables that are changed within --------------------------------------
    REAL, INTENT(inout) ::      &
         pc_gas(kbdim,klev,ngas),   & ! gas phase concentrations at each grid point [mol/m3]
         prv(kbdim,klev),           & ! Water vapour mixing ratio  [kg/kg]
         prs(kbdim,klev),           & ! Saturation mixing ratio    [kg/kg]
         prsi(kbdim,klev)             ! Saturation mixing ratio over ice   [kg/kg]

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &
         paero(kbdim,klev,fn2b),      &
         pprecp(kbdim,klev,nprc),     &
         pice(kbdim,klev,nice),       &
         psnow(kbdim,klev,nsnw)

    INTEGER, INTENT(in) :: level      ! thermodynamical level

    !-- Output statistics --------------------------------------
    REAL, DIMENSION(kbdim,klev), INTENT(out) :: &
            coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, coag_nprecp, &
            coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
            cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
            autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
            act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
            melt_vice, melt_nice, melt_vsnow, melt_nsnow

    ! Reset outputs
    coag_vaero=0.; coag_naero=0.; coag_vcloud=0.; coag_ncloud=0.; coag_vprecp=0.; coag_nprecp=0.
    coag_vice=0.; coag_nice=0.; coag_vsnow=0.; coag_nsnow=0.
    cond_vaero=0.; cond_vcloud=0.; cond_vprecp=0.; cond_vice=0.; cond_vsnow=0.
    autoc_vprecp=0.; autoc_nprecp=0.; autoc_vsnow=0.; autoc_nsnow=0.
    act_vcloud=0.; act_ncloud=0.; nucl_vice=0.; nucl_nice=0.;
    melt_vice=0.; melt_nice=0.; melt_vsnow=0.; melt_nsnow=0.

    ! Coagulation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    IF (lscoag) THEN
       coag_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)
       coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)
       coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
       coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
       coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
       coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
       coag_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
       coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
       coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
       coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
       CALL coagulation(kbdim,  klev,                           &
                        paero,  pcloud, pprecp, pice, psnow,    &
                        ptstep, ptemp,  ppres                   )
       coag_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)-coag_vaero(:,:)
       coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)-coag_naero(:,:)
       coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-coag_vcloud(:,:)
       coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-coag_ncloud(:,:)
       coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-coag_vprecp(:,:)
       coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-coag_nprecp(:,:)
       coag_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-coag_vice(:,:)
       coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-coag_nice(:,:)
       coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-coag_vsnow(:,:)
       coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-coag_nsnow(:,:)
    ENDIF

    ! Condensation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    IF (lscnd) THEN
        cond_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)
        cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
        cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
        cond_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        CALL condensation(kbdim,  klev,                         &
                          paero,  pcloud, pprecp, pice, psnow,  &
                          pc_gas, ngas,                         &
                          prv, prs, prsi, ptemp, ppres, ptstep  )
        cond_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)-cond_vaero(:,:)
        cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-cond_vcloud(:,:)
        cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-cond_vprecp(:,:)
        cond_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-cond_vice(:,:)
        cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-cond_vsnow(:,:)
    ENDIF

    ! Autoconversion (liquid)
    !   Statistics: change in total rain water volume (=change in cloud water) and rain drop number concentration
    IF (lsauto) THEN
         autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
         autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
         CALL autoconv2(kbdim,klev,pcloud, pprecp)
         !CALL autoconv_sb(kbdim,klev,ptstep,pcloud,pprecp)
         autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-autoc_vprecp(:,:)
         autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-autoc_nprecp(:,:)
    ENDIF

    ! Cloud activation
    !   Statistics: change in total cloud water volume (=change in cloud water) and cloud drop number concentration
    IF (lsactiv ) THEN
         act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
         act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
         CALL cloud_activation(kbdim,  klev,          &
                               ptemp,  ppres, prv,    &
                               prs,    paero, pcloud  )
         act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-act_vcloud(:,:)
         act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-act_ncloud(:,:)
    ENDIF

    ! Ice nucleation
    !   Statistics: change in total ice/snow* water volume and ice/snow number concentration
    !   * ice nucleation can also produce snow and in this case snow formation rate is saved
    !     to the autoconversion variables (no ice category; autoconversion disabled)
    IF (lsicenucl .AND. fixinc>=0.) THEN
        ! Fixed ice number concentration
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
        CALL fixed_ice_driver(kbdim, klev,             &
                             pcloud, pice,   psnow,    &
                             ptemp,  ppres,  prv,  prsi)
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-nucl_vice(:,:)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-nucl_nice(:,:)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
    ELSEIF (lsicenucl .AND. (ice_hom .OR. ice_imm .OR. ice_dep)) THEN
        ! Modelled ice nucleation
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
        CALL ice_nucl_driver(kbdim,klev,   &
                          paero,pcloud,pprecp,pice,psnow, &
                          ptemp,prv,prs,prsi,ptstep)
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-nucl_vice(:,:)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-nucl_nice(:,:)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
    ENDIF

    ! Melting of ice and snow
    !   Statistics: change in total ice and snow water volume and number concentrations
    IF (lsicmelt) THEN
         melt_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
         melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
         melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
         melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
         CALL ice_melt(kbdim,klev,pcloud,pice,pprecp,psnow,ptemp)
         melt_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-melt_vice(:,:)
         melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-melt_nice(:,:)
         melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-melt_vsnow(:,:)
         melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-melt_nsnow(:,:)
    ENDIF

    ! Snow formation ~ autoconversion from ice
    !   Statistics: change in total snow water volume (=change in ice water) and snow number concentration
    IF (lsautosnow) THEN
         autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
         autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
         CALL autosnow(kbdim,klev,pice,psnow)
         autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
         autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
    ENDIF

    ! Size distribution bin update
    IF (lsdistupdate ) &
         CALL distr_update(kbdim, klev,             &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )

  END SUBROUTINE salsa


END MODULE mo_salsa
