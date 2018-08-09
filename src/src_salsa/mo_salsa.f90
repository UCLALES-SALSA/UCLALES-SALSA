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

  SUBROUTINE salsa(kproma,   kbdim,    klev,    krow,       &
                   ppres,    prv, prs, prsi,    ptemp, ptstep,     &
                   pc_h2so4, pc_ocnv,  pc_ocsv, pc_hno3,    &
                   pc_nh3,   paero,    pcloud,  pprecp,     &
                   pice, psnow,                             &
                   pactd,    pw,    prtcl, level,           &
                   coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, &
                   coag_nprecp, coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
                   cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
                   autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
                   act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
                   melt_vice, melt_nice, melt_vsnow, melt_nsnow)

    USE mo_salsa_dynamics, only : coagulation, condensation
    USE mo_salsa_update, ONLY : distr_update
    USE mo_salsa_cloud, only : cloud_activation, autoconv2, &
            autosnow, ice_fixed_NC, ice_nucl_driver, ice_melt

    USE mo_submctl, ONLY :      &
         fn2b,               & ! size section and composition indices
         t_section,                 & ! For cloud bins
         ncld,                      &
         nprc,                      &
         nice,                      &
         nsnw,                      &
         lscoag,                    &
         lscnd,                     &
         lsauto,                    &
         lsautosnow,                &
         lsactiv,                   &
         lsicenucl,                 &
         lsicmelt,                  &
         lsdistupdate,              &
         fixinc, ice_hom, ice_imm, ice_dep

    USE class_componentIndex, ONLY : ComponentIndex

    IMPLICIT NONE

    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kproma,                    & ! number of horiz. grid points in a slab (='kproma')
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         krow                         ! local latitude index


    REAL, INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep                          ! time step [s]

    TYPE(ComponentIndex), INTENT(in) :: prtcl

    REAL, INTENT(in) ::            & ! Vertical velocity
         pw(kbdim,klev)

    !-- Input variables that are changed within --------------------------------------
    REAL, INTENT(inout) ::      & ! gas phase concentrations at each grid point [#/m3]
         pc_h2so4(kbdim,klev),      & ! sulphuric acid
         pc_hno3 (kbdim,klev),      & ! nitric acid
         pc_nh3  (kbdim,klev),      & ! ammonia
         pc_ocnv (kbdim,klev),      & ! nonvolatile organic compound
         pc_ocsv (kbdim,klev),      & ! semivolatile organic compound
         prv(kbdim,klev),           & ! Water vapour mixing ratio  [kg/kg]
         prs(kbdim,klev),           & ! Saturation mixing ratio    [kg/kg]
         prsi(kbdim,klev)             ! Saturation mixing ratio over ice   [kg/kg]

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &
         paero(kbdim,klev,fn2b),      &
         pprecp(kbdim,klev,nprc),     &
         pice(kbdim,klev,nice),       &
         psnow(kbdim,klev,nsnw)

    TYPE(t_section), INTENT(out) :: &
         pactd(kbdim,klev,ncld)

    INTEGER, INTENT(in) :: level                         ! thermodynamical level

    !-- Output statistics --------------------------------------
    REAL, DIMENSION(kbdim,klev), INTENT(out) :: &
            coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, coag_nprecp, &
            coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
            cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
            autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
            act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
            melt_vice, melt_nice, melt_vsnow, melt_nsnow

    !-- Local variables ------------------------------------------------------------------

    INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level

    zpbl(:)=1

    ! Coagulation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    coag_vaero(:,:)=SUM(paero(:,:,:)%volc(8),DIM=3)
    coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)
    coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)
    coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
    coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)
    coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
    coag_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)
    coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
    coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)
    coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
    IF (lscoag) &
       CALL coagulation(kproma, kbdim,  klev,                   &
                        paero,  pcloud, pprecp, pice, psnow,    &
                        ptstep, ptemp,  ppres                   )
    coag_vaero(:,:)=SUM(paero(:,:,:)%volc(8),DIM=3)-coag_vaero(:,:)
    coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)-coag_naero(:,:)
    coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)-coag_vcloud(:,:)
    coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-coag_ncloud(:,:)
    coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)-coag_vprecp(:,:)
    coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-coag_nprecp(:,:)
    coag_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)-coag_vice(:,:)
    coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-coag_nice(:,:)
    coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)-coag_vsnow(:,:)
    coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-coag_nsnow(:,:)

    ! Condensation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    cond_vaero(:,:)=SUM(paero(:,:,:)%volc(8),DIM=3)
    cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)
    cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)
    cond_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)
    cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)
    IF (lscnd) &
       CALL condensation(kproma,  kbdim,    klev,     krow,     &
                         paero,   pcloud,   pprecp,             &
                         pice,    psnow,                        &
                         pc_h2so4, pc_ocnv, pc_ocsv,  pc_hno3,  &
                         pc_nh3, prv, prs, prsi, ptemp, ppres,  &
                         ptstep, zpbl, prtcl                    )
    cond_vaero(:,:)=SUM(paero(:,:,:)%volc(8),DIM=3)-cond_vaero(:,:)
    cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)-cond_vcloud(:,:)
    cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)-cond_vprecp(:,:)
    cond_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)-cond_vice(:,:)
    cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)-cond_vsnow(:,:)

    ! Autoconversion (liquid)
    !   Statistics: change in total rain water volume (=change in cloud water) and rain drop number concentration
    autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)
    autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
    IF (lsauto) &
         CALL autoconv2(kproma,kbdim,klev, &
                        pcloud, pprecp )
    autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(8),DIM=3)-autoc_vprecp(:,:)
    autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-autoc_nprecp(:,:)

    ! Cloud activation
    !   Statistics: change in total cloud water volume (=change in cloud water) and cloud drop number concentration
    act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)
    act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
    IF (lsactiv )  &
         CALL cloud_activation(kproma, kbdim, klev,   &
                               ptemp,  ppres, prv,    &
                               prs,    pw,    paero,  &
                               pcloud, pactd          )
    act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(8),DIM=3)-act_vcloud(:,:)
    act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-act_ncloud(:,:)

    ! Ice nucleation
    !   Statistics: change in total ice water volume and ice number concentration
    nucl_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)
    nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
    IF (lsicenucl .AND. fixinc>=0.) THEN
        ! Fixed ice number concentration
        CALL  ice_fixed_NC(kproma, kbdim, klev,   &
                             pcloud,   pice,   &
                             ptemp,  ppres,  prv,  prsi)
    ELSEIF (lsicenucl .AND. (ice_hom .OR. ice_imm .OR. ice_dep)) THEN
        ! Modelled ice nucleation
        CALL ice_nucl_driver(kproma,kbdim,klev,   &
                          paero,pcloud,pprecp,pice, &
                          ptemp,prv,prs,prsi,ptstep)
    ENDIF
    nucl_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)-nucl_vice(:,:)
    nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-nucl_nice(:,:)

    ! Melting of ice and snow
    !   Statistics: change in total ice and snow water volume and number concentrations
    melt_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)
    melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
    melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)
    melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
    IF (lsicmelt) &
         CALL ice_melt(kproma,kbdim,klev,              &
                       pcloud,pice,pprecp,psnow,ptemp)
    melt_vice(:,:)=SUM(pice(:,:,:)%volc(8),DIM=3)-melt_vice(:,:)
    melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-melt_nice(:,:)
    melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)-melt_vsnow(:,:)
    melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-melt_nsnow(:,:)

    ! Snow formation ~ autoconversion from ice
    !   Statistics: change in total snow water volume (=change in ice water) and snow number concentration
    autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)
    autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
    IF (lsautosnow) &
         CALL autosnow(kproma,kbdim,klev, &
                       pice, psnow )
    autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(8),DIM=3)-autoc_vsnow(:,:)
    autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)

    ! Size distribution bin update
    IF (lsdistupdate ) &
         CALL distr_update(kproma, kbdim, klev,     &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )

  END SUBROUTINE salsa


END MODULE mo_salsa
