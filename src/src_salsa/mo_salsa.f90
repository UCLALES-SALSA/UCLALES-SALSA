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
                   pactd,    pw,    prtcl, level, pdn      )

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
         lsfixinc,                  &
         lsicmelt,                  &
         lsdistupdate

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
         ptstep,                       & ! time step [s]
         pdn(kbdim,klev)                 ! air density

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

    !-- Local variables ------------------------------------------------------------------

    INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level

    zpbl(:)=1

    ! Coagulation
    IF (lscoag) &
       CALL coagulation(kproma, kbdim,  klev,                   &
                        paero,  pcloud, pprecp, pice, psnow,    &
                        ptstep, ptemp,  ppres                   )

    ! Condensation
    IF (lscnd) &
       CALL condensation(kproma,  kbdim,    klev,     krow,     &
                         paero,   pcloud,   pprecp,             &
                         pice,    psnow,                        &
                         pc_h2so4, pc_ocnv, pc_ocsv,  pc_hno3,  &
                         pc_nh3, prv, prs, prsi, ptemp, ppres,  &
                         ptstep, zpbl, prtcl                    )

    ! Autoconversion (liquid)
    IF (lsauto) &
         CALL autoconv2(kproma,kbdim,klev, &
                        pcloud, pprecp, ptstep )

    ! Cloud activation
    IF (lsactiv )  &
         CALL cloud_activation(kproma, kbdim, klev,   &
                               ptemp,  ppres, prv,    &
                               prs,    pw,    paero,  &
                               pcloud, pactd          )

    ! Fixed ice number concentration
    IF (lsfixinc) &
          CALL  ice_fixed_NC(kproma, kbdim, klev,   &
                             pcloud,   pice,   &
                             prv,    prsi,    &
                             pdn     )

    ! Ice nucleation
    IF (lsicenucl) &
        CALL ice_nucl_driver(kproma,kbdim,klev,   &
                          paero,pcloud,pprecp,pice,psnow, &
                          ptemp,ppres,prv,prsi,ptstep)

    ! Melting of ice and snow
    IF (lsicmelt) &
         CALL ice_melt(kproma,kbdim,klev,              &
                       pcloud,pice,pprecp,psnow,ptemp)

    ! Snow formation ~ autoconversion from ice
    IF (lsautosnow) &
         CALL autosnow(kproma,kbdim,klev, &
                       pice, psnow, ptstep )

    ! Size distribution bin update
    IF (lsdistupdate ) &
         CALL distr_update(kproma, kbdim, klev,     &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )

  END SUBROUTINE salsa


END MODULE mo_salsa
