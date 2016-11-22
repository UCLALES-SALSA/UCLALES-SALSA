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
                   ppres,    prv, prs, prsi,    ptemp, ptt, ptstep,     &
                   pc_h2so4, pc_ocnv,  pc_ocsv, pc_hno3,    &
                   pc_nh3,   paero,    pcloud,  pprecp,     &
                   pice, psnow,                             &
                   pactd,    pw,    prtcl, time,level      )

    USE mo_salsa_dynamics, only : coagulation, condensation
    USE mo_salsa_update, ONLY : distr_update
    USE mo_salsa_cloud, only : cloud_activation, autoconv2, &
              ice_immers_nucl,ice_hom_nucl,ice_het_nucl,ice_melt, &
              autosnow

    USE mo_submctl, ONLY :      &
         fn2b,               & ! size section and composition indices
         t_section,                 & ! For cloud bins
         ncld,                      &
         nprc,                      &
         nice,                      & ! ice
         nsnw,                      & ! snow
         lscoag,                    &
         lscnd,                     &
         lsauto,                    &
         lsautosnow,                &
         lsactiv,                   &
         lsichom,                   &
         lsichet,                   &
         lsicimmers,                &
         lsicmelt,                  &
         lsdistupdate,              &
         debug

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
         ptt(kbdim,klev),              & ! temperature tendency
         ptstep,                       &   ! time step [s]
         time                             ! time

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
         prv(kbdim,klev),           & ! Water vapour mixing ratio  [kg/m3]
         prs(kbdim,klev),           & ! Saturation mixing ratio    [kg/m3]
         prsi(kbdim,klev)              ! Saturation mixing ratio over ice   [kg/m3]

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

    IF (debug) WRITE(*,*) 'start salsa'

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
                        pcloud, pprecp     )

    ! Cloud activation
    IF (lsactiv )  &
         CALL cloud_activation(kproma, kbdim, klev,   &
                               ptemp,  ppres, prv,    &
                               prs,    pw,    paero,  &
                               pcloud, pactd          )

    ! Immersion freezing
    IF (lsicimmers) &
         CALL ice_immers_nucl(kproma,kbdim,klev,            &
                              pcloud,pice,ppres,            &
                              ptemp,ptt,prv,prs,ptstep,time )

    ! Homogenous nucleation Morrison et al. 2005 eq. (25)
    IF (lsichom) &
        CALL ice_hom_nucl(kproma,kbdim,klev,       &
                          pcloud,pice,paero,ppres, &
                          ptemp,prv,prs,ptstep)

    !! heterogenous nucleation Morrison et al. 2005 eq. (27)
    IF (lsichet) &
        CALL ice_het_nucl(kproma,kbdim,klev,       &
                          pcloud,pice,paero,ppres, &
                          ptemp,prv,prs,ptstep)
    
    ! Melting of ice
    IF (lsicmelt) &
         CALL ice_melt(kproma,kbdim,klev,              &
                       pcloud,pice,pprecp,psnow,ppres, &
                       ptemp,prv,prs,ptstep)

    ! Snow formation ~ autoconversion for ice
    IF (lsautosnow) &
         CALL autosnow(kproma,kbdim,klev, &
                       pice, psnow        )

    ! Size distribution bin update
    IF (lsdistupdate ) &
         CALL distr_update(kproma, kbdim, klev,     &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )

  END SUBROUTINE salsa


END MODULE mo_salsa
