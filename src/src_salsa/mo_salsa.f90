MODULE mo_salsa

  ! --------------------------------------------------------------------------
  ! The SALSA subroutine:
  ! paerml has been removed alltogether because it is rather useless
  ! from this point. Replaced with pvols which is used for all the computation
  ! anyway. Apparently paerml is needed by HAM, so the conversion between pvols
  ! and paerml should be done there upon coupling.
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
                   pactd,    pw,       dbg2,    prtcl, time,level      )

    USE mo_salsa_properties
    USE mo_salsa_dynamics
    USE mo_salsa_nucleation
    USE mo_salsa_update
    USE mo_salsa_cloud ! included for autoconversion and cloud activation


    USE mo_submctl, ONLY :      &
         pi6,                       & ! pi/6
         avog,                      &
         in1a,  in2a,               & ! size section and composition indices
         in2b,  fn1a,               &
         fn2a,  fn2b,               &
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
         debug,                     &
         rhosu, msu, mvsu,          & ! properties of compounds
         rhono, mno,                &
         rhonh, mnh,                &
         rhooc, moc,                &
         rhobc, mbc,                &
         rhoss, mss,                &
         rhowa,                     &
         rhodu, mdu,                &
         nlim,                      &
         mvnh                       ! logical switch for recalculation of zdwet


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

    LOGICAL, INTENT(in) :: dbg2!, testswitch

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

                                ! aerosol distribution properties at each grid point
         !pm6rp (kbdim,klev,fn2b)   ,& ! mean mode actual radius (wet radius for soluble mode
                                ! and dry radius for insoluble modes) [cm]
         !pm6dry(kbdim,klev,fn2a)   ,& ! dry radius [cm]
         !pww   (kbdim,klev,fn2b)   ,& ! aerosol water content for each mode [kg(water) m-3(air)]
         !prhop (kbdim,klev,fn2b)   ,& ! mean mode particle density [g cm-3]
         !ppbl(kbdim)               ,&

    TYPE(t_section), INTENT(inout) :: &
         pcloud(kbdim,klev,ncld),     &
         paero(kbdim,klev,fn2b),      &
         pprecp(kbdim,klev,nprc),     &
         pice(kbdim,klev,nice),       &
         psnow(kbdim,klev,nsnw)

    TYPE(t_section), INTENT(out) :: &
         pactd(kbdim,klev,ncld)

    INTEGER, INTENT(in) :: level                         ! thermodynamical level

    INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level
    !REAL :: zrh(kbdim,klev)           ! Relative humidity


    !-- Output variables -----------------------------------------------------------------

    !   in each bin

    !-- Local variables ------------------------------------------------------------------

    LOGICAL::moving_center
    INTEGER :: ii, jj, kk, nn

    REAL :: zww(kbdim,klev,fn2b)
    REAL :: zrhop(kbdim,klev,fn2b)


    REAL ::                       &
         zcore   (kbdim,klev,fn2b),   & ! volume of the core in one size bin
         !zdwet   (kbdim,klev,fn2b),   &  ! EI TARVI??
         zvq     (kbdim,klev,fn2b),   &
         !zlwc    (kbdim,klev,fn2b),   & ! liquid water in a droplet in one size bin
         zddry                          !EI TARVI??

    IF (debug) WRITE(*,*) 'start salsa'

    moving_center=.false.

    zcore = 0.
    zddry = 0.
    zpbl = 1


    ! Nää tehdään jo LESsissä! Mutta tarvittaessa tyhjää binit tässäkin
    ! SALSA-ajurissa tapahtuvan advektiotendenssikikkailun takia
    DO kk = in1a,fn1a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(paero(ii,jj,kk)%numc < nlim) CYCLE
             zddry = (sum(paero(ii,jj,kk)%volc(1:2)) / &
                       paero(ii,jj,kk)%numc/pi6)**(1./3.)
             IF(zddry < 1.e-10) THEN
  WRITE(*,*) 'hep!',paero(ii,jj,kk)%numc,kk
                pc_h2so4(ii,jj)  = pc_h2so4(ii,jj) + paero(ii,jj,kk)%volc(1) * rhosu / msu * avog
                pc_ocsv(ii,jj)   = pc_ocsv(ii,jj) + paero(ii,jj,kk)%volc(2) * rhooc / moc * avog
                pc_hno3(ii,jj)   = pc_hno3(ii,jj) + paero(ii,jj,kk)%volc(6) * rhono / mno * avog
                pc_nh3(ii,jj)    = pc_nh3(ii,jj)  + paero(ii,jj,kk)%volc(7) * rhonh / mnh * avog
                paero(ii,jj,kk)%numc  = 0.0
                paero(ii,jj,kk)%volc(:) = 0.0

 !               STOP 'hep!'
             END IF
          END DO
       END DO
    END DO

    DO kk = in2a,fn2a      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(paero(ii,jj,kk)%numc < nlim) CYCLE
             zddry = (sum(paero(ii,jj,kk)%volc(:)) / &
                       paero(ii,jj,kk)%numc/pi6)**(1./3.)
             IF(zddry < 1.e-10) THEN
                WRITE(*,*) 'hep',paero(ii,jj,kk)%numc,kk,paero(ii,jj,kk)%volc(:), zddry
                pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + paero(ii,jj,kk)%volc(1) * rhosu / msu * avog
                pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + paero(ii,jj,kk)%volc(2) * rhooc / moc * avog
                pc_hno3(ii,jj)   = pc_hno3(ii,jj) + paero(ii,jj,kk)%volc(6) * rhono / mno * avog
                pc_nh3(ii,jj)    = pc_nh3(ii,jj)  + paero(ii,jj,kk)%volc(7) * rhonh / mnh * avog
                paero(ii,jj,kk)%numc  = 0.0
                paero(ii,jj,kk)%volc(:) = 0.0

!                STOP 'hep'
             END IF
          END DO
       END DO
    END DO

    DO kk = in2b,fn2b      ! size bin
       DO jj = 1,klev      ! vertical grid
          DO ii = 1,kproma ! horizontal grid
             IF(paero(ii,jj,kk)%numc < nlim) CYCLE
             zddry = (sum(paero(ii,jj,kk)%volc(:)) / &
                       paero(ii,jj,kk)%numc/pi6)**(1./3.)
             IF(zddry < 1.e-10) THEN
                pc_h2so4(ii,jj)   = pc_h2so4(ii,jj) + paero(ii,jj,kk)%volc(1) * rhosu / msu * avog
                pc_ocsv(ii,jj)    = pc_ocsv(ii,jj) + paero(ii,jj,kk)%volc(2) * rhooc / moc * avog
                pc_hno3(ii,jj)   = pc_hno3(ii,jj) + paero(ii,jj,kk)%volc(6) * rhono / mno * avog
                pc_nh3(ii,jj)    = pc_nh3(ii,jj)  + paero(ii,jj,kk)%volc(7) * rhonh / mnh * avog
                paero(ii,jj,kk)%numc  = 0.0
                paero(ii,jj,kk)%volc(:) = 0.0
                !WRITE(*,*) 'he'
                STOP 'he'
             END IF
          END DO
       END DO
    END DO

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

    zww    = (pi6*paero(:,:,:)%dwet**3 - paero(:,:,:)%core)*rhowa*paero(:,:,:)%numc ! Lisää nitraattien kontrib

    DO jj = 1,klev      ! vertical grid
       DO ii = 1,kproma ! horizontal kproma in the slab

          DO kk = in1a, fn1a
             zvq(ii,jj,kk) = sum(paero(ii,jj,kk)%volc(1:2)) + &
                  zww(ii,jj,kk)/ rhowa

             ! check if enough aerosol mass present

             IF(zvq(ii,jj,kk) > 1.e-30) THEN
                zrhop(ii,jj,in1a:fn1a) = (paero(ii,jj,kk)%volc(1)*rhosu     &
                     + paero(ii,jj,kk)%volc(2)*rhooc + zww(ii,jj,kk)) /     &
                     zvq(ii,jj,kk)                                    &
                     ! conversion from kg/m3 to g/cm3
                     /1000.
             ELSE

                ! if not enough mass, assume density of water in g/cm3 to avoid NaN
                zrhop(ii,jj,kk) = rhowa/1000.

             END IF

          END DO

          DO kk = in2a, fn2a
             zvq(ii,jj,kk) = sum(paero(ii,jj,kk)%volc(1:7))+                  &
                             zww(ii,jj,kk)/rhowa

             IF(zvq(ii,jj,kk) > 1.e-30) THEN

                zrhop(ii,jj,kk) = (paero(ii,jj,kk)%volc(1)*rhosu +          &
                                   paero(ii,jj,kk)%volc(2)*rhooc +          &
                                   paero(ii,jj,kk)%volc(6)*rhono +          &
                                   paero(ii,jj,kk)%volc(7)*rhonh +          &
                                   paero(ii,jj,kk)%volc(3)*rhobc +          &
                                   paero(ii,jj,kk)%volc(4)*rhodu +          &
                                   paero(ii,jj,kk)%volc(5)*rhoss +          &
                                   zww(ii,jj,kk))/                    &
                                   zvq(ii,jj,kk) / 1000.

             ELSE

                zrhop(ii,jj,kk) = rhowa/1000.

             END IF

          END DO

          DO kk = in2b, fn2b

             zvq(ii,jj,kk) = sum(paero(ii,jj,kk)%volc(1:4)) + SUM(paero(ii,jj,kk)%volc(6:7)) + &
                              zww(ii,jj,kk) / rhowa

             IF(zvq(ii,jj,kk) > 1.e-30) THEN
                zrhop(ii,jj,kk) = (paero(ii,jj,kk)%volc(1)*rhosu +          &
                                   paero(ii,jj,kk)%volc(2)*rhooc +          &
                                   paero(ii,jj,kk)%volc(6)*rhono +          &
                                   paero(ii,jj,kk)%volc(7)*rhonh +          &
                                   paero(ii,jj,kk)%volc(3)*rhobc +          &
                                   paero(ii,jj,kk)%volc(4)*rhodu +          &
                                   zww(ii,jj,kk))/                    &
                                   zvq(ii,jj,kk)/1000.
             ELSE
                zrhop(ii,jj,kk) = rhobc/1000.
             END IF

          END DO
       END DO
    END DO

  END SUBROUTINE salsa


END MODULE mo_salsa
