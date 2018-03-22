MODULE mo_salsa
   IMPLICIT NONE

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
                    ppres,    prv, prs, prsi,    ptemp,  ptstep,     &
                    pc_h2so4, pc_ocnv,  pc_ocsv, pc_hno3,    &
                    pc_nh3,   allSALSA,                         &
                    pactd,    pw, level )
     
     USE classSection
     USE mo_salsa_dynamics, only : coagulation, condensation
     USE mo_salsa_update, ONLY : distr_update
     USE mo_salsa_cloud, only : cloud_activation, autoconv2, &
                                autosnow,ice_fixed_NC, ice_nucl_driver,ice_melt
     
     USE mo_submctl, ONLY :      &
          fn2b,                      & ! size section and composition indices
          ncld,                      &
          nprc,                      &
          nice,                      &
          nsnw,                      & 
          ntotal,                    &
          fixinc,                    &
          lscoag,                    &
          lscnd,                     &
          lsauto,                    &
          lsautosnow,                &
          lsactiv,                   &
          lsicenucl,                 &
          lsicemelt,                  &
          lsdistupdate,              &
          lscheckarrays,             &
          ice_hom, ice_imm, ice_dep

     IMPLICIT NONE

     !-- Input parameters and variables --------------------------------------------------
     INTEGER, INTENT(in) ::      &
          kproma,                    & ! number of horiz. grid points in a slab (='kproma')
          kbdim,                     & ! dimension for arrays (='kbdim')
          klev,                      & ! number of vertical levels (='klev')
          krow                         ! local latitude index
     
     REAL, INTENT(in) ::            &
          ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
          ptemp(kbdim,klev),            & ! temperature at each grid point [K]
          ptstep                          ! time step [s]
     
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
     
     TYPE(Section), INTENT(inout) :: &
          allSALSA(kbdim,klev,ntotal)
     
     TYPE(Section), INTENT(out) :: &
          pactd(kbdim,klev,ncld)
     
     INTEGER, INTENT(in) :: level                         ! thermodynamical level
     
     !-- Local variables ------------------------------------------------------------------
     
     INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level
     
     zpbl(:) = 1
     
     ! Coagulation
     IF (lscoag%state) &
          CALL coagulation( kproma, kbdim,  klev,                   &
                            allSALSA,    &
                            ptstep, ptemp,  ppres   )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"COAG")

     ! Condensation
     IF (lscnd%state) &
          CALL condensation(kproma, kbdim,    klev,     krow,          &
                            level, allSALSA,                    &
                            pc_h2so4, pc_ocnv, pc_ocsv,  pc_hno3,  &
                            pc_nh3, prv, prs, prsi, ptemp, ppres,  &
                            ptstep, zpbl                    )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"CONDENSATION")

     ! Autoconversion (liquid)
     IF (lsauto%state) &
          CALL autoconv2(kproma,kbdim,klev, &
                         ptstep   )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"AUTOCONV")

     ! Cloud activation
     IF (lsactiv%state )  &
          CALL cloud_activation(kproma, kbdim, klev,   &
                                ptemp,  ppres, prv,    &
                                prs,    pw           , &
                                pactd          )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"ACTIVATION")

     ! Ice nucleation
     IF (lsicenucl%state) THEN
        IF (fixinc>0. .AND. .NOT. ANY([ice_hom,ice_imm,ice_dep])) THEN
           ! Fixed ice number concentration
           CALL  ice_fixed_NC(kproma, kbdim, klev,   &
                              ptemp,  ppres,  prv,  prsi)
        ELSE IF (ANY([ice_hom,ice_imm,ice_dep])) THEN
           ! Modelled ice nucleation
           CALL ice_nucl_driver(kproma,kbdim,klev,   &
                                ptemp,prv,prs,prsi,ptstep)
        END IF
     END IF

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"ICENUC")

     ! Melting of ice and snow
     IF (lsicemelt%state) &
          CALL ice_melt(kproma,kbdim,klev,              &
                        ptemp)

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"ICEMELT")

     ! Snow formation ~ autoconversion from ice
     IF (lsautosnow%state) &
          CALL autosnow(kproma,kbdim,klev)

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"AUTOSNOW")

     ! Size distribution bin update
     IF (lsdistupdate ) &
          CALL distr_update(kproma, kbdim, klev,     &
                            allSALSA, level       )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,ntotal,allSALSA,"DISTUPDATE")

   END SUBROUTINE salsa

   ! -------------------------------

   SUBROUTINE check_arrays(kbdim,klev,ntotal,array,position)
     USE classSection, ONLY : Section
     IMPLICIT NONE
     ! Check that particle arrays remain positive and
     ! check for NANs
     INTEGER, INTENT(in) :: kbdim,klev,ntotal
     CHARACTER(len=*), INTENT(in) :: position
     TYPE(Section), INTENT(in) :: array(kbdim,klev,ntotal)
     
     INTEGER :: ii,jj,nn

     DO jj = 1,klev
        DO ii = 1,kbdim
           DO nn = 1,ntotal
              IF ( array(ii,jj,nn)%numc < 0. .OR.  &
                   ANY( array(ii,jj,nn)%volc(:) < 0. ) ) THEN

                 WRITE(*,*) 'SALSA: NEGATIVE CONCENTRATIONS, BIN ',nn
                 WRITE(*,*) 'NUMC', array(ii,jj,nn)%numc
                 WRITE(*,*) 'VOLC', array(ii,jj,nn)%volc(:)
                 WRITE(*,*) 'At '//TRIM(position)

              END IF

              IF ( array(ii,jj,nn)%numc /= array(ii,jj,nn)%numc  .OR.  &
                   ANY( array(ii,jj,nn)%volc(:) /= array(ii,jj,nn)%volc(:) ) ) THEN

                   WRITE(*,*) 'SALSA: NAN CONCENTRATIONS, BIN ',nn
                   WRITE(*,*) 'NUMC', array(ii,jj,nn)%numc
                   WRITE(*,*) 'VOLC', array(ii,jj,nn)%volc(:)
                   WRITE(*,*) 'At '//TRIM(position)
                   
              END IF
           END DO
        END DO
     END DO

   END SUBROUTINE check_arrays


END MODULE mo_salsa
