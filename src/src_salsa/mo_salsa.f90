MODULE mo_salsa
  USE classSection, ONLY : Section
  USE mo_salsa_dynamics, only : coagulation, condensation
  USE mo_salsa_update, ONLY : distr_update
  USE mo_salsa_cloud, only : cloud_activation, autoconv2
  USE mo_salsa_cloud_ice, ONLY : &
       ice_fixed_NC, ice_nucl => ice_nucl_driver,ice_melt
  USE mo_salsa_cloud_ice_SE, ONLY : ice_nucl_SE => ice_nucl_driver
  USE mo_salsa_secondary_ice, ONLY: secondary_ice 
  USE mo_submctl, ONLY :      &
       spec, &
       ncld,                      &
       ntotal,                    &
       fixinc,                    &
       lscoag,                    &
       lscnd,                     &
       lsauto,                    &
       lsactiv,                   &
       lsicenucl,                 &
       lsicemelt,                 &
       lssecice,                  &
       lsdistupdate,              &
       lscheckarrays,             &
       lsicehom, lsiceimm, lsicedep, &
       ice_theta_dist
  USE mo_salsa_types, ONLY : allSALSA,ice
  
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

   SUBROUTINE salsa(kproma,   kbdim,    klev,    krow,             &
                    ppres,    prv, prs, prsi,    ptemp,  ptstep,   &
                    pc_h2so4, pc_ocnv,  pc_ocsv, pc_hno3,          & 
                    pc_nh3,   pactd,    pw,      level, lcharge    )
     
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
          
     TYPE(Section), INTENT(out) :: &
          pactd(kbdim,klev,ncld)
     
     INTEGER, INTENT(in) :: level                         ! thermodynamical level
     
     LOGICAL, INTENT(in) :: lcharge

     !-- Local variables ------------------------------------------------------------------
     
     INTEGER :: zpbl(kbdim)            ! Planetary boundary layer top level
     
     zpbl(:) = 1
     
     ! Coagulation
     IF (lscoag%state) &
          CALL coagulation( kproma, kbdim,  klev,    &
                            ptstep, ptemp,  ppres,   &                        
                            lcharge                  )
     
     IF (lscheckarrays) CALL check_arrays(kbdim,klev,"COAG")

     ! Secondary ice production
     IF (lssecice%state) CALL secondary_ice(kbdim, kproma, klev, ppres, ptemp, ptstep, prv, prsi)
     
     ! Condensation
     IF (lscnd%state) &
          CALL condensation(kproma,   kbdim,    klev,     krow,      &
                            pc_h2so4, pc_ocnv,  pc_ocsv,  pc_hno3,   &
                            pc_nh3,   prv,      prs,      prsi,      &
                            ptemp,    ppres,    ptstep,   zpbl       )

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,"CONDENSATION")

     ! Autoconversion (liquid)
     IF (lsauto%state .AND. lsauto%mode == 2) &
          CALL autoconv2(kproma,kbdim,klev, ptstep)

     ! Cloud activation
     IF (lsactiv%state )  &
          CALL cloud_activation(kproma, kbdim, klev,   &
                                ptemp,  ppres, prv,    &
                                prs,    pw,    pactd   )

     ! Ice nucleation
     IF (lsicenucl%state .OR. lsicenucl%mode == 2) THEN ! If mode=2, call even if state=false
        IF (fixinc>0. .AND. .NOT. ANY([lsicehom,lsiceimm,lsicedep])) THEN
           ! Fixed ice number concentration
           CALL  ice_fixed_NC(kproma, kbdim, klev,   &
                              ptemp,  ppres,  prv,  prsi)
        ELSE IF (ANY([lsicehom,lsiceimm,lsicedep])) THEN
           ! Modelled ice nucleation
           IF (ice_theta_dist) THEN
              CALL ice_nucl_SE(kproma,kbdim,klev,       &
                               ptemp,prv,prs,prsi,ptstep)
           ELSE
              CALL ice_nucl(kproma,kbdim,klev,        &
                            ptemp,prv,prs,prsi,ptstep )
           END IF              
        END IF
     END IF

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,"ICENUC")

     ! Melting of ice and snow
     IF (lsicemelt%state) &
          CALL ice_melt(kproma,kbdim,klev,ptemp)

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,"ICEMELT")
     
     ! Size distribution bin update
     IF (lsdistupdate ) &
          CALL distr_update(kbdim, klev, level) ! kproma

     IF (lscheckarrays) CALL check_arrays(kbdim,klev,"DISTUPDATE")

   END SUBROUTINE salsa

   ! -------------------------------

   SUBROUTINE check_arrays(kbdim,klev,position)
     IMPLICIT NONE
     ! Check that particle arrays remain positive and
     ! check for NANs
     INTEGER, INTENT(in) :: kbdim,klev
     CHARACTER(len=*), INTENT(in) :: position

     INTEGER :: ii,jj,nn

     DO jj = 1,klev
        DO ii = 1,kbdim
           DO nn = 1,ntotal
              IF ( allSALSA(ii,jj,nn)%numc < 0. .OR.  &
                   ANY( allSALSA(ii,jj,nn)%volc(:) < 0. ) ) THEN

                 WRITE(*,*) 'SALSA: NEGATIVE CONCENTRATIONS, BIN ',nn
                 WRITE(*,*) 'NUMC', allSALSA(ii,jj,nn)%numc
                 WRITE(*,*) 'VOLC', allSALSA(ii,jj,nn)%volc(:)
                 WRITE(*,*) 'At '//TRIM(position)

              END IF

              IF ( allSALSA(ii,jj,nn)%numc /= allSALSA(ii,jj,nn)%numc  .OR.  &
                   ANY( allSALSA(ii,jj,nn)%volc(:) /= allSALSA(ii,jj,nn)%volc(:) ) ) THEN

                   WRITE(*,*) 'SALSA: NAN CONCENTRATIONS, BIN ',nn
                   WRITE(*,*) 'NUMC', allSALSA(ii,jj,nn)%numc
                   WRITE(*,*) 'VOLC', allSALSA(ii,jj,nn)%volc(:)
                   WRITE(*,*) 'At '//TRIM(position)
                   
              END IF
           END DO
        END DO
     END DO

   END SUBROUTINE check_arrays


END MODULE mo_salsa
