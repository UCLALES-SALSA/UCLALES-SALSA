 
MODULE radiation_main
  USE grid, ONLY : iradtyp, nxp, nyp, nzp, dzt,       &
                   dn0, pi0, pi1, sst, cntlat, CCN,   &
                   a_pexnr, a_temp, a_tt,             &
                   a_rv, a_rp, a_rc, a_ri,            &
                   a_npp, a_rpp,                      &
                   a_ncloudp, a_nprecpp, a_mprecpp,   &
                   a_nicep,                           &
                   a_rflx, a_sflx,                    &
                   a_fus, a_fds,                      &
                   a_fuir, a_fdir,                    &
                   albedo, prtcl, level

  USE mo_submctl, ONLY : nprc,ira,fra
  USE class_ComponentIndex, ONLY : getNComp
  USE radiation, ONLY : d4stream
  IMPLICIT NONE
  
  REAL    :: sfc_albedo = 0.05
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadConstPress = .FALSE. ! Keep constant pressure levels
  INTEGER :: RadPrecipBins = 0 ! Add precipitation bins to cloud water (for level 3 and up)
  CHARACTER (len=50) :: radsounding = 'datafiles/dsrt.lay'  ! Juha: Added so the radiation background sounding can be given
                                                            ! from the NAMELIST
  
  CONTAINS

    SUBROUTINE rad_interface(time_in)
      IMPLICIT NONE
      
      REAL, INTENT(in) :: time_in
      
      INTEGER :: nspec
      REAL :: znc(nzp,nxp,nyp), zrc(nzp,nxp,nyp), zni(nzp,nxp,nyp), zri(nzp,nxp,nyp)
      
      nspec = getNcomp(prtcl)
      
      ! Radiation
      ! -------------
      
      !
      ! Level 1-3
      ! ---------
      IF (level <= 3) THEN
         znc(:,:,:) = CCN
         zrc(:,:,:) = a_rc(:,:,:) ! Cloud water only
         IF (level == 3 .AND. RadPrecipBins > 0) THEN ! Add precipitation (all or nothing)
            znc(:,:,:) = znc(:,:,:) + a_npp(:,:,:)
            zrc(:,:,:) = zrc(:,:,:) + a_rpp(:,:,:)
         END IF
         CALL d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
              dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rv, zrc, znc, a_tt,  &
              a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
              useMcICA=useMcICA, ConstPrs=RadConstPress)
         
      !
      ! Level 4
      ! -----------
      ELSE IF (level == 4) THEN
         znc(:,:,:) = SUM(a_ncloudp(:,:,:,:),DIM=4) ! Cloud droplets
         zrc(:,:,:) = a_rc(:,:,:) ! Cloud and aerosol water
         IF (RadPrecipBins > 0) THEN ! Add precipitation bins
            ! Water is the last species (nspec+1)
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,nspec*nprc+ira:nspec*nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         CALL d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
                       dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
                       a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
                       useMcICA=useMcICA, ConstPrs=RadConstPress)
         
      !
      ! Level 5
      ! ----------
      ELSE IF (level == 5) THEN
         znc(:,:,:) = SUM(a_ncloudp(:,:,:,:),DIM=4) ! Cloud droplets
         zrc(:,:,:) = a_rc(:,:,:) ! Cloud and aerosol water
         IF (RadPrecipBins > 0) THEN ! Add precipitation bins
            ! Water is the last species (nspec+1)
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,nspec*nprc+ira:nspec*nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         zni(:,:,:) = SUM(a_nicep(:,:,:,:),DIM=4) ! Ice
         zri(:,:,:) = a_ri(:,:,:) ! Ice (no aerosol ice?)
         CALL d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
                       dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
                       a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, ice=zri,nice=zni,radsounding=radsounding, &
                       useMcICA=useMcICA, ConstPrs=RadConstPress)
      END IF


    END SUBROUTINE rad_interface
     

END MODULE radiation_main
