
MODULE radiation_main
  USE grid, ONLY : iradtyp, nxp, nyp, nzp, dzt,       &
                   dn0, pi0, pi1, sst, cntlat, CCN,   &
                   a_pexnr, a_temp, a_tt,             &
                   a_rv, a_rp, a_rc, a_ri,            &
                   a_npp, a_rpp,                      &
                   a_maerop, a_naerop,                &
                   a_ncloudp, a_nprecpp, a_mprecpp,   &
                   a_nicep, a_nsnowp, a_msnowp,       &
                   a_rflx, a_sflx,                    &
                   a_fus, a_fds,                      &
                   a_fuir, a_fdir,                    &
                   albedo, level

  USE mo_submctl, ONLY : nprc,nsnw,ira,fra,isa,fsa,spec
  USE radiation, ONLY : d4stream
  IMPLICIT NONE
  
  REAL    :: sfc_albedo = 0.1
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadConstPress = .FALSE. ! Keep constant pressure levels
  INTEGER :: RadPrecipBins = 0 ! Add precipitation bins to cloud water (for level 3 and up)
  INTEGER :: RadSnowBins = 0 ! Add snow bins to cloud ice (for level 5 and up)
  CHARACTER (len=50) :: radsounding = 'datafiles/dsrt.lay'  ! Juha: Added so the radiation background sounding can be given
                                                            ! from the NAMELIST
  LOGICAL :: laerorad = .FALSE. ! Use binned aerosol for radiation? Only for level > 4

  CONTAINS

    SUBROUTINE rad_interface(time_in)
      IMPLICIT NONE
      
      REAL, INTENT(in) :: time_in  ! Time, decimal days
      
      INTEGER :: nspec
      REAL :: znc(nzp,nxp,nyp), zrc(nzp,nxp,nyp), zni(nzp,nxp,nyp), zri(nzp,nxp,nyp)
      
      nspec = spec%getNSpec()  
      nspec = MAX(1,nspec) ! This is to avoid some problems with LEV3 even though not actually even used. 
                           ! Avoids the need for more switches

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
         CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo, &
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
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,(nspec-1)*nprc+ira:(nspec-1)*nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         IF (laerorad) THEN
            CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo, &
                          dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
                          a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
                          useMcICA=useMcICA, ConstPrs=RadConstPress, maerop=a_maerop, naerop=a_naerop)
         ELSE
            CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo, &
                          dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
                          a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
                          useMcICA=useMcICA, ConstPrs=RadConstPress)
         END IF

      !
      ! Level 5
      ! ----------
      ELSE IF (level == 5) THEN
         znc(:,:,:) = SUM(a_ncloudp(:,:,:,:),DIM=4) ! Cloud droplets
         zrc(:,:,:) = a_rc(:,:,:) ! Cloud and aerosol water
         IF (RadPrecipBins > 0) THEN ! Add precipitation bins
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,(nspec-1)*nprc+ira:(nspec-1)*nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         zni(:,:,:) = SUM(a_nicep(:,:,:,:),DIM=4) ! Ice
         zri(:,:,:) = a_ri(:,:,:) ! Ice (no aerosol ice?)
         IF (RadSnowBins>0) THEN ! Add snow bins
            zri(:,:,:) = zri(:,:,:) + SUM(a_msnowp(:,:,:,(nspec-1)*nsnw+isa:(nspec-1)*nsnw+min(RadSnowBins,fsa)),DIM=4)
            zni(:,:,:) = zni(:,:,:) + SUM(a_nsnowp(:,:,:,isa:min(RadSnowBins,fsa)),DIM=4)
         END IF
         CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo, &
                       dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
                       a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, ice=zri,nice=zni,radsounding=radsounding, &
                       useMcICA=useMcICA, ConstPrs=RadConstPress, maerop=a_maerop, naerop=a_naerop)
      END IF


    END SUBROUTINE rad_interface
     

END MODULE radiation_main
