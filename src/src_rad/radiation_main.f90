
MODULE radiation_main
  USE mo_aux_state, ONLY : dzt, dn0, pi0, pi1
  USE mo_diag_state, ONLY : a_pexnr, a_temp, a_rv, a_rc, a_ri, a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo
  USE mo_progn_state, ONLY : a_tt, a_rp, a_npp, a_rpp, a_maerop, a_naerop, a_ncloudp, a_nprecpp, a_mprecpp, a_nicep
  USE grid, ONLY : iradtyp, nxp, nyp, nzp,       &
                   sst, cntlat, CCN, level

  USE mo_submctl, ONLY : nprc,ira,fra,spec
  USE radiation, ONLY : d4stream
  IMPLICIT NONE
  
  REAL    :: sfc_albedo = 0.1
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadConstPress = .FALSE. ! Keep constant pressure levels
  INTEGER :: RadPrecipBins = 0 ! Add precipitation bins to cloud water (for level 3 and up)
  CHARACTER (len=50) :: radsounding = 'datafiles/dsrt.lay'  ! Juha: Added so the radiation background sounding can be given
                                                            ! from the NAMELIST
  LOGICAL :: laerorad = .FALSE. ! Use binned aerosol for radiation? Only for level > 4

  CONTAINS

    SUBROUTINE rad_interface(time_in)
      IMPLICIT NONE
      
      REAL, INTENT(in) :: time_in  ! Time, decimal days
      
      INTEGER :: nspec
      REAL :: znc(nzp,nxp,nyp), zrc(nzp,nxp,nyp), zni(nzp,nxp,nyp), zri(nzp,nxp,nyp)
      
      nspec = spec%getNSpec(type="wet")  
      nspec = MAX(1,nspec) ! This is to avoid some problems with LEV3 even though not actually even used. 
                           ! Avoids the need for more switches
      
      ! Radiation
      ! -------------

      ! NOTE: Provide the arguments to radiation submodel as FLOATS for now, i.e. do not use the FloatArray datatypes!
      
      !
      ! Level 1-3
      ! ---------
      IF (level <= 3) THEN
         znc(:,:,:) = CCN
         zrc(:,:,:) = a_rc%d(:,:,:) ! Cloud water only
         IF (level == 3 .AND. RadPrecipBins > 0) THEN ! Add precipitation (all or nothing)
            znc(:,:,:) = znc(:,:,:) + a_npp%d(:,:,:)
            zrc(:,:,:) = zrc(:,:,:) + a_rpp%d(:,:,:)
         END IF
         CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo, &
              dn0%d, pi0%d, pi1%d, dzt%d, a_pexnr%d, a_temp%d, a_rv%d, zrc, znc, a_tt%d,  &
              a_rflx%d, a_sflx%d, a_fus%d, a_fds%d, a_fuir%d, a_fdir%d, albedo%d, radsounding=radsounding, &
              useMcICA=useMcICA, ConstPrs=RadConstPress)
         
      !
      ! Level 4
      ! -----------
      ELSE IF (level == 4) THEN
         znc(:,:,:) = SUM(a_ncloudp%d(:,:,:,:),DIM=4) ! Cloud droplets
         zrc(:,:,:) = a_rc%d(:,:,:) ! Cloud and aerosol water
         IF (RadPrecipBins > 0) THEN ! Add precipitation bins            
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp%d(:,:,:,(nspec-1)*nprc+ira:(nspec-1) *   &
                                          nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp%d(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         IF (laerorad) THEN
            CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo,   &
                          dn0%d, pi0%d, pi1%d, dzt%d, a_pexnr%d, a_temp%d, a_rp%d,  &
                          zrc, znc, a_tt%d, a_rflx%d, a_sflx%d, a_fus%d, a_fds%d,   &
                          a_fuir%d, a_fdir%d, albedo%d, radsounding=radsounding,    &
                          useMcICA=useMcICA, ConstPrs=RadConstPress,                &
                          maerop=a_maerop%d, naerop=a_naerop%d)
         ELSE
            CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo,   &
                          dn0%d, pi0%d, pi1%d, dzt%d, a_pexnr%d, a_temp%d, a_rp%d,  &
                          zrc, znc, a_tt%d, a_rflx%d, a_sflx%d, a_fus%d, a_fds%d,   &
                          a_fuir%d, a_fdir%d, albedo%d, radsounding=radsounding,    &
                          useMcICA=useMcICA, ConstPrs=RadConstPress)
         END IF

      !
      ! Level 5
      ! ----------
      ELSE IF (level == 5) THEN
         znc(:,:,:) = SUM(a_ncloudp%d(:,:,:,:),DIM=4) ! Cloud droplets
         zrc(:,:,:) = a_rc%d(:,:,:) ! Cloud and aerosol water
         IF (RadPrecipBins > 0) THEN ! Add precipitation bins
            zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp%d(:,:,:,(nspec-1)*nprc+ira:(nspec-1) *   &
                                          nprc+min(RadPrecipBins,fra)),DIM=4)
            znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp%d(:,:,:,ira:min(RadPrecipBins,fra)),DIM=4)
         END IF
         zni(:,:,:) = SUM(a_nicep%d(:,:,:,:),DIM=4) ! Ice
         zri(:,:,:) = a_ri%d(:,:,:) ! Ice (no aerosol ice?)
         CALL d4stream(nzp, nxp, nyp, nspec, cntlat, time_in, sst, sfc_albedo,   &
                       dn0%d, pi0%d, pi1%d, dzt%d, a_pexnr%d, a_temp%d, a_rp%d,  &
                       zrc, znc, a_tt%d, a_rflx%d, a_sflx%d, a_fus%d, a_fds%d,   &
                       a_fuir%d, a_fdir%d, albedo%d, ice=zri,nice=zni,           &
                       radsounding=radsounding,useMcICA=useMcICA,                &
                       ConstPrs=RadConstPress, maerop=a_maerop%d,                &
                       naerop=a_naerop%d)
      END IF


    END SUBROUTINE rad_interface
     

END MODULE radiation_main
