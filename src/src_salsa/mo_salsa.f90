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
                   ptemp,  pedr,   ptstep,              &
                   pc_gas, ngas,                        &
                   paero,  pcloud, pprecp, pice, psnow, &
                   level, sflg, nstat, sdata, slist)

    USE mo_salsa_dynamics, only : coagulation, condgas, gpparth2o
    USE mo_vbs_partition, ONLY : vbs_gas_phase_chem, vbs_condensation
    USE mo_salsa_update, ONLY : distr_update
    USE mo_salsa_cloud, only : cloud_activation, autoconv2, autoconv_sb, &
            autosnow, fixed_ice_driver,  param_ice_driver, ice_nucl_driver, ice_melt, &
            sip_hm, sip_iibr, sip_df

    USE mo_submctl, ONLY :      &
         fn2b,ncld,nprc,nice,nsnw,nvbs,    &
         t_section,part_h2so4,part_ocnv,   &
         lscoag,lscnd,lscndgas,nlcndh2oae, &
         nlcndh2ocl,nlcndh2oic,            &
         lsauto,auto_sb,lsautosnow,lsactiv,&
         lsicenucl,lsicmelt,lsdistupdate,  &
         fixinc, ice_diag, ice_hom, ice_imm, ice_dep, &
         nlsip_hm, nlsip_iibr, nlsip_df

    IMPLICIT NONE

    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         ngas                         ! number of gases


    REAL, INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         pedr(kbdim,klev),             & ! eddy dissipation rate
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

    LOGICAL, INTENT(IN) :: sflg       ! statistics sampling flag
    INTEGER, INTENT(IN) :: nstat      ! number of outputs
    REAL, INTENT(OUT) :: sdata(kbdim,klev,nstat) ! output data array
    CHARACTER(LEN=7), DIMENSION(:), INTENT(IN) :: slist ! names of the output variables

    ! Reset outputs
    sdata(:,:,:) = 0.

    ! Coagulation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    IF (lscoag) THEN
       IF (sflg) CALL salsa_var_stat('coag',0)
       CALL coagulation(kbdim,  klev,                           &
                        paero,  pcloud, pprecp, pice, psnow,    &
                        ptstep, ptemp,  ppres, pedr             )
       IF (sflg) CALL salsa_var_stat('coag',1)
    ENDIF

    ! Secondary ice production
    IF (lscoag .AND. nlsip_hm) THEN !  Hallett-Mossop
        IF (sflg) CALL salsa_var_stat('siph',0)
        CALL sip_hm(kbdim, klev, pice, psnow, ptemp)
        IF (sflg) CALL salsa_var_stat('siph',1)
    ENDIF
    IF (lscoag .AND. nlsip_iibr) THEN ! Ice-ice collisional breakup
        IF (sflg) CALL salsa_var_stat('sipi',0)
        CALL sip_iibr(kbdim, klev, pice, psnow, ptemp)
        IF (sflg) CALL salsa_var_stat('sipi',1)
    ENDIF
    IF (lscoag .AND. nlsip_df) THEN ! Droplet fragmentation during freezing
        IF (sflg) CALL salsa_var_stat('sipd',0)
        CALL sip_df(kbdim, klev, pcloud, pprecp, pice, psnow, ptemp)
        IF (sflg) CALL salsa_var_stat('sipd',1)
    ENDIF

    ! Condensation of H2SO4 and non-volatile organic vapor
    IF (lscnd .AND. lscndgas .AND. (part_h2so4 .OR. part_ocnv)) THEN
        CALL condgas(kbdim, klev, &
            paero, pcloud, pprecp, pice, psnow, &
            pc_gas, ngas, ptemp, ppres, ptstep)
    ENDIF
    ! Condensation of semivolatile organic species and gas phase chemistry
    IF (lscnd .AND. lscndgas .AND. nvbs>0) THEN
        ! Gas phase chemistry (oxidation)
        IF (sflg) CALL salsa_var_stat('oxid',0)
        CALL vbs_gas_phase_chem(kbdim, klev, &
            ptemp, ptstep, pc_gas, ngas)
        IF (sflg) CALL salsa_var_stat('oxid',1)

        ! Condensation
        IF (sflg) CALL salsa_var_stat('ocon',0)
        call vbs_condensation(kbdim, klev, &
            paero, pcloud, pprecp, pice, psnow, &
            pc_gas, ngas, ptemp, ppres, ptstep)
        IF (sflg) CALL salsa_var_stat('ocon',1)
    ENDIF

    ! Condensation of water vapor
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    IF (lscnd .AND. (nlcndh2ocl .OR. nlcndh2oae .OR. nlcndh2oic)) THEN
        IF (sflg) CALL salsa_var_stat('cond',0)
        CALL gpparth2o(kbdim, klev, &
            paero, pcloud, pprecp, pice, psnow, &
            ptemp, ppres, prs, prsi, prv, ptstep)
        IF (sflg) CALL salsa_var_stat('cond',1)
    ENDIF

    ! Autoconversion (liquid)
    !   Statistics: change in total rain water volume (=change in cloud water) and rain drop number concentration
    IF (lsauto) THEN
         IF (sflg) CALL salsa_var_stat('auto',0)
         IF (auto_sb) THEN
            CALL autoconv_sb(kbdim,klev,ptstep,pcloud,pprecp)
         ELSE
            CALL autoconv2(kbdim,klev,pcloud, pprecp)
         ENDIF
         IF (sflg) CALL salsa_var_stat('auto',1)
    ENDIF

    ! Cloud activation
    !   Statistics: change in total cloud water volume (=change in cloud water) and cloud drop number concentration
    IF (lsactiv ) THEN
         IF (sflg) CALL salsa_var_stat('cact',0)
         CALL cloud_activation(kbdim,  klev,          &
                               ptemp, prv, prs, paero, pcloud)
         IF (sflg) CALL salsa_var_stat('cact',1)
    ENDIF

    ! Ice nucleation
    !   Statistics: change in total ice/snow* water volume and ice/snow number concentration
    !   * ice nucleation can also produce snow and in this case snow formation rate is saved
    !     to the autoconversion variables (no ice category; autoconversion disabled)
    IF (lsicenucl) THEN
      IF (sflg) CALL salsa_var_stat('nucl',0) ! Total
      IF (fixinc>=0. .OR. ice_diag<0) THEN
        ! Fixed (fixinc>0.0) or diagnostic (ice_diag<0) ice number concentration
        IF (sflg) CALL salsa_var_stat('nucf',0) ! Fixed ice
        CALL fixed_ice_driver(kbdim, klev,             &
                             pcloud, pice,   psnow,    &
                             ptemp,  ppres,  prv, prs, prsi)
        IF (sflg) CALL salsa_var_stat('nucf',1)
      ENDIF
     IF (ice_diag>0) THEN
        ! Various other ice formation parameterizations
        IF (sflg) CALL salsa_var_stat('nucf',0) ! Use the same variable as for fixed ice
        CALL param_ice_driver(kbdim, klev,             &
                             paero, pcloud, pprecp, pice, psnow, &
                             ice_diag, ptemp, ppres, prv, prs, prsi)
        IF (sflg) CALL salsa_var_stat('nucf',1)
      ENDIF
      IF (ice_hom .OR. ice_imm .OR. ice_dep) THEN
        ! Modelled ice nucleation
        IF (sflg) CALL salsa_var_stat('nucm',0) ! Modelled ice
        CALL ice_nucl_driver(kbdim,klev,   &
                          paero,pcloud,pprecp,pice,psnow, &
                          ptemp,prv,prs,prsi,ptstep)
        IF (sflg) CALL salsa_var_stat('nucm',1)
      ENDIF
      IF (sflg) CALL salsa_var_stat('nucl',1)
    ENDIF

    ! Melting of ice and snow
    !   Statistics: change in total ice and snow water volume and number concentrations
    IF (lsicmelt) THEN
         IF (sflg) CALL salsa_var_stat('melt',0)
         CALL ice_melt(kbdim,klev,pcloud,pice,pprecp,psnow,ptemp)
         IF (sflg) CALL salsa_var_stat('melt',1)
    ENDIF

    ! Snow formation ~ autoconversion from ice
    !   Statistics: change in total snow water volume (=change in ice water) and snow number concentration
    IF (lsautosnow) THEN
         IF (sflg) CALL salsa_var_stat('auto',0) ! Note: the same output name for warm cloud autoconversion!
         CALL autosnow(kbdim,klev,pice,psnow)
         IF (sflg) CALL salsa_var_stat('auto',1)
    ENDIF

    ! Size distribution bin update
    IF (lsdistupdate ) THEN
        IF (sflg) CALL salsa_var_stat('dist',0)
         CALL distr_update(kbdim, klev,             &
                           paero,  pcloud, pprecp,  &
                           pice, psnow, level       )
        IF (sflg) CALL salsa_var_stat('dist',1)
    ENDIF


    CONTAINS

      ! Functions for generating requested raw output data for the LES model. These functions
      ! use parameters (e.g. nstat, sdata, slist, paero, pcloud,...) from the calling program.

      ! 1) Master routine
      SUBROUTINE salsa_var_stat(prefix,ncall)
        IMPLICIT NONE
        ! Input
        character (len=4), intent (in) :: prefix ! Process name
        INTEGER, INTENT(IN) :: ncall ! Before or after process call
        ! Local
        INTEGER :: i
        CHARACTER(LEN=7) :: nam
        !
        ! Examine outputs
        DO i=1,nstat
            ! The first four characters contain process name, which should match with the given prefix
            nam=slist(i)
            IF (nam(1:4)==prefix) THEN
                ! The 5th character is just '_', but the 6th character indicates the output type
                ! (n=number, r=water mixing ratio or an integer)
                !
                ! The last (7th) character is phase (a, c, r, i, s or g)
                SELECT CASE (nam(7:7))
                    CASE ('a')
                        CALL salsa_rate_stat(i,nam(6:6),fn2b,paero,ncall)
                    CASE ('c')
                        CALL salsa_rate_stat(i,nam(6:6),ncld,pcloud,ncall)
                    CASE ('r')
                        CALL salsa_rate_stat(i,nam(6:6),nprc,pprecp,ncall)
                    CASE ('i')
                        CALL salsa_rate_stat(i,nam(6:6),nice,pice,ncall)
                    CASE ('s')
                        CALL salsa_rate_stat(i,nam(6:6),nsnw,psnow,ncall)
                    CASE ('g')
                        CALL salsa_rate_stat_gas(i,nam(6:6),ngas,pc_gas,ncall)
                END SELECT
            ENDIF
        ENDDO
        !
      END SUBROUTINE salsa_var_stat

      ! 2a) The actual subroutine for doing the work for t_section types
      SUBROUTINE salsa_rate_stat(i,tchar,nb,curr,ncall)
        USE mo_submctl, ONLY : rhowa
        IMPLICIT NONE
        ! Inputs
        INTEGER, intent (in) :: i ! the i:th output
        CHARACTER, intent (in) :: tchar ! output type character (n, r, or an integer number)
        INTEGER, INTENT(IN) :: nb ! number of bins
        TYPE(t_section), INTENT(in) :: curr(kbdim,klev,nb) ! current state of a hydrometeor
        INTEGER, INTENT(IN) :: ncall ! calling before (0) or after (1) the main function call
        ! Local
        INTEGER :: j
        !
        ! Check the output type
        SELECT CASE (tchar)
            CASE ('n')
                IF (ncall==0) THEN
                    ! Set: sum of number concentrations over bins (#/m3)
                    sdata(:,:,i) = SUM(curr(:,:,:)%numc,DIM=3)
                ELSE
                    ! Difference
                    sdata(:,:,i) = SUM(curr(:,:,:)%numc,DIM=3) - sdata(:,:,i)
                ENDIF
            CASE('r')
                ! Sum of water volume mixing ratio over bins (m3/m3)
                IF (ncall==0) THEN
                    ! Set: sum of water volume mixing ratios over bins (m3/m3)
                    sdata(:,:,i) = SUM(curr(:,:,:)%volc(1),DIM=3)
                ELSE
                    ! Difference multiplied by water density (kg/m3)
                    sdata(:,:,i) = (SUM(curr(:,:,:)%volc(1),DIM=3) - sdata(:,:,i))*rhowa
                ENDIF
            CASE DEFAULT
                ! Species-dependent mixing ratio outputs
                ! Convert to integer - used as is
                READ(UNIT=tchar,FMT='(I1)') j
                IF (ncall==0) THEN
                    ! Set: sum of component j volume mixing ratio over bins (m3/m3)
                    sdata(:,:,i) = SUM(curr(:,:,:)%volc(j),DIM=3)
                ELSE
                    ! Difference multiplied by density (kg/m3)
                    sdata(:,:,i) = (SUM(curr(:,:,:)%volc(j),DIM=3) - sdata(:,:,i))*rhowa
                ENDIF
        END SELECT
        !
      END SUBROUTINE salsa_rate_stat

      ! 2b) The same for gases (not a t_section type)
      SUBROUTINE salsa_rate_stat_gas(i,tchar,ng,curr,ncall)
        USE mo_submctl, ONLY : mws_gas
        ! The list of requested outputs
        IMPLICIT NONE
        ! Inputs
        INTEGER, intent (in) :: i !  the i:th output
        CHARACTER, intent (in) :: tchar ! output type character (an integer number)
        INTEGER, INTENT(IN) :: ng ! number of gases
        REAL, INTENT(in) :: curr(kbdim,klev,ng) ! current state of the parameter
        INTEGER, INTENT(IN) :: ncall ! calling before (0) or after (1) the main function call
        ! Local
        INTEGER :: j
        !
        ! Convert to integer - used as is
        READ(UNIT=tchar,FMT='(I1)') j
        IF (j==0) j=10 ! Allow tenth gas
        ! Gas phase component j molar mixing ratio (mol/m3)
        IF (ncall==0) THEN
            sdata(:,:,i) = curr(:,:,j)
        ELSE
            sdata(:,:,i) = (curr(:,:,j) - sdata(:,:,i))*mws_gas(j) ! Convert to kg/m3
        ENDIF
        !
      END SUBROUTINE salsa_rate_stat_gas

  END SUBROUTINE salsa


END MODULE mo_salsa
