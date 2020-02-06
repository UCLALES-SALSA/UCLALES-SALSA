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
                   ptemp,  ptstep, petime,              &
                   pc_gas, ngas,                        &
                   paero,  pcloud, pprecp, pice, psnow, &
                   level, sflg, nstat, sdata, slist,    &
                   coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, &
                   coag_nprecp, coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
                   cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
                   autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
                   act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
                   melt_vice, melt_nice, melt_vsnow, melt_nsnow)

    USE mo_salsa_dynamics, only : coagulation, condgas, gpparth2o
    USE mo_vbs_partition, ONLY : vbs_gas_phase_chem, vbs_condensation
    USE mo_salsa_update, ONLY : distr_update
    USE mo_salsa_cloud, only : cloud_activation, autoconv2, autoconv_sb, &
            autosnow, fixed_ice_driver, ice_nucl_driver, ice_melt

    USE mo_submctl, ONLY :      &
         fn2b,ncld,nprc,nice,nsnw,nvbs,    &
         t_section,part_h2so4,part_ocnv,   &
         lscoag,lscnd,lscndgas,nlcndh2oae, &
         nlcndh2ocl,nlcndh2oic,            &
         lsauto,auto_sb,lsautosnow,lsactiv,&
         lsicenucl,lsicmelt,lsdistupdate,  &
         fixinc, ice_hom, ice_imm, ice_dep

    IMPLICIT NONE

    !-- Input parameters and variables --------------------------------------------------
    INTEGER, INTENT(in) ::          &
         kbdim,                     & ! dimension for arrays (='kbdim')
         klev,                      & ! number of vertical levels (='klev')
         ngas                         ! number of gases


    REAL, INTENT(in) ::            &
         ppres(kbdim,klev),            & ! atmospheric pressure at each grid point [Pa]
         ptemp(kbdim,klev),            & ! temperature at each grid point [K]
         ptstep,                       & ! time step [s]
         petime                          ! elapsed time [s]

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

    !-- Output statistics --------------------------------------
    REAL, DIMENSION(kbdim,klev), INTENT(out) :: &
            coag_vaero, coag_naero, coag_vcloud, coag_ncloud, coag_vprecp, coag_nprecp, &
            coag_vice, coag_nice, coag_vsnow, coag_nsnow, &
            cond_vaero, cond_vcloud, cond_vprecp, cond_vice, cond_vsnow, &
            autoc_vprecp, autoc_nprecp, autoc_vsnow, autoc_nsnow, &
            act_vcloud, act_ncloud, nucl_vice, nucl_nice, &
            melt_vice, melt_nice, melt_vsnow, melt_nsnow

    ! Reset outputs
    coag_vaero=0.; coag_naero=0.; coag_vcloud=0.; coag_ncloud=0.; coag_vprecp=0.; coag_nprecp=0.
    coag_vice=0.; coag_nice=0.; coag_vsnow=0.; coag_nsnow=0.
    cond_vaero=0.; cond_vcloud=0.; cond_vprecp=0.; cond_vice=0.; cond_vsnow=0.
    autoc_vprecp=0.; autoc_nprecp=0.; autoc_vsnow=0.; autoc_nsnow=0.
    act_vcloud=0.; act_ncloud=0.; nucl_vice=0.; nucl_nice=0.;
    melt_vice=0.; melt_nice=0.; melt_vsnow=0.; melt_nsnow=0.

    sdata(:,:,:) = 0.

    ! Coagulation
    !   Statistics: change in total water volume concentration for each hydrometeor (m^3/m^3)
    IF (lscoag) THEN
       coag_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)
       coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)
       coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
       coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
       coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
       coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
       coag_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
       coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
       coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
       coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)

       IF (sflg) CALL salsa_var_stat('coag',0)
       CALL coagulation(kbdim,  klev,                           &
                        paero,  pcloud, pprecp, pice, psnow,    &
                        ptstep, ptemp,  ppres                   )
       IF (sflg) CALL salsa_var_stat('coag',1)

       coag_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)-coag_vaero(:,:)
       coag_naero(:,:)=SUM(paero(:,:,:)%numc,DIM=3)-coag_naero(:,:)
       coag_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-coag_vcloud(:,:)
       coag_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-coag_ncloud(:,:)
       coag_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-coag_vprecp(:,:)
       coag_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-coag_nprecp(:,:)
       coag_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-coag_vice(:,:)
       coag_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-coag_nice(:,:)
       coag_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-coag_vsnow(:,:)
       coag_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-coag_nsnow(:,:)
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
            ptemp, ptstep, petime, pc_gas, ngas)
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
        cond_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)
        cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
        cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
        cond_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        IF (sflg) CALL salsa_var_stat('cond',0)
        CALL gpparth2o(kbdim, klev, &
            paero, pcloud, pprecp, pice, psnow, &
            ptemp, ppres, prs, prsi, prv, ptstep)
        IF (sflg) CALL salsa_var_stat('cond',1)
        cond_vaero(:,:)=SUM(paero(:,:,:)%volc(1),DIM=3)-cond_vaero(:,:)
        cond_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-cond_vcloud(:,:)
        cond_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-cond_vprecp(:,:)
        cond_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-cond_vice(:,:)
        cond_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-cond_vsnow(:,:)
    ENDIF

    ! Autoconversion (liquid)
    !   Statistics: change in total rain water volume (=change in cloud water) and rain drop number concentration
    IF (lsauto) THEN
         autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)
         autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)
         IF (sflg) CALL salsa_var_stat('auto',0)
         IF (auto_sb) THEN
            CALL autoconv_sb(kbdim,klev,ptstep,pcloud,pprecp)
         ELSE
            CALL autoconv2(kbdim,klev,pcloud, pprecp)
         ENDIF
         IF (sflg) CALL salsa_var_stat('auto',1)
         autoc_vprecp(:,:)=SUM(pprecp(:,:,:)%volc(1),DIM=3)-autoc_vprecp(:,:)
         autoc_nprecp(:,:)=SUM(pprecp(:,:,:)%numc,DIM=3)-autoc_nprecp(:,:)
    ENDIF

    ! Cloud activation
    !   Statistics: change in total cloud water volume (=change in cloud water) and cloud drop number concentration
    IF (lsactiv ) THEN
         act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)
         act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)
         IF (sflg) CALL salsa_var_stat('cact',0)
         CALL cloud_activation(kbdim,  klev,          &
                               ptemp,  ppres, prv,    &
                               prs,    paero, pcloud  )
         IF (sflg) CALL salsa_var_stat('cact',1)
         act_vcloud(:,:)=SUM(pcloud(:,:,:)%volc(1),DIM=3)-act_vcloud(:,:)
         act_ncloud(:,:)=SUM(pcloud(:,:,:)%numc,DIM=3)-act_ncloud(:,:)
    ENDIF

    ! Ice nucleation
    !   Statistics: change in total ice/snow* water volume and ice/snow number concentration
    !   * ice nucleation can also produce snow and in this case snow formation rate is saved
    !     to the autoconversion variables (no ice category; autoconversion disabled)
    IF (lsicenucl .AND. fixinc>=0.) THEN
        ! Fixed ice number concentration
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
        IF (sflg) CALL salsa_var_stat('nucl',0)
        CALL fixed_ice_driver(kbdim, klev,             &
                             pcloud, pice,   psnow,    &
                             ptemp,  ppres,  prv,  prsi)
        IF (sflg) CALL salsa_var_stat('nucl',1)
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-nucl_vice(:,:)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-nucl_nice(:,:)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
    ELSEIF (lsicenucl .AND. (ice_hom .OR. ice_imm .OR. ice_dep)) THEN
        ! Modelled ice nucleation
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
        IF (sflg) CALL salsa_var_stat('nucl',0)
        CALL ice_nucl_driver(kbdim,klev,   &
                          paero,pcloud,pprecp,pice,psnow, &
                          ptemp,prv,prs,prsi,ptstep)
        IF (sflg) CALL salsa_var_stat('nucl',1)
        nucl_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-nucl_vice(:,:)
        nucl_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-nucl_nice(:,:)
        autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
        autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
    ENDIF

    ! Melting of ice and snow
    !   Statistics: change in total ice and snow water volume and number concentrations
    IF (lsicmelt) THEN
         melt_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)
         melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)
         melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
         melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
         IF (sflg) CALL salsa_var_stat('melt',0)
         CALL ice_melt(kbdim,klev,pcloud,pice,pprecp,psnow,ptemp)
         IF (sflg) CALL salsa_var_stat('melt',1)
         melt_vice(:,:)=SUM(pice(:,:,:)%volc(1),DIM=3)-melt_vice(:,:)
         melt_nice(:,:)=SUM(pice(:,:,:)%numc,DIM=3)-melt_nice(:,:)
         melt_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-melt_vsnow(:,:)
         melt_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-melt_nsnow(:,:)
    ENDIF

    ! Snow formation ~ autoconversion from ice
    !   Statistics: change in total snow water volume (=change in ice water) and snow number concentration
    IF (lsautosnow) THEN
         autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)
         autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)
         IF (sflg) CALL salsa_var_stat('auto',0) ! Note: the same output name for warm cloud autoconversion!
         CALL autosnow(kbdim,klev,pice,psnow)
         IF (sflg) CALL salsa_var_stat('auto',1)
         autoc_vsnow(:,:)=SUM(psnow(:,:,:)%volc(1),DIM=3)-autoc_vsnow(:,:)
         autoc_nsnow(:,:)=SUM(psnow(:,:,:)%numc,DIM=3)-autoc_nsnow(:,:)
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
        USE mo_submctl, ONLY : rhowa, rhoic, rhosn
        IMPLICIT NONE
        ! Input
        character (len=4), intent (in) :: prefix ! Process name
        INTEGER, INTENT(IN) :: ncall ! Before or after process call
        !
        ! Level 4
        IF (level>=4) THEN
            CALL salsa_rate_stat(prefix,'a',fn2b,paero,ncall,rhowa)
            CALL salsa_rate_stat(prefix,'c',ncld,pcloud,ncall,rhowa)
            CALL salsa_rate_stat(prefix,'r',nprc,pprecp,ncall,rhowa)
        ENDIF
        ! Level 5
        IF (level>=5) THEN
            CALL salsa_rate_stat(prefix,'i',nice,pice,ncall,rhoic)
            CALL salsa_rate_stat(prefix,'s',nsnw,psnow,ncall,rhosn)
        ENDIF
        ! Gases
        IF (ngas>0) CALL salsa_rate_stat_gas(prefix,'g',ngas,pc_gas,ncall)
        !
      END SUBROUTINE salsa_var_stat

      ! 2a) The actual subroutine for doing the work for t_section types
      SUBROUTINE salsa_rate_stat(prefix,tchar,nb,curr,ncall,rhowa)
        USE mo_submctl, ONLY : nspec, dens
        IMPLICIT NONE
        ! Inputs
        character (len=4), intent (in) :: prefix ! variable name related to a process, e.g. 'coag'
        CHARACTER, intent (in) :: tchar ! character indicating the type of the input hydrometeor (a, c, r, i, s)
        INTEGER, INTENT(IN) :: nb ! number of bins
        TYPE(t_section), INTENT(in) :: curr(kbdim,klev,nb) ! current state of a hydrometeor
        INTEGER, INTENT(IN) :: ncall ! calling before (0) or after (1) the main function call
        REAL, INTENT(IN) :: rhowa ! Density of water (liquid, ice or snow)
        ! Local
        INTEGER :: i, j
        CHARACTER(LEN=7) :: nam
        !
        ! Is this output selected
        DO i=1,nstat
            IF ( prefix//'_n'//tchar == slist(i) ) THEN
                IF (ncall==0) THEN
                    ! Set: sum of number concentrations over bins (#/m3)
                    sdata(:,:,i) = SUM(curr(:,:,:)%numc,DIM=3)
                ELSE
                    ! Difference
                    sdata(:,:,i) = SUM(curr(:,:,:)%numc,DIM=3) - sdata(:,:,i)
                ENDIF
            ENDIF
            IF ( prefix//'_r'//tchar == slist(i) ) THEN
                ! Sum of water volume mixing ratio over bins (m3/m3)
                IF (ncall==0) THEN
                    ! Set: sum of water volume mixing ratios over bins (m3/m3)
                    sdata(:,:,i) = SUM(curr(:,:,:)%volc(1),DIM=3)
                ELSE
                    ! Difference multiplied by water density (kg/m3)
                    sdata(:,:,i) = (SUM(curr(:,:,:)%volc(1),DIM=3) - sdata(:,:,i))*rhowa
                ENDIF
            ENDIF
            ! Species and phase-dependent mixing ratio outputs
            DO j=1,nspec+1
                WRITE(nam,'(A4,A1,I1,A1)') prefix,'_',j,tchar ! e.g. 'cond_1a'
                IF ( nam == slist(i) ) THEN
                    IF (ncall==0) THEN
                        ! Set: sum of component j volume mixing ratio over bins (m3/m3)
                        sdata(:,:,i) = SUM(curr(:,:,:)%volc(j),DIM=3)
                    ELSE
                        ! Difference multiplied by density (kg/m3)
                        sdata(:,:,i) = (SUM(curr(:,:,:)%volc(j),DIM=3) - sdata(:,:,i))*dens(j)
                        IF (j==1) sdata(:,:,i)=sdata(:,:,i)*rhowa/dens(j)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
      END SUBROUTINE salsa_rate_stat

      ! 2b) The same for gases
      SUBROUTINE salsa_rate_stat_gas(prefix,tchar,ng,curr,ncall)
        ! The list of requested outputs
        IMPLICIT NONE
        ! Inputs
        character (len=4), intent (in) :: prefix ! variable name related to a process, e.g. 'coag'
        CHARACTER, intent (in) :: tchar ! character indicating the type of the input hydrometeor (here just g)
        INTEGER, INTENT(IN) :: ng ! number of gases
        REAL, INTENT(in) :: curr(kbdim,klev,ng) ! current state of the parameter
        INTEGER, INTENT(IN) :: ncall ! calling before (0) or after (1) the main function call
        ! Local
        INTEGER :: i, j
        CHARACTER(LEN=7) :: nam
        !
        ! Is this output selected
        DO i=1,nstat
            ! Species dependent
            DO j=1,ng
                WRITE(nam,'(A4,A1,I1,A1)') prefix,'_',j,tchar ! e.g. 'cond_1g'
                IF ( nam == slist(i) ) THEN
                    ! Gas phase component j molar mixing ratio (mol/m3)
                    IF (ncall==0) THEN
                        sdata(:,:,i) = curr(:,:,j)
                    ELSE
                        sdata(:,:,i) = curr(:,:,j) - sdata(:,:,i)
                    ENDIF
                ENDIF
            ENDDO
        ENDDO
      END SUBROUTINE salsa_rate_stat_gas

  END SUBROUTINE salsa


END MODULE mo_salsa
