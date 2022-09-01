!****************************************************************
! Module mo_vbs_init: initialization of the VOC/VBS part for UCLALES-SALSA
!
! To be used with UCLALES-SALSA only. Initializes local (i.e. within the VBS code)
! gas and aerosol species based on options given in the SALSA namelist. Exports
! properties of new species (VOCs and VBS species) into UCLALES-SALSA.
!
! Tomi Raatikainen (FMI)
!
!****************************************************************
!
MODULE mo_vbs_init

  IMPLICIT NONE

  CONTAINS

SUBROUTINE init_vbs(nvbs_setup, laqsoa)
    USE mo_vbs, ONLY : vbs_species, spec_moleweight, spec_density, spec_kappa
    USE mo_vbsctl, ONLY : vbs_voc_set, vbs_set, aqsoa_set,  &
            vbs_nvocs, vbs_ngroup, aqsoa_ngroup
    USE mo_submctl, ONLY : dens, mws, diss, nspec, zspec, moc, dissoc, rhooc, &
            nvocs, nvbs, naqsoa, mws_gas, ngases, zgas, &
            id_oh, id_no3, id_o3, conc_oh, conc_no3, conc_o3, ox_prescribed, ngases_diag
    ! Inputs
    LOGICAL, INTENT(IN) :: laqsoa
    INTEGER, INTENT(IN) :: nvbs_setup
    ! Local
    INTEGER :: i, j
    CHARACTER(LEN=3) :: tmp
    INTEGER :: subm_naerospec=0, subm_ngasspec=0
    REAL :: kappaoc

    ! 1) VBS initialization
    ! *********************
    ! VBS setup: define VOC, VBS and aqSOA sets and their physical properties
    kappaoc=dissoc*( mws(1)/dens(1) )/( moc/rhooc )
    CALL vbs_species(nvbs_setup, laqsoa, rhooc, moc, kappaoc)

    ! Set indexes: SALSA/LES species always before VBS species
    subm_naerospec = nspec+1 ! Water and other aerosol species
    subm_ngasspec = ngases ! SALSA gases and vapors like SO2 and oxidants
    ! a) VOCs - just gas phase
    DO i = 1, vbs_nvocs
        subm_ngasspec = subm_ngasspec + 1
        vbs_voc_set(i)%id_gas = subm_ngasspec ! Index to gas concentrations
    END DO
    ! b) VBS bins - both gas and particle phases
    DO i = 1, vbs_ngroup
        subm_ngasspec = subm_ngasspec + 1
        vbs_set(i)%id_gas = subm_ngasspec
        subm_naerospec = subm_naerospec + 1
        vbs_set(i)%id_vols = subm_naerospec
    END DO
    ! c) aqSOA - gas and particle phases
    DO i = 1, aqsoa_ngroup
        subm_ngasspec = subm_ngasspec + 1
        aqsoa_set(i)%id_gas = subm_ngasspec
        subm_naerospec = subm_naerospec + 1
        aqsoa_set(i)%id_vols = subm_naerospec
    END DO


    ! 2) Outputs
    ! **********
    ! From VBS to SALSA
    nvocs = vbs_nvocs ! Number of VOCs
    nvbs = vbs_ngroup ! Number of VBS bins
    naqsoa = aqsoa_ngroup ! Number of aqSOA species

    ! 1) Gases/vapors: count the total number and save molecular weights for SALSA/LES
    ! a) Oxidants (OH, O3 and NO) - either prognostic or diagnostic (constant concentration)
    IF (conc_oh>=0.) THEN
        IF (ox_prescribed) THEN
            ngases_diag = ngases_diag + 1 ! Diagnostic
            id_oh = vbs_nvocs + vbs_ngroup + aqsoa_ngroup + ngases_diag ! Append after the prognostic variables
        ELSE
            ngases = ngases + 1 ! Prognostic
            id_oh = ngases ! Index to gas phase
        ENDIF
        zgas(id_oh) = 'OH' ! Name
        mws_gas(id_oh) = 17e-3 ! OH(g)
    ENDIF
    IF (conc_o3>=0.) THEN
        IF (ox_prescribed) THEN
            ngases_diag = ngases_diag + 1
            id_o3 = vbs_nvocs + vbs_ngroup + aqsoa_ngroup + ngases_diag
        ELSE
            ngases = ngases + 1
            id_o3 = ngases
        ENDIF
        zgas(id_o3) = 'O3'
        mws_gas(id_o3) = 48e-3 ! O3(g)
    ENDIF
    IF (conc_no3>=0.) THEN
        IF (ox_prescribed) THEN
            ngases_diag = ngases_diag + 1
            id_no3 = vbs_nvocs + vbs_ngroup + aqsoa_ngroup + ngases_diag
        ELSE
            ngases = ngases + 1
            id_no3 = ngases
        ENDIF
        zgas(id_no3) = 'NO'
        mws_gas(id_no3) = 52e-3 ! NO3(g)
    ENDIF
    ! b) VOCs
    DO i=1, vbs_nvocs
        ngases = ngases + 1
        mws_gas(ngases)=spec_moleweight(vbs_voc_set(i)%spid)*1e-3
        WRITE(tmp,"('VO',I1)") i
        zgas(ngases)=tmp
    ENDDO
    ! c) VBS species
    DO i=1, vbs_ngroup
        ngases = ngases + 1
        mws_gas(ngases)=spec_moleweight(vbs_set(i)%spid)*1e-3
        WRITE(tmp,"('VB',I1)") i
        zgas(ngases)=tmp
    ENDDO
    ! d) aqSOA species
    DO i=1, aqsoa_ngroup
        ngases = ngases + 1
        mws_gas(ngases)=spec_moleweight(aqsoa_set(i)%spid)*1e-3
        WRITE(tmp,"('AQ',I1)") i
        zgas(ngases)=tmp
    ENDDO

    ! Add VBS and aqSOA to SALSA aerosol species list
    ! a) VBS species
    DO i=1,vbs_ngroup
        nspec=nspec + 1
        j=nspec+1 ! Water is the first, species but it is not included in nspec
        dens(j)=spec_density(vbs_set(i)%spid)
        mws(j)=spec_moleweight(vbs_set(i)%spid)*1e-3
        ! kappa-Kohler: k=diss_i*v_w/v_i => diss_i=k*v_i/v_w, where v=M/rho
        diss(j)=spec_kappa(vbs_set(i)%spid)*( mws(j)/dens(j) )/( mws(1)/dens(1) )
        ! SALSA/LES name
        WRITE(tmp,"('VB',I1)") i
        zspec(j)=tmp
    ENDDO
    ! b) aqSOA
    DO i=1,aqsoa_ngroup
        nspec=nspec + 1
        j=nspec+1
        dens(j)=spec_density(aqsoa_set(i)%spid)
        mws(j)=spec_moleweight(aqsoa_set(i)%spid)*1e-3
        diss(j)=spec_kappa(aqsoa_set(i)%spid)*( mws(j)/dens(j) )/( mws(1)/dens(1) )
        WRITE(tmp,"('AQ',I1)") i
        zspec(j)=tmp
    ENDDO

END SUBROUTINE init_vbs

END MODULE mo_vbs_init
