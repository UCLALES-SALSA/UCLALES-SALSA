!****************************************************************
! Module mo_vbs: defines parameters for the volatile and semi-volatile compounds
!
! To be used with UCLALES-SALSA only. Subroutine mo_vbs_species is used to
! define all vapor phase parameters for the specified setup.
!
! Code developers
!  Declan O'Donnell (MPI-Met) - soa code; used as guideline here
!  Zhang (MPI-Met) - modifications to soa code
!  Thomas Kuehn (UEF) - original code (2015)
!  Joonas Merikanto (FMI) - added aqueous phase SOA
!  Tomi Raatikainen (FMI) - modified for UCLALES-SALSA (2019)
!
!****************************************************************

MODULE mo_vbs

  IMPLICIT NONE
  PRIVATE

  ! MW, density and kappa arrays
  REAL, PUBLIC, SAVE :: spec_moleweight(100), spec_density(100), spec_kappa(100)

  ! Initialize species
  PUBLIC :: vbs_species

CONTAINS

  SUBROUTINE vbs_species(nvbs_setup,laqsoa,densoc,mwoc,kappaoc)

     USE mo_vbsctl, ONLY: &
       t_voc_prec,         &
       t_vbs_group,        &
       t_aq_soa,           &
       vbs_nvocs,          &
       vbs_voc_set,        &
       vbs_ngroup,         &
       vbs_set,            &
       aqsoa_ngroup,       &
       aqsoa_set

    ! -----------------------------------------------------------------------
    !
    ! SUBROUTINE vbs_species
    !
    ! depending on nvbs_setup:
    ! - creates the species for the VOC, VBS and aqSOA species
    ! - sets up the control structures vbs_voc_set and vbs_set
    !
    ! -----------------------------------------------------------------------

    INTEGER, INTENT(IN) :: nvbs_setup
    LOGICAL, INTENT(IN) :: laqsoa
    REAL, INTENT(IN) :: densoc, mwoc, kappaoc

    REAL, PARAMETER :: argas = 8.314472 ! [J/K/mol] molar/universal/ideal gas constant

    INTEGER :: jv
    INTEGER :: spid_temp

    ! -----------------------------------------------------------------------
    ! executable procedure
    ! -----------------------------------------------------------------------


    ! the additional aqueous phase soa species, if laqsoa is on
    IF (laqsoa) THEN
       aqsoa_ngroup = 2
    ELSE
       aqsoa_ngroup = 0
    ENDIF

    SELECT CASE(nvbs_setup)

    CASE(1) ! 3 class scheme with volatilities 0, 1, and 10 ug/m3

       ! Defining the VOC precursor species
       ! we use five VOC precursors, Monoterpenes and Isoprene from Megan,
       ! and Toluene, Xylene, and Benzene as anthropogenics
       vbs_nvocs = 5

       ! the VBS has three volatility groups
       vbs_ngroup = 3

       ! allocating memory for the precursor properties:
       IF (.NOT. ALLOCATED(vbs_voc_set)) THEN
          ALLOCATE(vbs_voc_set(vbs_nvocs))!
          DO jv=1,vbs_nvocs
             ! allocating memory for the stoichiometric coefficients
             IF (.NOT. ALLOCATED(vbs_voc_set(jv)%stoich_coeff_oh)) THEN
                ALLOCATE(vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+aqsoa_ngroup), &
                         vbs_voc_set(jv)%stoich_coeff_o3(vbs_ngroup+aqsoa_ngroup), &
                         vbs_voc_set(jv)%stoich_coeff_no3(vbs_ngroup+aqsoa_ngroup))
             END IF
          END DO
       END IF

       ! allocating memory for the basis set:
       IF (.NOT. ALLOCATED(vbs_set)) THEN
          ALLOCATE(vbs_set(vbs_ngroup))
       END IF

       ! Defining the VOC species

       ! Monoterpenes (APIN)
       CALL new_species(&
            mw          = 136.,                  &
            idx         = vbs_voc_set(1)%spid &
            )
       !oxidation rates and their temperature dependence (Arrhenius)
       ! OH (data from IUPAC)
       vbs_voc_set(1)%k_0_OH     = 1.2E-11 ! pre-factor [m3/(mol*s)]
       vbs_voc_set(1)%Eact_p_OH  = 440.    ! reduced activation energy [K]: Eact_p=Eact/R
       ! O3 (data from IUPAC)
       vbs_voc_set(1)%k_0_O3     = 6.3E-16 ! pre-factor [m3/(mol*s)]
       vbs_voc_set(1)%Eact_p_O3  = -580.   ! reduced activation energy [K]: Eact_p=Eact/R
       ! NO3
       vbs_voc_set(1)%k_0_NO3    = 1.2E-12 !pre-factor [m3/(mol*s)]
       vbs_voc_set(1)%Eact_p_NO3 = 490. !reduced activation energy [K]: Eact_p=Eact/R

       ! Stoichiometric Coefficients (different amount depending on laqsoa)
       IF (laqsoa) THEN
          ! from Harri (hi NOx)
          !                            (  OC   ,  VBS1   ,  VBS10  , IEPOX , Glyx)
          vbs_voc_set(1)%stoich_coeff_oh=(/0.1, 0.037, 0.088, 0.0, 0.0/)
          vbs_voc_set(1)%stoich_coeff_o3=(/0.1, 0.037, 0.088, 0.0, 0.0/)
          vbs_voc_set(1)%stoich_coeff_no3=(/0.1, 0.037, 0.088, 0.0, 0.0/)
       ELSE
          !                            (  OC      ,  VBS1   , VBS10    )
          vbs_voc_set(1)%stoich_coeff_oh=(/0.1, 0.037, 0.088/) ! from Harri (hi NOx)
          !vbs_voc_set(1)%stoich_coeff_oh=(/0.002, 0.003, 0.065/) ! from Harri (low NOx)
          vbs_voc_set(1)%stoich_coeff_o3=(/0.1, 0.037, 0.088/)
          vbs_voc_set(1)%stoich_coeff_no3=(/0.1, 0.037, 0.088/)
       ENDIF

       ! Isoprene (C5H8)
       CALL new_species(&
            mw          = 68.,                   &
            idx         = vbs_voc_set(2)%spid &
          )
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_set(2)%k_0_OH     = 2.7E-11 !pre-factor [m3/(mol*s)]
       vbs_voc_set(2)%Eact_p_OH  = 390. !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_set(2)%k_0_O3     = 1.03E-14 !pre-factor [m3/(mol*s)]
       vbs_voc_set(2)%Eact_p_O3  = -1995. !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_set(2)%k_0_NO3    = 3.15E-12 !pre-factor [m3/(mol*s)]
       vbs_voc_set(2)%Eact_p_NO3 = -450. !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC   ,  VBS1    ,   VBS10  ,   IEPOX , Glyx)
          vbs_voc_set(2)%stoich_coeff_oh=(/0.0, 0.0295, 0.0453, 0.525, 0.025/)
          vbs_voc_set(2)%stoich_coeff_o3=(/0.0, 0.0295, 0.0453, 0.0, 0.0/)
          vbs_voc_set(2)%stoich_coeff_no3=(/0.0, 0.0295, 0.0453, 0.0, 0.0/)
       ELSE
          !                            (  OC   ,  VBS1    , VBS10     )
          vbs_voc_set(2)%stoich_coeff_oh=(/0.0, 0.0295, 0.0453/)
          vbs_voc_set(2)%stoich_coeff_o3=(/0.0, 0.0295, 0.0453/)
          vbs_voc_set(2)%stoich_coeff_no3=(/0.0, 0.0295, 0.0453/)
       ENDIF

       ! Toluene (TOL)
       CALL new_species(&
            mw          = 92.,                   &
            idx         = vbs_voc_set(3)%spid &
            )
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_set(3)%k_0_OH     = 1.81E-12 !pre-factor [m3/(mol*s)]
       vbs_voc_set(3)%Eact_p_OH  = 338. !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_set(3)%k_0_O3     = 0.0 !pre-factor [m3/(mol*s)]
       vbs_voc_set(3)%Eact_p_O3  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_set(3)%k_0_NO3    = 0.0 !pre-factor [m3/(mol*s)]
       vbs_voc_set(3)%Eact_p_NO3 = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX , Glyx    )
          vbs_voc_set(3)%stoich_coeff_oh=(/0.36, 0.0, 0.0, 0.0, 0.24/) ! to be reviewed
          vbs_voc_set(3)%stoich_coeff_o3=(/0.36, 0.0, 0.0, 0.0, 0.0/)
          vbs_voc_set(3)%stoich_coeff_no3=(/0.36, 0.0, 0.0, 0.0, 0.0/)
       ELSE
          !                            (  OC    ,  VBS1 ,   VBS10)
          vbs_voc_set(3)%stoich_coeff_oh=(/0.36, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(3)%stoich_coeff_o3=(/0.36, 0.0, 0.0/)
          vbs_voc_set(3)%stoich_coeff_no3=(/0.36, 0.0, 0.0/)
       ENDIF

       ! Xylene (XYL)
       CALL new_species(&
            mw          = 106.,                  &
            idx         = vbs_voc_set(4)%spid &
            )
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_set(4)%k_0_OH     = 2.31E-11 !pre-factor [m3/(mol*s)]
       vbs_voc_set(4)%Eact_p_OH  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_set(4)%k_0_O3     = 0.0 !pre-factor [m3/(mol*s)]
       vbs_voc_set(4)%Eact_p_O3  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_set(4)%k_0_NO3    = 2.6E-16 !pre-factor [m3/(mol*s)]
       vbs_voc_set(4)%Eact_p_NO3 = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 ,  IEPOX,  Glyx   )
          vbs_voc_set(4)%stoich_coeff_oh=(/0.30, 0.0, 0.0, 0.0, 0.25/) ! to be reviewed
          vbs_voc_set(4)%stoich_coeff_o3=(/0.30, 0.0, 0.0, 0.0, 0.0/)
          vbs_voc_set(4)%stoich_coeff_no3=(/0.30, 0.0, 0.0, 0.0, 0.0/)
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_set(4)%stoich_coeff_oh=(/0.30, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(4)%stoich_coeff_o3=(/0.30, 0.0, 0.0/)
          vbs_voc_set(4)%stoich_coeff_no3=(/0.30, 0.0, 0.0/)
       ENDIF

       ! Benzene (BENZ)
       CALL new_species(&
            mw          = 66.,                   &
            idx         = vbs_voc_set(5)%spid &
            )
       !oxidation rates and their temperature dependence (Arrhenius)
       !OH
       vbs_voc_set(5)%k_0_OH     = 2.33E-12 !pre-factor [m3/(mol*s)]
       vbs_voc_set(5)%Eact_p_OH  = -193. !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       vbs_voc_set(5)%k_0_O3     = 0.0 !pre-factor [m3/(mol*s)]
       vbs_voc_set(5)%Eact_p_O3  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       vbs_voc_set(5)%k_0_NO3    = 0.0 !pre-factor [m3/(mol*s)]
       vbs_voc_set(5)%Eact_p_NO3 = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R

       ! Stoichiometric Coefficients
       IF (laqsoa) THEN
          !                            (  OC    ,  VBS1 , VBS10 , IEPOX ,  Glyx   )
          vbs_voc_set(5)%stoich_coeff_oh=(/0.37, 0.0, 0.0, 0.0, 0.35/) ! to be reviewed
          vbs_voc_set(5)%stoich_coeff_o3=(/0.37, 0.0, 0.0, 0.0, 0.0/)
          vbs_voc_set(5)%stoich_coeff_no3=(/0.37, 0.0, 0.0, 0.0, 0.0/)
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_set(5)%stoich_coeff_oh=(/0.37, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(5)%stoich_coeff_o3=(/0.37, 0.0, 0.0/)
          vbs_voc_set(5)%stoich_coeff_no3=(/0.37, 0.0, 0.0/)
       ENDIF

       ! Volatility Basis set C*=0 ug/m3
       CALL new_species(&
            mw          = 186.,                      & ! to be reviewed
            density     = 1320.,                     & ! to be reviewed
            kappa       = 0.037,                     & ! Petters and Kreidenweis (2007)
            idx         = spid_temp                     &
            )

       vbs_set(1)%spid        = spid_temp
       vbs_set(1)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            0.0e-6/spec_moleweight(spid_temp)     ! [M]=g/mol
       vbs_set(1)%T0          = 298               ! T0 [K]
       vbs_set(1)%Hvap_eff    = 30e3/argas        ! eff. evap. enthalpy [K]

       ! Volatility Basis set C*=1 ug/m3
       CALL new_species(&
            mw          = 186.,                      & ! to be reviewed
            density     = 1320.,                     & ! to be reviewed
            kappa       = 0.037,                     & ! Petters and Kreidenweis (2007)
            idx         = spid_temp                     &
            )

       vbs_set(2)%spid        = spid_temp
       vbs_set(2)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            1.0e-6/spec_moleweight(spid_temp)      ! [M]=g/mol
       vbs_set(2)%T0          = 298.               ! T0 [K]
       vbs_set(2)%Hvap_eff    = 30e3/argas        ! eff. evap. enthalpy [K]

       ! Volatility Basis set C*=10 ug/m3
       CALL new_species(&
            mw          = 186.,                      & ! to be reviewed
            density     = 1320.,                     & ! to be reviewed
            kappa       = 0.037,                     & !Petters and Kreidenweis (2007)
            idx         = spid_temp                     &
            )

       vbs_set(3)%spid        = spid_temp
       vbs_set(3)%C0          = &                    ! equ. vapor con. at T0 [mol/m3]
            10.0e-6/spec_moleweight(spid_temp)     ! [M]=g/mol
       vbs_set(3)%T0          = 298.               ! T0 [K]
       vbs_set(3)%Hvap_eff    = 30e3/argas        ! eff. evap. enthalpy [K]

    CASE(2) ! VBS-only schemes
       ! predetermined volatility bins
       ! physical properties (density, mw, kappa) taken from SALSA organics

       ! VBS scheme
       vbs_ngroup = 5 ! Four semivolatile bins and one non-volatile bin
       vbs_nvocs = 0  ! No VOCs or aqSOA

       ! Defining VBS species
       IF (.NOT. ALLOCATED(vbs_set)) ALLOCATE(vbs_set(vbs_ngroup))
       ! equ. vapor conc. at T0 (ug=1e-9 kg converted to mol)
       vbs_set(:)%C0 = (/0.0,1e-5,1e-3,1e-1,1e1/)*1e-9/mwoc ! C0 [mol/m3]
       vbs_set(:)%T0 = 298. ! T0 [K]
       vbs_set(:)%Hvap_eff = 30e3/argas ! eff. evap. enthalpy [K]
       DO jv=1, vbs_ngroup
            ! Physical properties from SALSA (mw from kg/mol to g/mol)
            CALL new_species(mw=mwoc*1e3, density=densoc, kappa=kappaoc, idx=vbs_set(jv)%spid)
       ENDDO

    CASE(3) ! VBS scheme with VOCs
       ! volatility bins based on VOC oxidation (the scheme from Farina et al., JGR, 2010)
       ! physical properties (density, mw, kappa) taken from SALSA organics

       ! VBS scheme
       vbs_ngroup = 5 ! One non-volatile bin and four semivolatile bins from Farina et al. (2010)
       vbs_nvocs = 2  ! Two VOCs (ALPH and ISOP), and aqSOA is possible

       ! Defining the VOC species
       IF (.NOT. ALLOCATED(vbs_voc_set)) THEN
            ALLOCATE(vbs_voc_set(vbs_nvocs))
            DO jv=1,vbs_nvocs
                ALLOCATE(vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+aqsoa_ngroup), &
                         vbs_voc_set(jv)%stoich_coeff_o3(vbs_ngroup+aqsoa_ngroup), &
                         vbs_voc_set(jv)%stoich_coeff_no3(vbs_ngroup+aqsoa_ngroup))
                vbs_voc_set(jv)%stoich_coeff_oh(:)=0.
                vbs_voc_set(jv)%stoich_coeff_o3(:)=0.
                vbs_voc_set(jv)%stoich_coeff_no3(:)=0.
            END DO
       END IF
       ! ALPH (MTP)
       ! Properties from Mielonen et al. (Atmosphere, 2018)
       CALL new_species(mw=136., idx=vbs_voc_set(1)%spid)
       vbs_voc_set(1)%k_0_OH     = 1.2E-11
       vbs_voc_set(1)%Eact_p_OH  = 440.
       vbs_voc_set(1)%k_0_O3     = 6.3E-16
       vbs_voc_set(1)%Eact_p_O3  = -580.
       vbs_voc_set(1)%k_0_NO3    = 1.2E-12
       vbs_voc_set(1)%Eact_p_NO3 = 490.
       ! Stoichiometric coefficients (Farina et al., Table 3):
       vbs_voc_set(1)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.07,0.06,0.24,0.41/) ! VOC + OH/O3 in low NOx
       vbs_voc_set(1)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.07,0.06,0.24,0.41/) ! VOC + OH/O3 in low NOx
       vbs_voc_set(1)%stoich_coeff_no3(2:vbs_ngroup) = (/0.07,0.06,0.24,0.41/) ! VOC + NO3
       ! ISOP (Isoprene)
       CALL new_species(mw=68., idx=vbs_voc_set(2)%spid)
       vbs_voc_set(2)%k_0_OH     = 2.7E-11
       vbs_voc_set(2)%Eact_p_OH  = 390.
       vbs_voc_set(2)%k_0_O3     = 1.03E-14
       vbs_voc_set(2)%Eact_p_O3  = -1995.
       vbs_voc_set(2)%k_0_NO3    = 3.15E-12
       vbs_voc_set(2)%Eact_p_NO3 = -450.
       ! Stoichiometric coefficients (Farina et al., Table 3):
       vbs_voc_set(2)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.02,0.02,0.0,0.0/)  ! VOC + OH/O3 in low NOx
       vbs_voc_set(2)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.02,0.02,0.0,0.0/)  ! VOC + OH/O3 in low NOx
       vbs_voc_set(2)%stoich_coeff_no3(2:vbs_ngroup) = (/0.01,0.02,0.01,0.0/) ! Isoprene + NO3

       ! aqSOA from VOCs (IEPOX and GLYX)
       IF (laqsoa) THEN
            ! Isoprene OH and O3 oxidation:
            !   ISOP+OH (->0.7 ISOPOO->0.75 IEPOX) => 0.525 IEPOX + ... (Bates et al., 2014; k=2.7E-11*exp(390/T))
            !   ISOP+OH => 0.025 GLYX (Jenkin et al., 2015)
            vbs_voc_set(2)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.525, 0.025/)
            !   ISOP+O3 => 0.01 GLYX + ... (McNeill et al., 2012; k=1.23e-14*exp(-2013/T))
            vbs_voc_set(2)%stoich_coeff_o3(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.01/)
            ! Other (hard coded): IEPOX+OH => 0.24 GLYX + ... (Jacobs et al, 2013; see below)
       ENDIF

       ! Defining VBS species
       IF (.NOT. ALLOCATED(vbs_set)) ALLOCATE(vbs_set(vbs_ngroup))
       ! equ. vapor conc. at T0 (ug=1e-9 kg converted to mol)
       vbs_set(:)%C0 = (/0.0,1.0,10.,100.,1000./)*1e-9/mwoc ! C0 [mol/m3]
       vbs_set(:)%T0 = 298. ! T0 [K]
       vbs_set(:)%Hvap_eff = 30e3/argas ! eff. evap. enthalpy [K]
       DO jv=1, vbs_ngroup
            ! Physical properties from SALSA (mw from kg/mol to g/mol)
            CALL new_species(mw=mwoc*1e3, density=densoc, kappa=kappaoc, idx=vbs_set(jv)%spid)
       ENDDO

    CASE DEFAULT
       WRITE (*,'(a,i0)') &
            'ERROR: No scheme is implemented for nvbs_setup =',&
            nvbs_setup

       STOP 'vbs_species: run terminated'

    END SELECT


    !>> thk: added aqueous phase soa directly
    IF (laqsoa) THEN

       ! allocating memory for the wet SOA set:
       IF (.NOT. ALLOCATED(aqsoa_set)) THEN
          ALLOCATE(aqsoa_set(aqsoa_ngroup))
       END IF

       ! Isoprene epoxide (IEPOX)
       CALL new_species(&
            mw          = 118.,                      & ! to be reviewed
            density     = 1320.,                     & ! to be reviewed
            kappa       = 0.037,                     & !Petters and Kreidenweis (2007)
            idx         = spid_temp                     &
            )

       ! Physical parameters
       aqsoa_set(1)%spid        = spid_temp
       aqsoa_set(1)%Eff_henry_cloud   = 1.0E5        ! Effective Henry's constant for cloud droplets
       aqsoa_set(1)%Eff_henry_aerosol = 1.0E8        ! Effective Henry's constant for aerosols
       !OH
       aqsoa_set(1)%k_0_OH     = 1.25E-11 !pre-factor [m3/(mol*s)] ((Bates et al.. Average between cis and trans isomers)
       !aqsoa_set(1)%k_0_OH     = 3.56E-11 !pre-factor [m3/(mol*s)] ((Jacobs et al.. average between IEPOX1 and IEPOX4)
       aqsoa_set(1)%Eact_p_OH  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       aqsoa_set(1)%k_0_O3     = 0.0 !pre-factor [m3/(mol*s)]
       aqsoa_set(1)%Eact_p_O3  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       aqsoa_set(1)%k_0_NO3    = 0.0 !pre-factor [m3/(mol*s)]
       aqsoa_set(1)%Eact_p_NO3 = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !Photodissociation
       aqsoa_set(1)%photodis   = 0.0 !desctruction rate [% s-1] with sunlight


       ! Glyoxal (Glyx)
       CALL new_species(&
            mw          = 58.,                      & ! to be reviewed
            density     = 1320.,                     & ! to be reviewed
            kappa       = 0.037,                     & !Petters and Kreidenweis (2007)
            idx         = spid_temp                     &
            )

       ! Physical parameters
       aqsoa_set(2)%spid        = spid_temp
       aqsoa_set(2)%Eff_henry_cloud   = 4.19E5  ! Effective Henry's constant for cloud droplets
       aqsoa_set(2)%Eff_henry_aerosol = 3.0E8  ! Effective Henry's constant for aerosols (Kampf et al. 2013)
       !OH
       aqsoa_set(2)%k_0_OH     = 3.1e-12 !pre-factor [m3/(mol*s)] (IUPAC)
       aqsoa_set(2)%Eact_p_OH  = 340.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !O3
       aqsoa_set(2)%k_0_O3     = 0.0 !pre-factor [m3/(mol*s)]
       aqsoa_set(2)%Eact_p_O3  = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !NO3
       !aqsoa_set(2)%k_0_NO3    = 6.E-13 !pre-factor [m3/(mol*s)] (McNeill et al., 2012)
       !aqsoa_set(2)%Eact_p_NO3 = -2058.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       aqsoa_set(2)%k_0_NO3    = 4.E-16 !pre-factor [m3/(mol*s)] (IUPAC)
       aqsoa_set(2)%Eact_p_NO3 = 0.0 !reduced activation energy [K*mol]: Eact_p = Eact/R
       !Photodissociation
       aqsoa_set(2)%photodis   = 3.141593/2.*8.21E-5 !desctruction rate [% s-1] with sunlight

    END IF !laqsoa

    !<< thk

  END SUBROUTINE vbs_species



  ! *****************************************
  ! Simplified version of the original new_species for UCLALES-SALSA:
  ! Take only the relevant parameters.
  SUBROUTINE new_species(mw,density,kappa,idx)

    REAL, INTENT(IN)           :: mw       ! Molecular weight (required for conversions)
    REAL, INTENT(IN), OPTIONAL :: density  ! Density
    REAL, INTENT(IN), OPTIONAL :: kappa    ! Kappa-Koehler coefficient
    INTEGER, INTENT(OUT)       :: idx      ! index to species list

    !--- Local counter for the total number of species defined
    INTEGER :: i = 0

    !--- executable procedure

    ! Increment number of species instances and store data
    i = i + 1

    spec_moleweight(i) = mw
    IF (PRESENT(density)) spec_density(i)=density
    IF (PRESENT(kappa)) spec_kappa(i)=kappa

    ! return current index
    idx = i

  END SUBROUTINE new_species


END MODULE mo_vbs
