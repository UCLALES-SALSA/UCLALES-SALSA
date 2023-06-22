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

  ! Isoprene and monoterpene indexes to VOC arrays
  INTEGER, PUBLIC :: voc_id_isop=-1, voc_id_mtp=-1

  ! Initialize species
  PUBLIC :: vbs_species

CONTAINS

  SUBROUTINE vbs_species(nvbs_setup,laqsoa,densoc,mwoc,kappaoc,conc_voc)

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
    REAL, INTENT(INOUT) :: conc_voc(:)

    REAL, PARAMETER :: argas = 8.314472 ! [J/K/mol] molar/universal/ideal gas constant
    REAL :: dens_corr

    INTEGER :: jv

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
       vbs_voc_set(1)%mw = 136.
       voc_id_mtp = 1
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
          vbs_voc_set(1)%stoich_coeff_o3=(/0.1, 0.037, 0.088/)
          vbs_voc_set(1)%stoich_coeff_no3=(/0.1, 0.037, 0.088/)
          !vbs_voc_set(1)%stoich_coeff_oh=(/0.002, 0.003, 0.065/) ! from Harri (low NOx)
       ENDIF

       ! Isoprene (C5H8)
       vbs_voc_set(2)%mw = 68.
       voc_id_isop = 2
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
          vbs_voc_set(2)%stoich_coeff_o3=(/0.0, 0.0295, 0.0453, 0.525, 0.025/)
          vbs_voc_set(2)%stoich_coeff_no3=(/0.0, 0.0295, 0.0453, 0.525, 0.025/)
       ELSE
          !                            (  OC   ,  VBS1    , VBS10     )
          vbs_voc_set(2)%stoich_coeff_oh=(/0.0, 0.0295, 0.0453/)
          vbs_voc_set(2)%stoich_coeff_o3=(/0.0, 0.0295, 0.0453/)
          vbs_voc_set(2)%stoich_coeff_no3=(/0.0, 0.0295, 0.0453/)
       ENDIF

       ! Toluene (TOL)
       vbs_voc_set(3)%mw = 92.
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
          vbs_voc_set(3)%stoich_coeff_o3=(/0.36, 0.0, 0.0, 0.0, 0.24/)
          vbs_voc_set(3)%stoich_coeff_no3=(/0.36, 0.0, 0.0, 0.0, 0.24/)
       ELSE
          !                            (  OC    ,  VBS1 ,   VBS10)
          vbs_voc_set(3)%stoich_coeff_oh=(/0.36, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(3)%stoich_coeff_o3=(/0.36, 0.0, 0.0/)
          vbs_voc_set(3)%stoich_coeff_no3=(/0.36, 0.0, 0.0/)
       ENDIF

       ! Xylene (XYL)
       vbs_voc_set(4)%mw = 106.
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
          vbs_voc_set(4)%stoich_coeff_o3=(/0.30, 0.0, 0.0, 0.0, 0.25/)
          vbs_voc_set(4)%stoich_coeff_no3=(/0.30, 0.0, 0.0, 0.0, 0.25/)
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_set(4)%stoich_coeff_oh=(/0.30, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(4)%stoich_coeff_o3=(/0.30, 0.0, 0.0/)
          vbs_voc_set(4)%stoich_coeff_no3=(/0.30, 0.0, 0.0/)
       ENDIF

       ! Benzene (BENZ)
       vbs_voc_set(5)%mw = 66.
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
          vbs_voc_set(5)%stoich_coeff_o3=(/0.37, 0.0, 0.0, 0.0, 0.35/)
          vbs_voc_set(5)%stoich_coeff_no3=(/0.37, 0.0, 0.0, 0.0, 0.35/)
       ELSE
          !                            (  OC    ,  VBS1 , VBS10  )
          vbs_voc_set(5)%stoich_coeff_oh=(/0.37, 0.0, 0.0/) ! to be reviewed
          vbs_voc_set(5)%stoich_coeff_o3=(/0.37, 0.0, 0.0/)
          vbs_voc_set(5)%stoich_coeff_no3=(/0.37, 0.0, 0.0/)
       ENDIF

       ! Volatility Basis set C*=0 ug/m3
       vbs_set(1)%mw    = 186.
       vbs_set(1)%dens  = 1320.
       vbs_set(1)%kappa = 0.1 ! 0.037

       vbs_set(1)%C0          = 0.0e-6/vbs_set(1)%mw  ! equ. vapor con. at T0 [mol/m3]
       vbs_set(1)%T0          = 298               ! T0 [K]
       vbs_set(1)%Hvap_eff    = 30e3/argas        ! eff. evap. enthalpy [K]

       ! Volatility Basis set C*=1 ug/m3
       vbs_set(2)%mw    = 186.
       vbs_set(2)%dens  = 1320.
       vbs_set(2)%kappa = 0.037

       vbs_set(2)%C0          = 1.0e-6/vbs_set(2)%mw  ! equ. vapor con. at T0 [mol/m3]
       vbs_set(2)%T0          = 298.               ! T0 [K]
       vbs_set(2)%Hvap_eff    = 30e3/argas        ! eff. evap. enthalpy [K]

       ! Volatility Basis set C*=10 ug/m3
       vbs_set(3)%mw    = 186.
       vbs_set(3)%dens  = 1320.
       vbs_set(3)%kappa = 0.037

       vbs_set(3)%C0          = 10.0e-6/vbs_set(3)%mw ! equ. vapor con. at T0 [mol/m3]
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
       ! Physical properties from SALSA (mw from kg/mol to g/mol)
       vbs_set(:)%mw = mwoc*1e3
       vbs_set(:)%dens = densoc
       vbs_set(:)%kappa = kappaoc

    CASE(3) ! VBS scheme with VOCs (Farina et al., J. Geophys. Res., 115, D09202, 2010)
        ! volatility bins based on VOC oxidation scheme from Farina et al. (JGR, 2010)
        ! one non-volatile bin added for non-volatile VOC oxidation products and for VBS(g) oxidation
        ! oxidation rates are from Prank et al. (ACP, 2022) and Mielonen et al. (Atmosphere, 2018)
        ! physical properties (density, mw, kappa) taken from SALSA organics
        ! active VOCs selected based on NAMELIST input conc_voc(:)>1e-40
        ! aqSOA is possible (either on or off)

        ! VBS scheme
        vbs_ngroup = 5 ! One non-volatile bin and four semivolatile bins from Farina et al. (2010)
        vbs_nvocs = COUNT(conc_voc(1:5)>1e-40)  ! VOCs: ALPH, ISOP, TOL/ARO1, XYL/ARO2, and BENZ/ARO1

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

        jv=0
        IF (conc_voc(1)>1e-40) THEN
            jv=jv+1
            ! ALPH (MTP)
            vbs_voc_set(jv)%mw = 136
            voc_id_mtp = jv
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Prank et al., 2022)
            vbs_voc_set(jv)%k_0_OH     = 1.2E-11
            vbs_voc_set(jv)%Eact_p_OH  = 440.
            vbs_voc_set(jv)%k_0_O3     = 6.3E-16
            vbs_voc_set(jv)%Eact_p_O3  = -580.
            vbs_voc_set(jv)%k_0_NO3    = 1.2E-12
            vbs_voc_set(jv)%Eact_p_NO3 = 490.
            ! Stoichiometric coefficients (Farina et al., Table 3):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.07,0.06,0.24,0.41/) ! VOC + OH/O3 in low NOx
            vbs_voc_set(jv)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.07,0.06,0.24,0.41/) ! VOC + OH/O3 in low NOx
            vbs_voc_set(jv)%stoich_coeff_no3(2:vbs_ngroup) = (/0.07,0.06,0.24,0.41/) ! VOC + NO3
            !
            ! aqSOA: not from monoterpenes
            !
            ! Update conc_voc
            IF (jv<1) conc_voc(jv)=conc_voc(1)
        ENDIF

        IF (conc_voc(2)>1e-40) THEN
            jv=jv+1
            ! ISOP (Isoprene)
            vbs_voc_set(jv)%mw = 68.
            voc_id_isop = jv
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Prank et al., 2022)
            vbs_voc_set(jv)%k_0_OH     = 2.7E-11
            vbs_voc_set(jv)%Eact_p_OH  = 390.
            vbs_voc_set(jv)%k_0_O3     = 1.03E-14
            vbs_voc_set(jv)%Eact_p_O3  = -1995.
            vbs_voc_set(jv)%k_0_NO3    = 3.15E-12
            vbs_voc_set(jv)%Eact_p_NO3 = -450.
            ! Stoichiometric coefficients (Farina et al., Table 3):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.02,0.02,0.0,0.0/)  ! VOC + OH/O3 in low NOx
            vbs_voc_set(jv)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.02,0.02,0.0,0.0/)  ! VOC + OH/O3 in low NOx
            vbs_voc_set(jv)%stoich_coeff_no3(2:vbs_ngroup) = (/0.01,0.02,0.01,0.0/) ! Isoprene + NO3
            !
            ! aqSOA (IEPOX and GLYX) from isoprene OH and O3 oxidation:
            IF (laqsoa) THEN
                ! ISOP+OH (->0.7 ISOPOO->0.75 IEPOX) => 0.525 IEPOX + ... (Bates et al., 2014; k=2.7E-11*exp(390/T))
                ! ISOP+OH => 0.025 GLYX (Jenkin et al., 2015)
                vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.525, 0.025/)
                ! ISOP+O3 => 0.01 GLYX + ... (McNeill et al., 2012; k=1.23e-14*exp(-2013/T))
                vbs_voc_set(jv)%stoich_coeff_o3(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.01/)
                ! Other (hard coded): IEPOX+OH => 0.24 GLYX + ... (Jacobs et al, 2013)
            ENDIF
            !
            ! Update conc_voc
            IF (jv<2) conc_voc(jv)=conc_voc(2)
        ENDIF

        IF (conc_voc(3)>1e-40) THEN
            jv=jv+1
            ! ARO1 or toluene (TOL)
            vbs_voc_set(jv)%mw = 92.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 1.81E-12
            vbs_voc_set(jv)%Eact_p_OH  = 338.
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Farina et al., Table 3):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.01,0.22,0.53,0.64/) ! VOC + OH/O3 in low NOx
            !
            ! aqSOA: TOL + OH -> ... + 0.238 GLYX +...   k=1.81e-12*exp(338/T) (McNeil et al., 2012)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.24/)
            !
            ! Update conc_voc
            IF (jv<3) conc_voc(jv)=conc_voc(3)
        ENDIF

        IF (conc_voc(4)>1e-40) THEN
            jv=jv+1
            ! ARO2 or xylene (XYL)
            vbs_voc_set(jv)%mw = 106.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 2.31E-11
            vbs_voc_set(jv)%Eact_p_OH  = 0.0
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Farina et al., Table 3):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.01,0.22,0.53,0.64/) ! VOC + OH/O3 in low NOx
            !
            ! aqSOA: XYL + OH -> ... + 0.247 GLYX + ...   k=1.43e-11 (McNeil et al., 2012)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.25/)
            !
            ! Update conc_voc
            IF (jv<4) conc_voc(jv)=conc_voc(4)
        ENDIF

        IF (conc_voc(5)>1e-40) THEN
            jv=jv+1
            ! ARO1 or benzene (BENZ)
            vbs_voc_set(jv)%mw = 78.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 2.33E-12
            vbs_voc_set(jv)%Eact_p_OH  = -193.
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Farina et al., Table 3):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.01,0.22,0.53,0.64/) ! VOC + OH/O3 in low NOx
            !
            ! aqSOA: BENZ + OH -> 0.35 GLYX (Mielonen et al., 2018)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.35/)
            !
            ! Update conc_voc
            IF (jv<5) conc_voc(jv)=conc_voc(5)
        ENDIF

        IF (jv /= vbs_nvocs) STOP 'vbs_species: error in setting up VOC species!'

        ! Defining VBS species
        IF (.NOT. ALLOCATED(vbs_set)) ALLOCATE(vbs_set(vbs_ngroup))
        ! equ. vapor conc. at T0 (ug=1e-9 kg converted to mol)
        vbs_set(:)%C0 = (/0.0,1.0,10.,100.,1000./)*1e-9/mwoc ! C0 [mol/m3]
        vbs_set(:)%T0 = 298. ! T0 [K]
        vbs_set(:)%Hvap_eff = 30e3/argas ! eff. evap. enthalpy [K]
        ! Physical properties from SALSA (mw from kg/mol to g/mol)
        vbs_set(:)%mw = mwoc*1e3
        vbs_set(:)%dens = densoc
        vbs_set(:)%kappa = kappaoc

        ! Density correction based on differences between LES and VBS scheme (assuming 1.4 g/cm3)
        dens_corr = densoc/1400.
        DO jv=1,vbs_nvocs
            vbs_voc_set(jv)%stoich_coeff_oh(:) = vbs_voc_set(jv)%stoich_coeff_oh(:) * dens_corr
            vbs_voc_set(jv)%stoich_coeff_o3(:) = vbs_voc_set(jv)%stoich_coeff_o3(:) * dens_corr
            vbs_voc_set(jv)%stoich_coeff_no3(:) = vbs_voc_set(jv)%stoich_coeff_no3(:) * dens_corr
        END DO

    CASE(4) ! VBS scheme with VOCs (Lane et al., Atmos. Environ., 42, 7539-7451, 2008)
        ! volatility bins based on VOC oxidation scheme from Lane et al. (AE, 2008)
        ! one non-volatile bin added for non-volatile VOC oxidation products and for VBS(g) oxidation
        ! oxidation rates are from Prank et al. (ACP, 2022) and Mielonen et al. (Atmosphere, 2018)
        ! physical properties (density, mw, kappa) taken from SALSA organics
        ! active VOCs selected based on NAMELIST input conc_voc(:)>1e-40
        ! aqSOA is possible (either on or off)

        ! VBS scheme
        vbs_ngroup = 5 ! One non-volatile bin and four semivolatile bins from Lane et al. (2008)
        vbs_nvocs = COUNT(conc_voc(1:5)>1e-40)  ! VOCs: ALPH, ISOP, TOL/ARO1, XYL/ARO2, and BENZ/ARO1

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

        jv=0
        IF (conc_voc(1)>1e-40) THEN
            jv=jv+1
            ! ALPH (MTP)
            vbs_voc_set(jv)%mw = 136
            voc_id_mtp = jv
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Prank et al., 2022)
            vbs_voc_set(jv)%k_0_OH     = 1.2E-11
            vbs_voc_set(jv)%Eact_p_OH  = 440.
            vbs_voc_set(jv)%k_0_O3     = 6.3E-16
            vbs_voc_set(jv)%Eact_p_O3  = -580.
            vbs_voc_set(jv)%k_0_NO3    = 1.2E-12
            vbs_voc_set(jv)%Eact_p_NO3 = 490.
            ! Stoichiometric coefficients (Lane et al., Table 2):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.072,0.061,0.239,0.405/) ! Base case
            vbs_voc_set(jv)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.072,0.061,0.239,0.405/) ! Base case
            vbs_voc_set(jv)%stoich_coeff_no3(2:vbs_ngroup) = (/0.072,0.061,0.239,0.405/) ! Base case
            !
            ! aqSOA: not from monoterpenes
            !
            ! Update conc_voc
            IF (jv<1) conc_voc(jv)=conc_voc(1)
        ENDIF

        IF (conc_voc(2)>1e-40) THEN
            jv=jv+1
            ! ISOP (Isoprene)
            vbs_voc_set(jv)%mw = 68.
            voc_id_isop = jv
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Prank et al., 2022)
            vbs_voc_set(jv)%k_0_OH     = 2.7E-11
            vbs_voc_set(jv)%Eact_p_OH  = 390.
            vbs_voc_set(jv)%k_0_O3     = 1.03E-14
            vbs_voc_set(jv)%Eact_p_O3  = -1995.
            vbs_voc_set(jv)%k_0_NO3    = 3.15E-12
            vbs_voc_set(jv)%Eact_p_NO3 = -450.
            ! Stoichiometric coefficients (Lane et al., Table 2):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup)  = (/0.006,0.02,0.01,0.0/)  ! Base case
            vbs_voc_set(jv)%stoich_coeff_o3(2:vbs_ngroup)  = (/0.006,0.02,0.01,0.0/)  ! Base case
            vbs_voc_set(jv)%stoich_coeff_no3(2:vbs_ngroup) = (/0.006,0.02,0.01,0.0/)  ! Base case
            !
            ! aqSOA (IEPOX and GLYX) from isoprene OH and O3 oxidation:
            IF (laqsoa) THEN
                ! ISOP+OH (->0.7 ISOPOO->0.75 IEPOX) => 0.525 IEPOX + ... (Bates et al., 2014; k=2.7E-11*exp(390/T))
                ! ISOP+OH => 0.025 GLYX (Jenkin et al., 2015)
                vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.525, 0.025/)
                ! ISOP+O3 => 0.01 GLYX + ... (McNeill et al., 2012; k=1.23e-14*exp(-2013/T))
                vbs_voc_set(jv)%stoich_coeff_o3(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.01/)
                ! Other (hard coded): IEPOX+OH => 0.24 GLYX + ... (Jacobs et al, 2013)
            ENDIF
            !
            ! Update conc_voc
            IF (jv<2) conc_voc(jv)=conc_voc(2)
        ENDIF

        IF (conc_voc(3)>1e-40) THEN
            jv=jv+1
            ! ARO1 or toluene (TOL)
            vbs_voc_set(jv)%mw = 92.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 1.81E-12
            vbs_voc_set(jv)%Eact_p_OH  = 338.
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Lane et al., Table 2):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.01,0.03,0.075,0.25/) ! Base case
            !
            ! aqSOA: TOL + OH -> ... + 0.238 GLYX +...   k=1.81e-12*exp(338/T) (McNeil et al., 2012)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.24/)
            !
            ! Update conc_voc
            IF (jv<3) conc_voc(jv)=conc_voc(3)
        ENDIF

        IF (conc_voc(4)>1e-40) THEN
            jv=jv+1
            ! ARO2 or xylene (XYL)
            vbs_voc_set(jv)%mw = 106.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 2.31E-11
            vbs_voc_set(jv)%Eact_p_OH  = 0.0
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Lane et al., Table 2):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.02,0.04,0.08,0.25/) ! Base case
            !
            ! aqSOA: XYL + OH -> ... + 0.247 GLYX + ...   k=1.43e-11 (McNeil et al., 2012)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.25/)
            !
            ! Update conc_voc
            IF (jv<4) conc_voc(jv)=conc_voc(4)
        ENDIF

        IF (conc_voc(5)>1e-40) THEN
            jv=jv+1
            ! ARO1 or benzene (BENZ)
            vbs_voc_set(jv)%mw = 78.
            ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Mielonen et al., 2018)
            vbs_voc_set(jv)%k_0_OH     = 2.33E-12
            vbs_voc_set(jv)%Eact_p_OH  = -193.
            vbs_voc_set(jv)%k_0_O3     = 0.0
            vbs_voc_set(jv)%Eact_p_O3  = 0.0
            vbs_voc_set(jv)%k_0_NO3    = 0.0
            vbs_voc_set(jv)%Eact_p_NO3 = 0.0
            ! Stoichiometric coefficients (Lane et al., Table 2):
            vbs_voc_set(jv)%stoich_coeff_oh(2:vbs_ngroup) = (/0.01,0.03,0.075,0.25/) ! Base case
            !
            ! aqSOA: BENZ + OH -> 0.35 GLYX (Mielonen et al., 2018)
            IF (laqsoa) vbs_voc_set(jv)%stoich_coeff_oh(vbs_ngroup+1:vbs_ngroup+2)=(/0.0, 0.35/)
            !
            ! Update conc_voc
            IF (jv<5) conc_voc(jv)=conc_voc(5)
        ENDIF

        IF (jv /= vbs_nvocs) STOP 'vbs_species: error in setting up VOC species!'

        ! Defining VBS species
        IF (.NOT. ALLOCATED(vbs_set)) ALLOCATE(vbs_set(vbs_ngroup))
        ! equ. vapor conc. at T0 (ug=1e-9 kg converted to mol)
        vbs_set(:)%C0 = (/0.0,1.0,10.,100.,1000./)*1e-9/mwoc ! C0 [mol/m3]
        vbs_set(:)%T0 = 298. ! T0 [K]
        vbs_set(:)%Hvap_eff = 30e3/argas ! eff. evap. enthalpy [K]
        ! Physical properties from SALSA (mw from kg/mol to g/mol)
        vbs_set(:)%mw = mwoc*1e3
        vbs_set(:)%dens = densoc
        vbs_set(:)%kappa = kappaoc

        ! Density correction based on differences between LES and VBS scheme (assuming 1.0 g/cm3)
        dens_corr = densoc/1000.
        DO jv=1,vbs_nvocs
            vbs_voc_set(jv)%stoich_coeff_oh(:) = vbs_voc_set(jv)%stoich_coeff_oh(:) * dens_corr
            vbs_voc_set(jv)%stoich_coeff_o3(:) = vbs_voc_set(jv)%stoich_coeff_o3(:) * dens_corr
            vbs_voc_set(jv)%stoich_coeff_no3(:) = vbs_voc_set(jv)%stoich_coeff_no3(:) * dens_corr
        END DO

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
       aqsoa_set(1)%mw    = 118.
       aqsoa_set(1)%dens  = 1320.
       aqsoa_set(1)%kappa = 0.037

       ! Effective Henry's law constants for cloud droplets and aerosols [M/atm]
       aqsoa_set(1)%Eff_henry_cloud   = 1.0E5
       aqsoa_set(1)%Eff_henry_aerosol = 1.0E8 ! Nguyen et al. (2014)
       ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Arrhenius)
       !OH
       aqsoa_set(1)%k_0_OH     = 1.25E-11 ! Bates et al. (2014): average between cis and trans isomers
       !aqsoa_set(1)%k_0_OH     = 3.56E-11 ! Jacobs et al. (2013): average between IEPOX1 and IEPOX4)
       aqsoa_set(1)%Eact_p_OH  = 0.0
       !O3
       aqsoa_set(1)%k_0_O3     = 0.0
       aqsoa_set(1)%Eact_p_O3  = 0.0
       !NO3
       aqsoa_set(1)%k_0_NO3    = 0.0
       aqsoa_set(1)%Eact_p_NO3 = 0.0
       ! Photodissociation rate [1/s]
       aqsoa_set(1)%photodis   = 0.0


       ! Glyoxal (Glyx)
       aqsoa_set(2)%mw    = 58.
       aqsoa_set(2)%dens  = 1320.
       aqsoa_set(2)%kappa = 0.037

       ! Effective Henry's law constants for cloud droplets and aerosols [M/atm]
       aqsoa_set(2)%Eff_henry_cloud   = 4.19E5
       aqsoa_set(2)%Eff_henry_aerosol = 3.0E8  ! Kampf et al. (2013)
       ! Oxidation rates [cm3/(molecule*s)] and their temperature dependence [T] (Arrhenius)
       !OH
       !aqsoa_set(2)%k_0_OH     = 11.4e-12 ! McNeill et al. (2012)
       !aqsoa_set(2)%Eact_p_OH  = 0.0
       aqsoa_set(2)%k_0_OH     = 3.1e-12 ! IUPAC
       aqsoa_set(2)%Eact_p_OH  = 340.0
       !O3
       aqsoa_set(2)%k_0_O3     = 0.0
       aqsoa_set(2)%Eact_p_O3  = 0.0
       !NO3
       !aqsoa_set(2)%k_0_NO3    = 6.E-13 ! McNeill et al. (2012)
       !aqsoa_set(2)%Eact_p_NO3 = -2058.0
       aqsoa_set(2)%k_0_NO3    = 4.E-16 ! IUPAC
       aqsoa_set(2)%Eact_p_NO3 = 0.0
       ! Photodissociation rate [1/s]
       aqsoa_set(2)%photodis   = 0.0 !  3.141593/2.*8.21E-5

    END IF !laqsoa

    !<< thk

  END SUBROUTINE vbs_species

END MODULE mo_vbs
