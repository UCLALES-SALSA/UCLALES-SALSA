!****************************************************************
! Module mo_vbsctl: data types for the gas phase species
!
! To be used with UCLALES-SALSA only. Data types contain parameters for
! the gas phase species including volatile organics compounds (VOC),
! volatility basis set (VBS) species, and aqSOA species.
!
! Thomas Kuehn (UEF)  - initial setup
! Tomi Raatikainen (FMI) - modified for UCLALES-SALSA (2019)
!
!****************************************************************

MODULE mo_vbsctl

  IMPLICIT NONE

  ! Collection of parameters for the volatility base precursers
  TYPE :: t_voc_prec
     !cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: id_gas      ! gas phase tracer index

     !oxidation rates and their temperature dependence (Arrhenius)
     !OH
     REAL :: k_0_OH     ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_OH  ! reduced activation energy [K]: Eact_p = Eact/R
     !O3
     REAL :: k_0_O3     ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_O3  ! reduced activation energy [K]: Eact_p = Eact/R
     !NO3
     REAL :: k_0_NO3    ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     ! define what fraction goes where
     ! (vector length depends on size of basis set)
     REAL, ALLOCATABLE :: stoich_coeff(:) !Stoichiometric Coefficients
  END TYPE t_voc_prec

  ! Collection of the volatility base parameters
  TYPE :: t_vbs_group
     ! cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: id_gas      ! gas phase tracer index
     INTEGER :: id_vols     ! index in vols structure (for SALSA)

     ! effective saturation concentration -- Farina et al (2010)
     REAL :: C0         ! saturation concentration at reference temperature [kg/m3]
     REAL :: T0         ! reference temperature [K]
     REAL :: Hvap_eff   ! effective heat of vaporization Hvap_eff=Hvap/R  [K]

     ! physical properties
     REAL :: mv         ! molecular volume [m3]
  END TYPE t_vbs_group


  ! Collection of the aqsoa parameters
  TYPE :: t_aq_soa
     ! cross-referencing
     INTEGER :: spid        ! species index
     INTEGER :: id_gas      ! gas phase tracer index
     INTEGER :: id_vols    ! index in vols structure (for SALSA)

     !oxidation rates and their temperature dependence (Arrhenius)
     !OH
     REAL :: k_0_OH     ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_OH  ! reduced activation energy [K]: Eact_p = Eact/R
     !O3
     REAL :: k_0_O3     ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_O3  ! reduced activation energy [K]: Eact_p = Eact/R
     !NO3
     REAL :: k_0_NO3    ! pre-factor [mol/(m3*s)]
     REAL :: Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     !Photodissociation
     REAL :: photodis   ! pre-factor [rate s-1]

     !Effective Henry coefficients
     REAL :: Eff_henry_aerosol   ! pre-factor [rate s-1]
     REAL :: Eff_henry_cloud   ! pre-factor [rate s-1]

     ! physical properties
     REAL :: mv         ! molecular volume
  END TYPE t_aq_soa


  ! the following four variable will be set in
  ! *vbs_species* in *mo_vbs*
  TYPE(t_voc_prec), ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_voc_set
  TYPE(t_vbs_group), ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_set
  TYPE(t_aq_soa), ALLOCATABLE, TARGET, DIMENSION(:) :: aqsoa_set

  INTEGER :: vbs_nvocs = 0    ! number of VBS precursors
  INTEGER :: vbs_ngroup = 0   ! number of VBS groups
  INTEGER :: aqsoa_ngroup = 0 ! number of wet SOA groups

END MODULE mo_vbsctl
