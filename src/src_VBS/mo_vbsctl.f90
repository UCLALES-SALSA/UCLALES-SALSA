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
     ! cross-referencing (for SALSA)
     INTEGER :: id_gas      ! gas phase tracer index

     ! OH, O3 and NO3 oxidation rates and their temperature dependence (Arrhenius)
     REAL :: k_0_OH, k_0_O3, k_0_NO3          ! pre-factor [cm3/(molecule*s)]
     REAL :: Eact_p_OH, Eact_p_O3, Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     ! physical properties
     REAL :: mw ! molecular weight (g/mol)

     ! define what fraction goes where (stoichiometric coefficients)
     REAL, ALLOCATABLE :: stoich_coeff_oh(:), stoich_coeff_o3(:), stoich_coeff_no3(:)
  END TYPE t_voc_prec

  ! Collection of the volatility base parameters
  TYPE :: t_vbs_group
     ! cross-referencing (for SALSA)
     INTEGER :: id_gas      ! gas phase tracer index
     INTEGER :: id_vols     ! index in vols structure

     ! effective saturation concentration
     REAL :: C0         ! saturation concentration at reference temperature [mol/m3]
     REAL :: T0         ! reference temperature [K]
     REAL :: Hvap_eff   ! effective heat of vaporization Hvap_eff=Hvap/R  [K]

     ! physical properties
     REAL :: mw, dens, kappa ! molecular weight (g/mol), density (kg/m3) and kappa (-)
  END TYPE t_vbs_group

  ! Collection of the aqsoa parameters
  TYPE :: t_aq_soa
     ! cross-referencing (for SALSA)
     INTEGER :: id_gas      ! gas phase tracer index
     INTEGER :: id_vols     ! index in vols structure

     ! OH, O3 and NO3 oxidation rates and their temperature dependence (Arrhenius)
     REAL :: k_0_OH, k_0_O3, k_0_NO3          ! pre-factor [cm3/(molecule*s)]
     REAL :: Eact_p_OH, Eact_p_O3, Eact_p_NO3 ! reduced activation energy [K]: Eact_p = Eact/R

     ! photodissociation rate [1/s]
     REAL :: photodis

     ! effective Henry's law constants [M/atm]
     REAL :: Eff_henry_aerosol
     REAL :: Eff_henry_cloud

     ! physical properties
     REAL :: mw, dens, kappa ! molecular weight (g/mol), density (kg/m3) and kappa (-)
  END TYPE t_aq_soa


  ! these will be set in mo_vbs
  TYPE(t_voc_prec), ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_voc_set
  TYPE(t_vbs_group), ALLOCATABLE, TARGET, DIMENSION(:) :: vbs_set
  TYPE(t_aq_soa), ALLOCATABLE, TARGET, DIMENSION(:) :: aqsoa_set

  INTEGER :: vbs_nvocs = 0    ! number of VBS precursors
  INTEGER :: vbs_ngroup = 0   ! number of VBS groups
  INTEGER :: aqsoa_ngroup = 0 ! number of wet SOA groups

END MODULE mo_vbsctl
