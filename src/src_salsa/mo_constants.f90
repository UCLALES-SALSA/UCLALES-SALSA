MODULE mo_constants

  IMPLICIT NONE

! Universal constants
  REAL, PARAMETER :: avo   = 6.02214e23    ! Avogadro constant in 1/mol

! Molar weights in g/mol
  REAL, PARAMETER :: amd   = 28.970        ! molecular weight of dry air

! Dry air and water vapour thermodynamic constants
  REAL, PARAMETER :: rd    = 287.05        ! gas constant for dry air in J/K/kg
  REAL, PARAMETER :: rv    = 461.51        ! gas constant for water vapour
                                                  ! in J/K/kg
  REAL, PARAMETER :: alv    = 2.5008e6      ! latent heat for vaporisation in J/kg
  REAL, PARAMETER :: als    = 2.8345e6      ! latent heat for sublimation in J/kg
  REAL, PARAMETER :: alf    = als-alv          ! latent heat for fusion in J/kg

! Earth and earth orbit parameters
  REAL, PARAMETER :: g     = 9.80665       ! gravity acceleration in m/s2
 
END MODULE mo_constants
