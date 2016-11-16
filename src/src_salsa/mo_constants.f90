MODULE mo_constants

  IMPLICIT NONE

! Universal constants
  REAL, PARAMETER :: api   = 3.14159265358979323846 ! pi
  REAL, PARAMETER :: asqrt2= 1.41421356237309504880
! REAL, PARAMETER :: ar    = 8.314e3       ! universal gas constant in J/K/kmol
  REAL, PARAMETER :: argas = 8.314409      ! universal gas constant (R) in J/K/mol
  REAL, PARAMETER :: avo   = 6.02214e23    ! Avogadro constant in 1/mol
  REAL, PARAMETER :: ak    = 1.380662e-23  ! Boltzmann constant in J/K
  REAL, PARAMETER :: stbo  = 5.67E-8       ! Stephan-Boltzmann constant in W/m2/K4
  REAL, PARAMETER :: ap0   = 101325.       ! Standard pressure in Pa

! Molar weights in g/mol
  REAL, PARAMETER :: amco2 = 44.011        ! molecular weight of carbon dioxide
  REAL, PARAMETER :: amch4 = 16.043        ! molecular weight of methane
  REAL, PARAMETER :: amo3  = 47.9982       ! molecular weight of ozone
  REAL, PARAMETER :: amn2o = 44.013        ! molecular weight of N2O
  REAL, PARAMETER :: amc11 =137.3686       ! molecular weight of CFC11
  REAL, PARAMETER :: amc12 =120.9140       ! molecular weight of CFC12
  REAL, PARAMETER :: amo2  = 31.9988       ! molecular weight of molecular oxygen
  REAL, PARAMETER :: amw   = 18.0154       ! molecular weight of water vapor
  REAL, PARAMETER :: amd   = 28.970        ! molecular weight of dry air

! Dry air and water vapour thermodynamic constants
  REAL, PARAMETER :: cpd   = 1005.46       ! specific heat of dry air at constant
                                                  ! pressure in J/K/kg
  REAL, PARAMETER :: cpv   = 1869.46       ! specific heat of water vapour at
                                                  ! constant pressure in J/K/kg
  REAL, PARAMETER :: rd    = 287.05        ! gas constant for dry air in J/K/kg
  REAL, PARAMETER :: rv    = 461.51        ! gas constant for water vapour
                                                  ! in J/K/kg
  REAL, PARAMETER :: rcpd  = 1.0/cpd       ! auxiliary constant in K*kg/J
  REAL            :: vtmpc1= rv/rd-1.0     ! dimensionless auxiliary constant
  REAL            :: vtmpc2= cpv/cpd-1.0   ! dimensionless auxiliary constant

! H2O related constants  (liquid, ice, snow), phase change constants
  REAL, PARAMETER :: rhoh2o = 1000.0        ! density of liquid water in kg/m3
  REAL, PARAMETER :: rhosea = 1025.0        ! density of sea water in kg/m3
  REAL, PARAMETER :: rhoice = 917.0         ! density of ice in kg/m3
  REAL, PARAMETER :: rhosno = 300.0         ! density of snow in kg/m3
  REAL, PARAMETER :: rhoiw  = rhoice/rhoh2o    ! density ratio (ice/water)
  REAL, PARAMETER :: alv    = 2.5008e6      ! latent heat for vaporisation in J/kg
  REAL, PARAMETER :: als    = 2.8345e6      ! latent heat for sublimation in J/kg
  REAL, PARAMETER :: alf    = als-alv          ! latent heat for fusion in J/kg
  REAL, PARAMETER :: cpliq  = 4218.         ! specific heat for liquid water J/K/kg
  REAL, PARAMETER :: cpsea  = 3994.         ! specific heat for sea water J/K/kg
  REAL, PARAMETER :: cpice  = 2106.         ! specific heat for ice J/K/kg
  REAL, PARAMETER :: cpsno  = 2090.         ! specific heat for snow J/K/kg
  REAL, PARAMETER :: alice  = 2.1656        ! thermal conductivity of ice in W/K/m
  REAL, PARAMETER :: alsno  = 0.31          ! thermal conductivity of snow in W/K/m
  REAL, PARAMETER :: tmelt  = 273.15        ! melting temperature of ice/snow in K

! Earth and earth orbit parameters
  REAL, PARAMETER :: a     = 6371000.0     ! radius of the earth in m
  REAL, PARAMETER :: rae   = 0.1277E-2     ! ratio of atmosphere to earth radius
  REAL, PARAMETER :: omega = .7292E-4      ! solid rotation velocity of the earth
                                                  ! in 1/s
  REAL, PARAMETER :: secperday = 86400.    ! seconds per day
  REAL, PARAMETER :: g     = 9.80665       ! gravity acceleration in m/s2
  REAL, PARAMETER :: qg    = 1.0/g         ! inverse of gravity acceleration
 
! Constants used for computation of saturation mixing ratio
! over liquid water (*c_les*) or ice(*c_ies*)
  REAL, PARAMETER :: c1es  = 610.78           !
  REAL, PARAMETER :: c2es  = c1es*rd/rv          !
  REAL, PARAMETER :: c3les = 17.269           !
  REAL, PARAMETER :: c3ies = 21.875           !
  REAL, PARAMETER :: c4les = 35.86            !
  REAL, PARAMETER :: c4ies = 7.66             !
  REAL, PARAMETER :: c5les = c3les*(tmelt-c4les) !
  REAL, PARAMETER :: c5ies = c3ies*(tmelt-c4ies) !
  REAL, PARAMETER :: c5alvcp = c5les*alv/cpd     !
  REAL, PARAMETER :: c5alscp = c5ies*als/cpd     !
  REAL, PARAMETER :: alvdcp  = alv/cpd           !
  REAL, PARAMETER :: alsdcp  = als/cpd           !

! Specifications, thresholds, and derived constants for the following subroutines:
! s_lake, s_licetemp, s_sicetemp, meltpond, meltpond_ice, update_albedo_ice_meltpond

  REAL, PARAMETER :: dmix     = 10.0   ! mixed-layer depth of lakes in m
  REAL, PARAMETER :: dmixsea  = 50.0   ! mixed-layer depth of ocean in m
  REAL, PARAMETER :: dice     = 0.05   ! minimum ice thickness in m
  REAL, PARAMETER :: dicepond = 0.01   ! minimum ice thickness of pond ice in m
  REAL, PARAMETER :: dicelim  = 0.10   ! threshold ice thickness for pond closing in m
  REAL, PARAMETER :: dpondmin = 0.01   ! minimum pond depth for pond fraction in m
  REAL, PARAMETER :: albpondi = 0.30   ! albedo of pond ice

  REAL, PARAMETER :: snicecond = alice/alsno * rhoh2o/rhosno
  REAL, PARAMETER :: hcapmix   = rhoh2o*cpliq*dmix     ! heat capacity of lake mixed layer in J/K/m2
  REAL, PARAMETER :: hcapice   = rhoice*cpice*dice     ! heat capacity of upper ice layer
  REAL, PARAMETER :: hcapicep  = rhoice*cpice*dicepond ! heat capacity of upper pond ice layer
  REAL, PARAMETER :: rhoilf    = rhoice*alf            ! [J/m3]
  REAL, PARAMETER :: rhowlf    = rhoh2o*alf            ! [J/m3]
  REAL, PARAMETER :: hcaprilf  = hcapmix/rhoilf        ! [m/K]
  REAL, PARAMETER :: rilfhcap  = rhoilf/hcapmix        ! [K/m]
  REAL, PARAMETER :: tfreez    = dice*rilfhcap         ! cooling below tmelt required to form dice
 
END MODULE mo_constants
