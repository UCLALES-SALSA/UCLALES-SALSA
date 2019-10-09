MODULE classSection
  IMPLICIT NONE

  TYPE Section
     REAL :: vhilim,     & ! bin volume at the high limit
             vlolim,     & ! - '' - at the low limit
             vratiohi,   & ! volume ratio between the center and high limit
             vratiolo,   & ! - '' - and the low limit
             dmid,       & ! bin middle diameter
             !******************************************************
             ! ^ Do NOT change the stuff above after initialization !
             !******************************************************
             dwet,       & ! Wet diameter or mean droplet diameter
             numc,       & ! Number concentration of particles/droplets
             core          ! Volume of dry particle
     
     REAL ::   volc(8)     ! Volume concentrations of aerosol species + water. Only used compouds are given values. For technical reasons this cannot be allocatable at least for now....
     INTEGER :: phase   ! Water phase 1: liquid, 2: ice, 3: snow

  END TYPE Section

  INTERFACE Section
     PROCEDURE :: cnstr
  END INTERFACE Section

  CONTAINS
    
    FUNCTION cnstr(iphase)
      IMPLICIT NONE
      TYPE(Section) :: cnstr
      INTEGER, INTENT(in) :: iphase ! Phase identifier, 1: liquid, 2: ice, 3: snow

      cnstr%vhilim = 0.; cnstr%vlolim = 0.
      cnstr%vratiohi = 0; cnstr%vratiolo = 0.
      cnstr%dmid = 0.; cnstr%dwet = 0.
      cnstr%numc = 0.; cnstr%core = 0.
      cnstr%phase = iphase
      cnstr%volc(:) = 0.
    END FUNCTION cnstr



END MODULE classSection
