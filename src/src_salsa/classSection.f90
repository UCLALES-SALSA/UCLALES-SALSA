MODULE classSection
  IMPLICIT NONE

  REAL, PARAMETER, PRIVATE :: pi6 = 0.5235988  
  ! Copied from mo_submctl to avoid cyclic dependencies. Some rearranging needed to avoid this situation
  ! ALSO/RELATED: classSection dependency of mo_submctl should be removed to stuff could be imported from there to here as well!!!!


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
             ddry,       & ! Dry diameter
             numc,       & ! Number concentration of particles/droplets
             core          ! Volume of dry particle
     
     REAL :: volc(8)     ! Volume concentrations of aerosol species + water. These are taken as bulk volume.
                         ! Only used compouds are given values. 
                         ! For technical reasons this cannot be allocatable at least for now....

     REAL :: rhoimean     ! The mean ice density for frozen particles. Takes into account only the bulk ice composition

     REAL :: rhoieff      ! Effective ice density. This should take into account non-spherical shape as well as the 
                          ! bulk ice composition



     ! These are only used for frozen hydrometeors
     !-----------------------------------------------------------
     !REAL :: virimebulk   ! Bulk volume concentration of rimed ice for frozen hydrometeors

     REAL :: vrime         ! The true volume concentration of rimed ice, taking into account particle shape
     !REAL :: viprist      ! The true volume concentration of pristine ice, taking into account the particle shape

     INTEGER :: phase    ! Phase identifier, 1: aerosol 2: cloud droplet, 3: precip, 4: ice, 5: snow
     REAL    :: nlim     ! Lower limit for number concentration used in many calculations, depends on particle type
     REAL    :: dlim     ! Category specific diameter limit used e.g. in coagulation calculations

     CONTAINS
       PROCEDURE :: updateDiameter
       PROCEDURE :: updateRhomean


  END TYPE Section
  INTERFACE Section
     PROCEDURE :: cnstr
  END INTERFACE Section

  CONTAINS
    
    FUNCTION cnstr(iphase, inlim, idlim)
      TYPE(Section) :: cnstr
      INTEGER, INTENT(in) :: iphase ! Phase identifier: 1: aerosol, 2: cloud droplets, 3: precip, 4: ice, 5: snow
      REAL, INTENT(in) :: inlim   ! Lower limit for number concentration depending on particle type
      REAL, INTENT(in) :: idlim

      cnstr%vhilim = 0.; cnstr%vlolim = 0.
      cnstr%vratiohi = 0; cnstr%vratiolo = 0.
      cnstr%dmid = 0.; cnstr%dwet = 0.
      cnstr%numc = 0.; cnstr%core = 0.
      cnstr%volc(:) = 0.
      !cnstr%virimebulk = 0.
      cnstr%vrime = 0.
      !cnstr%viprist = 0.
      cnstr%phase = iphase
      cnstr%nlim = inlim
      cnstr%dlim = idlim

    END FUNCTION cnstr

    ! 
    ! Function for calculating dimension (or wet diameter) for any particle type
    ! - Aerosol, cloud and rain are spherical
    ! - Snow and ice can be irregular and their densities can be size-dependent
    !
    ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
    ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
    !
    SUBROUTINE updateDiameter(SELF,limit)
      CLASS(Section), INTENT(inout) :: SELF
      LOGICAL,INTENT(in) :: limit  ! True -> constrain the result wet diameter by dlim

      SELF%dwet = 3.e-9
      IF (SELF%numc > SELF%nlim) &
           SELF%dwet=(SUM(SELF%volc(:))/SELF%numc/pi6)**(1./3.)

      IF (limit) SELF%dwet = MIN(SELF%dwet,SELF%dlim)
      
    END SUBROUTINE updateDiameter

    ! 
    ! Subroutine updateRhoeff
    ! Updates the effective density of the particle
    ! -------------------------------------------------
    !
    SUBROUTINE updateRhomean(SELF, rhoic, rhori, iwa)
      CLASS(Section), INTENT(inout) :: SELF
      REAL, INTENT(in) :: rhoic, rhori  ! THESE COULD BE REMOVED IF CLASSSECTION dependency in mo_submctl would be removed!! Could then just import spec
      INTEGER, INTENT(in) :: iwa   ! The same for this one

      REAL :: mass_p, mass_r, mass_t

      ! convert to masses -> get the mass mean density
      
      mass_p = rhoic*(SELF%volc(iwa) - SELF%vrime)
      mass_r = rhori*SELF%vrime
      mass_t = mass_p + mass_r

      SELF%rhoimean = (mass_p*rhoic + mass_r*rhori)/mass_t

    END SUBROUTINE updateRhomean




END MODULE classSection
