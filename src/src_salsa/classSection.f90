MODULE classSection
  USE mo_submctl, ONLY : pi6, spec
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
             ddry,       & ! Dry diameter (actual)
             dins,       & ! Diameter of the insoluble part (if present)
             numc,       & ! Number concentration of particles/droplets
             core          ! Volume of dry particle
     
     REAL :: volc(8)     ! Volume concentrations of aerosol species + water. These are taken as bulk volume.
                         ! Only used compouds are given values. 
                         ! For technical reasons this cannot be allocatable at least for now....

     REAL :: rhomean     ! The mean ice density for frozen particles. Takes into account only the bulk ice composition

     REAL :: rhoeff      ! Effective ice density. This should take into account non-spherical shape as well as the 
                          ! bulk ice composition



     ! These are only used for frozen hydrometeors
     !-----------------------------------------------------------
     !REAL :: virimebulk   ! Bulk volume concentration of rimed ice for frozen hydrometeors

     REAL :: vrime         ! The true volume concentration of rimed ice, taking into account particle shape
     !REAL :: viprist      ! The true volume concentration of pristine ice, taking into account the particle shape

     INTEGER :: phase    ! Phase identifier, 1: aerosol 2: cloud droplet, 3: precip, 4: ice
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
      INTEGER, INTENT(in) :: iphase ! Phase identifier: 1: aerosol, 2: cloud droplets, 3: precip, 4: ice
      REAL, INTENT(in) :: inlim   ! Lower limit for number concentration depending on particle type
      REAL, INTENT(in) :: idlim

      cnstr%vhilim = 0.; cnstr%vlolim = 0.
      cnstr%vratiohi = 0; cnstr%vratiolo = 0.
      cnstr%dmid = 0.; cnstr%dwet = 0.; cnstr%dins = 0.
      cnstr%numc = 0.; cnstr%core = 0.
      cnstr%volc(:) = 0.
      cnstr%vrime = 0.
      cnstr%phase = iphase
      cnstr%nlim = inlim
      cnstr%dlim = idlim
      cnstr%rhomean = 0.
      cnstr%rhoeff = 0.

    END FUNCTION cnstr

    ! 
    ! Function for calculating dimension (or wet diameter) for any particle type
    ! - Aerosol, cloud and rain are spherical
    ! - Snow and ice can be irregular and their densities can be size-dependent
    !
    ! Correct dimension is needed for irregular particles (e.g. ice and snow) for calculating fall speed (deposition and coagulation)
    ! and capacitance (condensation). Otherwise compact spherical structure can be expected,
    !
    ! This should be modified to account for the particle shape as well.
    !
    SUBROUTINE updateDiameter(SELF,limit,type)
      CLASS(Section), INTENT(inout) :: SELF
      LOGICAL,INTENT(in) :: limit  ! True -> constrain the result wet diameter by dlim
      CHARACTER(len=3), INTENT(in), OPTIONAL :: type  ! "dry" or "wet", "ins" (insoluble) or "all". Default is "wet"

      CHARACTER(len=3) :: swtyp
      INTEGER :: nwet, ndry

      nwet = spec%getNSpec(type="wet")
      ndry = spec%getNspec(type="dry")
      
      IF ( .NOT. PRESENT(type) ) THEN
         swtyp = "wet"
      ELSE
         swtyp = type
      END IF
       
      IF (SELF%numc > SELF%nlim) THEN
         IF (ANY(swtyp == ["wet","all"])) THEN
            SELF%dwet = 1.e-10
            ! Only ice will have vrime =/ 0
            SELF%dwet = ( (SUM(SELF%volc(1:nwet))+SELF%vrime)/SELF%numc/pi6 )**(1./3.)
         END IF

         IF (ANY(swtyp == ["dry","all"])) THEN
            SELF%ddry = 1.e-10
            SELF%ddry = (SUM(SELF%volc(1:ndry))/SELF%numc/pi6)**(1./3.)
         END IF

         IF (ANY(swtyp == ["ins","all"]) .AND. ALL( spec%ind_insoluble(:) > 0 )) THEN ! If there is any insoluble active, all the indices sohuld be > 0
            SELF%dins = 1.e-10
            SELF%dins = ( SUM(SELF%volc(spec%ind_insoluble))/SELF%numc/pi6 )**(1./3.)
         END IF
      ELSE
         IF (ANY(swtyp == ["wet","all"])) SELF%dwet = SELF%dmid

         IF (ANY(swtyp == ["dry","all"])) SELF%ddry = SELF%dmid

         IF (ANY(swtyp == ["ins","all"])) SELF%dins = SELF%dmid
      END IF

      IF (limit) THEN
         IF (ANY(swtyp == ["wet","all"])) &
              SELF%dwet = MIN(SELF%dwet,SELF%dlim)
         IF (ANY(swtyp == ["dry","all"])) &
              SELF%ddry = MIN(SELF%ddry,SELF%dlim)
         IF (ANY(swtyp == ["ins","all"])) &
              SELF%dins = MIN(SELF%dins,SELF%dlim)
      END IF

    END SUBROUTINE updateDiameter

    ! 
    ! Subroutine updateRhomean
    !
    ! updateRhomean just gets the bulk mass weighted mean ice density
    ! for partially rimed particles. Another subroutine should be added
    ! where the particle shape is also taken into account for getting
    ! the true effective density (once implemented, the effective
    ! density is the one that should be used for fall velocities etc)
    !
    ! -------------------------------------------------
    !
    SUBROUTINE updateRhomean(SELF)
      CLASS(Section), INTENT(inout) :: SELF
      REAL :: mass_p, mass_r, mass_t

      INTEGER :: iwa
      
      iwa = spec%getIndex("H2O")

      SELF%rhomean = spec%rhoic
      
      ! convert to masses -> get the mass mean density - this is needed for ice, for others it's always rhowa
      IF (SELF%phase > 3) THEN
         IF (SELF%numc > SELF%nlim) THEN
            mass_p = spec%rhoic*SELF%volc(iwa)
            mass_r = spec%rhori*SELF%vrime
            mass_t = mass_p + mass_r
            SELF%rhomean = (mass_p*spec%rhoic + mass_r*spec%rhori)/mass_t
         ELSE
            SELF%rhomean = spec%rhoic
         END IF
      ELSE
         SELF%rhomean = spec%rhowa
      END IF

    END SUBROUTINE updateRhomean

    !
    ! Subroutine addRimedVol
    !
    ! Addition of rimed ice volume to a particle e.g. from ice-precip collection
    ! is not straightforward because of it's effect on the bulk ice mean density.
    ! The contributions must first be converted into mass and go from there...
    !
    !
    !SUBROUTINE addRimedVol(SELF, addvol)
    !  CLASS(Section)

    !END SUBROUTINE addRimedVol




END MODULE classSection
