MODULE emission_main
  USE mo_seasalt_emission

  USE mo_submctl, ONLY : pi6, in1a, fn2a, in2b, fn2b, nbins, aerobins, spec, pi6
  USE mo_salsa_types, ONLY : aero

  USE mo_salsa_sizedist, ONLY : size_distribution  ! Could this be packaged somehow differently?

  USE grid, ONLY: deltax, deltay, dzt, zt,  & ! Note dzt is inverse of the distance
                  nxp,nyp,nzp,              &
                  a_up, a_vp, a_dn,      &
                  a_maerot, a_naerot
  USE util, ONLY: smaller, closest, getMassIndex
  USE exceptionHandling, ONLY: errorMessage

  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "emission_main"
  
  !
  ! ----------------------------------------------
  ! Contains arrays for emitted size distribution
  !
  TYPE EmitSizeDist
     REAL, ALLOCATABLE :: numc(:)
     REAL, ALLOCATABLE :: mass(:)
  END TYPE EmitSizeDist

  !
  ! --------------------------------------------
  ! Type for emission configuration options
  !
  TYPE  EmitConfig
     INTEGER          :: emitType = 1                 ! 1: Natural seasalt emissions, 2: custom artificial emissions 
     INTEGER          :: regime = 1                   ! Destination bin regime for emitted aerosol. 1: A, 2: B
     REAL             :: start_time = 0.,  &          ! Start time for emission (s)
                         end_time = 86400.            ! End time for emission (s)

     ! Parameters below valid for emitType > 1
     CHARACTER(len=3) :: species = 'SS '              ! Which aerosol species to emit (must conform to names in src_salsa/classSpecies)
     REAL             :: emitHeightMin = -999.,  &    ! Min height of airborne emissions (m)
                         emitHeightMax = -999.        ! Max height (m)
     INTEGER          :: emitLevMin = -999            ! Integer levels corresponding to the heights; If both heights and lev indices are specified,
     INTEGER          :: emitLevMax = -999            ! the level indices are preferred (so use one or the other...).

     INTEGER          :: emitSizeDistType = 1         ! 1: Monochromatic aerosol, 2: modal size disribution (lognormal)
     REAL             :: emitDiam = 10.e-6,    &      ! Assumed (dry )diameter of the particles (mode diameter for emitType=2).
                         emitNum  = 10000.            ! Number consentration of particles emitted per second #/m3/s (mode concentration for emitType=2)
     REAL             :: emitSigma = 2.0              ! Geometric standard deviation for emitSizeDist=2
  END TYPE EmitConfig
  
  INTEGER :: nEmissionModes = 0                             ! Number of emission modes (NAMELIST variable, max=5)
  INTEGER, PARAMETER :: maxEmissionModes = 5                ! Max number of emission modes
  TYPE(EmitConfig), TARGET :: emitModes(MaxEmissionModes)   ! Configuration instances for emission modes
  TYPE(EmitSizeDist), TARGET :: emitData(MaxEmissionModes)  ! Emission size distribution data for each setup/mode.

  ! The limit of 5 emission modes can be easily increased if needed. The maximum has to be hardcoded here
  ! as a PARAMETER in order to be able to give the configuration options directly from NAMELIST (the instances must be allocated).  

  CONTAINS

  SUBROUTINE init_emission
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "init_emission"

    INTEGER :: ibin
    INTEGER :: st,en
    INTEGER :: nprof
    INTEGER :: mb1,mb12,mb2,mb22,nc1,nc2
    REAL :: core(nbins), naero(1,1,nbins)

    DO nprof = 1,nEmissionModes

       ASSOCIATE(emd => emitModes(nprof), edt => emitData(nprof))

         WRITE(*,*) ''
         WRITE(*,*) '------------------------------------'
         WRITE(*,*) TRIM(name)//': INITIALIZING EMISSION PROFILE NO ',nprof
         
         ! The aerosol species specified for emission MUST also be active in SALSA configuration
         IF ( .NOT. spec%isUsed(emd%species) ) THEN
            CALL errorMessage(global_name, name, &
                 'Attempt to emit <'//TRIM(emd%species)//'> but the compound '// &
                 'is not set to be used in the SALSA namelist.')
            STOP
         END IF
         
         ! Initialize the emission size distribution arrays
         ALLOCATE( edt%numc(nbins), edt%mass(nbins*spec%getNSpec()) )
         edt%numc(:) = 0.; edt%mass(:) = 0.
         
         IF (emd%emitType == 1) THEN       ! Parameterized sea salt emissions
            
            WRITE(*,*) 'SEA SALT EMISSIONS; NO INITIALIZATION IMPLEMENTED (NOR NEEDED?)'
            edt%numc(:) = 0.; edt%mass(:) = 0. ! i.e. do nothing
            
         ELSE IF (emd%emitType == 2) THEN  ! Artificial emission
            
            ! Index limits for the regimes
            CALL regime_limits(emd%regime,st,en)
            
            IF (emd%emitSizeDistType == 1 ) THEN ! Monochromatic aerosol
               
               ! Determine the destination bin of the monochromatic emission from the
               ! specified diameter
               ibin = smaller(aerobins(st:en),emd%emitDiam) - 1
               
               ! Get bin indices for emitted mass tendencies (in addition to specified species, put also
               ! some water due to current limitations in condensation in mo_salsa_dynamics. Review this later!)
               CALL bin_indices(nprof,nbins,st,en,   &
                                nc1, nc2, mb1, mb2,  &
                                mb12, mb22, ibin=ibin)
               
               ! Place the emission number and mass concentrations to the emission data instance
               edt%numc(st+ibin) = emd%emitNum
               edt%mass(mb1) = emd%emitNum *   &
                    (pi6*emd%emitDiam**3)*spec%rholiq(nc1)
               ! Small amount of water due the condensation issue
               edt%mass(mb2) = 0.001 * edt%mass(mb1)
               
            ELSE IF (emd%emitSizeDistType == 2 ) THEN! Lognormal mode
               
               ! Consider specified emitNum and emitDiam as lognormal mode values. Also emitSigma needs to be specified here
               
               ! First, get the core volume of single particle in each aerosol bin based on the bin mean diameter
               core(:) = pi6 * aero(1,1,:)%dmid**3
               
               ! Get the binned size distribution (brackets used because size_distribution expects arrays for mode values)
               CALL size_distribution(1,1,1,1,[emd%emitNum],    &
                                              [emd%emitDiam],   &
                                              [emd%emitSigma],  &
                                               naero)
               !naero = 1.

               ! -- Get the indices for emitted species and water.
               CALL bin_indices(nprof,nbins,st,en,  &
                                nc1, nc2, mb1, mb2, &
                                mb12, mb22          )
               
               ! Set the emission number and mass concentrations to the emission data instance
               edt%numc(st:en) = naero(1,1,st:en)
               edt%mass(mb1:mb12) = naero(1,1,st:en)*core(st:en)*spec%rholiq(nc1)
               ! Small amount of water
               edt%mass(mb2:mb22) = 0.001 * edt%mass(mb1:mb12)
               
            END IF
            
         END IF
       
         ! Set up the emission levels. This will preferentially use the level indices, except if they're not given by the namelist
         CALL init_emission_heights(emd)
       
       END ASSOCIATE

    END DO

    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) ''

  END SUBROUTINE init_emission

  !
  ! ----------------------------------------------------------------------------
  ! Subroutine init_emission_heights: Sorts out the level and height information
  !                                   in the emission configuration
  !
  SUBROUTINE init_emission_heights(emd)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "init_emission_heights"

    TYPE(EmitConfig), INTENT(inout) :: emd
    
    INTEGER :: ilev
    INTEGER :: maxi
    REAL    :: maxf


    ASSOCIATE( emitLevMin => emd%emitLevMin,       &
               emitLevMax => emd%emitLevMax,       &
               emitHeightMin => emd%emitHeightMin, &
               emitHeightMax => emd%emitHeightMax  )
    
      ! At least one of the height or level definitions in the emitConfig must be specified
      IF ( ALL( [emitLevMin, emitLevMax] == -999 ) .AND.  &
           ALL( [emitHeightMin, emitHeightMax] == -999. ) ) THEN
      
          CALL errorMessage(global_name, name, &
               "At least of the the height or level definitions must be"// &
               " specified in emitConfig")
          STOP

      END IF

      IF ( emitLevMin == -999 .AND. emitLevMax == -999) THEN 

         ! EmitConfig must have emitHeightMin or emitHeightMax specified as a positive value
         maxf = MAX(emitHeightMin, emitHeightMax)
         IF (emitHeightMin == -999.) emitHeightMin = maxf
         IF (emitHeightMax == -999.) emitHeightMax = maxf

         ilev = closest(zt,emitHeightMin)
         emitLevMin = ilev
         ilev = closest(zt,emitHeightMax)
         emitLevMax = ilev
         
      ELSE IF (emitLevMin > 0 .OR. emitLevMax > 0) THEN
         
         maxi = MAX(emitLevMin, emitLevMax)
         IF (emitLevMin == -999) emitLevMin = maxi
         IF (emitLevMax == -999) emitLevMax = maxi
         
         ! Update the emission height levels according to the indices
         emitHeightMin = zt(emitLevMin)
         emitHeightMax = zt(emitLevMax)
         
      END IF
      
    END ASSOCIATE

  END SUBROUTINE init_emission_heights

  !
  ! --------------------------------------------------------------------------
  ! Subroutine regime_limits: Returns the limiting indices for A or B regimes
  !                           (reg=1 or reg=2, respectively)
  !
  SUBROUTINE regime_limits(reg,st,en)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "regime_limits"

    INTEGER, INTENT(in) :: reg
    INTEGER, INTENT(out) :: st,en

    IF ( reg == 1 ) THEN
       st = in1a
       en = fn2a
    ELSE IF ( reg == 2 ) THEN
       st = in2b
       en = fn2b
    END IF

  END SUBROUTINE regime_limits

  !
  ! -----------------------------------------------------------------------------
  ! Subroutine bin_indices: Returns bin and mass indices needed in init_emission
  !
  SUBROUTINE bin_indices(nprof,nbins,st,en,          &
                         nc_emit, nc_h2o,      &
                         mb_emit_1, mb_h2o_1,  &
                         mb_emit_2, mb_h2o_2, ibin)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "bin_indices"

    INTEGER, INTENT(in) :: nprof,nbins, st, en
    INTEGER, INTENT(in), OPTIONAL :: ibin

    INTEGER, INTENT(out) :: nc_emit, nc_h2o,      &  ! Mass indices for emitted species and water
                            mb_emit_1, mb_h2o_1,  &  ! First index of the regime for emission size distribution 
                                                     ! (if ibin present, this will be the single emission bin for monochromatic emissions)
                            mb_emit_2, mb_h2o_2      ! Last index of the regime for emission size distribution
                                                     ! (will be put out, but not used for monochromatic emissions)

    nc_emit  = spec%getIndex(emitModes(nprof)%species)
    nc_h2o  = spec%getIndex("H2O")

    WRITE(*,*) name,'---------------------------------'
    WRITE(*,*) nc_emit, nc_h2o
    WRITE(*,*) '--------------------------------------'

    IF ( PRESENT(ibin) ) THEN
       mb_emit_1 = getMassIndex(nbins,st+ibin,nc_emit)
       mb_h2o_1 = getMassIndex(nbins,st+ibin,nc_h2o)       
    ELSE
       mb_emit_1  = getMassIndex(nbins,st,nc_emit)
       mb_h2o_1  = getMassIndex(nbins,st,nc_h2o)
    END IF

    mb_emit_2 = getMassIndex(nbins,en,nc_emit)
    mb_h2o_2 = getMassIndex(nbins,en,nc_h2o)

  END SUBROUTINE bin_indices

  !
  ! -------------------------------------------------------------------
  ! subroutine aerosol_emission:  calls methods to calculate emitted
  !                               aerosols from ground/sea
  !  
  ! Adapted from the original code by Antti Kukkurainen
  ! Juha Tonttila, FMI, 2017
  !
  SUBROUTINE aerosol_emission(time_in)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "aerosol_emission"

    REAL, INTENT(in) :: time_in   ! time in seconds
    LOGICAL :: condition
    INTEGER :: pr

    ! Loop over all specified emission profiles
    DO pr = 1,nEmissionModes
       ASSOCIATE(emd => emitModes(pr), edt => emitData(pr))
         condition = getCondition(emd,time_in)
         IF (condition .AND. emd%emitType == 2) CALL custom_emission(edt,emd)
       END ASSOCIATE
    END DO

  END SUBROUTINE aerosol_emission

  !
  ! ---------------------------------------------------------------
  ! Simulates the emission of seasalt particles from
  ! an ocean surface as a function of the 10-m wind
  ! speed.
  !
  !SUBROUTINE surface_emission()
  !  IMPLICIT NONE

 !   CHARACTER(len=50), PARAMETER :: name = "surface_emission"
 !   
 !   REAL :: mass_flux(1,nbins) !mass flux at given radius
 !   REAL :: numb_flux(1,nbins)        !number flux at given radius
 !   REAL :: pseaice(1) = 0            !sea ice fraction
 !   REAL :: velo10m_salsa(1,1)        !wind speed

!    INTEGER :: nc, st, en, ii, jj
!    INTEGER :: in, fn
    
!    ! Surface seasalt emissions possible only if sea salt aerosol is used
!    IF (spec%isUsed(spec%nss)) THEN
!       nc=spec%getIndex(spec%nss)
       
!       IF (esrfc%regime == 1) THEN
!          in = in1a
!          fn = fn2a
!       ELSE IF (esrfc%regime == 2) THEN
!          in = in2b
!          fn = fn2b
!       END IF
!       st=(nc-1)*nbins+in
!       en=(nc-1)*nbins+fn
!       
!       DO jj=3,nyp-2
!          DO ii=3,nxp-2
!             
!             velo10m_salsa = SQRT(a_up(2,ii,jj)**2 + a_vp(2,ii,jj)**2)
!             
!             CALL seasalt_emissions_lsce_salsa(1, 1, 1, pseaice, velo10m_salsa, mass_flux, numb_flux)
!             
!             !number of particles + more particles per unit of time * scaling factor [#/kg]
!             a_naerot(2,ii,jj,in:fn) = a_naerot(2,ii,jj,in:fn) + numb_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             !mass + more mass per unit of time * scaling factor [kg/kg]
!             a_maerot(2,ii,jj,st:en) = a_maerot(2,ii,jj,st:en) + mass_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             
!          END DO
!       END DO
!       
!    END IF
!    
!  END SUBROUTINE surface_emission

  !
  ! ----------------------------------------------------------------
  ! Subroutine custom_emissio: "Customized" emission routine, mainly
  !                            for simulating atnhropogenic emissions,
  !                            such as ship or aircraft emissions etc.
  !                            Support for point sources will be included
  !                            soon. Now only does domain-wide emissions
  !                            at specified altitude and time.
  !
  SUBROUTINE custom_emission(edt,emd)
    IMPLICIT NONE
    
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding"
    
    TYPE(EmitSizeDist), INTENT(in) :: edt  ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd   ! Emission configuration instance

    INTEGER :: i,j,bb,ss,mm

    WRITE(*,*) '========================'
    WRITE(*,*) 'CALCULATING EMISSIONS'
    WRITE(*,*) '========================'
    

    ASSOCIATE( k1 => emd%emitLevMin, k2 => emd%emitLevMax )
    
      DO bb = 1,nbins
         DO j = 1,nyp
            DO i = 1,nxp
               a_naerot(k1:k2,i,j,bb) = a_naerot(k1:k2,i,j,bb) + edt%numc(bb)
               DO ss = 1,spec%getNSpec()
                  mm = getMassIndex(nbins,bb,ss)
                  a_maerot(k1:k2,i,j,mm) = a_maerot(k1:k2,i,j,mm) + edt%mass(mm)
               END DO
            END DO
         END DO
      END DO
      
    END ASSOCIATE

  END SUBROUTINE custom_emission
  
  ! ----------------------------------------------------------

  FUNCTION getCondition(emd,time)
    IMPLICIT NONE
    LOGICAL :: getCondition
    CLASS(EmitConfig), INTENT(in) :: emd
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "getCondition"

    getCondition = (                                &
                    emd%start_time <= time    .AND. &
                    emd%end_time > time            &
                   )   

  END FUNCTION getCondition
  


END MODULE emission_main
