MODULE emission_init
  USE mo_submctl, ONLY : spec, nbins, pi6, in1a, fn2a, in2b, fn2b
  USE emission_types, ONLY : emitConfig, emitModes, emitData, emitType3, nEmissionModes
  USE mo_lagrangian_tracker, ONLY : lagrangian_tracker
  USE mo_salsa_types, ONLY : aero
  USE mo_aux_state, ONLY : aetot,xt,yt,zt
  USE mo_salsa_sizedist, ONLY : size_distribution
  USE util, ONLY : smaller, getMassIndex, closest, arr_resize
  USE exceptionHandling, ONLY : errorMessage
  USE grid, ONLY : deltax, deltay,nxp,nyp
  USE mpi_interface, ONLY : myid
  IMPLICIT NONE

  CHARACTER(len=100), PARAMETER :: global_name = "emission_init"
  
  CONTAINS
    
    SUBROUTINE init_emission()
      IMPLICIT NONE
      CHARACTER(len=50), PARAMETER :: name = "init_emission"
      
      INTEGER :: ibin
      INTEGER :: st,en
      INTEGER :: nprof
      INTEGER :: mb1,mb12,mb2,mb22,nc1,nc2
      REAL :: core(nbins), naero(1,1,nbins)
      REAl, ALLOCATABLE ::  x(:), y(:), z(:)
      
      DO nprof = 1,nEmissionModes
         
         ASSOCIATE(emd => emitModes(nprof), edt => emitData(nprof))
           
           ! Use pointers for individual emission configuration and data instances to clean up the code
           IF (myid == 0) THEN
              WRITE(*,*) ''
              WRITE(*,*) '------------------------------------'
              WRITE(*,*) TRIM(name)//': INITIALIZING EMISSION PROFILE NO ',nprof
           END IF
           
           ! The aerosol species specified for emission MUST also be active in SALSA configuration
           IF ( .NOT. spec%isUsed(emd%species) ) THEN
              CALL errorMessage(global_name, name, &
                   'Attempt to emit <'//TRIM(emd%species)//'> but the compound '// &
                   'is not set to be used in the SALSA namelist.')
              STOP
           END IF
           
           ! Initialize the emission size distribution arrays
           ALLOCATE( edt%numc(nbins), edt%mass(nbins*spec%getNSpec(type="wet")) )
           edt%numc(:) = 0.; edt%mass(:) = 0.
           
           IF (emd%emitType == 1) THEN       ! Parameterized sea salt emissions
              IF (myid == 0) THEN   
                 WRITE(*,*) 'SEA SALT EMISSIONS; NO INITIALIZATION IMPLEMENTED (NOR NEEDED?)'
              END IF
              
              edt%numc(:) = 0.; edt%mass(:) = 0. ! i.e. do nothing
              
           ELSE IF (emd%emitType >= 2) THEN 
              
              ! Index limits for the regimes
              CALL regime_limits(emd%regime,st,en)
              
              IF (emd%emitSizeDistType == 1 ) THEN ! Monochromatic aerosol
                 
                 ! Determine the destination bin of the monochromatic emission from the
                 ! specified diameter
                 ibin = smaller(aetot%d(st:en),emd%emitDiam) - 1
                 
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
                 naero = 0.
                 CALL size_distribution( 1,1,1,1,st,en,   &
                                         [emd%emitNum],   &
                                         [emd%emitDiam],  &
                                         [emd%emitSigma], &
                                         naero            )
                 
                 ! -- Get the indices for emitted species and water.
                 CALL bin_indices( nprof,nbins,st,en,  &
                                   nc1, nc2, mb1, mb2, &
                                   mb12, mb22          )
                 
                 ! Set the emission number and mass concentrations to the emission data instance
                 edt%numc(st:en) = naero(1,1,st:en)
                 edt%mass(mb1:mb12) = naero(1,1,st:en)*core(st:en)*spec%rholiq(nc1)
                 ! Small amount of water
                 edt%mass(mb2:mb22) = 0.001 * edt%mass(mb1:mb12)
                 
              END IF
              
              ! Set up the emission levels. This will preferentially use the level indices, except if they're not given by the namelist
              !CALL init_emission_heights(emd)
              
           !END IF
           
             IF (ANY(emd%emitType == [3,5]))  THEN
                ASSOCIATE (emdT3 => emitType3(nprof))
                CALL init_emitType3_map(emd%emitMap,x,y,z,emdT3%np)
                CALL lagrangian_tracker(emdT3%ix,emdT3%iy,emdT3%iz,emdT3%t,emdT3%np, &
                                        emdT3%t_in,emdT3%t_out,emdT3%t_trac,         &
                                        emd%scS,emd%start_time,emd%end_time,x,y,z    )
              END ASSOCIATE           
             ELSE              
              ! Set up the emission levels. This will preferentially use the level indices, except if they're not given by the namelist
                CALL init_emission_heights(emd)
             END IF
          END IF 
           
         END ASSOCIATE
         
      END DO
      
      IF (myid == 0) THEN
         WRITE(*,*) '-------------------------------------'
         WRITE(*,*) ''
      END IF
      
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
      
      
      ASSOCIATE( emitLevMin => emd%emitLevMin,                        &
                 emitLevMax => emd%emitLevMax,                        &
                 emitHeightMin => emd%emitHeightMin,                  &
                 emitHeightMax => emd%emitHeightMax                   )
        
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
           
           ilev = closest(zt%d,emitHeightMin)
           emitLevMin = ilev
           ilev = closest(zt%d,emitHeightMax)
           emitLevMax = ilev
           
        ELSE IF (emitLevMin > 0 .OR. emitLevMax > 0) THEN
           
           maxi = MAX(emitLevMin, emitLevMax)
           IF (emitLevMin == -999) emitLevMin = maxi
           IF (emitLevMax == -999) emitLevMax = maxi
           
           ! Update the emission height levels according to the indices
           emitHeightMin = zt%d(emitLevMin)
           emitHeightMax = zt%d(emitLevMax)
           
        END IF

      END ASSOCIATE
      
  END SUBROUTINE init_emission_heights

  ! -----------------------------------------------------------------------------
  ! Subroutine init_emitType3_map: reads emitMap and calls lagrangian_tracker to return 
  ! cell indices and residence times of the emission source in the cells 
  SUBROUTINE init_emitType3_map(emitMap,x,y,z,np)
    CHARACTER(len=50), PARAMETER :: name = "init_emitType3_map"
    
    CHARACTER(len=40), INTENT(in) :: emitMap
    REAl, ALLOCATABLE, INTENT(out) ::  x(:), y(:), z(:)
    INTEGER, INTENT(out) :: np 
    
    REAl    :: deltax2, deltay2, xlim1, xlim2, ylim1, ylim2
    INTEGER :: i
    REAl    :: eps = 2e-2
    
    
    IF (myid == 0) THEN
       WRITE(*,*) ''
       WRITE(*,*) '------------------------------------'
       WRITE(*,*) TRIM(name)//': READING EMISSION MAP FROM ', TRIM(emitMap)
    END IF
    
    np = 8
    OPEN(1001,file=emitMap,status='old',form='formatted')
    
    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(y)) DEALLOCATE(y)
    IF (ALLOCATED(z)) DEALLOCATE(z)
    
    ALLOCATE (x(np), y(np), z(np))
        
    deltax2 = deltax/2
    deltay2 = deltay/2
    
    xlim1 = MINVAL(xt%d) + 3*deltax2
    xlim2 = MAXVAL(xt%d) - 3*deltax2
    ylim1 = MINVAL(yt%d) + 3*deltay2
    ylim2 = MAXVAL(yt%d) - 3*deltay2
    
    i = 1
    DO WHILE (.TRUE.)
       READ(1001,*,end=101) x(i), y(i), z(i)
       IF (x(i) == xlim1) x(i) = xlim1 - eps
       IF (x(i) == xlim2) x(i) = xlim2 - eps
       IF (y(i) == ylim1) y(i) = ylim1 + eps
       IF (y(i) == ylim2) y(i) = ylim2 + eps
       
       IF ( i >= np) THEN
          np = np*2
          CALL arr_resize(x,np)
          CALL arr_resize(y,np)
          CALL arr_resize(z,np)
       END IF
       i = i + 1
    END DO
    
101 CONTINUE
    CLOSE(1001)
    
    np = i - 1
    CALL arr_resize(x,np)
    CALL arr_resize(y,np)
    CALL arr_resize(z,np)
    
  END SUBROUTINE init_emitType3_map

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
    
    IF (myid == 0) THEN
      WRITE(*,*) name,'---------------------------------'
      WRITE(*,*) nc_emit, nc_h2o
      WRITE(*,*) '--------------------------------------'
    END IF
 
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
  
END MODULE emission_init
