!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
! 
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
! 
PROGRAM ucla_les

   USE mpi_interface, ONLY : myid

   IMPLICIT NONE

   REAL :: t1, t2
   
   CALL cpu_time(t1)
   CALL driver
   CALL cpu_time(t2)

   IF (myid == 0) THEN
      PRINT "(/,' ',49('-')/,' ',A16,F10.1,' s')", '  Execution time: ', t2-t1
      STOP ' ..... Normal termination'
   END IF

CONTAINS

   !----------------------------------------------------------------------
   ! Subroutine Driver:  This is the main program driver.  It calls routines
   ! to read the model initialization file, and configure memory and pointes.
   ! It also calls the routines which initialize the model and timestep it.
   !
   SUBROUTINE driver

      USE grid, ONLY          : define_grid, define_vars, level, nxp, nyp, nzp, nxpart
      USE init, ONLY          : initialize
      USE mo_field_init, ONLY : initialize_FieldArrays
      USE step, ONLY          : stepper
      USE mpi_interface, ONLY : init_mpi, define_decomp,                    &
                                init_alltoall_reorder, appl_finalize
      
      ! Added for SALSA
      USE mo_salsa_init, ONLY : define_salsa, salsa_initialize

      IMPLICIT NONE

      INTEGER :: ierror
      
      CALL init_mpi

      CALL define_parm

      IF (level >= 4) CALL define_salsa(level) ! Read SALSA namelist etc.

      IF (level >= 4) CALL salsa_initialize() ! All salsa variables and bin diameters are now initialized
      
      CALL define_decomp(nxp, nyp, nxpart)

      CALL define_grid

      CALL init_alltoall_reorder(nxp, nyp, nzp)

      CALL define_vars

      CALL initialize_FieldArrays ! the source scalar arrays must be allocated at this point (define_vars)

      CALL initialize ! Added initialization of aerosol size distributions here + a single call
                      ! for SALSA to set up cloud microphysics

      CALL stepper

      CALL appl_finalize(ierror)

      RETURN
   END SUBROUTINE driver

   !
   ! ----------------------------------------------------------------------
   ! Subroutine Read_nl: Driver for reading model namelist
   !
   SUBROUTINE define_parm

    USE util, ONLY              : fftinix,fftiniy
    USE sgsm, ONLY              : csx, prndtl
    USE srfc, ONLY              : isfctyp, zrough, ubmin, dthcon, drtcon, C_heat,          &
                                  deepSoilTemp, lConstSoilWater, lConstSoilHeatCap
    USE step, ONLY              : timmax, istpfl, corflg, outflg, frqhis,                    &
                                  strtim, radfrq
    USE nudg_defs, ONLY         : nudge_time, nudge_zmin, nudge_zmax,                                &
                                  ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero
    USE emission_types, ONLY    : emitModes, nEmissionModes
    USE grid, ONLY              : deltaz, deltay, deltax, nzp, nyp, nxp, nxpart,                     &
                                  dtlong, dzrat,dzmax, th00, umean, vmean, naddsc, level,            &
                                  filprf, expnme, isgstyp, igrdtyp, iradtyp, lnudging, lemission,    &
                                  lpback, pbncsrc, nfpt, distim, runtype, CCN,sst,W1,W2,W3, &
                                  cntlat, varlist_main, varlist_ps, varlist_ts
    USE init, ONLY              : us, vs, ts, rts, ps, hs, ipsflg, itsflg,iseed, hfilin,             &
                                  zrand, zrndamp, init_type
    USE init_warm_bubble, ONLY  : bubble_center, bubble_diameter, bubble_temp_ampl
    USE forc, ONLY              : div, case_name     ! Divergence, forcing case name
    USE radiation_main, ONLY    : radsounding,   &
                                  sfc_albedo,    &
                                  useMcICA,      &
                                  laerorad,      &
                                  RadConstPress, &
                                  RadPrecipBins
    USE mcrp, ONLY              : sed_aero, sed_cloud, sed_precp, sed_ice, init_mcrp_switches, &
                                  bulk_autoc, bulkScheme
    USE mpi_interface, ONLY     : myid, appl_abort, ver, author
    USE mo_output, ONLY         : ts_intvl, ps_intvl, main_intvl
    USE mo_check_state, ONLY    : breakUndefOutput
    USE mo_stats_parameters, ONLY : TH_rc, TH_rr, TH_ri, TH_rrate
    
    IMPLICIT NONE
    
    NAMELIST /model/     &
         expnme    ,       & ! experiment name
         nxpart    ,       & ! whether partition in x direction?
         naddsc    ,       & ! Number of additional scalars
         corflg , cntlat , & ! coriolis flag
         nfpt   , distim , & ! rayleigh friction points, dissipation time
         level  , CCN    , & ! Microphysical model Number of CCN per kg of air
         nxp    , nyp    , nzp   ,  & ! number of x, y, z points
         deltax , deltay , deltaz , & ! delta x, y, z (meters)
         dzrat  , dzmax  , igrdtyp, & ! stretched grid parameters
         timmax , dtlong , istpfl , & ! timestep control
         runtype, hfilin , filprf , & ! type of run (INITIAL or HISTORY)
         frqhis , outflg , &          ! freq of history writes, output flg
         strtim ,                   & ! Model start time
         iradtyp,                   & ! Radiation type
         isgstyp, csx    , prndtl , & ! SGS model type, parameters
         lnudging, lemission,       & ! master switch for nudging, aerosol emissions
         lpback, pbncsrc,       & ! Switch for piggybacking microphysics, switch for fixed or SALSA CDNC with PB
         div, case_name, &            ! divergence for LEVEL 4
         sed_aero, sed_cloud, sed_precp, sed_ice,  & ! Sedimentation (T/F)
         bulk_autoc,                & ! autoconversion (and accretion) switch for level < 4 
         bulkScheme                   ! 
    
    NAMELIST /initialization/      &
         init_type,                & ! Type of initialization: 1: random perturbations, 2: warm bubble
         bubble_center,            & ! Center coordinates for warm bubble (z,x,y) 
         bubble_diameter,          & ! Diameter for warm bubble (z,x,y)   
         bubble_temp_ampl,         & ! Temperature amplitude for the warm bubble, assume sinusoidal
         ipsflg, itsflg,           & ! sounding flags
         hs, ps, ts,               & ! sounding heights, pressure, temperature
         us, vs, rts,              & ! sounding E/W winds, water vapor
         umean, vmean, th00,       & ! gallilean E/W wind, basic state
         iseed, zrand, zrndamp       ! random seed
    
    NAMELIST /radiation/           &
         radfrq,                   & ! radiation type flag RADFRQ NOT USED ANYWHERE, VARIABLE DECLARED IN STEP.F90
         radsounding, sfc_albedo,  & ! Name of the radiation sounding file, surface albedo
         useMcICA,                 & ! Use the Monte Carlo Independent Column Approximation method (T/F)
         laerorad,                 & ! Use the binned aerosol data for radiation (with SALSA)
         RadConstPress,            & ! keep constant pressure levels (T/F) 
         RadPrecipBins               ! add precipitation bins cloud water (0, 1, 2, 3,...)
    
    NAMELIST /nudge/   &
         nudge_time,                       & ! Total nudging time (independent of spin-up)
         nudge_zmin, nudge_zmax,           & ! Altitude (m) range for nudging
         ndg_theta,                        & ! Temperature nudging
         ndg_rv,                           & ! Water vapor mixing ratio nudging
         ndg_u,                            & ! wind nudging
         ndg_v,                            & ! Horizontal wind nudging
         ndg_aero                             ! Aerosol number concentration nudging
    
    NAMELIST /emission/ &              
         nEmissionModes,      & ! Number of emission profiles to be used (max 5)
         emitModes              ! Emission configs
    
    NAMELIST /surface/      &
         isfctyp,           &   ! Surface parameterization type
         zrough,            &   ! roughness length
         ubmin,             &   !
         sst,               &   ! Surface temperature
         dthcon,            &   ! Sensible heat flux
         drtcon,            &   ! Latent heat flux 
         C_heat,            &   ! Soil heat capacity (Only with isfctyp=5)
         deepSoilTemp,      &   ! Deep soil temperature (Only with isfctyp=5)
         W1,                &   ! Soil water contents
         W2,                &  
         W3,                &
         lConstSoilWater,   &   ! Keep soil water content(s) constant? (Only with isfctyp=5)
         lConstSoilHeatCap      ! Keep soil heat capacity con
    
    NAMELIST /output/       &
         breakUndefOutput,  &   ! Stop the model if undefined output variables are entered to varlists (T/F)
         ts_intvl,          &   ! Output interval for TS statistics (s)
         ps_intvl,          &   ! Output interval for PS statistics (s)
         main_intvl,        &   ! Output interval for the main analysis files (s)
         varlist_main,      &   ! List of variables for the main analysis output
         varlist_ps,        &   ! List of variables for the PS statistics output
         varlist_ts,        &   ! List of variables for the TS statistics output
         TH_rc,             &   ! Threshold for cloud water mix rat for conditional averaging
         TH_rr,             &   ! Threshold for drizzle/rain mix rat for conditional averaging
         TH_ri,             &   ! Threshold for ice mix rat for conditional averaging
         TH_rrate               ! Threshold for precip rate for conditional averaging
         
    NAMELIST /version/  &
         ver, author        ! Information about UCLALES-SALSA version and author
    
    ps = 0.
    ts = th00
    !
    ! these are for initializing the temp variables used in ffts in x and y
    ! directions.
    !
    fftinix = 1
    fftiniy = 1

    ! Initialize output lists
    varlist_main = ''
    varlist_ps = ''
    varlist_ts = ''
    
    !
    ! Initialize some process switches (mcrp...etc). Need to be done here, before reading the NAMELIST!
    ! -------------------------------------------------------------------------------------------------
    CALL init_mcrp_switches()
    
    
    !
    ! read namelist from specified file
    !
    OPEN  (1,status='old',file='NAMELIST')
    REWIND(1)
    READ (1, nml=model)
    REWIND(1)
    READ  (1, nml=initialization)
    REWIND(1)
    READ (1, nml=radiation)
    REWIND(1)
    READ (1, nml=surface)
    IF (lnudging) THEN
       REWIND(1)
       READ  (1, nml=nudge)
    ENDIF
    IF (lemission) THEN
       REWIND(1)
       READ  (1, nml=emission)
    ENDIF
    REWIND(1)
    READ(1, nml=output)
    REWIND(1)
    READ  (1, nml=version) 
    CLOSE(1)
    
    !
    ! write file variable control to standard output
    !
    IF (myid == 0) THEN
       IF (runtype == 'HISTORY') THEN
          WRITE(*,601) expnme, hfilin, timmax
       ELSE
          WRITE(*,600) expnme, timmax
       END IF
       IF (outflg) WRITE(*,602) filprf, frqhis, main_intvl
       !
       ! Do some cursory error checking in namelist variables
       !
       IF (laerorad .AND. level < 4) THEN 
          WRITE(*,*) "WARNING: laerorad=TRUE valid only for level >= 4, setting to FALSE"
          laerorad = .FALSE.
       END IF
       
       IF (MIN(nxp,nyp) < 5) THEN
          PRINT *, '  ABORTING: min(nxp,nyp) must be > 4.'
          CALL appl_abort(0)
       END IF
       
       IF (nzp < 3 ) THEN
          PRINT *, '  ABORTING: nzp must be > 2 '
          CALL appl_abort(0)
       END IF
       
       IF (cntlat < -90. .OR. cntlat > 90.) THEN
          PRINT *, '  ABORTING: central latitude out of bounds.'
          CALL appl_abort(0)
       END IF
       
       IF (lpback .AND. level /= 4) THEN
          WRITE(*,*) "ABORTING: lpback=TRUE currently valid only with level=4"
          CALL appl_abort(0)
       END IF
    END IF
    
600 FORMAT(//' ',49('-')/,' ',/,'  Initial Experiment: ',A50 &
         /,'  Final Time:         ',F8.1,' s'              )
601 FORMAT(//' ',49('-')/,' ',/,'  Restart Experiment: ',A50 &
         /,'  Restart File: ',A30,                           &
         /,'  Final Time: ',F10.1,' s'              )
602 FORMAT('  Output File Stem:   ',A50                        &
         /,'  History Frequency:  ',F7.1,                    &
         /,'  Analysis Frequency: ',F7.1,                    &
         /,'  Model spinup period: ',F7.1)
    
    RETURN
  END SUBROUTINE define_parm
  
END PROGRAM ucla_les

