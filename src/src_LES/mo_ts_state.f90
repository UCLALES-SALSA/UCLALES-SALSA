MODULE mo_ts_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_ts_procedures
  IMPLICIT NONE

  ! Note some of these have to be stored in memory in order to calculate them in the physics routines,
  ! others are associated with and onDemand procedure to calculate the values upon output. The onDemand
  ! subroutines should be implemented in mo_ts_procedures and they must have identical IO interface (*name,*output).

  ! There are 3 ways to calculate the statistics:
  ! 1. For stats variables associated with an onDemand procedure but no storage, assume the stats variable to
  !    have IDENTICAL NAME WITH THE CORRESPONDING FULL FIELD from where they are calculated (assuming there is one).
  !    This allows you to use automatic fetching of data using mo_stats_finder in the subroutines associated with onDemand. 
  ! 2. Stats variables associated with an onDemand procedure, but a UNIQUE name: An explicit paring of the stats variable
  !    and the full field source variable must be defined in the onDemand procedure in mo_ts_procedures, or in other ways
  !    determine the calculation of the stats variable by its name.
  ! 3. Stats variables associated with storage space and no onDemand procedure: calculate the values during the timestep
  !    flagged for statistics output (tsflg in mo_output.f90) within the appropriate physics subroutine.
  !
  ! Regarding 3), an accumulation scheme would be easy to implement, but for now is left for future.
  ! Note that for all cases 1)-3), MPI reduction operations should be applied before writing output, since the statistics IO
  ! WILL ONLY BE PERFORMED BY RANK=MPI_ROOT (from mpi_interface.f90).

  ! Adding new variables:
  ! Declare your new TYPE(FloatArray0d) variable. Add initialization and registration of the variable to the TYPE(FieldArray) :: TS
  ! in the subroutine setTSVariables found below. Just check how other variables are defined and do the same. Only thing to consider
  ! are the 3 points listed above.
  !   - If you want to use onDemand, first see if mo_ts_procedures already has a suitable subroutine to calculate your variable.
  !   - If so, and your variable fits point 1), just associate the onDemand variable in the initialization and you're good to go. Nothing else needed
  !     (except adding the variable name to output list in the NAMELIST)
  !   - If instead you see a suitable subroutine but your variable fits point 2), you must add a new paring for your variable and a
  !     source variable name in mo_ts_procedures, or otherwise direct the calculation explicitly by your stats variable name. Again
  !     nothing else needed, except the NAMELIST entry
  !   - Instead write your own subroutine in mo_ts_procedures and associate onDemand with that. Again nothing else needed
  !   - If you don't want to use onDemand, associate the variable with the storage instead. Now you need to calculate the value during
  !     a statistics timestep using the TYPE(FloatArray0d) variable you declared. Just make sure the end result is the global representative
  !     value across all processes and defined for rank mpi_root !!
  
  TYPE(FloatArray0d), TARGET :: vtke,              &  ! 1:
                                sfcbflux,          &  ! 2:
                                wmax,              &  ! 3:
                                tsrf,              &  ! 4:
                                shf_bar,lhf_bar,   &  ! 5,6:
                                lwp_bar,rwp_bar,   &  ! 7,8:
                                iwp_bar,           &  ! 9:
                                ctop,cbase,cfrac,  &  ! 10,11,12:
                                sfc_rain              ! 13:

  REAL, ALLOCATABLE, TARGET :: a_ts0d(:)
  INTEGER, PARAMETER        :: nts0d = 13

  
  CONTAINS

    SUBROUTINE setTSVariables(TS,outputlist,level)
      TYPE(FieldArray), INTENT(inout) :: TS
      CHARACTER(len=*), INTENT(in)    :: outputlist(:)
      INTEGER, INTENT(in)             :: level
      CLASS(*), POINTER :: pipeline => NULL()

      INTEGER :: nts

      ALLOCATE(a_ts0d(nts0d))
      a_ts0d = 0.
      nts = 0
      
      nts = nts + 1
      vtke = FloatArray0d(a_ts0d(nts))
      pipeline => vtke
      CALL TS%newField("vtke", "Vertical integral of total TKE", "kg/s", "time",   &  ! NOT IMPLEMENTED
                       ANY(outputlist == "vtke"), pipeline)

      nts = nts + 1
      sfcbflux = FloatArray0d(a_ts0d(nts))
      pipeline => sfcbflux
      CALL TS%newField("sfcbflux", "Surface buoyancy flux", "m/s2", "time",   &  ! NOT IMPLEMENTED
                       ANY(outputlist == "sfcbflux"), pipeline)

      wmax = FloatArray0d()
      wmax%onDemand => globalMax
      pipeline => wmax
      CALL TS%newField("wmax", "Maximum vertical velocity", "m/s", "time",    &
                       ANY(outputlist == "wmax"), pipeline)

      nts = nts + 1
      tsrf = FloatArray0d(a_ts0d(nts))
      pipeline => tsrf
      CALL TS%newField("tsrf", "Mean surface temperature", "K", "time",   &    ! NOT IMPLEMENTED
                       ANY(outputlist == "tsrf"), pipeline)

      nts = nts + 1
      shf_bar = FloatArray0d(a_ts0d(nts))
      pipeline => shf_bar
      CALL TS%newField("shf_bar", "Mean surface sensible heat flux", "W/m2", "time",   &  ! NOT IMPLEMENTED
                       ANY(outputlist == "shf_bar"), pipeline)

      nts = nts + 1
      lhf_bar = FloatArray0d(a_ts0d(nts))
      pipeline => lhf_bar
      CALL TS%newField("lhf_bar", "Mean surface latent heat flux", "W/m2", "time",   &  ! NOT IMPLEMENTED
                       ANY(outputlist == "lhf_bar"), pipeline)

      lwp_bar = FloatArray0d()
      lwp_bar%onDemand => globalAvgWaterPaths
      pipeline => lwp_bar
      CALL TS%newField("lwp_bar", "Mean liquid water path (cloud water)", "kg/m2", "time",   & 
                       ANY(outputlist == "lwp_bar"), pipeline)

      rwp_bar = FloatArray0d()
      rwp_bar%onDemand => globalAvgWaterPaths
      pipeline => rwp_bar
      CALL TS%newField("rwp_bar", "Mean precipitating water path", "kg/m2", "time",   &
                       ANY(outputlist == "rwp_bar"), pipeline)

      IF (level == 5) THEN
         iwp_bar = FloatArray0d()
         iwp_bar%onDemand => globalAvgWaterPaths
         pipeline => iwp_bar
         CALL TS%newField("iwp_bar", "Mean frozen water path", "kg/m2", "time",   &
                          ANY(outputlist == "iwp_bar"), pipeline)
      END IF
         
      ctop = FloatArray0d()
      ctop%onDemand => globalAvgCloudBoundaries
      pipeline => ctop
      CALL TS%newField("ctop", "Mean cloud top height", "m", "time",   &
                       ANY(outputlist == "ctop"), pipeline)

      cbase = FloatArray0d()
      cbase%onDemand => globalAvgCloudBoundaries
      pipeline => cbase
      CALL TS%newField("cbase", "Mean cloud base height", "m", "time",   &
                       ANY(outputlist == "cbase"), pipeline)

      cfrac = FloatArray0d()
      cfrac%onDemand => globalAvgCloudFraction
      pipeline => cfrac
      CALL TS%newField("cfrac", "Mean total cloud fraction", "m", "time",   &
                       ANY(outputlist == "cfrac"), pipeline)

      nts = nts + 1
      sfc_rain = FloatArray0d(a_ts0d(nts))
      pipeline => sfc_rain
      CALL TS%newField("sfc_rain", "Mean surface rain rate", "mm/hour", "time",   &   ! NOT IMPLEMENTED
                       ANY(outputlist == "sfc_rain"), pipeline)

    END SUBROUTINE setTSVariables
        
END MODULE mo_ts_state
