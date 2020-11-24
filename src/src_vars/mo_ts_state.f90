MODULE mo_ts_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_ts_procedures
  IMPLICIT NONE

  ! Contains variable definitions for timeseries statistics outputs.
  ! STATISTICAL OUTPUT IS ONLY WRITTEN BY THE ROOT PROCESS, I.E. RANK=MPI_ROOT. THEREFORE MAKE SURE TO USE
  ! MPI REDUCTION OPERATIONS TO DEFINE THE FINAL RESULT IN THE ROOT PROCESS.
  
  ! Note some of these have to be stored in memory in order to calculate them in the physics routines,
  ! others are associated with and onDemand procedure to calculate the values upon output. The onDemand
  ! subroutines should be implemented in mo_ts_procedures and they must have identical IO interface (*name,*output).

  ! There are 3 ways to calculate the statistics:
  ! 1. For stats variables associated with an onDemand procedure but no storage, assume the stats variable to
  !    have IDENTICAL SHORT NAME WITH THE CORRESPONDING FULL FIELD from where they are calculated (assuming there is one).
  !    This allows you to use automatic fetching of data using mo_stats_finder in the subroutines found in mo_ts_procedures. 
  ! 2. Stats variables associated with an onDemand procedure, but a UNIQUE short name: An explicit paring of the stats variable
  !    and the full field source variable as well as the calculation of the statistic must be implemented in mo_ts_procedures
  !    and associated with the onDemand procedure.
  ! 3. Stats variables associated with storage space and no onDemand procedure: calculate the values during the timestep
  !    flagged for statistics output (tsflg in mo_output.f90) within the appropriate physics subroutine.
  !
  ! Regarding 3), the likely case would be accumulation statistics, which require memory. These are yet to be implemented.
  ! Note that for all cases 1)-3), MPI reduction operations should be applied before writing output, since the
  ! statistics io will only be performed by rank=mpi_root (from mpi_interface.f90)

  ! Adding new variables:
  ! Declare your new TYPE(FloatArray0d) variable. Add initialization and registration of the variable to the TYPE(FieldArray) :: TS
  ! in the subroutine setTSVariables found below. Just check how other variables are defined and do the same. Only thing to consider
  ! are the 3 ways to calculate the statistics listed above.
  !   - If you want to use onDemand, first see if mo_ts_procedures already has a suitable subroutine to calculate your variable
  !     (e.g. simple domain mean is already implemented).
  !   - If so, and your variable fits point 1) for calculating the statistic, just associate the onDemand variable in the initialization
  !     and you're good to go. Nothing else needed.
  !   - If instead you see a suitable subroutine but your variable fits point 2), you must add a new paring for your variable and a
  !     source variable name, and the calculation of the statistic in mo_ts_procedures.
  !   - If you don't want to use onDemand, associate the variable with the storage instead. Now you need to calculate the value during
  !     a statistics timestep, or accumulate the value over consecutive timesteps using the TYPE(FloatArray1d) variable you declared.

  ! On-demand statistics; Actually I'm not sure if it's really necessary to define these instances for each variable,
  ! since (at least currently) they're only used to relay the info to the FieldArray instance TS. On the other hand,
  ! the list of instances makes it easy to keep track on the implemented variables, and they could also be used as access
  ! points for the corresponding parameters not only in output, but also in model processes if they happen to need them.
  TYPE(FloatArray0d), TARGET :: ts_wmax,ts_wmin,                        &
                                ts_tsrf, ts_ustar,                      &  
                                ts_shf,ts_lhf,ts_lwp,ts_rwp,ts_iwp,     &  
                                tsic_lwp, tsic_iwp,                     &  
                                tspr_rwp,                               &  
                                ts_ctop,ts_cbase,ts_cfrac,              &  
                                ts_sfcrrate,tspr_sfcrrate,                &  
                                ts_SSmax, ts_SSimax

  ! Piggybacking slave microphysics diagnostics
  TYPE(FloatArray0d), TARGET :: ts_pb_sfcrrate, tspr_pb_sfcrrate

  ! NOT IMPLEMENTED
  ! Accumulated statistics
  TYPE(FloatArray0d), TARGET :: ts_shfacc,ts_lhfacc,                    &  
                                ts_lwpacc,ts_rwpacc,ts_iwpacc                               

                                
  ! Storage for accumulated statistics                                 
  REAL, ALLOCATABLE, TARGET :: a_ts0d(:)
  INTEGER, PARAMETER        :: nts0d = 5

  
  CONTAINS

    SUBROUTINE setTSVariables(TS,outputlist,level,lpback)
      TYPE(FieldArray), INTENT(inout) :: TS
      CHARACTER(len=*), INTENT(in)    :: outputlist(:)
      INTEGER, INTENT(in)             :: level
      LOGICAL, INTENT(in)             :: lpback
      CLASS(*), POINTER :: pipeline => NULL()

      INTEGER :: nts

      ! ----------------------------
      ! Accumulated statistics
      ! ----------------------------

      ! Accumulation procedures have not been implemented yet.
      
      ALLOCATE(a_ts0d(nts0d))
      a_ts0d = 0.
      nts = 0

      nts = nts + 1
      ts_shfacc = FloatArray0d("shfacc",trgt=a_ts0d(nts))
      pipeline => ts_shfacc
      CALL TS%newField(ts_shfacc%shortName, "Acc mean sensible heat flux", "W m-2",   &
                       "time", ANY(outputlist == ts_shfacc%shortName), pipeline, in_group=["tsacc"])

      nts = nts + 1
      ts_lhfacc = FloatArray0d("lhfacc",trgt=a_ts0d(nts))
      pipeline => ts_lhfacc
      CALL TS%newField(ts_lhfacc%shortName, "Acc mean latent heat flux", "W m-2",     &
                       "time", ANY(outputlist == ts_lhfacc%shortName), pipeline, in_group=["tsacc"])

      nts = nts + 1
      ts_lwpacc = FloatArray0d("lwpacc",trgt=a_ts0d(nts))
      pipeline => ts_lwpacc
      CALL TS%newField(ts_lwpacc%shortName, "Acc mean liquid water path", "kg m-2",   &
                       "time", ANY(outputlist == ts_lwpacc%shortName), pipeline, in_group=["tsacc"])

      IF (level >= 3) THEN
         nts = nts + 1
         ts_rwpacc = FloatArray0d("rwpacc",trgt=a_ts0d(nts))
         pipeline => ts_rwpacc
         CALL TS%newField(ts_rwpacc%shortName,"Acc mean rain water path", "kg m-2",   &
                          "time", ANY(outputlist == ts_rwpacc%shortName), pipeline, in_group=["tsacc"])
      END IF

      IF (level == 5) THEN
         nts = nts + 1
         ts_iwpacc = FloatArray0d("iwpacc",trgt=a_ts0d(nts))
         pipeline => ts_iwpacc
         CALL TS%newField(ts_iwpacc%shortName, "Acc mean ice water path", "kg m-2",   &
                          "time", ANY(outputlist == ts_iwpacc%shortName), pipeline, in_group=["tsacc"])
      END IF
      
      
      ! ----------------------
      ! On-demand statistics
      ! ----------------------
      ts_wmax = FloatArray0d("wmax",srcname="wwind")
      ts_wmax%onDemand => tsMax
      pipeline => ts_wmax
      CALL TS%newField(ts_wmax%shortName, "Maximum vertical velocity", "m s-1", "time",    &
                       ANY(outputlist == ts_wmax%shortName), pipeline)

      ts_wmin = FloatArray0d("wmin",srcname="wwind")
      ts_wmin%onDemand => tsMin
      pipeline => ts_wmin
      CALL TS%newField(ts_wmin%shortName, "Minimum vertical velocity", "m s-1", "time",   &   
                       ANY(outputlist == ts_wmin%shortName), pipeline)

      ts_ustar = FloatArray0d("ustar",srcname="ustar")
      ts_ustar%onDemand => tsSfcMean
      pipeline => ts_ustar
      CALL TS%newField(ts_ustar%shortName, "Surface friction velocity", "m s-1", "time",   &
                       ANY(outputlist == ts_ustar%shortName), pipeline)
      
      ts_tsrf = FloatArray0d("tsrf",srcname="temp")
      ts_tsrf%onDemand => tsSfcMean
      pipeline => ts_tsrf
      CALL TS%newField(ts_tsrf%shortName, "Surface temperature", "K", "time",    &
                       ANY(outputlist == ts_tsrf%shortName), pipeline)

      ts_shf = FloatArray0d("shf",srcname="shf")
      ts_shf%onDemand => tsSfcMean
      pipeline => ts_shf
      CALL TS%newField(ts_shf%shortName, "Surface sensible heat flux", "W m-2", "time",   &
                       ANY(outputlist == ts_shf%shortName), pipeline)

      ts_lhf = FloatArray0d("lhf",srcname="lhf")
      ts_lhf%onDemand => tsSfcMean
      pipeline => ts_lhf
      CALL TS%newField(ts_lhf%shortName, "Surface latent heat flux", "W m-2", "time",   &
                       ANY(outputlist == ts_lhf%shortName), pipeline)

      ts_sfcrrate = FloatArray0d("sfcrrate")
      ts_sfcrrate%onDemand => tsSfcMean
      pipeline => ts_sfcrrate
      CALL TS%newField(ts_sfcrrate%shortName, "Surface rain rate", "W m-2", "time",   &
                       ANY(outputlist == ts_sfcrrate%shortName), pipeline)

      tspr_sfcrrate = FloatArray0d("pr_sfcrrate",srcname="sfcrrate")
      tspr_sfcrrate%onDemand => tsInPrecipMean
      pipeline => tspr_sfcrrate
      CALL TS%newField(tspr_sfcrrate%shortName, "Surface rain rate in precipitating grid points", "W m-2", "time",   &
                       ANY(outputlist == tspr_sfcrrate%shortName), pipeline)
      
      ts_lwp = FloatArray0d("lwp",srcname="lwp")
      ts_lwp%onDemand => tsSfcMean
      pipeline => ts_lwp
      CALL TS%newField(ts_lwp%shortName, "Liquid water path", "kg m-2", "time",   &
                       ANY(outputlist == ts_lwp%shortName), pipeline)

      tsic_lwp = FloatArray0d("ic_lwp",srcname="lwp")
      tsic_lwp%onDemand => tsInLiqMean
      pipeline => tsic_lwp
      CALL TS%newField(tsic_lwp%shortName, "Liquid water path in cloudy columns", "kg m-2", "time",   &
                       ANY(outputlist == tsic_lwp%shortName), pipeline)
      
      IF (level >= 3) THEN
         ts_rwp = FloatArray0d("rwp",srcname="rwp")
         ts_rwp%onDemand => tsSfcMean
         pipeline => ts_rwp
         CALL TS%newField(ts_rwp%shortName,"Rain water path", "kg m-2", "time",    &
                          ANY(outputlist == ts_rwp%shortName), pipeline)

         tspr_rwp = FloatArray0d("pr_rwp",srcname="rwp")
         tspr_rwp%onDemand => tsInPrecipMean
         pipeline => tspr_rwp
         CALL TS%newField(tspr_rwp%shortName, "Rain water path in precipitating columns", "kg m-2", "time",   &
                          ANY(outputlist == tspr_rwp%shortName), pipeline)
         
      END IF
      
      IF (level == 5) THEN
         ts_iwp = FloatArray0d("iwp",srcname="iwp")
         ts_iwp%onDemand => tsSfcMean
         pipeline => ts_iwp
         CALL TS%newField(ts_iwp%shortName, "Ice water path", "kg m-2", "time",  &
                          ANY(outputlist == ts_iwp%shortName), pipeline)

         tsic_iwp = FloatArray0d("ic_iwp",srcname="iwp")
         tsic_iwp%onDemand => tsInIceMean
         pipeline => tsic_iwp
         CALL TS%newField(tsic_iwp%shortName, "Ice water path in cloudy columns", "kg m-2", "time",  &
                          ANY(outputlist == tsic_iwp%shortName), pipeline)
      END IF

      ts_ctop = FloatArray0d("ctop")
      ts_ctop%onDemand => tscloudBoundaries
      pipeline => ts_ctop
      CALL TS%newField(ts_ctop%shortName, "Cloud top height", "m", "time",   &
                       ANY(outputlist == ts_ctop%shortName), pipeline)

      ts_cbase = FloatArray0d("cbase")
      ts_cbase%onDemand => tscloudBoundaries
      pipeline => ts_cbase
      CALL TS%newField(ts_cbase%shortName, "Cloud base height", "m", "time",   &
                       ANY(outputlist == ts_cbase%shortName), pipeline)

      ts_cfrac = FloatArray0d("cfrac")
      ts_cfrac%onDemand => tscloudFraction
      pipeline => ts_cfrac
      CALL TS%newField(ts_cfrac%shortName, "Cloud fraction", "1", "time",   &
                       ANY(outputlist == ts_cfrac%shortName), pipeline)
      
      ts_SSmax = FloatArray0d("SSmax",srcname="rh")
      ts_SSmax%onDemand => tsMax
      pipeline => ts_SSmax
      CALL TS%newField(ts_SSmax%shortName, "Maximum supersaturation", "1", "time",  &
                       ANY(outputlist == ts_SSmax%shortName), pipeline)

      IF (level == 5) THEN
         ts_SSimax = FloatArray0d("SSimax",srcname="rhi")
         ts_SSimax%onDemand => tsMax
         pipeline => ts_SSimax
         CALL TS%newField(ts_SSimax%shortName, "Maximum supersaturation w.r.t. ice (for T < 273K)", "1", "time",   &
                          ANY(outputlist == ts_SSimax%shortName), pipeline)
      END IF

      IF (lpback) THEN

         ts_pb_sfcrrate = FloatArray0d("pb_sfcrrate")
         ts_pb_sfcrrate%onDemand => tsSfcMean
         pipeline => ts_pb_sfcrrate
         CALL TS%newField(ts_pb_sfcrrate%shortName, "Surface precip flux (bulk slave)", "W m-2", "time",   &
                          ANY(outputlist == ts_pb_sfcrrate%shortName), pipeline)
         
         tspr_pb_sfcrrate = FloatArray0d("pr_pb_sfcrrate",srcname="pb_sfcrrate")
         tspr_pb_sfcrrate%onDemand => tsInPrecipMean
         pipeline => tspr_pb_sfcrrate
         CALL TS%newField(tspr_pb_sfcrrate%shortName, "Conditional surface precip flux (bulk slave)", "W m-2", "time",  &
                          ANY(outputlist == tspr_pb_sfcrrate%shortName), pipeline)
         
      END IF

      
    END SUBROUTINE setTSVariables
        
END MODULE mo_ts_state
