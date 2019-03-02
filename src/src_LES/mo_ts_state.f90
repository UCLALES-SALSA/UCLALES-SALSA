MODULE mo_ts_state
  USE classFieldArray
  USE mo_structured_datatypes
  USE mo_check_state, ONLY : checkOutputs
  IMPLICIT NONE

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
      
      IF (ANY(outputlist == "vtke")) THEN
         nts = nts + 1
         vtke = FloatArray0d(a_ts0d(nts))
         pipeline => vtke
         CALL TS%newField("vtke", "Vertical integral of total TKE", "kg/s", "time",   &
                          ANY(outputlist == "vtke"), pipeline)
      END IF
      
      IF (ANY(outputlist == "sfcbflux")) THEN
         nts = nts + 1
         sfcbflux = FloatArray0d(a_ts0d(nts))
         pipeline => sfcbflux
         CALL TS%newField("sfcbflux", "Surface buoyancy flux", "m/s2", "time",   &
                          ANY(outputlist == "sfcbflux"), pipeline)
      END IF

      IF (ANY(outputlist == "wmax")) THEN
         nts = nts + 1
         wmax = FloatArray0d(a_ts0d(nts))
         pipeline => wmax
         CALL TS%newField("wmax", "Maximum vertical velocity", "m/s", "time",    &
                          ANY(outputlist == "wmax"), pipeline)
      END IF

      IF (ANY(outputlist == "tsrf")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => tsrf
         CALL TS%newField("tsrf", "Mean surface temperature", "K", "time",   &
                          ANY(outputlist == "tsrf"), pipeline)
      END IF

      IF (ANY(outputlist == "shf_bar")) THEN
         nts = nts + 1
         shf_bar = FloatArray0d(a_ts0d(nts))
         pipeline => shf_bar
         CALL TS%newField("shf_bar", "Mean surface sensible heat flux", "W/m2", "time",   &
                          ANY(outputlist == "shf_bar"), pipeline)
      END IF

      IF (ANY(outputlist == "lhf_bar")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => lhf_bar
         CALL TS%newField("lhf_bar", "Mean surface latent heat flux", "W/m2", "time",   &
                          ANY(outputlist == "lhf_bar"), pipeline)
      END IF

      IF (ANY(outputlist == "lwp_bar")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => lwp_bar
         CALL TS%newField("lwp_bar", "Mean liquid water path (cloud water)", "kg/m2", "time",   &
                          ANY(outputlist == "lwp_bar"), pipeline)
      END IF

      IF (ANY(outputlist == "rwp_bar")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => rwp_bar
         CALL TS%newField("rwp_bar", "Mean precipitating water path", "kg/m2", "time",   &
                          ANY(outputlist == "rwp_bar"), pipeline)
      END IF

      IF (ANY(outputlist == "iwp_bar") .AND. level == 5) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => iwp_bar
         CALL TS%newField("iwp_bar", "Mean frozen water path", "kg/m2", "time",   &
                          ANY(outputlist == "iwp_bar"), pipeline)
      END IF
      
      IF (ANY(outputlist == "ctop")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => ctop
         CALL TS%newField("ctop", "Mean cloud top height", "m", "time",   &
                          ANY(outputlist == "ctop"), pipeline)
      END IF

      IF (ANY(outputlist == "cbase")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => cbase
         CALL TS%newField("cbase", "Mean cloud base height", "m", "time",   &
                          ANY(outputlist == "cbase"), pipeline)
      END IF

      IF (ANY(outputlist == "cfrac")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => cfrac
         CALL TS%newField("cfrac", "Mean total cloud fraction", "m", "time",   &
                          ANY(outputlist == "cfrac"), pipeline)
      END IF
      
      IF (ANY(outputlist == "sfc_rain")) THEN
         nts = nts + 1
         tsrf = FloatArray0d(a_ts0d(nts))
         pipeline => sfc_rain
         CALL TS%newField("sfc_rain", "Mean surface rain rate", "mm/hour", "time",   &
                          ANY(outputlist == "sfc_rain"), pipeline)
      END IF

      ! Check the user specified output variable list for bad entries
      CALL checkOutputs(outputlist,TS)
      
    END SUBROUTINE setTSVariables
        
END MODULE mo_ts_state
