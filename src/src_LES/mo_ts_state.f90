MODULE mo_ts_state
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
         CALL TS%newField("tsrf", "Surface temperature", "K", "time",   &
                          ANY(outputlist == "tsrf"), pipeline)

      END IF
      

    END SUBROUTINE setTSVariables
    
  
  

END MODULE mo_ts_state
