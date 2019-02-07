MODULE mo_diag_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  IMPLICIT NONE

  SAVE
  
  !
  ! Mandatory diagnostic quantities
  !
  TYPE(FloatArray3D), TARGET :: a_theta  ! dry potential temp (k)
  TYPE(FloatArray3D), TARGET :: a_temp   ! Absolute temperature (K)
  TYPE(FloatArray3D), TARGET :: a_pexnr  ! perturbation exner func
  TYPE(FloatArray3D), TARGET :: a_press  ! pressure (hpa)
  TYPE(FloatArray3D), TARGET :: a_rc     ! Total cloud water +rain (level<=3) or aerosol+cloud (level>=4) water mixing ratio
  TYPE(FloatArray3D), TARGET :: a_ri     ! Unrimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_riri   ! Rimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rv     ! Water vapor mixing ratio  (only for levels < 4!)
  TYPE(FloatArray3D), TARGET :: a_srp    ! Bulk precipitation mixing ratio (levels >= 4)
  TYPE(FloatArray3D), TARGET :: a_snrp   ! Bulk precipitation number mixing ratio (levels >=4) 
  TYPE(FloatArray3D), TARGET :: a_rh     ! Relative humidity
  TYPE(FloatArray3D), TARGET :: a_rsl    ! Water vapor saturation mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rhi    ! Relative humidity over ice
  TYPE(FloatArray3D), TARGET :: a_rsi    ! Water vapor saturation mixing ratio over ice
  TYPE(FloatArray3D), TARGET :: a_dn     ! Air density
  TYPE(FloatArray3D), TARGET :: a_rflx, a_sflx,   &  ! Radiation fluxes
                                a_fus, a_fds,     &
                                a_fuir, a_fdir

  TYPE(FloatArray2D), TARGET :: albedo               ! Albedo, CHECK DEFINITION
  TYPE(FloatArray2D), TARGET :: a_ustar              ! Friction velocity
  TYPE(FloatArray2D), TARGET :: a_tstar              ! turbulent temperature scale 
  TYPE(FloatArray2D), TARGET :: a_rstar              ! turbulent moisture scale
  TYPE(FloatArray2D), TARGET :: uw_sfc               ! Surface fluxes
  TYPE(FloatArray2D), TARGET :: vw_sfc
  TYPE(FloatArray2D), TARGET :: ww_sfc
  TYPE(FloatArray2D), TARGET :: wt_sfc
  TYPE(FloatArray2D), TARGET :: wq_sfc
  TYPE(FloatArray3D), TARGET :: a_rrate               ! Precipitation flux
  TYPE(FloatArray3D), TARGET :: a_irate            ! Precipitation flux, frozen

  
  CONTAINS

    SUBROUTINE setDiagnosticVariables(Diag,outputlist,memsize,level,iradtyp,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Diag
      CHARACTER(len=10), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(inout) :: memsize
      INTEGER, INTENT(in) :: level, iradtyp, nzp,nxp,nyp      
      CLASS(*), POINTER :: pipeline
      REAL :: zeros3d(nzp,nxp,nyp), zeros2d(nxp,nyp)
      INTEGER :: nxyz, nxy

      nxyz = nxp*nyp*nzp
      nxy = nxp*nyp
      
      zeros3d = 0.   ! To initialize the arrays
      zeros2d = 0.

      memsize = memsize + nxyz
      a_theta = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_theta
      CALL Diag%NewField("theta", "Potential temperature", "K", "tttt",    &
                         ANY(outputlist == "theta"), pipeline)

      memsize = memsize + nxyz
      a_temp = FloatArray3D(zeros3d,store=.TRUE.)
      pipeline => a_temp
      CALL Diag%NewField("temp", "Abolute temperature", "K", "tttt",      &
                         ANY(outputlist == "temp"), pipeline)

      memsize = memsize + nxyz
      a_pexnr = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_pexnr
      CALL Diag%NewField("penxr", "Exner function", " ", "tttt",       &
                         ANY(outputlist == "pexnr"), pipeline)

      memsize = memsize + nxyz
      a_press = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_press
      CALL Diag%NewField("press", "Pressure", "Pa", "tttt",       &
                         ANY(outputlist == "press"), pipeline)

      IF (level > 1) THEN
         memsize = memsize + nxyz
         a_rc = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rc
         CALL Diag%NewField("rc", "Cloud mixing ratio", "kg/kg", "tttt",    &
                            ANY(outputlist == "rc"), pipeline)
      END IF

      IF (level == 5) THEN
         memsize = memsize + nxyz
         a_ri = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_ri
         CALL Diag%NewField("ri", "Unrimed ice mixing ratio", "g/kg", "tttt",  &
                            ANY(outputlist == "ri"), pipeline)

         memsize = memsize + nxyz
         a_riri = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_riri
         CALL Diag%NewField("riri", "Rimed ice mixing ratio", "kg/kg", "tttt",  &
                            ANY(outputlist == "riri"), pipeline)         
      END IF

      IF (level < 4) THEN
         memsize = memsize + nxyz
         a_rv = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rv
         CALL Diag%NewField("rv", "Water vapor mixing ratio", "kg/kg", "tttt",    &
                            ANY(outputlist == "rv"), pipeline)
      END IF

      IF (level >=4 ) THEN
         memsize = memsize + nxyz
         a_srp = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_srp
         CALL Diag%NewField("srp", "Precipitation mixing ratio", "kg/kg", "tttt",  &
                            ANY(outputlist == "srp"), pipeline)

         memsize = memsize + nxyz
         a_snrp = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_snrp
         CALL Diag%NewField("snrp", "Precipitation number mixing ratio", "#/kg", "tttt", &
                            ANY(outputlist == "snrp"), pipeline)
      END IF

      memsize = memsize + nxyz
      a_rh = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_rh
      CALL Diag%NewField("rh", "Relative humidity", "1", "tttt",    &
                         ANY(outputlist == "rh"), pipeline)

      memsize = memsize + nxyz
      a_rsl = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_rsl
      CALL Diag%NewField("rsl", "Liquid saturation mixing ratio", "kg/kg", "tttt",    &
                         ANY(outputlist == "rsl"), pipeline)

      IF (level == 5) THEN
         memsize = memsize + nxyz
         a_rhi = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rhi
         CALL Diag%NewField("rhi", "Relative humidity over ice", "1", "tttt",   &
                            ANY(outputlist == "rhi"), pipeline)

         memsize = memsize + nxyz
         a_rsi = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rsi
         CALL Diag%NewField("rsi", "Ice saturation mixing ratio", "kg/kg", "tttt",   &
                            ANY(outputlist == "rsi"), pipeline)
      END IF

      memsize = memsize + nxyz
      a_dn = FloatArray3d(zeros3d,store=.TRUE.)
      pipeline => a_dn
      CALL Diag%NewField("dn", "Air density", "kg/m3", "tttt",   &
                         ANY(outputlist == "dn"), pipeline)

      IF (iradtyp > 0) THEN
         memsize = memsize + nxyz
         a_rflx = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rflx
         CALL Diag%NewField("rflx", "Net longwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "rflx"), pipeline)           
      END IF

      IF (iradtyp >= 3) THEN
         memsize = memsize + nxyz
         a_sflx = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_sflx
         CALL Diag%NewField("sflx", "Net shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "sflx"), pipeline)

         memsize = memsize + nxy
         albedo = FloatArray2d(zeros2d,store=.TRUE.)   ! FIKSAA MAHDOLLISUUS 2D OUTPUTILLE
         pipeline => albedo
         CALL Diag%NewField("albedo", "Albedo", "1", "xtytt",     &
                            ANY(outputlist == "albedo"), pipeline)

         memsize = memsize + nxyz
         a_fus = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_fus
         CALL Diag%NewField("fus", "Upwelling shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fus"), pipeline)

         memsize = memsize + nxyz
         a_fds = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_fds
         CALL Diag%NewField("fds", "Downwelling shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fds"), pipeline)

         memsize = memsize + nxyz
         a_fuir = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_fuir
         CALL Diag%NewField("fuir", "Upwelling longwave flux", "W/m2", "tttt",    &
                            ANY(outputlist == "fuir"), pipeline)         

         memsize = memsize + nxyz
         a_fdir = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_fdir
         CALL Diag%NewField("fdir", "Downwelling longwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fdir"), pipeline)               
      END IF

      memsize = memsize + nxy
      a_ustar = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => a_ustar
      CALL Diag%NewField("ustar", "Turbulent friction velocity", "m/s", "xtytt",   &
                         ANY(outputlist == "ustar"), pipeline) 

      memsize = memsize + nxy
      a_tstar = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => a_tstar
      CALL Diag%NewField("tstar", "Turbulent scale temperature", "K", "xtytt",    &
                         ANY(outputlist == "tstar"), pipeline) 

      memsize = memsize + nxy
      a_rstar = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => a_rstar
      CALL Diag%NewField("rstar", "Turbulent scale CHECK", "CHECK", "xtytt",    &
                         ANY(outputlist == "rstar"), pipeline) 

      memsize = memsize + nxy
      uw_sfc = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => uw_sfc
      CALL Diag%NewField("uw_sfc", "Vertical momentum flux with u wind", "CHECK", "xtytt",   &
                         ANY(outputlist == "uw_sfc"), pipeline) 

      memsize = memsize + nxy
      vw_sfc = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => vw_sfc
      CALL Diag%NewField("vw_sfc", "Vertical momentum flux with v wind", "CHECK", "xtytt",    &
                         ANY(outputlist == "vw_sfc"), pipeline) 

      memsize = memsize + nxy
      ww_sfc = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => ww_sfc
      CALL Diag%NewField("ww_sfc", "Vertical wind covariance", "CHECK", "xtytt",    &
                         ANY(outputlist == "ww_sfc"), pipeline) 

      memsize = memsize + nxy
      wt_sfc = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => wt_sfc
      CALL Diag%NewField("wt_sfc", "Vertical temperature flux", "CHECK", "xtytt",   &
                         ANY(outputlist == "wt_sfc"), pipeline) 

      memsize = memsize + nxy
      wq_sfc = FloatArray2d(zeros2d,store=.TRUE.)
      pipeline => wq_sfc
      CALL Diag%NewField("wq_sfc", "Vertical moisture flux", "CHECK", "xtytt",    &
                         ANY(outputlist == "wq_sfc"), pipeline) 

      IF (level >= 3) THEN
         memsize = memsize + nxy
         a_rrate = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_rrate
         CALL Diag%NewField("rrate", "Liquid surface precipitation", "CHECK", "xtytt",   &
                            ANY(outputlist == "rrate"), pipeline) 
      END IF

      IF (level == 5) THEN
         memsize = memsize + nxy
         a_irate = FloatArray3d(zeros3d,store=.TRUE.)
         pipeline => a_irate
         CALL Diag%NewField("irate", "Frozen surface precipitation", "CHECK", "xtytt",    &
                            ANY(outputlist == "irate"), pipeline) 
      END IF
         
                 
    END SUBROUTINE

    !!
        
END MODULE mo_diag_state
