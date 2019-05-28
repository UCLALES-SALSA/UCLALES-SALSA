MODULE mo_diag_state
  USE util, ONLY : Extend_last
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  IMPLICIT NONE

  SAVE
  
  !
  ! Mandatory diagnostic quantities
  !
  TYPE(FloatArray3D), TARGET :: a_theta  ! 1: dry potential temp (k)
  TYPE(FloatArray3D), TARGET :: a_temp   ! 2: Absolute temperature (K)
  TYPE(FloatArray3D), TARGET :: a_pexnr  ! 3: perturbation exner func
  TYPE(FloatArray3D), TARGET :: a_press  ! 4: pressure (hpa)
  TYPE(FloatArray3D), TARGET :: a_rc     ! 5: Total cloud water +rain (level<=3) or aerosol+cloud (level>=4) water mixing ratio
  TYPE(FloatArray3D), TARGET :: a_ri     ! 6: Unrimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_riri   ! 7: Rimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rv     ! 8: Water vapor mixing ratio  (only for levels < 4!)
  TYPE(FloatArray3D), TARGET :: a_srp    ! 9: Bulk precipitation mixing ratio (levels >= 4)
  TYPE(FloatArray3D), TARGET :: a_snrp   ! 10: Bulk precipitation number mixing ratio (levels >=4) 
  TYPE(FloatArray3D), TARGET :: a_rh     ! 11: Relative humidity
  TYPE(FloatArray3D), TARGET :: a_rsl    ! 12: Water vapor saturation mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rhi    ! 13: Relative humidity over ice
  TYPE(FloatArray3D), TARGET :: a_rsi    ! 14: Water vapor saturation mixing ratio over ice
  TYPE(FloatArray3D), TARGET :: a_dn     ! 15: Air density
  TYPE(FloatArray3D), TARGET :: a_rflx, a_sflx,   &  ! 16, 17: Radiation fluxes
                                a_fus, a_fds,     &  ! 18, 19: 
                                a_fuir, a_fdir       ! 20, 21:
  TYPE(FloatArray3D), TARGET :: a_rrate              ! 22: Precipitation flux
  TYPE(FloatArray3D), TARGET :: a_irate              ! 23: Precipitation flux, frozen

  REAL, ALLOCATABLE, TARGET :: a_diag3d(:,:,:,:) 
  INTEGER, PARAMETER :: ndiag3d = 23   ! Remember to update if adding new variables!!
  
  !
  ! Two dimensional variables that need to be stored during the timestep
  !
  TYPE(FloatArray2D), TARGET :: albedo               ! 1: Albedo, CHECK DEFINITION
  TYPE(FloatArray2D), TARGET :: a_ustar              ! 2: Friction velocity
  TYPE(FloatArray2D), TARGET :: a_tstar              ! 3: turbulent temperature scale 
  TYPE(FloatArray2D), TARGET :: a_rstar              ! 4: turbulent moisture scale
  TYPE(FloatArray2D), TARGET :: uw_sfc               ! 5: Surface fluxes
  TYPE(FloatArray2D), TARGET :: vw_sfc               ! 6: 
  TYPE(FloatArray2D), TARGET :: ww_sfc               ! 7:
  TYPE(FloatArray2D), TARGET :: wt_sfc               ! 8: 
  TYPE(FloatArray2D), TARGET :: wq_sfc               ! 9:

  TYPE(FloatArray2D), TARGET :: a_sfcrrate           ! 10: Surface rain rate
  TYPE(FloatArray2D), TARGET :: a_sfcirate           ! 11: Surface frozen precipitation
  
  REAL, ALLOCATABLE, TARGET :: a_diag2d(:,:,:)
  INTEGER, PARAMETER :: ndiag2d = 11  ! Remember to update if adding new variables!!
  
  !
  ! Microphysical process rates -- they need to be stored during the timestep because they cannot be simply diagnosed afterwards
  ! For now, these are BULK process rates only for water/ice, except where indicated otherwise !! Number concentration rate given for particle formation processes,
  ! mass concetration rate for others
  !
  TYPE(FloatArray3d), TARGET :: m_autoconversion    ! 1: Bulk autoconversion rate (mass)
  TYPE(FloatArray3d), TARGET :: m_accretion         ! 2: Bulk accretion rate (mass)
  TYPE(FloatArray3d), TARGET :: m_ACcoll_dry        ! 3: Bulk cloud collection of aerosol (dry aerosol mass)
  TYPE(FloatArray3d), TARGET :: m_APcoll_dry        ! 4: Bulk precipitation collection of aerosol (dry aerosol mass)            
  TYPE(FloatArray3d), TARGET :: m_AIcoll_dry        ! 5: Bulk ice collection of aerosol (dry aerosol mass)
  TYPE(FloatArray3d), TARGET :: n_activation        ! 6: Bulk cloud activation rate (number)
  TYPE(FloatArray3d), TARGET :: n_icehom            ! 7: Bulk homogeneous freezing rate (number)
  TYPE(FloatArray3d), TARGET :: n_icedep            ! 8: Bulk deposition freezing rate (number)
  TYPE(Floatarray3d), TARGET :: n_iceimm            ! 9: Bulk immersion freezing rate (number)
  TYPE(FloatArray3d), TARGET :: m_conda             ! 10: Bulk condensation rate of water on aerosol (mass)
  TYPE(FloatArray3d), TARGET :: m_condc             ! 11: Bulk condensation rate of water on cloud droplets (mass)
  TYPE(FloatArray3d), TARGET :: m_condp             ! 12: Bulk condensation rate of water on precipitation (mass)
  TYPE(FloatArray3d), TARGET :: m_condi             ! 13: Bulk condensation (deposition) rate of water on ice (mass)

  REAL, ALLOCATABLE, TARGET :: a_rateDiag3d(:,:,:,:)
  INTEGER, PARAMETER :: nratediag3d = 13 ! Remember to update if adding new variables!!

  CONTAINS

    SUBROUTINE setDiagnosticVariables(Diag,outputlist,memsize,level,iradtyp,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Diag
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(inout) :: memsize
      INTEGER, INTENT(in) :: level, iradtyp, nzp,nxp,nyp      
      CLASS(*), POINTER :: pipeline => NULL()
      INTEGER :: nxyz, nxy, n2d,n3d,n3dr

      nxyz = nxp*nyp*nzp
      nxy = nxp*nyp

      n3d = 0
      n3dr = 0
      n2d = 0
      
      ALLOCATE(a_diag3d(nzp,nxp,nyp,ndiag3d))
      ALLOCATE(a_diag2d(nxp,nyp,ndiag2d))
      ALLOCATE(a_rateDiag3d(nzp,nxp,nyp,nratediag3d))

      a_diag3d = 0.
      a_diag2d = 0.
      a_rateDiag3d = 0.
      
      ! First entry for a_diag3d
      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_theta = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_theta
      CALL Diag%newField("theta", "Potential temperature", "K", "tttt",    &
                         ANY(outputlist == "theta"), pipeline)
      
      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_temp = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_temp
      CALL Diag%newField("temp", "Abolute temperature", "K", "tttt",      &
                         ANY(outputlist == "temp"), pipeline)
      
      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_pexnr = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_pexnr
      CALL Diag%newField("pexnr", "Exner function", "1", "tttt",       &
                         ANY(outputlist == "pexnr"), pipeline)

      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_press = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_press
      CALL Diag%newField("press", "Pressure", "Pa", "tttt",       &
                         ANY(outputlist == "press"), pipeline)

      IF (level > 1) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_rc = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rc
         CALL Diag%newField("rc", "Cloud mixing ratio", "kg/kg", "tttt",    &
                            ANY(outputlist == "rc"), pipeline)
      END IF

      IF (level == 5) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_ri = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_ri
         CALL Diag%newField("ri", "Unrimed ice mixing ratio", "g/kg", "tttt",  &
                            ANY(outputlist == "ri"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_riri = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_riri
         CALL Diag%newField("riri", "Rimed ice mixing ratio", "kg/kg", "tttt",  &
                            ANY(outputlist == "riri"), pipeline)         
      END IF

      IF (level < 4) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_rv = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rv
         CALL Diag%newField("rv", "Water vapor mixing ratio", "kg/kg", "tttt",    &
                            ANY(outputlist == "rv"), pipeline)
      END IF

      IF (level >=4 ) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_srp = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_srp
         CALL Diag%newField("srp", "Precipitation mixing ratio", "kg/kg", "tttt",  &
                            ANY(outputlist == "srp"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_snrp = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_snrp
         CALL Diag%newField("snrp", "Precipitation number mixing ratio", "#/kg", "tttt", &
                            ANY(outputlist == "snrp"), pipeline)
      END IF

      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_rh = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_rh
      CALL Diag%newField("rh", "Relative humidity", "1", "tttt",    &
                         ANY(outputlist == "rh"), pipeline)

      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_rsl = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_rsl
      CALL Diag%newField("rsl", "Liquid saturation mixing ratio", "kg/kg", "tttt",    &
                         ANY(outputlist == "rsl"), pipeline)

      IF (level == 5) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_rhi = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rhi
         CALL Diag%newField("rhi", "Relative humidity over ice", "1", "tttt",   &
                            ANY(outputlist == "rhi"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_rsi = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rsi
         CALL Diag%newField("rsi", "Ice saturation mixing ratio", "kg/kg", "tttt",   &
                            ANY(outputlist == "rsi"), pipeline)
      END IF

      memsize = memsize + nxyz
      n3d = n3d+1
      pipeline => NULL()
      a_dn = FloatArray3d(a_diag3d(:,:,:,n3d))
      pipeline => a_dn
      CALL Diag%newField("dn", "Air density", "kg/m3", "tttt",   &
                         ANY(outputlist == "dn"), pipeline)

      IF (iradtyp > 0) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_rflx = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rflx
         CALL Diag%newField("rflx", "Net longwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "rflx"), pipeline)           
      END IF

      IF (iradtyp >= 3) THEN
         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_sflx = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_sflx
         CALL Diag%newField("sflx", "Net shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "sflx"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_fus = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_fus
         CALL Diag%newField("fus", "Upwelling shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fus"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_fds = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_fds
         CALL Diag%newField("fds", "Downwelling shortwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fds"), pipeline)

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_fuir = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_fuir
         CALL Diag%newField("fuir", "Upwelling longwave flux", "W/m2", "tttt",    &
                            ANY(outputlist == "fuir"), pipeline)         

         memsize = memsize + nxyz
         n3d = n3d+1
         pipeline => NULL()
         a_fdir = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_fdir
         CALL Diag%newField("fdir", "Downwelling longwave flux", "W/m2", "tttt",   &
                            ANY(outputlist == "fdir"), pipeline)               
      END IF

      IF (level >= 3) THEN
         memsize = memsize + nxy
         n3d = n3d+1
         pipeline => NULL()
         a_rrate = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rrate
         CALL Diag%newField("rrate", "Liquid precipitation flux", "W/m2", "tttt",   &
                            ANY(outputlist == "rrate"), pipeline) 
      END IF

      IF (level == 5) THEN
         memsize = memsize + nxy
         n3d = n3d+1
         pipeline => NULL()
         a_irate = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_irate
         CALL Diag%newField("irate", "Frozen precipitation flux", "W/m2", "tttt",    &
                            ANY(outputlist == "irate"), pipeline) 
      END IF

      ! First diag2d entry
      IF (iradtyp >= 3) THEN
         memsize = memsize + nxy
         n2d = n2d+1
         pipeline => NULL()
         albedo = FloatArray2d(a_diag2d(:,:,n2d))  
         pipeline => albedo
         CALL Diag%newField("albedo", "Albedo", "1", "xtytt",     &
                            ANY(outputlist == "albedo"), pipeline)
      END IF
         
      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      a_ustar = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => a_ustar
      CALL Diag%newField("ustar", "Turbulent friction velocity", "m/s", "xtytt",   &
                         ANY(outputlist == "ustar"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      a_tstar = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => a_tstar
      CALL Diag%newField("tstar", "Turbulent scale temperature", "K", "xtytt",    &
                         ANY(outputlist == "tstar"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      a_rstar = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => a_rstar
      CALL Diag%newField("rstar", "Turbulent scale CHECK", "CHECK", "xtytt",    &
                         ANY(outputlist == "rstar"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      uw_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => uw_sfc
      CALL Diag%newField("uw_sfc", "Vertical momentum flux with u wind", "CHECK", "xtytt",   &
                         ANY(outputlist == "uw_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      vw_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => vw_sfc
      CALL Diag%newField("vw_sfc", "Vertical momentum flux with v wind", "CHECK", "xtytt",    &
                         ANY(outputlist == "vw_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      ww_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => ww_sfc
      CALL Diag%newField("ww_sfc", "Vertical wind covariance", "CHECK", "xtytt",    &
                         ANY(outputlist == "ww_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      wt_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => wt_sfc
      CALL Diag%newField("wt_sfc", "Vertical temperature flux", "CHECK", "xtytt",   &
                         ANY(outputlist == "wt_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      wq_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => wq_sfc
      CALL Diag%newField("wq_sfc", "Vertical moisture flux", "CHECK", "xtytt",    &
                         ANY(outputlist == "wq_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      a_sfcrrate = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => a_sfcrrate
      CALL Diag%newField("sfcrrate", "Surface rain rate", "mm/h", "xtytt",     &
                         ANY(outputlist == "sfcrrate"), pipeline)

      IF ( level == 5 ) THEN
         memsize = memsize + nxy
         n2d = n2d+1
         pipeline => NULL()
         a_sfcirate = FloatArray2d(a_diag2d(:,:,n2d))
         pipeline => a_sfcirate
         CALL Diag%newField("sfcirate", "Surface frozen precip (liquid eqv)", "mm/h", "xtytt",    &
                            ANY(outputlist == "sfcirate"), pipeline)
      END IF
      

      ! First rateDiag3d entry
      n3dr = n3dr + 1
      pipeline => NULL()
      m_autoconversion = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_autoconversion
      CALL Diag%newField("autoconversion", "Autoconversion rate, h2o mass", "kg/kgs", "tttt",   &
                         ANY(outputlist == "autoconversion"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_accretion = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_accretion
      CALL Diag%newField("accretion", "Accretion rate, h2o mass", "kg/kgs", "tttt",   &
                         ANY(outputlist == "autoconversion"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_ACcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_ACcoll_dry
      CALL Diag%newField("ACcoll", "Cloud collection of aerosol, dry mass", "kg/kgs", "tttt",  &
                         ANY(outputlist == "ACcoll"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_APcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_APcoll_dry
      CALL Diag%newField("APcoll","Rain collection of aerosol, dry mass", "kg/kgs", "tttt",   &
                         ANY(outputlist == "APcoll"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_AIcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_AIcoll_dry
      CALL Diag%newField("AIcoll","Ice collection of aerosol, dry mass", "kg/kgs", "tttt",   &
                         ANY(outputlist == "AIcoll"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      n_activation = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => n_activation
      CALL Diag%newField("activation","Cloud activation rate, number", "#/kgs","tttt",   &
                         ANY(outputlist == "activation"), pipeline)
      
      n3dr = n3dr + 1
      pipeline => NULL()
      n_icehom = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => n_icehom
      CALL Diag%newfield("icehom","Homogeneous freezing rate, number", "#/kgs","tttt",   &
                         ANY(outputlist == "icehom"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      n_icedep = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => n_icedep
      CALL Diag%newField("icedep","Deposition freezing rate, number", "#/kgs", "tttt",  &
                         ANY(outputlist == "icedep"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      n_iceimm = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => n_iceimm
      CALL Diag%newField("iceimm","Immersion freezing rate, number", "#/kgs", "tttt",  &
                         ANY(outputlist == "iceimm"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_conda = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_conda
      CALL Diag%newField("conda","H2O condensation rate on aerosol", "#/kgs", "tttt",  &
                         ANY(outputlist == "conda"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_condc = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_condc
      CALL Diag%newField("condc","H2O condensation rate on cloud droplets", "#/kgs", "tttt",  &
                         ANY(outputlist == "condc"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_condp = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_condp
      CALL Diag%newField("condp","H2O condensation rate on precipiation", "#/kgs", "tttt",    &
                         ANY(outputlist == "condp"), pipeline)

      n3dr = n3dr + 1
      pipeline => NULL()
      m_condi = FloatArray3d(a_rateDiag3d(:,:,:,n3dr))
      pipeline => m_condi
      CALL Diag%newField("condi","H2O condensation rate on ice", "#/kgs", "tttt",   &
                         ANY(outputlist == "condi"), pipeline)
      
      pipeline => NULL()
      
    END SUBROUTINE
        
END MODULE mo_diag_state
