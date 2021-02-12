MODULE mo_diag_state
  USE util, ONLY : Extend_last
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : nprc, nice
  IMPLICIT NONE

  SAVE
  
  ! -----------------------------------------------------------------------
  ! Mandatory diagnostic quantities (stored)
  !
  TYPE(FloatArray3D), TARGET :: a_theta  ! 1: dry potential temp (k)
  TYPE(FloatArray3D), TARGET :: a_temp   ! 2: Absolute temperature (K)
  TYPE(FloatArray3D), TARGET :: a_pexnr  ! 3: perturbation exner func
  TYPE(FloatArray3D), TARGET :: a_press  ! 4: pressure (hpa)
  TYPE(FloatArray3D), TARGET :: a_rtot   ! 5: Total water mix rat for level >= 4 (vapor + condensate)
  TYPE(FloatArray3D), TARGET :: a_rc     ! 6: Total cloud water +rain (level<=3) or aerosol+cloud (level>=4) water mixing ratio
  TYPE(FloatArray3D), TARGET :: a_ri     ! 7: Unrimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_riri   ! 8: Rimed ice mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rv     ! 9: Water vapor mixing ratio  (only for levels < 4!)
  TYPE(FloatArray3D), TARGET :: a_srp    ! 10: Bulk precipitation mixing ratio (levels >= 4)
  TYPE(FloatArray3D), TARGET :: a_snrp   ! 11: Bulk precipitation number mixing ratio (levels >=4) 
  TYPE(FloatArray3D), TARGET :: a_rh     ! 12: Relative humidity
  TYPE(FloatArray3D), TARGET :: a_rsl    ! 13: Water vapor saturation mixing ratio
  TYPE(FloatArray3D), TARGET :: a_rhi    ! 14: Relative humidity over ice
  TYPE(FloatArray3D), TARGET :: a_rsi    ! 15: Water vapor saturation mixing ratio over ice
  TYPE(FloatArray3D), TARGET :: a_dn     ! 16: Air density
  TYPE(FloatArray3D), TARGET :: a_rflx, a_sflx,   &  ! 17, 18: Radiation fluxes
                                a_fus, a_fds,     &  ! 19, 20: 
                                a_fuir, a_fdir       ! 21, 22:
  TYPE(FloatArray3D), TARGET :: a_rrate              ! 23: Precipitation flux
  TYPE(FloatArray3D), TARGET :: a_irate              ! 24: Precipitation flux, frozen

  REAL, ALLOCATABLE, TARGET :: a_diag3d(:,:,:,:) 
  INTEGER, PARAMETER :: ndiag3d = 24   ! Remember to update if adding new variables!!

  !-------------------------------------------------------------------
  ! Binned diagnostic variables mainly for output
  !
  TYPE(FloatArray4d), TARGET :: d_VtPrc, d_VtIce ! Precipitation and ice terminal fall velocities
  REAL, ALLOCATABLE, TARGET :: d_binned(:,:,:,:)
  
  !----------------------------------------------------------------------------
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

  ! -------------------------------------------------------------------------------
  ! Do these need to be stored? If not, move to mo_derived_state?
  TYPE(FloatArray2D), TARGET :: a_sfcrrate           ! 10: Surface rain rate
  TYPE(FloatArray2D), TARGET :: a_sfcirate           ! 11: Surface frozen precipitation
  
  REAL, ALLOCATABLE, TARGET :: a_diag2d(:,:,:)
  INTEGER, PARAMETER :: ndiag2d = 11  ! Remember to update if adding new variables!!
  
  ! --------------------------------------------------------------------------------------------------------------------------------------------------
  ! Microphysical process rates from SALSA -- they need to be stored during the timestep because they cannot be simply diagnosed afterwards
  ! For now, these are BULK process rates only for water/ice, except where indicated otherwise !! Number concentration rate given for particle formation processes,
  ! mass concetration rate for others
  !
  TYPE(FloatArray3d), TARGET :: s_m_autoc        ! 1: Bulk autoconversion rate (mass, following bin regime limits)
  TYPE(FloatArray3d), TARGET :: s_m_autoc80      ! 2: Bulk autoconversion rate (mass, for drops past 80um)
  TYPE(FloatArray3d), TARGET :: s_m_autoc50      ! 3: Bulk autoconversion rate (mass, for drops past 50 um)
  TYPE(FloatArray3d), TARGET :: s_m_accr         ! 4: Bulk accretion rate (mass)
  TYPE(FloatArray3d), TARGET :: s_m_accr80       ! 5: Bulk accretion rate (mass, by drops L.T. 80um)
  TYPE(FloatArray3d), TARGET :: s_m_accr50       ! 6: Bulk accretion rate (mass, by drops L.T. 50 um)
  TYPE(FloatArray3d), TARGET :: s_m_ACcoll_dry   ! 7: Bulk cloud collection of aerosol (dry aerosol mass)
  TYPE(FloatArray3d), TARGET :: s_m_APcoll_dry   ! 8: Bulk precipitation collection of aerosol (dry aerosol mass)            
  TYPE(FloatArray3d), TARGET :: s_m_AIcoll_dry   ! 9: Bulk ice collection of aerosol (dry aerosol mass)
  TYPE(FloatArray3d), TARGET :: s_n_activ        ! 10: Bulk cloud activation rate (number)
  TYPE(FloatArray3d), TARGET :: s_n_icehom       ! 11: Bulk homogeneous freezing rate (number)
  TYPE(FloatArray3d), TARGET :: s_n_icedep       ! 12: Bulk deposition freezing rate (number)
  TYPE(Floatarray3d), TARGET :: s_n_iceimm       ! 13: Bulk immersion freezing rate (number)
  TYPE(FloatArray3d), TARGET :: s_m_conda        ! 14: Bulk condensation rate of water on aerosol (mass)
  TYPE(FloatArray3d), TARGET :: s_m_condc        ! 15: Bulk condensation rate of water on cloud droplets (mass)
  TYPE(FloatArray3d), TARGET :: s_m_condp        ! 16: Bulk condensation rate of water on precipitation (mass)
  TYPE(FloatArray3d), TARGET :: s_m_condi        ! 17: Bulk condensation (deposition) rate of water on ice (mass)
  INTEGER, PARAMETER :: nratediag3d_salsa = 17
  
  ! Microphysical process rates from bulk microphysics.
  TYPE(FloatArray3d), TARGET :: b_m_autoc        ! 1: Autoconversion rate (mass)
  TYPE(FloatArray3d), TARGET :: b_n_autoc        ! 2: Autoconversion rate (number)
  TYPE(FloatArray3d), TARGET :: b_m_accr         ! 2: Accretion rate (mass)
  INTEGER, PARAMETER :: nratediag3d_bulk = 3
  
  REAL, ALLOCATABLE, TARGET :: a_rateDiag3d(:,:,:,:)

  ! --------------------------------------------------------------------
  ! Diagnostic variables for piggybacking, used in the "slave" bulk scheme
  !
  TYPE(FloatArray3D), TARGET :: pb_theta,       & ! 1: Potential temp
                                pb_temp,        & ! 2: Absolute temp
                                pb_rv,          & ! 3: Water vapor mix rat
                                pb_rc,          & ! 4: Cloud + precip mix rat (similar to level 3)
                                pb_rh,          & ! 5: Relative humidity
                                pb_rsl,         & ! 6: Saturation mixing ratio
                                pb_rrate          ! 7: Precip flux

  TYPE(FloatArray2d), TARGET :: pb_sfcrrate       ! 1: Surface precip flux; Could this be also as derived on-demand diagnostic?
  REAL, ALLOCATABLE, TARGET :: pb_diag3d(:,:,:,:)
  REAL, ALLOCATABLE, TARGET :: pb_diag2d(:,:,:)
  INTEGER, PARAMETER :: npbdiag3d = 7
  INTEGER, PARAMETER :: npbdiag2d = 1

  
  CONTAINS

    SUBROUTINE setDiagnosticVariables(Diag,outputlist,memsize,level,iradtyp,lpback,nzp,nxp,nyp)
      TYPE(FieldArray), INTENT(inout) :: Diag
      CHARACTER(len=*), INTENT(in) :: outputlist(:)
      INTEGER, INTENT(inout) :: memsize
      INTEGER, INTENT(in) :: level, iradtyp, nzp,nxp,nyp
      LOGICAL, INTENT(in) :: lpback
      CLASS(*), POINTER :: pipeline => NULL()
      INTEGER :: nxyz, nxy, n2d,n3d,nr3d,n4db, npb3d, npb2d
      INTEGER :: nbinned

      nxyz = nxp*nyp*nzp
      nxy = nxp*nyp
      
      ALLOCATE(a_diag3d(nzp,nxp,nyp,ndiag3d))
      ALLOCATE(a_diag2d(nxp,nyp,ndiag2d))
      a_diag3d = 0.
      a_diag2d = 0.
      n3d = 0
      n2d = 0

      IF (level < 4) THEN
         ALLOCATE(a_rateDiag3d(nzp,nxp,nyp,nratediag3d_bulk))
      ELSE IF (level >= 4 .AND. .NOT. lpback) THEN
         ALLOCATE(a_rateDiag3d(nzp,nxp,nyp,nratediag3d_salsa))
      ELSE IF (level >=4 .AND. lpback) THEN
         ALLOCATE(a_rateDiag3d(nzp,nxp,nyp,nratediag3d_bulk+nratediag3d_salsa))
      END IF         
      a_rateDiag3d = 0.
      nr3d = 0
      
      nbinned = 0
      IF ( level >= 4) THEN
         IF ( level >= 4) nbinned = nbinned + nprc
         IF ( level > 4) nbinned = nbinned + nice
         ALLOCATE(d_binned(nzp,nxp,nyp,nbinned))
         d_binned = 0.
         n4db = 0
      END IF
         
      IF (lpback) THEN
         ALLOCATE(pb_diag3d(nzp,nxp,nyp,npbdiag3d))
         pb_diag3d = 0.
         npb3d = 0
         ALLOCATE(pb_diag2d(nxp,nyp,npbdiag2d))
         pb_diag2d = 0.
         npb2d = 0
      END IF

      ! -----------------------------------
      ! General diagnostic variables
      
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

      IF (level >= 4) THEN
         memsize = memsize + nxyz
         n3d = n3d + 1
         pipeline => NULL()
         a_rtot = FloatArray3d(a_diag3d(:,:,:,n3d))
         pipeline => a_rtot
         CALL Diag%newField("rtot", "Total waer mixing ratio", "kg/kg", "tttt",  &
                            ANY(outputlist == "rtot"), pipeline)
      END IF
         
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

      IF (level >= 2) THEN
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
      CALL Diag%newField("uw_sfc", "Vertical momentum flux with u wind", "m2 s-2", "xtytt",   &
                         ANY(outputlist == "uw_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      vw_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => vw_sfc
      CALL Diag%newField("vw_sfc", "Vertical momentum flux with v wind", "m2 s-2", "xtytt",    &
                         ANY(outputlist == "vw_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      ww_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => ww_sfc
      CALL Diag%newField("ww_sfc", "Vertical wind covariance", "m2 s-2", "xtytt",    &
                         ANY(outputlist == "ww_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      wt_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => wt_sfc
      CALL Diag%newField("wt_sfc", "Vertical temperature flux", "m K s-1", "xtytt",   &
                         ANY(outputlist == "wt_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      wq_sfc = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => wq_sfc
      CALL Diag%newField("wq_sfc", "Vertical moisture flux", "m Kg Kg-1 s-1", "xtytt",    &
                         ANY(outputlist == "wq_sfc"), pipeline) 

      memsize = memsize + nxy
      n2d = n2d+1
      pipeline => NULL()
      a_sfcrrate = FloatArray2d(a_diag2d(:,:,n2d))
      pipeline => a_sfcrrate
      CALL Diag%newField("sfcrrate", "Surface rain rate", "W m-2", "xtytt",     &
                         ANY(outputlist == "sfcrrate"), pipeline)

      IF ( level == 5 ) THEN
         memsize = memsize + nxy
         n2d = n2d+1
         pipeline => NULL()
         a_sfcirate = FloatArray2d(a_diag2d(:,:,n2d))
         pipeline => a_sfcirate
         CALL Diag%newField("sfcirate", "Surface frozen precip", "W m-2", "xtytt",    &
                            ANY(outputlist == "sfcirate"), pipeline)
      END IF

      ! -----------------------------------
      ! Process rate diagnostics

      IF (level >= 4) THEN
         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_autoc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_autoc
         CALL Diag%newField("s_m_autoc", "Autoconversion rate, h2o mass", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_autoc"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_autoc50 = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_autoc50
         CALL Diag%newField("s_m_autoc50", "Autoconversion rate, h2o mass, over 50um", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_autoc50"), pipeline)
         
         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_autoc80 = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_autoc80
         CALL Diag%newField("s_m_autoc80", "Autoconversion rate, h2o mass, over 80um", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_autoc80"), pipeline)
         
         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_accr = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_accr
         CALL Diag%newField("s_m_accr", "Accretion rate, h2o mass", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_accr"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_accr50 = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_accr50
         CALL Diag%newField("s_m_accr50", "Accretion rate, h2o mass, over 50um", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_accr50"), pipeline)
         
         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_accr80 = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_accr80
         CALL Diag%newField("s_m_accr80", "Accretion rate, h2o mass, over 80um", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_accr80"), pipeline)
         
         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_ACcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_ACcoll_dry
         CALL Diag%newField("s_m_ACcoll_dry", "Cloud collection of aerosol, dry mass", "kg/kgs", "tttt",  &
                            ANY(outputlist == "s_m_ACcoll_dry"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_APcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_APcoll_dry
         CALL Diag%newField("s_m_APcoll_dry","Rain collection of aerosol, dry mass", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_APcoll_dry"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_AIcoll_dry = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_AIcoll_dry
         CALL Diag%newField("s_m_AIcoll_dry","Ice collection of aerosol, dry mass", "kg/kgs", "tttt",   &
                            ANY(outputlist == "s_m_AIcoll_dry"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_n_activ = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_n_activ
         CALL Diag%newField("s_n_activ","Cloud activation rate, number", "#/kgs","tttt",   &
                            ANY(outputlist == "s_n_activ"), pipeline)
      
         nr3d = nr3d + 1
         pipeline => NULL()
         s_n_icehom = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_n_icehom
         CALL Diag%newfield("s_n_icehom","Homogeneous freezing rate, number", "#/kgs","tttt",   &
                            ANY(outputlist == "s_n_icehom"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_n_icedep = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_n_icedep
         CALL Diag%newField("s_n_icedep","Deposition freezing rate, number", "#/kgs", "tttt",  &
                            ANY(outputlist == "s_n_icedep"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_n_iceimm = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_n_iceimm
         CALL Diag%newField("s_n_iceimm","Immersion freezing rate, number", "#/kgs", "tttt",  &
                            ANY(outputlist == "s_n_iceimm"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_conda = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_conda
         CALL Diag%newField("s_m_conda","H2O condensation rate on aerosol", "#/kgs", "tttt",  &
                            ANY(outputlist == "s_m_conda"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_condc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_condc
         CALL Diag%newField("s_m_condc","H2O condensation rate on cloud droplets", "#/kgs", "tttt",  &
                            ANY(outputlist == "s_m_condc"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_condp = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_condp
         CALL Diag%newField("s_m_condp","H2O condensation rate on precipiation", "#/kgs", "tttt",    &
                            ANY(outputlist == "s_m_condp"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         s_m_condi = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => s_m_condi
         CALL Diag%newField("s_m_condi","H2O condensation rate on ice", "#/kgs", "tttt",   &
                            ANY(outputlist == "s_m_condi"), pipeline)

         IF (lpback) THEN
            nr3d = nr3d + 1
            pipeline => NULL()
            b_m_autoc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
            pipeline => b_m_autoc
            CALL Diag%newField("b_m_autoc", "Autoconversion rate, bulk scheme (mass)", "kg/kg s",  &
                               "tttt", ANY(outputlist == "b_m_autoc"), pipeline)

            nr3d = nr3d + 1
            pipeline => NULL()
            b_n_autoc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
            pipeline => b_n_autoc
            CALL Diag%newField("b_n_autoc", "Autoconversion rate, bulk scheme (number)", "#/kg s",  &
                               "tttt", ANY(outputlist == "b_n_autoc"), pipeline)

            nr3d = nr3d + 1
            pipeline => NULL()
            b_m_accr = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
            pipeline => b_m_accr
            CALL Diag%newField("b_m_accr", "Accretion rate, bulk scheme (mass)", "kg/kg s",  &
                               "tttt", ANY(outputlist == "b_m_accr"), pipeline)
            
         END IF         

      ELSE IF (level <= 3) THEN
         nr3d = nr3d + 1
         pipeline => NULL()
         b_m_autoc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => b_m_autoc
         CALL Diag%newField("b_m_autoc", "Autoconversion rate, bulk scheme (mass)", "kg/kg s",  &
                            "tttt", ANY(outputlist == "b_m_autoc"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         b_n_autoc = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => b_n_autoc
         CALL Diag%newField("b_n_autoc", "Autoconversion rate, bulk scheme (number)", "#/kg s",  &
                            "tttt", ANY(outputlist == "b_n_autoc"), pipeline)

         nr3d = nr3d + 1
         pipeline => NULL()
         b_m_accr = FloatArray3d(a_rateDiag3d(:,:,:,nr3d))
         pipeline => b_m_accr
         CALL Diag%newField("b_m_accr", "Accretion rate, bulk scheme (mass)", "kg/kg s",  &
                            "tttt", ANY(outputlist == "b_m_accr"), pipeline)
            
      END IF      
        
      ! Binned variables
      n4db = 1
      IF ( level >= 4) THEN
         pipeline => NULL()
         d_VtPrc = FloatArray4d(d_binned(:,:,:,n4db:n4db+nprc-1))
         pipeline => d_VtPrc
         CALL Diag%newField("VtPrc", "Terminal fall speed of precip", "m/s", "ttttprc",   &
                            ANY(outputlist == "VtPrc"), pipeline)
         n4db = n4db + nprc
      END IF

      IF (level > 4) THEN
         pipeline => NULL()
         d_VtIce = FloatArray4d(d_binned(:,:,:,n4db:n4db+nice-1))
         pipeline => d_VtIce
         CALL Diag%newField("VtIce", "Terminal fall speed of ice", "m/s", "ttttice",   &
                            ANY(outputlist == "VtIce"), pipeline)
         n4db = n4db + nice
      END IF

      ! Piggybacking variables for "slave" microphysics
      IF (lpback) THEN
         npb3d = npb3d + 1
         pipeline => NULL()
         pb_theta = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_theta
         CALL Diag%newField("pb_theta", "Potential temperature (bulk slave)", "K", "tttt",   &
                            ANY(outputlist == "pb_theta"), pipeline)
         
         npb3d = npb3d + 1
         pipeline => NULL()
         pb_temp = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_temp
         CALL Diag%newField("pb_temp", "Absolute temperature (bulk slave)", "K", "tttt",    &
                            ANY(outputlist == "pb_temp"), pipeline)         

         npb3d = npb3d + 1
         pipeline => NULL()
         pb_rv = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_rv
         CALL Diag%newField("pb_rv", "Water vapor mixing ratio (bulk slave)", "kg/kg", "tttt",   &
                            ANY(outputlist == "pb_rv"), pipeline)

         npb3d = npb3d + 1
         pipeline => NULL()
         pb_rc = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_rc
         CALL Diag%newField("pb_rc", "Condensate mixing ratio (bulk slave)", "kg/kg", "tttt",    &
                            ANY(outputlist == "pb_rc"), pipeline)

         npb3d = npb3d + 1
         pipeline => NULL()
         pb_rh = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_rh
         CALL Diag%newField("pb_rh", "Relative humidity (bulk slave)", "1", "tttt",    &
                            ANY(outputlist == "pb_rh"), pipeline)
         
         npb3d = npb3d + 1
         pipeline => NULL()
         pb_rsl = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_rsl
         CALL Diag%newField("pb_rsl", "Saturation mixing ratio (bulk slave)", "kg/kg", "tttt",    &
                            ANY(outputlist == "pb_rsl"), pipeline)

         npb3d = npb3d + 1
         pipeline => NULL()
         pb_rrate = FloatArray3d(pb_diag3d(:,:,:,npb3d))
         pipeline => pb_rrate
         CALL Diag%newField("pb_rrate", "Rain flux (bulk slave)", "W m-2", "tttt",     &
                            ANY(outputlist == "pb_rrate"), pipeline)

         npb2d = 1
         pipeline => NULL()
         pb_sfcrrate = FloatArray2d(pb_diag2d(:,:,npb2d))
         pipeline => pb_sfcrrate
         CALL Diag%newField("pb_sfcrrate", "Surface rain flux (bulk slave)", "W m-2", "xtytt",   &
                            ANY(outputlist == "pb_sfcrrate"), pipeline)
         
      END IF
      
      pipeline => NULL()
      
    END SUBROUTINE setDiagnosticVariables
        
END MODULE mo_diag_state
