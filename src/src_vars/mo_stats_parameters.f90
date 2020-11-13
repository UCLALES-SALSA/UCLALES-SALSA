MODULE mo_stats_parameters
  IMPLICIT NONE

  ! Contains some parameters and constants needed to calculate statistical outputs
  ! for ps and ts files. Including threshold values for conditional averaging etc.
  ! These can be specified using the 'output' namelist.
  
  REAL :: TH_rc = 1.e-5,   &     ! Threshold for cloud water kg/kg
          TH_ri = 1.e-6,   &     ! Threshold for cloud ice kg/kg
          TH_rr = 1.e-6,   &     ! Threshold for drizzle/rain kg/kg
          TH_rrate = 30.          ! Threshold for precipitation flux W m-2 


END MODULE mo_stats_parameters
