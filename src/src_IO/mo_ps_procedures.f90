MODULE mo_ps_procedures
  USE grid, ONLY : nzp,nxp,nyp
  USE mo_stats_finder
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: globalAvgProfile
  
  CONTAINS

    SUBROUTINE globalAvgProfile(name,output)
      USE util, ONLY : get_avg3_root
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output(nzp)
      
      REAL :: fvar(nzp,nxp,nyp)
      
      fvar = 0.
      output = 0.
      CALL stats_get3d(name,fvar)
      
      CALL get_avg3_root(nzp,nxp,nyp,fvar,output)      
      
    END SUBROUTINE globalAvgProfile
    

END MODULE mo_ps_procedures
