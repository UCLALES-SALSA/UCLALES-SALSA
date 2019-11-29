MODULE mo_ps_procedures
  USE grid, ONLY : nzp,nxp,nyp
  USE mo_stats_finder
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: globalAvgProfile, globalAvgProfileBinned
  
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

    SUBROUTINE globalAvgProfileBinned(name,output,nstr,nend)
      USE util, ONLY : get_avg3_binned_root
      CHARACTER(len=*), INTENT(in) :: name
      INTEGER, INTENT(in) :: nstr,nend
      REAL, INTENT(out) :: output(nzp,nend-nstr+1) 

      REAL :: fvar(nzp,nxp,nyp,nend-nstr+1)
     
      fvar = 0.
      output = 0.
      CALL stats_get3d_binned(name,fvar,nstr,nend)
      
      CALL get_avg3_binned_root(nzp,nxp,nyp,nend-nstr+1,fvar,output)

    END SUBROUTINE globalAvgProfileBinned

    

END MODULE mo_ps_procedures
