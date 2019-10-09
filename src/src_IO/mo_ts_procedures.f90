MODULE mo_ts_procedures
  USE grid, ONLY : nzp,nxp,nyp,level
  USE mo_aux_state, ONLY : dzt,zt
  USE mo_stats_finder
  IMPLICIT NONE

  PRIVATE

  PUBLIC :: globalMax, globalMean3d, globalAvgWaterPaths,    &
            globalAvgCloudBoundaries, globalAvgCloudFraction

  CONTAINS

    ! Check if variable exists in the FieldArray instances, then check if the data is stored
    ! or if it is associated with an onDemand procedure. Finally calculate the statistics.

    !
    ! Global max values
    ! 
    SUBROUTINE globalMax(name,output)
      USE util, ONLY : get_gmax
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output
      REAL :: fvar(nzp,nxp,nyp)      
      fvar = 0.
      output = 0.
      CALL stats_get3d(name,fvar)      
      CALL get_gmax(nzp,nxp,nyp,fvar,output)            
    END SUBROUTINE globalMax

    !
    ! Global mean from 3d
    !
    SUBROUTINE globalMean3d(name,output)
      USE util, ONLY : get_avg4_root
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output
      REAL :: fvar(nzp,nxp,nyp)      
      fvar = 0.
      output = 0.
      CALL stats_get3d(name,fvar)      
      CALL get_avg4_root(nzp,nxp,nyp,fvar,output,dzt%d)            
    END SUBROUTINE globalMean3d

    !
    ! Global average water paths
    !
    SUBROUTINE globalAvgWaterPaths(name,output)
      USE util, ONLY : get_avg3_root
      USE mo_diag_state, ONLY : a_rc, a_srp, a_ri, a_riri, a_dn
      USE mo_progn_state, ONLY : a_rpp
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output
      REAL :: field(nzp,nxp,nyp)
      REAL :: prof(nzp)
      prof = 0.
      field = 0.
      output = 0.
      SELECT CASE(name)
      CASE('lwp_bar')
         field = a_rc%d*a_dn%d
      CASE('rwp_bar')
         IF (level < 4) THEN
            field = a_rpp%d*a_dn%d
         ELSE
            field = a_srp%d*a_dn%d
         END IF
      CASE('iwp_bar')
         field = (a_riri%d+a_ri%d)*a_dn%d
      END SELECT
      CALL get_avg3_root(nzp,nxp,nyp,field,prof)
      output = SUM(prof/dzt%d)           
    END SUBROUTINE globalAvgWaterPaths

    !
    ! Global average cloud boundaries
    ! Take the global mean profile of cloud water and diagnose from there
    ! 
    SUBROUTINE globalAvgCloudBoundaries(name,output)
      USE util, ONLY : get_avg3_root
      USE mo_diag_state, ONLY : a_rc,a_ri,a_riri      
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output
      REAL :: field(nzp,nxp,nyp)
      REAL :: prof(nzp)
      LOGICAL :: cldmask(nzp)
      REAL, PARAMETER :: TH = 1.e-6
      INTEGER :: k
      field = 0.
      prof = 0.
      output = 0.
      IF (level < 5) THEN
         field = a_rc%d
      ELSE
         field = a_rc%d + a_ri%d + a_riri%d
      END IF         
      CALL get_avg3_root(nzp,nxp,nyp,field,prof)
      cldmask = ( prof > TH  )
      IF (ANY(cldmask)) THEN ! any clouds found?
         SELECT CASE(name)
         CASE('ctop')  ! Catch the first cloud top starting from model top
            DO k = nzp,1,-1
               IF (cldmask(k)) THEN
                  output = zt%d(k)
                  EXIT
               END IF
            END DO
         CASE('cbase') ! Catch the first cloud base starting from surface
            DO k = 2,nzp
               IF (cldmask(k)) THEN
                  output = zt%d(k)
                  EXIT
               END IF                  
            END DO
         END SELECT
      END IF      
    END SUBROUTINE globalAvgCloudBoundaries

    !
    ! Global average cloud fraction
    ! Mask individual columns where there is any cloud as 1, the rest as 0 and take the global mean
    ! 
    SUBROUTINE globalAvgCloudFraction(name,output)
      USE util, ONLY : get_avg2_root
      USE mo_diag_state, ONLY : a_rc,a_ri,a_riri
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(out) :: output
      REAL :: field(nzp,nxp,nyp)
      REAL :: lateral(nxp,nyp)
      REAL, PARAMETER :: TH = 1.e-6
      field = 0.
      output = 0.
      lateral = 0.
      field = a_rc%d
      IF (level == 5) field = field + a_riri%d + a_riri%d      
      lateral(:,:) = MERGE(1.,0., ANY( field(:,:,:) > TH, DIM=1 ))      
      CALL get_avg2_root(nxp,nyp,lateral,output)      
    END SUBROUTINE globalAvgCloudFraction
      
END MODULE mo_ts_procedures
