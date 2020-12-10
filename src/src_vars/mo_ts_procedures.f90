MODULE mo_ts_procedures
  USE grid, ONLY : deltax,deltay,nzp,nxp,nyp,level
  USE mo_aux_state, ONLY : dzt,zt
  USE mo_stats_finder
  USE mo_structured_datatypes, ONLY : FloatArray0d
  USE mo_stats_parameters
  IMPLICIT NONE

  PRIVATE
  
  PUBLIC :: tsMin,tsMax,tsSfcMean,tsInCloudMean,               &
            tsInLiqMean,tsInIceMean,tsInPrecipMean,            &
            tsCloudBoundaries,tsCloudFraction,    &
            tsTotalMass

  CONTAINS

    ! Check if variable exists in the FieldArray instances, then check if the data is stored
    ! or if it is associated with an onDemand procedure. Finally calculate the statistics.

    !
    ! -------------------------------------------------
    ! SUBROUTINE tsMin:
    ! Gets a simple global minimum
    ! To use this procedure, the srcName must be defined in the
    ! FloatArray instance!
    !
    SUBROUTINE tsMin(SELF,output)
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      CHARACTER(len=50) :: dim
      CALL stats_getDim(SELF%srcName,dim)   ! Determine the dimensions of source variable
      SELECT CASE(dim)
      CASE ('tttt','mttt','tmtt','ttmt')
         CALL ts3dMin(SELF%srcName,output)
      CASE ('xtytt')
         CALL ts2dMin(SELF%srcName,output)
      END SELECT
    END SUBROUTINE tsMin
    ! -------------------------------------------
    SUBROUTINE ts3dMin(iname,output)
      USE util, ONLY : get_gmin
      CHARACTER(len=*), INTENT(in) :: iname  ! Name of the variable to be processed
      REAL, INTENT(out) :: output            ! The result min value          
      REAL :: fvar(nzp,nxp,nyp)              ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      CALL get_gmin(nzp,nxp,nyp,fvar,output)       
    END SUBROUTINE ts3dMin
    ! -------------------------------------------
    SUBROUTINE ts2dMin(iname,output)
      USE util, ONLY : get_gmin
      CHARACTER(len=*), INTENT(in) :: iname  ! Name of the variable to be processed
      REAL, INTENT(out) :: output            ! the result min value
      REAL :: fvar(nxp,nyp)                  ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      CALL get_gmin(nxp,nyp,fvar,output)       
    END SUBROUTINE ts2dMin
    
    !
    ! -------------------------------------------------
    ! SUBROUTINE tsMax:
    ! Gets a simple global maximum 
    ! To use this procedure, the srcName must be defined in the
    ! FloatArray instance!
    !
    SUBROUTINE tsMax(SELF,output)
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      CHARACTER(len=50) :: dim
      CALL stats_getDim(SELF%srcName,dim)   ! Determine the dimensions of the source variable
      SELECT CASE(dim)
      CASE ('tttt','mttt','tmtt','ttmt')
         CALL ts3dMax(SELF%srcName,output)
      CASE ('xtytt')
         CALL ts2dMax(SELF%srcName,output)
      END SELECT
    END SUBROUTINE tsMax
    ! -------------------------------------------
    SUBROUTINE ts3dMax(iname,output)
      USE util, ONLY : get_gmax
      CHARACTER(len=*), INTENT(in) :: iname    ! Name of the variable to be processed
      REAL, INTENT(out) :: output              ! The result max value
      REAL :: fvar(nzp,nxp,nyp)                ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      CALL get_gmax(nzp,nxp,nyp,fvar,output)       
    END SUBROUTINE ts3dMax
    ! -------------------------------------------
    SUBROUTINE ts2dMax(iname,output)
      USE util, ONLY : get_gmax
      CHARACTER(len=*), INTENT(in) :: iname    ! Name of the variable to be processed
      REAL, INTENT(out) :: output              ! The result max value
      REAL :: fvar(nxp,nyp)                    ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      CALL get_gmax(nxp,nyp,fvar,output)       
    END SUBROUTINE ts2dMax

    !
    ! ---------------------------------------------------------
    ! SUBROUTINE tsSfcMean:
    ! Gets surface mean values. If the source variable is 2d
    ! in the x-y plane, it is assumed to be a surface variable!!
    ! To use this procedure, the srcName must be defined in the
    ! FloatArray instance!
    !
    SUBROUTINE tsSfcMean(SELF,output)
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      CHARACTER(len=50) :: dim
      CALL stats_getDim(SELF%srcName,dim)      ! Determine the dimensions of the source variable
      SELECT CASE(dim)
      CASE ('tttt','mttt','tmtt','ttmt')
         CALL ts3dSfcMean(SELF%srcName,output)
      CASE ('xtytt')
         CALL ts2dSfcMean(SELF%srcName,output)
      END SELECT
    END SUBROUTINE
    ! --------------
    SUBROUTINE ts3dSfcMean(iname,output)
      USE util, ONLY : get_avg2dh
      CHARACTER(len=*), INTENT(in) :: iname     ! Name of the variable to be processed
      REAL, INTENT(out) :: output               ! The result mean value
      REAL :: fvar(nzp,nxp,nyp)                 ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      output = get_avg2dh(nxp,nyp,fvar(2,:,:))      
    END SUBROUTINE ts3dSfcMean
    ! ---------------
    SUBROUTINE ts2dSfcMean(iname,output)
      USE util, ONLY : get_avg2dh
      CHARACTER(len=*), INTENT(in) :: iname     ! Name of the variable to be processed
      REAL, INTENT(out) :: output               ! The result mean value
      REAL :: fvar(nxp,nyp)                     ! Values of the variable to be processed
      fvar = 0.
      output = 0.      
      CALL stats_get(iname,fvar)      
      output = get_avg2dh(nxp,nyp,fvar(:,:))      
    END SUBROUTINE ts2dSfcMean

    !
    ! ----------------------------------------------
    ! SUBROUTINE tsInCloudMean:
    ! Get conditionally averaged mean values inside a cloud volume.
    ! If the source variable is 2d in x-y plane, it is assumed to
    ! be a surface variable or e.g. LWP.
    ! 
    SUBROUTINE tsInCloudMean(SELF,output)
      USE mo_diag_state, ONLY : a_rc, a_ri, a_riri
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: cond(nzp,nxp,nyp)
      cond = ( a_rc%d > TH_rc .OR. a_ri%d+a_riri%d > TH_ri )
      CALL ts3dCondMean(SELF%srcName,cond,output)
    END SUBROUTINE tsInCloudMean
    !
    ! -----------------------------------------------
    ! SUBROUTINE tsInLiqMean
    ! Get conditionally averaged mean values inside a liquid
    ! cloud volume. If the source variable is 2d in x-y plane,
    ! it is assumed to be a surface variable or e.g. LWP.
    !
    SUBROUTINE tsInLiqMean(SELF,output)
      USE mo_diag_state, ONLY : a_rc, a_dn
      USE mo_derived_state, ONLY : lwp
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: cond(nzp,nxp,nyp)
      REAL :: zlwp(nxp,nyp)
      CHARACTER(len=50) :: dim
      cond = .FALSE.
      CALL stats_getDim(SELF%srcName,dim)
      SELECT CASE(dim)
      CASE('tttt','mttt','tmtt','ttmt')      
         cond = ( a_rc%d > TH_rc )
         CALL ts3dCondMean(SELF%srcName,cond,output)
      CASE('xtytt')
         CALL lwp%onDemand('lwp',zlwp)
         cond(2,:,:) = ( zlwp > TH_rc*a_dn%d(2,:,:)*MAXVAL(1./dzt%d) )
         CALL ts2dCondMean(SELF%srcName,cond(2,:,:),output)
      END SELECT         
    END SUBROUTINE tsInLiqMean
    !
    ! ------------------------------------------------
    ! SUBROUTINE tsInIceMean
    ! Get conditionally averaged values inside an ice
    ! cloud volume. If the source variable is 2d in x-y plane,
    ! it is assumed to be a surface variable or e.g. LWP.
    !
    SUBROUTINE tsInIceMean(SELF,output)
      USE mo_diag_state, ONLY : a_ri,a_riri, a_dn
      USE mo_derived_state, ONLY : iwp
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: cond(nzp,nxp,nyp)
      REAL :: ziwp(nxp,nyp)
      CHARACTER(len=50) :: dim
      cond = .FALSE.
      CALL stats_getDim(SELF%srcName,dim)
      SELECT CASE(dim)
      CASE('tttt','mttt','tmtt','ttmt')
         cond = ( a_ri%d+a_riri%d > TH_ri   )
         CALL ts3dCondMean(SELF%srcName,cond,output)
      CASE('xtytt')
         CALL iwp%onDemand('iwp',ziwp)
         cond(2,:,:) = (ziwp > TH_rc*a_dn%d(2,:,:)*MAXVAL(1./dzt%d))
         CALL ts2dCondMean(SELF%srcName,cond(2,:,:),output)
      END SELECT
    END SUBROUTINE tsInIceMean
    !
    ! --------------------------------------------------
    ! SUBROUTINE tsInPrecipMean
    ! Get conditionally averaged values in grid points
    ! containing drizzle or rain. If the source variable
    ! is 2d in x-y plane, it is assumed to be a surface
    ! variable or e.g. RWP
    !
    SUBROUTINE tsInPrecipMean(SELF,output)
      USE mo_diag_state, ONLY : a_rrate
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: cond(nzp,nxp,nyp)
      CHARACTER(len=50) :: dim
      cond(:,:,:) = .FALSE.
      cond = ( a_rrate%d > TH_rrate )
      CALL stats_getDim(SELF%srcName,dim)
      SELECT CASE(dim)
      CASE('tttt','mttt','tmtt','ttmt')
         CALL ts3dCondMean(SELF%srcName,cond,output)
      CASE('xtytt')
         CALL ts2dCondMean(SELF%srcName,cond(2,:,:),output)
      END SELECT      
    END SUBROUTINE tsInPrecipMean    
    ! ----------------
    SUBROUTINE ts2dCondMean(iname,cond,output)
      USE util, ONLY : get_avg2_root
      CHARACTER(len=*), INTENT(in) :: iname     ! Name of the variable to be processed
      LOGICAL, INTENT(in) :: cond(nxp,nyp)      ! Logical mask for conditional averaging
      REAL, INTENT(out) :: output               ! The result mean value
      REAL :: fvar(nxp,nyp)                     ! Values of the variable to be processed
      fvar = 0.
      output = 0.
      CALL stats_get(iname,fvar)
      CALL get_avg2_root(nxp,nyp,fvar,output,cond)
    END SUBROUTINE ts2dCondMean
    ! ---------------
    SUBROUTINE ts3dCondMean(iname,cond,output)
      USE util, ONLY : get_avg4_root
      CHARACTER(len=*), INTENT(in) :: iname     ! Name of the variable to be processed
      LOGICAL, INTENT(in) :: cond(nzp,nxp,nyp)  ! Logical mask for conditional averaging
      REAL, INTENT(out) :: output            ! the result mean value
      REAL :: fvar(nzp,nxp,nyp)                 ! Values of the variable to be processed
      fvar = 0.
      output = 0.
      CALL stats_get(iname,fvar)
      CALL get_avg4_root(nzp,nxp,nyp,fvar,output,dzt%d,cond)
    END SUBROUTINE ts3dCondMean

    !
    ! -------------------------------------------------
    ! SUBROUTINE tsCloudBoundaries
    ! Determines the cloud base and top heights.
    !
    SUBROUTINE tsCloudBoundaries(SELF,output)
      USE util, ONLY : get_avg3_root
      USE mo_diag_state, ONLY : a_rc, a_ri, a_riri
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: liqmask(nzp,nxp,nyp), icemask(nzp,nxp,nyp)
      LOGICAL :: profmask(nzp)
      REAL :: liqprof(nzp), iceprof(nzp)
      INTEGER :: k
      
      output = 0.; liqprof = 0.; iceprof = 0.
      liqmask = .FALSE.; icemask = .FALSE.; profmask = .FALSE.
      
      liqmask = ( a_rc%d > TH_rc )
      CALL get_avg3_root(nzp,nxp,nyp,a_rc%d,liqprof,liqmask) ! Avg in-cloud cloud condensate

      IF (level < 5) THEN
         profmask = ( liqprof > TH_rc )
      ELSE
         icemask = ( a_ri%d + a_riri%d > TH_ri )      
         CALL get_avg3_root(nzp,nxp,nyp,a_ri%d+a_riri%d,iceprof,icemask) ! Avg in-cloud cloud condensate
         profmask = ( liqprof > TH_rc .OR. iceprof > TH_ri )
      END IF

      IF ( .NOT. ANY(profmask) ) RETURN ! No clouds found
      
      SELECT CASE(SELF%shortName)
      CASE('ctop')
         DO k = nzp,1,-1
            IF (profmask(k)) THEN
               output = zt%d(k)
               EXIT
            END IF
         END DO
      CASE('cbase')
         DO k = 2,nzp
            IF (profmask(k)) THEN
               output = zt%d(k)
               EXIT
            END IF
         END DO
      END SELECT
      
    END SUBROUTINE tsCloudBoundaries

    !
    ! -------------------------------------------
    ! SUBROUTINE tsCloudFraction:
    ! Get the domain mean total cloud fraction.
    !
    SUBROUTINE tsCloudFraction(SELF,output)
      USE util, ONLY : get_avg2_root
      USE mo_diag_state, ONLY : a_rc,a_ri,a_riri
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      LOGICAL :: cldmask(nzp,nxp,nyp)
      REAL :: lateral(nxp,nyp)
      output = 0.
      lateral = 0.
      cldmask = .FALSE.
      cldmask =( a_rc%d > TH_rc )
      IF (level == 5) cldmask = ( cldmask .OR. a_ri%d + a_riri%d > TH_ri )     
      lateral(:,:) = MERGE(1.,0., ANY(cldmask, DIM=1) )      
      CALL get_avg2_root(nxp,nyp,lateral,output)      
    END SUBROUTINE tsCloudFraction

    !SUBROUTINE tsSurface(SELF,name,output)
    !  CLASS(FloatArray0d), INTENT(in) :: SELF
    !  CHARACTER(len=*), INTENT(in) :: name
    !  REAL, INTENT(out) :: output
    !  output = 0.
    !END SUBROUTINE tsSurface
    
    !
    ! -------------------------------------------
    ! SUBROUTINE tsTotalMass:
    ! Calculates the total mass of everything except air in the domain
    ! THIS ONLY WORK NOW FOR LEVEL >= 4
    !
    SUBROUTINE tsTotalMass(SELF,output)
      USE mo_progn_state, ONLY : a_maerop, a_mcloudp, a_mprecpp, a_micep, a_rp
      USE mo_diag_state, ONLY : a_dn
      CLASS(FloatArray0d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output
      REAL :: box(nzp)
      INTEGER :: i,j,k

      box(1:nzp) = deltax*deltay/dzt%d(1:nzp)
      
      output = 0.

      DO j = 1,nyp
         DO i = 1,nxp
            DO k = 1,nzp
               output = output + SUM(a_maerop%d(k,i,j,:))*box(k)*a_dn%d(k,i,j)
               output = output + SUM(a_mcloudp%d(k,i,j,:))*box(k)*a_dn%d(k,i,j)
               output = output + SUM(a_mprecpp%d(k,i,j,:))*box(k)*a_dn%d(k,i,j)
               IF (level == 5) output = output + SUM(a_micep%d(k,i,j,:))*box(k)*a_dn%d(k,i,j)
               output = output + a_rp%d(k,i,j)*box(k)*a_dn%d(k,i,j)
            END DO
         END DO
      END DO

      
    END SUBROUTINE tsTotalMass
    
      
END MODULE mo_ts_procedures
