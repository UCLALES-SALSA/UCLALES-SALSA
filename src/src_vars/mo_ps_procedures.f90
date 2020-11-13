MODULE mo_ps_procedures
  USE grid, ONLY : nzp,nxp,nyp, dtlt, level
  USE mo_stats_finder
  USE mo_structured_datatypes, ONLY : FloatArray1d
  USE defs
  IMPLICIT NONE

  ! Calculate profile statistics. The simple domain mean profile is automatically
  ! available for all 3d (including binned) fields defined during model runtime
  ! (found in the FieldArray instances Prog, Vector, Diag or Derived).
  
  ! Higher-order statistics are yet to be implemented, but will follow the exact same logic.
  ! Conditional statistics also require specific implementations.
  
  PRIVATE

  PUBLIC :: globalMinProfile, globalMaxProfile,                           & 
            globalMeanProfile, globalMeanProfileBinned, globalVarProfile, &
            incloudMeanProfile, precipMeanProfile,                        &
            meanScalarFlux, meanMomFlux, meanTKEres

  CONTAINS

    SUBROUTINE globalMinProfile(SELF,output)
      USE util, ONLY : get_gmin
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: fvar(nzp,nxp,nyp)
      fvar = 0.
      output = 0.
      CALL stats_get(SELF%srcName,fvar)
      CALL get_gmin(nzp,nxp,nyp,fvar,output)
    END SUBROUTINE globalMinProfile

    ! ----------------------------------------------------

    SUBROUTINE globalMaxProfile(SELF,output)
      USE util, ONLY : get_gmax
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: fvar(nzp,nxp,nyp)
      fvar = 0.
      output = 0.
      CALL stats_get(SELF%srcName,fvar)
      CALL get_gmax(nzp,nxp,nyp,fvar,output)
    END SUBROUTINE globalMaxProfile

    ! -------------------------------------------------------
    
    SUBROUTINE globalMeanProfile(SELF,output)
      USE util, ONLY : get_avg3_root
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)      
      REAL :: fvar(nzp,nxp,nyp)
      fvar = 0.
      output = 0.
      CALL stats_get(SELF%srcName,fvar)      
      CALL get_avg3_root(nzp,nxp,nyp,fvar,output)            
    END SUBROUTINE globalMeanProfile

    ! ----------------------------------------------------

    SUBROUTINE globalVarProfile(SELF,output)
      USE util, ONLY : get_var3_root
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: fvar(nzp,nxp,nyp)      
      fvar = 0.
      output = 0.
      CALL stats_get(SELF%srcName,fvar)
      CALL get_var3_root(nzp,nxp,nyp,fvar,output)
    END SUBROUTINE globalVarProfile
      
    ! ----------------------------------------------------
    
    SUBROUTINE globalMeanProfileBinned(name,output,nstr,nend)
      USE util, ONLY : get_avg3_binned_root
      CHARACTER(len=*), INTENT(in) :: name
      INTEGER, INTENT(in) :: nstr,nend
      REAL, INTENT(out) :: output(nzp,nend-nstr+1) 
      REAL :: fvar(nzp,nxp,nyp,nend-nstr+1)     
      fvar = 0.
      output = 0.
      CALL stats_get3d_binned(name,fvar,nstr,nend)      
      CALL get_avg3_binned_root(nzp,nxp,nyp,nend-nstr+1,fvar,output)
    END SUBROUTINE globalMeanProfileBinned

    ! -----------------------------------------------------
    
    SUBROUTINE incloudMeanProfile(SELF,output)
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: fvar(nzp,nxp,nyp)
      fvar = 0.
      output = 0.
    END SUBROUTINE incloudMeanProfile

    ! -----------------------------------------------------

    SUBROUTINE precipMeanProfile(SELF,output)
      CLASS(FloatArray1d), INTENT(in) :: SELF 
      REAL, INTENT(out) :: output(:)
      output = 0.
    END SUBROUTINE precipMeanProfile

    ! -----------------------------------------------------    

    SUBROUTINE meanScalarFlux(SELF,output)
      USE advf, ONLY : mamaos
      USE mo_aux_state, ONLY : dzt, dzm, dn0
      USE util, ONLY : get_avg3_root
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)

      REAL :: wp(nzp,nxp,nyp)
      REAL :: fvar(nzp,nxp,nyp)
      REAL :: tmp1(nzp,nxp,nyp)
      tmp1 = 0.  ! needed just as a dummy argument for the output diagnostic...
      output = 0.
      fvar = 0.
      
      CALL stats_get("wwind",wp)
      
      SELECT CASE(SELF%shortName)
      CASE("tw_res")
         CALL stats_get("theta",fvar)
         tmp1 = fvar
         CALL mamaos(nzp,nxp,nyp,wp,fvar,tmp1,dzt,dzm,dn0,dtlt,.FALSE.)        
         wp = wp * cp
      CASE("qw_res")
         IF (level < 4) THEN
            CALL stats_get("rv",fvar)
         ELSE
            CALL stats_get("rp",fvar)
         END IF
         tmp1 = fvar
         CALL mamaos(nzp,nxp,nyp,wp,fvar,tmp1,dzt,dzm,dn0,dtlt,.FALSE.)
         wp = wp * alvl
      CASE("lw_res")
         CALL stats_get("rc",fvar)
         IF (level < 4) THEN
            CALL stats_get("rpp",tmp1)
         ELSE
            CALL stats_get("srp",tmp1)
         END IF
         fvar = fvar + tmp1   ! Total liquid flux           
         tmp1 = fvar
         CALL mamaos(nzp,nxp,nyp,wp,fvar,tmp1,dzt,dzm,dn0,dtlt,.FALSE.)
         wp = wp * alvl  
      END SELECT
      
      CALL get_avg3_root(nzp,nxp,nyp,wp,output)
      
    END SUBROUTINE meanScalarFlux
    
    ! --------------------

    SUBROUTINE meanMomFlux(SELF,output)
      USE advl, ONLY : ladvzu, ladvzv, ladvzw
      USE mo_aux_state, ONLY : dn0
      USE util, ONLY : get_avg3_root
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: wpm(nzp,nxp,nyp)
      REAL :: tmp1(nzp),tmp2(nzp,nxp,nyp)
      REAL :: flx(nzp,nxp,nyp),fvar(nzp,nxp,nyp)
      INTEGER :: k
      tmp1 = 1.
      tmp2 = 0.
            
      CALL stats_get("wwind",wpm)
      wpm = 0.
      DO k = 2,nzp-1
         wpm(k,:,:) = wpm(k,:,:)*0.5*(dn0%d(k)+dn0%d(k+1))
      END DO
                     
      SELECT CASE(SELF%shortName)
      CASE("uw_res")   
         CALL stats_get("uwind",fvar)    ! This takes always the p-value although here it really should be the c value...
         CALL ladvzu(nzp,nxp,nyp,fvar,tmp2,wpm,flx,tmp1)
      CASE("vw_res")
         CALL stats_get("vwind",fvar)
         CALL ladvzv(nzp,nxp,nyp,fvar,tmp2,wpm,flx,tmp1)
      CASE("ww_res")
         CALL stats_get("wwind",fvar)
         CALL ladvzw(nzp,nxp,nyp,fvar,tmp2,wpm,flx,tmp1)
      END SELECT

      CALL get_avg3_root(nzp,nxp,nyp,flx,output)
         
    END SUBROUTINE meanMomFlux
    
    ! -------------------------------

    SUBROUTINE meanTKEres(SELF,output)
      USE mo_vector_state, ONLY : a_up, a_vp, a_wp
      USE util, ONLY : get_var3_root
      CLASS(FloatArray1d), INTENT(in) :: SELF
      REAL, INTENT(out) :: output(:)
      REAL :: uvar(nzp), vvar(nzp), wvar(nzp)      
      CALL get_var3_root(nzp,nxp,nyp,a_up%d,uvar)
      CALL get_var3_root(nzp,nxp,nyp,a_vp%d,vvar)
      CALL get_var3_root(nzp,nxp,nyp,a_wp%d,wvar)
      output = 0.5 * ( uvar + vvar + wvar )            
    END SUBROUTINE meanTKEres

    
END MODULE mo_ps_procedures
