
!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE step

  USE mo_submctl, ONLY : spec, nice
  USE util, ONLY : getMassIndex, calc_correlation
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray3d, FloatArray4d ! poista vika
  
  IMPLICIT NONE
  
  INTEGER :: istpfl = 1
  REAL    :: timmax = 18000.
  LOGICAL :: corflg = .FALSE.
  
  REAL    :: frqhis =  9000.
  REAL    :: radfrq =  0.
  
  REAL    :: time   =  0.
  REAL    :: strtim =  0.0    ! In decimal days, 0.5 mid-day
  LOGICAL :: outflg = .TRUE.
  
CONTAINS
   !
   ! ----------------------------------------------------------------------
   ! Subroutine model:  This is the main driver for the model's time
   ! integration.  It calls the routine tstep, which steps through the
   ! physical processes active on a time-step and updates variables.  It
   ! then checks to see whether or not different output options are
   ! satisfied.
   SUBROUTINE stepper

      USE mpi_interface, ONLY : myid, double_scalar_par_max
      USE mo_vector_state, ONLY : a_uc, a_vc, a_wc, a_up, a_vp, a_wp
      USE grid, ONLY : dtl, dtlt,  &
                       dtlv, dtlong, nzp, nyp, nxp, level
      USE thrm, ONLY : thermo
      USE mo_output, ONLY : write_main, close_main, write_ps, close_ps, sflg, ps_intvl, main_intvl
      USE mo_history, ONLY : write_hist
      
      LOGICAL, PARAMETER :: StopOnCFLViolation = .FALSE.
      REAL, PARAMETER :: cfl_upper = 0.50, cfl_lower = 0.30

      REAL         :: t1,t2,tplsdt,begtime
      REAL(kind=8) :: cflmax,gcflmax
      INTEGER      :: istp
      LOGICAL :: cflflg
      !
      ! Timestep loop for program
      !
      begtime = time
      istp = 0

      CALL cpu_time(t1)

      DO WHILE (time + 0.1*dtl < timmax)

         istp = istp+1
         tplsdt = time + dtl + 0.1*dtl
         sflg = ( mod(tplsdt,ps_intvl) < dtl  &
            .OR. tplsdt >= timmax  .OR. tplsdt < 2.*dtl)

         CALL t_step(cflflg,cflmax)

         time = time + dtl

         CALL double_scalar_par_max(cflmax,gcflmax)
         cflmax = gcflmax

         IF (cflmax > cfl_upper .OR. cflmax < cfl_lower) THEN
            CALL tstep_reset(nzp,nxp,nyp,a_up%d,a_vp%d,a_wp%d,a_uc%d,a_vc%d,a_wc%d,     &
                             dtl,dtlong,cflmax,cfl_upper,cfl_lower)
            dtlv = 2.*dtl
            dtlt = dtl
         END IF

         !
         ! output control
         !
         IF (mod(tplsdt,ps_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_ps(time)

         IF ((mod(tplsdt,frqhis) < dtl .OR. time >= timmax) .AND. outflg)   &
            CALL write_hist(2, time)

         IF ((mod(tplsdt,main_intvl) < dtl .OR. time >= timmax) .AND. outflg) THEN
            CALL thermo(level)
            CALL write_main(time)
         END IF

         IF (cflflg) THEN
            cflflg = .FALSE.
            IF (StopOnCFLViolation) CALL write_hist(-1,time)
         END IF

         IF(myid == 0) THEN
            CALL cpu_time(t2) ! t2-t1 is the actual time from output
            IF (mod(istp,istpfl) == 0 ) THEN
               PRINT "('   Timestep # ',i6," //     &
                  "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
                  istp, time, t2-t1
               CALL cpu_time(t1)
            END IF
         END IF

      END DO

      CALL write_hist(1, time)
      CALL close_main()
      CALL close_ps()

   END SUBROUTINE stepper
   !
   !----------------------------------------------------------------------
   ! Subroutine tstep_reset: Called to adjust current velocity and reset
   ! timestep based on cfl limits
   !
   SUBROUTINE tstep_reset(n1,n2,n3,up,vp,wp,uc,vc,wc,dtl,dtmx,cfl,c1,c2)

      INTEGER, INTENT (in)      :: n1,n2,n3
      REAL, INTENT (in)         :: up(n1,n2,n3),vp(n1,n2,n3),wp(n1,n2,n3),dtmx,c1,c2
      REAL(kind=8), INTENT (IN) :: cfl
      REAL, INTENT (inout)      :: uc(n1,n2,n3),vc(n1,n2,n3),wc(n1,n2,n3),dtl

      INTEGER :: i,j,k
      REAL    :: cbar, dtl_old

      cbar = (c1+c2)*0.5
      dtl_old = dtl

      IF (cfl > c1) dtl = min(dtmx,dtl*cbar/c1)
      IF (cfl < c2) dtl = min(dtmx,dtl*cbar/c2)

      DO j = 1, n3
         DO i = 1, n2
            DO k = 1, n1
               uc(k,i,j) = up(k,i,j) + (uc(k,i,j)-up(k,i,j))*dtl/dtl_old
               vc(k,i,j) = vp(k,i,j) + (vc(k,i,j)-vp(k,i,j))*dtl/dtl_old
               wc(k,i,j) = wp(k,i,j) + (wc(k,i,j)-wp(k,i,j))*dtl/dtl_old
            END DO
         END DO
      END DO

   END SUBROUTINE tstep_reset

   ! 
   !----------------------------------------------------------------------
   ! Subroutine set_LES_runtime: Set the status of process switches e.g.
   ! if they have a defined spinup time etc.
   !
   SUBROUTINE set_LES_runtime(time)
     USE mcrp, ONLY : sed_aero,   &
                      sed_cloud,  &
                      sed_precp,  &
                      sed_ice,    &
                      bulk_autoc
     USE grid, ONLY : level
     IMPLICIT NONE

     REAL, INTENT(in) :: time

     IF ( sed_aero%switch .AND. time > sed_aero%delay ) sed_aero%state = .TRUE.
     IF ( sed_cloud%switch .AND. time > sed_cloud%delay ) sed_cloud%state = .TRUE.
     IF ( sed_precp%switch .AND. time > sed_precp%delay ) sed_precp%state = .TRUE.
     IF ( sed_ice%switch .AND. time > sed_ice%delay ) sed_ice%state = .TRUE.
     IF ( bulk_autoc%switch .AND. time > bulk_autoc%delay ) bulk_autoc%state = .TRUE.
     IF (level < 5) THEN
        sed_ice%state = .FALSE.
     END IF

   END SUBROUTINE set_LES_runtime

   !
   !----------------------------------------------------------------------
   ! Subroutine t_step: Called by driver to timestep through the LES
   ! routines.  Within many subroutines, data is accumulated during
   ! the course of a timestep for the purposes of statistical analysis.
   !
   SUBROUTINE t_step(cflflg,cflmax)

      USE grid, ONLY : level,dtlt,      &
                       nxp,nyp,nzp,   &
                       a_nactd,  a_vactd
      USE mo_vector_state, ONLY : a_wp
      USE mo_field_types, ONLY : Diag, Prog
      USE sgsm, ONLY : diffuse
      USE srfc, ONLY : surface
      USE thrm, ONLY : thermo
      USE mcrp, ONLY : micro
      USE prss, ONLY : poisson
      USE advf, ONLY : fadvect, newdroplet
      USE advl, ONLY : ladvect
      USE forc, ONLY : forcings
      USE util, ONLY : maskactiv !Juha: Included for SALSA

      USE mo_salsa_driver, ONLY : run_SALSA

      USE constrain_SALSA, ONLY : SALSA_diagnostics, tend_constrain2

      LOGICAL, INTENT (out)      :: cflflg
      REAL(KIND=8), INTENT (out) :: cflmax

      REAL    :: zwp(nzp,nxp,nyp)  !! FOR SINGLE-COLUMN RUNS

      INTEGER :: nspec

      
      CALL set_LES_runtime(time)

      zwp = 0.5  ! single column run vertical velocity

      cflflg = .FALSE.

      ! Reset ALL tendencies here.
      !----------------------------------------------------------------
      ! "Scalar" timestep
      CALL tend0(.FALSE.)

      ! Put the newly activated to zero
      IF (level >= 4) THEN
         a_vactd = 0.
         a_nactd = 0.
      END IF

      CALL surface()

      CALL diffuse

      CALL sponge(0)

      IF (level >= 1) CALL forcings(time,strtim)

      CALL update_sclrs
      CALL tend0(.TRUE.)
      CALL thermo(level)

      ! SALSA timestep
      ! -----------------------
      IF (level >= 4) THEN

         nspec = spec%getNSpec(type="wet") ! Aerosol components + water
            
         IF ( nxp == 5 .AND. nyp == 5 ) THEN
            ! 1D -runs
            CALL run_SALSA(Diag,Prog,nzp,nxp,nyp,nspec,   &
                           zwp,a_nactd,a_vactd,dtlt,      &
                           time,level,.FALSE.             )
         ELSE
            !! for 2D or 3D runs
            CALL run_SALSA(Diag,Prog,nzp,nxp,nyp,nspec,   &
                           a_wp%d,a_nactd,a_vactd,dtlt,     &
                           time,level,.FALSE.             )
             
         END IF !nxp==5 and nyp == 5

         CALL tend_constrain2()
         CALL update_sclrs
         CALL tend0(.TRUE.)
         CALL SALSA_diagnostics()
         CALL thermo(level)

      END IF ! level >= 4

         
      !-------------------------------------------
      ! "Deposition" timestep
      ! Dont perform sedimentation or level 3 autoconversion during spinup (internal switches implemented)
      CALL micro(level)
      IF (level >= 4) CALL tend_constrain2()
      CALL update_sclrs
      CALL tend0(.TRUE.)
      IF (level >= 4) CALL SALSA_diagnostics()
      CALL thermo(level)

      !-------------------------------------------
      ! "Advection" timestep
      CALL fadvect
      
      IF (level >= 4) CALL tend_constrain2()
      CALL update_sclrs
      CALL tend0(.TRUE.)
      IF (level >= 4) CALL SALSA_diagnostics()
      CALL thermo(level)
      
      CALL corlos

      CALL ladvect

      CALL buoyancy

      CALL sponge(1)

      CALL poisson

      CALL cfl (cflflg, cflmax)

      IF (level >= 4) CALL SALSA_diagnostics()
      CALL thermo(level)

   END SUBROUTINE t_step
   !
   !----------------------------------------------------------------------
   ! Subroutine tend0: sets all tendency arrays to zero
   !
   SUBROUTINE tend0(sclonly)

     USE grid, ONLY : nscl, a_st, newsclr
     USE mo_vector_state, ONLY : a_ut, a_vt, a_wt

      LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

      INTEGER :: n

      IF( .NOT. sclonly) THEN
         a_ut%d = 0.; a_vt%d = 0.; a_wt%d = 0.
      END IF
      DO n = 1, nscl
         CALL newsclr(n)
         a_st = 0.
      END DO

   END SUBROUTINE tend0

   !
   !----------------------------------------------------------------------
   ! Subroutine cfl: Driver for calling CFL computation subroutine
   !
   SUBROUTINE cfl(cflflg,cflmax)

      USE grid, ONLY : nxp,nyp,nzp,dxi,dyi,dtlt
      USE mo_vector_state, ONLY : a_up, a_vp, a_wp
      USE mo_aux_state, ONLY : dzt
      !USE stat, ONLY : fill_scalar

      LOGICAL, INTENT(out) :: cflflg
      REAL(KIND=8), INTENT (out)   :: cflmax
      REAL, PARAMETER :: cflnum = 0.95

      cflmax =  cfll(nzp,nxp,nyp,a_up%d,a_vp%d,a_wp%d,dxi,dyi,dzt,dtlt)

      cflflg = (cflmax > cflnum)
      IF (cflflg) PRINT *, 'Warning CFL Violation :', cflmax

   END SUBROUTINE cfl
   !
   !----------------------------------------------------------------------
   ! Subroutine cfll: Checks CFL criteria, brings down the model if the
   ! maximum thershold is exceeded
   !
   REAL(KIND=8) FUNCTION cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w
      REAL, INTENT (in)    :: dxi,dyi,dtlt
      TYPE(FloatArray1d), INTENT(in) :: dzt

      INTEGER :: i, j, k
      cfll = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               cfll = max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                      abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt%d(k))))
            END DO
         END DO
      END DO

   END FUNCTION cfll
   !
   !----------------------------------------------------------------------
   ! Subroutine update_sclrs:  Updates scalars by applying tendency and
   ! boundary conditions
   !
   SUBROUTINE update_sclrs

      USE grid, ONLY : a_sp, a_st, nscl, nxyzp, nxp, nyp, nzp, &
                       dtlt, newsclr, isgstyp
      USE mo_aux_state, ONLY : dzt
      USE mo_progn_state, ONLY : a_qp
      USE sgsm, ONLY : tkeinit
      USE util, ONLY : sclrset

      INTEGER :: n

      DO n = 1, nscl
         CALL newsclr(n)
         CALL update(nzp,nxp,nyp,a_sp,a_st,dtlt)
         CALL sclrset('mixd',nzp,nxp,nyp,a_sp,dzt%d)
      END DO

      IF (isgstyp == 2) THEN
         CALL tkeinit(nxyzp,a_qp%d)
      END IF

   END SUBROUTINE update_sclrs
   !
   ! ----------------------------------------------------------------------
   ! Subroutine update:
   !
   SUBROUTINE update(n1,n2,n3,a,fa,dt)

      INTEGER, INTENT(in)   :: n1, n2, n3
      REAL, INTENT (in)     :: fa(n1,n2,n3),dt
      REAL, INTENT (inout)  :: a(n1,n2,n3)
      INTEGER :: i, j, k

      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 2, n1-1
               a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
            END DO
         END DO
      END DO

   END SUBROUTINE update
   !
   ! ----------------------------------------------------------------------
   ! Subroutine buoyancy:
   !
   SUBROUTINE buoyancy
     
     USE grid, ONLY : nxp, nyp, nzp, th00, level
     USE mo_vector_state, ONLY : a_wt
     !USE stat, ONLY : sflg, comp_tke
     USE mo_diag_state, ONLY : a_rv, a_rc, a_theta, a_srp, a_ri, a_riri
     USE mo_progn_state, ONLY : a_rp
     USE mo_aux_state, ONLY : pi1
     USE util, ONLY : ae1mm
     USE thrm, ONLY : update_pi1
     
     REAL :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)
     
     IF (level < 4) THEN
        rv = a_rv%d ! Water vapor
        rc = a_rp%d - a_rv%d ! Total condensate (cloud + precipitation)
     ELSE IF (level == 4) THEN
        rv = a_rp%d ! Water vapor
        rc = a_rc%d + a_srp%d 
     ELSE IF (level == 5) THEN
        rv = a_rp%d
        rc = a_rc%d + a_srp%d + a_ri%d + a_riri%d  ! Total condensed water (aerosol+cloud+precipitation+ice)
     END IF
     call boyanc(nzp,nxp,nyp,a_wt%d,a_theta,rv,th00,a_tmp1,rc)
     
     CALL ae1mm(nzp,nxp,nyp,a_wt%d,awtbar)
     CALL update_pi1(nzp,awtbar,pi1)
     
   END SUBROUTINE buoyancy
   !
   ! ----------------------------------------------------------------------
   ! Subroutine boyanc:
   !
   SUBROUTINE boyanc(n1,n2,n3,wt,th,rv,th00,scr,rc)

     USE defs, ONLY : g, ep2

     INTEGER, INTENT(in) :: n1,n2,n3
     REAL, INTENT(in)    :: th00
     TYPE(FloatArray3d), INTENT(in) :: th

     REAL, INTENT(in)    :: rv(n1,n2,n3)  ! water vapor                                      
     REAL, INTENT(in)    :: rc(n1,n2,n3)  ! Total condensed water (aerosol, cloud, rain, ice and snow) mixing ratio

     REAL, INTENT(inout) :: wt(n1,n2,n3)
     REAL, INTENT(out)   :: scr(n1,n2,n3)

     INTEGER :: k, i, j
     REAL    :: gover2
     
     gover2 = 0.5*g

     do j=3,n3-2
        do i=3,n2-2
           do k=1,n1
              scr(k,i,j)=gover2*((th%d(k,i,j)*(1.+ep2*rv(k,i,j))-th00)/th00-rc(k,i,j))
           end do
           
           do k=2,n1-2
              wt(k,i,j)=wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
           end do
        end do
     end do
     
   END SUBROUTINE boyanc
   !
   ! ----------------------------------------------------------------------
   ! Subroutine corlos:  This is the coriolis driver, its purpose is to
   ! from the coriolis accelerations for u and v and add them into the
   ! accumulated tendency arrays of ut and vt.
   !
   SUBROUTINE corlos

     USE defs, ONLY : omega
     USE mo_aux_state, ONLY : u0, v0
      USE grid, ONLY : nxp, nyp, nzp, cntlat
      USE mo_vector_state, ONLY : a_uc, a_vc, a_ut, a_vt
      
      LOGICAL, SAVE :: initialized = .FALSE.
      REAL, SAVE    :: fcor

      INTEGER :: i, j, k

      IF (corflg) THEN
         IF (.NOT. initialized) fcor = 2.*omega*sin(cntlat*0.01745329)
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 2, nzp
                  a_ut%d(k,i,j) = a_ut%d(k,i,j) - fcor*(v0%d(k)-0.25*                   &
                                (a_vc%d(k,i,j)+a_vc%d(k,i+1,j)+a_vc%d(k,i,j-1)+a_vc%d(k,i+1,j-1)))
                  a_vt%d(k,i,j) = a_vt%d(k,i,j) + fcor*(u0%d(k)-0.25*                   &
                                (a_uc%d(k,i,j)+a_uc%d(k,i-1,j)+a_uc%d(k,i,j+1)+a_uc%d(k,i-1,j+1)))
               END DO
            END DO
         END DO
         initialized = .TRUE.
      END IF

   END SUBROUTINE corlos
   !
   ! ----------------------------------------------------------------------
   ! Subroutine sponge: does the rayleigh friction for the momentum terms,
   ! and newtonian damping of thermal term the damping is accumulated with the
   ! other tendencies
   !
   SUBROUTINE sponge (isponge)

      USE grid, ONLY : nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th00
      USE mo_vector_state, ONLY : a_up, a_vp, a_wp, a_ut, a_vt, a_wt
      USE mo_progn_state, ONLY : a_tp, a_tt
      USE mo_aux_state, ONLY : u0, v0, th0
      
      INTEGER, INTENT (in) :: isponge

      INTEGER :: i, j, k, kk

      IF (maxval(spng_tfct) > epsilon(1.) .AND. nfpt > 1) THEN
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = nzp-nfpt, nzp-1
                  kk = k+1-(nzp-nfpt)
                  IF (isponge == 0) THEN
                     a_tt%d(k,i,j) = a_tt%d(k,i,j) - spng_tfct(kk)*                   &
                                   (a_tp%d(k,i,j)-th0%d(k)+th00)
                  ELSE
                     a_ut%d(k,i,j) = a_ut%d(k,i,j) - spng_tfct(kk)*(a_up%d(k,i,j)-u0%d(k))
                     a_vt%d(k,i,j) = a_vt%d(k,i,j) - spng_tfct(kk)*(a_vp%d(k,i,j)-v0%d(k))
                     a_wt%d(k,i,j) = a_wt%d(k,i,j) - spng_wfct(kk)*(a_wp%d(k,i,j))
                  END IF
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sponge

END MODULE step
