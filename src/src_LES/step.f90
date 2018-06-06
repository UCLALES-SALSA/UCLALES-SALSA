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

  USE mo_submctl, ONLY : spec
  USE util, ONLY : getMassIndex, calc_correlation 
  
  IMPLICIT NONE
  
  INTEGER :: istpfl = 1
  REAL    :: timmax = 18000.
  LOGICAL :: corflg = .FALSE.
  
  REAL    :: frqhis =  9000.
  REAL    :: frqanl =  3600.
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

      USE grid, ONLY : dtl, dzt, zt, zm, nzp, dn0, u0, v0, a_up, a_vp, a_wp, &
                       a_uc, a_vc, a_wc, write_hist, write_anal, close_anal, dtlt,  &
                       dtlv, dtlong, nzp, nyp, nxp, level,                          &
                       ! For mass budged
                       a_rp, a_rc, a_srp, a_dn

      USE stat, ONLY : sflg, savg_intvl, ssam_intvl, write_ps, close_stat, mcflg, acc_massbudged,  &
                       write_massbudged
      USE thrm, ONLY : thermo

      LOGICAL, PARAMETER :: StopOnCFLViolation = .FALSE.
      REAL, PARAMETER :: cfl_upper = 0.50, cfl_lower = 0.30

      REAL         :: t1,t2,tplsdt,begtime
      REAL(kind=8) :: cflmax,gcflmax
      INTEGER      :: istp, iret
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
         sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl  &
            .OR. tplsdt >= timmax  .OR. tplsdt < 2.*dtl)

         CALL t_step(cflflg,cflmax)

         time = time + dtl

         CALL double_scalar_par_max(cflmax,gcflmax)
         cflmax = gcflmax

         IF (cflmax > cfl_upper .OR. cflmax < cfl_lower) THEN
            CALL tstep_reset(nzp,nxp,nyp,a_up,a_vp,a_wp,a_uc,a_vc,a_wc,     &
                             dtl,dtlong,cflmax,cfl_upper,cfl_lower)
            dtlv = 2.*dtl
            dtlt = dtl
         END IF

         !
         ! output control
         !
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_ps(nzp,dn0,u0,v0,zm,zt,time)

         IF ((mod(tplsdt,frqhis) < dtl .OR. time >= timmax) .AND. outflg)   &
            CALL write_hist(2, time)
         IF (mod(tplsdt,savg_intvl) < dtl .OR. time >= timmax .OR. time == dtl)   &
            CALL write_hist(1, time)

         IF ((mod(tplsdt,frqanl) < dtl .OR. time >= timmax) .AND. outflg) THEN
            CALL thermo(level)
            CALL write_anal(time)
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

      IF (mcflg) THEN
         !
         ! Juha:
         ! Get the final statistics of atmospheric water for mass budged
         CALL acc_massbudged(nzp,nxp,nyp,1,dtlt,dzt,a_dn,    &
                             rv=a_rp,rc=a_rc,prc=a_srp)

         CALL write_massbudged

      END IF ! mcflg

      CALL write_hist(1, time)
      iret = close_anal()
      iret = close_stat()

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
   ! Subroutine t_step: Called by driver to timestep through the LES
   ! routines.  Within many subroutines, data is accumulated during
   ! the course of a timestep for the purposes of statistical analysis.
   !
   SUBROUTINE t_step(cflflg,cflmax)

      USE grid, ONLY : level, dtl, dtlt, Tspinup,                                         &
                       ! Added parameters for interfacing with SALSA
                       nxp, nyp, nzp, a_press, a_temp, a_rsl,                             &
                       a_rc, a_wp, a_rp, a_rt, a_rh,                                      &
                       a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                       a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                       a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                       a_gaerop, a_gaerot, a_dn,  a_nactd,  a_vactd,            &
                       a_rsi

      USE stat, ONLY : sflg, statistics
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

      LOGICAL, INTENT (out)      :: cflflg
      REAL(KIND=8), INTENT (out) :: cflmax

      LOGICAL :: zactmask(nzp,nxp,nyp)
      REAL    :: zwp(nzp,nxp,nyp)  !! FOR SINGLE-COLUMN RUNS
      INTEGER :: zrm

      INTEGER :: n4

      zwp = 0.5

      cflflg = .FALSE.

      ! The runmode parameter zrm is used by SALSA only
      zrm = 3
      IF ( time < Tspinup ) zrm = 2

      ! Reset ALL tendencies here.
      !----------------------------------------------------------------
      ! "Scalar" timestep
      CALL tend0(.FALSE.)

      ! Put the newly activated to zero
      IF (level >= 4) THEN
         a_vactd = 0.
         a_nactd = 0.
      END IF

      IF (level >= 4 .AND. time < 1.) THEN
         CALL thermo(level)
         CALL SALSA_diagnostics
      END IF

      CALL surface()

      CALL diffuse

      CALL sponge(0)

      IF (level >= 1) THEN

         CALL thermo(level)

         CALL forcings(time,strtim)

         IF (level >= 4) THEN

            n4 = spec%getNSpec() ! Aerosol components + water

            CALL tend_constrain(n4)
            CALL update_sclrs
            CALL tend0(.TRUE.)

            IF ( nxp == 5 .AND. nyp == 5 ) THEN
               ! 1D -runs
               CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,zwp,a_dn,  &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              zrm, dtlt, time, level  )
            ELSE
               !! for 2D or 3D runs
               CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_wp,a_dn,  &
                              a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                              a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                              a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                              a_nicep,   a_nicet,   a_micep,   a_micet,    &
                              a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                              a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                              zrm, dtlt, time, level  )
             
            END IF !nxp==5 and nyp == 5

            CALL tend_constrain(n4)
         END IF

      END IF ! level

      CALL update_sclrs

      !-------------------------------------------
      ! "Deposition" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Dont perform sedimentation or level 3 autoconversion during spinup
      IF (zrm == 3) CALL micro(level)

      IF (level >= 4) CALL tend_constrain(n4)
      CALL update_sclrs

      !-------------------------------------------
      ! "Advection" timestep
      ! -- Reset only scalar tendencies
      CALL tend0(.TRUE.)

      ! Mask for cloud base activation
      IF (level >= 4) CALL maskactiv(zactmask,nxp,nyp,nzp,2,a_rh,rc=a_rc,w=a_wp)
      ! Get tendencies from cloud base activation
      IF (level >= 4) CALL newdroplet(zactmask)

      CALL fadvect

      IF (level >= 4)  &
         CALL tend_constrain(n4)

      CALL update_sclrs

      CALL thermo(level)

      IF (level >= 4)  THEN
         CALL SALSA_diagnostics
         CALL thermo(level)
      END IF

      CALL corlos

      CALL ladvect

      CALL buoyancy

      CALL sponge(1)

      CALL poisson

      CALL cfl (cflflg, cflmax)

      CALL thermo(level)

      IF (level >= 4)  THEN
         CALL SALSA_diagnostics
         call thermo(level)
      ENDIF

      IF (sflg) THEN
         CALL statistics (time+dtl)
      END IF

   END SUBROUTINE t_step
   !
   !----------------------------------------------------------------------
   ! Subroutine tend0: sets all tendency arrays to zero
   !
   SUBROUTINE tend0(sclonly)

      USE grid, ONLY : a_ut, a_vt, a_wt, nscl, a_st, newsclr

      LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

      INTEGER :: n

      IF( .NOT. sclonly) THEN
         a_ut = 0.; a_vt = 0.; a_wt = 0.
      END IF
      DO n = 1, nscl
         CALL newsclr(n)
         a_st = 0.
      END DO

   END SUBROUTINE tend0
   !
   !----------------------------------------------------------------------
   ! In case of negative tendencies to SALSA arrays, put some constrains
   ! in order to avoid concentrations going negative. This will possibly
   ! slightly affect the conservation of mass - needs testing/revision
   ! Juha Tonttila, FMI, 2014
   !
   SUBROUTINE tend_constrain(nn)

      USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                       a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                       a_nicep,  a_nicet, a_nsnowp, a_nsnowt,                            & ! ice'n'snow
                       a_micep,  a_micet, a_msnowp, a_msnowt,                            & ! ice'n'snow
                       dtlt, nxp,nyp,nzp, level
      USE mo_submctl, ONLY : nbins, ncld, nprc, &
                             nice, nsnw          !ice'n'snow

      INTEGER, INTENT(in) :: nn

      INTEGER :: cc, ii,jj,kk,ni

      DO jj = 3, nyp-2

         DO ii = 3, nxp-2

            DO kk = 1, nzp

               ! Aerosols
               DO cc = 1, nbins

                  IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_naerot(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop(kk,ii,jj,cc))/dtlt,a_naerot(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( ((1.e-10-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtlt,  &
                                                                 a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                     END DO

                  END IF

               END DO

               ! Cloud droplets
               DO cc = 1, ncld

                  IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_ncloudt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp(kk,ii,jj,cc))/dtlt,a_ncloudt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtlt,  &
                                                                 a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                     END DO

                  END IF

               END DO

               ! Precipitation
               DO cc = 1, nprc

                  IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nprecpt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp(kk,ii,jj,cc))/dtlt,a_nprecpt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtlt,  &
                                                                 a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                     END DO

                  END IF

               END DO

             ! ice particles
             IF (level<5) CYCLE
             DO cc = 1,nice

                  IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nicet(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep(kk,ii,jj,cc))/dtlt,a_nicet(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_micet(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtlt,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                     END DO

                  END IF

               END DO

               ! Snow
               DO cc = 1, nsnw

                  IF ( a_nsnowp(kk,ii,jj,cc)+a_nsnowt(kk,ii,jj,cc)*dtlt < 0. ) THEN

                     a_nsnowt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp(kk,ii,jj,cc))/dtlt,a_nsnowt(kk,ii,jj,cc))
                     DO ni = 1, nn
                        a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) = MAX( ((1.e-10-1.0)*a_msnowp(kk,ii,jj,(ni-1)*nsnw+cc))/dtlt,  &
                                                                a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) )
                     END DO

                  END IF

               END DO

            END DO ! kk

         END DO ! ii

      END DO ! jj

   END SUBROUTINE tend_constrain
   !
   !----------------------------------------------------------------------
   ! Subroutine cfl: Driver for calling CFL computation subroutine
   !
   SUBROUTINE cfl(cflflg,cflmax)

      USE grid, ONLY : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtlt
      USE stat, ONLY : fill_scalar

      LOGICAL, INTENT(out) :: cflflg
      REAL(KIND=8), INTENT (out)   :: cflmax
      REAL, PARAMETER :: cflnum = 0.95

      cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtlt)

      cflflg = (cflmax > cflnum)
      IF (cflflg) PRINT *, 'Warning CFL Violation :', cflmax
      CALL fill_scalar(1,REAL(cflmax))

   END SUBROUTINE cfl
   !
   !----------------------------------------------------------------------
   ! Subroutine cfll: Checks CFL criteria, brings down the model if the
   ! maximum thershold is exceeded
   !
   REAL(KIND=8) FUNCTION cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, DIMENSION (n1,n2,n3), INTENT (in) :: u, v, w
      REAL, INTENT (in)    :: dxi,dyi,dzt(n1),dtlt

      INTEGER :: i, j, k
      cfll = 0.
      DO j = 3, n3-2
         DO i = 3, n2-2
            DO k = 1, n1
               cfll = max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                      abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
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

      USE grid, ONLY : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
                       dtlt, newsclr, isgstyp
      USE sgsm, ONLY : tkeinit
      USE util, ONLY : sclrset

      INTEGER :: n

      DO n = 1, nscl
         CALL newsclr(n)
         CALL update(nzp,nxp,nyp,a_sp,a_st,dtlt)
         CALL sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
      END DO

      IF (isgstyp == 2) THEN
         CALL tkeinit(nxyzp,a_qp)
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
     
     USE grid, ONLY : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, &
          a_rp, a_srp, a_ri, a_srs, nxp, nyp, nzp, dzm, th00, level, pi1
     USE stat, ONLY : sflg, comp_tke
     USE util, ONLY : ae1mm
     USE thrm, ONLY : update_pi1
     
     REAL :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)
     
     IF (level < 4) THEN
        rv = a_rv ! Water vapor
        rc = a_rp - a_rv ! Total condensate (cloud + precipitation)
     ELSE IF (level >= 4) THEN
        rv = a_rp ! Water vapor
        rc = a_rc + a_srp + a_ri + a_srs ! Total condensed water (aerosol+cloud+precipitation+ice+snow)
     END IF
     call boyanc(nzp,nxp,nyp,a_wt,a_theta,rv,th00,a_tmp1,rc)
     
     CALL ae1mm(nzp,nxp,nyp,a_wt,awtbar)
     CALL update_pi1(nzp,awtbar,pi1)
     
     IF (sflg)  CALL comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_tmp1)
     
   END SUBROUTINE buoyancy
   !
   ! ----------------------------------------------------------------------
   ! Subroutine boyanc:
   !
   SUBROUTINE boyanc(n1,n2,n3,wt,th,rv,th00,scr,rc)

      USE defs, ONLY : g, ep2

      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in)    :: th00,th(n1,n2,n3),  &
                             rv(n1,n2,n3)  ! water vapor
                                      
      REAL, INTENT(in)    :: rc(n1,n2,n3)  ! Total condensed water (aerosol, cloud, rain, ice and snow) mixing ratio

      REAL, INTENT(inout) :: wt(n1,n2,n3)
      REAL, INTENT(out)   :: scr(n1,n2,n3)

      INTEGER :: k, i, j
      REAL    :: gover2

      gover2 = 0.5*g

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rv(k,i,j))-th00)/th00-rc(k,i,j))
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
      USE grid, ONLY : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0, cntlat

      LOGICAL, SAVE :: initialized = .FALSE.
      REAL, SAVE    :: fcor

      INTEGER :: i, j, k

      IF (corflg) THEN
         IF (.NOT. initialized) fcor = 2.*omega*sin(cntlat*0.01745329)
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 2, nzp
                  a_ut(k,i,j) = a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                                (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                  a_vt(k,i,j) = a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                                (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
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

      USE grid, ONLY : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
                       nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

      INTEGER, INTENT (in) :: isponge

      INTEGER :: i, j, k, kk

      IF (maxval(spng_tfct) > epsilon(1.) .AND. nfpt > 1) THEN
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = nzp-nfpt, nzp-1
                  kk = k+1-(nzp-nfpt)
                  IF (isponge == 0) THEN
                     a_tt(k,i,j) = a_tt(k,i,j) - spng_tfct(kk)*                   &
                                   (a_tp(k,i,j)-th0(k)+th00)
                  ELSE
                     a_ut(k,i,j) = a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                     a_vt(k,i,j) = a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                     a_wt(k,i,j) = a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                  END IF
               END DO
            END DO
         END DO
      END IF

   END SUBROUTINE sponge

   !
   ! ---------------------------------------------------------------------
   ! SALSA_diagnostics: Update properties for the current timestep:
   !                    E.g. if enough water has evaporated from droplets,
   !                    deplete the cloud droplet bins and move CCN material
   !                    back to the aerosol regime.
   !                    In addition, update the diagnostic scalars for total grid-cell
   !                    liquid water contents.
   !
   ! Juha Tonttila, FMI, 2014
   ! Tomi Raatikainen, FMI, 2016

   SUBROUTINE SALSA_diagnostics
      USE grid, ONLY : nxp,nyp,nzp,    &
                       a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,a_gaerop, &
                       a_rc, a_srp,a_snrp,  &
                       a_rh, a_temp, a_ri,a_srs,a_snrs,a_rhi,                                      &
                       a_nicep,a_micep,a_nsnowp,a_msnowp, level,   &
                       binMixrat
      USE mo_submctl, ONLY : nbins,ncld,nprc,ica,fca,icb,fcb,ira,fra,              &
                             in1a,in2a,fn2a,fn2b,                        &
                             nice,nsnw,iia,fia,iib,fib,isa,fsa,                    &
                             spec, surfw0, rg, nlim, prlim, pi, &
                             lscndgas, pi6, avog,                                  &
                             aerobins
      IMPLICIT NONE

      INTEGER :: i,j,k,bb,bc,ba,s,sc,sa,str,end,nc,nn

      REAL :: cd

      REAL :: ra, rb

    REAL :: zvol
    REAL :: zdh2o
    REAL :: ns, zbb, zaa ! Number of moles, Raoult effect, Kelvin effect; For calculating the critical radius
    REAL :: cdcld,cdprc ! Critical diameter for cloud droplets and precipitation
    REAL :: zdrms, zwams

    ! Remove negative values
    a_naerop = MAX(0.,a_naerop)
    a_ncloudp = MAX(0.,a_ncloudp)
    a_nprecpp = MAX(0.,a_nprecpp)
    a_maerop = MAX(0.,a_maerop)
    a_mcloudp = MAX(0.,a_mcloudp)
    a_mprecpp = MAX(0.,a_mprecpp)
    
    a_nicep = MAX(0.,a_nicep)
    a_nsnowp = MAX(0.,a_nsnowp)
    a_micep = MAX(0.,a_micep)
    a_msnowp = MAX(0.,a_msnowp)
    
    nn = spec%getNSpec() ! total number of species
    
    ! Remove particles that have number but not mass
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp
             ! Aerosols
             DO bc = 1,nbins    
                IF (a_naerop(k,i,j,bc) > 0. .AND. SUM(a_maerop(k,i,j,bc:getMassIndex(nbins,bc,nn):nbins)) <= 0.) THEN
                   a_naerop(k,i,j,bc) = 0.
                   a_maerop(k,i,j,bc:getMassIndex(nbins,bc,nn):nbins) = 0.
                END IF
             END DO

             ! Clouds
             DO bc = 1,ncld
                IF (a_ncloudp(k,i,j,bc) > 0. .AND. SUM(a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld)) <= 0.) THEN
                   a_ncloudp(k,i,j,bc) = 0.
                   a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld) = 0.
                END IF
             END DO ! ncld

             ! Precipitation
             DO bc = 1,nprc
                IF (a_nprecpp(k,i,j,bc) > 0. .AND. a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn)) <= 0.) THEN
                   a_nprecpp(k,i,j,bc) = 0.
                   a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn):nprc) = 0.
                END IF
             END DO ! nprc

             ! Ice
             IF (level<5) CYCLE
             DO bc = 1,nice
                IF (a_nicep(k,i,j,bc) > 0. .AND. SUM(a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice)) <= 0.) THEN
                   a_nicep(k,i,j,bc) = 0.
                   a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice) = 0.
                END IF
             END DO ! ncld

             ! Snow
             DO bc = 1,nsnw
                IF (a_nsnowp(k,i,j,bc) > 0. .AND. a_msnowp(k,i,j,getMassIndex(nsnw,bc,nn)) <= 0.) THEN
                   a_nsnowp(k,i,j,bc) = 0.
                   a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn):nsnw) = 0.
                END IF
             END DO ! nsnw

          END DO !k
       END DO !i
    END DO !j

   ! Ghost species
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp

             ! Loop over cloud droplet bins
             DO bc = ica%cur,fcb%cur

                IF ( a_ncloudp(k,i,j,bc) > nlim .AND. a_rh(k,i,j)<0.999 .AND.   &
                     a_mcloudp(k,i,j,getMassIndex(ncld,bc,nn)) < 1.e-5  ) THEN

                   ! Critical diameter
                   ns = SUM( spec%diss(1:nn-1)*a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn-1):ncld)/spec%MM(1:nn-1) ) / &
                        a_ncloudp(k,i,j,bc)
                   zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                   zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                   cdcld = SQRT(3.*zbb/zaa)

                   ! Wet diameter
                   zvol = SUM( a_mcloudp(k,i,j,bc:getMassIndex(ncld,bc,nn):ncld)/spec%rholiq(1:nn) )/a_ncloudp(k,i,j,bc)
                   zdh2o = (zvol/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.2*(critical size) or 2 um
                   IF ( zdh2o < MAX(0.2*cdcld,2.e-6) .OR.     &
                        a_mcloudp(k,i,j,getMassIndex(ncld,bc,nn)) < 1.e-25*a_ncloudp(k,i,j,bc) )  THEN

                      IF (bc<=fca%cur) THEN
                          ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                      ELSE
                          ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                      END IF
                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(ncld,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO

                   END IF ! critical diameter

                END IF  ! blim

             END DO ! bc

             ! Loop over precipitation bins
             DO bc = ira,fra

                IF ( a_nprecpp(k,i,j,bc) > prlim .AND. a_rh(k,i,j)<0.999 .AND.  &
                     a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn)) < 1.e-6 ) THEN

                   ! Critical radius
                   ns = SUM( spec%diss(1:nn-1)*a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc)/spec%MM(1:nn-1) ) / &
                        a_nprecpp(k,i,j,bc)
                   zbb = 6.*spec%mwa*ns/(pi*spec%rhowa)
                   zaa = 4.*spec%mwa*surfw0/(rg*spec%rhowa*a_temp(k,i,j))
                   cdprc = SQRT(3.*zbb/zaa)

                   ! Wet radius
                   zvol = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn):nprc)/spec%rholiq(1:nn) )/a_nprecpp(k,i,j,bc)
                   zdh2o = (zvol/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.02*critical radius or 2 um
                   IF ( zdh2o < MAX(0.02*cdprc,2.e-6) .OR.   &
                        a_mprecpp(k,i,j,getMassIndex(nprc,bc,nn))<1e-25*a_nprecpp(k,i,j,bc) ) THEN
                      ! Move evaporating rain drops to a soluble aerosol bin with
                      ! the closest match in dry particle mass. Ain't perfect but
                      ! the bin update subroutine in SALSA will take care of the rest.
                      zvol = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc) )/a_nprecpp(k,i,j,bc) ! Dry mass

                      ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
                      cd = SUM( a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc)/spec%rholiq(1:nn-1) ) / &
                           a_nprecpp(k,i,j,bc) ! Dry radius
                      cd = (cd/pi6)**(1./3.)
                      ba=in2a !Ignore 1a and note that "aerobins" contains the lower limit of bin dry diameter
                      DO WHILE (cd>=aerobins(ba+1) .AND. ba<fn2a)
                         ba=ba+1
                      END DO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb) <= nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSE IF (a_naerop(k,i,j,ba) <= nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra = calc_correlation(a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins),  &
                                               a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc),   &
                                               nn-1)

                         rb = calc_correlation(a_maerop(k,i,j,bb:getMassIndex(nbins,bb,nn-1):nbins),  &
                                               a_mprecpp(k,i,j,bc:getMassIndex(nprc,bc,nn-1):nprc),   &
                                               nn-1)
                         IF (ra<rb) ba = bb
                      END IF                     

                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                      a_nprecpp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nprc,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                         a_mprecpp(k,i,j,sc) = 0.
                      END DO

                   END IF ! Critical radius

                END IF ! prlim

             END DO ! bc

             ! Loop over ice bins
             DO bc = iia%cur,fib%cur

                IF ( a_nicep(k,i,j,bc) > prlim .AND. a_rhi(k,i,j)<0.999 .AND.   &
                     a_micep(k,i,j,getMassIndex(nice,bc,nn)) < 1.e-15 ) THEN

                   ! Ice and snow don't have a critical size, but lose particles when water content becomes low enough or size is < 2e-6m
                   
                   ! Diameter
                   cd = SUM( a_micep(k,i,j,bc:getMassIndex(nice,bc,nn):nice)/spec%rhoice(1:nn) ) / &
                        a_nicep(k,i,j,bc)
                   cd = (cd/pi6)**(1./3.)

                   ! Lose ice when dry to total mass ratio is more than 0.5
                   CALL binMixrat("ice","dry",bc,i,j,k,zdrms)
                   CALL binMixrat("ice","wet",bc,i,j,k,zwams)
                   zvol = zdrms/zwams

                   IF ( zvol>0.5 .OR. cd < 2.e-6 ) THEN
                      IF (bc<=fia%cur) THEN
                         ba = iia%par + (bc-iia%cur) ! Index for parallel aerosol bin
                      ELSE
                         ba = iib%par + (bc-iib%cur) ! Index for parallel aerosol bin
                      END IF

                      ! Move the number of particles from ice to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move mass material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nice,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF  ! prlim

             END DO ! bc

             ! Loop over snow bins
             DO bc = isa,fsa

                IF ( a_nsnowp(k,i,j,bc) > prlim .AND. a_rhi(k,i,j)<0.999 .AND.  &
                     a_msnowp(k,i,j,getMassIndex(nsnw,bc,nn)) < 1.e-20 ) THEN

                   cd = SUM( a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn):nsnw)/spec%rhosnow(1:nn) ) / &
                        a_nsnowp(k,i,j,bc)
                   cd = (cd/pi6)**(1./3.)

                   ! Lose snow when dry to total mass ratio is more than 0.5 or diameter < 2.e-6m
                   CALL binMixrat("snow","dry",bc,i,j,k,zdrms)
                   CALL binMixrat("snow","wet",bc,i,j,k,zwams)
                   zvol = zdrms/zwams

                   IF ( zvol>0.5 .OR. cd < 2.e-6) THEN
 
                      ! Move evaporating snow to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle diameter (a and b bins)
                      cd = SUM( a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw)/spec%rhosnow(1:nn-1) ) / &
                           a_nsnowp(k,i,j,bc)
                      cd = (cd/pi6)**(1./3.)
                      
                      ba=in2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd>=aerobins(ba+1) .AND. ba<fn2a)
                         ba=ba+1
                      END DO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb) <= nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSE IF (a_naerop(k,i,j,ba) <= nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra = calc_correlation(a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins),  &
                                               a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw),    &
                                               nn-1)
                         rb = calc_correlation(a_maerop(k,i,j,bb:getMassIndex(nbins,bb,nn-1):nbins),  &
                                               a_msnowp(k,i,j,bc:getMassIndex(nsnw,bc,nn-1):nsnw),    &
                                               nn-1)
                         IF (ra<rb) ba = bb
                      END IF


                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nsnowp(k,i,j,bc)
                      a_nsnowp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = getMassIndex(nsnw,bc,s)
                         sa = getMassIndex(nbins,ba,s)
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_msnowp(k,i,j,sc)
                         a_msnowp(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF ! prlim

             END DO ! bc

             ! Loop over aerosol bins
             DO ba = 1,nbins
                IF (a_naerop(k,i,j,ba) > nlim) THEN
                   zvol = SUM( a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn-1):nbins)/spec%rholiq(1:nn-1) )/ &
                          a_naerop(k,i,j,ba) ! Dry volume

                   ! Particles smaller then 0.1 nm diameter are set to zero 
                   IF ( zvol < pi6*1.e-10**3 ) THEN
                      ! Volatile species to the gas phase
                      IF (spec%isUsed('SO4') .AND. lscndgas) THEN
                         nc = spec%getIndex('SO4')
                         s = getMassIndex(nbins,ba,nc)
                         a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / spec%msu * avog
                      END IF
                      IF (spec%isUsed('OC') .AND. lscndgas) THEN
                         nc = spec%getIndex('OC')
                         s = getMassIndex(nbins,ba,nc)
                         a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / spec%moc * avog
                      END IF
                      IF (spec%isUsed('NO') .AND. lscndgas) THEN
                         nc = spec%getIndex('NO')
                         s = getMassIndex(nbins,ba,nc)
                         a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / spec%mno * avog
                      END IF
                      IF (spec%isUsed('NH') .AND. lscndgas) THEN
                         nc = spec%getIndex('NH')
                         s = getMassIndex(nbins,ba,nc)
                         a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / spec%mnh * avog
                      END IF

                      ! Mass and number to zero (insolube species and water are lost)
                      a_maerop(k,i,j,ba:getMassIndex(nbins,ba,nn):nbins) = 0.
                      a_naerop(k,i,j,ba) = 0.
                   END IF
                END IF
             END DO

          END DO   ! k
       END DO   ! i
    END DO   ! j

      !!!!!!!!!!!!!!!!!!!!!!!
      ! Update diagnostic tracers
      !!!!!!!!!!!!!!!!!!!!!!!

      ! Liquid water content
      nc = spec%getIndex('H2O')
      ! Aerosols, regimes a and b
      str = getMassIndex(nbins,in1a,nc)
      end = getMassIndex(nbins,fn2b,nc)
      a_rc(:,:,:) = SUM(a_maerop(:,:,:,str:end),DIM=4)
      ! Clouds, regime a and b
      str = getMassIndex(ncld,ica%cur,nc)
      end = getMassIndex(ncld,fcb%cur,nc)
      a_rc(:,:,:) = a_rc(:,:,:) + SUM(a_mcloudp(:,:,:,str:end),DIM=4)
      ! Precipitation
      str = getMassIndex(nprc,ira,nc)
      end = getMassIndex(nprc,fra,nc)
      a_srp(:,:,:) = SUM(a_mprecpp(:,:,:,str:end),DIM=4)
      a_snrp(:,:,:) = SUM(a_nprecpp(:,:,:,ira:fra),DIM=4)

      ! ice, regimes a and b
      str = getMassIndex(nice,iia%cur,nc)
      end = getMassIndex(nice,fib%cur,nc)
      a_ri(:,:,:) = SUM(a_micep(:,:,:,str:end),DIM=4)
      ! Snow
      str = getMassIndex(nsnw,isa,nc)
      end = getMassIndex(nsnw,fsa,nc)
      a_srs(:,:,:) = SUM(a_msnowp(:,:,:,str:end),DIM=4)
      a_snrs(:,:,:) = SUM(a_nsnowp(:,:,:,isa:fsa),DIM=4)

   END SUBROUTINE SALSA_diagnostics


END MODULE step
