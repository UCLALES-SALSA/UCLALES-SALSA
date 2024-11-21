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
module step

  implicit none

  integer :: istpfl = 1
  real    :: timmax = 18000.
  logical :: corflg = .false.

  real    :: frqhis =  9000.
  real    :: frqrst =  3600.
  real    :: frqanl =  3600.
  real    :: anl_start = -1.

  real    :: time   =  0.
  real    :: strtim =  0.0
  real    :: cntlat =  31.5 ! 30.0
  logical :: outflg = .true.


contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine model:  This is the main driver for the model's time
  ! integration.  It calls the routine tstep, which steps through the
  ! physical processes active on a time-step and updates variables.  It
  ! then checks to see whether or not different output options are
  ! satisfied.
  subroutine stepper

    use mpi_interface, only : myid, double_scalar_par_max

    use grid, only : dtl, zt, zm, nzp, dn0, u0, v0, &
         write_hist, write_anal, close_anal, &
         dtlong, nzp
    use stat, only : sflg, csflg, cswrite, savg_intvl, ssam_intvl, write_ps, close_stat

    real, parameter :: cfl_upper = 0.5

    real    :: t1,t2,tplsdt
    REAL(kind=8) :: cflmax,gcflmax
    integer :: istp, iret
    !
    ! Timestep loop for program
    !
    istp = 0

    call cpu_time(t1)

    do while (time < timmax)
       ! Limit time step based on the Courant-Friedrichs-Lewy condition
       call cfl(cflmax)
       call double_scalar_par_max(cflmax,gcflmax)
       cflmax = gcflmax
       dtl = min(dtlong,dtl*cfl_upper/(cflmax+epsilon(1.)))

       ! Determine when to compute statistics
       !    - When a given output or profile sampling time (n*tstep) will be reached or exceeded for the first time
       !    - After the first call (time=0)
       tplsdt = time + dtl ! Time after t_step
       sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl .or. tplsdt < 1.1*dtl)
       cswrite = csflg .and. (mod(tplsdt,savg_intvl) < dtl .or. tplsdt < 1.1*dtl)

       call t_step
       time = time + dtl

       ! Write profiles
       if (mod(tplsdt,savg_intvl) < dtl .or. tplsdt < 1.1*dtl)   &
            call write_ps(nzp,dn0,u0,v0,zm,zt,time)

       ! Write restarts (*.<time>s and *.rst)
       if (mod(tplsdt,frqhis) < dtl .and. outflg)   &
            call write_hist(2, time)

       if (mod(tplsdt,frqrst) < dtl .and. outflg)   &
            call write_hist(1, time)

       ! Write analysis files
       if (mod(tplsdt,frqanl) < dtl .and. outflg .and. time >= anl_start)   &
            call write_anal(time)

       if(myid == 0) then
          istp = istp+1
          if (mod(istp,istpfl) == 0 ) THEN
              call cpu_time(t2) ! t2-t1 is the actual CPU time from the previous output
              print "('   Timestep # ',i6," //     &
                 "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
                 istp, time, t2-t1
              call cpu_time(t1)
          ENDIF
       endif

    enddo

    call write_hist(1, time)
    iret = close_anal()
    iret = close_stat()

  end subroutine stepper
  !
  !----------------------------------------------------------------------
  ! subroutine t_step: Called by driver to timestep through the LES
  ! routines.  Within many subroutines, data is accumulated during
  ! the course of a timestep for the purposes of statistical analysis.
  !
  subroutine t_step()

    use grid, only : level, dtl, Tspinup, sst, nxp, nyp, nzp, a_press,    &
                     a_temp,  a_rp, a_rt, a_rsl, a_rsi, a_dn, a_edr,      &
                     a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                     a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                     a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                     a_gaerop, a_gaerot, nspt, nbins, ncld, nprc, nice, nsnw,           &
                     nudge_theta, nudge_rv, nudge_u, nudge_v, nudge_ccn,                &
                     ifSeaSpray, ifSeaVOC, sea_tspinup
    use stat, only : sflg, statistics, les_rate_stats, out_mcrp_data, out_mcrp_list, &
                     out_mcrp_nout, mcrp_var_save
    use sgsm, only : diffuse
    use srfc, only : surface, marine_aero_flux, marine_gas_flux, srfc_gas_flux
    use thrm, only : thermo
    use mcrp, only : micro
    use prss, only : poisson
    use advf, only : fadvect
    use advl, only : ladvect
    use forc, only : forcings

    USE mo_salsa_driver, ONLY : run_SALSA
    USE mo_submctl, ONLY : nvbs_setup

    real :: xtime
    LOGICAL :: zrm

    xtime = time/86400. + strtim

    ! Reset ALL tendencies here.
    !----------------------------------------------------------------
    ! "Scalar" timestep
    CALL tend0(.FALSE.)

    call surface(sst)
    IF (level > 3 .AND. time >= sea_tspinup) THEN
        if(ifSeaSpray)then
          CALL marine_aero_flux(MAX(271.15,sst)) ! T>-2 C
        endif
        if(nvbs_setup>=0 .and. ifSeaVOC) then
            call marine_gas_flux(sst)
        elseif (nvbs_setup>=0) THEN
            CALL srfc_gas_flux()
        endif
    END IF

    IF (sflg) CALL les_rate_stats('srfc')
    call update_sclrs(.TRUE.) ! Update edges for diffusion
    CALL tend0(.TRUE.)

    call diffuse

    IF (sflg) CALL les_rate_stats('diff')
    call update_sclrs
    CALL tend0(.TRUE.)

    call sponge(0)

    call thermo(level)

    call forcings(xtime,cntlat,sst)

    IF (sflg) CALL les_rate_stats('forc')

    IF (level >= 4) THEN
        call update_sclrs
        CALL tend0(.TRUE.)

        ! The runmode parameter zrm is used by SALSA only
        zrm = time < Tspinup

        CALL run_SALSA(nxp,nyp,nzp,nspt,nbins,ncld,nprc,nice,nsnw, &
                  a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_dn,a_edr, &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_gaerop,  a_gaerot,  zrm, dtl, time, level, &
                  sflg, out_mcrp_nout, out_mcrp_list, out_mcrp_data)

        IF (sflg) CALL les_rate_stats('mcrp')
        call update_sclrs
        CALL tend0(.TRUE.)

        ! Some diagnostic tests
        CALL SALSA_diagnostics

        ! Save user-selected details about SALSA microphysics
        IF (sflg) CALL mcrp_var_save()

        ! Sedimentation (not during spinup)
        IF (time >= Tspinup) CALL sedim_SALSA
        IF (sflg) CALL les_rate_stats('sedi')
    ELSEIF (time >= Tspinup) THEN
        ! Don't perform level 3 microphysics during spinup
        call update_sclrs
        CALL tend0(.TRUE.)

        CALL micro(level)

        IF (sflg) CALL les_rate_stats('mcrp')

        ! Save user-selected details about Seifert and Beheng microphysics
        IF (sflg) CALL mcrp_var_save()

    END IF

    CALL update_sclrs (.TRUE.) ! Update edges for advection
    call tend0(.TRUE.)

    !-------------------------------------------
    ! "Advection" timestep
    CALL fadvect

    IF (sflg) CALL les_rate_stats('advf')

    CALL update_sclrs

    CALL thermo(level)

    ! Nudging
     IF (nudge_theta/=0 .OR. nudge_rv/=0 .OR. nudge_u/=0 .OR. &
            nudge_v/=0 .OR. (level>3 .AND. nudge_ccn/=0) ) THEN

        ! Reset tendencies
        call tend0(.TRUE.)

        CALL nudging(time)

        IF (sflg) CALL les_rate_stats('nudg')
        CALL update_sclrs

        CALL thermo(level)
    ENDIF

    call corlos

    call ladvect

    call buoyancy

    call sponge(1)

    call poisson

    CALL thermo(level)

    if (sflg) call statistics (time+dtl)

  end subroutine t_step
  !
  !----------------------------------------------------------------------
  !
  ! Nudging towards the initial state (temperature, water vapor,
  ! horizontal winds and aerosol and/or cloud droplets).
  !
  ! TR 22.3.2017
  !
  SUBROUTINE nudging(time)

    use grid, only : level, dtl, nxp, nyp, nzp, nbins, ncld, nice, &
                zt, a_rp, a_rt, a_rc, a_ri, &
                a_naerop, a_naerot, a_ncloudp, a_nicep, &
                nspt, a_dn, a_maerop, a_maerot, &
                a_tp, a_tt, a_up, a_ut, a_vp, a_vt, &
                !th0, th00, rt0, u0, v0, &
                nudge_theta, nudge_theta_time, nudge_theta_zmin, nudge_theta_zmax, nudge_theta_tau, &
                nudge_rv, nudge_rv_time, nudge_rv_zmin, nudge_rv_zmax, nudge_rv_tau, &
                nudge_u, nudge_u_time, nudge_u_zmin, nudge_u_zmax, nudge_u_tau, &
                nudge_v, nudge_v_time, nudge_v_zmin, nudge_v_zmax, nudge_v_tau, &
                nudge_ccn, nudge_ccn_time, nudge_ccn_zmin, nudge_ccn_zmax, nudge_ccn_tau, &
                theta_ref, rv_ref, u_ref, v_ref, aero_ref, nudge_init
    USE mo_submctl, ONLY : in2a, nlim

    IMPLICIT NONE
    REAL, INTENT(IN) :: time
    REAL :: aero_target(nzp,nxp,nyp,nbins)

    ! Initialization
    IF (nudge_init) THEN
        ! Note: the first temperature and humidity values can include random
        ! perturbations, so could take the target values from soundings (th0, rt0).
        ! There are no wind perturbations, but can still could use u0 and v0.
        !
        ! (Liquid water) potential temperature: nudge towards initial theta
        IF (nudge_theta/=0) THEN
            ALLOCATE(theta_ref(nzp))
            theta_ref(:)=a_tp(:,3,3)
            !theta_ref(:)=th0(:)-th00 ! Initial state from soundings
        ENDIF
        !
        ! Water vapor mixing ratio based on total water
        !   Levels 0-3: total = total water (a_rp)
        !   Levels 4-5: total = water vapor (a_rp) + aerosol, cloud and rain water (a_rc) + ice and snow (a_ri)
        IF (nudge_rv/=0)  THEN
            ALLOCATE(rv_ref(nzp))
            IF (level>3) THEN
                rv_ref(:)=a_rp(:,3,3)+a_rc(:,3,3)+a_ri(:,3,3)
            ELSE ! Levels 0-3
                rv_ref(:)=a_rp(:,3,3) ! This includes all
            ENDIF
            !rv_ref(:)=rt0(:) ! Initial state from soundings
        ENDIF
        !
        ! Horizontal winds
        IF (nudge_u/=0) THEN
            ALLOCATE(u_ref(nzp))
            u_ref(:)=a_up(:,3,3)
            !u_ref(:)=u0(:) ! Initial state from soundings
        ENDIF
        IF (nudge_v/=0) THEN
            ALLOCATE(v_ref(nzp))
            v_ref(:)=a_vp(:,3,3)
            !v_ref(:)=v0(:) ! Initial state from soundings
        ENDIF
        !
        ! Nudge level 4 and 5 aerosol concentration based on total CCN = aerosol + cloud droplets + ice.
        ! Precipitation and snow are not included, because these cannot be related to a specific aerosol bin.
        IF (level>3 .AND. nudge_ccn/=0) THEN
            ! Nudge aerosol based on the total number (aerosol+cloud+ice)
            ALLOCATE(aero_ref(nzp,nbins))
            aero_ref(:,:)=a_naerop(:,3,3,:)
            aero_ref(:,in2a:nbins)=aero_ref(:,in2a:nbins)+a_ncloudp(:,3,3,1:ncld)
            IF (level==5) aero_ref(:,in2a:nbins)=aero_ref(:,in2a:nbins)+a_nicep(:,3,3,1:nice)
        ENDIF
        !
        ! Initialized
        nudge_init=.FALSE.
    ENDIF

    ! (Liquid water) potential temperature:
    IF (nudge_theta>0) &
        CALL nudge_any(nxp,nyp,nzp,zt,a_tp,a_tt,theta_ref,dtl,time,nudge_theta, &
            nudge_theta_time,nudge_theta_zmin,nudge_theta_zmax,nudge_theta_tau)

    ! Water vapor
    IF (nudge_rv>0) THEN
        IF (level>3) THEN
            ! Nudge water vapor (a_rp) based on total (vapor + cloud + rain [+ ice + snow])
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rc+a_ri,a_rt,rv_ref,dtl,time,nudge_rv, &
                nudge_rv_time,nudge_rv_zmin,nudge_rv_zmax,nudge_rv_tau)
        ELSE
            ! Nudge total water (a_rp) based on total water
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp,a_rt,rv_ref,dtl,time,nudge_rv, &
                nudge_rv_time,nudge_rv_zmin,nudge_rv_zmax,nudge_rv_tau)
        ENDIF
    ENDIF

    ! Horizontal winds
    IF (nudge_u>0) &
         CALL nudge_any(nxp,nyp,nzp,zt,a_up,a_ut,u_ref,dtl,time,nudge_u, &
            nudge_u_time,nudge_u_zmin,nudge_u_zmax,nudge_u_tau)
    IF (nudge_v>0) &
        CALL nudge_any(nxp,nyp,nzp,zt,a_vp,a_vt,v_ref,dtl,time,nudge_v, &
            nudge_v_time,nudge_v_zmin,nudge_v_zmax,nudge_v_tau)

    ! Aerosol
    IF (level>3 .AND. nudge_ccn/=0) THEN
        ! Target aerosol concentration = aerosol(t)+cloud(t)+ice(t)
        aero_target(:,:,:,:)=a_naerop(:,:,:,:)
        aero_target(:,:,:,in2a:nbins)=aero_target(:,:,:,in2a:nbins)+a_ncloudp(:,:,:,1:ncld)
        IF (level==5) aero_target(:,:,:,in2a:nbins)=aero_target(:,:,:,in2a:nbins)+a_nicep(:,:,:,1:nice)
        ! Apply to sectional data
        CALL nudge_any_2d(nxp,nyp,nzp,nbins,zt,aero_target,a_naerot,aero_ref,dtl,time,nudge_ccn, &
            nudge_ccn_time,nudge_ccn_zmin,nudge_ccn_zmax,nudge_ccn_tau)
        ! Change aerosol mass so that the mean dry and wet sizes are unchanged
        CALL adj_salsa_maerot(nxp,nyp,nzp,nbins,nspt-1,nlim,a_dn,a_naerop,a_naerot,a_maerop,a_maerot)
    ENDIF

  END SUBROUTINE nudging
  !
  ! Nudging for any 3D field based on 1D target
  SUBROUTINE nudge_any(nxp,nyp,nzp,zt,ap,at,trgt,dt,time,iopt,tref,zmin,zmax,tau)
    USE util, ONLY : get_avg3
    IMPLICIT NONE
    INTEGER :: nxp,nyp,nzp
    REAL :: zt(nzp), ap(nzp,nxp,nyp), at(nzp,nxp,nyp)
    REAL :: dt, time
    REAL :: trgt(nzp)
    REAL :: tref,zmin,zmax,tau
    INTEGER :: iopt
    INTEGER :: kk
    REAL :: avg(nzp)
    !
    IF (iopt==1 .AND. time<=tref) THEN
        ! Soft nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        CALL get_avg3(nzp,nxp,nyp,ap,avg)
        DO kk = 1,nzp
            IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                at(kk,:,:)=at(kk,:,:)-(avg(kk)-trgt(kk))/max(tau,dt)
        ENDDO
    ELSEIF (iopt==2 .AND. time<=tref) THEN
        ! Hard nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO kk = 1,nzp
            IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                at(kk,:,:)=at(kk,:,:)-(ap(kk,:,:)-trgt(kk))/max(tau,dt)
        ENDDO
    ENDIF
    !
  END SUBROUTINE nudge_any
  !
  ! Nudging for any 4D field based on 2D target
  SUBROUTINE nudge_any_2d(nxp,nyp,nzp,nb,zt,ap,at,trgt,dt,time,iopt,tref,zmin,zmax,tau)
    USE util, ONLY : get_avg3
    IMPLICIT NONE
    INTEGER :: nxp,nyp,nzp,nb
    REAL :: zt(nzp), ap(nzp,nxp,nyp,nb), at(nzp,nxp,nyp,nb)
    REAL :: dt, time
    REAL :: trgt(nzp,nb)
    REAL :: tref,zmin,zmax,tau
    INTEGER :: iopt
    INTEGER :: ii, kk
    REAL :: avg(nzp)
    !
    IF (iopt==1 .AND. time<=tref) THEN
        ! Soft nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO ii=1,nb
            CALL get_avg3(nzp,nxp,nyp,ap(:,:,:,ii),avg)
            DO kk = 1,nzp
                IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                    at(kk,:,:,ii)=at(kk,:,:,ii)-(avg(kk)-trgt(kk,ii))/max(tau,dt)
            ENDDO
        ENDDO
    ELSEIF (iopt==2 .AND. time<=tref) THEN
        ! Hard nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO ii=1,nb
            DO kk = 1,nzp
                IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                    at(kk,:,:,ii)=at(kk,:,:,ii)-(ap(kk,:,:,ii)-trgt(kk,ii))/max(tau,dt)
            ENDDO
        ENDDO
    ENDIF
    !
  END SUBROUTINE nudge_any_2d
  !
  ! Calculate SALSA aerosol mass tendency so that number nudging doesn't change the mean size or composition
  SUBROUTINE adj_salsa_maerot(nxp,nyp,nzp,nb,ns,nlim,adn,anp,ant,amp,amt)
    IMPLICIT NONE
    INTEGER :: nxp,nyp,nzp,nb,ns ! Dimensions
    REAL, INTENT(IN) :: nlim ! Aerosol number concentration limit (#/m3)
    REAL, INTENT(IN) :: adn(nzp,nxp,nyp) ! Air density (kg/m3)
    REAL, INTENT(IN) :: anp(nzp,nxp,nyp,nb), ant(nzp,nxp,nyp,nb) ! Aerosol number and number tendency
    REAL, INTENT(IN) :: amp(nzp,nxp,nyp,nb*(ns+1)) ! Aerosol mass
    REAL, INTENT(INOUT) :: amt(nzp,nxp,nyp,nb*(ns+1)) ! Aerosol mass tendency
    INTEGER :: i, j, k, bc
    !
    DO j = 3,nyp-2
      DO i = 3,nxp-2
        DO k = 1,nzp
          DO bc = 1,nb
            IF (anp(k,i,j,bc)*adn(k,i,j) > nlim) THEN
              amt(k,i,j,bc:ns*nb+bc:nb) = ant(k,i,j,bc)/anp(k,i,j,bc)*amp(k,i,j,bc:ns*nb+bc:nb)
            END IF
          END DO
        END DO
      END DO
    END DO
    !
  END SUBROUTINE adj_salsa_maerot
  !
  !----------------------------------------------------------------------
  ! subroutine tend0: sets all tendency arrays to zero
  !
  subroutine tend0(sclonly)

    use grid, only : a_ut, a_vt, a_wt, a_sclrt

    LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

    IF( .NOT. sclonly) THEN
       a_ut=0.; a_vt=0.; a_wt=0.
    ENDIF
    a_sclrt=0.

  end subroutine tend0
  !
  !----------------------------------------------------------------------
  ! Subroutine cfl: Driver for calling CFL computation subroutine
  !
  subroutine cfl(cflmax)

    use grid, only : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtl
    use stat, only : fill_scalar

    real(KIND=8), intent (out)   :: cflmax
    real, parameter :: cflnum=0.95

    cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtl)

    if (cflmax > cflnum) print *, 'Warning CFL Violation :', cflmax
    call fill_scalar(REAL(cflmax),'cfl    ','max')

  end subroutine cfl
  !
  !----------------------------------------------------------------------
  ! Subroutine cfll: Gets the peak CFL number
  !
  real(KIND=8) function cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (in)    :: dxi,dyi,dzt(n1),dtlt

    integer :: i, j, k
    cfll=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             cfll=max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                  abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
          end do
       end do
    end do

  end function cfll
  !
  !----------------------------------------------------------------------
  ! subroutine update_sclrs:  Updates scalars by applying tendency and
  ! boundary conditions
  !
  subroutine update_sclrs(doedges)

    use grid, only : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
         dtl, newsclr, isgstyp
    use sgsm, only : tkeinit
    use util, only : sclrset

    logical, optional, intent(in) :: doedges

    integer :: n

    do n=1,nscl
       call newsclr(n)
       call update(nzp,nxp,nyp,a_sp,a_st,dtl)
       call sclrset('cnst',nzp,nxp,nyp,a_sp,dzt,doedges)
    end do

    if (isgstyp == 2) then
       call tkeinit(nxyzp,a_qp)
    end if

  end subroutine update_sclrs
  !
  ! ----------------------------------------------------------------------
  ! subroutine update:
  !
  subroutine update(n1,n2,n3,a,fa,dt)

    integer, intent(in)   :: n1, n2, n3
    real, intent (in)     :: fa(n1,n2,n3),dt
    real, intent (in out) :: a(n1,n2,n3)
    integer :: i, j, k

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
          end do
       end do
    end do

  end subroutine update
  !
  ! ----------------------------------------------------------------------
  ! subroutine buoyancy:
  !
  subroutine buoyancy

    use grid, only : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, &
         a_rp, a_rpp, a_ri, nxp, nyp, nzp, dzm, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1

    real :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)

    IF (level<4) THEN
       rv = a_rv ! Water vapor
       rc = a_rc + a_rpp + a_ri ! Total condensate (cloud + precipitation + total ice)
    ELSE
       rv = a_rp ! Water vapor
       rc = a_rc + a_ri ! Total condensed water (aerosol+cloud+precipitation+ice+snow)
    END IF
    call boyanc(nzp,nxp,nyp,a_wt,a_theta,rv,th00,a_tmp1,rc)

    call ae1mm(nzp,nxp,nyp,a_wt,awtbar)
    call update_pi1(nzp,awtbar,pi1)

    if (sflg)  call comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_tmp1)

  end subroutine buoyancy
  !
  ! ----------------------------------------------------------------------
  ! subroutine boyanc:
  !
  subroutine boyanc(n1,n2,n3,wt,th,rv,th00,scr,rc)

    use defs, only: g, ep2

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: th00,th(n1,n2,n3),  &
                           rv(n1,n2,n3), &  ! Total water vapour mixing ratio
                           rc(n1,n2,n3)     ! Total condensed water (aerosol, cloud, rain, ice and snow) mixing ratio
    real, intent(inout) :: wt(n1,n2,n3)
    real, intent(out)   :: scr(n1,n2,n3)

    integer :: k, i, j
    real :: gover2

    gover2  = 0.5*g

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

  end subroutine boyanc
  !
  ! ----------------------------------------------------------------------
  ! subroutine corlos:  This is the coriolis driver, its purpose is to
  ! from the coriolis accelerations for u and v and add them into the
  ! accumulated tendency arrays of ut and vt.
  !
  subroutine corlos

    use defs, only : omega
    use grid, only : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0

    logical, save :: initialized = .False.
    real, save    :: fcor

    integer :: i, j, k

    if (corflg) then
       if (.not.initialized) fcor=2.*omega*sin(cntlat*0.01745329)
       do j=3,nyp-2
          do i=3,nxp-2
             do k=2,nzp
                a_ut(k,i,j)=a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                     (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                a_vt(k,i,j)=a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                     (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
             end do
          end do
       end do
       initialized = .True.
    end if

  end subroutine corlos
!
! ----------------------------------------------------------------------
! subroutine sponge: does the rayleigh friction for the momentum terms,
! and newtonian damping of thermal term the damping is accumulated with the
! other tendencies
!
  subroutine sponge (isponge)

    use grid, only : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
         nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00, spongeinit
    use util, only : get_pustat_vector

    integer, intent (in) :: isponge

    integer :: i, j, k, kk
    real :: tbar(nfpt),ubar(nfpt),vbar(nfpt), fact

    if (maxval(spng_tfct) > epsilon(1.) .and. nfpt > 1) then
     if(spongeinit) then ! Nudge sponge layer back to initial profile (default)
       do j=3,nyp-2
          do i=3,nxp-2
             do k=nzp-nfpt,nzp-1
                kk = k+1-(nzp-nfpt)
                if (isponge == 0) then
                   a_tt(k,i,j)=a_tt(k,i,j) - spng_tfct(kk)*                   &
                        (a_tp(k,i,j)-th0(k)+th00)
                else
                   a_ut(k,i,j)=a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                   a_vt(k,i,j)=a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                   a_wt(k,i,j)=a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                end if
             end do
          end do
       end do
     else                 ! Nudge sponge layer to bulk value
       fact = 1./float((nxp-4)*(nyp-4))
       do k = nzp-nfpt,nzp-1
          kk = k+1-(nzp-nfpt)
          tbar(kk) = sum(a_tp(k,3:nxp-2,3:nyp-2))*fact
          ubar(kk) = sum(a_up(k,3:nxp-2,3:nyp-2))*fact
          vbar(kk) = sum(a_vp(k,3:nxp-2,3:nyp-2))*fact
       end do

       CALL get_pustat_vector('avg', nfpt, tbar)
       CALL get_pustat_vector('avg', nfpt, ubar)
       CALL get_pustat_vector('avg', nfpt, vbar)

       do j=3,nyp-2
          do i=3,nxp-2
             do k=nzp-nfpt,nzp-1
                kk = k+1-(nzp-nfpt)
                if (isponge == 0) then
                   a_tt(k,i,j)=a_tt(k,i,j) - spng_tfct(kk)*(a_tp(k,i,j)-tbar(kk))
                else
                   a_ut(k,i,j)=a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-ubar(kk))
                   a_vt(k,i,j)=a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-vbar(kk))
                   a_wt(k,i,j)=a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                end if
             end do
          end do
       end do
     end if
    end if

  end subroutine sponge

  !
  ! ---------------------------------------------------------------------
  ! SALSA_diagnostics
  !
  SUBROUTINE SALSA_diagnostics()
    USE grid, ONLY : tmp_prcp, tmp_icep, tmp_snwp, tmp_gasp
    USE mo_submctl, ONLY : prlim

    ! Check that ignored species or bins are not used
    IF (ALLOCATED(tmp_prcp)) THEN
        IF (ANY(tmp_prcp>prlim)) STOP 'Non-zero rain even when disabled!'
    ENDIF
    IF (ALLOCATED(tmp_icep)) THEN
        IF (ANY(tmp_icep>prlim)) STOP 'Non-zero ice even when disabled!'
    ENDIF
    IF (ALLOCATED(tmp_snwp)) THEN
        IF (ANY(tmp_snwp>prlim)) STOP 'Non-zero snow even when disabled!'
    ENDIF
    IF (ALLOCATED(tmp_gasp)) THEN
        IF (ANY(tmp_gasp>prlim)) STOP 'Non-zero gas even when disabled!'
    ENDIF

  END SUBROUTINE SALSA_diagnostics
  !
  ! ---------------------------------------------------------------------
  ! sedim_SALSA: calculates sedimentation for all SALSA species
  ! Juha: The code below is a modified version of the original one by Zubair
  ! Jaakko: Modified for the use of ice and snow bins
  !
  SUBROUTINE sedim_SALSA()
    USE mo_submctl, ONLY : nlim,prlim
    use defs, only : alvl, alvi
    use grid, only : nzp, nxp, nyp, nspt, nbins, ncld, nprc, nice, nsnw,   &
                        sed_aero, sed_cloud, sed_precp, sed_ice, sed_snow, &
                        level, dtl, dzt, a_dn, a_temp, a_theta,            &
                        a_naerop,  a_naerot,  a_maerop,  a_maerot,         &
                        a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,        &
                        a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,        &
                        a_nicep,   a_nicet,   a_micep,   a_micet,          &
                        a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,         &
                        a_ustar, a_tt
    IMPLICIT NONE

    ! Ice only when level=5
    sed_snow = sed_snow .AND. level==5
    sed_ice = sed_ice .AND. level==5

    IF (sed_aero ) CALL DepositionAny(nzp,nxp,nyp,nspt,nbins,a_temp,a_theta,a_dn,a_ustar,dzt, &
                            a_naerop,a_maerop,a_naerot,a_maerot,a_tt,dtl,nlim,alvl,1)

    IF (sed_cloud) CALL DepositionAny(nzp,nxp,nyp,nspt,ncld,a_temp,a_theta,a_dn,a_ustar,dzt, &
                            a_ncloudp,a_mcloudp,a_ncloudt,a_mcloudt,a_tt,dtl,nlim,alvl,2)

    IF (sed_precp) CALL DepositionAny(nzp,nxp,nyp,nspt,nprc,a_temp,a_theta,a_dn,a_ustar,dzt, &
                            a_nprecpp,a_mprecpp,a_nprecpt,a_mprecpt,a_tt,dtl,prlim,alvl,3)

    IF (sed_ice  ) CALL DepositionAny(nzp,nxp,nyp,nspt,nice,a_temp,a_theta,a_dn,a_ustar,dzt, &
                            a_nicep,a_micep,a_nicet,a_micet,a_tt,dtl,prlim,alvi,4)

    IF (sed_snow ) CALL DepositionAny(nzp,nxp,nyp,nspt,nsnw,a_temp,a_theta,a_dn,a_ustar,dzt, &
                            a_nsnowp,a_msnowp,a_nsnowt,a_msnowt,a_tt,dtl,prlim,alvi,5)

  END SUBROUTINE !sedim_SALSA

  ! -----------------------------------------------------------------
  SUBROUTINE DepositionAny(n1,n2,n3,n4,nn,tk,th,adn,ustar,dzt,numc,mass,numct,masst,tlt,dt,clim,alv,flag)
    USE mo_submctl, ONLY : terminal_vel, calc_eff_radius
    use defs, only : pi, cp, kb, g
    use stat, only : sflg, SALSA_precip_stats
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,n4       ! Grid numbers, number of chemical species
    INTEGER, INTENT(in) :: nn                ! Number of bins
    REAL, INTENT(in) :: tk(n1,n2,n3)         ! Absolute temprature
    REAL, INTENT(in) :: th(n1,n2,n3)         ! Potential temprature
    REAL, INTENT(in) :: adn(n1,n2,n3)        ! Air density
    REAL, INTENT(in) :: ustar(n2,n3)         ! Friction velocity
    REAL, INTENT(in) :: dzt(n1)              ! Inverse of grid level thickness
    REAL, INTENT(in) :: numc(n1,n2,n3,nn)    ! Particle number concentration
    REAL, INTENT(in) :: mass(n1,n2,n3,nn*n4) ! Particle mass mixing ratio
    REAL, INTENT(in) :: dt                   ! timestep
    REAL, INTENT(IN) :: clim                ! Concentration limit (#/m^3)
    REAL, INTENT(IN) :: alv                 ! Specific latent heat of vaporization(J/kg)
    INTEGER, INTENT(IN) :: flag         ! An option for identifying aerosol, cloud, precipitation, ice and snow
    REAL, INTENT(out) :: numct(n1,n2,n3,nn)    ! Particle number concentration tendency
    REAL, INTENT(out) :: masst(n1,n2,n3,nn*n4) ! Particle mass mixing ratio tendency
    REAL, INTENT(out) :: tlt(n1,n2,n3) ! Ice-liquid water potential temperature tendency

    REAL :: flxdivm(n1,n2,n3,nn*n4), flxdivn(n1,n2,n3,nn) ! Mass and number divergence
    REAL :: depflxm(n1,n2,n3,nn*n4), depflxn(n1,n2,n3,nn) ! Mass and number deposition fluxes
    INTEGER :: i,j,k,kk
    INTEGER :: bin,bs

    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]
    real, parameter :: A = 1.249      ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87

    REAL :: avis, kvis           ! Viscosity of air, kinematic viscosity
    REAL :: lambda              ! Mean free path
    REAL :: va                    ! Thermal speed of air molecule
    REAL :: Kn, GG               ! Knudsen number, slip correction
    REAL :: vc                    ! Critical fall speed
    REAL :: mdiff                ! Particle diffusivity
    REAL :: rt, Sc, St

    REAL :: pmass(n4), rwet, fd, frac, pdn
    flxdivm = 0.
    flxdivn = 0.
    depflxm = 0.
    depflxn = 0.

    DO j = 3,n3-2
       DO i = 3,n2-2
          DO k=n1-1,2,-1

             ! atm modelling Eq.4.54
             avis=1.8325e-5*(416.16/(tk(k,i,j)+120.0))*(tk(k,i,j)/296.16)**1.5
             kvis =  avis/adn(k,i,j)
             va = sqrt(8*kb*tk(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2*avis/(adn(k,i,j)*va) !mean free path

             ! Fluxes
             !------------------
             DO bin = 1,nn
                IF (numc(k,i,j,bin)*adn(k,i,j)<clim) CYCLE

                ! Calculate wet size
                !   n4 = number of active species
                !   bin = size bin
                pmass(:)=mass(k,i,j,bin:(n4-1)*nn+bin:nn)/numc(k,i,j,bin)
                rwet=calc_eff_radius(n4,pmass,flag)

                ! Calculate effective density
                pdn = sum(pmass(:))/(4./3.*pi*rwet**3.)

                ! Terminal velocity
                Kn = lambda/rwet
                GG = 1.+ Kn*(A+B*exp(-C/Kn))
                vc = terminal_vel(rwet,pdn,adn(k,i,j),avis,GG,flag)

                IF (k==2) THEN ! The level just above surface
                    ! Particle diffusitivity  (15.29) in jacobson book
                    mdiff = (kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet)
                    Sc = kvis/mdiff
                    St = vc*ustar(i,j)**2.0/g*kvis
                    if (St<0.01) St=0.01
                    rt = 1.0/MAX(epsilon(1.0),(ustar(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18
                    vc = (1./rt) + vc
                ENDIF

                kk = k ! current interface - start from the first below current bin
                fd = vc*dt ! distance from the current interface to the leading edge of the falling particles
                DO WHILE ( fd > 0. .AND. kk > 1 )
                    ! Fraction of grid cell k passing through interface kk is fd*dzt(k), where nagtive values mean that the
                    ! interface is not reached (fd<0 ignored) and values larger than one mean that all particles fall through.
                    ! Absolute particle number flux per unit area and time is obtained by multiplying concentrations by
                    ! rho_air*dz/dt.
                    frac = MIN(1., fd*dzt(k) )*adn(k,i,j)/dzt(k)/dt ! kg/m^2/s

                    ! Flux for the particle mass (positive down)
                    DO bs = bin, (n4-1)*nn + bin, nn
                        depflxm(kk,i,j,bs) = depflxm(kk,i,j,bs) + mass(k,i,j,bs)*frac ! kg/m^2/s
                    END DO
                    ! Flux for the particle number
                    depflxn(kk,i,j,bin) = depflxn(kk,i,j,bin) + numc(k,i,j,bin)*frac ! #/m^2/s

                    ! The next interface
                    kk=kk-1
                    fd = fd - 1./dzt(kk)
                ENDDO
             END DO

             ! Divergency (here positive values mean decreasing mass)
             !  dc/dt=(dn/dt)/(rho*A*dz)=(dn_out/A/dt-dn_in/A/dt)/rho/dz
             flxdivm(k,i,j,:) = (depflxm(k,i,j,:)-depflxm(k+1,i,j,:))*dzt(k)/adn(k,i,j) ! kg/kg/s
             flxdivn(k,i,j,:) = (depflxn(k,i,j,:)-depflxn(k+1,i,j,:))*dzt(k)/adn(k,i,j) ! #/kg/s

          END DO ! k

          ! Assume constant flux at the domain top (level n1-2 is calculated correctly, so use that)
          flxdivm(n1,i,j,:) = flxdivm(n1-2,i,j,:)
          flxdivm(n1-1,i,j,:) = flxdivm(n1-2,i,j,:)
          flxdivn(n1,i,j,:) = flxdivn(n1-2,i,j,:)
          flxdivn(n1-1,i,j,:) = flxdivn(n1-2,i,j,:)

          ! Account for changes in ice-liquid water potential temperature
          DO k = 2,n1
             tlt(k,i,j) = tlt(k,i,j) + SUM(flxdivm(k,i,j,1:nn))*(alv/cp)*th(k,i,j)/tk(k,i,j)
          END DO ! k

       END DO ! i
    END DO ! j

    ! Number and mass tendencies
    numct = numct - flxdivn
    masst = masst - flxdivm

    ! Statistics
    IF (sflg) CALL SALSA_precip_stats(n1,n2,n3,n4,nn,depflxm,flag)

  END SUBROUTINE DepositionAny

end module step
