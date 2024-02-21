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

  logical :: lsvarflg = .false.

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
    use thrm, only : thermo

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

    use grid, only : level, dtl, Tspinup, sst, u0, v0, umean, vmean,                    &
                     ! Added parameters for interfacing with SALSA
                     nxp, nyp, nzp, a_press, a_temp, a_rp, a_rt, a_rsl, a_rsi, a_dn,    &
                     a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                     a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                     a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                     a_gaerop, a_gaerot,                                                &
                     nspec, nbins, ncld, nprc, nice, nsnw, &
                     nudge_theta, nudge_rv, nudge_u, nudge_v, nudge_ccn, &
                     ifSeaSpray, ifSeaVOC, sea_tspinup, a_edr

    use stat, only : sflg, statistics, les_rate_stats, out_mcrp_data, out_mcrp_list, out_mcrp_nout, mcrp_var_save
    use sgsm, only : diffuse
    use srfc, only : surface, marine_aero_flux, marine_gas_flux, srfc_gas_flux
    use thrm, only : thermo
    use mcrp, only : micro
    use prss, only : poisson
    use advf, only : fadvect
    use advl, only : ladvect
    use forc, only : forcings
    use lsvar, only : varlscale
    USE mo_salsa_driver, ONLY : run_SALSA
    USE mo_submctl, ONLY : nvbs_setup

    real :: xtime
    INTEGER :: zrm
    INTEGER :: n4

    xtime = time/86400. + strtim

    n4 = nspec + 1 ! Aerosol components + water

    ! The runmode parameter zrm is used by SALSA only
    zrm = 3
    IF ( time < Tspinup ) zrm = 2

    ! Reset ALL tendencies here.
    !----------------------------------------------------------------
    ! "Scalar" timestep
    CALL tend0(.FALSE.)

    ! Large scale forcing based on specified SST and geostrophic winds
    if (lsvarflg) THEN
        call varlscale(time,sst,u0,v0)
        u0(:) = u0(:) - umean
        v0(:) = v0(:) - vmean
    ENDIF

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
        CALL tend_constrain(n4)
    END IF

    IF (sflg) CALL les_rate_stats('srfc')
    call update_sclrs
    CALL tend0(.TRUE.)

    call diffuse

    IF (level>3) CALL tend_constrain(n4)
    IF (sflg) CALL les_rate_stats('diff')
    call update_sclrs
    CALL tend0(.TRUE.)

    call sponge(0)

    if (level >= 0) then

       call thermo(level)

       call forcings(xtime,cntlat,sst)

       IF (level>3) CALL tend_constrain(n4)
       IF (sflg) CALL les_rate_stats('forc')
       call update_sclrs

       IF (level >= 4) THEN

          CALL tend0(.TRUE.)

          CALL run_SALSA(nxp,nyp,nzp,n4,nbins,ncld,nprc,nice,nsnw, &
                  a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_dn,a_edr, &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_gaerop,  a_gaerot,  zrm, dtl, time, level, &
                  sflg, out_mcrp_nout, out_mcrp_list, out_mcrp_data)

          CALL tend_constrain(n4)
          IF (sflg) CALL les_rate_stats('mcrp')
          call update_sclrs

          ! Save user-selected details about SALSA microphysics
          IF (sflg) CALL mcrp_var_save()
       END IF

    end if ! level

    !-------------------------------------------
    ! "Deposition" timestep
    ! -- Reset only scalar tendencies
    CALL tend0(.TRUE.)

    ! Dont perform sedimentation or level 3 autoconversion during spinup
    IF (zrm == 3) THEN
        CALL micro(level)

        IF (level >= 4) CALL tend_constrain(n4)

        IF (sflg .AND. level < 4 ) CALL les_rate_stats('mcrp')
        IF (sflg .AND. level >= 4) CALL les_rate_stats('sedi')

        ! Save user-selected details about Seifert and Beheng microphysics
        IF (sflg .AND. level < 4) CALL mcrp_var_save()

        CALL update_sclrs
    END IF

    !-------------------------------------------
    ! "Advection" timestep
    ! -- Reset only scalar tendencies
    call tend0(.TRUE.)

    CALL fadvect

    IF (level >= 4) CALL tend_constrain(n4)

    IF (sflg) CALL les_rate_stats('advf')

    CALL update_sclrs

    CALL thermo(level)

    ! Nudging
     IF (nudge_theta/=0 .OR. nudge_rv/=0 .OR. nudge_u/=0 .OR. &
            nudge_v/=0 .OR. (level>3 .AND. nudge_ccn/=0) ) THEN

        ! Reset tendencies
        call tend0(.TRUE.)

        ! Update diagnostic tracers
        IF (level >= 4)  THEN
             CALL SALSA_diag_update
             call thermo(level)
        ENDIF

        CALL nudging(time)

        IF (level >= 4) CALL tend_constrain(n4)

        IF (sflg) CALL les_rate_stats('nudg')

        CALL update_sclrs

        CALL thermo(level)
    ENDIF

    IF (level >= 4)  THEN
         CALL SALSA_diagnostics(.true.)
         call thermo(level)
    ENDIF

    call corlos

    call ladvect

    call buoyancy

    call sponge(1)

    call poisson

    CALL thermo(level)

    IF (level >= 4)  THEN
         CALL SALSA_diagnostics(.false.)
         call thermo(level)
         IF (sflg) THEN
            CALL les_rate_stats('diag') ! This includes both calls
            call tend0(.TRUE.) ! Tendencies were calculated for this purpose only, so set to zero now
         ENDIF
    ENDIF

    if (sflg) then
       call statistics (time+dtl)
    end if

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
                zt, a_rp, a_rt, a_rc, a_srp, a_ri, a_srs, &
                a_naerop, a_naerot, a_ncloudp, a_nicep, &
                nspec, a_dn, a_maerop, a_maerot, &
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
        !   Levels 4-5: total = water vapor (a_rp) + aerosol and cloud water (a_rc) + rain water (a_srp)
        !                            + ice (a_ri) + snow (a_srs)
        IF (nudge_rv/=0)  THEN
            ALLOCATE(rv_ref(nzp))
            IF (level>3) THEN
                rv_ref(:)=a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)+a_ri(:,3,3)+a_srs(:,3,3)
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
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rc+a_srp+a_ri+a_srs,a_rt,rv_ref,dtl,time,nudge_rv, &
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
        IF (level==5) aero_target(:,:,:,in2a:nbins)=aero_target(:,:,:,in2a:nbins)-a_nicep(:,:,:,1:nice)
        ! Apply to sectional data
        CALL nudge_any_2d(nxp,nyp,nzp,nbins,zt,aero_target,a_naerot,aero_ref,dtl,time,nudge_ccn, &
            nudge_ccn_time,nudge_ccn_zmin,nudge_ccn_zmax,nudge_ccn_tau)
        ! Change aerosol mass so that the mean dry and wet sizes are unchanged
        CALL adj_salsa_maerot(nxp,nyp,nzp,nbins,nspec,nlim,a_dn,a_naerop,a_naerot,a_maerop,a_maerot)
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

    use grid, only : a_ut, a_vt, a_wt, nscl, a_st, newsclr

    LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

    integer :: n

    IF( .NOT. sclonly) THEN
       a_ut=0.; a_vt=0.; a_wt=0.
    ENDIF
    do n=1,nscl
       call newsclr(n)
       a_st=0.
    end do

  end subroutine tend0
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
                     a_nicep, a_nicet, a_nsnowp, a_nsnowt, a_micep, a_micet, a_msnowp, a_msnowt,  &
                     dtl, nxp,nyp,nzp,nbins,ncld,nprc,nice,nsnw,level

    INTEGER, INTENT(in) :: nn

    INTEGER :: cc, ii,jj,kk,ni

    DO jj = 3,nyp-2

       DO ii = 3,nxp-2

          DO kk = 1,nzp

             ! Aerosols
             DO cc = 1,nbins

                IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_naerot(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop(kk,ii,jj,cc))/dtl,a_naerot(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( ((1.e-10-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtl,  &
                                                               a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                   END DO

                END IF

             END DO

             ! Cloud droplets
             DO cc = 1,ncld

                IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_ncloudt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp(kk,ii,jj,cc))/dtl,a_ncloudt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtl,  &
                                                               a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                   END DO

                END IF

             END DO

             ! Precipitation
             DO cc = 1,nprc

                IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nprecpt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp(kk,ii,jj,cc))/dtl,a_nprecpt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtl,  &
                                                               a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                   END DO

                END IF

             END DO

             ! ice particles
             IF (level<5) CYCLE
             DO cc = 1,nice

                IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nicet(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep(kk,ii,jj,cc))/dtl,a_nicet(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_micet(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtl,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                   END DO

                END IF

             END DO

             ! Snow
             DO cc = 1,nsnw

                IF ( a_nsnowp(kk,ii,jj,cc)+a_nsnowt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nsnowt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp(kk,ii,jj,cc))/dtl,a_nsnowt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) = MAX( ((1.e-10-1.0)*a_msnowp(kk,ii,jj,(ni-1)*nsnw+cc))/dtl,  &
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
  subroutine update_sclrs

    use grid, only : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
         dtl, newsclr, isgstyp
    use sgsm, only : tkeinit
    use util, only : sclrset

    integer :: n

    do n=1,nscl
       call newsclr(n)
       call update(nzp,nxp,nyp,a_sp,a_st,dtl)
       call sclrset('cnst',nzp,nxp,nyp,a_sp,dzt)
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
         a_rp, a_rpp, a_srp, a_ri, a_srs, nxp, nyp, nzp, dzm, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1

    real :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)

    IF (level<4) THEN
       rv = a_rv ! Water vapor
       rc = a_rc + a_rpp + a_ri ! Total condensate (cloud + precipitation + total ice)
    ELSE
       rv = a_rp ! Water vapor
       rc = a_rc + a_srp + a_ri + a_srs ! Total condensed water (aerosol+cloud+precipitation+ice+snow)
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
         nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

    integer, intent (in) :: isponge

    integer :: i, j, k, kk

    if (maxval(spng_tfct) > epsilon(1.) .and. nfpt > 1) then
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
    end if

  end subroutine sponge

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

  SUBROUTINE SALSA_diagnostics(reset_stats)
    USE grid, ONLY : nxp,nyp,nzp,nbins,ncld,nprc,nice,nsnw,nspec, &
                     a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,      &
                     a_nicep,a_micep,a_nsnowp,a_msnowp,a_gaerop, &
                     a_naerot,a_maerot,a_ncloudt,a_mcloudt,a_nprecpt,a_mprecpt,      &
                     a_nicet,a_micet,a_nsnowt,a_msnowt,a_gaerot, &
                     a_rp, a_rsl, a_rsi, a_temp, a_dn, level, dtl, &
                     tmp_prcp, tmp_icep, tmp_snwp, tmp_gasp
    USE mo_submctl, ONLY : fn1a,in2a,fn2a, &
                     diss, mws, dens, &
                     surfw0, rg, nlim, prlim, pi, pi6, &
                     lscndgas, part_h2so4, part_ocnv, iso, ioc, isog, iocg, &
                     aerobins, calc_correlation
    USE stat, ONLY : sflg
    IMPLICIT NONE

    LOGICAL :: reset_stats

    INTEGER :: i,j,k,bc,ba,bb,s,sc,sa,nn

    REAL :: zvol, ra, rb, a_rh(nzp,nxp,nyp), a_rhi(nzp,nxp,nyp)
    REAL :: ns, cd

    a_rh(:,3:nxp-2,3:nyp-2) = a_rp(:,3:nxp-2,3:nyp-2)/a_rsl(:,3:nxp-2,3:nyp-2)
    a_rhi(:,3:nxp-2,3:nyp-2) = a_rp(:,3:nxp-2,3:nyp-2)/a_rsi(:,3:nxp-2,3:nyp-2)

    ! Change in concentrations due to diagnostics (e.g. release of cloud/rain/ice/snow into aerosol,
    ! cleaning particles without mass or negligible number concentration)
    IF (sflg .AND. reset_stats) THEN
        ! Statistics output step: need to calculate tendencies over
        ! both calls, so use tendencies as temporary variables!
        a_naerot=a_naerop; a_maerot=a_maerop
        a_ncloudt=a_ncloudp; a_mcloudt=a_mcloudp
        a_nprecpt=a_nprecpp; a_mprecpt=a_mprecpp
        a_nicet=a_nicep; a_micet=a_micep
        a_nsnowt=a_nsnowp; a_msnowt=a_msnowp
        a_gaerot=a_gaerop
    ENDIF

    nn = nspec + 1 ! Aerosol components + water

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

    ! Remove particles that have number but no mass.
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp
             ! Aerosols
             DO bc = 1,nbins
                IF (a_naerop(k,i,j,bc) > 0. .AND. SUM(a_maerop(k,i,j,bc:nspec*nbins+bc:nbins)) <= 0.) THEN
                   a_naerop(k,i,j,bc) = 0.
                   a_maerop(k,i,j,bc:nspec*nbins+bc:nbins) = 0.
                END IF
             END DO

             ! Clouds
             DO bc = 1,ncld
                IF (a_ncloudp(k,i,j,bc) > 0. .AND. SUM(a_mcloudp(k,i,j,bc:nspec*ncld+bc:ncld)) <= 0.) THEN
                   a_ncloudp(k,i,j,bc) = 0.
                   a_mcloudp(k,i,j,bc:nspec*ncld+bc:ncld) = 0.
                END IF
             END DO ! ncld

             ! Precipitation
             DO bc = 1,nprc
                IF (a_nprecpp(k,i,j,bc) > 0. .AND. a_mprecpp(k,i,j,bc) <= 0.) THEN
                   a_nprecpp(k,i,j,bc) = 0.
                   a_mprecpp(k,i,j,bc:nspec*nprc+bc:nprc) = 0.
                END IF
             END DO ! nprc

             ! Ice
             IF (level<5) CYCLE
             DO bc = 1,nice
                IF (a_nicep(k,i,j,bc) > 0. .AND. SUM(a_micep(k,i,j,bc:nspec*nice+bc:nice)) <= 0.) THEN
                   a_nicep(k,i,j,bc) = 0.
                   a_micep(k,i,j,bc:nspec*nice+bc:nice) = 0.
                END IF
             END DO ! ncld

             ! Snow
             DO bc = 1,nsnw
                IF (a_nsnowp(k,i,j,bc) > 0. .AND. a_msnowp(k,i,j,bc) <= 0.) THEN
                   a_nsnowp(k,i,j,bc) = 0.
                   a_msnowp(k,i,j,bc:nspec*nsnw+bc:nsnw) = 0.
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
             DO bc = 1,ncld
                IF ( a_ncloudp(k,i,j,bc)*a_dn(k,i,j) > nlim .AND. a_rh(k,i,j)<0.999 .AND. &
                        a_mcloudp(k,i,j,bc)<1e-5 ) THEN
                   ! Critical diameter (assuming soluble CCN)
                   ns = SUM( diss(2:nn)*a_mcloudp(k,i,j,ncld+bc:nspec*ncld+bc:ncld)/mws(2:nn) )/a_ncloudp(k,i,j,bc)
                   cd = 3.*SQRT(ns*rg*a_temp(k,i,j)/(2.*pi*surfw0))

                   ! Wet diameter
                   zvol = (SUM( a_mcloudp(k,i,j,bc:nspec*ncld+bc:ncld)/dens(1:nn) )/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.2*critical diameter or 2 um or if there is no water
                   IF ( zvol < MAX(0.2*cd,2.e-6) .OR. a_mcloudp(k,i,j,bc)<1e-25*a_ncloudp(k,i,j,bc) ) THEN
                      ba = fn1a + bc ! Index for parallel aerosol bin
                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = (s-1)*ncld + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO
                   END IF ! critical diameter
                END IF  ! blim
             END DO ! bc

             ! Loop over precipitation bins
             DO bc = 1,nprc
                IF ( a_nprecpp(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rh(k,i,j)<0.999 .AND. &
                        a_mprecpp(k,i,j,bc)<1e-6 ) THEN
                   ! Critical diameter
                   ns = SUM( diss(2:nn)*a_mprecpp(k,i,j,nprc+bc:nspec*nprc+bc:nprc)/mws(2:nn) )/a_nprecpp(k,i,j,bc)
                   cd = 3.*SQRT(ns*rg*a_temp(k,i,j)/(2.*pi*surfw0))

                   ! Wet diameter
                   zvol = (SUM( a_mprecpp(k,i,j,bc:nspec*nprc+bc:nprc)/dens(1:nn) )/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.02*critical diameter or 2 um or if there is no water
                   IF ( zvol < MAX(0.02*cd,2.e-6) .OR. a_mprecpp(k,i,j,bc)<1e-25*a_nprecpp(k,i,j,bc) ) THEN

                      ! Move evaporating precipitation to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle radius (a and b bins)
                      cd = 0.5*(SUM( a_mprecpp(k,i,j,nprc+bc:nspec*nprc+bc:nprc)/dens(2:nn) )/a_nprecpp(k,i,j,bc)/pi6)**(1./3.) ! Dry radius
                      ba=fn2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd<aerobins(ba) .AND. ba>in2a)
                         ba=ba-1
                      ENDDO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (bb>nbins) THEN
                         ! b bins not used so select a
                         !ba = ba
                      ELSEIF (a_naerop(k,i,j,bb)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSEIF (a_naerop(k,i,j,ba)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra=calc_correlation(a_maerop(k,i,j,nbins+ba:nspec*nbins+ba:nbins), &
                                a_mprecpp(k,i,j,nprc+bc:nspec*nprc+bc:nprc),nspec)
                         rb=calc_correlation(a_maerop(k,i,j,nbins+bb:nspec*nbins+bb:nbins), &
                                a_mprecpp(k,i,j,nprc+bc:nspec*nprc+bc:nprc),nspec)
                         IF (ra<rb) ba = bb
                      ENDIF

                      ! Move the number of particles from precipitation to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                      a_nprecpp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = (s-1)*nprc + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                         a_mprecpp(k,i,j,sc) = 0.
                      END DO
                   END IF ! Critical diameter
                END IF ! prlim
             END DO ! bc

             ! Loop over ice bins
             DO bc = 1,nice
                IF ( a_nicep(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rhi(k,i,j)<0.999 .AND. &
                        a_micep(k,i,j,bc)<1e-8 ) THEN
                   ! Diameter (assuming water density for ice)
                   cd = ( SUM(a_micep(k,i,j,bc:nspec*nice+bc:nice)/dens(1:nn))/a_nicep(k,i,j,bc)/pi6)**(1./3.)

                   ! Dry to total mass ratio
                   zvol = SUM( a_micep(k,i,j,nice+bc:nspec*nice+bc:nice) )/SUM( a_micep(k,i,j,bc:nspec*nice+bc:nice) )

                   ! Ice and snow don't have a critical size, but lose particles smaller than 2e-6 m and particles which dry to total mass ratio is more than 0.5
                   IF ( zvol>0.5 .OR. cd<2e-6 ) THEN
                      ba = fn1a + bc ! Index for parallel aerosol bin
                      ! Move the number of particles from ice to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move mass to aerosol (including water)
                      DO s = 1,nn
                         sc = (s-1)*nice + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO
                   END IF
                END IF  ! prlim
             END DO ! bc

             ! Loop over snow bins
             DO bc = 1,nsnw
                IF ( a_nsnowp(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rhi(k,i,j)<0.999 .AND. &
                        a_msnowp(k,i,j,bc)<1e-8 ) THEN
                   ! Diameter (assuming water density for snow)
                   cd = ( SUM(a_msnowp(k,i,j,bc:nspec*nsnw+bc:nsnw)/dens(1:nn))/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)

                   ! Dry to total mass ratio
                   zvol = SUM( a_msnowp(k,i,j,nsnw+bc:nspec*nsnw+bc:nsnw) )/SUM( a_msnowp(k,i,j,bc:nspec*nsnw+bc:nsnw) )

                   ! Lose particles smaller than 2e-6 m and particles which dry to total mass ratio is more than 0.5
                   IF ( zvol>0.5  .OR. cd<2.e-6 ) THEN

                      ! Move evaporating snow to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle radius (a and b bins)
                      cd = 0.5*(SUM( a_msnowp(k,i,j,nsnw+bc:nspec*nsnw+bc:nsnw)/dens(2:nn) )/a_nsnowp(k,i,j,bc)/pi6)**(1./3.) ! Dry radius
                      ba=fn2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd<aerobins(ba) .AND. ba>in2a)
                         ba=ba-1
                      ENDDO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (bb>nbins) THEN
                         ! b bins are not used so select a
                         !ba = ba
                     ELSEIF (a_naerop(k,i,j,bb)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSEIF (a_naerop(k,i,j,ba)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra=calc_correlation(a_maerop(k,i,j,nbins+ba:nspec*nbins+ba:nbins), &
                                a_msnowp(k,i,j,nsnw+bc:nspec*nsnw+bc:nsnw),nspec)
                         rb=calc_correlation(a_maerop(k,i,j,nbins+bb:nspec*nbins+bb:nbins), &
                                a_msnowp(k,i,j,nsnw+bc:nspec*nsnw+bc:nsnw),nspec)
                         IF (ra<rb) ba = bb
                      ENDIF

                      ! Move the number of particles from snow to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nsnowp(k,i,j,bc)
                      a_nsnowp(k,i,j,bc) = 0.

                      ! Move mass to aerosol (including water)
                      DO s = 1,nn
                         sc = (s-1)*nsnw + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_msnowp(k,i,j,sc)
                         a_msnowp(k,i,j,sc) = 0.
                      END DO
                   END IF
                END IF ! prlim
             END DO ! bc

             ! Loop over aerosol bins
             DO ba = 1,nbins
                IF (a_naerop(k,i,j,ba)*a_dn(k,i,j) > nlim) THEN
                   zvol = SUM( a_maerop(k,i,j,nbins+ba:nspec*nbins+ba:nbins)/dens(2:nn) )/a_naerop(k,i,j,ba) ! Dry volume

                   ! Particles smaller than 0.1 nm diameter are set to zero
                   IF ( zvol < pi6*1.e-10**3 ) THEN
                      ! Volatile species to the gas phase
                      IF (lscndgas .AND. part_h2so4) THEN
                         s = (iso-1)*nbins + ba
                         a_gaerop(k,i,j,isog) = a_gaerop(k,i,j,isog) + a_maerop(k,i,j,s)
                      END IF
                      IF (lscndgas .AND. part_ocnv) THEN
                         s = (ioc-1)*nbins + ba
                         a_gaerop(k,i,j,iocg) = a_gaerop(k,i,j,iocg) + a_maerop(k,i,j,s)
                      END IF

                      ! Mass and number to zero (insolube species and water are lost)
                      a_maerop(k,i,j,ba:nspec*nbins+ba:nbins) = 0.
                      a_naerop(k,i,j,ba) = 0.
                   END IF
                END IF
             END DO

          END DO   ! k
       END DO   ! i
    END DO   ! j

    IF (sflg .AND. .NOT.reset_stats) THEN
        a_naerot=(a_naerop-a_naerot)/dtl
        a_maerot=(a_maerop-a_maerot)/dtl
        a_ncloudt=(a_ncloudp-a_ncloudt)/dtl
        a_mcloudt=(a_mcloudp-a_mcloudt)/dtl
        a_nprecpt=(a_nprecpp-a_nprecpt)/dtl
        a_mprecpt=(a_mprecpp-a_mprecpt)/dtl
        a_nicet=(a_nicep-a_nicet)/dtl
        a_micet=(a_micep-a_micet)/dtl
        a_nsnowt=(a_nsnowp-a_nsnowt)/dtl
        a_msnowt=(a_msnowp-a_msnowt)/dtl
        a_gaerot=(a_gaerop-a_gaerot)/dtl
    ENDIF

    !!!!!!!!!!!!!!!!!!!!!!!
    ! Update diagnostic tracers
    !!!!!!!!!!!!!!!!!!!!!!!

    CALL SALSA_diag_update


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
  ! SALSA_diag_update: Update diagnostic concentration tracers

  SUBROUTINE SALSA_diag_update
    USE grid, ONLY : a_maerop,a_mcloudp,a_mprecpp,a_nprecpp, &
            a_micep,a_msnowp,a_nsnowp,nbins,ncld,nprc,nice,nsnw,&
            a_rc,a_srp,a_snrp,a_ri,a_srs,a_snrs

    ! Liquid water content (water is the first species)
    ! Aerosol and cloud droplets, regimes a and b
    a_rc(:,:,:) = SUM(a_maerop(:,:,:,1:nbins),DIM=4) + &
                  SUM(a_mcloudp(:,:,:,1:ncld),DIM=4)
    ! Precipitation
    a_srp(:,:,:) = SUM(a_mprecpp(:,:,:,1:nprc),DIM=4)
    a_snrp(:,:,:) = SUM(a_nprecpp(:,:,:,:),DIM=4)

    ! ice, regimes a and b
    a_ri(:,:,:) = SUM(a_micep(:,:,:,1:nice),DIM=4)
    ! Snow
    a_srs(:,:,:) = SUM(a_msnowp(:,:,:,1:nsnw),DIM=4)
    a_snrs(:,:,:) = SUM(a_nsnowp(:,:,:,:),DIM=4)

  END SUBROUTINE SALSA_diag_update


end module step
