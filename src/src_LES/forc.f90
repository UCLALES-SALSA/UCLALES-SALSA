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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE forc

  USE grid, ONLY: nxp,nyp,nzp,iradtyp,lnudging,lemission,  &
                  CCN,level,dtlt
  USE mo_vector_state, ONLY : a_ut, a_up, a_vt, a_vp
  USE mo_field_state, ONLY : SALSA_tracers_4d                  
  USE mo_aux_state, ONLY : zm,zt,dzt,dzm,dn0,pi0,pi1
  USE mo_diag_state, ONLY : a_pexnr,a_temp,a_rv,a_rc,a_rflx,a_sflx,   &
                            a_fus,a_fds,a_fuir,a_fdir,albedo
  USE mo_progn_state, ONLY : a_tt,a_tp,a_rt,a_rp,a_rpp, a_npp,        &
                             a_ncloudp,a_nprecpp,a_mprecpp,a_nicep,   &
                             a_chargeTimet, a_chargeTimep  
  USE mpi_interface, ONLY : myid, appl_abort
  USE util, ONLY : get_avg2dh
  USE defs, ONLY      : cp
  !USE stat, ONLY      : sflg
  USE radiation_main, ONLY : rad_interface, useMcICA, iradtyp
  USE nudg, ONLY : nudging
  USE emission_main, ONLY : aerosol_emission
  USE emission_types, ONLY : emitModes
  USE mo_structured_datatypes
  IMPLICIT NONE

  ! these are now all namelist parameters
  CHARACTER (len=10) :: case_name = 'none'               

  REAL    :: div = 0.

CONTAINS
  !
  ! -------------------------------------------------------------------
  ! Subroutine forcings:  calls the appropriate large-scale forcings.
  !
  ! Contains: explicit and parameterized radiation calculations, nudging,
  ! large-scale divergence effects and other case-specific forcings.
  !
  SUBROUTINE forcings(time,strtim)
    
    REAL,  INTENT (in) :: time, strtim  ! time in seconds (since model start), strtim in decimal days
    REAL :: xka, fr0, fr1
    REAL :: time_decday

    ! NOT FINISHED; PUT LARGE-SCALE FORCINGS/CASE-SPECIFIC STUFF IN THEIR OWN PACKAGES

    time_decday = time/86400. + strtim

    ! DIVERGENCE GIVEN FROM NAMELIST
    IF (trim(case_name) == 'atex') THEN
       xka = 130.
       fr0 = 74.
       fr1 = 0.
       div = 0.
    ELSE
       xka = 85.
       fr0 = 70.
       fr1 = 22.
       !div = 3.75e-6
    END IF

    !IF (trim(case_name) == 'ascos') THEN
    !    ! Full radiation calculations when saving data (stat/sflg=.TRUE. when saving)
    !    useMcICA = .NOT. sflg
    !END IF

    ! 
    ! Nudging
    ! ------------
    !
    IF (lnudging) CALL nudging(time)

    !
    ! Aerosol emissions
    ! --------------------
    !
    IF (lemission .AND. level >= 4) THEN
       CALL aerosol_emission(time)

       IF (ANY(emitModes(:)%emitType > 0)) THEN
          ! For particle charging, reduce the time tracer by one timestep where it is > 0
          a_chargeTimet%d = a_chargeTimet%d - MIN(dtlt, a_chargeTimep%d)
       END IF
    END IF
    
    SELECT CASE(iradtyp)
    CASE (1)
       ! No radiation, just large-scale forcing. 
       ! Note, there's a slight discrepancy between lev 1-3 and lev 4 with a_rp
       ! (total water vs vapour): it perhaps doesn't make sence to change the
       ! tendency of condesated water due to subsidence in level4 and for level
       ! 1-3 total mixing ratio is the only prognostic water variable.
       ! -------------------------------------------------
       IF ( case_name /= 'none' ) THEN
          CALL case_forcing(nzp,nxp,nyp,zt,dzt,dzm,div,a_tp,a_rp,a_tt,a_rt)
       END IF

    CASE (2)
       ! Some original special CASEs....
       ! ---------------------------------
       IF ( level >= 4) THEN
          IF(myid == 0) WRITE(*,*) 'FORCING: selection not implemented for level 4 or 5, iradtyp = ',iradtyp
          CALL appl_abort(0)
       END IF

       SELECT CASE(level)
       CASE(1) 
          CALL smoke_rad(nzp, nxp, nyp, dn0, a_rflx, zm, dzt,a_tt,a_rp)
       CASE(2)
          CALL gcss_rad(nzp, nxp, nyp, xka, fr0, fr1, div, a_rc, dn0,     &
                        a_rflx, zt, zm, dzt, a_tt, a_tp, a_rt, a_rp)

       END SELECT
       IF (trim(case_name) == 'atex') CALL case_forcing(nzp, nxp, nyp,    &
                                                        zt, dzt, dzm, div, a_tp, a_rp, a_tt, a_rt)
    CASE (3)
       ! Radiation + large scale forcing
       !----------------------------------
       CALL rad_interface(time_decday)
       
       IF ( case_name /= 'none') THEN
          CALL case_forcing(nzp,nxp,nyp,zt,dzt,dzm,div,a_tp,a_rp,a_tt,a_rt)
       END IF

    CASE (4)
       ! ??
       !---
       CALL bellon(nzp, nxp, nyp, a_rflx, a_sflx, zt, dzt, dzm, a_tt, a_tp,&
                   a_rt, a_rp, a_ut%d, a_up%d, a_vt%d, a_vp%d)
    END SELECT 

  END SUBROUTINE forcings

  !
  ! -------------------------------------------------------------------
  ! Subroutine gcss_rad:  call simple radiative parameterization and
  ! simultaneously update fields due to vertical motion as given by div
  !
  SUBROUTINE gcss_rad(n1,n2,n3,xka,fr0,fr1,div,rc,dn0,flx,zt,zm,dzt,   &
                      tt,tl,rtt,rt)

    INTEGER, INTENT (in)             :: n1, n2, n3
    REAL, INTENT (in)                :: xka, fr0, fr1, div
    TYPE(FloatArray1d), INTENT(in)   :: zt, zm, dzt, dn0
    TYPE(FloatArray3d), INTENT(in)   :: rc, tl, rt
    TYPE(FloatArray3d), INTENT(inout) :: tt, rtt
    TYPE(FloatArray3d), INTENT (out)  :: flx

    INTEGER :: i, j, k, km1, kp1, ki
    REAL    :: lwp(n2,n3), fact

    lwp = 0.
    DO j = 3, n3-2
       DO i = 3, n2-2
          ki = n1
          DO k = 1, n1
             km1 = max(1,k-1)
             lwp(i,j) = lwp(i,j)+max(0.,rc%d(k,i,j)*dn0%d(k)*(zm%d(k)-zm%d(km1)))
             flx%d(k,i,j) = fr1*exp(-1.*xka*lwp(i,j))
             IF ( (rc%d(k,i,j) > 0.01e-3) .AND. (rt%d(k,i,j) >= 0.008) ) ki = k
          END DO

          fact = div*cp*dn0%d(ki)
          DO k = 2, n1
             km1 = max(2,k-1)
             lwp(i,j) = lwp(i,j)-max(0.,rc%d(k,i,j)*dn0%d(k)*(zm%d(k)-zm%d(k-1)))
             flx%d(k,i,j) = flx%d(k,i,j)+fr0*exp(-1.*xka*lwp(i,j))
             IF (zm%d(k) > zm%d(ki) .AND. ki > 1 .AND. fact > 0.) THEN
                flx%d(k,i,j) = flx%d(k,i,j) + fact*(0.25*(zm%d(k)-zm%d(ki))**1.333 + &
                             zm%d(ki)*(zm%d(k)-zm%d(ki))**(1./3.))
             END IF
             tt%d(k,i,j) = tt%d(k,i,j)-(flx%d(k,i,j)-flx%d(km1,i,j))*dzt%d(k)/(dn0%d(k)*cp)
          END DO
          !
          ! subsidence
          !
          IF (div /= 0.) THEN
             DO k = 2, n1-2
                kp1 = k+1
                tt%d(k,i,j)  = tt%d(k,i,j) + &
                             div*zt%d(k)*(tl%d(kp1,i,j)-tl%d(k,i,j))*dzt%d(k)
                rtt%d(k,i,j) = rtt%d(k,i,j) + &
                             div*zt%d(k)*(rt%d(kp1,i,j)-rt%d(k,i,j))*dzt%d(k)
             END DO
          END IF
       END DO
    END DO

  END SUBROUTINE gcss_rad
  !
  ! -------------------------------------------------------------------
  ! Subroutine smoke_rad:  call simple radiative parameterization for
  ! the smoke cloud
  !
  SUBROUTINE smoke_rad(n1,n2,n3,dn0,flx,zm,dzt,tt,rt)

    INTEGER, INTENT (in)              :: n1,n2, n3
    TYPE(FloatArray1d), INTENT(in)    :: zm,dzt,dn0
    TYPE(FloatArray3d), INTENT(in)    :: rt
    TYPE(FloatArray3d), INTENT(inout) :: tt
    TYPE(FloatArray3d), INTENT(out)   :: flx
    REAL, PARAMETER      :: xka = 50.0, fr0 = 60.0

    INTEGER :: i,j,k, km1, ki
    REAL    :: smoke(n2,n3)

    smoke = 0.
    DO j = 3, n3-2
       DO i = 3, n2-2
          ki = n1
          DO k = 1, n1
             km1 = max(1,k-1)
             smoke(i,j) = smoke(i,j)+max(0.,rt%d(k,i,j)*dn0%d(k)*(zm%d(k)-zm%d(km1)))
          END DO

          DO k = 2, n1
             km1 = max(2,k-1)
             smoke(i,j) = smoke(i,j)-max(0.,rt%d(k,i,j)*dn0%d(k)*(zm%d(k)-zm%d(k-1)))
             flx%d(k,i,j) = fr0*exp(-1.*xka*smoke(i,j))
             tt%d(k,i,j) = tt%d(k,i,j)-(flx%d(k,i,j)-flx%d(km1,i,j))*dzt%d(k)/(dn0%d(k)*cp)
          END DO
       END DO
    END DO

  END SUBROUTINE smoke_rad
  !
  ! -------------------------------------------------------------------
  ! Subroutine case_forcing: adjusts tendencies according to a specified
  ! large scale forcing.  Normally CASE (run) specific.
  ! 
  ! BTW: None of the arguments in the subroutine call are really necessary, since they are all imported 
  !      globally to this module.
  SUBROUTINE case_forcing(n1,n2,n3,zt,dzt,dzm,zdiv,tl,rt,tt,rtt)

    USE mpi_interface, ONLY : pecount, double_scalar_par_sum,myid, appl_abort
    !USE stat, ONLY : get_zi

    INTEGER, INTENT (in) :: n1,n2, n3
    TYPE(FloatArray1d), INTENT (in)      :: zt, dzt, dzm
    REAL, INTENT(in)                     :: zdiv
    TYPE(FloatArray3d), INTENT (in)      :: tl, rt
    TYPE(FloatArray3d), INTENT (inout)   :: tt, rtt

    TYPE(FloatArray4d), POINTER :: varp => NULL(), vart => NULL()
    
    INTEGER :: i,j,k,kp1,b,c
    REAL, DIMENSION (n1) :: sf
    REAL, PARAMETER :: zmx_sub = 2260. ! originally 2260.

    REAL (kind=8) :: zig, zil
    REAL          :: zibar

    zig = 0.0; zil = 0.0; zibar = 0.0
    kp1 = 0
    sf(:) = -zdiv*zt%d(:)*dzt%d(:)
    SELECT CASE (trim(case_name))
    CASE('default')
       !
       ! User specified divergence used as a simple large scle forcing for moisture and temperature fields
       ! + also aerosol for level > 4
       ! -------------------------------------------------------------------------------------------------
       !
       DO j = 3, n3-2
          DO i = 3, n2-2
             DO k = 2, n1-1
                kp1 = k+1
                tt%d(k,i,j) = tt%d(k,i,j) - (tl%d(kp1,i,j)-tl%d(k,i,j))*sf(k)
                rtt%d(k,i,j) = rtt%d(k,i,j) - (rt%d(kp1,i,j)-rt%d(k,i,j))*sf(k)
             END DO
          END DO
       END DO

       ! Some additional stuff needed for SALSA. a_salsa array has all the necessary tracers
       ! and its association depends already on level, so no need to any extra checks here.
       IF (level >= 4) THEN

          DO b = 1,SALSA_tracers_4d%count
             CALL SALSA_tracers_4d%getData(1,varp,index=b)
             CALL SALSA_tracers_4d%getData(2,vart,index=b)
             DO c = 1,SIZE(varp%d,DIM=4)
                DO j = 3,n3-2
                   DO i = 3,n2-2
                      vart%d(2:n1-1,i,j,c) = vart%d(2:n1-1,i,j,c) -   &
                           (varp%d(3:n1,i,j,c) - varp%d(2:n1-1,i,j,c))*sf(2:n1-1)
                   END DO
                END DO
             END DO
          END DO

          varp => NULL(); vart => NULL()
          
       END IF

    CASE('rico')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       DO k = 2, n1-2
          IF (zt%d(k) < zmx_sub) THEN
             sf(k) = -0.005*zt%d(k)/zmx_sub
          ELSE
             sf(k) = -0.005
          END IF
          sf(k) = sf(k)*dzt%d(k)
       END DO

       DO j = 3, n3-2
          DO i = 3, n2-2
             DO k = 2, n1-2
                !
                ! subsidence
                ! 
                kp1 = k+1
                tt%d(k,i,j)  = tt%d(k,i,j) - ( tl%d(kp1,i,j) - tl%d(k,i,j) )*sf(k)
                rtt%d(k,i,j) = rtt%d(k,i,j) - ( rt%d(kp1,i,j) - rt%d(k,i,j) )*sf(k)
                !
                ! temperature advection and radiative cooling
                !
                tt%d(k,i,j) = tt%d(k,i,j)  - 2.5/86400.
                !
                ! moisture advection
                !
                IF (zt%d(k) <= 2980.) THEN
                   rtt%d(k,i,j) = rtt%d(k,i,j) - (1. - 1.3456*zt%d(k)/2980.)/8.64e7
                ELSE
                   rtt%d(k,i,j) = rtt%d(k,i,j) + .3456/8.64e7
                END IF
             END DO
          END DO
       END DO

    CASE ('bomex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       DO k = 2, n1-2
          IF (zt%d(k) < 1500.) THEN
             sf(k) = -0.0065*zt%d(k)/1500.
          ELSE
             sf(k) = min(0.,-0.0065  + 0.0065*(zt%d(k)-1500.)/600.)
          END IF
          sf(k) = sf(k)*dzt%d(k)
       END DO

       DO j = 3, n3-2
          DO i = 3, n2-2
             DO k = 2, n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                IF (zt%d(k) < 1500.) THEN
                   tt%d(k,i,j) = tt%d(k,i,j) - ( tl%d(kp1,i,j)-tl%d(k,i,j) )*sf(k) &
                              - 2.315e-5
                ELSE IF (zt%d(k) < 2000.) THEN
                   tt%d(k,i,j) = tt%d(k,i,j) - ( tl%d(kp1,i,j)-tl%d(k,i,j) )*sf(k) &
                              - 2.315e-5*(1.- (zt%d(k)-1500.)*1.e-3)
                END IF
                !
                ! moisture advection
                !
                rtt%d(k,i,j) = rtt%d(k,i,j) - ( rt%d(kp1,i,j) - rt%d(k,i,j) )*sf(k)
                IF (zt%d(k) < 300.) THEN
                   rtt%d(k,i,j) = rtt%d(k,i,j) - 1.2e-8
                ELSE IF (zt%d(k) < 500.) THEN
                   rtt%d(k,i,j) = rtt%d(k,i,j) - 1.2e-8*(1.- (zt%d(k)-300.)/200.)
                END IF
             END DO
          END DO
       END DO
    CASE ('atex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       zil = get_zi (n1, n2, n3, 2, rt, dzm, zt, 6.5e-3)
       CALL double_scalar_par_sum(zil,zig)
       zibar = REAL(zig/pecount)

       DO k = 2, n1-2
          IF (zt%d(k) < zibar) THEN
             sf(k) = -0.0065*zt%d(k)/1500.
          ELSE
             sf(k) = min(0.,-0.0065*(1 - (zt%d(k)-zibar)/300.))
          END IF
          sf(k) = sf(k)*dzt%d(k)
       END DO

       DO j = 3, n3-2
          DO i = 3, n2-2
             DO k = 2, n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                IF (zt%d(k) < zibar) THEN
                   tt%d(k,i,j) = tt%d(k,i,j) - ( tl%d(kp1,i,j)-tl%d(k,i,j) )*sf(k) &
                              - 2.315e-5*(1. + (1.- zt%d(k)/zibar)/2.)
                ELSE IF (zt%d(k) < zibar+300.) THEN
                   tt%d(k,i,j) = tt%d(k,i,j) - ( tl%d(kp1,i,j)-tl%d(k,i,j) )*sf(k) &
                              - 2.315e-5*(1.- (zt%d(k)-zibar)/300.)
                END IF
                !
                ! moisture advection
                !
                rtt%d(k,i,j) = rtt%d(k,i,j) - ( rt%d(kp1,i,j) - rt%d(k,i,j) )*sf(k)
                IF (zt%d(k) < zibar) rtt%d(k,i,j) = rtt%d(k,i,j)  - 1.5e-8
             END DO
          END DO
       END DO
        !
    CASE ('ascos')
        ! ASCOS
        ! ---------
        !
        DO k = 2, n1-2
            ! calculate subsidence factor (wsub / dz)
            sf(k) = -5.0e-6*min(2000.0,zt%d(k))*dzt%d(k)
        END DO
        !
        DO j = 3, n3-2
            DO i = 3, n2-2
                DO k = 2, n1-2
                    !
                    ! Temperature and humidity advection due to subsidence
                    !
                    kp1 = k+1
                    tt%d(k,i,j)  =  tt%d(k,i,j) - ( tl%d(kp1,i,j) - tl%d(k,i,j) )*sf(k)
                    rtt%d(k,i,j) = rtt%d(k,i,j) - ( rt%d(kp1,i,j) - rt%d(k,i,j) )*sf(k)
                END DO
            END DO
        END DO
        !
    CASE ('barba')
        ! Barbados
        ! -----------
        ! Large scale subsidence: w(z)=w0*(1-exp(z/H)), where w0=7.5 mm/s and H=1000 m.
        ! Radiative cooling rate: 2.5 K/day
        ! No temperature or humidity advection
        !
        ! calculate subsidence factor (wsub / dz)
        DO k = 2, n1-2
            sf(k) = -7.5e-3*(1.0-exp(-zt%d(k)/1000.0))*dzt%d(k)
        END DO
        !
        DO j = 3, n3-2
            DO i = 3, n2-2
                DO k = 2, n1-2
                    ! Subsidence
                    kp1 = k+1
                    tt%d(k,i,j)  =  tt%d(k,i,j) - ( tl%d(kp1,i,j) - tl%d(k,i,j) )*sf(k)
                    rtt%d(k,i,j) = rtt%d(k,i,j) - ( rt%d(kp1,i,j) - rt%d(k,i,j) )*sf(k)
                    !
                    ! Radiative cooling: 2.5 K/day
                    tt%d(k,i,j) = tt%d(k,i,j)  - 2.5/86400.
                END DO
            END DO
        END DO
        !
    CASE ('amazon')
        ! Amazon
        ! --------
        ! - to be added -
        !
    CASE DEFAULT
       IF (myid == 0) PRINT *, '  ABORTING: inproper CALL to radiation'
       CALL appl_abort(0)
    END SELECT

  END SUBROUTINE case_forcing
  !
  ! -------------------------------------------------------------------
  ! Subroutine bellon_rad:  call simple radiative parameterization
  !
  SUBROUTINE bellon(n1,n2,n3,flx,sflx,zt,dzt,dzm,tt,tl,rtt,rt, ut,u,vt,v)

    INTEGER, INTENT (in) :: n1,n2, n3

    TYPE(FloatArray1d), INTENT (in)              :: zt, dzt, dzm
    TYPE(FloatArray3d), INTENT(inout)            :: tt, tl, rtt, rt  
    REAL, DIMENSION (n1, n2, n3), INTENT (inout) :: ut,u,vt,v
    TYPE(FloatArray3d), INTENT (inout)             :: flx, sflx
    REAL, PARAMETER      :: w0 = 7.5e-3, H = 1000., Qrate = 2.5/86400.

    INTEGER :: i,j,k,kp1
    REAL    :: grad,wk

    DO j = 3, n3-2
       DO i = 3, n2-2
          !
          ! subsidence
          !
          flx%d(1,i,j)  = 0.
          sflx%d(1,i,j) = 0.
          DO k = 2, n1-2
             kp1 = k+1
             wk = w0*(1.-exp(-zt%d(k)/H))
             grad = Qrate/wk
             flx%d(k,i,j)  = wk*((tl%d(kp1,i,j)-tl%d(k,i,j))*dzt%d(k)-grad)
             sflx%d(k,i,j) = wk*((rt%d(kp1,i,j)-rt%d(k,i,j))*dzt%d(k)-grad)
             tt%d(k,i,j) = tt%d(k,i,j) + flx%d(k,i,j)
             rtt%d(k,i,j)= rtt%d(k,i,j) + &
                         wk*(rt%d(kp1,i,j)-rt%d(k,i,j))*dzt%d(k)
             ut(k,i,j) =  ut(k,i,j) + &
                          wk*(u(kp1,i,j)-u(k,i,j))*dzm%d(k)
             vt(k,i,j) =  vt(k,i,j) + &
                          wk*(v(kp1,i,j)-v(k,i,j))*dzm%d(k)
          END DO
          flx%d(n1,  i,j)  = 0.
          flx%d(n1-1,i,j)  = 0.
          sflx%d(n1,  i,j) = 0.
          sflx%d(n1-1,i,j) = 0.
       END DO
    END DO

  END SUBROUTINE bellon

  ! POISTA:: ::
  ! THIS WAS IN STAT.F90
   REAL FUNCTION get_zi (n1, n2, n3, itype, sx, xx, z, threshold)

      INTEGER, INTENT (in) :: n1, n2, n3, itype
      TYPE(FloatArray1d), INTENT (in)    :: xx, z
      REAL, INTENT(in) :: threshold
      TYPE(FloatArray3d), INTENT(in)     :: sx
      
      INTEGER :: i, j, k, kk
      REAL    :: zibar, sval, dmy, scr(n2,n3)

      get_zi = -999.
      SELECT CASE(itype)
         CASE (1)
            !
            ! find level at which sx=threshold (xx is one over grid spacing)
            !
            zibar = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  k = 2
                  DO WHILE (k < n1-2 .AND. sx%d(k,i,j) > threshold)
                     k = k+1
                  END DO
                  IF (k == n1-2) zibar = -999.
                  IF (zibar /= -999.) zibar = zibar + z%d(k-1) +  &
                                             (threshold - sx%d(k-1,i,j))/xx%d(k-1)     /  &
                                             (sx%d(k,i,j) - sx%d(k-1,i,j) + epsilon(1.))
               END DO
            END DO
            IF (zibar /= -999.) get_zi = zibar/REAL((n3-4)*(n2-4))

         CASE(2)
            !
            ! find level of maximum gradient (xx is one over grid spacing)
            !
            scr = 0.
            DO j = 3, n3-2
               DO i = 3, n2-2
                  sval = 0.
                  DO k = 2, n1-5
                     dmy = (sx%d(k+1,i,j)-sx%d(k,i,j))*xx%d(k)
                     IF (dmy > sval) THEN
                        sval = dmy
                        scr(i,j) = z%d(k)
                     END IF
                  END DO
               END DO
            END DO
            get_zi = get_avg2dh(n2,n3,scr)

         CASE(3)
            !
            ! find level where xx is a maximum
            !
            sval = -huge(1.)
            kk = 1
            DO k = 2, n1
               IF (xx%d(k) > sval) THEN
                  kk = k
                  sval = xx%d(k)
               END IF
            END DO
            get_zi = z%d(kk)

         CASE(4)
            !
            ! find level where xx is a minimum
            !
            sval = huge(1.)
            kk = 1
            DO k = 2, n1-2
               IF (xx%d(k) < sval) THEN
                  kk = k
                  sval = xx%d(k)
               END IF
            END DO
            get_zi = z%d(kk)
      END SELECT

   END FUNCTION get_zi


  
END MODULE forc
