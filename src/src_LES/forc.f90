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
module forc

  use defs, only      : cp
  use radiation, only : d4stream
  use stat, only : sflg
  implicit none

  ! these are now all namelist parameters
  character (len=10) :: case_name = 'none'               
  character (len=50) :: radsounding = 'datafiles/dsrt.lay'  ! Juha: Added so the radiation background sounding can be given
                                                            ! from the NAMELIST
  REAL :: sfc_albedo = 0.05
  REAL :: div = 0., zmaxdiv = 1e6
  REAL :: xka = 119., fr0 = 96.2, fr1 = 61.2, alpha = 1.0
  REAL :: rc_limit = 0.01e-3, rt_limit = 8e-3
  LOGICAL :: useMcICA = .TRUE.
  LOGICAL :: RadConstPress = .FALSE. ! Keep constant pressure levels
  INTEGER :: RadPrecipBins = 0 ! Add precipitation bins to cloud water (for level 3 and up)
  INTEGER :: RadSnowBins = 0 ! Add snow bins to cloud ice (for level 5 and up)

contains
  !
  SUBROUTINE surface_naerot(fluksi)
    use grid, only : nxp, nyp, a_naerot, a_maerot
    USE mo_submctl, ONLY : nbins, aerobins, in2a, fn2a, in2b, pi6, iss, ih2o, rhowa, rhoss
    IMPLICIT NONE
    REAL :: fluksi(nxp,nyp,fn2a) ! Rate of change in number concentration (#/kg/s)
    INTEGER :: i, j
    REAL :: mdry, mwat
    !
    IF (iss<1) STOP 'Sea salt not included!'
    !
    ! Ignore 1a, because there is no sea salt
    DO i=in2a,fn2a
        ! Dry particle mass (assuming SS): use bin GMD
        !   Note: bin center is volume mean - using other than that will cause problems!
        mdry=4.*pi6*(aerobins(i)**3+aerobins(i+1)**3)*rhoss
        ! Water mass: assume volume growth factor of 10
        mwat=10.*mdry*rhowa/rhoss
        !
        ! Apply to 2b bins
        j = in2b + i - in2a
        a_naerot(2,:,:,j) = a_naerot(2,:,:,j) + fluksi(:,:,i)
        ! ... and specifically to SS
        j = (iss-1)*nbins + in2b + i - in2a
        a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + fluksi(:,:,i)*mdry
        !  ... and just add water
        j = (ih2o-1)*nbins + in2b + i - in2a
        a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + fluksi(:,:,i)*mwat
    ENDDO
  END SUBROUTINE surface_naerot

  !
  ! -------------------------------------------------------------------
  ! subroutine forcings:  calls the appropriate large-scale forcings
  !
  subroutine forcings(time_in, cntlat, sst)

    use grid, only: nxp, nyp, nzp, zm, zt, dzt, dzm, dn0, iradtyp, a_rc     &
         , a_rflx, a_sflx, albedo, a_tt, a_tp, a_rt, a_rp, a_pexnr, a_temp  &
         , a_rv, a_rpp, a_npp, CCN, pi0, pi1, level, a_maerop, &
         a_ncloudp, a_mcloudp, a_nprecpp, a_mprecpp, a_nicep, a_micep, a_nsnowp, a_msnowp, &
         a_fus, a_fds, a_fuir, a_fdir
    USE mo_submctl, ONLY : nbins,ncld,nice,nprc,nsnw
    use mpi_interface, only : myid, appl_abort

    real, intent (in) :: time_in, cntlat, sst
    REAL :: znc(nzp,nxp,nyp), zrc(nzp,nxp,nyp), zni(nzp,nxp,nyp), zri(nzp,nxp,nyp)

    select case(iradtyp)
    case (0)
       ! No radiation or large-scale forcing

    case (1)
       ! No radiation, just case-dependent large-scale forcing
       !
       IF ( case_name /= 'none' ) THEN
          call case_forcing(nzp,nxp,nyp,zt,dzt,dzm,div,a_tp,a_rp,a_tt,a_rt)
       END IF

    case (2)
       ! Parameterized radiation (GCSS) and subsidence (constant div up to z=zmaxdiv)
       !
       IF (level <= 3) THEN
          zrc(:,:,:) = a_rc(:,:,:) + a_rpp(:,:,:) ! Liquid water mixing ratio - radiative effects
          znc(:,:,:) = a_rp(:,:,:) + a_rpp(:,:,:) ! Total water mixing ratio - for determining inversion height
       ELSE
          zrc(:,:,:) = SUM(a_maerop(:,:,:,1:nbins),DIM=4) + &
                       SUM(a_mcloudp(:,:,:,1:ncld),DIM=4) + &
                       SUM(a_mprecpp(:,:,:,1:nprc),DIM=4) ! Aerosol, cloud and rain water
          znc(:,:,:) = a_rp(:,:,:) + zrc(:,:,:) ! Water vapor and liquid water, but no ice or snow
       ENDIF
       call new_gcss_rad(nzp, nxp, nyp, zrc, znc, a_rflx)

    case (3)
       ! Fu and Liou (1993) radiation code and case-dependent large-scale forcing
       !
       IF (level <= 3) THEN
          znc(:,:,:) = CCN
          zrc(:,:,:) = a_rc(:,:,:) ! Cloud water
          IF (level == 3 .AND. RadPrecipBins > 0) THEN
             ! Include clouds and rain - number is a mass-mean for cloud and rain species
             zrc(:,:,:) = a_rc(:,:,:) + a_rpp(:,:,:)
             WHERE (zrc>1e-10) znc = (max(0.,a_rpp*a_npp)+max(0.,a_rc*CCN))/zrc
          ENDIF
          call d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
               dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rv, zrc, znc, a_tt,  &
               a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
               useMcICA=useMcICA, ConstPrs=RadConstPress)

       ELSE IF (level == 4) THEN
          ! Water is the first SALSA species
          !zrc(:,:,:) = a_rc(:,:,:) ! Cloud and aerosol water
          zrc(:,:,:) = SUM(a_mcloudp(:,:,:,1:ncld),DIM=4) ! Cloud droplets
          znc(:,:,:) = SUM(a_ncloudp(:,:,:,:),DIM=4)
          IF (RadPrecipBins>0) THEN ! Add precipitation bins
             zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,1:min(RadPrecipBins,nprc)),DIM=4)
             znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,1:min(RadPrecipBins,nprc)),DIM=4)
          ENDIF
          CALL d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
               dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
               a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, radsounding=radsounding, &
               useMcICA=useMcICA, ConstPrs=RadConstPress)

       ELSE IF (level == 5) THEN
          ! Water is the first SALSA species
          !zrc(:,:,:) = a_rc(:,:,:) ! Cloud and aerosol water
          zrc(:,:,:) = SUM(a_mcloudp(:,:,:,1:ncld),DIM=4) ! Cloud droplets
          znc(:,:,:) = SUM(a_ncloudp(:,:,:,:),DIM=4)
          IF (RadPrecipBins>0) THEN ! Add precipitation bins
             zrc(:,:,:) = zrc(:,:,:) + SUM(a_mprecpp(:,:,:,1:min(RadPrecipBins,nprc)),DIM=4)
             znc(:,:,:) = znc(:,:,:) + SUM(a_nprecpp(:,:,:,1:min(RadPrecipBins,nprc)),DIM=4)
          ENDIF
          zri(:,:,:) = SUM(a_micep(:,:,:,1:nice),DIM=4) ! Ice
          zni(:,:,:) = SUM(a_nicep(:,:,:,:),DIM=4)
          IF (RadSnowBins>0) THEN ! Add snow bins
             zri(:,:,:) = zri(:,:,:) + SUM(a_msnowp(:,:,:,1:min(RadSnowBins,nsnw)),DIM=4)
             zni(:,:,:) = zni(:,:,:) + SUM(a_nsnowp(:,:,:,1:min(RadSnowBins,nsnw)),DIM=4)
          ENDIF
          CALL d4stream(nzp, nxp, nyp, cntlat, time_in, sst, sfc_albedo, &
               dn0, pi0, pi1, dzt, a_pexnr, a_temp, a_rp, zrc, znc, a_tt,  &
               a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, albedo, ice=zri,nice=zni,radsounding=radsounding, &
               useMcICA=useMcICA, ConstPrs=RadConstPress)

       END IF

       ! Case-dependent large-scale forcing
       IF ( case_name /= 'none') THEN
          CALL case_forcing(nzp,nxp,nyp,zt,dzt,dzm,div,a_tp,a_rp,a_tt,a_rt)
       END IF

    case default
       if (myid == 0) print *, '  ABORTING: improper call to forcing'
       call appl_abort(0)

    end select

  end subroutine forcings

  !
  ! -------------------------------------------------------------------
  ! subroutine gcss_rad:  call simple radiative parameterization and 
  ! simultaneously update fields due to vertical motion as given by div
  !
  subroutine gcss_rad(n1,n2,n3,xka,fr0,fr1,div,rc,dn0,flx,zt,zm,dzt,   &
       tt,tl,rtt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: xka, fr0, fr1, div
    real, intent (in)   :: zt(n1),zm(n1),dzt(n1),dn0(n1),rc(n1,n2,n3),   &
         tl(n1,n2,n3),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3),rtt(n1,n2,n3)
    real, intent (out)  :: flx(n1,n2,n3)

    integer :: i, j, k, km1, kp1, ki
    real    :: lwp(n2,n3), fact

    lwp=0.
    do j=3,n3-2
       do i=3,n2-2
          ki = n1
          do k=1,n1
             km1=max(1,k-1)
             lwp(i,j)=lwp(i,j)+max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
             flx(k,i,j)=fr1*exp(-1.*xka*lwp(i,j))
             if ( (rc(k,i,j) > 0.01e-3) .and. (rt(k,i,j) >= 0.008) ) ki=k
          enddo

          fact = div*cp*dn0(ki)
          do k=2,n1
             km1=max(2,k-1)
             lwp(i,j)=lwp(i,j)-max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=flx(k,i,j)+fr0*exp(-1.*xka*lwp(i,j))
             if (zm(k) > zm(ki) .and. ki > 1 .and. fact > 0.) then
                flx(k,i,j)=flx(k,i,j) + fact*(0.25*(zm(k)-zm(ki))**1.333 + &
                  zm(ki)*(zm(k)-zm(ki))**0.333333)
             end if
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzt(k)/(dn0(k)*cp)
          enddo
          !
          ! subsidence
          !
          if (div /= 0.) then
             do k=2,n1-2
                kp1 = k+1
                tt(k,i,j) = tt(k,i,j) + div*zt(k)*(tl(kp1,i,j)-tl(k,i,j))*dzt(k)
                rtt(k,i,j)=rtt(k,i,j) + div*zt(k)*(rt(kp1,i,j)-rt(k,i,j))*dzt(k)
             end do
          end if
       enddo
    enddo

  end subroutine gcss_rad


  subroutine new_gcss_rad(n1,n2,n3,rc,rt,flx)
    USE grid, ONLY : a_ncloudp, a_nprecpp, a_mprecpp, a_nicep, a_nsnowp, a_msnowp, &
         a_naerop, a_naerot, a_ncloudt, a_nicet, a_nsnowt, a_maerop, a_mcloudp, a_micep,  &
         a_maerot, a_mcloudt, a_micet, a_msnowt, a_nprecpt, a_mprecpt, level, a_temp, a_theta, &
         zt, zm, dzt, dzm, dn0, a_tt, a_tp, a_rt, a_rp
    implicit none
    integer, intent (in)::  n1,n2, n3
    real, intent (in)   ::  rc(n1,n2,n3),rt(n1,n2,n3)
    real, intent (out)  ::  flx(n1,n2,n3)

    integer :: i, j, k, kp1, ki
    real    :: lwp, fact
    real, dimension (n1) :: sf

    ! a) Radiation
    flx=0.
    do j=3,n3-2
       do i=3,n2-2
          ki=n1
          lwp=0.
          ! No cloud water at level k=1
          flx(1,i,j)=fr1*exp(-1.*xka*lwp)
          do k=2,n1
             lwp=lwp+max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=fr1*exp(-1.*xka*lwp)
             if ( rc(k,i,j) >= rc_limit .and. rt(k,i,j) >= rt_limit) ki=k
          enddo
          !
          fact=dn0(ki)*cp*div*alpha
          ! Level k=1
          flx(1,i,j)=flx(1,i,j)+fr0*exp(-1.*xka*lwp)
          do k=2,n1
             lwp=lwp-max(0.,rc(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=flx(k,i,j)+fr0*exp(-1.*xka*lwp)
             if (k > ki .and. fact > 0.) then
                flx(k,i,j)=flx(k,i,j) + fact*(0.25*(zm(k)-zm(ki))**1.333 + &
                    zm(ki)*(zm(k)-zm(ki))**0.333333)
             end if
             ! dtheta_il/dT=dtheta/dT=1/pi=theta/T => dtheta_il = theta/T*dT
             a_tt(k,i,j)=a_tt(k,i,j)-(flx(k,i,j)-flx(k-1,i,j))*dzt(k)/(dn0(k)*cp)*a_theta(k,i,j)/a_temp(k,i,j)
          enddo
      enddo
    enddo
    !
    ! b) Subsidence
    IF (ABS(div)<1e-10) RETURN
    !
    sf=0.
    do k=2,n1-2
        ! calculate subsidence factor (wsub / dz)
        sf(k) = -div*min( zmaxdiv,zt(k) )*dzt(k)
    end do
    !
    DO j=3,n3-2
        DO i=3,n2-2
            DO k=2,n1-1
                kp1 = k+1
                a_tt(k,i,j) = a_tt(k,i,j) - ( a_tp(kp1,i,j) - a_tp(k,i,j) )*sf(k)
                a_rt(k,i,j) = a_rt(k,i,j) - ( a_rp(kp1,i,j) - a_rp(k,i,j) )*sf(k)

                IF (level>=4) THEN
                  a_maerot(k,i,j,:) = a_maerot(k,i,j,:) - ( a_maerop(kp1,i,j,:) - a_maerop(k,i,j,:) )*sf(k)
                  a_mcloudt(k,i,j,:) = a_mcloudt(k,i,j,:) - ( a_mcloudp(kp1,i,j,:) - a_mcloudp(k,i,j,:) )*sf(k)
                  a_mprecpt(k,i,j,:) = a_mprecpt(k,i,j,:) - ( a_mprecpp(kp1,i,j,:) - a_mprecpp(k,i,j,:) )*sf(k)
                  a_naerot(k,i,j,:) = a_naerot(k,i,j,:) - ( a_naerop(kp1,i,j,:) - a_naerop(k,i,j,:) )*sf(k)
                  a_ncloudt(k,i,j,:) = a_ncloudt(k,i,j,:) - ( a_ncloudp(kp1,i,j,:) - a_ncloudp(k,i,j,:) )*sf(k)
                  a_nprecpt(k,i,j,:) = a_nprecpt(k,i,j,:) - ( a_nprecpp(kp1,i,j,:) - a_nprecpp(k,i,j,:) )*sf(k)
                ENDIF
                IF (level>=5) THEN
                  a_micet(k,i,j,:) = a_micet(k,i,j,:) - ( a_micep(kp1,i,j,:) - a_micep(k,i,j,:) )*sf(k)
                  a_msnowt(k,i,j,:) = a_msnowt(k,i,j,:) - ( a_msnowp(kp1,i,j,:) - a_msnowp(k,i,j,:) )*sf(k)
                  a_nicet(k,i,j,:) = a_nicet(k,i,j,:) - ( a_nicep(kp1,i,j,:) - a_nicep(k,i,j,:) )*sf(k)
                  a_nsnowt(k,i,j,:) = a_nsnowt(k,i,j,:) - ( a_nsnowp(kp1,i,j,:) - a_nsnowp(k,i,j,:) )*sf(k)
                ENDIF
            END DO
        END DO
    END DO
  end subroutine new_gcss_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine smoke_rad:  call simple radiative parameterization for 
  ! the smoke cloud
  !
  subroutine smoke_rad(n1,n2,n3,dn0,flx,zm,dzt,tt,rt)

    integer, intent (in):: n1,n2, n3
    real, intent (in)   :: zm(n1),dzt(n1),dn0(n1),rt(n1,n2,n3)
    real, intent (inout):: tt(n1,n2,n3)
    real, intent (out)  :: flx(n1,n2,n3)
    real, parameter     :: xka= 50.0, fr0=60.0

    integer :: i,j,k, km1, ki
    real    :: smoke(n2,n3)

    smoke=0.
    do j=3,n3-2
       do i=3,n2-2
          ki = n1
          do k=1,n1
             km1=max(1,k-1)
             smoke(i,j)=smoke(i,j)+max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(km1)))
          enddo

          do k=2,n1
             km1=max(2,k-1)
             smoke(i,j)=smoke(i,j)-max(0.,rt(k,i,j)*dn0(k)*(zm(k)-zm(k-1)))
             flx(k,i,j)=fr0*exp(-1.*xka*smoke(i,j))
             tt(k,i,j) =tt(k,i,j)-(flx(k,i,j)-flx(km1,i,j))*dzt(k)/(dn0(k)*cp)
          enddo
       enddo
    enddo

  end subroutine smoke_rad
  !
  ! -------------------------------------------------------------------
  ! subroutine case_forcing: adjusts tendencies according to a specified
  ! large scale forcing.  Normally case (run) specific.
  !
  subroutine case_forcing(n1,n2,n3,zt,dzt,dzm,zdiv,tl,rt,tt,rtt)

    use mpi_interface, only : pecount, double_scalar_par_sum,myid, appl_abort
    use stat, only : get_zi
    USE grid, ONLY : a_ncloudp, a_nprecpp, a_mprecpp, a_nicep, a_nsnowp, a_msnowp, &
         a_naerop, a_naerot, a_ncloudt, a_nicet, a_nsnowt, a_maerop, a_mcloudp, a_micep,  &
         a_maerot, a_mcloudt, a_micet, a_msnowt, a_nprecpt, a_mprecpt, level

    integer, intent (in):: n1,n2, n3
    real, dimension (n1), intent (in)          :: zt, dzt, dzm
    real, intent(in)                           :: zdiv
    real, dimension (n1,n2,n3), intent (in)    :: tl, rt
    real, dimension (n1,n2,n3), intent (inout) :: tt, rtt

    integer :: i,j,k,kp1
    real, dimension (n1) :: sf
    real, parameter :: zmx_sub = 2260. ! originally 2260.

    real (kind=8) :: zig, zil
    real          :: zibar

    zig = 0.0; zil = 0.0; zibar = 0.0
    kp1= 0
    select case (trim(case_name))
    case('default')
       !
       ! User specified divergence used as a simple large scle forcing for moisture and temperature fields
       ! -------------------------------------------------------------------------------------------------
       !
       DO k=2,n1-1
          sf(k)=-zdiv*MIN(zt(k),zmaxdiv)*dzt(k)
       END DO
       DO j=3,n3-2
          DO i=3,n2-2
             DO k=2,n1-1
                kp1 = k+1
                tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)

                IF (level>=4) THEN
                  a_maerot(k,i,j,:) = a_maerot(k,i,j,:) - ( a_maerop(kp1,i,j,:) - a_maerop(k,i,j,:) )*sf(k)
                  a_mcloudt(k,i,j,:) = a_mcloudt(k,i,j,:) - ( a_mcloudp(kp1,i,j,:) - a_mcloudp(k,i,j,:) )*sf(k)
                  a_mprecpt(k,i,j,:) = a_mprecpt(k,i,j,:) - ( a_mprecpp(kp1,i,j,:) - a_mprecpp(k,i,j,:) )*sf(k)
                  a_naerot(k,i,j,:) = a_naerot(k,i,j,:) - ( a_naerop(kp1,i,j,:) - a_naerop(k,i,j,:) )*sf(k)
                  a_ncloudt(k,i,j,:) = a_ncloudt(k,i,j,:) - ( a_ncloudp(kp1,i,j,:) - a_ncloudp(k,i,j,:) )*sf(k)
                  a_nprecpt(k,i,j,:) = a_nprecpt(k,i,j,:) - ( a_nprecpp(kp1,i,j,:) - a_nprecpp(k,i,j,:) )*sf(k)
                ENDIF
                IF (level>=5) THEN
                  a_micet(k,i,j,:) = a_micet(k,i,j,:) - ( a_micep(kp1,i,j,:) - a_micep(k,i,j,:) )*sf(k)
                  a_msnowt(k,i,j,:) = a_msnowt(k,i,j,:) - ( a_msnowp(kp1,i,j,:) - a_msnowp(k,i,j,:) )*sf(k)
                  a_nicet(k,i,j,:) = a_nicet(k,i,j,:) - ( a_nicep(kp1,i,j,:) - a_nicep(k,i,j,:) )*sf(k)
                  a_nsnowt(k,i,j,:) = a_nsnowt(k,i,j,:) - ( a_nsnowp(kp1,i,j,:) - a_nsnowp(k,i,j,:) )*sf(k)
                ENDIF
             END DO
          END DO
       END DO

    case('rico')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < zmx_sub) then
             sf(k) =  -0.005*zt(k)/zmx_sub
          else
             sf(k) =  -0.005 
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! subsidence
                ! 
                kp1 = k+1
                tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                !
                ! temperature advection and radiative cooling
                !
                tt(k,i,j) = tt(k,i,j)  - 2.5/86400.
                !
                ! moisture advection
                !
                if (zt(k) <= 2980.) then
                   rtt(k,i,j) = rtt(k,i,j)  - (1. -  1.3456*zt(k)/2980.)/8.64e7
                else
                   rtt(k,i,j) = rtt(k,i,j)  + .3456/8.64e7
                end if
             enddo
          enddo
       enddo

    case ('bomex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       do k=2,n1-2
          if (zt(k) < 1500.) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065  + 0.0065*(zt(k)-1500.)/600.)
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < 1500.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5
                else if (zt(k) < 2000.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-1500.)*1.e-3)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < 300.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8
                elseif (zt(k) < 500.) then
                   rtt(k,i,j) = rtt(k,i,j)  - 1.2e-8*(1.- (zt(k)-300.)/200.)
                end if
             enddo
          enddo
       enddo
    case ('atex')
       !
       ! calculate subsidence factor (wsub / dz)
       !
       zil = get_zi (n1, n2, n3, 2, rt, dzm, zt, 6.5e-3)
       call double_scalar_par_sum(zil,zig)
       zibar = real(zig/pecount)

       do k=2,n1-2
          if (zt(k) < zibar) then
             sf(k) =  -0.0065*zt(k)/1500.
          else
             sf(k) =  min(0.,-0.0065*(1 - (zt(k)-zibar)/300.))
          end if
          sf(k) = sf(k)*dzt(k)
       end do

       do j=3,n3-2
          do i=3,n2-2
             do k=2,n1-2
                !
                ! temperature advection and radiative cooling
                !
                kp1 = k+1
                if (zt(k) < zibar) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1. + (1.- zt(k)/zibar)/2.)
                else if (zt(k) < zibar+300.) then
                   tt(k,i,j) = tt(k,i,j) - ( tl(kp1,i,j)-tl(k,i,j) )*sf(k) &
                        - 2.315e-5*(1.- (zt(k)-zibar)/300.)
                end if
                !
                ! moisture advection
                !
                rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                if (zt(k) < zibar) rtt(k,i,j) = rtt(k,i,j)  - 1.5e-8
             enddo
          enddo
       enddo
        !
    case ('ascos')
        ! ASCOS
        ! ---------
        !
        do k=2,n1-2
            ! calculate subsidence factor (wsub / dz)
            sf(k) = -5.0e-6*min(2000.0,zt(k))*dzt(k)
        end do
        !
        do j=3,n3-2
            do i=3,n2-2
                do k=2,n1-2
                    !
                    ! Temperature and humidity advection due to subsidence
                    !
                    kp1 = k+1
                    tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                    rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                enddo
            enddo
        enddo
        !
    case ('barba')
        ! Barbados
        ! -----------
        ! Large scale subsidence: w(z)=w0*(1-exp(z/H)), where w0=7.5 mm/s and H=1000 m.
        ! Radiative cooling rate: 2.5 K/day
        ! No temperature or humidity advection
        !
        ! calculate subsidence factor (wsub / dz)
        do k=2,n1-2
            sf(k) = -7.5e-3*(1.0-exp(-zt(k)/1000.0))*dzt(k)
        end do
        !
        do j=3,n3-2
            do i=3,n2-2
                do k=2,n1-2
                    ! Subsidence
                    kp1 = k+1
                    tt(k,i,j)  =  tt(k,i,j) - ( tl(kp1,i,j) - tl(k,i,j) )*sf(k)
                    rtt(k,i,j) = rtt(k,i,j) - ( rt(kp1,i,j) - rt(k,i,j) )*sf(k)
                    !
                    ! Radiative cooling: 2.5 K/day
                    tt(k,i,j) = tt(k,i,j)  - 2.5/86400.
                enddo
            enddo
        enddo
        !
    CASE ('amazon')
        ! Amazon
        ! --------
        ! - to be added -
        !
    case default
       if (myid == 0) print *, '  ABORTING: inproper call to radiation'
       call appl_abort(0)
    end select

  end subroutine case_forcing
  !
  ! -------------------------------------------------------------------
  ! subroutine bellon_rad:  call simple radiative parameterization
  !
  subroutine bellon(n1,n2,n3,flx,sflx,zt,dzt,dzm,tt,tl,rtt,rt, ut,u,vt,v)

    integer, intent (in) :: n1,n2, n3

    real, dimension (n1), intent (in)            :: zt, dzt, dzm
    real, dimension (n1, n2, n3), intent (inout) :: tt, tl, rtt, rt, ut,u,vt,v
    real,  dimension (n1, n2, n3), intent (out)  :: flx, sflx
    real, parameter      :: w0= 7.5e-3, H=1000., Qrate = 2.5/86400.

    integer :: i,j,k,kp1
    real    :: grad,wk

    do j=3,n3-2
       do i=3,n2-2
          !
          ! subsidence
          !
          flx(1,i,j)  = 0.
          sflx(1,i,j) = 0.
          do k=2,n1-2
             kp1 = k+1
             wk = w0*(1.-exp(-zt(k)/H))
             grad = Qrate/wk
             flx(k,i,j)  = wk*((tl(kp1,i,j)-tl(k,i,j))*dzt(k)-grad)
             sflx(k,i,j) = wk*((rt(kp1,i,j)-rt(k,i,j))*dzt(k)-grad)
             tt(k,i,j) = tt(k,i,j) + flx(k,i,j)
             rtt(k,i,j)=rtt(k,i,j) + &
                  wk*(rt(kp1,i,j)-rt(k,i,j))*dzt(k)
             ut(k,i,j) =  ut(k,i,j) + &
                  wk*(u(kp1,i,j)-u(k,i,j))*dzm(k)
             vt(k,i,j) =  vt(k,i,j) + &
                  wk*(v(kp1,i,j)-v(k,i,j))*dzm(k)
          end do
          flx(n1,  i,j)  = 0.
          flx(n1-1,i,j)  = 0.
          sflx(n1,  i,j) = 0.
          sflx(n1-1,i,j) = 0.
       enddo
    enddo

  end subroutine bellon
 
end module forc
