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
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module srfc

  integer :: isfctyp = 0
  real    :: zrough =  0.001
  real    :: ubmin  =  0.20
  real    :: dthcon = 100.0
  real    :: drtcon = 0.0



! Sami added ----->
! Initial values for surface properties
 ! real :: W1 = 0.9   !Water content      ... Definition in grid now because of restrat files
 ! real :: W2 = 0.9
 ! real :: W3 = 0.9

  real ::  B1 = 6.5
  real ::  B3 = 7.26
  real ::  K_s1 = 5.8e-5
  real ::  K_s3 = 0.89e-5
  real ::  fii_s1 = -0.036
  real ::  fii_s3 = -0.085
  real ::  thetaS1 = 0.98 ! Soil porosity
  real ::  thetaS2 = 0.98
  real ::  thetaS3 = 0.488
  real ::  D1 = 0.1  !Depth of different layers !
  real ::  D2 = 0.3
  real ::  D3 = 0.6

! <--- Sami added


  ! Chlorophyll alpha concentration for marine organic emissions
  real :: wtrChlA = -1.   ! kg/m3
  logical :: ifPOCadd = .false.
  ! Isoprene and monoterpene concentrations in the ocean surface layer
  real :: wtrIsop = -1.   ! mol/m3
  real :: wtrMtrp = -1.   ! mol/m3

contains

  ! --------------------------------------------------------------------------
  ! Size-resolved marine aerosol production rates as a function of sea surface temperature (SST, K)
  ! and friction velocity (u_star, m/s). Parameterization from:
  !  Mortensson et al., Laboratory simulations and parameterization of the primary marine
  !  aerosol production, J. Geophys. Res., 108, 4297, doi:10.1029/2002JD002263, 2003
  ! The result is the rate of change in particle number concentration (#/kg/s), or tendency, for
  ! each size bin at the first level above sea surface. This is then used to compute tendencies
  ! for mass concentrations.
  !
  SUBROUTINE marine_aero_flux(sst,ovf)
    ! Parameters
    use defs, only: vonk, g
    use grid, only: nxp, nyp, a_ustar, a_dn, zm, nbins, a_naerot, a_maerot
    use mo_submctl, only: rhooc, rhoss, nvbs_setup, aerobins, in2a, fn2a, &
                    in2b, pi6, iss, ioc, ih2o, rhowa, rhoss, rhooc, nvbs_setup
    use mo_vbsctl, only: vbs_set
    use mo_vbs, only: spec_density
    use stat, only: flux_stat, sflg

    IMPLICIT NONE
    REAL, INTENT(IN) :: sst     ! Sea surface temperature (K)
    REAL, INTENT(INOUT) :: ovf(nxp,nyp,fn2a) ! volume fraction of organic matter
    ! Local
    INTEGER :: i, j, k
    REAL :: usum, zs, dia, Ak, Bk, rhorho, dia80, omf, vdry
    REAL :: w(nxp,nyp), u10(nxp,nyp),omf_max(nxp,nyp) ! whitecap cover, 10 m wind speed
    REAL :: flx(fn2a+1) ! Production rate for each size bin
    REAL :: dcdt(nxp,nyp,fn2a) ! Particle concentration tendency (#/kg/s)

    IF (iss<1 .and. ioc<1 .and. nvbs_setup<1) STOP 'No sea spray species included!'

    ! Calculate particle flux for each size bin limit
    DO k=1,fn2a+1
        dia=aerobins(k)*2.
        ! Parameters for Eq. 5 from Table 1
        IF (dia<0.020e-6) THEN
            Ak=0.
            Bk=0.
        ELSEIF (dia<0.145e-6) THEN
            Ak=-2.576e35*dia**4+5.932e28*dia**3-2.867e21*dia**2-3.003e13*dia-2.881e6
            Bk= 7.188e37*dia**4-1.616e31*dia**3+6.791e23*dia**2+1.829e16*dia+7.609e8
        ELSEIF (dia<0.419e-6) THEN
            Ak=-2.452e33*dia**4+2.404e27*dia**3-8.148e20*dia**2+1.183e14*dia-6.743e6
            Bk= 7.368e35*dia**4-7.310e29*dia**3+2.528e23*dia**2-3.787e16*dia+2.279e9
        ELSEIF (dia<2.8e-6) THEN
            Ak= 1.085e29*dia**4-9.841e23*dia**3+3.132e18*dia**2-4.165e12*dia+2.181e6
            Bk=-2.859e31*dia**4+2.601e26*dia**3-8.297e20*dia**2+1.105e15*dia-5.800e8
        ELSE
            Ak=0.
            Bk=0.
            ! Extrapolation: Ak*sst+Bk vs D is almost linear in the log-log scale
            Ak=139934 ! Ak(2.8e-6)
            Bk=-3.84343e7 ! Bk(2.8e-6)
            ! Cheat a bit by using variable Bk
            Bk=(Ak*sst+Bk)*(dia/2.8e-6)**(-3.5)
            Ak=0.
        ENDIF
        ! Equation 6: particle flux per whitecap area, dFp/dlog(Dp) [#/m^2/s]
        flx(k)=MAX(0.,Ak*sst+Bk)
    ENDDO

    ! Roughness height is needed for the 10 m wind speeds
    zs = zrough
    IF (zrough <= 0.) THEN ! Calculate
        usum = 0.
        DO j=3,nyp-2
            DO i=3,nxp-2
                usum = usum + a_ustar(i,j)
            END DO
        ENDDO
        usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
        zs = max(0.0001,(0.016/g)*usum**2)
    ENDIF

    ! Whitecap cover (% => fraction) based on 10 m wind speeds (eq. 2)
    u10(:,:) = a_ustar(:,:)/vonk*log(10.0/zs)
    w(:,:)=0.01*3.84e-4*u10**3.41
    ! Organic mass fraction grom Gantt ea (2011)
    ! conversion to mg from kg
    if(wtrChlA>0.)then
        if(nvbs_setup > 0)then
            rhorho = spec_density(vbs_set(1)%spid)/rhoss
        else
            rhorho = rhooc/rhoss ! ratio of densities of marine organic matter and salt
        endif
        omf_max = 1./(1.+exp(-2.63e6*wtrChlA+0.18*U10))
    endif

    ! Particle concentration tendency for each size bin
    dcdt(:,:,:)=0.
    DO k=in2a,fn2a ! Ignore 1a, because there is no sea salt
        ! Mean flux: 0.5*(flx(k)+flx(k+1))
        ! From dFp/dlog(Dp) to dFp: multiply by log10(rad(k+1)/rad(k))
        ! Multiply by the whitecap cover w
        ! Convert #/m^2/s to the rate of change in concetration (#/kg/s): multily by 1/rho/dz
        DO j=3,nyp-2
          DO i=3,nxp-2
            dcdt(i,j,k)=0.5*(flx(k)+flx(k+1))*log10(aerobins(k+1)/aerobins(k))*w(i,j)/a_dn(2,i,j)/(zm(3)-zm(2))

            ! Size dependent organic mass fraction from Gantt ea, 2011 (eq. 3)
            if(wtrChlA>0.)then
                dia = (4.*(aerobins(k)**3+aerobins(k+1)**3))**(1./3.)*1.e6 ! bin volume-mean diameter in micrometers
                ! Organic fraction is parameterized for particle diameter at 80% RH
                ! which depends on composition, so technically would need to iterate.
                ! However, hopefully it converges after limited nr of timesteps.
                ! Mean diameter at RH 80% ~2 x dry diameter for sea salt,
                ! organic part is assumed insoluble

                if( Ovf(i,j,k) < 0.)then ! first timestep
                        Ovf(i,j,k) = omf_max(i,j)/(omf_max(i,j)*(1.-rhorho)+rhorho)
                endif
                dia80 = (dia**3.*(Ovf(i,j,k)+8.*(1.-Ovf(i,j,k))))**(1./3.) !
                omf = omf_max(i,j)/(1.+0.03*exp(6.81*dia80))+0.03*omf_max(i,j)
                Ovf(i,j,k) = omf/(omf*(1.-rhorho)+rhorho)
            else
                Ovf(i,j,k) = 0.
            endif

            ! Different assumptions can be made about whether the organic fraction replaces the salt or is additive to it.
            ! If assumed additive, increase the emission flux accordingly
            if(ifPOCadd) dcdt(i,j,k) = dcdt(i,j,k) / (1. - Ovf(i,j,k))

          enddo
        enddo


        ! Apply to 2b bins, if possible
        IF (nbins<in2b) THEN
            i = k
        ELSE
            i = in2b + k - in2a
        ENDIF

        ! Note: bin center is volume mean - using other than that will cause problems!
        vdry = pi6*4.*(aerobins(k)**3+aerobins(k+1)**3)

        a_naerot(2,:,:,i) = a_naerot(2,:,:,i) + dcdt(:,:,k)
        ! ... and specifically to SS
        if(iss>0)then
            j = (iss-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*rhoss*(1.-ovf(:,:,k))
            !  ... and just add water at 80% RH (D80 = 2 x Ddry)
            j = (ih2o-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*rhowa*(1.-ovf(:,:,k))*7.
        endif
        ! .. organic fraction
        if (nvbs_setup>=0) then
            j = (vbs_set(1)%id_vols-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*spec_density(vbs_set(1)%spid)*(ovf(:,:,k))
        elseif(ioc>0)then
            j = (ioc-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*rhooc*(ovf(:,:,k))
        endif

    ENDDO

    if (sflg) then
        ! Convert dcdt [#/kg/s] to flux [#/m^2/s] and take sum over bins
        w(:,:) = SUM(dcdt(:,:,:),DIM=3)*a_dn(2,:,:)*(zm(3)-zm(2))
        call flux_stat(nxp,nyp,w,1)
    endif

  END SUBROUTINE marine_aero_flux

  ! --------------------------------------------------------------------------
  ! Gas emissions from ocean surface
  !
  SUBROUTINE marine_gas_flux(sst)
    ! Parameters
    use defs, only: vonk, g
    use grid, only: nxp, nyp, a_ustar, a_dn, zm, a_gaerot
    use mo_vbsctl, only: vbs_voc_set
    use mo_vbs, only: spec_moleweight
    use stat, only: flux_stat, sflg
    IMPLICIT NONE
    REAL, INTENT(IN) :: sst     ! Sea surface temperature (K)

    ! Local
    INTEGER :: i, j, k
    REAL :: usum, zs
    REAL :: u10(nxp,nyp)! whitecap cover, 10 m wind speed
    REAL :: flxIsop(nxp,nyp), flxMtrp(nxp,nyp) ! Gas flux (kg/m2/s)

    real, parameter :: schmidt_ref = 660.0, & ! CO2 in 293 K
 & per_sec = 1.0 / 3600.0, to_m = 0.01, p23 = 2.0 / 3.0, p12 = 0.5
    real :: schmidt_isoprene, schmidt_monoterp, sc_ratio, temp_c
    integer :: iGas

    flxIsop(:,:) = 0.
    flxMtrp(:,:) = 0.
    if(.not. (wtrMtrp > 0. .or. wtrIsop > 0.))return

    ! Roughness height is needed for the 10 m wind speeds
    zs = zrough
    IF (zrough <= 0.) THEN ! Calculate
        usum = 0.
        DO j=3,nyp-2
            DO i=3,nxp-2
                usum = usum + a_ustar(i,j)
            END DO
        ENDDO
        usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
        zs = max(0.0001,(0.016/g)*usum**2)
    ENDIF

    ! Whitecap cover (% => fraction) based on 10 m wind speeds (eq. 2)
    u10(:,:) = a_ustar(:,:)/vonk*log(10.0/zs)

    temp_c = sst - 273.15

    ! Saltzman et al., 1993 (might need in the future)
    !schmidt_dms = 2674.0 - 147.12 * temp_c + 3.726 * temp_c**2 + 0.038*temp_c**3

    ! 1. Isoprene
    ! Palmer & Shaw, 2005
    schmidt_isoprene = 3913.15 - 162.13*temp_c + 2.67*temp_c**2 - 0.012*temp_c**3
    if(wtrIsop > 0.)then
      sc_ratio = schmidt_ref / schmidt_isoprene
      iGas = vbs_voc_set(2)%id_gas
      DO j=3,nyp-2
        DO i=3,nxp-2
          ! First compute the transfer velocity (m/s)
!          select case (transfer_velocity_type)
!          case(Liss_Merlivat) ! Liss & Merlivat 1986
!            if (windspeed <= 3.6) then
!              flxIsop = 0.17 * windspeed * sc_ratio ** p23
!            else if (windspeed <= 13.0) then
!              flxIsop = 0.612 * sc_ratio**p23 + (2.85*windspeed - 10.26) * sc_ratio**p12
!            else
!              flxIsop = 0.612 * sc_ratio**p23 + (5.90*windspeed - 49.90) * sc_ratio**p12
!            end if
!            flxIsop = flxIsop * per_sec*to_m
!          case (Wannikhof) ! Wannikhof (2014)
            flxIsop(i,j) = 0.251*u10(i,j)*u10(i,j)*sc_ratio**p12 * per_sec*to_m
!          case default ! No flux
!            flxIsop = 0.
!          end select
          ! Flux
          flxIsop(i,j) = flxIsop(i,j) * wtrIsop * spec_moleweight(iGas) /1000.

          a_gaerot(2,i,j,iGas) = a_gaerot(2,i,j,iGas) + flxIsop(i,j) / a_dn(2,i,j) / (zm(3)-zm(2))
        enddo
      enddo
    endif

    ! 2. Monoterpenes
    if(wtrMtrp > 0.)then
      iGas = vbs_voc_set(1)%id_gas
      ! Temperature dependent Schmidt numbers for monoterpenes do not seem to be available.
      ! Following Moore & Grozsko (1999), assume that the ratio of the diffusivities
      ! is inversely proportional to the ratio of the molar volumes to the power 0.6
      ! [Wilke and Chang, 1955].
      ! Data from chemspider.com:
      ! Molar volume of DMS:      75.5±3.0 cm3
      ! Molar volume of isoprene: 101.1±3.0 cm3
      ! Molar Volume of a-pinene: 154.9±3.0 cm3
      ! Molar Volume of limonene: 163.3±3.0 cm3
      ! Take average of a-pinene and limonene - 159.1
      schmidt_monoterp = schmidt_isoprene*(159.1/101.1)**0.6
      sc_ratio = schmidt_ref / schmidt_monoterp
      DO j=3,nyp-2
        DO i=3,nxp-2
          ! Transfer velocity (m/s)
          flxMtrp(i,j) = 0.251*u10(i,j)*u10(i,j)*sc_ratio**p12 * per_sec*to_m
          ! Flux
          flxMtrp(i,j) = flxMtrp(i,j) * wtrMtrp * spec_moleweight(iGas) /1000.
          a_gaerot(2,i,j,iGas) = a_gaerot(2,i,j,iGas) + flxMtrp(i,j) / a_dn(2,i,j) / (zm(3)-zm(2))
        enddo
      enddo
    endif

    if (sflg) then
        call flux_stat(nxp,nyp,flxIsop,2)
        call flux_stat(nxp,nyp,flxMtrp,3)
    endif

  END SUBROUTINE marine_gas_flux

  !
  ! --------------------------------------------------------------------------
  ! SURFACE: Calculates surface fluxes using an algorithm chosen by ISFCTYP
  ! and fills the appropriate 2D arrays
  !
  !     default: specified thermo-fluxes (drtcon, dthcon)
  !     isfctyp=1: specified surface layer gradients (drtcon, dthcon)
  !     isfctyp=2: fixed lower boundary of water at certain sst
  !     isfctyp=3: bulk aerodynamic law with coefficients (drtcon, dthcon)
  !     isfctyp=4: regulate surface temperature to yield a constant surface buoyancy flux
  !
  ! Modified for level 4
  ! Juha Tonttila, FMI, 2014
  !

  subroutine surface(sst)

    use defs, only: vonk, p00, rcp, g, cp, alvl, ep2
    use grid, only: nzp, nxp, nyp, a_up, a_vp, a_theta, a_rv, a_rp, zt, psrf, th00  &
         , umean, vmean, a_ustar, a_tstar, a_rstar, uw_sfc, vw_sfc, ww_sfc    &
         , wt_sfc, wq_sfc, obl, dn0, level,dtl, a_sflx, a_rflx, precip, W1,W2,W3
    use thrm, only: rslf
    use stat, only: sfc_stat, sflg
    use mpi_interface, only : nypg, nxpg, double_array_par_sum


    implicit none
    real, optional, intent (inout) :: sst
    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp)       &
         ,wspd(nxp,nyp), bfct(nxp,nyp)
    real :: rx(nzp,nxp,nyp)


    real :: total_sw, total_rw, total_la, total_se, total_pre  ! Sami added
    real :: C_heat,lambda ! Sami added
    real :: K1,K2,K3,Kmean1,Kmean2,fii_1,fii_2,fii_3,Q3,Q12,Q23,ff1  ! Sami added



    integer :: i, j, iterate
    real    :: zs, bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
    real (kind=8) :: bfl(2), bfg(2)

    ! Added by Juha
    SELECT CASE(level)
       CASE(1,2,3)
          rx = a_rv
       CASE(4,5)
          rx = a_rp
    END SELECT

    select case(isfctyp)
    !
    ! use prescribed surface gradients dthcon, drton from NAMELIST
    ! use then similarity theory to compute the fluxes
    !
    case(1)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j)=dthcon
             drdz(i,j)=drtcon
             bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
          end do
       end do
       zs = zrough
       call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar     &
            ,a_rstar,obl)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
    !
    ! use prescribed SST and assume qsurf=qsat (i.e. ocean) to compute
    ! gradients. Then use similarity theory to predict the fluxes.
    !
    case(2)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       usum = 0.
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
             bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
             usum = usum + a_ustar(i,j)
          end do
       end do
       usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
       zs = zrough
       if (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
       call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar     &
            ,a_rstar,obl)
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
    !
    ! drtcon (wq=Ch*u*dth, Garrat p.55) and using prescribed sst
    ! and qsurf=qsat (ocean); note that here zrough is not the roughness
    !
    case(3)
       call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
       do j=3,nyp-2
          do i=3,nxp-2
             dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
             drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
             if (ubmin > 0.) then
                a_ustar(i,j) = sqrt(zrough)* wspd(i,j)
             else
                a_ustar(i,j) = abs(ubmin)
             end if
             a_tstar(i,j) =  dthcon * wspd(i,j)*dtdz(i,j)/a_ustar(i,j)
             a_rstar(i,j) =  drtcon * wspd(i,j)*drdz(i,j)/a_ustar(i,j)
             bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
          end do
       end do
       call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
            ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
    !
    ! fix surface temperature to yield a constant surface buoyancy flux dthcon
    !
    case(4)

       Vzt   = 10.* (log(zt(2)/zrough)/log(10./zrough))
       Vbulk = Vzt * (vonk/log(zt(2)/zrough))**2

       bfl(:) = 0.
       do j=3,nyp-2
          do i=3,nxp-2
             bfl(1) = bfl(1)+a_theta(2,i,j)
             bfl(2) = bfl(2)+rx(2,i,j) ! Juha: rx
          end do
       end do

       call double_array_par_sum(bfl,bfg,2)

       bfg(2) = bfg(2)/real((nxpg-4)*(nypg-4))
       bfg(1) = bfg(1)/real((nxpg-4)*(nypg-4))

       do iterate=1,5
          bflx  = ((sst -bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst) -bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
          sst1 = sst + 0.1
          bflx1 = ((sst1-bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst1)-bfg(2))) &
               * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
          sst  = sst + 0.1* (dthcon - bflx) / (bflx1-bflx)
       end do

       do j=3,nyp-2
          do i=3,nxp-2
             wt_sfc(i,j) = Vbulk * (sst -a_theta(2,i,j))
             wq_sfc(i,j) = Vbulk * (rslf(psrf,sst) - rx(2,i,j)) ! Juha: rx
             wspd(i,j)    = max(0.1,                                    &
                  sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
             bflx         = wt_sfc(i,j)*g/bfg(1) + g*ep2*wq_sfc(i,j)
             a_ustar(i,j) = diag_ustar(zt(2),zrough,bflx,wspd(i,j))
             uw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_up(2,i,j)+umean)/wspd(i,j)
             vw_sfc(i,j)  = -a_ustar(i,j)*a_ustar(i,j)                  &
                  *(a_vp(2,i,j)+vmean)/wspd(i,j)
             ww_sfc(i,j)  = 0.
             a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
          end do
       end do
    !
    !Sami addition: Calculate surface fluxes from energy budget
    ! This is based on Acs et al: A coupled soil moisture and surface
    ! temperature prediction model, Journal of applied meteorology, 30, 1991
    !
    case(5)

       total_sw = 0.0
       total_rw = 0.0
       total_la = 0.0
       total_se = 0.0
       total_pre = 0.0
       ffact = 1.

       !
       !   Calculate mean energy fluxes(Mean for each proceccors)
       !

        do j=3,nyp-2
           do i=3,nxp-2
              total_sw = total_sw+a_sflx(2,i,j)
              total_rw = total_rw + a_rflx(2,i,j)
              total_la = total_la + wq_sfc(i,j)*(0.5*(dn0(1)+dn0(2))*alvl)/ffact
              total_se = total_se + wt_sfc(i,j)*(0.5*(dn0(1)+dn0(2))*cp)/ffact
              total_pre=   total_pre +  precip(2,i,j)
           end do
        end do
        total_sw =  total_sw/real((nxp-4)*(nyp-4))
        total_rw =  total_rw/real((nxp-4)*(nyp-4))
        total_la =  total_la/real((nxp-4)*(nyp-4))
        total_se = total_se/real((nxp-4)*(nyp-4))
        total_pre = total_pre/real((nxp-4)*(nyp-4))


        ! From energy fluxes calculate new sirface temperature
        sst1 =sst
        C_heat = 1.29e3*(840+4187*thetaS1*W1) ! Eq 33
        lambda=1.5e-7*C_heat

        ! Determine moisture at different depths in the ground
        !

        K1=K_s1*W1**(2.*B1+3.)
        K2=K_s1*W2**(2.*B1+3.)
        K3=K_s3*W3**(2.*B3+3.)
        Kmean1=(D1*K1+D2*K2)/(D1+D2)
        Kmean2=(D2*K2+D3*K3)/(D2+D3)

        fii_1 = fii_s1*W1**(-B1)
        fii_2 = fii_s1*W2**(-B1)
        fii_3 = fii_s3*W3**(-B3)

        Q3 = K_s3*W3**(2.*B3+3.)*sin(0.05236) ! 3 degrees  Eq 8

        Q12 = Kmean1*2.*( (fii_1-fii_2)/(D1+D2)+1.0)
        Q23 = Kmean2*2.*( (fii_2-fii_3)/(D2+D3)+1.0)

        W1 = W1+1./(thetaS1*D1)*(-Q12-((total_la+total_pre)/((0.5*(dn0(1)+dn0(2))*alvl)/ffact))/(thetaS1*1000))*dtl
        W2 = W2+1./(thetaS2*D2)*(Q12-Q23)*dtl
        W3 = W3+1./(thetaS3*D3)*(Q23-Q3)*dtl
        !
        !  Following is copied from case (2). No idea if this is valid or not..
        !

        call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
        usum = 0.
        do j=3,nyp-2
           do i=3,nxp-2

              dtdz(i,j) = a_theta(2,i,j) - sst1*(p00/psrf)**rcp

              ff1=1.0
              IF(W1<0.75) ff1=W1/0.75
              ! Flux of moisture is limited by water content.
              drdz(i,j) = rx(2,i,j) - ff1*rslf(psrf,min(sst1,280.))
              !
              bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
              usum = usum + a_ustar(i,j)
           end do
        end do
        usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
        zs = zrough
        if (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
        call srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar     &
             ,a_rstar,obl)
        call sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar  &
             ,uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

           sst1 = sst1-(total_rw+total_la+total_se+((lambda*C_heat*7.27e-5/(2.0))**0.5*(SST1-280.0)))&
                /(2.0e-2*C_heat+(lambda*C_heat/(2.0*7.27e-5))**0.5)*dtl

           sst = sst1
    !
    ! fix thermodynamic fluxes at surface given values in energetic
    ! units and calculate momentum fluxes from similarity theory
    !
    case default
       ffact = 1.
       wt_sfc(1,1)  = ffact* dthcon/(0.5*(dn0(1)+dn0(2))*cp)
       wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)

       if (zrough <= 0.) then
          usum = 0.
          do j=3,nyp-2
             do i=3,nxp-2
                usum = usum + a_ustar(i,j)
             end do
          end do
          usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
          zs = max(0.0001,(0.016/g)*usum**2)
       else
          zs = zrough
       end if

       do j=3,nyp-2
          do i=3,nxp-2
             wt_sfc(i,j)=wt_sfc(1,1)
             wq_sfc(i,j)=wq_sfc(1,1)

             wspd(i,j)    = max(0.1,                                    &
                  sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
             if (ubmin > 0.) then
                bflx = g*wt_sfc(1,1)/th00
                if (level >= 2) bflx = bflx + g*ep2*wq_sfc(i,j)
                a_ustar(i,j) = diag_ustar(zt(2),zs,bflx,wspd(i,j))
             else
                a_ustar(i,j) = abs(ubmin)
             end if

             ffact = a_ustar(i,j)*a_ustar(i,j)/wspd(i,j)
             uw_sfc(i,j)  = -ffact*(a_up(2,i,j)+umean)
             vw_sfc(i,j)  = -ffact*(a_vp(2,i,j)+vmean)
             ww_sfc(i,j)  = 0.
             a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
             a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
          end do
       end do

    end select

    if (sflg) call sfc_stat(nxp,nyp,wt_sfc,wq_sfc,a_ustar,sst)

    return
  end subroutine surface
  !
  ! -------------------------------------------------------------------
  ! GET_SWNDS: returns surface winds valid at cell centers
  !
  subroutine get_swnds(n1,n2,n3,usfc,vsfc,wspd,up,vp,umean,vmean)

    implicit none

    integer, intent (in) :: n1, n2, n3
    real, intent (in)    :: up(n1,n2,n3), vp(n1,n2,n3), umean, vmean
    real, intent (out)   :: usfc(n2,n3), vsfc(n2,n3), wspd(n2,n3)

    integer :: i, j, ii, jj

    do j=3,n3-2
       jj = j-1
       do i=3,n2-2
          ii = i-1
          usfc(i,j) = (up(2,i,j)+up(2,ii,j))*0.5+umean
          vsfc(i,j) = (vp(2,i,j)+vp(2,i,jj))*0.5+vmean
          wspd(i,j) = max(abs(ubmin),sqrt(usfc(i,j)**2+vsfc(i,j)**2))
       enddo
    enddo

  end subroutine get_swnds
  !
  ! ----------------------------------------------------------------------
  ! FUNCTION GET_USTAR:  returns value of ustar using the below
  ! similarity functions and a specified buoyancy flux (bflx) given in
  ! kinematic units
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  real function diag_ustar(z,z0,bflx,wnd)

    use defs, only: vonk

    implicit none

    real, parameter      :: am   =  4.8   !   "          "         "
    real, parameter      :: bm   = 19.3   !   "          "         "
    real, parameter      :: eps  = 1.e-10 ! non-zero, small number

    real, intent (in)    :: z             ! height where u locates
    real, intent (in)    :: z0            ! momentum roughness height
    real, intent (in)    :: bflx          ! surface buoyancy flux (m^2/s^3)
    real, intent (in)    :: wnd           ! wind speed at z

    integer :: iterate
    real    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

    lnz   = log(z/z0)
    klnz  = vonk/lnz
    c1    = 3.14159/2. - 3.*log(2.)

    ustar =  wnd*klnz
    if (bflx /= 0.0) then
       do iterate=1,4
          lmo   = -(ustar**3)/(bflx*vonk + eps)
          zeta  = z/lmo
          if (zeta > 0.) then
             ustar =  vonk*wnd  /(lnz + am*zeta)
          else
             x     = sqrt( sqrt( 1.0 - bm*zeta ) )
             psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
             ustar = wnd*vonk/(lnz - psi1)
          end if
       end do
    end if

    diag_ustar = ustar

    return
  end function diag_ustar
  !
  ! ----------------------------------------------------------------------
  ! Subroutine srfcscls:  returns scale values based on Businger/Dye
  ! similarity functions.
  !
  ! phi_h (zeta > 0) =  Pr * (1 + ah * zeta)
  ! phi_h (zeta < 0) =  Pr * (1 - bh * zeta)^(-1/2)
  !
  ! phi_m (zeta > 0) =  (1 + am * zeta)
  ! phi_m (zeta < 0) =  (1 - bm * zeta)^(-1/4)
  !
  ! where zeta = z/lmo and lmo = (theta_rev/g*vonk) * (ustar^2/tstar)
  !
  ! Ref: Businger, 1973, Turbulent Transfer in the Atmospheric Surface
  ! Layer, in Workshop on Micormeteorology, pages 67-100.
  !
  ! Code writen March, 1999 by Bjorn Stevens
  !
  subroutine srfcscls(n2,n3,z,z0,th00,u,dth,drt,ustar,tstar,rstar,obl)

    use defs, only : vonk, g, ep2
    use grid, only: runtype

    implicit none

    ! MALTE
    !real, parameter     :: ah   =  4.7   ! stability function parameter
    !real, parameter     :: bh   = 15.0   !   "          "         "
    !real, parameter     :: am   =  4.7   !   "          "         "
    !real, parameter     :: bm   = 15.0   !   "          "         "
    ! Original
    real, parameter     :: ah   =  7.8   ! stability function parameter
    real, parameter     :: bh   = 12.0   !   "          "         "
    real, parameter     :: am   =  4.8   !   "          "         "
    real, parameter     :: bm   = 19.3   !   "          "         "
    !
    real, parameter     :: pr   = 0.74   ! prandlt number
    real, parameter     :: eps  = 1.e-10 ! non-zero, small number

    integer, intent(in) :: n2,n3         ! span of indicies covering plane
    real, intent(in)    :: z             ! height where u & T locate
    real, intent(in)    :: z0            ! momentum roughness height
    real, intent(in)    :: th00          ! reference temperature
    real, intent(in)    :: u(n2,n3)      ! velocities at z
    real, intent(in)    :: dth(n2,n3)    ! theta (th(z) - th(z0))
    real, intent(in)    :: drt(n2,n3)    ! qt(z) - qt(z0)
    real, intent(inout) :: ustar(n2,n3)  ! scale velocity
    real, intent(inout) :: tstar(n2,n3)  ! scale temperature
    real, intent(inout) :: rstar(n2,n3)  ! scale value of qt
    real, intent(inout) :: obl(n2,n3)    ! Obukhov Length

    logical, save :: first_call=.True.
    integer :: i,j,iterate
    real    :: lnz, klnz, betg
    real    :: x, psi1, psi2, Lold, Ldif, zeff, zeta, lmo, dtv
    logical    :: exititer

    lnz   = log(z/z0)
    klnz  = vonk/lnz
    betg  = th00/g

    do j=3,n3-2
       do i=3,n2-2
          dtv = dth(i,j) + ep2*th00*drt(i,j)

          !
          ! Neutral case
          !
          if (dtv == 0.) then
             ustar(i,j) = u(i,j)*klnz
             tstar(i,j) = dtv*klnz/pr
             lmo        = -1.e10

          !
          ! start iterations from values at previous tstep,
          ! unless the sign has changed or if it is the first call, then
          ! use neutral values.
          !
          else
             if ((runtype=='INITIAL' .and. first_call) .or. (tstar(i,j)*dtv <= 0.)) then
                ustar(i,j) = u(i,j)*klnz
                tstar(i,j) = dtv*klnz/pr
                lmo        = -1.e10
             end if

             if(ustar(i,j) == 0) ustar(i,j) = 0.1

             Lold  = 1e9
             Ldif  = 1e9
             iterate  = 0
             exititer = .false.

             do while(abs(Ldif)>0.1)
                lmo   = betg*ustar(i,j)**2/(vonk*tstar(i,j))
                Ldif       = lmo - Lold
                Lold       = lmo

                if ((dtv < 0) .and. (lmo > -0.001)) lmo = -0.001999
                if ((dtv > 0) .and. (lmo < +0.001)) lmo = +0.001777

                ! BvS : Following ECMWF, limit z/L for very stable conditions
                if(z/lmo > 5.) then
                   zeff = lmo * 5.
                   exititer = .true.
                else
                   zeff = z
                end if

                zeta  = zeff/lmo
                if(zeta <= 0) then
                    x     = (1.0 - bm*zeta )**0.25
                    psi1  = 3.14159265/2.- 2.*atan(x) + 2.*log((1.+x)/2.) + log((1.+x**2)/2.)
                    x     =  (1.0 - bh*zeta)**0.25
                    psi2  = 2.*log( (1. + x**2)/2. )
                ELSE
                    psi1  = - am*zeta
                    psi2  = - ah*zeta
                END IF

                ustar(i,j) = u(i,j)*vonk/(log(zeff/z0) - psi1)
                if(ustar(i,j)<0.) ustar(i,j) = 0.1
                tstar(i,j) = (dtv*vonk/pr)/(log(zeff/z0) - psi2)

                if(exititer) then
                   lmo        = zeff/5.
                   exit
                end if

                iterate = iterate + 1

                ! Limit L for day/night transitions
                if(lmo > 1e6)  lmo = 1e6
                if(lmo < -1e6) lmo = -1e6

                if(iterate>10000) stop 'Obukh. length not converged!'
             end do
          end if

          obl(i,j) = lmo
          rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
          tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)
       end do
    end do

    first_call = .False.

    return
  end subroutine srfcscls
  !
  ! ----------------------------------------------------------------------
  ! subroutine: sfcflxs:  this routine returns the surface fluxes based
  ! on manton-cotton algebraic surface layer equations.
  !
  subroutine sfcflxs(n2,n3,vk,ubar,u,v,xx,us,ts,rs,uw,vw,tw,rw,ww)
    implicit none
    real, parameter      :: cc=4.7,eps=1.e-20

    integer, intent(in)  :: n2,n3
    real, intent(in)     :: ubar(n2,n3),u(n2,n3),v(n2,n3),xx(n2,n3),vk
    real, intent(in)     :: us(n2,n3),ts(n2,n3),rs(n2,n3)
    real, intent(out)    :: uw(n2,n3),vw(n2,n3),tw(n2,n3),rw(n2,n3),ww(n2,n3)

    real    :: x(n2,n3),y(n2,n3)
    integer i,j

    do j=3,n3-2
       do i=3,n2-2

          uw(i,j)=-(u(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          vw(i,j)=-(v(i,j)/(ubar(i,j)+eps))*us(i,j)**2
          tw(i,j)=-ts(i,j)*us(i,j)
          rw(i,j)=-rs(i,j)*us(i,j)

          x(i,j) = xx(i,j)*vk*ts(i,j)*(ubar(i,j)/us(i,j))**2
          x(i,j) = x(i,j)*sqrt(sqrt(1.-15.*min(0.,x(i,j)))) &
               /(1.0+cc*max(0.,x(i,j)))
          y(i,j) =sqrt((1.-2.86*x(i,j))/(1.+x(i,j)* &
               (-5.39+x(i,j)*6.998 )))
          ww(i,j)=(0.27*max(6.25*(1.-x(i,j))*y(i,j),eps)-&
               1.18*x(i,j)*y(i,j))*us(i,j)**2
       enddo
    enddo
    return
  end subroutine sfcflxs

end module srfc



