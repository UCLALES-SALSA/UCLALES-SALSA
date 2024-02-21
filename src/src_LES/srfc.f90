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
  real    :: zrough_t = -999. ! Optional roughness height for temperature
  real    :: ubmin  =  0.20
  real    :: dthcon = 100.0
  real    :: drtcon = 0.0



! Sami added ----->
! Initial values for surface properties
  real ::  W1 = 0.9   !Water content
  real ::  W2 = 0.9
  real ::  W3 = 0.9
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
  ! Option for sea-spray aerosol source function parameterization
  integer :: ssa_param = 0
  ! Volume fraction of organic matter in sea spray
  real, allocatable, save :: ovf(:)

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
  SUBROUTINE marine_aero_flux(sst)
    ! Parameters
    use defs, only: vonk, g
    use grid, only: nxp, nyp, a_ustar, a_dn, zm, nbins, a_naerot, a_maerot
    use mo_submctl, only: aerobins, in2a, fn2a, &
                    in2b, pi6, iss, ioc, ih2o, dens, nspec
    use stat, only: fill_scalar, sflg
    use util, only : get_avg2dh

    IMPLICIT NONE
    REAL, INTENT(IN) :: sst     ! Sea surface temperature (K)
    ! Local
    INTEGER :: i, j, k
    REAL :: usum, zs, dia, rhorho, dia80, omf, omf_max, vdry, u10_bar
    REAL :: u10(nxp,nyp) ! 10 m wind speed
    REAL :: flx(fn2a+1) ! Production rate for each size bin
    REAL :: dcdt(nxp,nyp,fn2a) ! Particle concentration tendency (#/kg/s)

    IF (wtrChlA>0. .and. ioc<1) THEN
        STOP 'No sea spray species (OC) included!'
    ELSEIF (iss<1 .and. nspec==1) THEN
        ! If single component aerosol, then use the only one available
        iss=2
    ELSEIF (iss<1) THEN
        STOP 'No sea spray species (SS) included!'
    ENDIF

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
    u10(:,:) = a_ustar(:,:)/vonk*log(10.0/zs)
    ! Mean 10 m wind speed (for the domain)
    u10_bar = get_avg2dh(nxp,nyp,u10)

    ! Calculate particle flux, dF/dlogDp, for each size bin limit
    DO k=1,fn2a+1
        dia=aerobins(k)*2.
        select case (ssa_param)
        case(0)
            ! Combining Martensson et al. (2003), Monahan et al. (1986), and Smith and Harrison (1998)
            flx(k) = flux_mmsh(dia,sst,u10_bar)
        case(1)
            ! Gong (2003)
            flx(k) = flux_gong(dia,sst,u10_bar)
        case(2)
            ! Grythe et al. (2014)
            flx(k) = flux_grythe(dia,sst,u10_bar)
        case(3)
            ! Sofiev et al. (2011)
            flx(k) = flux_sofiev(dia,sst,u10_bar)
        case default
            STOP 'Sea-spray source function not supported!'
        end select
    ENDDO

    ! Organic mass fraction grom Gantt ea (2011)
    ! conversion to mg from kg
    if(wtrChlA>0.)then
        rhorho = dens(ioc)/dens(iss) ! ratio of densities of marine organic matter and salt
        omf_max = 1./(1.+exp(-2.63e6*wtrChlA+0.18*u10_bar))
    else
        omf_max = 0.
        rhorho = 1.
    endif
    ! Growth factor depends on organic volume fraction
    if (.not. allocated(ovf)) THEN
        ! first call
        ALLOCATE(ovf(fn2a+1))
        ovf(:) = omf_max/(omf_max*(1.-rhorho)+rhorho)
    endif

    if(wtrChlA>0.)then
        ! Size dependent organic mass fraction from Gantt ea, 2011 (eq. 3)
        DO k=in2a,fn2a+1 ! Ignore 1a, because there is no sea salt
            dia=aerobins(k)*2.e6  !compute for bin borders, average later
            ! Organic fraction is parameterized for particle aerodynamic diameter at 80% RH
            ! which depends on composition, so technically would need to iterate.
            ! However, hopefully it converges after limited nr of timesteps.
            ! Mean diameter at RH 80% ~2 x dry diameter for sea salt,
            ! organic part is assumed insoluble
            ! Correction factor from stokes diameter to aerodynamic diameter
            ! is 1.07 for sea salt at 80% RH and 1.15 for dry OC, so on average 1.1
            dia80 = (dia**3*(ovf(k)+8.*(1.-ovf(k))))**(1./3.)* 1.1
            omf = omf_max/(1.+0.03*exp(6.81*dia80))+0.03*omf_max
            ovf(k) = omf/(omf*(1.-rhorho)+rhorho)
        ENDDO
    else
        ovf(:) = 0.
    endif

    ! Particle concentration tendency for each size bin
    dcdt(:,:,:)=0.
    DO k=in2a,fn2a ! Ignore 1a, because there is no sea salt
        ! Mean flux: 0.5*(flx(k)+flx(k+1))
        ! From dFp/dlog(Dp) to dFp: multiply by log10(rad(k+1)/rad(k))
        ! Convert #/m^2/s to the rate of change in concentration (#/kg/s): multiply by 1/rho/dz
        DO j=3,nyp-2
          DO i=3,nxp-2
            dcdt(i,j,k)=0.5*(flx(k)+flx(k+1))*log10(aerobins(k+1)/aerobins(k))/a_dn(2,i,j)/(zm(3)-zm(2))

            ! Different assumptions can be made about whether the organic fraction replaces the salt or is additive to it.
            ! If assumed additive, increase the emission flux accordingly
            if(ifPOCadd) dcdt(i,j,k) = dcdt(i,j,k) / (1. - 0.5*(ovf(k)+ovf(k+1)))

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
                                              & vdry*dens(iss)*(1.-0.5*(ovf(k)+ovf(k+1)))
            !  ... and add water proportionally to the salt volume
            j = (ih2o-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*dens(ih2o)*(1.-0.5*(ovf(k)+ovf(k+1)))*7.
        endif
        ! .. organic fraction
        if(ioc>0)then
            j = (ioc-1)*nbins + i
            a_maerot(2,:,:,j) = a_maerot(2,:,:,j) + dcdt(:,:,k)* &
                                              & vdry*dens(ioc)*0.5*(ovf(k)+ovf(k+1))
        endif

    ENDDO

    if (sflg) then
        ! Convert dcdt [#/kg/s] to flux [#/m^2/s] and take sum over bins
        omf = SUM ( SUM( SUM(dcdt(3:nxp-2,3:nyp-2,:),DIM=3)*a_dn(2,3:nxp-2,3:nyp-2),DIM=2 ) )*(zm(3)-zm(2))/REAL((nxp-4)*(nyp-4))
        call fill_scalar(omf,'flx_aer')
        ! 10 m wind speed
        call fill_scalar(u10_bar,'u10    ')
    endif

  END SUBROUTINE marine_aero_flux


  ! -------------------------------------------------------------------
  ! Sea spray emissions parameterizations
  !   Inputs: dry diameter (Dp, m), sea surface temperature (SST, K) and 10 m wind speed (u10, m/s)
  !   Output: particle number flux, dF/dlog10(Dp), #/m^2/s
  !
  ! References
  !  Gong, S. L.: A parameterization of sea-salt aerosol source function for sub- and super-micron
  !     particles, Global Biogeochem. Cycles, 17, 1097, doi:10.1029/2003GB002079, 2003
  !  Grythe, H., Strom, J., Krejci, R., Quinn, P., and Stohl, A.: A review of sea-spray aerosol
  !     source functions using a large global set of sea salt aerosol concentration measurements,
  !     Atmos. Chem. Phys., 14, 1277-1297, https://doi.org/10.5194/acp-14-1277-2014, 2014
  !  Jaegle, L., Quinn, P. K., Bates, T. S., Alexander, B., and Lin, J.-T.: Global distribution of
  !     sea salt aerosols: new constraints from in situ and remote sensing observations,
  !     Atmos. Chem. Phys., 11, 3137-3157, https://doi.org/10.5194/acp-11-3137-2011, 2011.
  ! Martensson, E. M., Nilsson, E. D., de Leeuw, G., Cohen, L. H., and Hansson, H.-C.:
  !     Laboratory simulations and parameterization of the primary marine aerosol production,
  !     J. Geophys. Res., 108, 4297, doi:10.1029/2002JD002263, 2003
  ! Monahan, E. C., D. E. Spiel, and K. L. Davidson: A model of marine aerosol generation via
  !    whitecaps and wave disruption, in Oceanic Whitecaps, edited by E. Monahan, and G. M. Niocaill,
  !    pp. 167-174, D. Reidel, Norwell, Mass., 1986
  ! Smith, M.H., and Harrison,  N.M.: The sea spray generation function, J. Aerosol Sci.,
  !    Volume 29, Supplement 1, S189-S190, https://doi.org/10.1016/S0021-8502(98)00280-8, 1998
  ! Sofiev, M., Soares, J., Prank, M., de Leeuw, G., and Kukkonen, J.: A regional-to-global
  !    mode of emission and transport of sea salt particles in the atmosphere,
  !    J. Geophys. Res., 116, D21302, doi:10.1029/2010JD014713, 2011

  REAL FUNCTION flux_mmsh(dia,sst,u10)
    ! Flux parameterization combining those from Martensson et al. (2003), Monahan et al. (1986), and Smith and Harrison (1998)
    ! Valid from 0.02e-6 to 300e-6 m
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    !
    if (dia<=1e-6) then
        flux_mmsh = flux_martensson(dia,sst,u10)
    elseif (dia<=10e-6) then
        flux_mmsh = flux_monahan(dia,sst,u10)
    else
        flux_mmsh = flux_smith(dia,sst,u10)
    endif
    !
  END FUNCTION flux_mmsh

  REAL FUNCTION flux_martensson(dia,sst,u10)
    ! Flux parameterization from Martensson et al. (2003)
    ! Valid from 0.02e-6 to 2.8e-6 m
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL :: Ak, Bk, flx, w
    !
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
    ENDIF
    ! Eq. 6: particle flux per whitecap area, dFp/dlog(Dp) [#/m^2/s]
    flx=MAX(0.,Ak*sst+Bk)
    !
    ! Whitecap cover (% => fraction) based on 10 m wind speeds (eq. 2)
    w=3.84e-6*u10**3.41
    !
    ! Multiply by white cap cover
    flux_martensson=flx*w
  END FUNCTION flux_martensson

  REAL FUNCTION flux_monahan(dia,sst,u10)
    ! Flux parameterization from Monahan et al. (1986)
    ! Valid from 0.8e-6 to 8e-6 m, but reasonable values from about 0.1e-6 m up to 50e-6 m
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL :: rad80, b, flx, tw, twt
    !
    ! Particle radius at 80 RH% is about the same as the dry diameter (assuming GF=2)
    rad80=dia*1e6 ! Radius in microns
    !
    ! Flux parameterization, dF/dr, Eq. 1 in Gong (2003) or Eq. A2 in Grythe et al. (2014)
    b = (0.38-log10(rad80))/0.65
    flx = 1.373*u10**3.41/rad80**3*(1.+0.057*rad80**1.05)*10.**(1.19*exp(-b**2))
    !
    ! Additional: temperature dependency from Jaegle et al. (2011)
    tw = sst-273.15 ! Temperature in C
    twt = 0.3+0.1*tw-0.0076*tw**2+0.00021*tw**3
    !
    ! Convert from dF/dDp to dF/dlog10Dp: multiply by Dp*ln(10)
    flux_monahan = rad80*log(10.)*flx*twt
  END FUNCTION flux_monahan

  REAL FUNCTION flux_smith(dia,sst,u10)
    ! Flux parameterization from Smith and Harrison (1998)
    ! Valid from 1e-6 to 300e-6 m
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL :: flx, tw, twt
    !
    ! Partial particle flux, dF/dr80 = dF/dDp
    flx = 0.2*u10**3.5*exp(-1.5*log(dia/3e-6)**2)+6.8e-3*u10**3*exp(-3.3*log(dia/30e-6)**2)
    !
    ! Additional: temperature dependency from Jaegle et al. (2011)
    tw = sst-273.15 ! Temperature in C
    twt = 0.3+0.1*tw-0.0076*tw**2+0.00021*tw**3
    !
    ! Convert from dF/dDp to dF/dlog10Dp: multiply by Dp*ln(10)
    flux_smith = dia*1e6*log(10.)*flx*twt
  END FUNCTION flux_smith

  REAL FUNCTION flux_gong(dia,sst,u10)
    ! Flux parameterization from Gong (2003)
    ! Valid from 0.07e.6 to 20e-6 m, but reasonable values from about 0.01e-6 m up to 50e-6 m.
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL :: rad80, a, b, flx, tw, twt
    !
    ! Particle radius at 80 RH% is about the same as the dry diameter (assuming GF=2)
    rad80=dia*1e6 ! Radius in microns
    !
    ! Flux parameterization, dF/dr=dF/dDp, Eq. 2
    a = 4.7*(1.+30.*rad80)**(-0.017*rad80**(-1.44)) ! here theta=30
    b = (0.433-log10(rad80))/0.433
    flx = 1.373*u10**3.41/rad80**a*(1.+0.057*rad80**3.45)*10.**(1.607*exp(-b**2))
    !
    ! Additional: temperature dependency from Jaegle et al. (2011)
    tw = sst-273.15 ! Temperature in C
    twt = 0.3+0.1*tw-0.0076*tw**2+0.00021*tw**3
    !
    ! Convert from dF/dDp to dF/dlog10Dp: multiply by Dp*ln(10)
    flux_gong = rad80*log(10.)*flx*twt
  END FUNCTION flux_gong

  REAL FUNCTION flux_grythe(dia,sst,u10)
    ! The new flux parameterization from Grythe et al. (2014)
    ! Valid from 0.07e.6 to 10e-6 m, but reasonable values from about 0.01e-6 m up to 50e-6 m.
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL :: tw, twt, flx
    !
    ! Partial particle flux, dF/dDp, Eq. 7
    flx = 235.*u10**3.5*exp(-0.55*log(dia/0.1e-6)**2)+ &
          0.2*u10**3.5*exp(-1.5*log(dia/3e-6)**2) + &
          6.8e-3*u10**3*exp(-1.*log(dia/30e-6)**2)
    ! Note: the last coefficient is 6.8e-3 and not 6.8 (Smith and Harrison, 1998)
    !
    ! Temperature dependency, Eq. A7 (originally from Jaegle et al., 2011)
    tw = sst-273.15 ! Temperature in C
    twt=0.3+0.1*tw-0.0076*tw**2+0.00021*tw**3
    !
    ! Convert from dF/dDp to dF/dlog10(Dp): multiply by Dp*ln(10)
    flux_grythe = dia*1e6*log(10.)*flx*twt
  END FUNCTION flux_grythe

  REAL FUNCTION flux_sofiev(dia,sst,u10)
    ! Flux parameterization from Sofiev et al. (2011)
    ! Valid from 0.01e-6 to 10e-6 m, but can be extrapolated up to 100e-6 m
    REAL, INTENT(IN) :: dia,sst,u10 ! Dry diameter (m), sea surface temperature (K) and 10 m wind speed (m/s)
    REAL, PARAMETER :: sw=0.033 ! Salinity
    REAL :: diam, ai, bi, tw, ftw, fsw, w, flx
    !
    ! Diameter in micro meters
    diam=dia*1e6
    !
    ! Temperature
    tw = sst-273.15 ! Temperature in C
    if (tw<=-2.) then
        ai=0.092
        bi=-0.96
    elseif (tw<=5.) then
        ai=0.092+(0.15-0.092)/(5.+2.)*(tw+2.)
        bi=-0.96+(-0.88+0.96)/(5.+2.)*(tw+2.)
    elseif (tw<=15.) then
        ai=0.15+(0.48-0.15)/(15.-5.)*(tw-5.)
        bi=-0.88+(-0.36+0.88)/(15.-5.)*(tw-5.)
    elseif (tw<25.) then
        ai=0.48+(1.-0.48)/(25.-15.)*(tw-15.)
        bi=-0.36+(0.+0.36)/(25.-15.)*(tw-15.)
    else
        ai=1.
        bi=0.
    endif
    ftw=ai*diam**bi
    !
    ! Salinity
    if (sw<=0.) then
        ai=0.12
        bi=-0.71
    elseif (sw<=0.0092) then
        ai=0.12+(0.12-5.85e-5)/(0.0092-0.)*(sw-0.0092)
        bi=-0.71+(-0.71+1.7)/(0.0092-0.)*(sw-0.0092)
    else
        ai=0.12+(0.12-1.)/(0.0092-0.033)*(sw-0.0092)
        bi=-0.71+(-0.71+0.)/(0.0092-0.033)*(sw-0.0092)
    endif
    fsw=ai*diam**bi
    !
    ! Whitecap cover (% => fraction) based on 10 m wind speeds
    w=3.84e-6*u10**3.41
    !
    ! Partial particle flux, dF/dDp(0.033,25C), Eq. 6
    flx = 1e6*exp(-0.09/(diam+3e-3))/(2.+exp(-5./diam))*(1.+0.05*diam**1.05)/diam**3* &
        10.**(1.05*exp(-((0.27-log10(diam))/1.1)**2))
    !
    ! Convert from dF/dDp to dF/dlog10(Dp): multiply by Dp*ln(10)
    flux_sofiev=diam*log(10.)*flx*ftw*fsw*w
  END FUNCTION flux_sofiev

  ! --------------------------------------------------------------------------
  ! Gas emissions from ocean surface
  !
  SUBROUTINE marine_gas_flux(sst)
    ! Parameters
    use defs, only: vonk, g
    use grid, only: nxp, nyp, a_ustar, a_dn, zm, a_gaerot
    use mo_submctl, only: id_mtp, id_isop, mws_gas
    use stat, only: fill_scalar, sflg
    use util, only : get_avg2dh
    IMPLICIT NONE
    REAL, INTENT(IN) :: sst     ! Sea surface temperature (K)

    ! Local
    INTEGER :: i, j
    REAL :: usum, zs, u10_bar
    REAL :: u10(nxp,nyp)! whitecap cover, 10 m wind speed
    REAL :: flxIsop, flxMtrp ! Gas flux (kg/m2/s)

    real, parameter :: schmidt_ref = 660.0, & ! CO2 in 293 K
 & per_sec = 1.0 / 3600.0, to_m = 0.01, p23 = 2.0 / 3.0, p12 = 0.5
    real :: schmidt_isoprene, schmidt_monoterp, sc_ratio, temp_c
    integer :: iGas

    flxIsop = 0.
    flxMtrp = 0.
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
    ! Mean 10 m wind speed (for the domain)
    u10_bar = get_avg2dh(nxp,nyp,u10)

    temp_c = sst - 273.15

    ! Saltzman et al., 1993 (might need in the future)
    !schmidt_dms = 2674.0 - 147.12 * temp_c + 3.726 * temp_c**2 + 0.038*temp_c**3

    ! 1. Isoprene
    ! Palmer & Shaw, 2005
    schmidt_isoprene = 3913.15 - 162.13*temp_c + 2.67*temp_c**2 - 0.012*temp_c**3
    if(wtrIsop > 0. .and. id_isop>0)then
      sc_ratio = schmidt_ref / schmidt_isoprene
      iGas = id_isop

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
            flxIsop = 0.251*u10_bar*u10_bar*sc_ratio**p12 * per_sec*to_m
!          case default ! No flux
!            flxIsop = 0.
!          end select
      ! Flux
      flxIsop = flxIsop * wtrIsop * mws_gas(iGas)
      DO j=3,nyp-2
        DO i=3,nxp-2
          a_gaerot(2,i,j,iGas) = a_gaerot(2,i,j,iGas) + flxIsop / a_dn(2,i,j) / (zm(3)-zm(2))
        enddo
      enddo
    endif

    ! 2. Monoterpenes
    if(wtrMtrp > 0. .and. id_mtp>0)then
      iGas = id_mtp
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
      ! Transfer velocity (m/s)
      flxMtrp = 0.251*u10_bar*u10_bar*sc_ratio**p12 * per_sec*to_m
      ! Flux
      flxMtrp = flxMtrp * wtrMtrp * mws_gas(iGas)
      DO j=3,nyp-2
        DO i=3,nxp-2
          a_gaerot(2,i,j,iGas) = a_gaerot(2,i,j,iGas) + flxMtrp / a_dn(2,i,j) / (zm(3)-zm(2))
        enddo
      enddo
    endif

    if (sflg) then
        call fill_scalar(flxIsop,'flx_iso')
        call fill_scalar(flxMtrp,'flx_mt ')
        call fill_scalar(u10_bar,'u10    ')
    endif

  END SUBROUTINE marine_gas_flux

  ! --------------------------------------------------------------------------
  ! Constant surface fluxes (kg/m2/s) for active VOCs and VBS(g) and aqSOA(g) species
  !
  SUBROUTINE srfc_gas_flux()
    ! Parameters
    use grid, only: nxp, nyp, a_dn, zm, a_gaerot
    use mo_submctl, only: nvocs, nvbs, naqsoa, gas_srfc_flx
    IMPLICIT NONE
    INTEGER :: i, j, k

    DO k=1,nvocs+nvbs+naqsoa ! Gases always in this order
      IF (gas_srfc_flx(k)>0.) THEN
        DO j=3,nyp-2
          DO i=3,nxp-2
            ! Mass mixing ratio tendency (kg/kg/s) due to surface flux (kg/m2/s)
            a_gaerot(2,i,j,k) = a_gaerot(2,i,j,k) + gas_srfc_flx(k)/(a_dn(2,i,j)*(zm(3)-zm(2)))
          ENDDO
        ENDDO
      ENDIF
    ENDDO

  END SUBROUTINE srfc_gas_flux

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
         , wt_sfc, wq_sfc, obl, dn0, level,dtl, a_sflx, a_rflx, precip
    use thrm, only: rslf
    use stat, only: fill_scalar, sflg
    use util, only : get_avg2dh

    implicit none
    real, intent (inout) :: sst
    real :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp)       &
         ,wspd(nxp,nyp), bfct(nxp,nyp)
    real :: rx(nzp,nxp,nyp)

    real :: total_sw, total_rw, total_la, total_se, total_pre  ! Sami added
    real :: C_heat,lambda ! Sami added
    real :: K1,K2,K3,Kmean1,Kmean2,fii_1,fii_2,fii_3,Q3,Q12,Q23,ff1  ! Sami added

    integer :: i, j, iterate
    real    :: zs, bflx, ffact, sst1, bflx1, Vbulk, Vzt, bfg(2), usum

    ! Added by Juha
    SELECT CASE(level)
       CASE(0,1,2,3)
          rx = a_rv
       CASE(4,5)
          rx = a_rp
    END SELECT

    zs = -999. ! currently not set

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

       bfg(1) = get_avg2dh(nxp,nyp,a_theta(2,:,:))
       bfg(2) = get_avg2dh(nxp,nyp,rx(2,:,:)) ! Juha: rx

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

        ! From energy fluxes calculate new surface temperature
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
        ff1=1.0
        IF(W1<0.75) ff1=W1/0.75
        !
        !  Following is copied from case (2). No idea if this is valid or not..
        !
        call get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
        usum = 0.
        do j=3,nyp-2
           do i=3,nxp-2
              dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
              ! Flux of moisture is limited by water content.
              drdz(i,j) = rx(2,i,j) - ff1*rslf(psrf,min(sst,280.))
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

        sst = sst-(total_rw+total_la+total_se+((lambda*C_heat*7.27e-5/(2.0))**0.5*(sst-280.0)))&
                /(2.0e-2*C_heat+(lambda*C_heat/(2.0*7.27e-5))**0.5)*dtl

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
                bflx = g*wt_sfc(1,1)/th00 + g*ep2*wq_sfc(i,j)
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

    if (sflg) then
        ! Surface temperature (K)
        CALL fill_scalar(sst,'tsrf   ')
        ! Friction velocity
        usum = SUM(SUM(a_ustar(3:nxp-2,3:nyp-2),DIM=2))/float((nxp-4)*(nyp-4))
        CALL fill_scalar(usum,'ustar  ')
        ! Sensible heat flux
        usum = SUM(SUM(wt_sfc(3:nxp-2,3:nyp-2),DIM=2))/float((nxp-4)*(nyp-4))*cp*(dn0(1)+dn0(2))*0.5
        CALL fill_scalar(usum,'shf_bar')
        ! Latent heat flux
        usum = SUM(SUM(wq_sfc(3:nxp-2,3:nyp-2),DIM=2))/float((nxp-4)*(nyp-4))*alvl*(dn0(1)+dn0(2))*0.5
        CALL fill_scalar(usum,'lhf_bar')
        ! *** optional ts-outputs ***
        ! Obukhov lenght
        usum = SUM(SUM(obl(3:nxp-2,3:nyp-2),DIM=2))/float((nxp-4)*(nyp-4))
        CALL fill_scalar(usum,'obl    ')
        ! Surface pressure
        CALL fill_scalar(psrf,'psrf   ')
        ! Length scales
        CALL fill_scalar(zs,'z0m    ') ! Momentum zrough or the actually calculated value
        CALL fill_scalar(zrough_t,'z0t    ') ! Temperature roughness height
        ! Surface winds
        usum = SUM(SUM(wspd(3:nxp-2,3:nyp-2),DIM=2))/float((nxp-4)*(nyp-4))
        CALL fill_scalar(usum,'wspd   ')
    endif

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
    real    :: betg, z0t
    real    :: x, psi1, psi2, Lold, Ldif, zeff, zeta, lmo, dtv
    logical    :: exititer

    ! Optional temperature roughness height
    z0t=zrough_t
    IF (zrough_t<0.) z0t=z0

    betg  = th00/g

    do j=3,n3-2
       do i=3,n2-2
          dtv = dth(i,j) + ep2*th00*drt(i,j)

          !
          ! Neutral case
          !
          if (dtv == 0.) then
             ustar(i,j) = u(i,j)*vonk/log(z/z0)
             tstar(i,j) = dtv*vonk/log(z/z0t)/pr
             lmo        = -1.e10

          !
          ! start iterations from values at previous tstep,
          ! unless the sign has changed or if it is the first call, then
          ! use neutral values.
          !
          else
             if ((runtype=='INITIAL' .and. first_call) .or. (tstar(i,j)*dtv <= 0.)) then
                ustar(i,j) = u(i,j)*vonk/log(z/z0)
                tstar(i,j) = dtv*vonk/log(z/z0t)/pr
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
                tstar(i,j) = (dtv*vonk/pr)/(log(zeff/z0t) - psi2)

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



