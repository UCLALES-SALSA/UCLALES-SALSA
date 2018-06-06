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
MODULE fuliou

  USE defs,   ONLY : nv, nv1, mb, pi, totalpower, g, R, ep2
  USE cldwtr, ONLY : init_cldwtr, cloud_water, init_cldice, cloud_ice,  &
                     init_cldgrp, cloud_grp, init_aerorad, aero_rad  
  USE solver, ONLY : qft
  USE ckd
  USE mo_submctl, ONLY : nbins,spec

  IMPLICIT NONE

  LOGICAL, SAVE :: Initialized = .FALSE.
  INTEGER :: iseed = 0
  REAL, PARAMETER :: minSolarZenithCosForVis = 1.e-4

CONTAINS
  !
  !---------------------------------------------------------------------------
  ! Subroutine rad_init initialize data arrays for gases, ice model and water
  ! model on first call
  ! Juha: added initialization of aerosol model
  !
  SUBROUTINE rad_init()
    INTEGER, DIMENSION (:), ALLOCATABLE :: seed
    INTEGER :: isize

    IF (.NOT. Initialized) THEN
       ! Initialize random numbers for McICA
       CALL random_seed(size=isize)
       ALLOCATE (seed(isize))
       seed(:) = iseed
       CALL random_seed(put=seed)
       DEALLOCATE (seed)

       CALL init_ckd
       CALL init_cldwtr
       CALL init_cldice
       CALL init_cldgrp
       CALL init_aerorad
       Initialized = .TRUE.
    END IF

  END SUBROUTINE rad_init
  !
  ! ----------------------------------------------------------------------
  ! Subroutine set_random_offset is needed to generate true pseudorandom numbers for parallel runs
  SUBROUTINE set_random_offset(ioffset)
    IMPLICIT NONE
    INTEGER :: ioffset, i
    REAL    :: randomNumber

    ! Initialize random number generator, if not yet initialized
    IF (.NOT. Initialized) CALL rad_init()

    ! Call randon mumbers
    IF (ioffset > 0) THEN
        DO i = 1, ioffset
            CALL random_number(randomNumber)
        END DO
    END IF

  END SUBROUTINE set_random_offset
  !
  ! ----------------------------------------------------------------------
  ! Subroutine rad: Computes radiative fluxes using a band structure 
  ! defined by input ckd file
  !
  SUBROUTINE rad (as, u0, ss, pts, ee, pp, pt, ph, po, fds, fus, fdir, fuir, &
                  McICA, nspec, plwc, pre, piwc, pde, prwc, pgwc, maerobin, naerobin)

    INTEGER, INTENT(in) :: nspec
    REAL, INTENT(in) :: pp (nv1) ! pressure at interfaces

    REAL, DIMENSION(nv), INTENT (in) :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    REAL, OPTIONAL, DIMENSION(nv), INTENT (in) :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre,  & ! effective radius of cloud droplets [microns]
         piwc, & ! cloud ice water content [g/m^3]
         pde,  & ! effective diameter of ice particles [microns]
         prwc, & ! rain water content [g/m^3]
         pgwc    ! graupel water content

    REAL, OPTIONAL, INTENT(in) :: maerobin(nv,nspec*nbins),  & !  maerobin(:,:), naerobin(:,:)
                                  naerobin(nv,nbins)     !

    REAL, INTENT(in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         ee, & ! broadband surface emissivity (all IR bands given this value)
         u0, & ! cosine of solar zenith angle
         ss, & ! Solar constant
         pts   ! Surface skin temperature

    LOGICAL, INTENT(in) :: McICA
    
    REAL, DIMENSION(nv1), INTENT (out) ::  &
         fds, fus,  & ! downward and upward solar flux
         fdir, fuir   ! downward and upward ir flux

    CALL rad_ir(nspec,pts, ee, pp, pt, ph, po, fdir, fuir, McICA, &
                 plwc, pre, piwc, pde, prwc, pgwc, maerobin, naerobin)

    CALL rad_vis(nspec,as, u0, ss, pp, pt, ph, po, fds, fus, McICA, &
                 plwc, pre, piwc, pde, prwc, pgwc, maerobin, naerobin)

  END SUBROUTINE rad

  ! ----------------------------------------------------------------------
  ! Subroutine rad_ir 
  ! Computes IR radiative fluxes using a band structure 
  ! defined by input ckd file
  !
  SUBROUTINE rad_ir (nspec,pts, ee, pp, pt, ph, po, fdir, fuir, McICA, &
                     plwc, pre, piwc, pde, prwc, pgwc, maerobin, naerobin)

    INTEGER, INTENT(in) :: nspec
    REAL, INTENT(in) :: pp (nv1) ! pressure at interfaces

    REAL, DIMENSION(nv), INTENT(in) :: &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    REAL, OPTIONAL, DIMENSION(nv), INTENT(in) :: &
         plwc, & ! cloud liquid water content [g/m^3]
         pre,  & ! effective radius of cloud droplets [microns]
         piwc, & ! cloud ice water content [g/m^3]
         pde,  & ! effective diameter of ice particles [microns]
         prwc, & ! rain water content [g/m^3]
         pgwc    ! graupel water content

    REAL, OPTIONAL, INTENT(in) :: maerobin(nv,nspec*nbins), naerobin(nv,nbins) !maerobin(:,:), naerobin(:,:)
    REAL, INTENT(in) :: &
         ee, & ! broadband surface emissivity (all IR bands given this value)
         pts   ! Surface skin temperature
         
    LOGICAL, INTENT(in) :: McICA

    REAL, DIMENSION(nv1), INTENT (out) :: &
         fdir, fuir   ! downward and upward ir flux

    ! ----------------------------------------
    LOGICAL, PARAMETER :: irWeighted = .FALSE. 

    REAL, DIMENSION (nv)   :: tw,ww,tg,dz,tauNoGas, wNoGas, Tau, w
    REAL, DIMENSION (nv)   :: taer,waer
    REAL, DIMENSION (nv)   :: ti,wi,tgr,wgr
    REAL, DIMENSION (nv1)  :: fu1, fd1, bf
    REAL, DIMENSION (nv,4) :: www, pfNoGas, pf
    REAL, DIMENSION (nv,4) :: wwi,wwgr
    REAL, DIMENSION (nv,4) :: wwaer

    INTEGER :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    REAL    :: fuq2, xir_norm
    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bandWeights
    REAL    :: randomNumber
    ! ----------------------------------------

    IF (.NOT. Initialized) CALL rad_init()

    IF(.NOT. allocated(bandweights)) THEN 
      ALLOCATE(bandweights(size(ir_bands)))
      CALL computeIRBandWeights(ir_bands, irWeighted, bandWeights)
    END IF
    
    fdir(:) = 0.0; fuir(:) = 0.0

    CALL thicks(pp, pt, ph, dz) 

    IF (McICA) THEN
       !
       ! Select a single band and g-point (ib, ig1) and use these as the limits
       !   in the loop through the spectrum below. 
       !
       CALL random_number(randomNumber)
       CALL select_bandg(ir_bands, bandweights, randomNumber, ib, ig1) 
       ig2 = ig1
       iblimit = 1
    ELSE
       iblimit = size(ir_bands)
    END IF

    bandLoop: DO ibandloop = 1, iblimit
      IF (.NOT. McICA) THEN
         ib  = ibandloop
         ig1 = 1
         ig2 = kg(ir_bands(ib))
      END IF
      !WRITE(*,*) 'HEP IR', 1./center(band(ib+size(solar_bands)))
      !
      ! Water vapor continuum optical depth
      !
      CALL gascon ( center(ir_bands(ib)), pp, pt, ph, TauNoGas )
      wNoGas = 0.; pfNoGas  = 0. 
      IF (present(plwc)) THEN
        CALL cloud_water(ib + size(solar_bands), pre, plwc, dz, tw, ww, www)
        CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw, ww, www)
      END IF
      IF (present(piwc)) THEN
        CALL cloud_ice(ib + size(solar_bands), pde, piwc, dz, ti, wi, wwi)
        CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, ti, wi, wwi)
      END IF
      IF (present(pgwc)) THEN
        CALL cloud_grp(ib + size(solar_bands), pgwc, dz, tgr, wgr, wwgr)
        CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgr, wgr, wwgr)
      END IF 
      IF ( PRESENT(maerobin) .AND. PRESENT(naerobin) ) THEN
         CALL aero_rad(ib + size(solar_bands), nbins, nspec, maerobin, naerobin, &
                       dz, taer, waer, wwaer)
         CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, taer, waer, wwaer)
      END IF

      CALL planck(pt, pts, llimit(ir_bands(ib)), rlimit(ir_bands(ib)), bf)

      gPointLoop: DO ig = ig1, ig2
         tau = TauNoGas; w = wNoGas; pf = pfNoGas
         CALL gases (ir_bands(ib), ig, pp, pt, ph, po, tg )
         CALL combineOpticalProperties(tau, w, pf, tg)

         !
         ! Solver expects cumulative optical depth
         !
         DO k = 2, nv
           tau(k) = tau(k) + tau(k - 1)
         END DO
         CALL qft (.FALSE., ee, 0., 0., bf, tau, w, pf(:, 1), pf(:, 2),      &
                   pf(:, 3), pf(:, 4), fu1, fd1)

         IF (McICA) THEN 
            xir_norm = 1./bandweights(ib)
         ELSE
            xir_norm = gPointWeight(ir_bands(ib), ig)
         END IF

         fdir(:) = fdir(:) + fd1(:) * xir_norm
         fuir(:) = fuir(:) + fu1(:) * xir_norm
      END DO gPointLoop
    END DO bandLoop
    !
    ! fuq2 is the surface emitted flux in the band 0 - 280 cm**-1 with a
    ! hk of 0.03.
    !
    fuq2 = bf(nv1) * 0.03 * pi * ee
    fuir(:) = fuir(:) + fuq2
  END SUBROUTINE rad_ir
  ! ----------------------------------------------------------------------
  ! Subroutine rad_vis: Computes radiative fluxes using a band structure 
  ! defined by input ckd file
  !

  SUBROUTINE rad_vis (nspec,as, u0, ss, pp, pt, ph, po, fds, fus, McICA,  &
                      plwc, pre, piwc, pde, prwc, pgwc, maerobin, naerobin )

    INTEGER, INTENT(in) :: nspec
    REAL, INTENT(in) :: pp (nv1) ! pressure at interfaces

    REAL, DIMENSION(nv), INTENT (in) ::  &
         pt,   & ! temperature [K] at mid points
         ph,   & ! humidity mixing ratio in kg/kg
         po      ! ozone mixing ratio

    REAL, OPTIONAL, DIMENSION(nv), INTENT (in) ::  &
         plwc, & ! cloud liquid water content [g/m^3]
         pre,  & ! effective radius of cloud droplets [microns]
         piwc, & ! cloud ice water content [g/m^3]
         pde,  & ! effective diameter of ice particles [microns]
         prwc, & ! rain water content [g/m^3]
         pgwc    ! graupel water content

    REAL, OPTIONAL, INTENT(in) :: maerobin(nv,nspec*nbins), naerobin(nv,nbins) !maerobin(:,:), naerobin(:,:)

    REAL, INTENT(in) :: &
         as, & ! broadband albedo (all visible bands given this value)
         u0, & ! cosine of solar zenith angle
         ss    ! Solar constant
         
    LOGICAL, INTENT(in) :: McICA

    REAL, DIMENSION(nv1), INTENT (out)::  &
         fds, fus    ! downward and upward solar flux

    ! ----------------------------------------
    LOGICAL, PARAMETER :: solarWeighted = .FALSE. ! Could be .TRUE.?

    REAL, DIMENSION(nv)   :: tw,ww,tg,tgm,dz, tauNoGas, wNoGas, tau, w
    REAL, DIMENSION(nv)   :: ti,wi
    REAL, DIMENSION (nv)  :: taer,waer
    REAL, DIMENSION(nv)   :: tgr,wgr
    REAL, DIMENSION(nv1)  :: fu1, fd1, bf
    REAL, DIMENSION(nv,4) :: www, pfNoGas, pf
    REAL, DIMENSION(nv,4) :: wwi
    REAL, DIMENSION(nv,4) :: wwgr
    REAL, DIMENSION (nv,4):: wwaer

    REAL, DIMENSION(:), ALLOCATABLE, SAVE :: bandWeights

    INTEGER :: ib, ig, k, ig1, ig2, ibandloop, iblimit
    REAL    :: fuq1, xs_norm
    REAL    :: randomNumber
    ! ----------------------------------------

    IF (.NOT.Initialized) CALL rad_init()

    IF (.NOT. allocated(bandweights)) THEN 
      ALLOCATE(bandweights(size(solar_bands)))
      CALL computeSolarBandWeights(solar_bands, solarWeighted, bandWeights)
    END IF
    
    fds(:) = 0.0
    fus(:) = 0.0
    bf(:)  = 0.0
    
    !WRITE(*,*) 'HEP'
    !WRITE(*,*) '1', 1./center(solar_bands)
    !WRITE(*,*) '2',1./center(ir_bands)
    !WRITE(*,*) '3',1./center(band)


    IF(u0 > minSolarZenithCosForVis) THEN
      CALL thicks(pp, pt, ph, dz) 
  
      IF (McICA) THEN
         !
         ! Select a single band and g-point (ib, ig1) and use these as the
         ! limits in the loop through the spectrum below. 
         !
         CALL random_number(randomNumber)
         CALL select_bandg(solar_bands, bandweights, randomNumber, ib, ig1) 
         ig2 = ig1
         iblimit = 1 
      ELSE
         iblimit = size(solar_bands)
      END IF
  
      bandLoop: DO ibandloop = 1, iblimit
         !
         ! Select g points either all, or one depending on McICA
         !
         IF (.NOT. McICA) THEN
           ib  = ibandloop
           ig1 = 1
           ig2 = kg(solar_bands(ib))
         END IF
  
         !WRITE(*,*) 'HEP SW', 1./center(band(ib))

         !
         ! Rayleigh scattering
         !
         CALL rayle ( ib, u0, power(solar_bands(ib)), pp, pt, dz, tauNoGas, &
                      wNoGas, pfNoGas)
         !
         ! Water vapor continuum
         !
         CALL gascon ( center(solar_bands(ib)), pp, pt, ph, tgm )
         IF(any(tgm > 0.)) &
           CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgm)
         !
         ! Cloud water
         !
         IF (present(plwc)) THEN
           CALL cloud_water(ib, pre, plwc, dz, tw, ww, www)
           CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tw,ww,www)
         END IF
         IF (present(piwc)) THEN
           CALL cloud_ice(ib, pde, piwc, dz, ti, wi, wwi)
           CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, ti,wi,wwi)
         END IF
         IF (present(pgwc)) THEN
           CALL cloud_grp(ib,pgwc, dz, tgr, wgr, wwgr)
           CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, tgr, wgr,wwgr)
         END IF 
         IF ( PRESENT(maerobin) .AND. PRESENT(naerobin) ) THEN
            CALL aero_rad(ib, nbins, nspec, maerobin, naerobin, dz, taer, waer, wwaer)
            CALL combineOpticalProperties(TauNoGas, wNoGas, pfNoGas, taer, waer, wwaer)
         END IF
 
         gPointLoop: DO ig = ig1, ig2
            tau = TauNoGas; w = wNoGas; pf = pfNoGas
            CALL gases (solar_bands(ib), ig, pp, pt, ph, po, tg )
            CALL combineOpticalProperties(tau, w, pf, tg)
            
            !WRITE(*,*) SUM(taer(:))
            !
            ! Solver expects cumulative optical depth
            !
            DO k = 2, nv
               tau(k) = tau(k) + tau(k - 1)
               !WRITE(*,*) 'VIS',k,taer(k),waer(k),wwaer(k,:)
            END DO
            CALL qft (.TRUE., 0., as, u0, bf, tau, w, pf(:, 1), pf(:, 2),    &
                      pf(:, 3), pf(:, 4), fu1, fd1)
            IF (McICA) THEN 
               xs_norm = power(solar_bands(ib))/ bandweights(ib)
            ELSE
               xs_norm = gPointWeight(solar_bands(ib), ig)*power(solar_bands(ib))
            END IF
            fds(:) = fds(:) + fd1(:) * xs_norm
            fus(:) = fus(:) + fu1(:) * xs_norm
         END DO gPointLoop
      END DO bandLoop
      !
      ! In this model, we used the solar spectral irradiance determined by
      ! Thekaekara (1973), and 1340.0 W/m**2 is the solar energy contained 
      ! in the spectral region 0.2 - 4.0 um., thus scale solar fluxes by
      ! fuq1
      !
      !WRITE(*,*) "EKA", fds(90), fus(90), tau(90)

      fuq1 = ss / totalpower
      fds(:)  = fds(:)*fuq1
      fus(:)  = fus(:)*fuq1
      !WRITE(*,*) "TOKA", fds(90), fus(90), tau(90)
    END IF 
  END SUBROUTINE rad_vis
  ! ----------------------------------------------------------------------
  ! Subroutine select_bandg
  !
  ! Selects the band (i) and the g point (j) based on the probability a
  ! photon would be found in the wavelengths covered by the band and g
  ! point, which is given by g_prob.  Note g_prob sums to unity for both
  ! the solar bands (power > 0.) and the infrared bands respectively.
  ! 
  SUBROUTINE select_bandg(bands, bandweights, randomNumber, i, j)
    TYPE(band_properties), &
          DIMENSION(:),    &
             INTENT(in)  :: bands
    REAL, DIMENSION(:), &
             INTENT(in)  :: bandweights
    REAL,    INTENT(in)  :: randomNumber
    INTEGER, INTENT(out) :: i, j

    REAL :: cumulative
    
    i = 1; j = 1
    ! The probability contained in the first g point of the first band
    cumulative = gPointWeight(bands(i), j) * bandweights(i)

    DO WHILE (randomNumber > cumulative .AND. cumulative < 1.0)
       j = j+1
       IF (j > kg(bands(i)) ) THEN
          i = i+1
          j = 1
       END IF
       cumulative = cumulative + gPointWeight(bands(i), j) * bandweights(i)
    END DO
  END SUBROUTINE select_bandg

  ! ----------------------------------------------------------------------
  ! Subroutine thicks: Integrates the hydrostatic equation to provide 
  ! layer thicknesses
  ! 
  SUBROUTINE thicks(pp, pt, ph, dz) 

    REAL, INTENT (in)  :: pp(nv1), pt(nv), ph(nv)
    REAL, INTENT (out) :: dz(nv)

    INTEGER :: i
    REAL    :: tv

    DO  i = 1, nv
       tv = pt(i)*(1+0. + ep2*ph(i) )
       dz(i) = (R/g) * tv * alog( pp(i+1) / pp(i) )
    END DO
    
  END SUBROUTINE thicks
  ! ----------------------------------------------------------------------
  !
  SUBROUTINE combineOpticalProperties(tau,      ssa,      pF, &
                                      tauToAdd, ssaToAdd, pFtoAdd)
    REAL, DIMENSION(:),    INTENT(inout) :: tau, ssa
    REAL, DIMENSION(:, :), INTENT(inout) :: pF   ! Phs function (level, moment)
    REAL, DIMENSION(:),    INTENT(in)    :: tautoAdd
    REAL, DIMENSION(:),    OPTIONAL, INTENT(in) :: ssaToAdd
    REAL, DIMENSION(:, :), OPTIONAL, INTENT(in) :: pFToAdd ! Phs function

    INTEGER :: j
    !
    ! Adds optical properties to running sum
    !   If ssa and/or w[1-4] are not present we assume the new medium is 
    ! strictly absorbring
    ! 
    IF(present(ssaToAdd) .AND. present(pfToAdd)) THEN
       DO j = 1, size(pF, 2) 
          WHERE (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:) > 0.)
             pf(:, j) = (ssa(:)*tau(:)*pf(:, j) + ssaToAdd(:)*tauToAdd(:)     &
                        * pfToAdd(:, j))/(ssa(:)*tau(:) + ssaToAdd(:) * tauToAdd(:))
          ELSE WHERE
             pf(:, j) = 0.  
          END WHERE
       END DO
       WHERE (tau(:) + tauToAdd(:) > 0.)
          ssa(:) = (ssa(:) * tau(:) + ssaToAdd(:) * tauToAdd(:)) /            &
                   (tau(:) + tauToAdd(:))
       ELSE WHERE
          ssa(:) = 0. 
       END WHERE
       tau(:) = tau(:) + tauToAdd(:) 
    ELSE
      !
      ! New medium is absorbing - phase function doesn't change
      !
       ssa(:) = (ssa(:) * tau(:)) / (tau(:) + tauToAdd(:))
       tau(:) = tau(:) + tauToAdd(:) 
    END IF
 
  END SUBROUTINE combineOpticalProperties
  ! ----------------------------------------------------------------------
  ! Subroutine rayle:  computes optical properties associated with rayleigh
  ! scattering 
  !
  ! ri is the coefficient in Eq.(4.8) of Fu (1991) to compute the optical
  ! depth due to Rayleigh scattering in the solar bands.
  !
  ! tr, wr, and wwr are the optical depth, single scattering albedo,
  ! and expansion coefficients of the phase function ( 1, 2, 3, and
  ! 4 ) due to the Rayleigh scattering for a given layer.
  ! 
  SUBROUTINE rayle ( ib, u0, power, pp, pt, dz, tr, wr, wwr)
    INTEGER, INTENT(in) :: ib
    REAL, INTENT(in)    :: u0, power, pp(nv1), pt(nv), dz(nv)
    REAL, INTENT(out)   :: tr(nv), wr(nv), wwr(nv,4)

    REAL, PARAMETER :: ri(6) = (/ 0.9022e-5, 0.5282e-6, 0.5722e-7, &
                                  0.1433e-7, 0.4526e-8, 0.1529e-8 /)

    INTEGER :: i
    REAL    :: x

    IF ( ib == 1 ) THEN
       x = -3.902860e-6*u0*u0 + 6.120070e-6*u0 + 4.177440e-6
    ELSE
       x = ri(ib)
    END IF
    
    IF(power > 0.) THEN
      DO  i = 1, nv
        tr(i) = x * ( pp(i) + pp(i+1) ) * dz(i) * 0.5 / pt(i)
      END DO
      wr(:) = 1.0
      wwr(:, :) = 0. 
      wwr(:, 2) = 0.5
    ELSE
      tr(:)     = 0. 
      wr(:)     = 0.
      wwr(:, :) = 0. 
    END IF 
  END SUBROUTINE rayle

  ! *********************************************************************
  ! tgm(nv) are the optical depthes due to water vapor continuum absorp-
  ! tion in nv layers for a given band ib. We include continuum absorp-
  ! tion in the 280 to 1250 cm**-1 region. vv(11)-vv(17) are the central
  ! wavenumbers of each band in this region. 
  ! *********************************************************************
  SUBROUTINE gascon ( center, pp, pt, ph, tgm)
    REAL, INTENT(in)  :: center, pp(nv1), pt(nv), ph(nv)
    REAL, INTENT(out) :: tgm(nv)

    INTEGER :: k
    REAL    :: ff, pe, s, pmid

    IF ( center >= 280 .AND. center <= 1250.) THEN
       s = ( 4.18 +  5577.8 * exp ( - 0.00787 * center ) ) / 1013.25
       DO k = 1, nv
          pmid   = (pp(k) + pp(k+1))*0.5
          pe     = pmid * ph(k) / ( 0.622 + 0.378 * ph(k) )
          ff     = s*(pe +  0.002*pmid ) *exp (1800.0/pt(k) -6.08108)
          tgm(k) = ff * ph(k) * ( pp(k+1) - pp(k) ) * 1.019767
       END DO
    ELSE
       tgm(:) = 0.
    END IF
  END SUBROUTINE gascon
  ! ----------------------------------------------------------------------
  ! Subroutine planck:  Integrates planck function over band in xk sub
  ! intervals.  The temperatures at the interfaces are taken as the 
  ! average of the mid-point temperatures, the temperature at the lowest
  ! pressure interface is set to the temperature at the first mid point, 
  ! and surface temperatures are taken as the skin temperature.
  ! 
  SUBROUTINE planck ( pt, tskin, llimit, rlimit, bf)
    REAL, INTENT(in)  :: pt(nv), tskin, llimit, rlimit
    REAL, INTENT(out) :: bf(nv1) ! intensity [W/m^2/Sr]

    REAL, PARAMETER :: xk = 10.

    INTEGER :: k
    REAL    :: v1, v2, vmid, fq1, fq2, tk

    DO k = 1, nv1
       bf(k) = 0.0
    END DO

    v1 = llimit
    DO WHILE (v1 > rlimit+epsilon(rlimit))
       v2 = max(v1 - xk, rlimit)
       vmid = ( v1 + v2 ) * 0.5
       fq1 = 1.19107e-8 * vmid * vmid * vmid
       fq2 = 1.43884 * vmid
       DO k = 2, nv
          tk = (pt(k)+pt(k-1))*0.5
          IF (tk <= 0.) THEN
             PRINT*,'tk wrong',tk,v1,v2,rlimit,llimit
             stop
          END IF
          bf(k) = bf(k) + (fq1/(exp(fq2/tk) - 1.0))*(v1-v2)
       END DO
       IF (pt(1) <= 0.) THEN
          PRINT*,'pt wrong',pt(1),v1,v2,rlimit,llimit
          stop
       END IF
       bf(1) = bf(1) + (fq1/(exp(fq2/pt(1)) - 1.0))*(v1-v2)
       bf(nv1) = bf(nv1) + (fq1/(exp(fq2/tskin) - 1.0))*(v1-v2)
       v1 = v2
    END DO

  END SUBROUTINE planck
  !
  ! ---------------------------------------------------------------------------
  SUBROUTINE computeIRBandWeights(bands, weighted, bandweights)
    TYPE(band_properties), &
          DIMENSION(:), INTENT(in)  :: bands
    LOGICAL,            INTENT(in)  :: weighted
    REAL, DIMENSION(:), INTENT(out) :: bandweights
    
    INTEGER :: ib
    !
    ! find the weighting for band points so that the probability of a photon
    ! existing in the g-point range of a band can be calculated, and used for
    ! McICA calculations.  This is the relative band width for the IR bands
    !
    
    IF(size(bands) /= size(bandweights)) &
      STOP "Didn't provide the right amount of storage for band weights" 
    IF(any(isSolar(bands))) STOP "Can't compute IR band weights for solar bands." 
    
    IF (weighted) THEN
       DO ib = 1, size(bands)
         bandweights(ib) = (llimit(bands(ib)) - rlimit(bands(ib)))/(bllmx-brlmn)
       END DO 
    ELSE
       bandweights(:) = 1./(REAL(size(bands)))
    END IF
  END SUBROUTINE computeIRBandWeights
  ! ---------------------------------------------------------------------------
  SUBROUTINE computeSolarBandWeights(bands, weighted, bandweights)
    TYPE(band_properties), &
          DIMENSION(:), INTENT(in)  :: bands
    LOGICAL,            INTENT(in)  :: weighted
    REAL, DIMENSION(:), INTENT(out) :: bandweights
    
    INTEGER :: i
    !
    ! find the weighting for band points so that the probability of a photon
    ! existing in the g-point range of a band can be calculated, and used for
    ! McICA calculations.  This is the solar power for the solar bands
    !
        IF(size(bands) /= size(bandweights)) &
         STOP "Didn't provide the right amount of storage for band weights" 
    IF(any(.NOT. isSolar(bands))) STOP "Can't compute solar band weights in IR"

    IF(weighted) THEN 
       bandweights(:) = (/ (power(bands(i))/totalpower, i = 1, size(bands)) /)
    ELSE 
       bandweights(:) = 1./(REAL(size(bands)))
    END IF
  END SUBROUTINE computeSolarBandWeights

END MODULE fuliou
