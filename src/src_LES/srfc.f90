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
MODULE srfc
  USE grid, ONLY : nzp,nxp,nyp,dtl,deltaz,level, W1,W2,W3,cmbcnst,sst,psrf, &
                   deltax,deltay
  USE mo_structured_datatypes
  
  USE mo_salsa_types, ONLY : aero
  USE mo_salsa_sizedist, ONLY : size_distribution
  USE mo_progn_state, ONLY : a_maerot, a_naerot, a_naerop, a_indefp, a_indeft  
  USE util, ONLY: getMassIndex
  USE mo_submctl, ONLY : pi6, in1a, fn2a, in2b, fn2b, nbins, nliquid, spec, prlim, &
                         ice_theta_dist, ica, fca, icb, fcb, ncld, nprc, ice_theta_dist
  USE mpi_interface, ONLY : myid
  USE mo_diag_state, ONLY: a_tskin, a_qskin, a_fgi, a_weight, a_fcz0, a_fuelmcg, &
  			   a_ignitiontime, a_fuelburnt, a_firespread, a_areaburnt, &
  			   a_phiwc, a_phiwb, a_tcrit, a_R0
  USE ncio, ONLY : open_surf_nc, read_surf_nc_2d, close_nc
  USE thrm, ONLY: rslf
  USE emission_init, ONLY: regime_limits
  USE mo_mpi_io, ONLY: write_hist_field_2d
  
  
  IMPLICIT NONE

  ! NAMELIST parameters
  ! ------------------------------------------------
  INTEGER :: isfctyp = 0       ! Surface model type
  REAL    :: zrough =  0.1e-3  ! Roughness length (meters, check!)
  REAL    :: ubmin  =  0.20
  
  ! Specified surface fluxes for isfctyp=0/default. These are used for other purposes isfctyp > 0.
  ! See comments for subroutine SURFACE.
  REAL :: dthcon = 100.0  ! Sensible heat flux Wm-2
  REAL :: drtcon = 0.0    ! Latent heat flux Wm-2
  ! -----------------------------------------------

  
  ! Other module variables
  ! -----------------------------------------------------------
  ! Additional perturbation values for surface fluxes with isfctyp=0/default
  REAL, ALLOCATABLE :: lh_flx(:,:)  ! Latent heat Wm-2
  REAL, ALLOCATABLE :: sh_flx(:,:)  ! Sensible heat Wm-2
  
  ! <--- Sami added
  ! Sami added ----->
  ! for isfctyp == 5:
  ! Initial values for surface properties
  ! REAL :: W1 = 0.9   !Water content      ... Definition in grid now because of restrat files
  ! REAL :: W2 = 0.9
  ! REAL :: W3 = 0.9
  REAL ::  B1 = 6.5
  REAL ::  B3 = 7.26
  REAL ::  K_s1 = 0.89e-5  ! 0-30 cm Hydraulic constant at saturation App. Acs-1991 m/s
  REAL ::  K_s3 = 0.93e-5  ! 30-60 cm Hydraulic constant at saturation App. Acs-1991
  REAL ::  fii_s1 = -0.036 ! 0-30 cm soil moisture potential at saturation Table 2 Acs-1991 m
  REAL ::  fii_s3 = -0.085 ! 30-60 cm soil moisture potential at saturation Table 2 Acs-1991 m
  REAL ::  thetaS1 = 0.5490! volumetric soil moisture at saturation App. Acs-1991 m3/m3
  REAL ::  thetaS2 = 0.5490!  
  REAL ::  thetaS3 = 0.488  
  REAL ::  D1 = 0.1        ! thickness of the soil layers in m Acs-1991 Fig1
  REAL ::  D2 = 0.3
  REAL ::  D3 = 0.6
  ! Juha moved/added
  ! for isfctyp == 5:
  REAL :: C_heat = 2.e6                 ! Surface heat capacity
  REAL :: deepSoilTemp = 280.            ! Assumed deep soil layer temperature
  LOGICAL :: lConstSoilWater = .FALSE.   ! Keep the value(s) of surface water content constant (as specified in NAMELIST)
  LOGICAL :: lConstSoilHeatCap = .FALSE. ! Keep the value of surface heat capacity constant (as specified in NAMELIST)
 
  PUBLIC 
  REAL, ALLOCATABLE, SAVE :: firespreadg(:,:), ignitiontimeg(:,:), areaburntg(:,:)
  REAL, ALLOCATABLE, SAVE :: ignitiontimecell(:,:,:), areaignitedcell(:,:,:)
 
  
CONTAINS
  
  ! Called in init.f90
  SUBROUTINE surface_initialize()
    USE grid, ONLY : sst,psrf
    
    IMPLICIT NONE
    
    ALLOCATE(lh_flx(nxp,nyp),sh_flx(nxp,nyp))
    lh_flx = 0.
    sh_flx = 0. 
    
    WRITE(*,*) 'Initial', nxp,nyp
       
    IF (isfctyp==6) CALL surface_state()
            
    
  END SUBROUTINE surface_initialize
  
   ! --------------------------------------------------------------------------
   ! SURFACE: Calculates surface fluxes using an algorithm chosen by ISFCLYR
   ! and fills the appropriate 2D arrays
   !
   !   default: specified thermo-fluxes (drtcon, dthcon)
   !   isfctyp=1: specified surface layer gradients (drtcon, dthcon)
   !   isfctyp=2: fixed lower boundary of water at certain sst
   !   isfctyp=3: bulk aerodynamic law with coefficients (drtcon, dthcon)
   !   isfctyp=4: regulate surface temperature to yield a constant surface buoyancy flux
   !
   ! Modified for level 4: a_rv replaced by a local variable rx, which has
   ! values a_rv if level < 4, and a_rp if level == 4 (i.e. water vapour mixrat in both cases)
   !
   ! Juha Tonttila, FMI, 2014
   ! 
   !   isfctyp=5: calculate surface fluxes from energy budget based on Acs et al: 
   !              A coupled soil moisture and surface
   !              temperature prediction model, Journal of applied meteorology, 30, 1991
   !              Sami Romakkaniemi, FMI
   !   isfctyp=6: calculate surface fluxes as in default but add fluxes from 
   !              vegetation fire. Silvia Calderón FMI 2025
   
   !
   SUBROUTINE surface(time, dtl)

     USE defs, ONLY: vonk, p00, rcp, g, cp, alvl, ep2,ep, rowt, ep, omega, pi
     USE mo_diag_state, ONLY : a_theta, a_rv, a_ustar, a_tstar, a_rstar,  &
                               uw_sfc, vw_sfc, ww_sfc, wt_sfc, wq_sfc,    &
                               a_sflx, a_rflx, a_rrate, a_temp, & 
                               a_tskin, a_qskin, &
                               a_fgi, a_weight, a_fcz0, a_fuelmcg, & 
                               a_ignitiontime, a_fuelburnt, a_firespread, &
                               a_areaburnt, a_phiwc, a_phiwb, a_tcrit, a_R0                    
     USE mo_aux_state, ONLY : zt, dn0, xt, yt
     USE mo_progn_state, ONLY : a_rp
     USE mo_vector_state, ONLY : a_up, a_vp
     USE grid, ONLY: th00, umean, vmean, &
                     level, psrf, sst, W1,W2,W3,&
                      cmbcnst!,  &
                     !mc_ApVdom, deltaz
      !USE stat, ONLY: sfc_stat, sflg, mcflg, acc_massbudged
     USE mpi_interface, ONLY : nypg, nxpg, double_array_par_sum, & 
                               cyclics2d, cyclicc2d
      

      IMPLICIT NONE
      REAL :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp),       &
              wspd(nxp,nyp), bfct(nxp,nyp)
      REAL :: rx(nzp,nxp,nyp)
      
      REAL :: total_rw, total_la, total_se, total_pre !, total_sw  ! Sami added
      REAL :: lambda ! Sami added
      REAL :: K1,K2,K3,Kmean1,Kmean2,fii_1,fii_2,fii_3,Q3,Q12,Q23,ff1  ! Sami added
      
     
      INTEGER :: i, j, iterate, req(8)
      REAL    :: zs, bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
      REAL (kind=8) :: bfl(2), bfg(2)
      

      REAL :: mctmp(nxp,nyp) ! Helper for mass conservation statistics    
      
      REAL :: tt                  ! time in seconds after ignition time
      REAL,  INTENT (in) :: time  ! time in seconds (since model start)
      REAL,  INTENT (in) :: dtl   ! dtlt seconds since previous timestep
      REAL :: fraction_burnt      ! Fraction of fuel(i.e. vegetation) burnt (time, time+dtl)
      REAL :: total_la_fire       ! latent heat coming from moisture evaporation 
      REAL :: total_se_fire       ! sensible heat coming from combustion
      REAL :: total_aer_fire      ! mass of aerosol produced by fire
      REAL :: area_burnt_cell     ! area burnt between time and time+dtl
      REAL :: fuel_burnt          ! fuel burnt per time step in a cell
      REAL :: fbcell              ! fraction of the cell burning
      INTEGER :: kk
            
      mctmp = 0.
      fraction_burnt = 0.
      tt = 0.
      fuel_burnt = 0.
      fbcell= 0.
      
      
      ! Added by Juha
      SELECT CASE(level)
         CASE(1,2,3)
            rx = a_rv%d
         CASE(4,5)
            rx = a_rp%d
      END SELECT

      SELECT CASE(isfctyp)
            !
            ! set surface gradients
            !
         CASE(1)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up%d,a_vp%d,umean,vmean)
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = dthcon
                  drdz(i,j) = drtcon
                  bfct(i,j) = g*zt%d(2)/(a_theta%d(2,i,j)*wspd(i,j)**2)
               END DO
            END DO
            zs = zt%d(2)/zrough
            CALL srfcscls(nxp,nyp,zt%d(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

            !
            ! get fluxes from profiles
            !
         CASE(2)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up%d,a_vp%d,umean,vmean)
            usum = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = a_theta%d(2,i,j) - sst*(p00/psrf)**rcp
                  drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
                  bfct(i,j) = g*zt%d(2)/(a_theta%d(2,i,j)*wspd(i,j)**2)
                  usum = usum + a_ustar%d(i,j)
               END DO
            END DO
            usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
            zs = (zrough/100.)
            IF (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
            CALL srfcscls(nxp,nyp,zt%d(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
             !
             ! get fluxes from bulk formulae with coefficients given by dthcon and
             ! drtcon
             !
         CASE(3)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up%d,a_vp%d,umean,vmean)
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = a_theta%d(2,i,j) - sst*(p00/psrf)**rcp
                  drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
                  IF (ubmin > 0.) THEN
                     a_ustar%d(i,j) = sqrt(zrough)* wspd(i,j)
                  ELSE
                     a_ustar%d(i,j) = abs(ubmin)
                  END IF
                  a_tstar%d(i,j) =  dthcon * wspd(i,j)*dtdz(i,j)/a_ustar%d(i,j)
                  a_rstar%d(i,j) =  drtcon * wspd(i,j)*drdz(i,j)/a_ustar%d(i,j)
                  bfct(i,j) = g*zt%d(2)/(a_theta%d(2,i,j)*wspd(i,j)**2)
               END DO
            END DO
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
             !
             ! fix surface temperature to yield a constant surface buoyancy flux
             !
         CASE(4)

            Vzt   = 10.* (log(zt%d(2)/zrough)/log(10./zrough))
            Vbulk = Vzt * (vonk/log(zt%d(2)/zrough))**2

            bfl(:) = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  bfl(1) = bfl(1)+a_theta%d(2,i,j)
                  bfl(2) = bfl(2)+rx(2,i,j) ! Juha: rx
               END DO
            END DO

            CALL double_array_par_sum(bfl,bfg,2)

            bfg(2) = bfg(2)/REAL((nxpg-4)*(nypg-4))
            bfg(1) = bfg(1)/REAL((nxpg-4)*(nypg-4))
            

            DO iterate = 1, 5
               bflx  = ((sst -bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst) -bfg(2))) &
                        * 0.5*(dn0%d(1)+dn0%d(2))*cp*Vbulk
               sst1 = sst + 0.1
               bflx1 = ((sst1-bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst1)-bfg(2))) &
                        * 0.5*(dn0%d(1)+dn0%d(2))*cp*Vbulk
               sst  = sst + 0.1* (dthcon - bflx) / (bflx1-bflx)
            END DO

            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  wt_sfc%d(i,j) = Vbulk * (sst -a_theta%d(2,i,j))
                  wq_sfc%d(i,j) = Vbulk * (rslf(psrf,sst) - rx(2,i,j)) ! Juha: rx
                  wspd(i,j)    = max(0.1,                                    &
                                 sqrt((a_up%d(2,i,j)+umean)**2+(a_vp%d(2,i,j)+vmean)**2))
                  bflx         = wt_sfc%d(i,j)*g/bfg(1) + g*ep2*wq_sfc%d(i,j)
                  a_ustar%d(i,j) = diag_ustar(zt%d(2),zrough,bflx,wspd(i,j))
                  uw_sfc%d(i,j)  = -a_ustar%d(i,j)*a_ustar%d(i,j)                  &
                                 *(a_up%d(2,i,j)+umean)/wspd(i,j)
                  vw_sfc%d(i,j)  = -a_ustar%d(i,j)*a_ustar%d(i,j)                  &
                                 *(a_vp%d(2,i,j)+vmean)/wspd(i,j)
                  ww_sfc%d(i,j)  = 0.
                  a_rstar%d(i,j) = wq_sfc%d(i,j)/a_ustar%d(i,j)
                  a_tstar%d(i,j) = wt_sfc%d(i,j)/a_ustar%d(i,j)
               END DO
            END DO

         CASE(5)

            !
            !Sami addition: Calculate surface fluxes from energy budget
            ! This is based on Acs et al: A coupled soil moisture and surface
            ! temperature prediction model, Journal of applied meteorology, 30, 1991
            !
            !total_sw = 0.0
            total_rw = 0.0
            total_la = 0.0
            total_se = 0.0
            total_pre = 0.0
            ffact = 1.
            !, a_rflx, precip
            !
            !   Calculate mean energy fluxes(Mean for each proceccors)
            !
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  !total_sw = total_sw+a_sflx%d(2,i,j)
                  !       WRITE(*,*) a_rflx(:,:,:)
                  total_rw = total_rw + a_rflx%d(2,i,j)
                  total_la = total_la + wq_sfc%d(i,j)*(0.5*(dn0%d(1)+dn0%d(2))*alvl)/ffact
                  total_se = total_se + wt_sfc%d(i,j)*(0.5*(dn0%d(1)+dn0%d(2))*cp)/ffact
                  total_pre=   total_pre +  a_rrate%d(2,i,j)
               END DO
            END DO
            !total_sw  = total_sw/REAL((nxp-4)*(nyp-4))
            total_rw  = total_rw/REAL((nxp-4)*(nyp-4))
            total_la  = total_la/REAL((nxp-4)*(nyp-4))
            total_se  = total_se/REAL((nxp-4)*(nyp-4))
            total_pre = total_pre/REAL((nxp-4)*(nyp-4))

        ! From energy fluxes calculate new sirface temperature
        sst1 =sst
        IF (.NOT. lConstSoilHeatCap) THEN
           C_heat = 1.29e3*(840.+4187.*thetaS1*W1) ! Eq 33
        END IF

        lambda=1.5e-7*C_heat

        ! Determine moisture at different depths in the ground
        !

            K1 = K_s1*W1**(2.*B1+3.)
            K2 = K_s1*W2**(2.*B1+3.)
            K3 = K_s3*W3**(2.*B3+3.)
            Kmean1 = (D1*K1+D2*K2)/(D1+D2)
            Kmean2 = (D2*K2+D3*K3)/(D2+D3)

            fii_1 = fii_s1*W1**(-B1)
            fii_2 = fii_s1*W2**(-B1)
            fii_3 = fii_s3*W3**(-B3)

            Q3 = K_s3*W3**(2.*B3+3.)*sin(0.05236) ! 3 degrees  Eq 8

            Q12 = Kmean1*( 2.*(fii_1-fii_2)/(D1+D2)+1.0)
            Q23 = Kmean2*( 2.*(fii_2-fii_3)/(D2+D3)+1.0)

            IF (.NOT. lConstSoilWater) THEN
               W1 = W1 + (1./(thetaS1*D1)) * ( - Q12 + (total_pre/(0.5*(dn0%d(1)+dn0%d(2))*alvl)/ffact)         &
                                               - (total_la/(0.5*(dn0%d(1)+dn0%d(2))*alvl)/ffact)/rowt ) *dtl
               W2 = W2 + (1./(thetaS2*D2)) * (Q12-Q23) * dtl
               W3 = W3 + (1./(thetaS3*D3)) * (Q23-Q3) * dtl
            END IF
            !
            !  Following is copied from CASE (2). No idea if this is valid or not..
            !

            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up%d,a_vp%d,umean,vmean)
            usum = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2

                  dtdz(i,j) = a_theta%d(2,i,j) - sst1*(p00/psrf)**rcp

                  ff1 = 1.0
                  IF(W1 < 0.75) ff1 = W1/0.75
                  ! Flux of moisture is limited by water content.
                  drdz(i,j) = a_rp%d(2,i,j) - ff1*rslf(psrf,sst1)   !rslf(psrf,min(sst1,280.))  !  a_rv changed to a_rp (by Zubair)
                  !
                  bfct(i,j) = g*zt%d(2)/(a_theta%d(2,i,j)*wspd(i,j)**2)
                  usum = usum + a_ustar%d(i,j)
               END DO
            END DO
            usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
            zs = (zrough/100.)
            IF (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
            CALL srfcscls(nxp,nyp,zt%d(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

           sst1 = sst1-(total_rw+total_la+total_se+( SQRT(lambda*C_heat*7.27e-5/(2.0)) *(SST1-deepSoilTemp)))&
                /(2.0e-2*C_heat+ SQRT( lambda*C_heat/(2.0*7.27e-5)) )*dtl

            sst = sst1
         

	CASE(6)
        
         !
         ! Fix thermodynamic fluxes at surface given in energetic units 
         ! Calculate momentum fluxes from winds as in the DEFAULT METHOD
         ! Fire adds latent and sensible heat increasing buoyancy flux
         !
         ! The fire spreads using the spread rate calculated with Rothermel's model (1972)
         ! The fuel burned and heat released is calculated as in Mandel et al. (2011),
         ! the speed of burning (combustion reaction velocity as kg fuel consumed per time)
         ! is assumed to be independent of the wind speed and the fuel moisture.
         ! This speed is given as Tf or fuel burnt time that is equal to the 
         ! fuel weight (w) divided by 0.8514 (Eq.3 in Mandel-2011)
         ! 
         ! Mandel, J., Beezley, J. D., and Kochanski, A. K.: 
         ! Coupled atmosphere-wildland fire modeling with WRF 3.3 and SFIRE 2011, 
         ! Geosci. Model Dev., 4, 591–610, 
         ! https://doi.org/10.5194/gmd-4-591-2011, 2011. 
         ! Rothermel, R.C. (1972). A mathematical model for predicting fire spread 
         ! in wildland fuels. USDA Forest Service, 
         ! Intermountain Forest and Range Experiment Station, Research Paper INT–115, 40 p.
         ! https://research.fs.usda.gov/treesearch/32533
         
         ! Background surface fluxes from the runles
            ffact = 1.  
            lh_flx(:,:) = lh_flx(:,:) + drtcon
            sh_flx(:,:) = sh_flx(:,:) + dthcon
            total_la_fire = 0.
	    total_se_fire = 0.            
          
          ! Field of fire spread rates affected by horizontal wind velocity R = Ro(1+phiw)                      
            DO j = 3, nyp-2
                  DO i = 3, nxp-2
                     wspd(i,j) = max(0.1, &
                                sqrt((a_up%d(2,i,j)+umean)**2+(a_vp%d(2,i,j)+vmean)**2))                    	
                     a_firespread%d(i,j) = a_R0%d(i,j) + a_phiwc%d(i,j)*(MIN(wspd(i,j),6.0)*60/0.3048)**a_phiwb%d(i,j)  
                     ! Coen et al. (2013) wspd capped at 6 m/s         
                  END DO
            END DO
	    
	  ! Looping through the surface
            DO j = 3, nyp-2
               DO i = 3, nxp-2  
                  ! The cell has fuel (vegetation) and it is ignited
                  ! Each cell is divided in 20 fractions and fire spreads progressively
                  ! If the cell does not have fuel(i.e. vegetation) fgi =0                         
                  IF (a_fgi%d(i,j)>0 .AND. time > a_ignitiontime%d(i,j)) THEN   
                       ! Fraction of the cell ignited
                       fbcell = a_areaburnt%d(i,j)/(deltax*deltay)
                       IF      (fbcell<=0.05 .AND. time > a_ignitiontime%d(i,j)) THEN 
			  		ignitiontimecell(i,j,1) = a_ignitiontime%d(i,j)  
		       ELSEIF  ((fbcell>0.05 .AND. fbcell <= 0.1) .AND. (time > ignitiontimecell(i,j,1) .AND. &  
		       			ignitiontimecell(i,j,2) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,2) = time		       		                                          
		       ELSEIF  ((fbcell>0.1 .AND. fbcell <= 0.15) .AND. (time > ignitiontimecell(i,j,2) .AND. &
			  		ignitiontimecell(i,j,3) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,3) = time
		       ELSEIF  ((fbcell>0.15 .AND. fbcell <= 0.2) .AND. (time > ignitiontimecell(i,j,3) .AND. &
			  		ignitiontimecell(i,j,4) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,4) = time           
		       ELSEIF  ((fbcell>0.2 .AND. fbcell <= 0.25) .AND. (time > ignitiontimecell(i,j,4).AND. & 
		          		ignitiontimecell(i,j,5) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,5) = time
		       ELSEIF  ((fbcell>0.25 .AND. fbcell <= 0.30) .AND. (time > ignitiontimecell(i,j,5).AND. &
		          		ignitiontimecell(i,j,6) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,6) = time
		       ELSEIF  ((fbcell>0.30 .AND. fbcell <= 0.35) .AND. (time > ignitiontimecell(i,j,6).AND. &
		          		ignitiontimecell(i,j,7) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,7) = time 
		       ELSEIF  ((fbcell>0.35 .AND. fbcell <= 0.4) .AND. (time > ignitiontimecell(i,j,7).AND. &
		          		ignitiontimecell(i,j,8) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,8) = time
		       ELSEIF  ((fbcell>0.4 .AND. fbcell <= 0.45) .AND. (time > ignitiontimecell(i,j,8).AND. &
		          	 	ignitiontimecell(i,j,9) > 1.E12)) THEN 
		                  ignitiontimecell(i,j,9) = time  
		       ELSEIF  ((fbcell>0.45 .AND. fbcell <= 0.5) .AND. (time > ignitiontimecell(i,j,9).AND. &
		          		ignitiontimecell(i,j,10) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,10) = time
		       ELSEIF  ((fbcell>0.5 .AND. fbcell <= 0.55) .AND. (time > ignitiontimecell(i,j,10).AND. &
		          		ignitiontimecell(i,j,11) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,11) = time  
		       ELSEIF  ((fbcell>0.55 .AND. fbcell < 0.60) .AND. (time > ignitiontimecell(i,j,11).AND. &
		          		ignitiontimecell(i,j,12) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,12) = time
		       ELSEIF  ((fbcell>0.6 .AND. fbcell < 0.65) .AND. (time > ignitiontimecell(i,j,12).AND. &
		          		ignitiontimecell(i,j,13) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,13) = time
		       ELSEIF  ((fbcell>0.65 .AND. fbcell < 0.7) .AND. (time > ignitiontimecell(i,j,13).AND. &
		          		ignitiontimecell(i,j,14) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,14) = time                      
		       ELSEIF  ((fbcell>0.7 .AND. fbcell < 0.75) .AND. (time > ignitiontimecell(i,j,14).AND. &
		          		ignitiontimecell(i,j,15) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,15) = time             
		       ELSEIF  ((fbcell>0.75 .AND. fbcell < 0.8) .AND. (time > ignitiontimecell(i,j,15).AND. &
		          		ignitiontimecell(i,j,16) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,16) = time            
		       ELSEIF  ((fbcell>0.8 .AND. fbcell < 0.85) .AND. (time > ignitiontimecell(i,j,16).AND. &
		          		ignitiontimecell(i,j,17) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,17) = time    
		       ELSEIF  ((fbcell>0.85 .AND. fbcell < 0.9) .AND. (time > ignitiontimecell(i,j,17).AND. &
		          		ignitiontimecell(i,j,18) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,18) = time  
		       ELSEIF  ((fbcell>0.9 .AND. fbcell < 0.95) .AND. (time > ignitiontimecell(i,j,18).AND. &
		          		ignitiontimecell(i,j,19) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,19) = time 
		       ELSEIF  ((fbcell>0.95 .AND. fbcell < 0.999) .AND. (time > ignitiontimecell(i,j,19).AND. &
		          		ignitiontimecell(i,j,20) > 1.E12)) THEN   
		                  ignitiontimecell(i,j,20) = time                             
		       ELSEIF  (fbcell >= 0.999 .AND. ignitiontimecell(i,j,20)>1.0E12) THEN		          		 
		                  ignitiontimecell(i,j,20) = time     
		       END IF
		              	       	    
          	       fuel_burnt = 0.
          	       total_se_fire = 0.
          	       total_la_fire = 0.
          	       area_burnt_cell = 0.
          	       DO kk=1,20          	          
          	          IF (time > ignitiontimecell(i,j,kk)) THEN           	                
			  	tt = time -ignitiontimecell(i,j,kk)	
			  	! Similar to Mandel et al. 2011 using the fraction of remaining fuel		  	
			  	fraction_burnt = (EXP(-tt/a_weight%d(i,j)) - EXP(-(tt+dtl)/a_weight%d(i,j))) 
			  	! The fire moves at the same velocity in x and y directions
                       		! as a growing circle
			  	IF (areaignitedcell(i,j,kk) < deltax*deltay/20) THEN	  	
			  	    areaignitedcell(i,j,kk) = areaignitedcell(i,j,kk) + & 
			  	              MAX(pi*((tt+dtl)*a_firespread%d(i,j))**2 - &
                  	                          pi*(tt*a_firespread%d(i,j))**2, 0.)
			  	ELSE 			  				  	              
			  	    areaignitedcell(i,j,kk) = deltax*deltay/20	  				  				  	
			  	END IF
			  	!
			  	fuel_burnt = fuel_burnt + areaignitedcell(i,j,kk)*a_fgi%d(i,j)*fraction_burnt       	                                                               
                          	! for checking purposes
                          	!WRITE(*,*) 'kk, ignitiontimecell', kk, a_ignitiontime%d(i,j), ignitiontimecell(i,j,kk)        
                          	!WRITE(*,*) 'area_burnt_cell, area_burnt', areaignitedcell(i,j,kk), area_burnt_cell
                          	!WRITE(*,*) 'fuelburnt',fuel_burnt                                                           	   			     
		          	! Calculating the surface fluxes coming from combustion		          	
		          	! Average sensible heat released in time interval (t, t+Deltat) Mandel-2011-Eq.4 in W/m2
		          	total_se_fire = total_se_fire + a_fgi%d(i,j)*fraction_burnt/dtl * &
                                                1/(1+a_fuelmcg%d(i,j))*cmbcnst * &
                                                areaignitedcell(i,j,kk)/(deltax*deltay/20)
		          	! Average latent heat released in time interval (t, t+Deltat) Mandel-2011-Eq.5 in W/m2
		          	total_la_fire = total_la_fire + a_fgi%d(i,j)*fraction_burnt/dtl* &
                                                (a_fuelmcg%d(i,j)+0.56)/(1+a_fuelmcg%d(i,j))*alvl* &
                                                areaignitedcell(i,j,kk)/(deltax*deltay/20)
		          	tt = 0.
		          	fraction_burnt = 0.
		          	area_burnt_cell = area_burnt_cell + areaignitedcell(i,j,kk) 		          	
		          END IF
		       END DO
		       
		       !WRITE(*,*) 'area_burnt_cell, area_burnt', areaignitedcell(i,j,:)
		       a_fuelburnt%d(i,j) = fuel_burnt
		       a_areaburnt%d(i,j) = MIN(area_burnt_cell, deltax*deltay)
		       a_tcrit%d(i,j) = ignitiontimecell(i,j,20)
		       sh_flx(i,j) = sh_flx(i,j) + total_se_fire 
		       lh_flx(i,j) = lh_flx(i,j) + total_la_fire		
		       !WRITE(*,*) 'sh_fire, lh_fire', total_se_fire, total_la_fire 		                 
		 END IF
		 		 	  		   
                 wt_sfc%d(i,j) = sh_flx(i,j)/(0.5*(dn0%d(1)+dn0%d(2))*cp)
                 wq_sfc%d(i,j) = lh_flx(i,j)/(0.5*(dn0%d(1)+dn0%d(2))*alvl)                  
                  
                 IF (ubmin > 0.) THEN
                     ! It was th00, but it must be reference potential temperature within the ABL  
                     bflx = g*wt_sfc%d(i,j)/a_theta%d(2,i,j)    
                     IF (level >= 2) bflx = bflx + g*ep2*wq_sfc%d(i,j)
                     a_ustar%d(i,j) = diag_ustar(zt%d(2),a_fcz0%d(i,j),bflx,wspd(i,j))
                 ELSE
                     a_ustar%d(i,j) = abs(ubmin)
                 END IF
                 ffact = a_ustar%d(i,j)*a_ustar%d(i,j)/wspd(i,j)
                 uw_sfc%d(i,j)  = -ffact*(a_up%d(2,i,j)+umean)
                 vw_sfc%d(i,j)  = -ffact*(a_vp%d(2,i,j)+vmean)
                 ww_sfc%d(i,j)  = 0.
                 a_rstar%d(i,j) = wq_sfc%d(i,j)/a_ustar%d(i,j)
                 a_tstar%d(i,j) = wt_sfc%d(i,j)/a_ustar%d(i,j)
            END DO
          END DO
          
                 
          ! Updating area burnt
          CALL cyclics2d(nxp,nyp,a_areaburnt%d,req)                
          CALL cyclicc2d(nxp,nyp,a_areaburnt%d,req) 
          
          ! Checking if fire had spread across cells
          CALL update_ignition(time)   	  
	  
	  CALL cyclics2d(nxp,nyp,a_ignitiontime%d,req)     
          CALL cyclicc2d(nxp,nyp,a_ignitiontime%d,req)             
	  
                   
	  ! Reset for next timestep
	  lh_flx = 0.
	  sh_flx = 0.
         !
         ! fix thermodynamic fluxes at surface given values in energetic
         ! units and calculate  momentum fluxes from winds
         !
         CASE DEFAULT
            ffact = 1.  
            lh_flx(:,:) = lh_flx(:,:) + drtcon
            sh_flx(:,:) = sh_flx(:,:) + dthcon
                        
            IF (zrough <= 0.) THEN
               usum = 0.
               DO j = 3, nyp-2
                  DO i = 3, nxp-2
                     usum = usum + a_ustar%d(i,j)
                  END DO
               END DO
               usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
               zs = max(0.0001,(0.016/g)*usum**2)
            ELSE
               zs = zrough
            END IF

            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  wt_sfc%d(i,j) = sh_flx(i,j)/(0.5*(dn0%d(1)+dn0%d(2))*cp)
                  wq_sfc%d(i,j) = lh_flx(i,j)/(0.5*(dn0%d(1)+dn0%d(2))*alvl)
                  
                  wspd(i,j)   = max(0.1,                                    &
                                sqrt((a_up%d(2,i,j)+umean)**2+(a_vp%d(2,i,j)+vmean)**2))
                  IF (ubmin > 0.) THEN
                     bflx = g*wt_sfc%d(i,j)/th00   ! Why was the index 1,1??
                     IF (level >= 2) bflx = bflx + g*ep2*wq_sfc%d(i,j)
                     a_ustar%d(i,j) = diag_ustar(zt%d(2),zs,bflx,wspd(i,j))
                  ELSE
                     a_ustar%d(i,j) = abs(ubmin)
                  END IF
                  ffact = a_ustar%d(i,j)*a_ustar%d(i,j)/wspd(i,j)
                  uw_sfc%d(i,j)  = -ffact*(a_up%d(2,i,j)+umean)
                  vw_sfc%d(i,j)  = -ffact*(a_vp%d(2,i,j)+vmean)
                  ww_sfc%d(i,j)  = 0.
                  a_rstar%d(i,j) = wq_sfc%d(i,j)/a_ustar%d(i,j)
                  a_tstar%d(i,j) = wt_sfc%d(i,j)/a_ustar%d(i,j)
               END DO
            END DO
            
            ! Reset for next timestep
            lh_flx = 0.
            sh_flx = 0.
            
      END SELECT


      RETURN
   END SUBROUTINE surface
   !
   ! -------------------------------------------------------------------
   ! GET_SWNDS: returns surface winds valid at cell centers
   !
   SUBROUTINE get_swnds(n1,n2,n3,usfc,vsfc,wspd,up,vp,umean,vmean)

      IMPLICIT NONE

      INTEGER, INTENT (in) :: n1, n2, n3
      REAL, INTENT (in)    :: up(n1,n2,n3), vp(n1,n2,n3), umean, vmean
      REAL, INTENT (out)   :: usfc(n2,n3), vsfc(n2,n3), wspd(n2,n3)

      INTEGER :: i, j, ii, jj

      DO j = 3, n3-2
         jj = j-1
         DO i = 3, n2-2
            ii = i-1
            usfc(i,j) = (up(2,i,j)+up(2,ii,j))*0.5+umean
            vsfc(i,j) = (vp(2,i,j)+vp(2,i,jj))*0.5+vmean
            wspd(i,j) = max(abs(ubmin),sqrt(usfc(i,j)**2+vsfc(i,j)**2))
         END DO
      END DO

   END SUBROUTINE get_swnds
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
   REAL FUNCTION diag_ustar(z,z0,bflx,wnd)

      USE defs, ONLY : vonk

      IMPLICIT NONE

      REAL, PARAMETER   :: am  = 4.8    !   "          "         "
      REAL, PARAMETER   :: bm  = 19.3   !   "          "         "
      REAL, PARAMETER   :: eps = 1.e-10 ! non-zero, small number

      REAL, INTENT (in) :: z             ! height where u locates
      REAL, INTENT (in) :: z0            ! momentum roughness height
      REAL, INTENT (in) :: bflx          ! surface buoyancy flux (m^2/s^3)
      REAL, INTENT (in) :: wnd           ! wind speed at z

      INTEGER :: iterate
      REAL    :: lnz, klnz, c1, x, psi1, zeta, lmo, ustar

      lnz  = log(z/z0)
      klnz = vonk/lnz
      c1   = 3.14159/2. - 3.*log(2.)

      ustar = wnd*klnz
      IF (bflx /= 0.0) THEN
         DO iterate = 1, 4
            lmo  = -(ustar**3)/(bflx*vonk + eps)
            zeta = z/lmo
            IF (zeta > 0.) THEN
               ustar = vonk*wnd / (lnz + am*zeta)
            ELSE
               x     = sqrt( sqrt( 1.0 - bm*zeta ) )
               psi1  = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + c1
               ustar = wnd*vonk/(lnz - psi1)
            END IF
         END DO
      END IF

      diag_ustar = ustar

      RETURN
   END FUNCTION diag_ustar
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
   ! Layer, in Workshop on Micrometeorology, pages 67-100.
   !
   ! Code writen March, 1999 by Bjorn Stevens
   !
   SUBROUTINE srfcscls(n2,n3,z,z0,th00,u,dth,drt,ustar,tstar,rstar)
   

      USE defs, ONLY : vonk, g, ep2

      IMPLICIT NONE

      REAL, PARAMETER     :: ah  =  7.8   ! stability function parameter
      REAL, PARAMETER     :: bh  = 12.0   !   "          "         "
      REAL, PARAMETER     :: am  =  4.8   !   "          "         "
      REAL, PARAMETER     :: bm  = 19.3   !   "          "         "
      REAL, PARAMETER     :: pr  = 0.74   ! prandlt number
      REAL, PARAMETER     :: eps = 1.e-10 ! non-zero, small number

      INTEGER, INTENT(in) :: n2,n3         ! span of indicies covering plane
      REAL, INTENT(in)    :: z             ! height where u & T locate
      REAL, INTENT(in)    :: z0            ! momentum roughness height
      REAL, INTENT(in)    :: th00          ! reference temperature
      REAL, INTENT(in)    :: u(n2,n3)      ! velocities at z
      REAL, INTENT(in)    :: dth(n2,n3)    ! theta (th(z) - th(z0))
      REAL, INTENT(in)    :: drt(n2,n3)    ! qt(z) - qt(z0)
      TYPE(FloatArray2d), INTENT(inout) :: ustar  ! scale velocity
      TYPE(FloatArray2d), INTENT(inout) :: tstar  ! scale temperature
      TYPE(FloatArray2d), INTENT(inout) :: rstar  ! scale value of qt

      LOGICAL, SAVE :: first_call = .TRUE.
      INTEGER :: i,j,iterate
      REAL    :: lnz, klnz, betg, cnst1, cnst2
      REAL    :: x, y, psi1, psi2, zeta, lmo, dtv

      lnz   = log(z/z0)
      klnz  = vonk/lnz
      betg  = th00/g
      cnst2 = -log(2.)
      cnst1 = 3.14159/2. + 3.*cnst2

      DO j = 3, n3-2
         DO i = 3, n2-2
            dtv = dth(i,j) + ep2*th00*drt(i,j)
            !
            ! stable CASE
            !
            IF (dtv > 0.) THEN
               x     = (betg*u(i,j)**2)/dtv
               y     = (am - 0.5*x)/lnz
               x     = (x*ah - am**2)/(lnz**2)
               lmo   = -y + sqrt(x+y**2)
               zeta  = z/lmo
               ustar%d(i,j) =  vonk*u(i,j)  /(lnz + am*zeta)
               tstar%d(i,j) = (vonk*dtv/(lnz + ah*zeta))/pr
               !
               ! Neutral case
               !
            ELSE IF (dtv == 0.) THEN
               ustar%d = vonk*u(i,j)  /lnz
               tstar%d = vonk*dtv/(pr*lnz)
               !
               ! ustable case, start iterations from values at previous tstep,
               ! unless the sign has changed or if it is the first call, then
               ! use neutral values.
               !
            ELSE
               IF (first_call .OR. tstar%d(i,j)*dtv <= 0.) THEN
                  ustar%d(i,j) = u(i,j)*klnz
                  tstar%d(i,j) = (dtv*klnz/pr)
               END IF

               DO iterate = 1, 3
                  lmo  = betg*ustar%d(i,j)**2/(vonk*tstar%d(i,j))
                  zeta = z/lmo
                  x    = sqrt(sqrt( 1.0 - bm*zeta ) )
                  psi1 = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + cnst1
                  y    = sqrt(1.0 - bh*zeta)
                  psi2 = log(1.0 + y) + cnst2
                  ustar%d(i,j) = u(i,j)*vonk/(lnz - psi1)
                  tstar%d(i,j) = (dtv*vonk/pr)/(lnz - psi2)
               END DO
            END IF

            rstar%d(i,j) = tstar%d(i,j)*drt(i,j)/(dtv + eps)
            tstar%d(i,j) = tstar%d(i,j)*dth(i,j)/(dtv + eps)
         END DO
      END DO

      first_call = .FALSE.

      RETURN
   END SUBROUTINE srfcscls
   !
   ! ----------------------------------------------------------------------
   ! Subroutine: sfcflxs:  this routine returns the surface fluxes based
   ! on manton-cotton algebraic surface layer equations.
   !
   SUBROUTINE sfcflxs(n2,n3,vk,ubar,u,v,xx,us,ts,rs,uw,vw,tw,rw,ww)
      IMPLICIT NONE
      REAL, PARAMETER     :: cc = 4.7,eps = 1.e-20

      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(in)    :: ubar(n2,n3),u(n2,n3),v(n2,n3),xx(n2,n3),vk
      TYPE(FloatArray2d), INTENT(in)    :: us, ts, rs
      TYPE(FloatArray2d), INTENT(inout)   :: uw, vw, tw, rw, ww

      REAL :: x(n2,n3),y(n2,n3)
      INTEGER i,j

      DO j = 3, n3-2
         DO i = 3, n2-2

            uw%d(i,j) = -(u(i,j)/(ubar(i,j)+eps))*us%d(i,j)**2
            vw%d(i,j) = -(v(i,j)/(ubar(i,j)+eps))*us%d(i,j)**2
            tw%d(i,j) = -ts%d(i,j)*us%d(i,j)
            rw%d(i,j) = -rs%d(i,j)*us%d(i,j)

            x(i,j) = xx(i,j)*vk*ts%d(i,j)*(ubar(i,j)/us%d(i,j))**2
            x(i,j) = x(i,j)*sqrt(sqrt(1.-15.*min(0.,x(i,j)))) &
                     /(1.0+cc*max(0.,x(i,j)))
            y(i,j) = sqrt((1.-2.86*x(i,j))/(1.+x(i,j)* &
                     (-5.39+x(i,j)*6.998 )))
            ww%d(i,j)= (0.27*max(6.25*(1.-x(i,j))*y(i,j),eps)-&
                      1.18*x(i,j)*y(i,j))*us%d(i,j)**2
         END DO
      END DO
      RETURN
    END SUBROUTINE sfcflxs

! ------------------------------------------------------------------------------------
! Subroutine: update_ignition_time
! If a cell has burnt completely, the fire spreads to surrounding cells 
! as long as there is a positive gradient in the spread rate and
! the ignition time is below the non-ignited flag of 1.0E15 s

SUBROUTINE update_ignition(time)              
   USE mo_diag_state, ONLY: a_ignitiontime, a_areaburnt, a_firespread,a_fgi
   IMPLICIT NONE
   INTEGER :: i, j
   REAL, INTENT(IN) :: time
    
   DO j = 2, nyp-1
     DO i = 2, nxp-1 	                              
	IF (a_areaburnt%d(i,j)>=deltax*deltay) THEN
   	   IF (a_firespread%d(i,j+1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i,j+1)>1.E12) THEN   	    	   	   
   	      a_ignitiontime%d(i,j+1)   = time
   	   ELSE IF (a_firespread%d(i,j-1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i,j-1)>1.E12) THEN
   	      a_ignitiontime%d(i,j-1)   = time
   	   ELSE IF (a_firespread%d(i+1,j)>a_firespread%d(i,j).AND. a_ignitiontime%d(i+1,j)>1.E12) THEN   	    	   	   
   	      a_ignitiontime%d(i+1,j)   = time
   	   ELSE IF (a_firespread%d(i-1,j)>a_firespread%d(i,j).AND. a_ignitiontime%d(i-1,j)>1.E12) THEN
   	      a_ignitiontime%d(i-1,j)   = time
   	   ELSE IF (a_firespread%d(i+1,j+1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i+1,j+1)>1.E12) THEN   	    	   	   
   	      a_ignitiontime%d(i+1,j+1) = time
   	   ELSE IF (a_firespread%d(i+1,j-1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i+1,j-1)>1.E12) THEN
   	      a_ignitiontime%d(i+1,j-1) = time     
   	   ELSE IF (a_firespread%d(i-1,j+1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i-1,j+1)>1.E12) THEN   	    	   	   
   	      a_ignitiontime%d(i-1,j+1) = time
   	   ELSE IF (a_firespread%d(i-1,j-1)>a_firespread%d(i,j).AND. a_ignitiontime%d(i-1,j-1)>1.E12) THEN
   	      a_ignitiontime%d(i-1,j-1) = time  
   	   END IF 
       END IF
     END DO
   END DO 
   
   a_ignitiontime%d(1,:) = a_ignitiontime%d(2,:) 
   a_ignitiontime%d(nxp,:) = a_ignitiontime%d(nxp-1,:)
   a_ignitiontime%d(:,1) = a_ignitiontime%d(:,nyp-1) 
   a_ignitiontime%d(:,nyp) = a_ignitiontime%d(:,nyp)
   
END SUBROUTINE update_ignition


! ---------------------------------------------------------------------------
! Subroutine surface_state:  This routine reads properties
! needed to simulate burning vegetation on the surface
! 
SUBROUTINE surface_state()
  USE ncio, ONLY : open_surf_nc, read_surf_nc_2d, close_nc  
  USE mpi_interface, ONLY : xoffset, yoffset, wrxid, wryid, nxpg, nypg,   &
                            myid, nyprocs, nxprocs, ranktable
  USE mo_diag_state, ONLY: a_tskin, a_qskin, a_fgi, a_weight, a_fcz0, & 
  			   a_fuelmcg, a_ignitiontime, a_fuelburnt, &
  			   a_firespread, a_areaburnt, a_phiwc, a_phiwb, & 
  			   a_tcrit, a_R0
  USE grid, ONLY: sst, psrf,nxp,nyp
  
  IMPLICIT NONE
  LOGICAL :: READ_NC
  INTEGER :: ncid, nvar,nxp_global,nyp_global
  INTEGER :: istart, iend, jstart, jend
  REAL :: rskin
	
  REAL, ALLOCATABLE :: fgig(:,:), weightg(:,:), fcz0g(:,:), R0g(:,:)
  REAL, ALLOCATABLE :: fuelmcgg(:,:), phiwcg(:,:), phiwbg(:,:) 
  
  ! Read the NetCDF input when it is available
  INQUIRE(FILE='datafiles/surface_in.nc', EXIST=READ_NC)

  ! Open the input file
  IF (READ_NC) CALL open_surf_nc(ncid, nxp_global, nyp_global) 
  
  ! For checking purposes
  !WRITE(*,*) 'Internal', nxp_global,nyp_global 
  !WRITE(*,*) 'ranktable',ranktable
  !WRITE(*,*) 'wrxid,wryid', wrxid, wryid
  !WRITE(*,*) 'xoffset, yoffset',xoffset,yoffset
  
 ! Allocate input variables
  ALLOCATE(fgig(nxp_global,nyp_global), weightg(nxp_global,nyp_global), fcz0g(nxp_global,nyp_global))
  ALLOCATE(fuelmcgg(nxp_global,nyp_global), phiwcg(nxp_global,nyp_global), phiwbg(nxp_global,nyp_global))
  ALLOCATE(R0g(nxp_global,nyp_global))
  
  ALLOCATE(firespreadg(nxp_global,nyp_global),ignitiontimeg(nxp_global,nyp_global),areaburntg(nxp_global,nyp_global)) 
  ALLOCATE(ignitiontimecell(nxp,nyp,20))
  ALLOCATE(areaignitedcell(nxp,nyp,20))
 
  areaburntg = 0.
  ignitiontimecell = 1.0E15
  areaignitedcell = 0.
  
  IF (READ_NC) THEN
     ! Read the surface properties
     CALL read_surf_nc_2d(ncid, 'fgi', nxp_global, nyp_global, fgig)
     CALL read_surf_nc_2d(ncid, 'weight', nxp_global, nyp_global, weightg)
     CALL read_surf_nc_2d(ncid, 'fcz0',  nxp_global, nyp_global,fcz0g)
     CALL read_surf_nc_2d(ncid, 'fuelmcg',  nxp_global, nyp_global, fuelmcgg)
     CALL read_surf_nc_2d(ncid, 'ignitiontime', nxp_global, nyp_global, ignitiontimeg)	     
     CALL read_surf_nc_2d(ncid, 'R0', nxp_global, nyp_global, R0g)
     CALL read_surf_nc_2d(ncid, 'phiwc', nxp_global, nyp_global, phiwcg)
     CALL read_surf_nc_2d(ncid, 'phiwb', nxp_global, nyp_global, phiwbg)
     CALL close_nc(ncid)
     WRITE(*,*) 'Surface properties read successfully from datafiles/surface_in.nc'
  ELSE
     WRITE(*,*) 'No datafiles/surface_in.nc was read'
     WRITE(*,*) 'No fuel or vegetation in the model domain'
     WRITE(*,*) 'No surface_in.nc found — using defaults.'
     fgig = 0.0
     weightg = 7.
     fcz0g = 0.1
     fuelmcgg = 0.0
     ignitiontimeg =1.0E15 ! No fire because time<ignitiontime, then no fire
     R0g = 0.0
     phiwcg = 0.0
     phiwbg= 1.0	  
  END IF
 
  firespreadg = R0g
  
  istart = MAX(wrxid * (nxp_global-2)/nxprocs ,1)
  iend   = MIN((wrxid+1)*(nxp_global-2)/nxprocs+ 3, nxp_global)
  jstart = MAX(wryid * (nyp_global-2)/nyprocs, 1) 
  jend   = MIN((wryid+1)*(nyp_global-2)/nyprocs+3, nyp_global)
  
  ! for checking purposes
  !WRITE(*,*) 'istart,iend, jstart,jend', istart,iend, jstart,jend
  
  a_fgi%d    = fgig(istart:iend,jstart:jend)
  a_weight%d = weightg(istart:iend,jstart:jend)/0.8514 !Mandel 2011 Eq.3
  a_fcz0%d   = fcz0g(istart:iend,jstart:jend)
  a_fuelmcg%d = fuelmcgg(istart:iend,jstart:jend)
  a_ignitiontime%d = ignitiontimeg(istart:iend,jstart:jend)
  a_tskin%d = sst
  ! Assuming saturated surface
  rskin = rslf(psrf,sst) 
  a_qskin%d  = rskin / (1 + rskin) 
  a_fuelburnt%d =0.
  a_R0%d = R0g(istart:iend,jstart:jend)
  a_areaburnt%d  = 0.
  a_phiwc%d = phiwcg(istart:iend,jstart:jend)
  a_phiwb%d = phiwbg(istart:iend,jstart:jend)  
  a_tcrit%d  = 1.0E15 ! Same as no fire 
  a_firespread%d = R0g(istart:iend,jstart:jend)
  ignitiontimecell(:,:,1) = ignitiontimeg(istart:iend,jstart:jend)
  
END SUBROUTINE surface_state
   
 
END MODULE srfc
