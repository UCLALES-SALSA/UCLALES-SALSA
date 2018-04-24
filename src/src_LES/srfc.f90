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

   INTEGER :: isfctyp = 0
   REAL    :: zrough =  0.1e-3   !!!!! CHECK THAT IN METERS EVERYWHERE!!
   REAL    :: ubmin  =  0.20
   REAL    :: dthcon = 100.0
   REAL    :: drtcon = 0.0



   ! Sami added ----->
   ! Initial values for surface properties
   ! REAL :: W1 = 0.9   !Water content      ... Definition in grid now because of restrat files
   ! REAL :: W2 = 0.9
   ! REAL :: W3 = 0.9

   REAL ::  B1 = 6.5
   REAL ::  B3 = 7.26
   REAL ::  K_s1 = 5.8e-5
   REAL ::  K_s3 = 0.89e-5
   REAL ::  fii_s1 = -0.036
   REAL ::  fii_s3 = -0.085
   REAL ::  thetaS1 = 0.98 ! Soil porosity
   REAL ::  thetaS2 = 0.98
   REAL ::  thetaS3 = 0.488
   REAL ::  D1 = 0.1  !Depth of different layers !
   REAL ::  D2 = 0.3
   REAL ::  D3 = 0.6

! <--- Sami added

! Juha moved/added
  ! for isfctyp == 5:
  REAL :: C_heat = 2.e6                 ! Surface heat capacity
  REAL :: deepSoilTemp = 280.            ! Assumed deep soil layer temperature
  LOGICAL :: lConstSoilWater = .FALSE.   ! Keep the value(s) of surface water content constant (as specified in NAMELIST)
  LOGICAL :: lConstSoilHeatCap = .FALSE. ! Keep the value of surface heat capacity constant (as specified in NAMELIST)

CONTAINS
   !
   ! --------------------------------------------------------------------------
   ! SURFACE: Calcualtes surface fluxes using an algorithm chosen by ISFCLYR
   ! and fills the appropriate 2D arrays
   !
   !     default: specified thermo-fluxes (drtcon, dthcon)
   !     isfclyr=1: specified surface layer gradients (drtcon, dthcon)
   !     isfclyr=2: fixed lower boundary of water at certain sst
   !     isfclyr=3: bulk aerodynamic law with coefficeints (drtcon, dthcon)
   !
   ! Modified for level 4: a_rv replaced by a local variable rx, which has
   ! values a_rv if level < 4, and a_rp if level == 4 (i.e. water vapour mixrat in both cases)
   !
   ! Juha Tonttila, FMI, 2014
   !

   SUBROUTINE surface()

      USE defs, ONLY: vonk, p00, rcp, g, cp, alvl, ep2
      USE grid, ONLY: nzp, nxp, nyp, a_up, a_vp, a_theta, a_rv, a_rp, zt, dzt, psrf, th00,  &
                      umean, vmean, a_ustar, a_tstar, a_rstar, uw_sfc, vw_sfc, ww_sfc,      &
                      wt_sfc, wq_sfc, dn0, level,dtl, a_sflx, a_rflx, precip, a_dn,         &
                      W1,W2,W3, mc_ApVdom, dtlt,sst
      USE thrm, ONLY: rslf
      USE stat, ONLY: sfc_stat, sflg, mcflg, acc_massbudged
      USE mpi_interface, ONLY : nypg, nxpg, double_array_par_sum


      IMPLICIT NONE
      REAL :: dtdz(nxp,nyp), drdz(nxp,nyp), usfc(nxp,nyp), vsfc(nxp,nyp),       &
              wspd(nxp,nyp), bfct(nxp,nyp)
      REAL :: rx(nzp,nxp,nyp)


      REAL :: total_sw, total_rw, total_la, total_se, total_pre  ! Sami added
      REAL :: lambda ! Sami added
      REAL :: K1,K2,K3,Kmean1,Kmean2,fii_1,fii_2,fii_3,Q3,Q12,Q23,ff1  ! Sami added



      INTEGER :: i, j, iterate
      REAL    :: zs, bflx, ffact, sst1, bflx1, Vbulk, Vzt, usum
      REAL (kind=8) :: bfl(2), bfg(2)

      REAL :: mctmp(nxp,nyp) ! Helper for mass concenrvation statistics

      ! Added by Juha
      SELECT CASE(level)
         CASE(1,2,3)
            rx = a_rv
         CASE(4,5)
            rx = a_rp
      END SELECT

      SELECT CASE(isfctyp)
            !
            ! set surface gradients
            !
         CASE(1)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = dthcon
                  drdz(i,j) = drtcon
                  bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
               END DO
            END DO
            zs = zt(2)/zrough
            CALL srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

            !
            ! get fluxes from profiles
            !
         CASE(2)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
            usum = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
                  drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
                  bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
                  usum = usum + a_ustar(i,j)
               END DO
            END DO
            usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
            zs = (zrough/100.)
            IF (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
            CALL srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
             !
             ! get fluxes from bulk formulae with coefficients given by dthcon and
             ! drtcon
             !
         CASE(3)
            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  dtdz(i,j) = a_theta(2,i,j) - sst*(p00/psrf)**rcp
                  drdz(i,j) = rx(2,i,j) - rslf(psrf,sst) ! Juha: rx
                  IF (ubmin > 0.) THEN
                     a_ustar(i,j) = sqrt(zrough)* wspd(i,j)
                  ELSE
                     a_ustar(i,j) = abs(ubmin)
                  END IF
                  a_tstar(i,j) =  dthcon * wspd(i,j)*dtdz(i,j)/a_ustar(i,j)
                  a_rstar(i,j) =  drtcon * wspd(i,j)*drdz(i,j)/a_ustar(i,j)
                  bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
               END DO
            END DO
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)
             !
             ! fix surface temperature to yield a constant surface buoyancy flux
             !
         CASE(4)

            Vzt   = 10.* (log(zt(2)/zrough)/log(10./zrough))
            Vbulk = Vzt * (vonk/log(zt(2)/zrough))**2

            bfl(:) = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  bfl(1) = bfl(1)+a_theta(2,i,j)
                  bfl(2) = bfl(2)+rx(2,i,j) ! Juha: rx
               END DO
            END DO

            CALL double_array_par_sum(bfl,bfg,2)

            bfg(2) = bfg(2)/REAL((nxpg-4)*(nypg-4))
            bfg(1) = bfg(1)/REAL((nxpg-4)*(nypg-4))

            DO iterate = 1, 5
               bflx  = ((sst -bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst) -bfg(2))) &
                        * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
               sst1 = sst + 0.1
               bflx1 = ((sst1-bfg(1)) + bfg(1)*ep2*(rslf(psrf,sst1)-bfg(2))) &
                        * 0.5*(dn0(1)+dn0(2))*cp*Vbulk
               sst  = sst + 0.1* (dthcon - bflx) / (bflx1-bflx)
            END DO

            DO j = 3, nyp-2
               DO i = 3, nxp-2
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
               END DO
            END DO

         CASE(5)

            !
            !Sami addition: Calculate surface fluxes from energy budget
            !
            !
            ! This is based on Acs et al: A coupled soil moisture and surface
            ! temperature prediction model, Journal of applied meteorology, 30, 1991
            !
            total_sw = 0.0
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
                  total_sw = total_sw+a_sflx(2,i,j)
                  !       WRITE(*,*) a_rflx(:,:,:)
                  total_rw = total_rw + a_rflx(2,i,j)
                  total_la = total_la + wq_sfc(i,j)*(0.5*(dn0(1)+dn0(2))*alvl)/ffact
                  total_se = total_se + wt_sfc(i,j)*(0.5*(dn0(1)+dn0(2))*cp)/ffact
                  total_pre=   total_pre +  precip(2,i,j)
               END DO
            END DO
            total_sw  = total_sw/REAL((nxp-4)*(nyp-4))
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

            Q12 = Kmean1*2.*( (fii_1-fii_2)/(D1+D2)+1.0)
            Q23 = Kmean2*2.*( (fii_2-fii_3)/(D2+D3)+1.0)

            IF (.NOT. lConstSoilWater) THEN
               W1 = W1+1./(thetaS1*D1)*(-Q12-((total_la+total_pre)/((0.5*(dn0(1)+dn0(2))*alvl)/ffact))/(thetaS1*1000))*dtl
               W2 = W2+1./(thetaS2*D2)*(Q12-Q23)*dtl
               W3 = W3+1./(thetaS3*D3)*(Q23-Q3)*dtl
            END IF
            !
            !  Following is copied from CASE (2). No idea if this is valid or not..
            !

            CALL get_swnds(nzp,nxp,nyp,usfc,vsfc,wspd,a_up,a_vp,umean,vmean)
            usum = 0.
            DO j = 3, nyp-2
               DO i = 3, nxp-2

                  dtdz(i,j) = a_theta(2,i,j) - sst1*(p00/psrf)**rcp

                  ff1 = 1.0
                  IF(W1 < 0.75) ff1 = W1/0.75
                  ! Flux of moisture is limited by water content.
                  drdz(i,j) = a_rp(2,i,j) - ff1*rslf(psrf,min(sst1,280.))  !  a_rv changed to a_rp (by Zubair)
                  !
                  bfct(i,j) = g*zt(2)/(a_theta(2,i,j)*wspd(i,j)**2)
                  usum = usum + a_ustar(i,j)
               END DO
            END DO
            usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
            zs = (zrough/100.)
            IF (zrough <= 0.) zs = max(0.0001,(0.016/g)*usum**2)
            CALL srfcscls(nxp,nyp,zt(2),zs,th00,wspd,dtdz,drdz,a_ustar,a_tstar,     &
                          a_rstar)
            CALL sfcflxs(nxp,nyp,vonk,wspd,usfc,vsfc,bfct,a_ustar,a_tstar,a_rstar,  &
                         uw_sfc,vw_sfc,wt_sfc,wq_sfc,ww_sfc)

           sst1 = sst1-(total_rw+total_la+total_se+( SQRT(lambda*C_heat*7.27e-5/(2.0)) *(SST1-deepSoilTemp)))&
                /(2.0e-2*C_heat+ SQRT( lambda*C_heat/(2.0*7.27e-5)) )*dtl

            sst = sst1
         !
         ! fix thermodynamic fluxes at surface given values in energetic
         ! units and calculate  momentum fluxes from winds
         !
         CASE DEFAULT
            ffact = 1.
            wt_sfc(1,1)  = ffact* dthcon/(0.5*(dn0(1)+dn0(2))*cp)
            wq_sfc(1,1)  = ffact* drtcon/(0.5*(dn0(1)+dn0(2))*alvl)

            IF (zrough <= 0.) THEN
               usum = 0.
               DO j = 3, nyp-2
                  DO i = 3, nxp-2
                     usum = usum + a_ustar(i,j)
                  END DO
               END DO
               usum = max(ubmin,usum/float((nxp-4)*(nyp-4)))
               zs = max(0.0001,(0.016/g)*usum**2)
            ELSE
               zs = zrough
            END IF

            DO j = 3, nyp-2
               DO i = 3, nxp-2
                  wt_sfc(i,j) = wt_sfc(1,1)
                  wq_sfc(i,j) = wq_sfc(1,1)

                  wspd(i,j)   = max(0.1,                                    &
                                sqrt((a_up(2,i,j)+umean)**2+(a_vp(2,i,j)+vmean)**2))
                  IF (ubmin > 0.) THEN
                     bflx = g*wt_sfc(1,1)/th00
                     IF (level >= 2) bflx = bflx + g*ep2*wq_sfc(i,j)
                     a_ustar(i,j) = diag_ustar(zt(2),zs,bflx,wspd(i,j))
                  ELSE
                     a_ustar(i,j) = abs(ubmin)
                  END IF

                  ffact = a_ustar(i,j)*a_ustar(i,j)/wspd(i,j)
                  uw_sfc(i,j)  = -ffact*(a_up(2,i,j)+umean)
                  vw_sfc(i,j)  = -ffact*(a_vp(2,i,j)+vmean)
                  ww_sfc(i,j)  = 0.
                  a_rstar(i,j) = wq_sfc(i,j)/a_ustar(i,j)
                  a_tstar(i,j) = wt_sfc(i,j)/a_ustar(i,j)
               END DO
            END DO

      END SELECT

      IF ( mcflg ) THEN
         !
         ! Juha: Take moisture flux to mass budged statistics
         mctmp(:,:) = wq_sfc(:,:)*(0.5*(a_dn(1,:,:)+a_dn(2,:,:)))
         CALL acc_massbudged(nzp,nxp,nyp,2,dtlt,dzt,a_dn,       &
                             revap=mctmp,ApVdom=mc_ApVdom)
         !
         !
      END IF !mcflg
      IF (sflg) CALL sfc_stat(nxp,nyp,wt_sfc,wq_sfc,a_ustar,sst)

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
   ! Layer, in Workshop on Micormeteorology, pages 67-100.
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
      REAL, INTENT(inout) :: ustar(n2,n3)  ! scale velocity
      REAL, INTENT(inout) :: tstar(n2,n3)  ! scale temperature
      REAL, INTENT(inout) :: rstar(n2,n3)  ! scale value of qt

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
               ustar(i,j) =  vonk*u(i,j)  /(lnz + am*zeta)
               tstar(i,j) = (vonk*dtv/(lnz + ah*zeta))/pr
               !
               ! Neutral case
               !
            ELSE IF (dtv == 0.) THEN
               ustar = vonk*u(i,j)  /lnz
               tstar = vonk*dtv/(pr*lnz)
               !
               ! ustable case, start iterations from values at previous tstep,
               ! unless the sign has changed or if it is the first call, then
               ! use neutral values.
               !
            ELSE
               IF (first_call .OR. tstar(i,j)*dtv <= 0.) THEN
                  ustar(i,j) = u(i,j)*klnz
                  tstar(i,j) = (dtv*klnz/pr)
               END IF

               DO iterate = 1, 3
                  lmo  = betg*ustar(i,j)**2/(vonk*tstar(i,j))
                  zeta = z/lmo
                  x    = sqrt(sqrt( 1.0 - bm*zeta ) )
                  psi1 = 2.*log(1.0+x) + log(1.0+x*x) - 2.*atan(x) + cnst1
                  y    = sqrt(1.0 - bh*zeta)
                  psi2 = log(1.0 + y) + cnst2
                  ustar(i,j) = u(i,j)*vonk/(lnz - psi1)
                  tstar(i,j) = (dtv*vonk/pr)/(lnz - psi2)
               END DO
            END IF

            rstar(i,j) = tstar(i,j)*drt(i,j)/(dtv + eps)
            tstar(i,j) = tstar(i,j)*dth(i,j)/(dtv + eps)
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
      REAL, INTENT(in)    :: us(n2,n3),ts(n2,n3),rs(n2,n3)
      REAL, INTENT(out)   :: uw(n2,n3),vw(n2,n3),tw(n2,n3),rw(n2,n3),ww(n2,n3)

      REAL :: x(n2,n3),y(n2,n3)
      INTEGER i,j

      DO j = 3, n3-2
         DO i = 3, n2-2

            uw(i,j) = -(u(i,j)/(ubar(i,j)+eps))*us(i,j)**2
            vw(i,j) = -(v(i,j)/(ubar(i,j)+eps))*us(i,j)**2
            tw(i,j) = -ts(i,j)*us(i,j)
            rw(i,j) = -rs(i,j)*us(i,j)

            x(i,j) = xx(i,j)*vk*ts(i,j)*(ubar(i,j)/us(i,j))**2
            x(i,j) = x(i,j)*sqrt(sqrt(1.-15.*min(0.,x(i,j)))) &
                     /(1.0+cc*max(0.,x(i,j)))
            y(i,j) = sqrt((1.-2.86*x(i,j))/(1.+x(i,j)* &
                     (-5.39+x(i,j)*6.998 )))
            ww(i,j)= (0.27*max(6.25*(1.-x(i,j))*y(i,j),eps)-&
                      1.18*x(i,j)*y(i,j))*us(i,j)**2
         END DO
      END DO
      RETURN
   END SUBROUTINE sfcflxs

END MODULE srfc



