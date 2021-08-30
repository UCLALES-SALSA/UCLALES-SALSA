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
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE mcrp

  ! Juha 2020 11 19: these imports should be sorted out, as well as the input/output parameters
  !                  in the subroutine calls... Comments and some changes below. Generally, if something
  !                  is given here in global scope, there is no need to provide them as inputs! 
  USE defs, ONLY : alvl,alvi,rowt,pi,Rm,cp,kb,g,vonk   ! These should be in module scope since used wherever
  USE mo_aux_state, ONLY : dzt                         ! Same for this
  USE grid, ONLY : dtlt,nxp,nyp,nzp,CCN,pbncsrc    ! Should be in module scope AND therefore no need to provide them as inputs!

  USE mo_submctl, ONLY : spec       ! Theres some silly overlaps with this that I did not yet resolve.
 
  !USE stat, ONLY : sflg, updtst, acc_removal, mcflg, acc_massbudged, cs_rem_set
  USE classProcessSwitch, ONLY : ProcessSwitch
  USE mo_structured_datatypes
  USE omp_lib
  IMPLICIT NONE
   
  INTEGER :: bulkScheme = 1   ! Select bulk microphysics parameterizations:
                              ! 1: Seifert and Beheng
                              ! 2: Khairoutdinov and Kogan 2000 (mean radius formulation for autoconv)
                              ! 3: Khairoutdinol and Kogan 2000 (rc nc exponential formulation)
   
  TYPE(ProcessSwitch) :: sed_aero,  &
                         sed_cloud, &
                         sed_precp, &
                         sed_ice,   &
                         bulk_autoc
  !
  ! drop sizes definition is based on vanZanten (2005)
  ! cloud droplets' diameter: 2-50 e-6 m
  ! drizzle drops' diameter: 50-1000 e-6 m
  !
  REAL, PARAMETER :: eps0 = 1e-20       ! small number
  REAL, PARAMETER :: eps1 = 1e-9        ! small number
  REAL, PARAMETER :: rho_0 = 1.21       ! air density at surface
  
  REAL, PARAMETER :: D_min = 2.e-6      ! minimum diameter of cloud droplets
  REAL, PARAMETER :: D_bnd = 80.e-6     ! precip/cloud diameter boundary
  REAL, PARAMETER :: D_max = 1.e-3      ! maximum diameter of prcp drops
  
  REAL, PARAMETER :: X_min = (D_min**3)*rowt*pi/6. ! min cld mass
  REAL, PARAMETER :: X_bnd = (D_bnd**3)*rowt*pi/6. ! precip/cld bound mass
  REAL, PARAMETER :: X_max = (D_max**3)*rowt*pi/6. ! max prcp mass
  
  REAL, PARAMETER :: prw = pi * rowt / 6.
  
  CONTAINS
  
    !
    ! ---------------------------------------------------------------------
    ! init_micro: initialize the sedimentation switches
    ! 
    SUBROUTINE init_mcrp_switches
      IMPLICIT NONE
      
      ! Defaults are FALSE
      sed_aero = ProcessSwitch()
      sed_cloud = ProcessSwitch()
      sed_precp = ProcessSwitch()
      sed_ice = ProcessSwitch()
      bulk_autoc = ProcessSwitch()
      
    END SUBROUTINE init_mcrp_switches
    
    !
    ! ---------------------------------------------------------------------
    ! MICRO: sets up call to microphysics
    !
    SUBROUTINE micro(level)
      USE mo_diag_state, ONLY : a_rv,a_rc,a_theta,     &
                                a_temp,a_rsl,a_dn,a_ustar,             &
                                a_rrate, a_irate, a_sfcrrate, a_sfcirate,   &
                                d_VtPrc, d_VtIce
      USE mo_progn_state, ONLY : a_rp,a_tp,a_rt,a_tt,a_rpp,a_rpt,a_npp,a_npt
      USE mo_aux_state, ONLY : dn0
      INTEGER, INTENT(in) :: level
      INTEGER :: nspec
      
      SELECT CASE (level)
      CASE(2)
         IF (sed_cloud%state)  &
              CALL sedim_cd(a_theta,a_temp,a_rc,dn0,a_rrate,a_rt,a_tt)
      CASE(3)
         CALL mcrph(dn0,a_theta,a_temp,a_rv,a_rsl,a_rc,a_rpp,       &
                    a_npp,a_rrate,a_sfcrrate,a_rt,a_tt,a_rpt,a_npt  )
      CASE(4,5)
         nspec = spec%getNSpec(type="wet")
         ! Import tracers directly and not via arguments because theres so many...
         CALL sedim_SALSA(nspec,level,a_ustar,a_temp,a_theta,a_dn,a_rrate,   &
                          a_sfcrrate,a_irate,a_sfcirate,d_VtPrc,d_VtIce,a_tt )
                         
      CASE(0) ! For piggybacking call to level 3 microphysics. pb_mcrph just wraps the necessary calls to thermo and mcrph
         CALL pb_mcrph(dn0)
      END SELECT
      
    END SUBROUTINE micro

    !
    !------------------------------------------------------------------------------------
    ! PB_MCRPH: Wraps the bulk microphysics in "piggybacking" setup,
    !           i.e. the bulk microphysics is calculated alongside SALSA, but
    !           arising tendencies in dynamical/thermodynamical variables are
    !           not coupled with the model core. The piggybacking microphysical
    !           properties and process rates will have a dedicated set of
    !           diagnostics
    !
    SUBROUTINE pb_mcrph(dn0)
      USE thrm, ONLY : thermo
      USE mo_progn_state, ONLY : pb_rpp, pb_npp,  &  ! Import the slave microphysics variables here locally
                                 pb_rpt, pb_npt
      USE mo_diag_state, ONLY : pb_rv, pb_rc, pb_rh,   &
                                pb_theta, pb_temp,     &
                                pb_rsl, pb_rrate,      &
                                pb_sfcrrate
      USE mo_derived_state, ONLY : CDNC
      USE mo_ts_state, ONLY : tsic_CDNC
      TYPE(FloatArray1d), INTENT(in) :: dn0
      REAL :: strg(nzp,nxp,nyp,2)  ! dummy storage array
      REAL :: pbcdnc(nzp,nxp,nyp)
      TYPE(FloatArray3d) :: rttdummy, tltdummy

      ! These are here just because mcrph requires these tendencies as output arguments,
      ! but in the piggybacking mode we will want to neglect them. Check later if possible
      ! to avoid this mess. Associating and nullifying the pointers every timestep may be
      ! a slight performance issue.
      strg = 0.
      rttdummy = FloatArray3d(strg(:,:,:,1))
      tltdummy = FloatArray3d(strg(:,:,:,2))
      
      ! Call with level = 0 updates the diagnostic pb-variables using saturation adjustment to be
      ! consistent with the level 3 microphysics used as the "slave" scheme
      CALL thermo(0)

      ! Diagnose cdnc from SALSA
      pbcdnc = 0.
      SELECT CASE(pbncsrc)
      CASE(1)
         ! Use master CDNC gridpoint by gridpoint
         CALL CDNC%onDemand("CDNC",pbcdnc)
      CASE(2)
         ! Use the domain mean of the master CDNC.
         ! For this, modify the parameter CCN.
         ! (note that CCN is also used in radiation, but only for levels 1-3,
         ! so this should not have any side effects)
         CALL tsic_CDNC%onDemand(CCN,root=.FALSE.)
      END SELECT
         
      ! A regular call to level 3 microphysics, except now using the pb-variables and neglecting
      ! the contributions to total water mix rat and temperature tendencies.
      IF (pbncsrc == 1) THEN
         ! Use the master CDNC grid-point by grid-point
         CALL mcrph(dn0,pb_theta,pb_temp,pb_rv,pb_rsl,pb_rc,pb_rpp,    &      
                    pb_npp,pb_rrate,pb_sfcrrate,rttdummy,tltdummy,     &
                    pb_rpt,pb_npt,pbcdnc)
      ELSE
         ! Use constant CCN or the mean of master CDNC (assigned to the variable CCN)
         CALL mcrph(dn0,pb_theta,pb_temp,pb_rv,pb_rsl,pb_rc,pb_rpp,    &      
                    pb_npp,pb_rrate,pb_sfcrrate,rttdummy,tltdummy,     &
                    pb_rpt,pb_npt)         
      END IF
      
      rttdummy%d => NULL()
      tltdummy%d => NULL()
    END SUBROUTINE pb_mcrph    
    !
    ! ---------------------------------------------------------------------
    ! MCRPH: calls microphysical parameterization
    !
    SUBROUTINE mcrph(dn0,th,tk,rv,rs,rc,rp,np,rrate,         &
                     sfcrrate,rtt,tlt,rpt,npt,pcdnc)
      
      TYPE(FloatArray3d), INTENT (in)    :: th, tk, rv, rs
      TYPE(FloatArray1d), INTENT (in)    :: dn0
      TYPE(FloatArray3d), INTENT (inout) :: rc, rtt, tlt, rpt, npt, np, rp
      TYPE(FloatArray3d), INTENT (inout) :: rrate
      TYPE(FloatArray2d), INTENT (inout) :: sfcrrate
      REAL, INTENT(in), OPTIONAL         :: pcdnc(nzp,nxp,nyp) ! For piggybacking runs, thus optional.
                                                              ! Diagnosed locally from SALSA, so datatype
                                                              ! is a simple REAL array
      
      INTEGER :: i, j, k
      
      !
      ! Microphysics following Seifert Beheng (2001, 2005)
      ! note that the order below is important as the rc array is
      ! redefined in cld_dgn below and is assumed to be cloud water
      ! after that and total condensate priort to that
      !
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 1, nzp
               rp%d(k,i,j) = max(0., rp%d(k,i,j))
               np%d(k,i,j) = max(min(rp%d(k,i,j)/X_bnd,np%d(k,i,j)),rp%d(k,i,j)/X_max)
            END DO
         END DO
      END DO
      
      CALL wtr_dff_SB(dn0,rp,np,rc,rs,rv,tk,rpt,npt)
      
      IF (bulk_autoc%state) THEN
         SELECT CASE(bulkScheme)
            CASE(1)
               CALL auto_SB(dn0,rc,rp,rpt,npt,pcdnc)
               CALL accr_SB(dn0,rc,rp,np,rpt,npt)
            CASE(2)
               CALL auto_KK_1(rc,rpt,npt,pcdnc)
               CALL accr_KK(dn0,rc,rp,np,rpt,npt)
            CASE(3)
               CALL auto_KK_2(rc,rpt,npt,pcdnc)
               CALL accr_KK(dn0,rc,rp,np,rpt,npt)
            END SELECT
      END IF
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               rp%d(k,i,j)  = rp%d(k,i,j) + max(-rp%d(k,i,j)/dtlt,rpt%d(k,i,j))*dtlt
               np%d(k,i,j)  = np%d(k,i,j) + max(-np%d(k,i,j)/dtlt,npt%d(k,i,j))*dtlt
               rpt%d(k,i,j) = 0.
               npt%d(k,i,j) = 0.
               rp%d(k,i,j)  = max(0., rp%d(k,i,j))
               np%d(k,i,j)  = max(min(rp%d(k,i,j)/X_bnd,np%d(k,i,j)),rp%d(k,i,j)/X_max)
            END DO
         END DO
      END DO

      ! zero the rate diagnostic for current timestep
      rrate%d = 0.
      sfcrrate%d = 0.
      
      IF (sed_precp%state) CALL sedim_rd(dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)
      
      IF (sed_cloud%state) CALL sedim_cd(th,tk,rc,dn0,rrate,rtt,tlt,pcdnc)

      sfcrrate%d(:,:) = rrate%d(2,:,:)
      
    END SUBROUTINE mcrph
    !
    ! ---------------------------------------------------------------------
    ! WTR_DFF_SB: calculates the evolution of the both number- and
    ! mass mixing ratio large drops due to evaporation in the absence of
    ! cloud water.
    !
    SUBROUTINE wtr_dff_SB(dn0,rp,np,rl,rs,rv,tk,rpt,npt)
      TYPE(FloatArray3d), INTENT(in)    :: tk, np, rp, rs, rv
      TYPE(FloatArray1d), INTENT(in)    :: dn0
      TYPE(FloatArray3d), INTENT (inout):: rpt, npt, rl
      
      REAL, PARAMETER    :: Kt = 2.5e-2    ! conductivity of heat [J/(sKm)]
      REAL, PARAMETER    :: Dv = 3.e-5     ! diffusivity of water vapor [m2/s]
      
      INTEGER             :: i, j, k
      REAL                :: Xp, Dp, G, S, cerpt, cenpt
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp
               IF (rp%d(k,i,j) > rl%d(k,i,j)) THEN   
                  Xp = rp%d(k,i,j)/ (np%d(k,i,j)+eps0)
                  Xp = MIN(MAX(Xp,X_bnd),X_max)
                  Dp = ( Xp / prw )**(1./3.)
                  
                  G = 1. / (1. / (dn0%d(k)*rs%d(k,i,j)*Dv) + &
                       alvl*(alvl/(Rm*tk%d(k,i,j))-1.) / (Kt*tk%d(k,i,j)))
                  S = rv%d(k,i,j)/rs%d(k,i,j) - 1.
                  
                  IF (S < 0) THEN
                     cerpt = 2. * pi * Dp * G * S * np%d(k,i,j)
                     cenpt = cerpt * np%d(k,i,j) / rp%d(k,i,j)
                     rpt%d(k,i,j) = rpt%d(k,i,j) + cerpt
                     npt%d(k,i,j) = npt%d(k,i,j) + cenpt
                  END IF
               END IF
               rl%d(k,i,j) = max(0.,rl%d(k,i,j) - rp%d(k,i,j))  ! What 
            END DO
         END DO
      END DO
            
    END SUBROUTINE wtr_dff_SB
    !
    ! ---------------------------------------------------------------------
    ! AUTO_SB:  calculates the evolution of mass- and number mxg-ratio for
    ! drizzle drops due to autoconversion. The autoconversion rate assumes
    ! f(x)=A*x**(nu_c)*exp(-Bx), an exponential in drop MASS x. It can
    ! be reformulated for f(x)=A*x**(nu_c)*exp(-Bx**(mu)), where formu=1/3
    ! one would get a gamma dist in drop diam -> faster rain formation.
    !
    SUBROUTINE auto_SB(dn0,rc,rp,rpt,npt,pcdnc)
      USE mo_diag_state, ONLY : b_m_autoc, b_n_autoc  ! for rate diagnostics
      
      TYPE(FloatArray1d), INTENT(in)    :: dn0
      TYPE(FloatArray3d), INTENT(in)    :: rc, rp
      TYPE(FloatArray3d), INTENT(inout) :: rpt, npt
      REAL, INTENT(in), OPTIONAL        :: pcdnc(nzp,nxp,nyp) ! For piggybacking runs, thus optional.
                                                             ! Diagnosed locally from SALSA, so datatype
                                                             ! is a simple REAL array
      
      REAL, PARAMETER :: nu_c = 0            ! width parameter of cloud DSD
      REAL, PARAMETER :: k_c  = 9.44e+9      ! Long-Kernel
      REAL, PARAMETER :: k_1  = 6.e+2        ! Parameter for phi function
      REAL, PARAMETER :: k_2  = 0.68         ! Parameter for phi function

      REAL, PARAMETER :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
      REAL, PARAMETER :: Eau = 5.67    ! autoconv. exponent in KK param.
      REAL, PARAMETER :: mmt = 1.e+6   ! transformation from m to \mu m

      INTEGER :: i, j, k
      REAL    :: k_au, Xc, Dc, au, tau, phi

      b_m_autoc%d = 0.
      b_n_autoc%d = 0.
      
      k_au = k_c / (20.*X_bnd) * (nu_c+2.)*(nu_c+4.)/(nu_c+1.)**2

      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               IF (PRESENT(pcdnc)) THEN
                  Xc = rc%d(k,i,j)/(pcdnc(k,i,j)+eps0)
               ELSE
                  Xc = rc%d(k,i,j)/(CCN+eps0)
               END IF
               IF (Xc > 0.) THEN
                  Xc = MIN(MAX(Xc,X_min),X_bnd)
                  au = k_au * dn0%d(k) * rc%d(k,i,j)**2 * Xc**2
                  !
                  ! small threshold that should not influence the result
                  !
                  IF (rc%d(k,i,j) > 1.e-6) THEN
                     tau = 1.0-rc%d(k,i,j)/(rc%d(k,i,j)+rp%d(k,i,j)+eps0)
                     tau = MIN(MAX(tau,eps0),0.9)
                     phi = k_1 * tau**k_2 * (1.0 - tau**k_2)**3
                     au  = au * (1.0 + phi/(1.0 - tau)**2)
                  END IF

                  rpt%d(k,i,j) = rpt%d(k,i,j) + au
                  npt%d(k,i,j) = npt%d(k,i,j) + au/X_bnd                  
                  !
                  ! rate diagnostics
                  b_m_autoc%d(k,i,j) = au
                  b_n_autoc%d(k,i,j) = au/X_bnd
                  !
               END IF
            END DO
         END DO
      END DO
      
    END SUBROUTINE auto_SB

    !
    !---------------------------------
    !auto_KK_1: The Khairoutdinov and Kogan (2000) autoconversion
    !           parameterization (the volume mean radius formulation).
    !
    SUBROUTINE auto_KK_1(rc,rpt,npt,pcdnc)
      USE mo_diag_state, ONLY : b_m_autoc, b_n_autoc  ! for rate diagnostics

      TYPE(FloatArray3d), INTENT(in)    :: rc
      TYPE(FloatArray3d), INTENT(inout) :: rpt, npt
      REAL, INTENT(in), OPTIONAL        :: pcdnc(nzp,nxp,nyp) ! For piggybacking runs, thus optional.
                                                             ! Diagnosed locally from SALSA, so datatype
                                                             ! is a simple REAL array

      REAL, PARAMETER :: Cau = 4.1e-15 ! autoconv. coefficient in KK param.
      REAL, PARAMETER :: Eau = 5.67    ! autoconv. exponent in KK param.
      REAL, PARAMETER :: mmt = 1.e+6   ! transformation from m to \mu m

      INTEGER :: i, j, k
      REAL    :: Xc, Dc, au
      
      b_m_autoc%d = 0.
      b_n_autoc%d = 0.
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               IF (PRESENT(pcdnc)) THEN
                  Xc = rc%d(k,i,j)/(pcdnc(k,i,j)+eps0)
               ELSE
                  Xc = rc%d(k,i,j)/(CCN+eps0)
               END IF
               IF (Xc > 0.) THEN
                  Dc = ( Xc / prw )**(1./3.)
                  au = Cau * (Dc * mmt / 2.)**Eau
                  !
                  rpt%d(k,i,j) = rpt%d(k,i,j) + au
                  npt%d(k,i,j) = npt%d(k,i,j) + au/X_bnd                  
                  !
                  ! rate diagnostics
                  b_m_autoc%d(k,i,j) = au
                  b_n_autoc%d(k,i,j) = au/X_bnd

               END IF
            END DO
         END DO
      END DO
      
    END SUBROUTINE auto_KK_1
    !
    !---------------------------------
    !auto_KK_2: The Khairoutdinov and Kogan (2000) autoconversion
    !           parameterization (the rc nc formulation).
    !
    SUBROUTINE auto_KK_2(rc,rpt,npt,pcdnc)
      USE mo_diag_state, ONLY : b_m_autoc, b_n_autoc  ! for rate diagnostics

      TYPE(FloatArray3d), INTENT(in)    :: rc
      TYPE(FloatArray3d), INTENT(inout) :: rpt, npt
      REAL, INTENT(in), OPTIONAL        :: pcdnc(nzp,nxp,nyp) ! For piggybacking runs, thus optional.
                                                             ! Diagnosed locally from SALSA, so datatype
                                                             ! is a simple REAL array

      REAL, PARAMETER :: Cau = 1350. ! autoconv. coefficient in KK param.
      REAL, PARAMETER :: E1 = 2.47    ! autoconv. exponent in KK param.
      REAL, PARAMETER :: E2 = -1.79   ! autoconv. exponent in KK param
      REAL, PARAMETER :: pcm = 1.e-6

      INTEGER :: i, j, k
      REAL    :: nc, Dc, au
      
      b_m_autoc%d = 0.
      b_n_autoc%d = 0.
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               IF (PRESENT(pcdnc)) THEN
                  nc = pcdnc(k,i,j)
               ELSE
                  nc = CCN
               END IF
              
               IF (nc > 0. .AND. rc%d(k,i,j) > 1.e-6) THEN
                  au = Cau * (rc%d(k,i,j)**E1) * ((nc*pcm)**E2) 
                  !
                  rpt%d(k,i,j) = rpt%d(k,i,j) + au
                  npt%d(k,i,j) = npt%d(k,i,j) + au/X_bnd                  
                  !
                  ! rate diagnostics
                  b_m_autoc%d(k,i,j) = au
                  b_n_autoc%d(k,i,j) = au/X_bnd

               END IF
            END DO
         END DO
      END DO
      
    END SUBROUTINE auto_KK_2
    
    !
    ! ---------------------------------------------------------------------
    ! ACCR_SB calculates the evolution of mass mxng-ratio due to accretion
    ! and self collection following Seifert & Beheng (2001).  Included is
    ! an alternative formulation for accretion only, following
    ! Khairoutdinov and Kogan
    !
    SUBROUTINE accr_SB(dn0,rc,rp,np,rpt,npt)
      USE mo_diag_state, ONLY : b_m_accr
      TYPE(FloatArray3d), INTENT(in)    :: rc, rp, np
      TYPE(FloatArray1d), INTENT(in)    :: dn0
      TYPE(FloatArray3d), INTENT(inout) :: rpt, npt

      REAL, PARAMETER :: k_r = 5.78
      REAL, PARAMETER :: k_1 = 5.e-4
      REAL, PARAMETER :: Cac = 67.     ! accretion coefficient in KK param.
      REAL, PARAMETER :: Eac = 1.15    ! accretion exponent in KK param.

      INTEGER :: i, j, k
      REAL    :: tau, phi, ac, sc

      b_m_accr%d = 0.
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               IF (rc%d(k,i,j) > 0. .AND. rp%d(k,i,j) > 0.) THEN
                  tau = 1.0-rc%d(k,i,j)/(rc%d(k,i,j)+rp%d(k,i,j)+eps0)
                  tau = MIN(MAX(tau,eps0),1.)
                  phi = (tau/(tau+k_1))**4
                  ac  = k_r * rc%d(k,i,j) * rp%d(k,i,j) * phi * sqrt(rho_0*dn0%d(k))

                  rpt%d(k,i,j) = rpt%d(k,i,j) + ac
                  !
                  ! Rate diagnostic
                  b_m_accr%d(k,i,j) = ac
                  !
               END IF
               sc = k_r * np%d(k,i,j) * rp%d(k,i,j) * sqrt(rho_0*dn0%d(k))
               npt%d(k,i,j) = npt%d(k,i,j) - sc
            END DO
         END DO
      END DO
      
    END SUBROUTINE accr_SB
    
    !----------------------------------------------------------------------
    ! ACCR_KK: accretion rate according to Khairoutdinov and Kogan (2000)
    !
    SUBROUTINE accr_KK(dn0,rc,rp,np,rpt,npt)
      USE mo_diag_state, ONLY : b_m_accr
      TYPE(FloatArray3d), INTENT(in)    :: rc, rp, np
      TYPE(FloatArray1d), INTENT(in)    :: dn0
      TYPE(FloatArray3d), INTENT(inout) :: rpt, npt

      REAL, PARAMETER :: k_r = 5.78
      REAL, PARAMETER :: k_1 = 5.e-4
      REAL, PARAMETER :: Cac = 67.     ! accretion coefficient in KK param.
      REAL, PARAMETER :: Eac = 1.15    ! accretion exponent in KK param.

      INTEGER :: i, j, k
      REAL    :: tau, phi, ac, sc

      b_m_accr%d = 0.
      
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            DO k = 2, nzp-1
               IF (rc%d(k,i,j) > 0. .AND. rp%d(k,i,j) > 0.) THEN

                  ac = Cac * (rc%d(k,i,j) * rp%d(k,i,j))**Eac 
                  !
                  rpt%d(k,i,j) = rpt%d(k,i,j) + ac
                  !
                  ! Rate diagnostic
                  b_m_accr%d(k,i,j) = ac
                  !
               END IF
               sc = k_r * np%d(k,i,j) * rp%d(k,i,j) * sqrt(rho_0*dn0%d(k))
               npt%d(k,i,j) = npt%d(k,i,j) - sc
            END DO
         END DO
      END DO
      
    END SUBROUTINE accr_KK
    
    !
    ! ---------------------------------------------------------------------
    ! SEDIM_RD: calculates the sedimentation of the rain drops and its
    ! effect on the evolution of theta_l and r_t.  This is expressed in
    ! terms of Dp the mean diameter, not the mass weighted mean diameter
    ! as is used elsewhere.  This is just 1/lambda in the exponential
    ! distribution
    !
    SUBROUTINE sedim_rd(dn0,rp,np,tk,th,rrate,rtt,tlt,rpt,npt)
      TYPE(FloatArray1d), INTENT (in)           :: dn0
      TYPE(FloatArray3d), INTENT (in)           :: rp, np, th, tk
      TYPE(FloatArray3d), INTENT (inout)        :: rrate
      TYPE(FloatArray3d), INTENT (inout)        :: rtt, tlt, rpt, npt

      REAL, PARAMETER :: a2 = 9.65       ! in SI [m/s]
      REAL, PARAMETER :: c2 = 6e2        ! in SI [1/m]
      REAL, PARAMETER :: Dv = 25.0e-6    ! in SI [m/s]
      REAL, PARAMETER :: cmur1 = 10.0    ! mu-Dm-relation for rain following
      REAL, PARAMETER :: cmur2 = 1.20e+3 ! Milbrandt&Yau 2005, JAS, but with
      REAL, PARAMETER :: cmur3 = 1.5e-3  ! revised constants
      REAL, PARAMETER :: aq = 6.0e3
      REAL, PARAMETER :: bq = -0.2
      REAL, PARAMETER :: an = 3.5e3
      REAL, PARAMETER :: bn = -0.1

      INTEGER :: i, j, k, kp1, kk, km1
      REAL    :: b2, Xp, Dp, Dm, mu, flxdiv, tot,sk, mini, maxi, cc, zz
      REAL, DIMENSION(nzp) :: nslope,rslope,dn,dr, rfl, nfl, vn, vr, cn, cr

      b2 = a2*exp(c2*Dv)

      DO j = 3, nyp-2
         DO i = 3, nxp-2

            nfl(nzp) = 0.
            rfl(nzp) = 0.
            DO k = nzp-1, 2, -1
               Xp = rp%d(k,i,j) / (np%d(k,i,j)+eps0)
               Xp = MIN(MAX(Xp,X_bnd),X_max)
               !
               ! Adjust Dm and mu-Dm and Dp=1/lambda following Milbrandt & Yau
               !
               Dm = ( 6. / (rowt*pi) * Xp )**(1./3.)
               mu = cmur1*(1.+tanh(cmur2*(Dm-cmur3)))
               Dp = (Dm**3/((mu+3.)*(mu+2.)*(mu+1.)))**(1./3.)

               
               SELECT CASE(bulkScheme)
               CASE(1)
                  ! Seifert and Beheng
                  vn(k) = sqrt(dn0%d(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(1.+mu)))
                  vr(k) = sqrt(dn0%d(k)/1.2)*(a2 - b2*(1.+c2*Dp)**(-(4.+mu)))
               CASE(2,3)
                  !Khairoutdinov and Kogan 2000
                  vn(k) = max(0.,an * Dp + bn)
                  vr(k) = max(0.,aq * Dp + bq)
               END SELECT

            END DO

            DO k = 2, nzp-1
               kp1 = min(k+1,nzp-1)
               km1 = max(k,2)
               cn(k) = 0.25*(vn(kp1)+2.*vn(k)+vn(km1))*dzt%d(k)*dtlt
               cr(k) = 0.25*(vr(kp1)+2.*vr(k)+vr(km1))*dzt%d(k)*dtlt
            END DO

            !...piecewise linear method: get slopes
            DO k = nzp-1, 2, -1
               dn(k) = np%d(k+1,i,j)-np%d(k,i,j)
               dr(k) = rp%d(k+1,i,j)-rp%d(k,i,j)
            END DO
            dn(1)  = dn(2)
            dn(nzp) = dn(nzp-1)
            dr(1)  = dr(2)
            dr(nzp) = dr(nzp-1)
            DO k = nzp-1, 2, -1
               !...slope with monotone limiter for np
               sk = 0.5 * (dn(k-1) + dn(k))
               mini = min(np%d(k-1,i,j),np%d(k,i,j),np%d(k+1,i,j))
               maxi = max(np%d(k-1,i,j),np%d(k,i,j),np%d(k+1,i,j))
               nslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(np%d(k,i,j)-mini), &
                  &                              2.*(maxi-np%d(k,i,j)))
               !...slope with monotone limiter for rp
               sk = 0.5 * (dr(k-1) + dr(k))
               mini = min(rp%d(k-1,i,j),rp%d(k,i,j),rp%d(k+1,i,j))
               maxi = max(rp%d(k-1,i,j),rp%d(k,i,j),rp%d(k+1,i,j))
               rslope(k) = 0.5 * sign(1.,sk)*min(abs(sk), 2.*(rp%d(k,i,j)-mini), &
                  &                              2.*(maxi-rp%d(k,i,j)))
            END DO

            rfl(nzp-1) = 0.
            nfl(nzp-1) = 0.
            DO k = nzp-2, 2, -1

               kk = k
               tot = 0.0
               zz  = 0.0
               cc  = min(1.,cn(k))
               DO WHILE (cc > 0 .AND. kk <= nzp-1)
                  tot = tot + dn0%d(kk)*(np%d(kk,i,j)+nslope(kk)*(1.-cc))*cc/dzt%d(kk)
                  zz  = zz + 1./dzt%d(kk)
                  kk  = kk + 1
                  cc  = min(1.,cn(kk) - zz*dzt%d(kk))
               END DO
               nfl(k) = -tot /dtlt

               kk = k
               tot = 0.0
               zz  = 0.0
               cc  = min(1.,cr(k))
               DO WHILE (cc > 0 .AND. kk <= nzp-1)
                  tot = tot + dn0%d(kk)*(rp%d(kk,i,j)+rslope(kk)*(1.-cc))*cc/dzt%d(kk)
                  zz  = zz + 1./dzt%d(kk)
                  kk  = kk + 1
                  cc  = min(1.,cr(kk) - zz*dzt%d(kk))
               END DO
               rfl(k) = -tot /dtlt ! kg/m2/s?

               kp1 = k+1
               flxdiv = (rfl(kp1)-rfl(k))*dzt%d(k)/dn0%d(k)  ! kg/kg/s
               rpt%d(k,i,j) =rpt%d(k,i,j)-flxdiv
               rtt%d(k,i,j) =rtt%d(k,i,j)-flxdiv
               tlt%d(k,i,j) =tlt%d(k,i,j)+flxdiv*(alvl/cp)*th%d(k,i,j)/tk%d(k,i,j)

               npt%d(k,i,j) = npt%d(k,i,j)-(nfl(kp1)-nfl(k))*dzt%d(k)/dn0%d(k)

               rrate%d(k,i,j) = rrate%d(k,i,j) - rfl(k)/dn0%d(k) * alvl*0.5*(dn0%d(k)+dn0%d(kp1))
            END DO
         END DO
      END DO
      
    END SUBROUTINE sedim_rd
    ! 
    ! ---------------------------------------------------------------------
    ! SEDIM_CD: calculates the cloud-droplet sedimentation flux and its effect
    ! on the evolution of r_t and theta_l assuming a log-normal distribution
    ! 
    SUBROUTINE sedim_cd(th,tk,rc,dn0,rrate,rtt,tlt,pcdnc)      
      TYPE(FloatArray3d), INTENT (in)         :: th,tk,rc
      TYPE(FloatArray1d), INTENT (in)         :: dn0
      TYPE(FloatArray3d), INTENT (inout)      :: rrate
      TYPE(FloatArray3d), INTENT (inout)      :: rtt,tlt
      REAL, INTENT(in), OPTIONAL              :: pcdnc(nzp,nxp,nyp) ! For piggybacking runs, thus optional.
                                                                   ! Diagnosed locally from SALSA, so datatype
                                                                   ! is a simple REAL array
      
      REAL, PARAMETER :: c = 1.19e8 ! Stokes fall velocity coef [m^-1 s^-1]
      REAL, PARAMETER :: sgg = 1.2  ! geometric standard dev of cloud droplets

      INTEGER :: i, j, k, kp1
      REAL    :: Dc, Xc, vc, flxdiv
      REAL    :: rfl(nzp)

      ! 
      ! calculate the precipitation flux and its effect on r_t and theta_l
      ! 
      
      !!$omp parallel
      !!$omp do firstprivate(Xc, Dc, rfl, kp1,flxdiv)
      DO j = 3, nyp-2
         DO i = 3, nxp-2
            rfl(nzp) = 0.
            DO k = nzp-1, 2, -1
               IF (PRESENT(pcdnc)) THEN
                  Xc = rc%d(k,i,j) / (pcdnc(k,i,j)+eps0)
               ELSE
                  Xc = rc%d(k,i,j) / (CCN+eps0)
               END IF
               Dc = ( Xc / prw )**((1./3.))
               Dc = MIN(MAX(Dc,D_min),D_bnd)
               vc = min(c*(Dc*0.5)**2 * exp(4.5*(log(sgg))**2),1./(dzt%d(k)*dtlt))
               rfl(k) = - rc%d(k,i,j) * vc  ! kg/kg m/s
               
               kp1 = k+1
               flxdiv = (rfl(kp1)-rfl(k))*dzt%d(k)  ! kg/kg/s
               rtt%d(k,i,j) = rtt%d(k,i,j)-flxdiv
               tlt%d(k,i,j) = tlt%d(k,i,j)+flxdiv*(alvl/cp)*th%d(k,i,j)/tk%d(k,i,j)
               rrate%d(k,i,j) = rrate%d(k,i,j) - rfl(k) * alvl * 0.5*(dn0%d(k) + dn0%d(k+1))               
            END DO
         END DO
      END DO
         
            
    END SUBROUTINE sedim_cd
    
    ! ---------------------------------------------------------------------
    ! SEDIM_AERO: calculates the salsa particles sedimentation and dry deposition flux  (.. Zubair) !
    !
    ! Juha: The code below is a modified version of the original one by Zubair
    !
    ! Juha: Rain is now treated completely separately (20151013)
    !
    ! Jaakko: Modified for the use of ice and snow bins
    SUBROUTINE sedim_SALSA(nspec,level,ustar,tk,th,adn,rrate,sfcrrate,    &
                           irate,sfcirate,VtPrc,VtIce,tlt                 )
      USE util, ONLY : getMassIndex
      USE mo_progn_state, ONLY :  a_naerop,  a_naerot,  a_maerop,  a_maerot,         &
                                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,        &
                                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,        &
                                  a_nicep,   a_nicet,   a_micep,   a_micet
      !USE grid, ONLY : mc_ApVdom
      USE mo_submctl, ONLY : nbins, ncld, nprc,           &
                             nice

      IMPLICIT NONE

      INTEGER, INTENT(in) :: nspec,level
      TYPE(FloatArray3d), INTENT(in) :: tk, th, adn
      TYPE(FloatArray2d), INTENT(in) :: ustar
      TYPE(FloatArray3d), INTENT(inout) :: tlt
      TYPE(FloatArray3d), INTENT(inout) :: rrate,irate
      TYPE(FloatArray2d), INTENT(inout) :: sfcrrate, sfcirate
      TYPE(FloatArray4d), INTENT(inout) :: VtPrc,VtIce
      
      INTEGER :: i,j,k,nc,istr,iend

      REAL :: prnt(nzp,nxp,nyp,nprc), prmt(nzp,nxp,nyp,nspec*nprc)     ! Number and mass tendencies due to fallout, precip
      REAL :: irnt(nzp,nxp,nyp,nice), irmt(nzp,nxp,nyp,(nspec+1)*nice)     ! Number and mass tendencies due to fallout, ice

      ! PArticle removal arrays, given in kg/(m2 s)
      REAL :: remaer(nxp,nyp,nspec*nbins),   &
              remcld(nxp,nyp,nspec*ncld),    &
              remprc(nxp,nyp,nspec*nprc),    &
              remice(nxp,nyp,(nspec+1)*nice)   ! nspec + 1 because of RIME

      ! Particle number removal arrays
      REAL :: andep(nxp,nyp,nbins),     &
              cndep(nxp,nyp,ncld)

      !REAL :: mctmp(n2,n3) ! Helper for mass conservation calculations

      ! Divergence fields
      REAL :: amdiv(nzp,nxp,nyp,nspec*nbins),    &
              cmdiv(nzp,nxp,nyp,nspec*ncld)
      REAL :: andiv(nzp,nxp,nyp,nbins),       &
              cndiv(nzp,nxp,nyp,ncld)

      remaer = 0.; remcld = 0.; remprc = 0.; remice = 0.
      prnt = 0.; prmt = 0.; irnt = 0.; irmt = 0.

      rrate%d = 0.; sfcrrate%d = 0.
      IF (level == 5) THEN
         irate%d = 0.; sfcirate%d = 0.
      END IF
      
      ! Sedimentation for slow (non-precipitating) particles
      !-------------------------------------------------------
      IF (sed_aero%state) THEN

         CALL DepositionSlow(nbins,nspec,tk,adn,ustar,a_naerop,a_maerop, &
                             andiv,amdiv,andep,remaer,rrate,sfcrrate,1     )
         
         a_naerot%d = a_naerot%d - andiv
         a_maerot%d = a_maerot%d - amdiv
         
         ! Account for changes in liquid water pot temperature
         nc = spec%getIndex('H2O')
         istr = getMassIndex(nbins,1,nc)
         iend = getMassIndex(nbins,nbins,nc)
         DO j = 3,nyp-2
            DO i = 3,nxp-2
               DO k = 2,nzp
                  tlt%d(k,i,j) = tlt%d(k,i,j) +      &
                       SUM(amdiv(k,i,j,istr:iend))*(alvl/cp)*th%d(k,i,j)/tk%d(k,i,j)
               END DO
            END DO
         END DO
         
      END IF ! sed_aero
      
      IF (sed_cloud%state) THEN
         
         CALL DepositionSlow(ncld,nspec,tk,adn,ustar,a_ncloudp,a_mcloudp, &
                             cndiv,cmdiv,cndep,remcld,rrate,sfcrrate,2     )
         
         a_ncloudt%d = a_ncloudt%d - cndiv
         a_mcloudt%d = a_mcloudt%d - cmdiv
         
         ! Account for changes in liquid water pot temperature
         nc = spec%getIndex('H2O')
         istr = getMassIndex(ncld,1,nc)
         iend = getMassIndex(ncld,ncld,nc)
         DO j = 3,nyp-2
            DO i = 3,nxp-2
               DO k = 2,nzp
                  tlt%d(k,i,j) = tlt%d(k,i,j) +      &
                       SUM(cmdiv(k,i,j,istr:iend))*(alvl/cp)*th%d(k,i,j)/tk%d(k,i,j)
               END DO
            END DO
         END DO
                  
      END IF ! sed_cloud
      
      ! ---------------------------------------------------------
      ! SEDIMENTATION/DEPOSITION OF FAST PRECIPITATING PARTICLES
      IF (sed_precp%state) THEN
         
         CALL DepositionFast(nprc,nspec,tk,adn,a_nprecpp,a_mprecpp,   &
                             prnt,prmt,remprc,rrate,sfcrrate,VtPrc,3  )
         
         a_nprecpt%d(:,:,:,:) = a_nprecpt%d(:,:,:,:) + prnt(:,:,:,:)/dtlt
         a_mprecpt%d(:,:,:,:) = a_mprecpt%d(:,:,:,:) + prmt(:,:,:,:)/dtlt
         nc = spec%getIndex('H2O')
         ! Surface precipitation W m-2
         !istr = getMassIndex(nprc,1,nc); iend = getMassIndex(nprc,nprc,nc)
         !sfcrrate%d(:,:) = SUM(remprc(:,:,istr:iend),DIM=3)*alvl
         
         ! Convert mass flux to heat flux (W/m^2) for output
         !rrate%d(:,:,:) = rrate%d(:,:,:)*alvl

         istr = getMassIndex(nprc,1,nc)
         iend = getMassIndex(nprc,nprc,nc)
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 1, nzp-1
                  tlt%d(k,i,j) = tlt%d(k,i,j) -           &
                       (SUM(prmt(k,i,j,istr:iend))/dtlt)*(alvl/cp)*th%d(k,i,j)/tk%d(k,i,j)
               END DO
            END DO
         END DO
      END IF
      
      IF (sed_ice%state .AND. level == 5) THEN                          
         CALL DepositionFast(nice,nspec+1,tk,adn,a_nicep,a_micep,     &
                             irnt,irmt,remice,irate,sfcirate,VtIce,4  )
         
         a_nicet%d(:,:,:,:) = a_nicet%d(:,:,:,:) + irnt(:,:,:,:)/dtlt
         a_micet%d(:,:,:,:) = a_micet%d(:,:,:,:) + irmt(:,:,:,:)/dtlt

         nc = spec%getIndex('H2O')
         ! Surface frozen precipitation as mm/h (liquid equivalent)
         !istr = getMassIndex(nice,1,nc); iend = getMassIndex(nice,nice,nc)
         !sfcirate%d(:,:) = SUM(remice(:,:,istr:iend),DIM=3)*3600.
         !istr = getMassIndex(nice,1,nc+1); iend = getMassIndex(nice,nice,nc+1)
         !sfcirate%d(:,:) = sfcirate%d(:,:) +  &
         !     SUM(remice(:,:,istr:iend),DIM=3)*3600.
         
         ! Convert mass flux to heat flux (W/m^2)
         !irate%d(:,:,:) = irate%d(:,:,:)*alvi

         istr = getMassIndex(nice,1,nc)
         iend = getMassIndex(nice,nice,nc+1) ! Include both pristine and rime
         DO j = 3, nyp-2
            DO i = 3, nxp-2
               DO k = 1, nzp-1
                  tlt%d(k,i,j) = tlt%d(k,i,j) -        &
                       (SUM(irmt(k,i,j,istr:iend))/dtlt)*(alvi/cp)*th%d(k,i,j)/tk%d(k,i,j)
               END DO
            END DO
         END DO
      END IF
      
      !IF (mcflg) THEN
      !   ! For mass conservation statistics
      !   mctmp(:,:) = 0.
      !   ss = spec%getIndex('H2O')
      !   istr = getMassIndex(nbins,1,ss); iend = getMassIndex(nbins,nbins,ss)
      !   mctmp(:,:) = mctmp(:,:) + SUM(remaer(:,:,istr:iend),dim=3)
      !   istr = getMassIndex(ncld,1,ss); iend = getMassIndex(ncld,ncld,ss)
      !   mctmp(:,:) = mctmp(:,:) + SUM(remcld(:,:,istr:iend),dim=3)
      !   istr = getMassIndex(nprc,1,ss); iend = getMassIndex(nprc,nprc,ss)
      !   mctmp(:,:) = mctmp(:,:) + SUM(remprc(:,:,istr:iend),dim=3)
      !   istr = getMassIndex(nice,1,ss); iend = getMassIndex(nice,nice,ss)
      !   mctmp(:,:) = mctmp(:,:) + SUM(remice(:,:,istr:iend),dim=3)
      !   CALL acc_massbudged(n1,n2,n3,3,tstep,dzt,a_dn,rdep=mctmp,ApVdom=mc_ApVdom)
      !END IF !mcflg
      ! Aerosol removal statistics
      !IF (sflg) CALL acc_removal(n2,n3,nspec,remaer,remcld,remprc,remice(:,:,1:nspec*nice)) ! Dont include rime here yet!!
      !IF (sflg) CALL cs_rem_set(n2,n3,nspec,remaer,remcld,remprc,remice(:,:,1:nspec*nice))

   END SUBROUTINE sedim_SALSA
 ! -----------------------------------------------------------------



   SUBROUTINE DepositionSlow(nb,ns,tk,adn,ustar,numc,mass,flxdivn,flxdivm,    &
                             depflxn,depflxm,rate,srate,flag)
    USE mo_particle_external_properties, ONLY : calcDiamLES, terminal_vel     
    USE util, ONLY : getBinMassArray
    USE mo_submctl, ONLY : nlim,prlim,pi6
    IMPLICIT NONE

    INTEGER, INTENT(in)            :: ns                ! number of chemical species (note that with ice the latter should contain also rime)
    INTEGER, INTENT(in)            :: nb                ! Number of bins
    TYPE(FloatArray3d), INTENT(in) :: tk                ! Absolute temprature
    TYPE(FloatArray3d), INTENT(in) :: adn               ! Air density
    TYPE(FloatArray2d), INTENT(in) :: ustar             !
    TYPE(FloatArray4d), INTENT(in) :: numc              ! Particle number concentration
    TYPE(FloatArray4d), INTENT(in) :: mass              ! Particle mass mixing ratio
    INTEGER, INTENT(IN)            :: flag         ! An option for identifying aerosol, cloud, precp and ice (1,2,3,4)
    REAL, INTENT(OUT)              :: flxdivm(nzp,nxp,nyp,nb*ns), flxdivn(nzp,nxp,nyp,nb) ! Mass and number divergency
    REAL, INTENT(OUT)              :: depflxn(nxp,nyp,nb), depflxm(nxp,nyp,nb*ns) ! Mass and number deposition fluxes to the surface
    TYPE(FloatArray3d), INTENT(inout) :: rate     ! Rain rate and surface rain rate diagnostics (cumulative, W/m^2) 
    TYPE(FloatArray2d), INTENT(inout) :: srate
    
    INTEGER :: i,j,k,kp1
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

    REAL :: rflm(nzp,nb*ns), rfln(nzp,nb), pmass(ns), dwet
    
    REAL :: zpm(nb*ns)  ! Bin mass array to clean things up
    REAL :: zpn(nb)     ! Bin number array to clean things up
    REAL :: zdn         ! Particle density
    
    REAL :: clim ! concentration limit

    clim = nlim
    IF (ANY(flag == [3,4])) clim = prlim
    
    flxdivm = 0.
    flxdivn = 0.
    depflxm = 0.
    depflxn = 0.
    zpm(:) = 0.
    zpn(:) = 0.
    !$omp parallel firstprivate(rflm,rfln,zpm,zpn,pmass)
    !$omp do private(k,bs,avis,kvis,lambda,dwet,mdiff,Sc,St,rt,Kn,zdn,GG,vc,kp1)
    DO j = 3,nyp-2
       DO i = 3,nxp-2

          rflm = 0.
          rfln = 0.
          
          DO k=nzp-1,2,-1
             kp1 = k+1

             ! atm modelling Eq.4.54
             avis=1.8325e-5*(416.16/(tk%d(k,i,j)+120.0))*(tk%d(k,i,j)/296.16)**1.5
             kvis =  avis/adn%d(k,i,j)
             va = sqrt(8*kb*tk%d(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2*avis/(adn%d(k,i,j)*va) !mean free path

             ! Fluxes
             !------------------
             ! -- Calculate the *corrections* for small particles
             zpm(:) = mass%d(k,i,j,:)
             zpn(:) = numc%d(k,i,j,:)

             pmass = 0.
             DO bin = 1,nb
                IF (zpn(bin)<clim) CYCLE

                ! Calculate wet size
                CALL getBinMassArray(nb,ns,bin,zpm,pmass)
                IF (flag < 4) THEN
                   dwet = calcDiamLES(ns,zpn(bin),pmass,flag,sph=.TRUE.)   
                ELSE
                   dwet = calcDiamLES(ns,zpn(bin),pmass,flag,sph=.FALSE.) ! For ice, this will be the max diameter for non-spherical particles
                END IF
                   
                ! Calculate the particle density based on dwet; note that for ice this will therefore be the "effective" density
                ! for a sphere with dwet, and may thus be very low for non-spherical particles
                zdn = SUM(pmass)/zpn(bin)/(pi6*dwet**3)

                ! Terminal velocity
                Kn = 2.*lambda/dwet      !lambda/rwet
                GG = 1.+ Kn*(A+B*exp(-C/Kn))
                vc = terminal_vel(dwet,zdn,adn%d(k,i,j),avis,GG,flag)

                ! This algorithm breaks down if the fall velocity is large enough to make the fall 
                ! distance greater than the grid level thickness (mainly an issue with very high 
                ! vertical resolutions). As a simple solution, limit the fall velocity based on grid 
                ! box thickness and timestep (should not cause any further issues).
                vc = MIN( vc, MIN(0.5*(1./dzt%d(k))/dtlt, 2.0) )

                ! POISTA
                IF (vc > 10. .OR. vc < 0.) WRITE(*,*) 'DEP SLOW ', vc
                
                IF (k==2) THEN ! The level just above surface
                    ! Particle diffusitivity  (15.29) in jacobson book
                    mdiff = (kb*tk%d(k,i,j)*GG)/(3.0*pi*avis*dwet) !(kb*tk(k,i,j)*GG)/(6.0*pi*avis*rwet)
                    Sc = kvis/mdiff
                    St = vc*ustar%d(i,j)**2.0/g*kvis
                    if (St<0.01) St=0.01
                    rt = 1.0/MAX(epsilon(1.0),(ustar%d(i,j)*(Sc**(-2.0/3.0)+10**(-3.0/St)))) ! atm chem&phy eq19.18
                    vc = (1./rt) + vc
                ENDIF

                ! Flux for the particle mass
                DO bs = bin, (ns-1)*nb + bin, nb
                   rflm(k,bs) = -zpm(bs)*vc   ! kg/kg m/s
                END DO

                ! Flux for the particle number
                 rfln(k,bin) = -zpn(bin)*vc
             END DO

             flxdivm(k,i,j,:) = (rflm(kp1,:)*adn%d(kp1,i,j)-rflm(k,:)*adn%d(k,i,j))*dzt%d(k)/adn%d(k,i,j)  ! kg/kg/s
             flxdivn(k,i,j,:) = (rfln(kp1,:)*adn%d(kp1,i,j)-rfln(k,:)*adn%d(k,i,j))*dzt%d(k)/adn%d(k,i,j)

                                               ! Sum of water flux from all bins
             rate%d(k,i,j) = rate%d(k,i,j) - SUM(rflm(k,(ns-1)*nb:ns*nb)) * adn%d(k,i,j) * alvl

          END DO ! k

          ! Deposition flux to surface
          k=2
          depflxm(i,j,:) = -rflm(k,:)*dzt%d(k) ! kg/kg/s
          depflxn(i,j,:) = -rfln(k,:)*dzt%d(k)

          srate%d(i,j) = srate%d(i,j) - SUM(rflm(k,(ns-1)*nb:ns*nb)) * adn%d(k,i,j) * alvl
          
       END DO ! i
    END DO ! j
    !$omp end do
    !$omp end parallel
    
  END SUBROUTINE DepositionSlow


  !------------------------------------------------------------------
  SUBROUTINE DepositionFast(nb,ns,tk,adn,numc,mass,prnt,prvt,remprc,rate,srate,Vt,flag)
    USE mo_particle_external_properties, ONLY : calcDiamLES, terminal_vel
    USE util, ONLY : getBinMassArray
    USE mo_submctl, ONLY : nlim,prlim,pi6
    USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ns,nb
    TYPE(FloatArray3d), INTENT(in) :: tk
    TYPE(FloatArray3d), INTENT(in) :: adn
    TYPE(FloatArray4d), INTENT(in) :: numc
    TYPE(FloatArray4d), INTENT(in) :: mass
    INTEGER, INTENT(IN) :: flag         ! An option for identifying liquid and ice (1,2,3,4)
    REAL, INTENT(out) :: prnt(nzp,nxp,nyp,nb), prvt(nzp,nxp,nyp,nb*ns)     ! Number and mass tendencies due to fallout
    REAL, INTENT(out) :: remprc(nxp,nyp,nb*ns)
    TYPE(FloatArray3d), INTENT(inout) :: rate ! Precip rate (W/m^2)
    TYPE(FloatArray2d), INTENT(inout) :: srate ! Surface precip rate (W/m^2)
    TYPE(FloatArray4d), INTENT(inout) :: Vt   ! Binned particle terminal velocity
    
    INTEGER :: k,i,j,bin
    INTEGER :: istr,iend

    real, parameter :: A = 1.249 ! fundamentals of atm. modelling pg509
    real, parameter :: B = 0.42
    real, parameter :: C = 0.87
    real, parameter :: M = 4.8096e-26 ! average mass of one air molecule, eq2.3 fundamentals of atm.
                                      ! modelling [kg molec-1]
    REAL :: GG
    REAL :: Kn
    REAL :: lambda      ! Mean free path
    REAL :: avis,kvis   ! Air viscosity, kinematic viscosity
    REAL :: va          ! Thermal speed of air molecule
    REAL :: vc

    ! For precipitation:
    REAL :: fd,fdmax,fdos ! Fall distance for rain drops, max fall distance, overshoot from nearest grid level
    REAL :: prnchg(nzp,nb), prvchg(nzp,nb,ns) ! Instantaneous changes in precipitation number and mass (volume)
    REAL :: dwet  ! Spherical wet diameter
    REAL :: dnsp  ! Diameter for non-spherical ice
    
    REAL :: prnumc, pmass(ns) ! Instantaneous source number and mass
    INTEGER :: kf, ni,fi
    LOGICAL :: prcdep  ! Deposition flag

    REAL :: clim  ! concentration limit
    
    REAL :: zpm(nb*ns), zpn(nb)  ! Bin mass and number arrays to clean things up
    REAL :: zdn                  ! Particle density
    REAL :: zdneff               ! Effective density for non-spherical

    TYPE(t_shape_coeffs) :: shape ! Used for ice
    
    clim = nlim
    IF (ANY(flag == [3,4])) clim = prlim

    ! Zero the output diagnostics for terminal velocity
    Vt%d = 0.
    
    remprc(:,:,:) = 0.
    prnt(:,:,:,:) = 0.
    prvt(:,:,:,:) = 0.
    zpm(:) = 0.
    zpn(:) = 0.
    !$omp parallel firstprivate(prnchg,prvchg,zpm,zpn,pmass,shape)
    !$omp do private(bin,ni,kvis,va,lambda,dwet,dnsp,vc,prcdep,kf,fdos,prnumc,fdmax,fd,zdn,Kn,GG,fi,avis)
    DO j = 3,nyp-2
       DO i = 3,nxp-2

          prnchg = 0.
          prvchg = 0.
          
          DO k=nzp-1,2,-1
          
             ! atm modelling Eq.4.54
             avis = 1.8325e-5*(416.16/(tk%d(k,i,j)+120.0))*(tk%d(k,i,j)/296.16)**1.5
             kvis = avis/adn%d(k,i,j) 
             va = sqrt(8.*kb*tk%d(k,i,j)/(pi*M)) ! thermal speed of air molecule
             lambda = 2.*avis/(adn%d(k,i,j)*va) !mean free path

             zpm(:) = mass%d(k,i,j,:)
             zpn(:) = numc%d(k,i,j,:)

             pmass = 0.
             ! Precipitation bin loop
             DO bin = 1,nb
                IF (zpn(bin) < clim) CYCLE

                ! Calculate wet size
                CALL getBinMassArray(nb,ns,bin,zpm,pmass)
                IF (flag < 4) THEN
                   dwet=calcDiamLES(ns,zpn(bin),pmass,flag,sph=.TRUE.)   
                ELSE
                   dwet=calcDiamLES(ns,zpn(bin),pmass,flag,sph=.TRUE.) 
                   dnsp=calcDiamLES(ns,zpn(bin),pmass,flag,sph=.FALSE.) ! For non-spherical ice, this is the max diameter of the crystal
                   CALL getShapeCoefficients(shape,SUM(pmass(1:ns-1)),pmass(ns),zpn(bin))
                END IF
                
                ! Calculate particle density based on dwet; for non-spherical ice this will get the "effective" density, which amy be
                ! quite low for large non-spherical crystals
                zdn = SUM(pmass)/zpn(bin)/(pi6*dwet**3)
                
                ! Terminal velocity
                Kn = 2.*lambda/dwet   !lambda/rwet
                GG = 1.+ Kn*(A+B*exp(-C/Kn))
                
                IF (flag < 4) THEN
                   vc = terminal_vel(dwet,zdn,adn%d(k,i,j),avis,GG,flag)
                ELSE
                   vc = terminal_vel(dwet,zdn,adn%d(k,i,j),avis,GG,flag,shape,dnsp)
                END IF
                ! Diagnostics
                Vt%d(k,i,j,bin) = vc
                
                ! Determine output flux for current level: Find the closest level to which the
                ! current drop parcel can fall within 1 timestep. If the lowest atmospheric level
                ! is reached, the drops are sedimented.
                
                ! Maximum fall distance:
                fdmax = dtlt*vc
                
                fd = 0.
                fi = 0
                prcdep = .FALSE. ! deposition flag
                DO WHILE ( fd < fdmax )
                   fd = fd + ( 1./dzt%d(k-fi) )
                   fi = fi + 1
                   ! Check if sedimentation occurs for current parcel
                   IF (k-fi <= 1) THEN
                      prcdep = .TRUE.
                      EXIT
                   END IF
                END DO
                ! Back up one level
                fi = fi-1
                kf = k - fi
                fd = fd - ( 1./dzt%d(kf) )
             
                ! How much the actual fall distance overshoots below the layer kf
                fdos = MIN(MAX(fdmax-fd,0.),1./dzt%d(kf))
             
                ! Remove the drops from the original level
                pmass = 0.
                prnumc = zpn(bin)
                prnchg(k,bin) = prnchg(k,bin) - prnumc
                DO ni = 1,ns
                   pmass(ni) = zpm((ni-1)*nb+bin)
                   prvchg(k,bin,ni) = prvchg(k,bin,ni) - pmass(ni)
                END DO ! ni

                ! Removal statistics
                IF (prcdep) THEN
                   DO ni=1,ns
                      remprc(i,j,(ni-1)*nb+bin) = remprc(i,j,(ni-1)*nb+bin) +    &
                           pmass(ni)*adn%d(k,i,j)*vc
                   END DO
                ENDIF ! prcdep

                ! precip rate statistics: removal of water from the current bin is accounted for
                ! Water is the last (nspec) species and precip rate is given here W/m^2
                ! Juha: For ice, have to sum up the pristine and rimed ice
                IF (flag == 4) THEN
                   rate%d(k,i,j)=rate%d(k,i,j)+(zpm((ns-2)*nb+bin)+zpm((ns-1)*nb+bin))*adn%d(k,i,j)*vc*alvi
                ELSE
                   rate%d(k,i,j)=rate%d(k,i,j)+zpm((ns-1)*nb+bin)*adn%d(k,i,j)*vc*alvl
                END IF
                                   
                ! Put the drops to new positions (may be partially the original grid cell as well!)
                IF (fdos*dzt%d(kf) > 0.5) THEN  ! Reduce numerical diffusion
                   prnchg(kf-1,bin) = prnchg(kf-1,bin) + prnumc
                   DO ni = 1,ns
                      prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + pmass(ni)
                   END DO
                ELSE
                   prnchg(kf-1,bin) = prnchg(kf-1,bin) + ( fdos*dzt%d(kf) )*prnumc
                   prnchg(kf,bin) = prnchg(kf,bin) + ( 1. - fdos*dzt%d(kf) )*prnumc
                   DO ni = 1,ns
                      prvchg(kf-1,bin,ni) = prvchg(kf-1,bin,ni) + ( fdos*dzt%d(kf) )*pmass(ni)
                      prvchg(kf,bin,ni) = prvchg(kf,bin,ni) + ( 1. - fdos*dzt%d(kf) )*pmass(ni)
                   END DO
                END IF ! diffusion
                
             END DO !bin
          END DO ! k

          ! Surface precip rate
          IF (flag == 4) THEN
             srate%d(i,j) = srate%d(i,j) + SUM(remprc(i,j,(ns-2)*nb:ns*nb)) * alvi
          ELSE
             srate%d(i,j) = srate%d(i,j) + SUM(remprc(i,j,(ns-1)*nb:ns*nb)) * alvl
          END IF
          
          prnt(:,i,j,:) = prnt(:,i,j,:) + prnchg(:,:)
          DO ni = 1,ns
             istr = (ni-1)*nb+1
             iend = ni*nb
             prvt(:,i,j,istr:iend) = prvt(:,i,j,istr:iend) + prvchg(:,:,ni)
          
          END DO !ni
       END DO ! i
    END DO ! j
    !$omp end do
    !$omp end parallel
    
  END SUBROUTINE DepositionFast



END MODULE mcrp
