MODULE nudg
  USE grid, ONLY : level, dtlt, nxp, nyp, nzp, &
                   zt, a_rp, a_rt, a_rpp, a_rc, a_srp, a_ri, a_srs, &
                   a_naerop, a_naerot, a_ncloudp, a_nicep, &
                   a_tp, a_tt, a_up, a_ut, a_vp, a_vt
  USE mo_submctl, ONLY : nbins, ncld, nice, in2a, fn2b
   IMPLICIT NONE

   ! Nudging options (nudge_*: 0=disabled, 1=soft, 2=hard), total nudging time (s), altitude range, and time constants (tau [s])
   LOGICAL :: useNudge = .FALSE.           ! Master switch for nudging
   INTEGER :: nudge_theta = 0, &           ! (liquid water) potential temperature, depending on the microphysical level
              nudge_rv = 0, &              ! Water vapor mixing ratio (maintain total water)
              nudge_u = 0, nudge_v = 0, &  ! Horizontal winds
              nudge_ccn = 0                ! Sectional aerosol for levels 4 and 5 (maintain aerosol+cloud+ice)
   REAL    :: nudge_time = 3600., nudge_zmin = -1.e10, nudge_zmax = 1.e10
   REAL    :: tau_theta = 300., tau_rv = 300., tau_u = 300., tau_v = 300., tau_ccn = 300.
   REAL, SAVE, ALLOCATABLE :: theta_ref(:), rv_ref(:), u_ref(:), v_ref(:), aero_ref(:,:)
   REAL, SAVE, ALLOCATABLE :: aero_target(:,:)

CONTAINS

  SUBROUTINE init_nudg
    IMPLICIT NONE

         ! (Liquid water) potential temperature: nudge towards th0(:)-th00 ! Juha: ??
         IF (nudge_theta /= 0) THEN
            ALLOCATE(theta_ref(nzp))
            theta_ref(:) = a_tp(:,3,3)
         END IF
         !
         ! Water vapor mixing ratio based on total water
         !   Levels 0-3: total = cloud water and water vapor (a_rp) + rain water (a_rpp)
         !   Levels 4-5: total = water vapor (a_rp) + condensate (a_rc) + rain water (a_srp)
         !                            + ice (a_ri) + snow (a_srs)
         IF (nudge_rv /= 0)  THEN
            ALLOCATE(rv_ref(nzp))
            IF (level == 5) THEN
               rv_ref(:) = a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)+a_ri(:,3,3)+a_srs(:,3,3)
            ELSE IF (level == 4) THEN
               rv_ref(:) = a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)
            ELSE ! Levels 0-3
               rv_ref(:) = a_rp(:,3,3)+a_rpp(:,3,3)
            END IF
         END IF
         !
         ! Horizontal winds
         IF (nudge_u /= 0) THEN
            ALLOCATE(u_ref(nzp))
            u_ref(:) = a_up(:,3,3)
         END IF
         IF (nudge_v /= 0) THEN
            ALLOCATE(v_ref(nzp))
            v_ref(:) = a_vp(:,3,3)
         END IF
         !
         ! Aerosol concentration for level 4. Nudge aerosol concentration based on
         ! total CCN = aerosol + cloud droplets + ice (a_nicep). Precipitation and snow
         ! are not included, because these cannot be related to a specific aerosol bin
         ! and their concentrations are low.
         IF (level > 3 .AND. nudge_ccn /= 0) THEN
            ! Nudge aerosol based on the total number (aerosol+cloud+ice)
            ALLOCATE(aero_ref(nzp,nbins),aero_target(nzp,nbins))
            aero_ref(:,:) = a_naerop(:,3,3,:)
            aero_ref(:,in2a:fn2b) = aero_ref(:,in2a:fn2b)+a_ncloudp(:,3,3,1:ncld)
            IF (level == 5) aero_ref(:,in2a:fn2b) = aero_ref(:,in2a:fn2b)+a_nicep(:,3,3,1:nice)
         END IF

  END SUBROUTINE init_nudg


   !
   !----------------------------------------------------------------------
   !
   ! Nudging towards the initial state (temperature, water vapor,
   ! horizontal winds and aerosol and/or cloud droplets).
   !
   ! TR 22.3.2017
   !

   SUBROUTINE nudging(time)

      IMPLICIT NONE
      REAL, INTENT(IN) :: time

      ! Nudging time is independent of the spin-up
      IF (time > nudge_time) RETURN

      ! (Liquid water) potential temperature:
      IF (nudge_theta > 0) &
         CALL nudge_any(nxp,nyp,nzp,zt,a_tp,a_tt,theta_ref,dtlt,tau_theta,nudge_theta)

      ! Water vapor
      IF (nudge_rv > 0) THEN
         IF (level > 3) THEN
            ! Nudge water vapor (a_rp) based on total (vapor + cloud + rain)
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rc+a_srp,a_rt,rv_ref,dtlt,tau_rv,nudge_rv)
         ELSE
            ! Nudge total water (a_rp) based on total + rain
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rpp,a_rt,rv_ref,dtlt,tau_rv,nudge_rv)
         END IF
      END IF

      ! Horizontal winds
      IF (nudge_u > 0) &
         CALL nudge_any(nxp,nyp,nzp,zt,a_up,a_ut,u_ref,dtlt,tau_u,nudge_u)
      IF (nudge_v > 0) &
         CALL nudge_any(nxp,nyp,nzp,zt,a_vp,a_vt,v_ref,dtlt,tau_v,nudge_v)

      ! Aerosol
      IF (level > 3 .AND. nudge_ccn /= 0) THEN
         ! Target aerosol concentration = total(t=0)-cloud(t)-ice(t)
         aero_target(:,:) = aero_ref(:,:)
         aero_target(:,in2a:fn2b) = aero_target(:,in2a:fn2b)-a_ncloudp(:,3,3,1:ncld)
         IF (level == 5) aero_target(:,in2a:fn2b) = aero_target(:,in2a:fn2b)-a_nicep(:,3,3,1:nice)
         ! Apply to sectional data
         CALL nudge_any_2d(nxp,nyp,nzp,nbins,zt,a_naerop,a_naerot,aero_target,dtlt,tau_ccn,nudge_ccn)
      END IF

   END SUBROUTINE nudging


   SUBROUTINE nudge_any(nxp,nyp,nzp,zt,ap,at,trgt,dt,tau,iopt)
      USE util, ONLY : get_avg3
      IMPLICIT NONE
      INTEGER :: nxp,nyp,nzp
      REAL    :: zt(nzp), ap(nzp,nxp,nyp), at(nzp,nxp,nyp)
      REAL    :: dt
      REAL    :: trgt(nzp)
      REAL    :: tau
      INTEGER :: iopt
      INTEGER :: kk
      REAL    :: avg(nzp)
      !
      IF (iopt == 1) THEN
         ! Soft nudging
         CALL get_avg3(nzp,nxp,nyp,ap,avg)
         DO kk = 1, nzp
            IF (nudge_zmin <= zt(kk) .AND. zt(kk) <= nudge_zmax) &
               at(kk,:,:) = at(kk,:,:)-(avg(kk)-trgt(kk))/max(tau,dt)
         END DO

      ELSE IF (iopt == 2) THEN
         ! Hard nudging
         DO kk = 1, nzp
            IF (nudge_zmin <= zt(kk) .AND. zt(kk) <= nudge_zmax) &
               at(kk,:,:) = at(kk,:,:)-(ap(kk,:,:)-trgt(kk))/max(tau,dt)
         END DO

      ELSE
         ! Unknown
         WRITE(*,*)'Unknown nudging option!'
         STOP
      END IF
     !
   END SUBROUTINE nudge_any

   !
   ! Nudging for any 4D field based on 2D target
   SUBROUTINE nudge_any_2d(nxp,nyp,nzp,nb,zt,ap,at,trgt,dt,tau,iopt)
      USE util, ONLY : get_avg3
      IMPLICIT NONE
      INTEGER :: nxp,nyp,nzp,nb
      REAL    :: zt(nzp), ap(nzp,nxp,nyp,nb), at(nzp,nxp,nyp,nb)
      REAL    :: dt
      REAL    :: trgt(nzp,nb)
      REAL    :: tau
      INTEGER :: iopt
      INTEGER :: ii, kk
      REAL    :: avg(nzp)
      !
      IF (iopt == 1) THEN
         ! Soft nudging
         DO ii = 1, nb
            CALL get_avg3(nzp,nxp,nyp,ap(:,:,:,ii),avg)
            DO kk = 1, nzp
               IF (nudge_zmin <= zt(kk) .AND. zt(kk) <= nudge_zmax) &
                  at(kk,:,:,ii) = at(kk,:,:,ii)-(avg(kk)-trgt(kk,ii))/max(tau,dt)
            END DO
         END DO

      ELSE IF (iopt == 2) THEN
         ! Hard nudging
         DO ii = 1, nb
            DO kk = 1, nzp
               IF (nudge_zmin <= zt(kk) .AND. zt(kk) <= nudge_zmax) &
                  at(kk,:,:,ii) = at(kk,:,:,ii)-(ap(kk,:,:,ii)-trgt(kk,ii))/max(tau,dt)
            END DO
         END DO

      ELSE
         ! Unknown
         WRITE(*,*)'Unknown nudging option!'
         STOP
      END IF
     !
   END SUBROUTINE nudge_any_2d


END MODULE nudg
