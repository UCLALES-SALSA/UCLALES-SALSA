MODULE nudg
  USE grid, ONLY : level, dtlt, nxp, nyp, nzp, &
                   zt, a_rp, a_rt, a_rpp, a_rc, a_srp, a_ri, a_srs, &
                   a_naerop, a_naerot, a_ncloudp, a_nicep, &
                   a_tp, a_tt, a_up, a_ut, a_vp, a_vt
  USE mo_submctl, ONLY : nbins, ncld, nice, in2a, fn2b
  USE nudg_defs 
  
   IMPLICIT NONE

   CHARACTER(len=50), PARAMETER :: global_name = "nudg"


  CONTAINS
  !
  ! Main nudging procedures
  ! --------------------------
  !
  SUBROUTINE init_nudg
    IMPLICIT NONE

    CHARACTER(len=50), PARAMETER :: name = "init_nudg"
    REAL :: hlp1d(nzp), hlp2d(nzp,nbins)

    hlp1d = 0.
    hlp2d = 0.
   
    ! (Liquid water) potential temperature: nudge towards th0(:)-th00 ! Juha: ??
    IF (ndg_theta%nudgetype > 0) THEN
       !Ali, if it is not allocated during reading a history file
       IF (.NOT.ALLOCATED(theta_ref) )  ALLOCATE(theta_ref(nzp))
       theta_ref(:) = a_tp(:,3,3)
    END IF
    !
    ! Water vapor mixing ratio based on total water
    !   Levels 0-3: total = cloud water and water vapor (a_rp) + rain water (a_rpp)
    !   Levels 4-5: total = water vapor (a_rp) + condensate (a_rc) + rain water (a_srp)
    !                            + ice (a_ri) + snow (a_srs)
    IF (ndg_rv%nudgetype > 0)  THEN
      !Ali, if it is not allocated during reading a history file 
      IF (.NOT.ALLOCATED(rv_ref) )  ALLOCATE(rv_ref(nzp))
       IF (level == 5) THEN
          hlp1d = a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)+a_ri(:,3,3)+a_srs(:,3,3)
          rv_ref(:) = hlp1d(:)
       ELSE IF (level == 4) THEN
          hlp1d = a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)
          rv_ref(:) = hlp1d(:)
       ELSE ! Levels 0-3
          hlp1d = a_rp(:,3,3)+a_rpp(:,3,3)
          rv_ref(:) = hlp1d(:)
       END IF
    END IF
    !
    ! Horizontal winds
    IF (ndg_u%nudgetype > 0) THEN
      !Ali, if it is not allocated during reading a history file 
      IF (.NOT.ALLOCATED(u_ref) ) ALLOCATE(u_ref(nzp))
       u_ref(:) = a_up(:,3,3)
    END IF
    IF (ndg_v%nudgetype > 0) THEN
      !Ali, if it is not allocated during reading a history file 
      IF (.NOT.ALLOCATED(v_ref) ) ALLOCATE(v_ref(nzp))
       v_ref(:) = a_vp(:,3,3)
    END IF
    !
    ! Aerosol concentration for level 4. Nudge aerosol concentration based on
    ! total CCN = aerosol + cloud droplets + ice (a_nicep). Precipitation and snow
    ! are not included, because these cannot be related to a specific aerosol bin
    ! and their concentrations are low.
    IF (level > 3 .AND. ndg_aero%nudgetype > 0) THEN
      !Ali, if it is not allocated during reading a history file 
      IF (.NOT.ALLOCATED(aero_ref) ) ALLOCATE(aero_ref(nzp,nbins))
      ALLOCATE(aero_target(nzp,nbins))
       ! Nudge aerosol based on the total number (aerosol+cloud+ice)
       hlp2d = a_naerop(:,3,3,:)
       hlp2d(:,in2a:fn2b) = hlp2d(:,in2a:fn2b) + a_ncloudp(:,3,3,1:ncld)
       IF (level == 5) hlp2d(:,in2a:fn2b) = hlp2d(:,in2a:fn2b)+a_nicep(:,3,3,1:nice)
       
       aero_ref(:,:) = hlp2d(:,:)
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
      CHARACTER(len=50), PARAMETER :: name = "nudging"
      
      ! (Liquid water) potential temperature:
      IF ( ndg_theta%nudgetype > 0 ) &
           CALL nudge_any(nxp,nyp,nzp,zt,a_tp,dtlt,time,   &
           ndg_theta,theta_ref,a_tt)
      
      ! Water vapor
      IF ( ndg_rv%nudgetype > 0 ) THEN
         IF (level > 3) THEN
            ! Nudge water vapor (a_rp) based on total (vapor + cloud + rain)
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rc+a_srp,dtlt,time,   &
                 ndg_rv,rv_ref,a_rt)
         ELSE
            ! Nudge total water (a_rp) based on total + rain
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rpp,dtlt,time,   &
                 ndg_rv,rv_ref,a_rt)
         END IF
      END IF
      
      ! Horizontal winds
      IF ( ndg_u%nudgetype > 0 ) &
           CALL nudge_any(nxp,nyp,nzp,zt,a_up,dtlt,time,    &
           ndg_u,u_ref,a_ut)
      IF ( ndg_v%nudgetype > 0 ) &
           CALL nudge_any(nxp,nyp,nzp,zt,a_vp,dtlt,time,    &
           ndg_v,v_ref,a_vt)
      
      ! Aerosol 
      IF (level > 3 .AND. ndg_aero%nudgetype > 0 ) THEN
         ! Target aerosol concentration = total(t=0)-cloud(t)-ice(t)
         aero_target(:,:) = aero_ref(:,:)
         aero_target(:,in2a:fn2b) = aero_target(:,in2a:fn2b) - a_ncloudp(:,3,3,1:ncld)
         IF (level == 5) aero_target(:,in2a:fn2b) = aero_target(:,in2a:fn2b) - a_nicep(:,3,3,1:nice)
         ! Apply to sectional data
         CALL nudge_any_2d(nxp,nyp,nzp,nbins,zt,a_naerop,dtlt,time,   &
              ndg_aero,aero_target,a_naerot)
      END IF
      
   END SUBROUTINE nudging


   SUBROUTINE nudge_any(nx,ny,nz,zt,ap,dt,time,ndg_var,trgt,at)
      USE util, ONLY : get_avg3
      IMPLICIT NONE
      INTEGER, INTENT(in)               :: nx,ny,nz
      REAL, INTENT(in)                  :: zt(nz), ap(nz,nx,ny) 
      REAL, INTENT(in)                  :: dt, time
      TYPE(t_nudge), INTENT(in)         :: ndg_var
      REAL, INTENT(in)                  :: trgt(nzp)
      REAL, INTENT(inout)               :: at(nz,nx,ny)
      CHARACTER(len=50), PARAMETER      :: name = "nudge_any"
      INTEGER :: ii,jj
      REAL    :: avg(nz)
      REAL    :: tauloc
      REAL    :: diff1d(nz)
      
      LOGICAL :: nudgelev(nz)
      LOGICAL :: master_condition

      master_condition = ( time < nudge_time .OR.   &
                           ndg_var%tau_max_continue )

      IF ( master_condition ) THEN

         nudgelev(:) = ( nudge_zmin <= zt(:) .AND. zt(:) <= nudge_zmax )
         tauloc = ndg_var%f_tau(time)
         !
         IF (ndg_var%nudgetype == 1) THEN
            ! Soft nudging
            CALL get_avg3(nz,nx,ny,ap,avg)
            
            diff1d(:) = MERGE( (avg(:)-trgt(:))/max(tauloc,dt), 0., nudgelev )
            DO jj = 1,ny
               DO ii = 1,nx
                  at(:,ii,jj) = at(:,ii,jj) - diff1d(:)
               END DO
            END DO
            
         ELSE IF (ndg_var%nudgetype == 2) THEN
            ! Hard nudging
            DO jj = 1,ny
               DO ii = 1,nx
                  diff1d(:) = MERGE( (ap(:,ii,jj)-trgt(:))/max(tauloc,dt), 0., nudgelev )
                  at(:,ii,jj) = at(:,ii,jj) - diff1d(:)
               END DO
            END DO
            
         ELSE
            ! Unknown
            WRITE(*,*)'Unknown nudging option!'
            STOP
         END IF
         !
      END IF  ! master

   END SUBROUTINE nudge_any

   !
   ! Nudging for any 4D field based on 2D target
   SUBROUTINE nudge_any_2d(nx,ny,nz,nb,zt,ap,dt,time,ndg_var,trgt,at)
     USE util, ONLY : get_avg3
     IMPLICIT NONE
     INTEGER, INTENT(in)               :: nx,ny,nz,nb
     REAL, INTENT(in)                  :: zt(nz), ap(nz,nx,ny,nb)
     REAL, INTENT(in)                  :: dt, time
     TYPE(t_nudge), INTENT(in)         :: ndg_var
     REAL, INTENT(in)                  :: trgt(nzp,nbins)
     REAL, INTENT(inout)               :: at(nz,nx,ny,nb)
     CHARACTER(len=50), PARAMETER      :: name = "nudge_any_2d"
     INTEGER :: ii, jj,nn
     REAL    :: avg(nz)
     
     REAL :: diff1d(nz)
     REAL :: tauloc
     
     LOGICAL :: nudgelev(nz)
     LOGICAL :: master_condition
     
     master_condition = ( time < nudge_time .OR.   &
                          ndg_var%tau_max_continue )
     
     IF ( master_condition ) THEN
        
        nudgelev(:) = ( nudge_zmin <= zt(:) .AND. zt(:) <= nudge_zmax )
        tauloc = ndg_var%f_tau(time)
        !
        IF (ndg_var%nudgetype == 1) THEN
           ! Soft nudging
           DO nn = 1, nb
              CALL get_avg3(nz,nx,ny,ap(:,:,:,nn),avg)
              DO jj = 1,ny
                 DO ii = 1,nx
                    diff1d(:) = (avg(:)-trgt(:,nn))/max(tauloc,dt)
                    at(:,ii,jj,nn) = at(:,ii,jj,nn) - diff1d(:)
                 END DO
              END DO
           END DO
           
        ELSE IF (ndg_var%nudgetype == 2) THEN
           ! Hard nudging
           DO nn = 1, nb            
              DO jj = 1,ny
                 DO ii = 1,nx
                    diff1d(:) = (ap(:,ii,jj,nn)-trgt(:,nn))/max(tauloc,dt)
                    at(:,ii,jj,nn) = at(:,ii,jj,nn) - diff1d(:)
                 END DO
              END DO
           END DO
           
        ELSE
           ! Unknown
           WRITE(*,*)'Unknown nudging option!'
           STOP
        END IF
        
     END IF ! master
     
   END SUBROUTINE nudge_any_2d
   
END MODULE nudg
