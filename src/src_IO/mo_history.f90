MODULE mo_history
  USE grid
  !Ali
  !These are must be written and read from history file
  !for consistent nudging initialization
  USE nudg_defs, ONLY : theta_ref, rv_ref, u_ref, v_ref, aero_ref, &
                        ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero
  USE mpi_interface, ONLY : appl_abort, myid, wrxid, wryid
  USE mo_vector_state
  USE mo_progn_state
  USE mo_diag_state
  USE mo_aux_state
  IMPLICIT NONE


  CONTAINS
    
    !
    ! ----------------------------------------------------------------------
    ! Subroutine write_hist:  This subroutine writes a binary history file
    !
    SUBROUTINE write_hist(htype, time)      
      INTEGER :: errcode = -17
      
      INTEGER, INTENT (in) :: htype
      REAL, INTENT (in)    :: time
      
      CHARACTER(len=20), PARAMETER :: name = "write_hist"
      
      CHARACTER (len=80) :: hname
      
      INTEGER :: n, iblank, ii, jj, kk,nn
      !
      ! create and open a new output file.
      !
      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(filprf)
      
      SELECT CASE(htype)
      CASE DEFAULT
         hname = trim(hname)//'.iflg'
      CASE(0)
         hname = trim(hname)//'.R'
      CASE(1)
         hname = trim(hname)//'.rst'
      CASE(2)
         iblank=index(hname,' ')
         WRITE(hname(iblank:iblank+7),'(a1,i6.6,a1)') '.', int(time), 's'
      END SELECT
      !
      ! Write fields
      !
      IF (myid == 0) PRINT "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
           ,hname
      OPEN(10,file=trim(hname), form='unformatted')
      
      WRITE(10) time,th00,umean,vmean,dtl,level,isgstyp,iradtyp,nzp,nxp,nyp,nscl
      WRITE(10) xt%d, xm%d, yt%d, ym%d, zt%d, zm%d, dn0%d, th0%d, u0%d, v0%d, pi0%d, &
                pi1%d, rt0%d, psrf,sst,W1,W2,W3 ! added by Zubair
      
      WRITE(10) a_ustar%d, a_tstar%d, a_rstar%d
      
      WRITE(10) a_pexnr%d
      WRITE(10) a_press%d
      WRITE(10) a_theta%d
      
      WRITE(10) a_up%d
      WRITE(10) a_vp%d
      WRITE(10) a_wp%d
      WRITE(10) a_uc%d
      WRITE(10) a_vc%d
      WRITE(10) a_wc%d
      
      DO n = 1, nscl
         CALL newsclr(n)  
         WRITE(10) a_sp
      END DO
      
      IF (ndg_theta%nudgetype > 0) THEN
         DO n = 1, nzp
            WRITE(10) theta_ref(n)
         END DO
      END IF
      
      IF (ndg_rv%nudgetype > 0) THEN
         DO n = 1, nzp
            WRITE(10) rv_ref(n)
         END DO
      END IF
      
      IF (ndg_u%nudgetype > 0) THEN
         DO n = 1, nzp
            WRITE(10) u_ref(n)
         END DO
      END IF
      
      IF (ndg_v%nudgetype > 0) THEN
         DO n = 1, nzp
            WRITE(10) v_ref(n)
         END DO
      END IF
      
      IF (level > 3 .AND. ndg_aero%nudgetype > 0) THEN
         WRITE(10) nbins
         DO n = 1, nzp
            DO nn = 1, nbins
               WRITE(10) aero_ref(n,nn)    
            END DO
         END DO
      END IF
      
      IF ( ASSOCIATED(a_rv%d)   ) WRITE(10) a_rv%d
      IF ( ASSOCIATED(a_rc%d)   ) WRITE(10) a_rc%d
      IF ( ASSOCIATED(a_rflx%d) ) WRITE(10) a_rflx%d
      CLOSE(10)
      
      IF (myid == 0 .AND. htype < 0) THEN
         PRINT *, 'CFL Violation'
         CALL appl_abort(errcode)
      END IF
      
      RETURN
    END SUBROUTINE write_hist
    !
    ! ----------------------------------------------------------------------
    ! Subroutine read_hist:  This subroutine reads a binary history file
    !
    !                        Modified for level 4
    !                Juha Tonttila, FMI, 20140828
    !
    
    SUBROUTINE read_hist(time, hfilin)      
      CHARACTER(len=80), INTENT(in) :: hfilin
      REAL, INTENT(out)             :: time
      
      CHARACTER(len=20), PARAMETER :: name = "read_hist"
      
      CHARACTER (len=80) :: hname
      INTEGER :: n, nxpx, nypx, nzpx, nsclx, iradx, isgsx, lvlx, ii, jj, kk
      LOGICAL :: exans
      REAL    :: umx, vmx, thx
      INTEGER :: nn, nnbins
      !
      ! open input file.
      !
      
      WRITE(hname,'(i4.4,a1,i4.4)') wrxid,'_',wryid
      hname = trim(hname)//'.'//trim(hfilin)
      
      inquire(file=trim(hname),exist=exans)
      IF (.NOT. exans) THEN
         PRINT *,'ABORTING: History file', trim(hname),' not found'
         CALL appl_abort(0)
      ELSE
         OPEN(10,file=trim(hname),status='old',form='unformatted')
         READ(10) time,thx,umx,vmx,dtl,lvlx,isgsx,iradx,nzpx,nxpx,nypx,nsclx
         
         IF (nxpx /= nxp .OR. nypx /= nyp .OR. nzpx /= nzp)  THEN
            IF (myid == 0) PRINT *, nxp, nyp, nzp, nxpx, nypx, nzpx
            CALL appl_abort(-1)
         END IF
         
         READ(10) xt%d, xm%d, yt%d, ym%d, zt%d, zm%d, dn0%d, th0%d, u0%d, v0%d, pi0%d, pi1%d, rt0%d, psrf,sst,W1,W2,W3
         
         READ(10) a_ustar%d, a_tstar%d, a_rstar%d
         
         READ(10) a_pexnr%d
         READ(10) a_press%d
         READ(10) a_theta%d
         
         READ(10) a_up%d
         READ(10) a_vp%d
         READ(10) a_wp%d
         READ(10) a_uc%d
         READ(10) a_vc%d
         READ(10) a_wc%d
         
         DO n = 1, nscl
            CALL newsclr(n)
            IF (n <= nsclx) READ(10) a_sp
         END DO

         IF (ndg_theta%nudgetype > 0) THEN
           ALLOCATE(theta_ref(nzp))
           DO n = 1, nzp
              READ(10) theta_ref(n)
           END DO
         END IF

         IF (ndg_rv%nudgetype > 0) THEN
           ALLOCATE(rv_ref(nzp))
           DO n = 1, nzp
              READ(10) rv_ref(n)
           END DO
         END IF

         IF (ndg_u%nudgetype > 0) THEN
           ALLOCATE(u_ref(nzp))
           DO n = 1, nzp
              READ(10) u_ref(n)
           END DO
         END IF

         IF (ndg_v%nudgetype > 0) THEN
           ALLOCATE(v_ref(nzp))
           DO n = 1, nzp
              READ(10) v_ref(n)
           END DO
         END IF

         IF (level > 3 .AND. ndg_aero%nudgetype > 0) THEN
           READ(10) nnbins
           ALLOCATE(aero_ref(nzp,nnbins))
           DO n = 1, nzp
             DO nn = 1, nbins
               READ(10) aero_ref(n,nn)    
             END DO
           END DO
         END IF
         
         DO n = nscl+1, nsclx
            READ(10)
         END DO

         IF (lvlx > 0 .AND. lvlx < 4) THEN
            IF (level > 0 .AND. lvlx < 4) THEN
               READ(10) a_rv%d
            ELSE
               READ(10)
            END IF
         END IF
         IF (lvlx > 1) THEN
            IF (level > 1) THEN
               READ(10) a_rc%d
            ELSE
               READ(10)
            END IF
         END IF
         IF (iradx > 0) THEN
            IF (iradtyp > 0) THEN
               READ(10) a_rflx%d
            ELSE
               READ(10)
            END IF
         END IF

         CLOSE(10)
         !
         ! adjust namelist and basic state appropriately
         !
         IF (thx /= th00) THEN
            IF (myid == 0) PRINT "('  th00 changed  -  ',2f8.2)",th00,thx
            a_tp%d(:,:,:) = a_tp%d(:,:,:) + thx - th00
         END IF
         IF (umx /= umean) THEN
            IF (myid == 0) PRINT "('  umean changed  -  ',2f8.2)",umean,umx
            a_up%d = a_up%d + umx - umean
         END IF
         IF (vmx /= vmean) THEN
            IF (myid == 0) PRINT "('  vmean changed  -  ',2f8.2)",vmean,vmx
            a_vp%d = a_vp%d + vmx - vmean
         END IF
         dtlv = 2.*dtl
         dtlt = dtl

      END IF

   END SUBROUTINE read_hist

END MODULE mo_history
