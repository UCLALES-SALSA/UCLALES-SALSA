MODULE mo_history
  USE grid
  !Ali
  !These are must be written and read from history file
  !for consistent nudging initialization
  USE nudg_defs, ONLY : theta_ref, rv_ref, u_ref, v_ref, aero_ref, &
                        ndg_theta, ndg_rv, ndg_u, ndg_v, ndg_aero
  USE mpi_interface, ONLY : appl_abort,myid
  USE mo_mpi_io
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
      INTEGER :: ichar

      TYPE(mpi_file_parameters) :: fhist
      
      INTEGER :: n, iblank,nn
      INTEGER :: globalInts(4)  ! just to reduce the number of separate mpi write calls
      REAL    :: globalFloats(5)
      INTEGER :: localInts(3)
      REAL    :: localFloats(5)
      INTEGER :: nudgetypes(5)  ! Aerosol missing - needs some special stuff
      !
      ! create and open a new output file.
      !

      hname = trim(filprf)
      ichar = LEN(TRIM(hname))
      WRITE(hname(ichar+1:ichar+12),'(a1,i6.6,a1,a1,a3)')  &
           '.',int(time),'s','.','rst'
      !
      ! Write fields
      !
      IF (myid == 0) PRINT "(//' ',49('-')/,' ',/,'   History write to: ',A30)" &
           ,hname

      globalInts = [level,isgstyp,iradtyp,nscl]
      globalFloats = [time,th00,umean,vmean,dtl]
      localInts = [nzp,nxp,nyp]
      localFloats = [psrf,sst,W1,W2,W3]
      nudgetypes = [ndg_theta%nudgetype,ndg_rv%nudgetype,ndg_u%nudgetype,   &
                    ndg_v%nudgetype,ndg_aero%nudgetype]
      
      CALL create_mpi_hist(trim(hname),fhist)

      ! These values are identical for all processes -> write only from root
      CALL write_hist_mpi(4,globalInts,.TRUE.,fhist)
      CALL write_hist_mpi(5,globalFloats,.TRUE.,fhist)

      CALL write_hist_mpi(5,nudgetypes,.TRUE.,fhist)

      ! Process specific parameters
      CALL write_hist_mpi(3,localInts,.FALSE.,fhist)
      CALL write_hist_mpi(5,localFloats,.FALSE.,fhist)
      
      ! Basic state arrays - identical for all processses -> write only from root
      CALL write_hist_mpi(nzp,dn0%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,th0%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,u0%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,v0%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,pi0%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,pi1%d,.TRUE.,fhist)
      CALL write_hist_mpi(nzp,rt0%d,.TRUE.,fhist)
      
      ! Grid displacements
      CALL write_hist_mpi(nxp,xt%d,.FALSE.,fhist)
      CALL write_hist_mpi(nxp,xm%d,.FALSE.,fhist)
      CALL write_hist_mpi(nyp,yt%d,.FALSE.,fhist)
      CALL write_hist_mpi(nyp,ym%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,zm%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,zt%d,.FALSE.,fhist)

      ! 2d fields
      CALL write_hist_mpi(nxp,nyp,a_ustar%d,.FALSE.,fhist)
      CALL write_hist_mpi(nxp,nyp,a_tstar%d,.FALSE.,fhist)
      CALL write_hist_mpi(nxp,nyp,a_rstar%d,.FALSE.,fhist)

      ! 3d fields
      CALL write_hist_mpi(nzp,nxp,nyp,a_pexnr%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_press%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_theta%d,.FALSE.,fhist)

      CALL write_hist_mpi(nzp,nxp,nyp,a_up%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_vp%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_wp%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_uc%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_vc%d,.FALSE.,fhist)
      CALL write_hist_mpi(nzp,nxp,nyp,a_wc%d,.FALSE.,fhist)

      ! Prognostic scalars
      DO n = 1, nscl
         CALL newsclr(n)  
         CALL write_hist_mpi(nzp,nxp,nyp,a_sp,.FALSE.,fhist)
      END DO

      IF (ndg_theta%nudgetype > 0)  &          
           CALL write_hist_mpi(nzp,theta_ref,.FALSE.,fhist)
      
      IF (ndg_rv%nudgetype > 0)  &
           CALL write_hist_mpi(nzp,rv_ref,.FALSE.,fhist)
      
      IF (ndg_u%nudgetype > 0)  &
           CALL write_hist_mpi(nzp,u_ref,.FALSE.,fhist)
      
      IF (ndg_v%nudgetype > 0)  &
           CALL write_hist_mpi(nzp,v_ref,.FALSE.,fhist)
      
      ! AEROSOL NUDGE REF MISSING FOR NOW !!!!
      
      IF ( ASSOCIATED(a_rv%d)   ) CALL write_hist_mpi(nzp,nxp,nyp,a_rv%d,.FALSE.,fhist)
      IF ( ASSOCIATED(a_rc%d)   ) CALL write_hist_mpi(nzp,nxp,nyp,a_rc%d,.FALSE.,fhist)
      IF ( ASSOCIATED(a_rflx%d) ) CALL write_hist_mpi(nzp,nxp,nyp,a_rflx%d,.FALSE.,fhist)

      CALL close_mpi_hist(fhist)
      
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
      INTEGER :: n, nxpx, nypx, nzpx, nsclx, iradx, isgsx, lvlx 
      LOGICAL :: exans
      REAL    :: umx, vmx, thx
      INTEGER :: nn, nnbins

      TYPE(mpi_file_parameters) :: fhist
      
      INTEGER :: globalInts(4)  ! just to reduce the number of separate mpi write calls
      REAL    :: globalFloats(5)
      INTEGER :: localInts(3)
      REAL    :: localFloats(5)
      INTEGER :: nudgetypes(5)

      globalInts = 0
      globalFloats = 0.
      localInts = 0
      localFloats = 0.
      nudgetypes = 0
      
      !
      ! open input file.
      !
      hname = trim(hfilin)//'.rst'
      
      inquire(file=trim(hname),exist=exans)
      IF (.NOT. exans) THEN
         PRINT *,'ABORTING: History file', trim(hname),' not found'
         CALL appl_abort(0)
      ELSE

         ! For parallel reading functions, the argument onlyroot=TRUE means that
         ! the same memory address is read by all processes, so no need to broadcast
         ! those variables afterwards. All reading should take place in exactly the same
         ! order as they were written.
         
         CALL open_mpi_hist(trim(hname),fhist)

         ! Some parameters.
         CALL read_hist_mpi(4,globalInts,.TRUE.,fhist)
         CALL read_hist_mpi(5,globalFloats,.TRUE.,fhist)
         CALL read_hist_mpi(5,nudgetypes,.TRUE.,fhist)
         CALL read_hist_mpi(3,localInts,.FALSE.,fhist)
         CALL read_hist_mpi(5,localFloats,.FALSE.,fhist)
         
         ! Decompose the input arrays and do some checking of the configuration
         lvlx = globalInts(1); isgsx = globalInts(2); iradx = globalInts(3); nsclx = globalInts(4)
         time = globalFloats(1); thx = globalFloats(2); umx = globalFloats(3)
         vmx = globalFloats(4); dtl = globalFloats(5)
         ndg_theta%nudgetype = nudgetypes(1); ndg_rv%nudgetype = nudgetypes(2)
         ndg_u%nudgetype = nudgetypes(3); ndg_v%nudgetype = nudgetypes(4)
         ndg_aero%nudgetype = nudgetypes(5)
         nzpx = localInts(1); nxpx = localInts(2); nypx = localInts(3)
         psrf = localFloats(1); sst = localFloats(2); W1 = localFloats(3)
         W2 = localFloats(4); W3 = localFloats(5)
         
         IF (lvlx /= level) THEN
            IF (myid == 0) WRITE(*,*) 'RESTART ERROR: inconsistent LEVEL: ',lvlx, level
            CALL appl_abort(-1)
         END IF
         IF (nsclx /= nscl) THEN
            IF (myid == 0) WRITE(*,*)    &
                 'RESTART ERROR: inconsistent configuration; check prognostic scalar variables. ',  &
                 nsclx, nscl
            CALL appl_abort(-1)
         END IF
         IF (nxpx /= nxp .OR. nypx /= nyp .OR. nzpx /= nzp)  THEN
            IF (myid == 0) WRITE(*,*) 'RESTART ERROR: inconsistent grid: ', nxp, nyp, nzp, nxpx, nypx, nzpx
            CALL appl_abort(-1)
         END IF
         
         CALL read_hist_mpi(nzp,dn0%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,th0%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,u0%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,v0%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,pi0%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,pi1%d,.TRUE.,fhist)
         CALL read_hist_mpi(nzp,rt0%d,.TRUE.,fhist)         

         CALL read_hist_mpi(nxp,xt%d,.FALSE.,fhist)
         CALL read_hist_mpi(nxp,xm%d,.FALSE.,fhist)
         CALL read_hist_mpi(nyp,yt%d,.FALSE.,fhist)
         CALL read_hist_mpi(nyp,ym%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,zm%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,zt%d,.FALSE.,fhist)

         CALL read_hist_mpi(nxp,nyp,a_ustar%d,.FALSE.,fhist)
         CALL read_hist_mpi(nxp,nyp,a_tstar%d,.FALSE.,fhist)
         CALL read_hist_mpi(nxp,nyp,a_rstar%d,.FALSE.,fhist)

         ! 3d fields
         CALL read_hist_mpi(nzp,nxp,nyp,a_pexnr%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_press%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_theta%d,.FALSE.,fhist)
         
         CALL read_hist_mpi(nzp,nxp,nyp,a_up%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_vp%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_wp%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_uc%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_vc%d,.FALSE.,fhist)
         CALL read_hist_mpi(nzp,nxp,nyp,a_wc%d,.FALSE.,fhist)
         
         ! Prognostic scalars
         DO n = 1, nscl
            CALL newsclr(n)  
            CALL read_hist_mpi(nzp,nxp,nyp,a_sp,.FALSE.,fhist)
         END DO
                
         IF (ndg_theta%nudgetype > 0) THEN 
            ALLOCATE(theta_ref(nzp))
            CALL read_hist_mpi(nzp,theta_ref,.FALSE.,fhist)
         END IF

         IF (ndg_rv%nudgetype > 0) THEN
           ALLOCATE(rv_ref(nzp))
           CALL read_hist_mpi(nzp,rv_ref,.FALSE.,fhist)
         END IF

         IF (ndg_u%nudgetype > 0) THEN
           ALLOCATE(u_ref(nzp))
           CALL read_hist_mpi(nzp,u_ref,.FALSE.,fhist)
         END IF

         IF (ndg_v%nudgetype > 0) THEN
           ALLOCATE(v_ref(nzp))
           CALL read_hist_mpi(nzp,v_ref,.FALSE.,fhist)
         END IF

         ! AEROREF STILL MISSING

         IF (ASSOCIATED(a_rv%d)) &
              CALL read_hist_mpi(nzp,nxp,nyp,a_rv%d,.FALSE.,fhist)
              
         IF (ASSOCIATED(a_rc%d)) &
              CALL read_hist_mpi(nzp,nxp,nyp,a_rc%d,.FALSE.,fhist)

         IF (ASSOCIATED(a_rflx%d) .AND. iradx > 0) &
              CALL read_hist_mpi(nzp,nxp,nyp,a_rflx%d,.FALSE.,fhist)
         
         CALL close_mpi_hist(fhist)
         
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
