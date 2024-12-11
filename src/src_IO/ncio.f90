MODULE ncio

  USE netcdf
  USE mpi_interface, ONLY : appl_abort, myid, pecount
  USE classFieldArray, ONLY : FieldArray
  
  IMPLICIT NONE

!  INTERFACE write_nc
!     MODULE PROCEDURE write_nc_0d, write_nc_1d,   &
!                      write_nc_2d, write_nc_3d,   &
!                      write_nc_4d
!  END INTERFACE write_nc
  
  PRIVATE
  
  PUBLIC :: close_nc, sync_nc, &  
            open_aero_nc, read_aero_nc_1d, read_aero_nc_2d,  &
            StreamDef


  TYPE StreamDef
     ! Includes key parameters and IDs for output streams. Which dimension
     ! IDs and environment arrays are defined depends on the output stream
     ! configuration. The user should make sure that the dimensions required
     ! by the variables in the FieldArray instances conform with the dimensions
     ! defined for each output stream!

     ! Includes also type bound procedures for creating the output stream.
     ! The procedures for writing output are currently NOT bound to this type,
     ! but they could be. However, the only benefit would be to make the code
     ! a little bit more streamlined.
     
     INTEGER :: ncid = 0, nrec = 0, nvar = 0 ! NetCDF file handle, current number of records,
                                             ! number of variables staged for output in the current stream
     CHARACTER(len=150) :: fname             ! File name of the current output stream
     
     ! Dimension IDs (regular grid)
     INTEGER :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0
     ! Dimensions IDs (SALSA bin axes)
     INTEGER :: aeaID=0, claID=0, aebID=0, clbID=0, prcID=0, iceID=0  

     ! Dimension environment arrays used for variable definitions
     INTEGER :: dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0,  &
                dim_ztt(2)  = 0, dim_zmt(2)  = 0, dim_xtytt(3) = 0,                  &
                dim_ttttaea(5) = 0, dim_ttttcla(5) = 0,    &
                dim_ttttaeb(5) = 0, dim_ttttclb(5) = 0,    &
                dim_ttttprc(5) = 0, dim_ttttice(5) = 0,    &
                dim_ttaea(2) = 0, dim_ttcla(2) = 0,        &
                dim_ttaeb(2) = 0, dim_ttclb(2) = 0,        &
                dim_ttprc(2) = 0, dim_ttice(2) = 0,        &
                dim_zttaea(3) = 0, dim_zttcla(3) = 0,      &
                dim_zttaeb(3) = 0, dim_zttclb(3) = 0,      &
                dim_zttprc(3) = 0, dim_zttice(3) = 0,      &
                dim_taea(2), dim_taeb(2),                  &
                dim_tcla(2), dim_tclb(2),                  &
                dim_tprc(2), dim_tice(2)
     
     CONTAINS

       PROCEDURE :: open_nc
       ! These two are currently also used elsewhere; this discrepancy needs to be solved
       !PROCEDURE :: close_nc
       !PROCEDURE :: sync_nc
       PROCEDURE :: define_nc_dims
       PROCEDURE :: define_nc_vars
       PROCEDURE :: write_nc_0d,write_nc_1d,write_nc_2d,   &
                    write_nc_3d,write_nc_4d
       GENERIC :: write_nc => write_nc_0d, write_nc_1d, write_nc_2d, &
                              write_nc_3d, write_nc_4d
       
       PROCEDURE, PRIVATE :: defvar_loop
       PROCEDURE, PRIVATE :: define_variable
             
  END TYPE StreamDef

  
  CONTAINS

    ! Procedures bound to the TYPE StreamDef
    !
    ! ----------------------------------------------------------------------
    ! Subroutine Open_NC: Opens a NetCDF File and identifies starting record
    !
    SUBROUTINE open_nc (SELF, fname, ename, time, npts, version, author, info)
      CLASS(StreamDef), INTENT(inout) :: SELF
      INTEGER, INTENT(in)             :: npts
      REAL, INTENT (in)               :: time
      CHARACTER (len=150), INTENT (in) :: fname, ename
      CHARACTER(LEN=80), INTENT(in) :: version, author, info
      
      REAL, ALLOCATABLE :: xtimes(:)
      
      CHARACTER (len=8)  :: date
      INTEGER :: iret, ncall, VarID, RecordDimID
      LOGICAL :: exans
      
      inquire(file=trim(fname),exist=exans)
      SELF%fname=fname
      
      ncall = 0
      IF (.NOT. exans) THEN
         CALL date_and_time(date)
         iret = nf90_create(fname,NF90_SHARE,SELF%ncid)
         
         iret = nf90_put_att(SELF%ncid,NF90_GLOBAL,'title',ename)
         iret = nf90_put_att(SELF%ncid,NF90_GLOBAL,'history','Created on '//date)
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'Source','UCLALES-SALSA '//trim(version))
         IF (len(author) > 0) iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'Author',trim(author)) ! Optional
         IF (len(info) > 0) iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'Info',trim(info)) ! Optional
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, '_FillValue',-999.)
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'NPTS',npts)
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'NPROCS',pecount)
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'PROCID',myid)
         iret = nf90_put_att(SELF%ncid, NF90_GLOBAL, 'IO_version',1.1)
      ELSE
         iret  = nf90_open(trim(fname), NF90_WRITE, SELF%ncid)
         iret  = nf90_inquire(SELF%ncid, unlimitedDimId = RecordDimID)
         iret  = nf90_inquire_dimension(SELF%ncid, RecordDimID, len=SELF%nrec)
         ncall = 1
         iret  = nf90_inq_varid(SELF%ncid,'time',VarID)
         ALLOCATE (xtimes(SELF%nrec+1))
         iret  = nf90_get_var(SELF%ncid, VarId, xtimes(1:SELF%nrec))
         ncall = 1
         DO WHILE(ncall <= SELF%nrec .AND. xtimes(ncall) < time - spacing(1.))
            ncall = ncall+1
         END DO
         DEALLOCATE(xtimes)
      END IF
      SELF%nrec = ncall
      CALL sync_nc(SELF%ncid)
      
    END SUBROUTINE open_nc

    ! -------------------------------------
    ! Define dimensions
    !
    SUBROUTINE define_nc_dims(SELF,n1,n2,n3, &
                              inae_a,incld_a,  &
                              inprc,inae_b,    &
                              incld_b,inice    )

      CLASS(StreamDef), INTENT(inout) :: SELF
      INTEGER, OPTIONAL, INTENT (in) :: n1, n2, n3    
      INTEGER, OPTIONAL, INTENT(in)  :: inae_a,incld_a,inprc, &
                                        inae_b,incld_b,       &
                                        inice            
      INTEGER :: iret, VarID

      
      ! Every output stream shall have Time as the unlimited record axis 
      iret = nf90_def_dim(SELF%ncid, 'time', NF90_UNLIMITED, SELF%timeID)
      
      IF (present(n1)) THEN
         iret = nf90_def_dim(SELF%ncid, 'zt', n1, SELF%ztID)
         iret = nf90_def_dim(SELF%ncid, 'zm', n1, SELF%zmID)
      END IF
      IF (present(n2)) THEN
         iret = nf90_def_dim(SELF%ncid, 'xt', n2, SELF%xtID)
         iret = nf90_def_dim(SELF%ncid, 'xm', n2, SELF%xmID)
      END IF
      IF (present(n3)) THEN
         iret = nf90_def_dim(SELF%ncid, 'yt', n3, SELF%ytID)
         iret = nf90_def_dim(SELF%ncid, 'ym', n3, SELF%ymID)
      END IF
      IF (present(inae_a)) THEN
         iret = nf90_def_dim(SELF%ncid, 'aea', inae_a, SELF%aeaID)
      END IF
      IF (present(inae_b)) THEN
         iret = nf90_def_dim(SELF%ncid, 'aeb', inae_b, SELF%aebID)
      END IF
      IF (present(incld_a)) THEN
         iret = nf90_def_dim(SELF%ncid, 'cla', incld_a, SELF%claID)
      END IF
      IF (present(incld_b)) THEN
         iret = nf90_def_dim(SELF%ncid, 'clb', incld_b, SELF%clbID)
      END IF
      IF (present(inprc)) THEN
         iret = nf90_def_dim(SELF%ncid, 'prc', inprc, SELF%prcID)
      END IF
      IF (present(inice)) THEN
         iret = nf90_def_dim(SELF%ncid, 'ice', inice, SELF%iceID)
      END IF
      
      SELF%dim_xtytt = [SELF%xtID,SELF%ytID,SELF%timeID]
      SELF%dim_ztt = [SELF%ztID,SELF%timeID]
      SELF%dim_zmt = [SELF%zmID,SELF%timeID]
      SELF%dim_tttt= [SELF%ztID,SELF%xtID,SELF%ytID,SELF%timeID]  ! thermo point
      SELF%dim_mttt= [SELF%zmID,SELF%xtID,SELF%ytID,SELF%timeID]  ! zpoint
      SELF%dim_tmtt= [SELF%ztID,SELF%xmID,SELF%ytID,SELF%timeID]  ! upoint
      SELF%dim_ttmt= [SELF%ztID,SELF%xtID,SELF%ymID,SELF%timeID]  ! ypoint
      
      ! Juha: dimension environments for size distribution variables
      SELF%dim_ttttaea = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%aeaID,SELF%timeID]
      SELF%dim_ttttaeb = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%aebID,SELF%timeID]
      SELF%dim_ttttcla = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%claID,SELF%timeID]
      SELF%dim_ttttclb = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%clbID,SELF%timeID]
      SELF%dim_ttttprc = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%prcID,SELF%timeID]
      SELF%dim_ttttice = [SELF%ztID,SELF%xtID,SELF%ytID,SELF%iceID,SELF%timeID]  ! Jaakko
      ! ---
      ! Zubair: dimension environments for avegare size distribution variables per bin - ts files
      SELF%dim_taea = [SELF%aeaID,SELF%timeID]
      SELF%dim_taeb = [SELF%aebID,SELF%timeID]
      SELF%dim_tcla = [SELF%claID,SELF%timeID]
      SELF%dim_tclb = [SELF%clbID,SELF%timeID]
      SELF%dim_tprc = [SELF%prcID,SELF%timeID]
      SELF%dim_tice = [SELF%iceID,SELF%timeID] ! Jaakko
      ! Zubair: dimension environments for avegare size distribution variables per bin - ps files
      SELF%dim_zttaea = [SELF%ztID,SELF%aeaID,SELF%timeID]
      SELF%dim_zttaeb = [SELF%ztID,SELF%aebID,SELF%timeID]
      SELF%dim_zttcla = [SELF%ztID,SELF%claID,SELF%timeID]
      SELF%dim_zttclb = [SELF%ztID,SELF%clbId,SELF%timeID]
      SELF%dim_zttprc = [SELF%ztID,SELF%prcID,SELF%timeID]
      SELF%dim_zttice = [SELF%ztID,SELF%iceID,SELF%timeID] ! Jaakko      

      ! Assume that every file has at least the Time variable, so define that here
      iret = nf90_def_var(SELF%ncid,'time',NF90_FLOAT,SELF%timeID,VarID)
      iret = nf90_put_att(SELF%ncid,VarID,'longname','Time')
      iret = nf90_put_att(SELF%ncid,VarID,'units'   ,'s')      

    END SUBROUTINE define_nc_dims

    ! ---------------------------
    ! Define variables
    !
    SUBROUTINE define_nc_vars(SELF,outInst)
      CLASS(StreamDef), INTENT(inout) :: SELF
      TYPE(FieldArray), INTENT(in) :: outInst(:) ! This shall be an array of instances to facilitate
                                                 ! multiple variable sets
      INTEGER :: iret, i, N, Nexist

      N = SIZE(outInst)
      
      IF ( SELF%nrec == 0) THEN
         DO i = 1,N
            CALL SELF%defvar_loop(outInst(i))
         END DO

         iret = nf90_enddef(SELF%ncid)
         CALL sync_nc(SELF%ncid)
         SELF%nrec = 1
         
      ELSE  ! nrec /= 0

         SELF%nvar = 1 ! Time is always present
         DO i = 1,N
            SELF%nvar = SELF%nvar + outInst(i)%count
         END DO
         
         iret = nf90_inquire(SELF%ncid, nVariables=Nexist)
         !SELF%nvar = Nexist
         IF (Nexist /= SELF%nvar) THEN
            CALL close_nc(SELF%ncid)
            IF (myid == 0) PRINT *, '  ABORTING: Incompatible Netcdf File',n,SELF%nvar
            CALL appl_abort(0)
         ELSE
            CALL sync_nc(SELF%ncid)
         END IF
            
      END IF

    END SUBROUTINE define_nc_vars
    !
    ! ------------------------------------------------------------------------------------
    ! Type bound procedure defvar_loop: loops over the output variables and performs the variable
    !                         definition cycle
    SUBROUTINE defvar_loop(SELF,varInst)
      CLASS(StreamDef), INTENT(inout) :: SELF
      TYPE(FieldArray), INTENT(in) :: varInst
      
      INTEGER :: lnvar, n, VarID, iret
      CHARACTER(len=50) :: name
      CHARACTER(len=100) :: long_name
      CHARACTER(len=50) :: unit,dim

      ! Check that the variable array is initialized
      IF (varInst%Initialized) THEN      
         lnvar = varInst%count
         DO n = 1,lnvar
            ! Could also use ASSOCIATE here
            name = varInst%list(n)%name
            long_name = varInst%list(n)%long_name
            unit = varInst%list(n)%unit
            dim = varInst%list(n)%dimension
            
            ! Determine dimension and define variable
            CALL SELF%define_variable(dim,name,VarID)
            
            iret = nf90_put_att(SELF%ncid,VarID,'longname',long_name)
            iret = nf90_put_att(SELF%ncid,VarID,'units'   ,unit)
         END DO
         SELF%nvar = SELF%nvar + lnvar
      END IF
         
    END SUBROUTINE defvar_loop
      
    !
    ! ------------------------------------------------------------------------------------
    ! Subroutine define_variable: defines the output variables for netcdf
    !
    SUBROUTINE define_variable(SELF,dim,name,VarID)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=50), INTENT(in) :: dim
      CHARACTER(len=50), INTENT(in) :: name
      INTEGER, INTENT(out) :: VarID
      
      INTEGER :: iret

      ! Should add check and error messages if the dimensions required by the
      ! variable are defined for the current output stream.
      
      SELECT CASE(dim)
      ! Grid axis dimensions + for simple time series
      CASE ('time')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%timeID,VarID)
      CASE ('zt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%ztID,VarID)
      CASE ('zm')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%zmID,VarID)
      CASE ('xt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%xtID,VarID)
      CASE ('xm')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%xmID,VarID)
      CASE ('yt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%ytID,VarID)
      CASE ('ym')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%ymID,VarID)
      CASE ('aea')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%aeaID,VarID)
      CASE ('aeb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%aebID,VarID)
      CASE ('cla')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%claID,VarID)
      CASE ('clb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%clbID,VarID)
      CASE ('prc')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%prcID,VarID)
      CASE ('ice')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%iceID,VarID)
      ! //
      ! Binned 3d output   
      CASE ('ttttaea')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttaea,VarID)
      CASE ('ttttaeb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttaeb,VarID)
      CASE ('ttttcla')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttcla,VarID)
      CASE ('ttttclb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttclb,VarID)
      CASE ('ttttprc')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttprc,VarID)
      CASE ('ttttice')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttttice,VarID)
      ! //
      ! Regular 3d output   
      CASE ('tttt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tttt,VarID)
      CASE ('mttt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_mttt,VarID)
      CASE ('tmtt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tmtt,VarID)
      CASE ('ttmt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_ttmt,VarID)
      ! //
      ! Regular 2d output (x-y)
      CASE ('xtytt')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_xtytt,VarID)
      ! //
      ! Binned time series
      CASE ('taea')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_taea,VarID)
      CASE ('taeb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_taeb,VarID)
      CASE ('tcla')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tcla,VarID)
      CASE ('tclb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tclb,VarID)
      CASE ('tprc')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tprc,VarID)
      CASE ('tice')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_tice,VarID)
      !//
      ! Binned time-height profiles   
      CASE ('zttaea')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttaea,VarID)
      CASE ('zttaeb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttaeb,VarID)
      CASE ('zttcla')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttcla,VarID)
      CASE ('zttclb')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttclb,VarID)
      CASE ('zttprc')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttprc,VarID)
      CASE ('zttice')
         iret = nf90_def_var(SELF%ncid,name,NF90_FLOAT,SELF%dim_zttice,VarID)
      ! //
      ! Regular time-height profiles
      CASE ('ztt')
         iret = nf90_def_var(SELF%ncID,name,NF90_FLOAT,SELF%dim_ztt,VarID)
      CASE ('zmt')
         iret = nf90_def_var(SELF%ncID,name,NF90_FLOAT,SELF%dim_zmt,VarID)
         
      CASE DEFAULT
         IF (myid == 0) PRINT *, '  ABORTING: NCIO: Bad dimensional information ',trim(dim)
         CALL appl_abort(0)
      END SELECT
            
    END SUBROUTINE define_variable

    ! ---------------------------------
    ! Writing output (type bound)
    !
    SUBROUTINE write_nc_0d(SELF,name,var,beg)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var
      INTEGER, INTENT(in) :: beg(1)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(SELF%ncid,name,VarID)
      iret = nf90_put_var(SELF%ncid,VarID,var,start=beg)
    END SUBROUTINE write_nc_0d
    ! --------------------------------
    SUBROUTINE write_nc_1d(SELF,name,var,beg,icnt)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:)
      INTEGER, INTENT(in) :: beg(:)             ! Use automatic arrays to facilitate both axis variables and time-height profile variables
      INTEGER, OPTIONAL, INTENT(in) :: icnt(:)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(SELF%ncid,name,VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg,count=icnt)     
      ELSE
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg)
      END IF
    END SUBROUTINE write_nc_1d
    ! -------------------------------
    SUBROUTINE write_nc_2d(SELF,name,var,beg,icnt)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:)
      INTEGER, INTENT(in) :: beg(3)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(3)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(SELF%ncid,name,VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg,count=icnt)     
      ELSE
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg)
      END IF
    END SUBROUTINE write_nc_2d
    ! -------------------------------
    SUBROUTINE write_nc_3d(SELF,name,var,beg,icnt)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:,:)
      INTEGER, INTENT(in) :: beg(4)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(4)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(SELF%ncid,name,VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg,count=icnt)     
      ELSE
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg)
      END IF
    END SUBROUTINE write_nc_3d
    ! --------------------------------
    SUBROUTINE write_nc_4d(SELF,name,var,beg,icnt)
      CLASS(StreamDef), INTENT(in) :: SELF
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:,:,:)
      INTEGER, INTENT(in) :: beg(5)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(5)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(SELF%ncid,name,VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg,count=icnt)     
      ELSE
         iret = nf90_put_var(SELF%ncid,VarID,var,start=beg)
      END IF
    END SUBROUTINE write_nc_4d




    
    ! Other subroutines (unbound)
    ! ----------------------------------
    
    !
    !
    !
    !

 ! FUNCTIONS FOR READING ERA5 DATA FROM A NETCDF FILE
 !


 !

 SUBROUTINE open_ERA5_nc(ncid,nc_levs,nc_times)
    IMPLICIT NONE

    INTEGER, INTENT(out) :: ncid,nc_levs,nc_times
    INTEGER :: iret, did

    ! Open file
    iret = nf90_open('datafiles/UCLALES_era5.nc',NF90_NOWRITE,ncid)

    ! Inquire the number of input levels
    iret = nf90_inq_dimid(ncid,'z',did)
    iret = nf90_inquire_dimension(ncid,did,len=nc_levs)

    iret = nf90_inq_dimid(ncid,'time',did)
    iret = nf90_inquire_dimension(ncid,did,len=nc_times)

    !iret = nf90_inq_dimid(ncid,'nmod',did)
    !iret = nf90_inquire_dimension(ncid,did,len=nc_nmod)

 END SUBROUTINE open_ERA5_nc

 ! ----------------------------------------------------------------------
 ! FUNCTIONS FOR READING AEROSOL SIZE DISTRIBUTIONS FROM A NETCDF FILE
 !
 SUBROUTINE open_aero_nc(ncid,nc_levs,nc_nspec,nc_nmod)
    IMPLICIT NONE

    INTEGER, INTENT(out) :: ncid,nc_levs,nc_nspec,nc_nmod
    INTEGER :: iret, did

    ! Open file
    iret = nf90_open('datafiles/aerosol_in.nc',NF90_NOWRITE,ncid)

    ! Inquire the number of input levels
    iret = nf90_inq_dimid(ncid,'levs',did)
    iret = nf90_inquire_dimension(ncid,did,len=nc_levs)

    iret = nf90_inq_dimid(ncid,'nspec',did)
    iret = nf90_inquire_dimension(ncid,did,len=nc_nspec)

    iret = nf90_inq_dimid(ncid,'nmod',did)
    iret = nf90_inquire_dimension(ncid,did,len=nc_nmod)

 END SUBROUTINE open_aero_nc
 !
 ! ---------------------------------------------------
 !
 SUBROUTINE read_aero_nc_1d(ncid,name,d1,var)
    IMPLICIT NONE

    INTEGER, INTENT(in)           :: ncid, d1
    CHARACTER(len=*), INTENT(in) :: name
    REAL, INTENT(out)             :: var(d1)

    INTEGER :: iret,vid

    iret = nf90_inq_varid(ncid,name,vid)
    iret = nf90_get_var(ncid,vid,var)

 END SUBROUTINE read_aero_nc_1d
 !
 ! ---------------------------------------------------
 !
 SUBROUTINE read_aero_nc_2d(ncid,name,d1,d2,var)
   IMPLICIT NONE
   
   INTEGER, INTENT(in)           :: ncid, d1,d2
   CHARACTER(len=*), INTENT(in) :: name
   REAL, INTENT(out)             :: var(d1,d2)
   
   INTEGER :: iret, vid
   
   iret = nf90_inq_varid(ncid,name,vid)
   iret = nf90_get_var(ncid,vid,var)
   
 END SUBROUTINE read_aero_nc_2d
 !
 ! -----------------------------------------------------
 !

 ! -------------------------------------------
 ! Closig and syncing with files on disk
 !
 SUBROUTINE close_nc(fid)
   IMPLICIT NONE
   
   INTEGER, INTENT(in) :: fid
   
   INTEGER :: iret
   
   iret = nf90_close(fid)
   
 END SUBROUTINE close_nc
 
 SUBROUTINE sync_nc(fid)
   IMPLICIT NONE

   INTEGER, INTENT(in) :: fid
   INTEGER :: iret

   iret = nf90_sync(fid)
   
 END SUBROUTINE sync_nc


 

 END MODULE ncio
