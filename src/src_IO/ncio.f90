MODULE ncio

  USE netcdf
  USE mpi_interface, ONLY : appl_abort, myid, pecount, wrxid, wryid
  USE classFieldArray, ONLY : FieldArray
  
  IMPLICIT NONE

  INTERFACE write_nc
     MODULE PROCEDURE write_nc_0d, write_nc_1d,   &
                      write_nc_2d, write_nc_3d,   &
                      write_nc_4d
  END INTERFACE write_nc
  
  PRIVATE
  
  PUBLIC :: open_nc, define_nc, define_nc_cs, &
       open_aero_nc, read_aero_nc_1d, read_aero_nc_2d,  &
       close_nc, sync_nc, write_nc 

  ! Save dimension IDs as global module variables
  INTEGER, SAVE :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,            &
             dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0,  &
             dim_ztt(2)  = 0, dim_zmt(2)  = 0, dim_xtytt(3)
      
  INTEGER, SAVE :: aeaID=0, claID=0, aebID=0, clbID=0, prcID=0,   &
             iceID=0,                     &
             dim_ttttaea(5) = 0, dim_ttttcla(5) = 0,    &
             dim_ttttaeb(5) = 0, dim_ttttclb(5) = 0,    &
             dim_ttttprc(5) = 0,                        &
             dim_ttttice(5) = 0,  &
             dim_ttaea(2) = 0, dim_ttcla(2) = 0,        &
             dim_ttaeb(2) = 0, dim_ttclb(2) = 0,        &
             dim_ttprc(2) = 0,                          &
             dim_ttice(2) = 0, &
             dim_zttaea(3) = 0, dim_zttcla(3) = 0,    &
             dim_zttaeb(3) = 0, dim_zttclb(3) = 0,    &
             dim_zttprc(3) = 0,                        &
             dim_zttice(3) = 0,                        &
             dim_taea(2), dim_taeb(2), dim_tcla(2), dim_tclb(2),  &
             dim_tprc(2), dim_tice(2)

  
  CONTAINS
    !
    ! ----------------------------------------------------------------------
    ! Subroutine Open_NC: Opens a NetCDF File and identifies starting record
    !
    SUBROUTINE open_nc (fname, ename, time, npts, ncid, nrec, version, author, info)
      
      INTEGER, INTENT(in)             :: npts
      INTEGER, INTENT(out)            :: ncid
      INTEGER, INTENT(out)            :: nrec
      REAL, INTENT (in)               :: time
      CHARACTER (len=150), INTENT (in) :: fname, ename
      CHARACTER(LEN=80) :: version, author, info
      
      REAL, ALLOCATABLE :: xtimes(:)
      
      CHARACTER (len=8)  :: date
      CHARACTER (len=200) :: lfname
      INTEGER :: iret, ncall, VarID, RecordDimID
      LOGICAL :: exans

      IF (pecount > 1) THEN
         WRITE(lfname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',wrxid,wryid,'.nc'
      ELSE
         WRITE(lfname,'(a,a3)') trim(fname),'.nc'
      END IF
      
      inquire(file=trim(lfname),exist=exans)
      
      ncall = 0
      IF (.NOT. exans) THEN
         CALL date_and_time(date)
         iret = nf90_create(lfname,NF90_SHARE,ncid)
         
         iret = nf90_put_att(ncid,NF90_GLOBAL,'title',ename)
         iret = nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'Source','UCLALES-SALSA '//trim(version))
         IF (len(author) > 0) iret = nf90_put_att(ncid, NF90_GLOBAL, 'Author',trim(author)) ! Optional
         IF (len(info) > 0) iret = nf90_put_att(ncid, NF90_GLOBAL, 'Info',trim(info)) ! Optional
         iret = nf90_put_att(ncid, NF90_GLOBAL, '_FillValue',-999.)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPTS',npts)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPROCS',pecount)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'PROCID',myid)
         iret = nf90_put_att(ncid, NF90_GLOBAL, 'IO_version',1.1)
      ELSE
         iret  = nf90_open(trim(lfname), NF90_WRITE, ncid)
         iret  = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
         iret  = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
         ncall = 1
         iret  = nf90_inq_varid(ncid,'time',VarID)
         ALLOCATE (xtimes(nrec+1))
         iret  = nf90_get_var(ncid, VarId, xtimes(1:nrec))
         ncall = 1
         DO WHILE(ncall <= nrec .AND. xtimes(ncall) < time - spacing(1.))
            ncall = ncall+1
         END DO
         DEALLOCATE(xtimes)
      END IF
      nrec = ncall
      iret = nf90_sync(ncid)
      
    END SUBROUTINE open_nc
    !
    ! ----------------------------------------------------------------------
    ! Subroutine Define_NC: Defines the structure of the nc file (if not
    ! already open)
    !
    ! Juha: Added more dimensions to represent bins for aerosol, cloud and
    ! precipitation particles.
    !
    SUBROUTINE define_nc(ncID,nRec,out_nvar,     &
                         outProg,outVector,      &
                         outDiag,outDerived,     &
                         outPS,outAxes,          &
                         n1,n2,n3,               &
                         inae_a,incld_a,         &
                         inprc,inae_b,           &
                         incld_b,inice           )
      
      INTEGER, INTENT (in)           :: ncID
      INTEGER, INTENT(inout) :: nRec
      INTEGER, INTENT(out) :: out_nvar
      TYPE(FieldArray), OPTIONAL, INTENT(in) :: outProg, outDiag, outDerived, &  ! these should be the subsets of FieldArray instances,
                                                outVector, outPS, outAxes        ! whose variables are marked as outputs, so the output
                                                                                 ! status does not have to be tested. New ones can be added
                                                                                 ! e.g. for statistical outputs
      INTEGER, OPTIONAL, INTENT (in) :: n1, n2, n3    
      INTEGER, OPTIONAL, INTENT(in)  :: inae_a,incld_a,inprc, &
                                        inae_b,incld_b,       &
                                        inice            

      CHARACTER (len=50) :: xnm
      INTEGER :: iret, n, VarID

      out_nvar = 0
      
      
      IF (nRec == 0) THEN
         iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
         IF (present(n1)) THEN
            iret = nf90_def_dim(ncID, 'zt', n1, ztID)
            iret = nf90_def_dim(ncID, 'zm', n1, zmID)
         END IF
         IF (present(n2)) THEN
            iret = nf90_def_dim(ncID, 'xt', n2, xtID)
            iret = nf90_def_dim(ncID, 'xm', n2, xmID)
         END IF
         IF (present(n3)) THEN
            iret = nf90_def_dim(ncID, 'yt', n3, ytID)
            iret = nf90_def_dim(ncID, 'ym', n3, ymID)
         END IF
         
         IF (present(inae_a)) THEN
            iret = nf90_def_dim(ncID, 'aea', inae_a, aeaID)
         END IF
         IF (present(inae_b)) THEN
            iret = nf90_def_dim(ncID, 'aeb', inae_b, aebID)
         END IF
         IF (present(incld_a)) THEN
            iret = nf90_def_dim(ncID, 'cla', incld_a, claID)
         END IF
         IF (present(incld_b)) THEN
            iret = nf90_def_dim(ncID, 'clb', incld_b, clbID)
         END IF
         IF (present(inprc)) THEN
            iret = nf90_def_dim(ncID, 'prc', inprc, prcID)
         END IF
         IF (present(inice)) THEN
            iret = nf90_def_dim(ncID, 'ice', inice, iceID)
         END IF

         dim_xtytt = [xtID,ytID,timeID]
         dim_ztt = [ztID,timeID]
         dim_zmt = [zmID,timeID]
         dim_tttt= [ztID,xtID,ytID,timeID]  ! thermo point
         dim_mttt= [zmID,xtID,ytID,timeID]  ! zpoint
         dim_tmtt= [ztID,xmID,ytID,timeID]  ! upoint
         dim_ttmt= [ztID,xtID,ymID,timeID]  ! ypoint
         
         ! Juha: dimension environments for size distribution variables
         dim_ttttaea = [ztID,xtID,ytID,aeaID,timeID]
         dim_ttttaeb = [ztID,xtID,ytID,aebID,timeID]
         dim_ttttcla = [ztID,xtID,ytID,claID,timeID]
         dim_ttttclb = [ztID,xtID,ytID,clbID,timeID]
         dim_ttttprc = [ztID,xtID,ytID,prcID,timeID]
         dim_ttttice = [ztID,xtID,ytID,iceID,timeID]  ! Jaakko
         ! ---
         ! Zubair: dimension environments for avegare size distribution variables per bin - ts files
         dim_taea = [aeaID,timeID]
         dim_taeb = [aebID,timeID]
         dim_tcla = [claID,timeID]
         dim_tclb = [clbID,timeID]
         dim_tprc = [prcID,timeID]
         dim_tice = [iceID,timeID] ! Jaakko
         ! Zubair: dimension environments for avegare size distribution variables per bin - ps files
         dim_zttaea = [ztID,aeaID,timeID]
         dim_zttaeb = [ztID,aebID,timeID]
         dim_zttcla = [ztID,claID,timeID]
         dim_zttclb = [ztID,clbId,timeID]
         dim_zttprc = [ztID,prcID,timeID]
         dim_zttice = [ztID,iceID,timeID] ! Jaakko
         
         ! Define variables
         ! First time
         iret = nf90_def_var(ncID,'time',NF90_FLOAT,timeID,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','Time')
         iret = nf90_put_att(ncID,VarID,'units'   ,'s')
         
         IF (PRESENT(outProg))  &
              CALL defvar_loop(ncID,outProg,out_nvar)
         IF (PRESENT(outVector)) &
              CALL defvar_loop(ncID,outVector,out_nvar)
         IF (PRESENT(outDiag))  &
              CALL defvar_loop(ncID,outDiag,out_nvar)
         IF (PRESENT(outAxes))  &
              CALL defvar_loop(ncID,outAxes,out_nvar)
         IF (PRESENT(outDerived)) &
              CALL defvar_loop(ncID,outDerived,out_nvar)
         IF (PRESENT(outPS)) &
              CALL defvar_loop(ncID,outPS,out_nvar)
         
         iret = nf90_enddef(ncID)
         iret = nf90_sync(ncID)
         nRec = 1
         
      ELSE  ! nrec /= 0
         
         iret = nf90_inquire(ncID, nVariables=n)
         IF (n /= out_nvar) THEN
            iret = nf90_close(ncID)
            IF (myid == 0) PRINT *, '  ABORTING: Incompatible Netcdf File',n,out_nvar
            CALL appl_abort(0)
         ELSE
            DO n = 1, outProg%count
               xnm = outProg%list(n)%name
               iret = nf90_inquire_variable(ncID, n, name=xnm)
            END DO
            DO n = 1, outVector%count
               xnm = outVector%list(n)%name
               iret = nf90_inquire_variable(ncId, n, name=xnm)
            END DO
            DO n = 1,outDiag%count
               xnm = outDiag%list(n)%name
               iret = nf90_inquire_variable(ncID, n, name=xnm)
            END DO
            DO n = 1, outAxes%count
               xnm = outAxes%list(n)%name
               iret = nf90_inquire_variable(ncID, n, name=xnm)
            END DO
            DO n = 1, outDerived%count
               xnm = outDerived%list(n)%name
               iret = nf90_inquire_variable(ncID, n, name=xnm)
            END DO
            DO n = 1, outPS%count
               xnm = outPS%list(n)%name
               iret = nf90_inquire_variable(ncID, n, name=xnm)
            END DO
            iret = nf90_sync(ncID)
         END IF

      END IF ! nrec
      
    END SUBROUTINE define_nc

    !
    ! ------------------------------------------------------------------------------------
    ! Subroutine defvar_loop: loops over the output variables and performs the variable
    !                         definition cycle
    !
    SUBROUTINE defvar_loop(ncID,varInst,out_nvar)
      INTEGER, INTENT(in) :: ncID
      TYPE(FieldArray), INTENT(in) :: varInst
      INTEGER, INTENT(inout) :: out_nvar
      
      INTEGER :: nvar, n, VarID, iret
      CHARACTER(len=50) :: name
      CHARACTER(len=100) :: long_name
      CHARACTER(len=10) :: unit,dim

      ! Check that the variable array is initialized
      IF (varInst%Initialized) THEN      
         nvar = varInst%count
         DO n = 1,nvar            
            name = varInst%list(n)%name
            long_name = varInst%list(n)%long_name
            unit = varInst%list(n)%unit
            dim = varInst%list(n)%dimension
            
            ! Determine dimension and define variable
            CALL define_variable(ncID,dim,name,VarID)
            
            iret = nf90_put_att(ncID,VarID,'longname',long_name)
            iret = nf90_put_att(ncID,VarID,'units'   ,unit)
         END DO
         out_nvar = out_nvar + nvar
      END IF
         
    END SUBROUTINE defvar_loop
      
    !
    ! ------------------------------------------------------------------------------------
    ! Subroutine define_variable: defines the output variables for netcdf
    !
    SUBROUTINE define_variable(ncID,dim,name,VarID)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=10), INTENT(in) :: dim
      CHARACTER(len=50), INTENT(in) :: name
      INTEGER, INTENT(out) :: VarID
      
      INTEGER :: iret
      
      SELECT CASE(dim)
      ! Grid axis dimensions + for simple time series
      CASE ('time')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,timeID  ,VarID)
      CASE ('zt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,ztID    ,VarID)
      CASE ('zm')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,zmID    ,VarID)
      CASE ('xt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,xtID    ,VarID)
      CASE ('xm')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,xmID    ,VarID)
      CASE ('yt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,ytID    ,VarID)
      CASE ('ym')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,ymID    ,VarID)
      CASE ('aea')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,aeaID   ,VarID)
      CASE ('aeb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,aebID   ,VarID)
      CASE ('cla')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,claID   ,VarID)
      CASE ('clb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,clbID   ,VarID)
      CASE ('prc')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,prcID   ,VarID)
      CASE ('ice')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,iceID   ,VarID)
      ! //
      ! Binned 3d output   
      CASE ('ttttaea')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttaea,VarID)
      CASE ('ttttaeb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttaeb,VarID)
      CASE ('ttttcla')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttcla,VarID)
      CASE ('ttttclb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttclb,VarID)
      CASE ('ttttprc')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttprc,VarID)
      CASE ('ttttice')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttttice,VarID)
      ! //
      ! Regular 3d output   
      CASE ('tttt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tttt,VarID)
      CASE ('mttt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_mttt,VarID)
      CASE ('tmtt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tmtt,VarID)
      CASE ('ttmt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ttmt,VarID)
      ! //
      ! Binned time series
      CASE ('taea')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_taea,VarID)
      CASE ('taeb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_taeb,VarID)
      CASE ('tcla')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tcla,VarID)
      CASE ('tclb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tclb,VarID)
      CASE ('tprc')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tprc,VarID)
      CASE ('tice')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_tice,VarID)
      !//
      ! Binned time-height profiles   
      CASE ('zttaea')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttaea,VarID)
      CASE ('zttaeb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttaeb,VarID)
      CASE ('zttcla')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttcla,VarID)
      CASE ('zttclb')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttclb,VarID)
      CASE ('zttprc')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttprc,VarID)
      CASE ('zttice')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zttice,VarID)
      ! //
      ! Regular time-height profiles
      CASE ('ztt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_ztt,VarID)
      CASE ('zmt')
         iret = nf90_def_var(ncID,name,NF90_FLOAT,dim_zmt,VarID)
         
      CASE DEFAULT
         IF (myid == 0) PRINT *, '  ABORTING: NCIO: Bad dimensional information ',trim(dim)
         CALL appl_abort(0)
      END SELECT
            
    END SUBROUTINE define_variable

    !
    !
    !
    !
    SUBROUTINE write_nc_0d(ncID,name,var,beg)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var
      INTEGER, INTENT(in) :: beg(1)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(ncID, name, VarID)
      iret = nf90_put_var(ncID, VarID, var, start=beg)     
    END SUBROUTINE write_nc_0d

    SUBROUTINE write_nc_1d(ncID,name,var,beg,icnt)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:)
      INTEGER, INTENT(in) :: beg(:)             ! Use automatic arrays to facilitate both axis variables and time-height profile variables
      INTEGER, OPTIONAL, INTENT(in) :: icnt(:)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(ncID, name, VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(ncID, VarID, var, start=beg, count=icnt)     
      ELSE
         iret = nf90_put_var(ncID, VarID, var, start=beg)
      END IF
    END SUBROUTINE write_nc_1d
    
    SUBROUTINE write_nc_2d(ncID,name,var,beg,icnt)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:)
      INTEGER, INTENT(in) :: beg(3)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(3)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(ncID, name, VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(ncID, VarID, var, start=beg, count=icnt)     
      ELSE
         iret = nf90_put_var(ncID, VarID, var, start=beg)
      END IF
    END SUBROUTINE write_nc_2d

    SUBROUTINE write_nc_3d(ncID,name,var,beg,icnt)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:,:)
      INTEGER, INTENT(in) :: beg(4)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(4)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(ncID, name, VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(ncID, VarID, var, start=beg, count=icnt)     
      ELSE
         iret = nf90_put_var(ncID, VarID, var, start=beg)
      END IF
    END SUBROUTINE write_nc_3d

    SUBROUTINE write_nc_4d(ncID,name,var,beg,icnt)
      INTEGER, INTENT(in) :: ncID
      CHARACTER(len=*), INTENT(in) :: name
      REAL, INTENT(in) :: var(:,:,:,:)
      INTEGER, INTENT(in) :: beg(5)
      INTEGER, OPTIONAL, INTENT(in) :: icnt(5)
      INTEGER :: iret, VarID
      iret = nf90_inq_varid(ncID, name, VarID)
      IF (PRESENT(icnt)) THEN
         iret = nf90_put_var(ncID, VarID, var, start=beg, count=icnt)     
      ELSE
         iret = nf90_put_var(ncID, VarID, var, start=beg)
      END IF
    END SUBROUTINE write_nc_4d
    
    
    !
    ! ----------------------------------------------------------------------
    ! Subroutine define_nc_cs: Defines the structure of a column statistics nc file
    !
    SUBROUTINE define_nc_cs(ncID, nRec, n2, n3, level, rad_level, spec_list, nspec )
      INTEGER, INTENT (in)    :: ncID, n2, n3, level, rad_level, nspec
      INTEGER, INTENT (inout) :: nRec ! nRec=0 means new files
      CHARACTER(LEN=3), INTENT (in) :: spec_list(nspec) ! SALSA species (e.g. SO4, Org,...)
      
      INTEGER, SAVE :: timeID=0, xtID=0, ytID=0
      INTEGER, SAVE :: dim_ttt(3)
      CHARACTER(LEN=7) nam
      INTEGER :: iret, VarID, ss
      
      IF (nRec == 0) THEN
         iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
         iret = nf90_def_dim(ncID, 'xt', n2, xtID)
         iret = nf90_def_dim(ncID, 'yt', n3, ytID)
         
         dim_ttt = [xtID,ytID,timeID]
         
         iret = nf90_def_var(ncID,'time',NF90_FLOAT,timeID  ,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','Time')
         iret = nf90_put_att(ncID,VarID,'units','s')
         
         iret = nf90_def_var(ncID,'xt',NF90_FLOAT,xtID    ,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','xt')
         iret = nf90_put_att(ncID,VarID,'units','m')
         
         iret = nf90_def_var(ncID,'yt',NF90_FLOAT,ytID    ,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','yt')
         iret = nf90_put_att(ncID,VarID,'units','m')
         
         iret = nf90_def_var(ncID,'lwp',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','LWP')
         iret = nf90_put_att(ncID,VarID,'units','kg/m2')
         
         iret = nf90_def_var(ncID,'rwp',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','RWP')
         iret = nf90_put_att(ncID,VarID,'units','kg/m2')
         
         iret = nf90_def_var(ncID,'Nc',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','Number of cloud droplets')
         iret = nf90_put_att(ncID,VarID,'units','m-3')
         
         iret = nf90_def_var(ncID,'Nr',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','Number of rain drops')
         iret = nf90_put_att(ncID,VarID,'units','m-3')
         
         iret = nf90_def_var(ncID,'nccnt',NF90_INT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','NCCNT')
         iret = nf90_put_att(ncID,VarID,'units','nccnt')
         
         iret = nf90_def_var(ncID,'nrcnt',NF90_INT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','nrcnt')
         iret = nf90_put_att(ncID,VarID,'units','nrcnt')
         
         iret = nf90_def_var(ncID,'zb',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','zb')
         iret = nf90_put_att(ncID,VarID,'units','m')
         
         iret = nf90_def_var(ncID,'zc',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','zc')
         iret = nf90_put_att(ncID,VarID,'units','m')
         
         iret = nf90_def_var(ncID,'zi1',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','zi1_bar')
         iret = nf90_put_att(ncID,VarID,'units','m')
         
         iret = nf90_def_var(ncID,'lmax',NF90_FLOAT,dim_ttt,VarID)
         iret = nf90_put_att(ncID,VarID,'longname','lmax')
         iret = nf90_put_att(ncID,VarID,'units','lmax')
         
         ! Can add: maximum/minimum vertical velocities and their variances,
         ! surface heat and humidity fluxes, buoyancy statistics,...
         
         IF (rad_level == 3) THEN
            iret = nf90_def_var(ncID,'albedo',NF90_FLOAT,dim_ttt,VarID)
            iret = nf90_put_att(ncID,VarID,'longname','Albedo')
            iret = nf90_put_att(ncID,VarID,'units','1')
         END IF
         
         IF (level >= 4) THEN
            ! Aerosol and water removal
            DO ss = 1, nspec
               nam = 'rm'//trim(spec_list(ss))//'dr'
               iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
               iret = nf90_put_att(ncID,VarID,'longname',nam)
               iret = nf90_put_att(ncID,VarID,'units','kg/m2s')
               
               nam = 'rm'//trim(spec_list(ss))//'cl'
               iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
               iret = nf90_put_att(ncID,VarID,'longname',nam)
               iret = nf90_put_att(ncID,VarID,'units','kg/m2s')
               
               nam = 'rm'//trim(spec_list(ss))//'pr'
               iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
               iret = nf90_put_att(ncID,VarID,'longname',nam)
               iret = nf90_put_att(ncID,VarID,'units','kg/m2s')
            END DO
         END IF
         
         IF (level>=5) THEN
            ! Ice 
            iret = nf90_def_var(ncID,'iwp',NF90_FLOAT,dim_ttt,VarID)
            iret = nf90_put_att(ncID,VarID,'longname','IWP')
            iret = nf90_put_att(ncID,VarID,'units','kg/ms')
            
            iret = nf90_def_var(ncID,'Ni',NF90_FLOAT,dim_ttt,VarID)
            iret = nf90_put_att(ncID,VarID,'longname','Number of ice particles')
            iret = nf90_put_att(ncID,VarID,'units','m-3')
            
            iret = nf90_def_var(ncID,'nicnt',NF90_INT,dim_ttt,VarID)
            iret = nf90_put_att(ncID,VarID,'longname','nicnt')
            iret = nf90_put_att(ncID,VarID,'units','nicnt')
            
            ! Aerosol and water removal with ice
            DO ss = 1,nspec
               nam = 'rm'//trim(spec_list(ss))//'ic'
               iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
               iret = nf90_put_att(ncID,VarID,'longname',nam)
               iret = nf90_put_att(ncID,VarID,'units','kg/m2s')
            END DO
         END IF
         
         IF (level==3) THEN
            ! Surface precipitation for levels 3
            iret = nf90_def_var(ncID,'prcp',NF90_FLOAT,dim_ttt,VarID)
            iret = nf90_put_att(ncID,VarID,'longname','prcp')
            iret = nf90_put_att(ncID,VarID,'units','kg/m2s')
         END IF
         
         iret = nf90_enddef(ncID)
         iret = nf90_sync(ncID)
         nRec = 1
      END IF
      
    END SUBROUTINE define_nc_cs

 !
 ! ----------------------------------------------------------------------
 ! FUNCTIONS FOR READING AEROSOL SIZE DISTRIBUTIONS FROM A NETCDF FILE
 !
 SUBROUTINE open_aero_nc(ncid,nc_levs,nc_nspec,nc_nmod)
    IMPLICIT NONE

    INTEGER, INTENT(out) :: ncid,nc_levs,nc_nspec,nc_nmod
    INTEGER :: iret, did

    ! Open file
    iret = nf90_open('aerosol_in.nc',NF90_NOWRITE,ncid)

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
