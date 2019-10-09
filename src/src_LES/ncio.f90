MODULE ncio

  USE netcdf
  USE mpi_interface, ONLY : appl_abort, myid, pecount, wrxid, wryid

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: open_nc, define_nc, define_nc_cs, &
            open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_aero_nc

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
    CHARACTER (len=*), INTENT (in) :: fname, ename
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
  SUBROUTINE define_nc(ncID, nRec, nVar, sx, n1, n2, n3, &
                       inae_a,incld_a,inprc,             &
                       inae_b,incld_b,inice_a,inice_b,insnw         )

    INTEGER, INTENT (in)           :: nVar, ncID
    INTEGER, OPTIONAL, INTENT (in) :: n1, n2, n3
    ! Juha: Added
    INTEGER, OPTIONAL, INTENT(in)  :: inae_a,incld_a,inprc, &
                                      inae_b,incld_b,       &
                                      inice_a,inice_b,insnw
    ! --
    INTEGER, INTENT (inout)        :: nRec
    CHARACTER (len=7), INTENT (in) :: sx(nVar)

    INTEGER, SAVE :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,            &
                     dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0,  &
                     dim_tt(2)  = 0, dim_mt(2)  = 0

    ! Juha: added
    INTEGER, SAVE :: aeaID=0, claID=0, aebID=0, clbID=0, prcID=0,   &
                     icaID=0, icbID=0,snowID=0,                     &
                     dim_ttttaea(5) = 0, dim_ttttcla(5) = 0,    &
                     dim_ttttaeb(5) = 0, dim_ttttclb(5) = 0,    &
                     dim_ttttprc(5) = 0,                        &
                     dim_ttttica(5) = 0, dim_tttticb(5) = 0,    &
                     dim_ttttsnw(5) = 0,                        &
                     dim_ttaea(2) = 0, dim_ttcla(2) = 0,        &
                     dim_ttaeb(2) = 0, dim_ttclb(2) = 0,        &
                     dim_ttprc(2) = 0,                          &
                     dim_ttica(2) = 0, dim_tticb(2) = 0,        &
                     dim_ttsnw(2) = 0,                          &
                     dim_ttztaea(3) = 0, dim_ttztcla(3) = 0,    &
                     dim_ttztaeb(3) = 0, dim_ttztclb(3) = 0,    &
                     dim_ttztprc(3) = 0,                        &
                     dim_ttztica(3) = 0, dim_ttzticb(3) = 0,    &
                     dim_ttztsnw(3) = 0
    !--

    CHARACTER (len=7) :: xnm
    INTEGER :: iret, n, VarID


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
       ! IF this is analysis file, dont write binned output by default!
       ! --------------------------------------------------------------
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
       IF (present(inice_a)) THEN
          iret = nf90_def_dim(ncID, 'ica', inice_a, icaID)
       END IF
       IF (present(inice_b)) THEN
          iret = nf90_def_dim(ncID, 'icb', inice_b, icbID)
       END IF
       IF (present(insnw)) THEN
          iret = nf90_def_dim(ncID, 'snw', insnw, snowID)
       END IF

       dim_tt = (/ztID,timeID/)
       dim_mt = (/zmID,timeID/)
       dim_tttt= (/ztID,xtID,ytID,timeID/)  ! thermo point
       dim_mttt= (/zmID,xtID,ytID,timeID/)  ! zpoint
       dim_tmtt= (/ztID,xmID,ytID,timeID/)  ! upoint
       dim_ttmt= (/ztID,xtID,ymID,timeID/)  ! ypoint

       ! Juha: dimension environments for size distribution variables
       dim_ttttaea = (/ztID,xtID,ytID,aeaID,timeID/)
       dim_ttttaeb = (/ztID,xtID,ytID,aebID,timeID/)
       dim_ttttcla = (/ztID,xtID,ytID,claID,timeID/)
       dim_ttttclb = (/ztID,xtID,ytID,clbID,timeID/)
       dim_ttttprc = (/ztID,xtID,ytID,prcID,timeID/)

       ! Jaakko: ice & snow
       dim_ttttica = (/ztID,xtID,ytID,icaID,timeID/)
       dim_tttticb = (/ztID,xtID,ytID,icbID,timeID/)
       dim_ttttsnw = (/ztID,xtID,ytID,snowID,timeID/)
       ! ---
       ! Zubair: dimension environments for avegare size distribution variables per bin - ts files
       dim_ttaea = (/aeaID,timeID/)
       dim_ttaeb = (/aebID,timeID/)
       dim_ttcla = (/claID,timeID/)
       dim_ttclb = (/clbID,timeID/)
       dim_ttprc = (/prcID,timeID/)
       ! Jaakko:
       dim_ttica = (/icaID,timeID/)
       dim_tticb = (/icbID,timeID/)
       dim_ttsnw = (/snowID,timeID/)
       ! Zubair: dimension environments for avegare size distribution variables per bin - ps files
       dim_ttztaea = (/ztID,aeaID,timeID/)
       dim_ttztaeb = (/ztID,aebID,timeID/)
       dim_ttztcla = (/ztID,claID,timeID/)
       dim_ttztclb = (/ztID,clbId,timeID/)
       dim_ttztprc = (/ztID,prcID,timeID/)
       ! Jaakko
       dim_ttztica = (/ztID,icaID,timeID/)
       dim_ttzticb = (/ztID,icbId,timeID/)
       dim_ttztsnw = (/ztID,snowID,timeID/)

       DO n = 1, nVar
          SELECT CASE(trim(ncinfo(2,sx(n))))
          CASE ('time')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,timeID  ,VarID)
          CASE ('zt')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,ztID    ,VarID)
          CASE ('zm')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,zmID    ,VarID)
          CASE ('xt')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,xtID    ,VarID)
          CASE ('xm')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,xmID    ,VarID)
          CASE ('yt')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,ytID    ,VarID)
          CASE ('ym')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,ymID    ,VarID)
          ! Juha: added for size distributions
          CASE ('aea')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,aeaID   ,VarID)
          CASE ('aeb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,aebID   ,VarID)
          CASE ('cla')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,claID   ,VarID)
          CASE ('clb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,clbID   ,VarID)
          CASE ('prc')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,prcID   ,VarID)
          !Jaakko added for ice and snow
          CASE ('ica')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,icaID   ,VarID)
          CASE ('icb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,icbID   ,VarID)
          CASE ('snow')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,snowID   ,VarID)
          !Juha added
          CASE ('ttttaea')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaea,VarID)
          CASE ('ttttaeb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaeb,VarID)
          CASE ('ttttcla')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttcla,VarID)
          CASE ('ttttclb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttclb,VarID)
          CASE ('ttttprc')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttprc,VarID)
          !Jaakko added
          CASE ('ttttica')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttica,VarID)
          CASE ('tttticb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttticb,VarID)
          CASE ('ttttsnw')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttsnw,VarID)
          ! ---
          CASE ('tttt')
             IF (present(n2) .AND. present(n3)) THEN
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttt,VarID)
             ELSE
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             END IF
          CASE ('mttt')
             IF (present(n2) .AND. present(n3)) THEN
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mttt,VarID)
             ELSE
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             END IF
          CASE ('tmtt')
             IF (present(n2) .AND. present(n3)) THEN
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tmtt,VarID)
             ELSE
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             END IF
          CASE ('ttmt')
             IF (present(n2) .AND. present(n3)) THEN
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttmt,VarID)
             ELSE
                iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mt,VarID)
             END IF
          CASE ('ttaea')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttaea,VarID)
          CASE ('ttaeb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttaeb,VarID)
          CASE ('ttcla')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttcla,VarID)
          CASE ('ttclb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttclb,VarID)
          CASE ('ttprc')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttprc,VarID)
          CASE ('ttica')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttica,VarID)
          CASE ('tticb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tticb,VarID)
          CASE ('ttsnw')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttsnw,VarID)
          CASE ('ttztaea')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaea,VarID)
          CASE ('ttztaeb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaeb,VarID)
          CASE ('ttztcla')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztcla,VarID)
          CASE ('ttztclb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztclb,VarID)
          CASE ('ttztprc')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztprc,VarID)
          CASE ('ttztica')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztica,VarID)
          CASE ('ttzticb')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttzticb,VarID)
          CASE ('ttztsnw')
             iret = nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztsnw,VarID)
          CASE DEFAULT
             IF (myid == 0) PRINT *, '  ABORTING: NCIO: Bad dimensional information ',trim(ncinfo(2,sx(n)))
             CALL appl_abort(0)
       END SELECT
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,sx(n)))
       iret = nf90_put_att(ncID,VarID,'units'   ,ncinfo(1,sx(n)))
    END DO
    iret = nf90_enddef(ncID)
    iret = nf90_sync(ncID)
    nRec = 1
 ELSE
    iret = nf90_inquire(ncID, nVariables=n)
    IF (n /= nVar) THEN
       iret = nf90_close(ncID)
       IF (myid == 0) PRINT *, '  ABORTING: Incompatible Netcdf File',n,nVar
       CALL appl_abort(0)
    ELSE
       DO n = 1, nVar
          xnm = sx(n)
          iret = nf90_inquire_variable(ncID, n, name=xnm)
       END DO
       iret = nf90_sync(ncID)
    END IF
 END IF

 END SUBROUTINE define_nc
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

       dim_ttt = (/xtID,ytID,timeID/)

       iret = nf90_def_var(ncID,'time',NF90_FLOAT,timeID  ,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'time'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'time'))

       iret = nf90_def_var(ncID,'xt',NF90_FLOAT,xtID    ,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'xt'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'xt'))

       iret = nf90_def_var(ncID,'yt',NF90_FLOAT,ytID    ,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'yt'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'yt'))

       iret = nf90_def_var(ncID,'lwp',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'lwp'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'lwp'))

       iret = nf90_def_var(ncID,'rwp',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'rwp'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'rwp'))

       iret = nf90_def_var(ncID,'Nc',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Nc'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'Nc'))

       iret = nf90_def_var(ncID,'Nr',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Nr'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'Nr'))

       iret = nf90_def_var(ncID,'nccnt',NF90_INT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nccnt'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'nccnt'))

       iret = nf90_def_var(ncID,'nrcnt',NF90_INT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nrcnt'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'nrcnt'))

       iret = nf90_def_var(ncID,'zb',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zb'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'zb'))

       iret = nf90_def_var(ncID,'zc',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zc'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'zc'))

       iret = nf90_def_var(ncID,'zi1',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zi1_bar'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'zi1_bar'))

       iret = nf90_def_var(ncID,'lmax',NF90_FLOAT,dim_ttt,VarID)
       iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'lmax'))
       iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'lmax'))

       ! Can add: maximum/minimum vertical velocities and their variances,
       ! surface heat and humidity fluxes, buoyancy statistics,...

       IF (rad_level == 3) THEN
          iret = nf90_def_var(ncID,'albedo',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'albedo'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'albedo'))
       END IF

       IF (level >= 4) THEN
          ! Aerosol and water removal
          DO ss = 1, nspec
             nam = 'rm'//trim(spec_list(ss))//'dr'
             iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
             iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
             iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))

             nam = 'rm'//trim(spec_list(ss))//'cl'
             iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
             iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
             iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))

             nam = 'rm'//trim(spec_list(ss))//'pr'
             iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
             iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
             iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))
          END DO
       END IF

       IF (level>=5) THEN
          ! Ice and snow
          iret = nf90_def_var(ncID,'iwp',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'iwp'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'iwp'))

          iret = nf90_def_var(ncID,'swp',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'swp'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'swp'))

          iret = nf90_def_var(ncID,'Ni',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Ni'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'Ni'))

          iret = nf90_def_var(ncID,'Ns',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Ns'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'Ns'))

          iret = nf90_def_var(ncID,'nicnt',NF90_INT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nicnt'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'nicnt'))

          iret = nf90_def_var(ncID,'nscnt',NF90_INT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nscnt'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'nscnt'))

          ! Aerosol and water removal with ice and snow
          DO ss = 1,nspec
             nam = 'rm'//trim(spec_list(ss))//'ic'
             iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
             iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
             iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))

             nam = 'rm'//trim(spec_list(ss))//'sn'
             iret = nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
             iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
             iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))
          END DO
       END IF

       IF (level==3) THEN
          ! Surface precipitation for levels 3
          iret = nf90_def_var(ncID,'prcp',NF90_FLOAT,dim_ttt,VarID)
          iret = nf90_put_att(ncID,VarID,'longname',ncinfo(0,'prcp'))
          iret = nf90_put_att(ncID,VarID,'units',ncinfo(1,'prcp'))
       END IF

       iret = nf90_enddef(ncID)
       iret = nf90_sync(ncID)
       nRec = 1
    END IF

 END SUBROUTINE define_nc_cs
 !
 ! ----------------------------------------------------------------------
 ! Subroutine nc_info: Gets long_name, units and dimension info given a
 ! short name.
 !
 CHARACTER (len=80) FUNCTION ncinfo(itype,short_name)

    CHARACTER (len=40) :: v_lnm ='scalar xx mixing ratio                  '

    INTEGER, INTENT (in) :: itype
    CHARACTER (len=*), INTENT (in) :: short_name

    INTEGER :: scalar_number

    SELECT CASE (trim(short_name))
       CASE ('sxx')
          READ(short_name(2:3),'(i2.2)') scalar_number
          WRITE(v_lnm(8:9),'(i2.2)') scalar_number
          IF (itype == 0) ncinfo = v_lnm
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('time')
          IF (itype == 0) ncinfo = 'Time'
          IF (itype == 1) ncinfo = 's'
          IF (itype == 2) ncinfo = 'time'
       CASE('zt')
          IF (itype == 0) ncinfo = 'Vertical displacement of cell centers'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'zt'
       CASE('zm')
          IF (itype == 0) ncinfo = 'Vertical displacement of cell edges'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'zm'
       CASE('xt')
          IF (itype == 0) ncinfo = 'East-west displacement of cell centers'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'xt'
       CASE('xm')
          IF (itype == 0) ncinfo = 'East-west displacement of cell edges'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'xm'
       CASE('yt')
          IF (itype == 0) ncinfo = 'North-south displacement of cell centers'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'yt'
       CASE('ym')
          IF (itype == 0) ncinfo = 'North-south displacement of cell edges'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ym'
       ! Juha: added for SALSA
       CASE('aea')
          IF (itype == 0) ncinfo = 'Aerosol size bins, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'aea'
       CASE('aeb')
          IF (itype == 0) ncinfo = 'Aerosol size bins, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'aeb'
       CASE('cla')
          IF (itype == 0) ncinfo = 'Cloud droplet size bins, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'cla'
       CASE('clb')
          IF (itype == 0) ncinfo = 'Cloud droplet size bins, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'clb'
       CASE('prc')
          IF (itype == 0) ncinfo = 'Precipitation size bins'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'prc'
       CASE('ica')
          IF (itype == 0) ncinfo = 'Ice cloud droplet size bins, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ica'
       CASE('icb')
          IF (itype == 0) ncinfo = 'Ice cloud droplet size bins, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'icb'
       CASE('snw')
          IF (itype == 0) ncinfo = 'Snow size bins'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'snow'
       !----
       CASE('u0')
          IF (itype == 0) ncinfo = 'Geostrophic zonal wind'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'zt'
       CASE('v0')
          IF (itype == 0) ncinfo = 'Geostrophic meridional wind'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'zt'
       CASE('dn0')
          IF (itype == 0) ncinfo = 'Base-state density'
          IF (itype == 1) ncinfo = 'kg/m^3'
          IF (itype == 2) ncinfo = 'zt'
       CASE('u')
          IF (itype == 0) ncinfo = 'Zonal wind'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'mttt'
       CASE('v')
          IF (itype == 0) ncinfo = 'Meridional wind'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'tmtt'
       CASE('w')
          IF (itype == 0) ncinfo = 'Vertical velocity'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('thl')
          IF (itype == 0) ncinfo = 'Liquid water potential temperature'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('thil')
          IF (itype == 0) ncinfo = 'Ice-Liquid water potential temperature'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('theta')
          IF (itype == 0) ncinfo = 'Potential temperature'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('p')
          IF (itype == 0) ncinfo = 'Pressure'
          IF (itype == 1) ncinfo = 'Pa'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('q')
          IF (itype == 0) ncinfo = 'Total water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('l')
          IF (itype == 0) ncinfo = 'Liquid water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('r')
          IF (itype == 0) ncinfo = 'Rain-water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('i')
          IF (itype == 0) ncinfo = 'Ice mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('s')
          IF (itype == 0) ncinfo = 'Snow mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('n')
          IF (itype == 0) ncinfo = 'Rain-drop number mixing ratio'
          IF (itype == 1) ncinfo = '#/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('stke')
          IF (itype == 0) ncinfo = 'Sub-filter scale TKE'
          IF (itype == 1) ncinfo = 'J/kg'
          IF (itype == 2) ncinfo = 'mttt'
       CASE('cfl')
          IF (itype == 0) ncinfo = 'Courant number'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'time'
       CASE('maxdiv')
          IF (itype == 0) ncinfo = 'Maximum divergence'
          IF (itype == 1) ncinfo = '1/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('zi1_bar')
          IF (itype == 0) ncinfo = 'Height of maximum theta gradient'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('zi2_bar')
          IF (itype == 0) ncinfo = 'Height of maximum theta variance'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('zi3_bar')
          IF (itype == 0) ncinfo = 'Height of minimum buoyancy flux'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('vtke')
          IF (itype == 0) ncinfo = 'Vertical integral of total TKE'
          IF (itype == 1) ncinfo = 'kg/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('sfcbflx')
          IF (itype == 0) ncinfo = 'Surface Buoyancy Flux'
          IF (itype == 1) ncinfo = 'm/s^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('wmax')
          IF (itype == 0) ncinfo = 'Maximum vertical velocity'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('tsrf')
          IF (itype == 0) ncinfo = 'Surface temperature'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'time'
       CASE('ustar')
          IF (itype == 0) ncinfo = 'Surface friction velocity'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('shf_bar')
          IF (itype == 0) ncinfo = 'Sensible heat flux'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('lhf_bar')
          IF (itype == 0) ncinfo = 'Latent heat flux'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('zi_bar')
          IF (itype == 0) ncinfo = 'Height of maximum total water mixing ratio gradient'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('lwp_bar','lwp')
          IF (itype == 0) ncinfo = 'Liquid-water path'
          IF (itype == 1) ncinfo = 'kg/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('lwp_var')
          IF (itype == 0) ncinfo = 'Liquid-water path variance'
          IF (itype == 1) ncinfo = 'kg^2/m^4'
          IF (itype == 2) ncinfo = 'time'
       CASE('iwp_bar','iwp')
          IF (itype == 0) ncinfo = 'Ice-water path'
          IF (itype == 1) ncinfo = 'kg/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('iwp_var')
          IF (itype == 0) ncinfo = 'Ice-water path variance'
          IF (itype == 1) ncinfo = 'kg/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('zc')
          IF (itype == 0) ncinfo = 'Cloud-top height'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('zb')
          IF (itype == 0) ncinfo = 'Cloud-base height'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('cfrac')
          IF (itype == 0) ncinfo = 'Cloud fraction'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'time'
       CASE('lmax')
          IF (itype == 0) ncinfo = 'Maximum liquid water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       case('imax')
          if (itype == 0) ncinfo = 'Maximum ice water mixing ratio'
          if (itype == 1) ncinfo = 'kg/kg'
          if (itype == 2) ncinfo = 'time'
       case('smax')
          if (itype == 0) ncinfo = 'Maximum snow water mixing ratio'
          if (itype == 1) ncinfo = 'kg/kg'
          if (itype == 2) ncinfo = 'time'
       CASE('albedo')
          IF (itype == 0) ncinfo = 'Reflected (TOA) shortwave radiation'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'time'
       CASE('rwp_bar','rwp')
          IF (itype == 0) ncinfo = 'Rain-water path'
          IF (itype == 1) ncinfo = 'kg/m^2'
          IF (itype == 2) ncinfo = 'time'
       case('swp_bar','swp')
          if (itype == 0) ncinfo = 'Snow-water path'
          if (itype == 1) ncinfo = 'kg/m^2'
          if (itype == 2) ncinfo = 'time'
       CASE('prcp')
          IF (itype == 0) ncinfo = 'Surface precipitation rate'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'time'
       case('sprcp')
          if (itype == 0) ncinfo = 'Surface snow precipitation rate'
          if (itype == 1) ncinfo = 'W/m^2'
          if (itype == 2) ncinfo = 'time'
       CASE('prcp_bc')
          IF (itype == 0) ncinfo = 'Below cloud precipitation rate'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'time'
       CASE('pfrac')
          IF (itype == 0) ncinfo = 'Surface precipitation fraction'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'time'
       CASE('sfrac')
          IF (itype == 0) ncinfo = 'Surface snow precipitation fraction'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'time'
       CASE('CCN')
          IF (itype == 0) ncinfo = 'Cloud condensation nuclei'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('nrain')
          IF (itype == 0) ncinfo = 'Conditionally sampled rain number mixing ratio'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('nrcnt')
          IF (itype == 0) ncinfo = 'Rain cell counts'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'time'
      case('nscnt')
          if (itype == 0) ncinfo = 'Snow cell counts'
          if (itype == 1) ncinfo = '#'
          if (itype == 2) ncinfo = 'time'
       CASE('nccnt')
          IF (itype == 0) ncinfo = 'Cloud cell counts'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'time'
       CASE('nicnt')
          IF (itype == 0) ncinfo = 'Ice cell counts'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'time'
       CASE('SS_max')
          IF (itype==0) ncinfo = 'Maximum supersaturation'
          IF (itype==1) ncinfo = '%'
          IF (itype==2) ncinfo = 'time'
       CASE('SSi_max')
          IF (itype==0) ncinfo = 'Maximum supersaturation over ice'
          IF (itype==1) ncinfo = '%'
          IF (itype==2) ncinfo = 'time'
       !
       !
       ! SALSA temporal statistics
       CASE('Nc_ic')
          IF (itype == 0) ncinfo = 'In-cloud CDNC'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Nca_ica')
          IF (itype==0) ncinfo = 'In-cloud CDNC (a bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Ncb_icb')
          IF (itype==0) ncinfo = 'In-cloud CDNC (b bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Na_oc')
          IF (itype == 0) ncinfo = 'Aerosol number concentration outside clouds'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Na_int')
          IF (itype == 0) ncinfo = 'In-cloud interstitial aerosol number concentration'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Naa_int')
          IF (itype==0) ncinfo = 'In-cloud interstitial aerosol number concentration (a bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Nab_int')
          IF (itype==0) ncinfo = 'In-cloud interstitial aerosol number concentration (b bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Ni_ic')
          IF (itype == 0) ncinfo = 'Ice number concentration in liquid clouds'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ni_ii')
          IF (itype == 0) ncinfo = 'Ice number concentration in icy regions'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Nia_iia')
          IF (itype==0) ncinfo = 'Ice number concentration in icy regions (a bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Nib_iib')
          IF (itype==0) ncinfo = 'Ice number concentration in icy regions (b bins)'
          IF (itype==1) ncinfo = 'kg^-1'
          IF (itype==2) ncinfo = 'time'
       CASE('Ni_is')
          IF (itype == 0) ncinfo = 'Ice number concentration in snowy regions'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ns_ic')
          IF (itype == 0) ncinfo = 'Snow number concentration in liquid clouds'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ns_ii')
          IF (itype == 0) ncinfo = 'Snow number concentration in icy regions'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ns_is')
          IF (itype == 0) ncinfo = 'Snow number concentration in snowy regions'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ra_int')
          IF (itype==0) ncinfo = 'Mean interstitial aerosol wet radius'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Raa_int')
          IF (itype==0) ncinfo = 'Mean interstitial aerosol wet radius (a bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rab_int')
          IF (itype==0) ncinfo = 'Mean interstitial aerosol wet radius (b bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rc_ic')
          IF (itype==0) ncinfo = 'Mean cloud droplet radius'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rca_ica')
          IF (itype==0) ncinfo = 'Mean cloud droplet radius (a bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rcb_icb')
          IF (itype==0) ncinfo = 'Mean cloud droplet radius (b bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Ri_ii')
          IF (itype == 0) ncinfo = 'Mean ice radius in icy regions'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('Ria_iia')
          IF (itype==0) ncinfo = 'Mean ice radius (a bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rib_iib')
          IF (itype==0) ncinfo = 'Mean ice radius (b bins)'
          IF (itype==1) ncinfo = 'm'
          IF (itype==2) ncinfo = 'time'
       CASE('Rs_is')
          IF (itype == 0) ncinfo = 'Mean snow radius in snowy regions'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'time'
       CASE('SO4_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet SO4 mass mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SO4_oc')
          IF (itype == 0) ncinfo = 'Aerosol SO4 mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SO4_int')
          IF (itype == 0) ncinfo = 'SO4 mass mixing ratio in intersitial aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SO4_ii')
          IF (itype==0) ncinfo = 'SO4 mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('SO4_is')
          IF (itype==0) ncinfo = 'SO4 mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('OC_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet OC mass mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('OC_oc')
          IF (itype == 0) ncinfo = 'Aerosol OC mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('OC_int')
          IF (itype == 0) ncinfo = 'OC mass mixing ratio in interstitial aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('OC_ii')
          IF (itype==0) ncinfo = 'OC mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('OC_is')
          IF (itype==0) ncinfo = 'OC mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('BC_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet BC mass mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('BC_oc')
          IF (itype == 0) ncinfo = 'Aerosol BC mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('BC_int')
          IF (itype == 0) ncinfo = 'BC mass mixing ratio in interstitial aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('BC_ii')
          IF (itype==0) ncinfo = 'BC mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('BC_is')
          IF (itype==0) ncinfo = 'BC mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('DU_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet DU mass mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('DU_oc')
          IF (itype == 0) ncinfo = 'Aerosol DU mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('DU_int')
          IF (itype == 0) ncinfo = 'DU mass mixing ration in interstitial aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('DU_ii')
          IF (itype==0) ncinfo = 'DU mass mixing ration in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('DU_is')
          IF (itype==0) ncinfo = 'DU mass mixing ration in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('SS_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet mass mixing ratio of SS'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SS_oc')
          IF (itype == 0) ncinfo = 'Aerosol SS mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SS_int')
          IF (itype == 0) ncinfo = 'SS mass mixing ratio in interstitial particles'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('SS_ii')
          IF (itype==0) ncinfo = 'SS mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('SS_is')
          IF (itype==0) ncinfo = 'SS mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('NH_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet mass mixing ratio of NH3'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NH_oc')
          IF (itype == 0) ncinfo = 'Aerosol NH3 mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NH_int')
          IF (itype == 0) ncinfo = 'NH3 mass mixing ratio in interstitial particles'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NH_ii')
          IF (itype==0) ncinfo = 'NH3 mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('NH_is')
          IF (itype==0) ncinfo = 'NH3 mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('NO_ic')
          IF (itype == 0) ncinfo = 'Cloud droplet mass mixing ratio of NO3'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NO_oc')
          IF (itype == 0) ncinfo = 'Aerosol NO3 mass mixing ratio outside clouds'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NO_int')
          IF (itype == 0) ncinfo = 'NO3 mass mixing ratio in interstitial particles'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'time'
       CASE('NO_ii')
          IF (itype==0) ncinfo = 'NO3 mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('NO_is')
          IF (itype==0) ncinfo = 'NO3 mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('H2O_ic')
          IF (itype==0) ncinfo = 'Cloud droplet mass mixing ratio of H2O'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('H2O_int')
          IF (itype==0) ncinfo = 'H2O mass mixing ratio in interstitial particles'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('H2O_ii')
          IF (itype==0) ncinfo = 'H2O mass mixing ratio in ice'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('H2O_is')
          IF (itype==0) ncinfo = 'H2O mass mixing ratio in snow'
          IF (itype==1) ncinfo = 'kg/kg'
          IF (itype==2) ncinfo = 'time'
       CASE('rmH2Oae','rmH2Odr')
          IF (itype == 0) ncinfo = 'Deposition of H2O with aerosols'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmH2Ocl')
          IF (itype == 0) ncinfo = 'Deposition of H2O with cloud droplets'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmH2Opr')
          IF (itype == 0) ncinfo = 'Deposition of water with rain'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmH2Oic')
          IF (itype == 0) ncinfo = 'Deposition of water with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
        CASE('rmH2Osn')
          IF (itype == 0) ncinfo = 'Deposition of water with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4dr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of SO4'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4cl')
          IF (itype == 0) ncinfo = 'Cloud deposition of SO4'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4pr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of SO4'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4ic')
          IF (itype == 0) ncinfo = 'Deposition of SO4 with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4sn')
          IF (itype == 0) ncinfo = 'Deposition of SO4 with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4wt')
          IF (itype == 0) ncinfo = 'Total wet deposition of SO4'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSO4tt')
          IF (itype == 0) ncinfo = 'Total deposition of SO4'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of OC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCcl')
          IF (itype == 0) ncinfo = 'Cloud deposition of OC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of OC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCic')
          IF (itype == 0) ncinfo = 'Deposition of OC with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCsn')
          IF (itype == 0) ncinfo = 'Deposition of OC with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of OC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmOCtt')
          IF (itype == 0) ncinfo = 'Total deposition of OC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of BC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCcl')
          IF (itype == 0) ncinfo = 'Cloud deposition of BC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of BC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCic')
          IF (itype == 0) ncinfo = 'Deposition of BC with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCsn')
          IF (itype == 0) ncinfo = 'Deposition of BC with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of BC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmBCtt')
          IF (itype == 0) ncinfo = 'Total deposition of BC'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of DU'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUcl')
          IF (itype == 0) ncinfo = 'Cloud deposition of DU'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of DU'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUic')
          IF (itype == 0) ncinfo = 'Deposition of DU with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUsn')
          IF (itype == 0) ncinfo = 'Deposition of DU with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of DU'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmDUtt')
          IF (itype == 0) ncinfo = 'Total deposition of DU'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSSdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of SS'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSScl')
          IF (itype == 0) ncinfo = 'Cloud deposition of SS'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSSpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of SS'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSSic')
          IF (itype == 0) ncinfo = 'Deposition of SS with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSSsn')
          IF (itype == 0) ncinfo = 'Deposition of SS with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSSwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of SS'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmSStt')
          IF (itype == 0) ncinfo = 'Total deposition of SS'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of NH3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHcl')
          IF (itype == 0) ncinfo = 'Cloud deposition of NH3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of NH3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHic')
          IF (itype == 0) ncinfo = 'Deposition of NH3 with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHsn')
          IF (itype == 0) ncinfo = 'Deposition of NH3 with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of NH3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNHtt')
          IF (itype == 0) ncinfo = 'Total deposition of NH3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOdr')
          IF (itype == 0) ncinfo = 'Aerosol deposition of NO3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOcl')
          IF (itype == 0) ncinfo = 'Cloud deposition of NO3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOpr')
          IF (itype == 0) ncinfo = 'Precipitation deposition of NO3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOic')
          IF (itype == 0) ncinfo = 'Deposition of NO3 with ice'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOsn')
          IF (itype == 0) ncinfo = 'Deposition of NO3 with snow'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOwt')
          IF (itype == 0) ncinfo = 'Total wet deposition of NO3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       CASE('rmNOtt')
          IF (itype == 0) ncinfo = 'Total deposition of NO3'
          IF (itype == 1) ncinfo = 'kg/m^2/s'
          IF (itype == 2) ncinfo = 'time'
       ! // SALSA temporal
       CASE('fsttm')
          IF (itype == 0) ncinfo = 'First sample time'
          IF (itype == 1) ncinfo = 's'
          IF (itype == 2) ncinfo = 'time'
       CASE('lsttm')
          IF (itype == 0) ncinfo = 'Last sample time'
          IF (itype == 1) ncinfo = 's'
          IF (itype == 2) ncinfo = 'time'
       CASE('nsmp')
          IF (itype == 0) ncinfo = 'Number of samples'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'time'
       CASE('u_2')
          IF (itype == 0) ncinfo = 'Variance of u wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('v_2')
          IF (itype == 0) ncinfo = 'Variance of v wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('w_2')
          IF (itype == 0) ncinfo = 'Second raw moment of w wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       case('theta_2')
          if (itype == 0) ncinfo = 'Variance of theta'
          if (itype == 1) ncinfo = 'K^2'
          if (itype == 2) ncinfo = 'tttt'
       CASE('w_3')
          IF (itype == 0) ncinfo = 'Third raw moment of w wind'
          IF (itype == 1) ncinfo = 'm^3/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       case('theta_3')
          if (itype == 0) ncinfo = 'Third moment of theta'
          if (itype == 1) ncinfo = 'K^3'
          if (itype == 2) ncinfo = 'tttt'
       CASE('tot_tw')
          IF (itype == 0) ncinfo = 'Total vertical flux of theta'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_tw')
          IF (itype == 0) ncinfo = 'Sub-filter scale vertical flux of theta'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('tot_uw')
          IF (itype == 0) ncinfo = 'Total vertical flux of u-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_uw')
          IF (itype == 0) ncinfo = 'Sub-filter scale vertical flux of u-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('tot_vw')
          IF (itype == 0) ncinfo = 'Total vertical flux of v-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_vw')
          IF (itype == 0) ncinfo = 'SGS vertical flux of v-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('tot_ww')
          IF (itype == 0) ncinfo = 'Total vertical flux of v-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_ww')
          IF (itype == 0) ncinfo = 'SGS vertical flux of w-wind'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('km')
          IF (itype == 0) ncinfo = 'Eddy viscosity'
          IF (itype == 1) ncinfo = 'm^2/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('kh')
          IF (itype == 0) ncinfo = 'Eddy dIFfusivity'
          IF (itype == 1) ncinfo = 'm^2/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('lmbd')
          IF (itype == 0) ncinfo = 'Mixing lengthscale'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('lmbde')
          IF (itype == 0) ncinfo = 'Dissipation lengthscale'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_tke')
          IF (itype == 0) ncinfo = 'Sub-filter scale TKE'
          IF (itype == 1) ncinfo = 'm^2/s^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_boy')
          IF (itype == 0) ncinfo = 'Subfilter Buoyancy production of TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_shr')
          IF (itype == 0) ncinfo = 'Shear production of SGS TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('boy_prd')
          IF (itype == 0) ncinfo = 'Buoyancy production of resolved TKE '
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('shr_prd')
          IF (itype == 0) ncinfo = 'Shear production of resolved TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('trans')
          IF (itype == 0) ncinfo = 'Net transport of resolved TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('diss')
          IF (itype == 0) ncinfo = 'Dissipation rate of resolved TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('dff_u')
          IF (itype == 0) ncinfo = 'u(du/dt) from dIFfusion'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('dff_v')
          IF (itype == 0) ncinfo = 'v(dv/dt) from dIFfusion'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('dff_w')
          IF (itype == 0) ncinfo = 'w(dw/dt) from dIFfusion'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('adv_u')
          IF (itype == 0) ncinfo = 'u(du/dt) from advection'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('adv_v')
          IF (itype == 0) ncinfo = 'v(dv/dt) from advection'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('adv_w')
          IF (itype == 0) ncinfo = 'w(dw/dt) from advection'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('prs_u')
          IF (itype == 0) ncinfo = 'u(du/dt) from pressure'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('prs_v')
          IF (itype == 0) ncinfo = 'v(dv/dt) from pressure'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('prs_w')
          IF (itype == 0) ncinfo = 'w(dw/dt) from pressure'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('prd_uw')
          IF (itype == 0) ncinfo = 'uw shear production'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('storage')
          IF (itype == 0) ncinfo = 'Rate of increase of resolved TKE'
          IF (itype == 1) ncinfo = 'm^2/s^3'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('q_2')
          IF (itype == 0) ncinfo = 'Variance of total water'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('q_3')
          IF (itype == 0) ncinfo = 'Third moment of total water'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('tot_qw')
          IF (itype == 0) ncinfo = 'Total vertical flux of q'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sfs_qw')
          IF (itype == 0) ncinfo = 'Sub-filter scale vertical flux of q'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('rflx')
          IF (itype == 0) ncinfo =  'Total Radiative flux'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('rflx2')
          IF (itype == 0) ncinfo = 'Variance of total radiative flux'
          IF (itype == 1) ncinfo = 'W^2/m^4'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sflx')
          IF (itype == 0) ncinfo = 'Shortwave radiative flux'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sflx2')
          IF (itype == 0) ncinfo = 'Variance of shortwave radiative flux'
          IF (itype == 1) ncinfo = 'W^2/m^4'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sw_up')
          IF (itype == 0) ncinfo = 'Upwelling shortwave radiation'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sw_down')
          IF (itype == 0) ncinfo = 'Downwelling shortwave radiation'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('lw_up')
          IF (itype == 0) ncinfo = 'Upwelling longwave radiation'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('lw_down')
          IF (itype == 0) ncinfo = 'Downwelling longwave radiation'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('l_2')
          IF (itype == 0) ncinfo = 'Variance of liquid water mixing ratio'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('l_3')
          IF (itype == 0) ncinfo = 'Third moment of liquid water mixing ratio'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('tot_lw')
          IF (itype == 0) ncinfo = 'Resolved turbulent flux of liquid water mixing ratio'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('sed_lw')
          IF (itype == 0) ncinfo = 'Sedimentation flux of r_l'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('cs1')
          IF (itype == 0) ncinfo = 'Fraction of cloudy columns (cs1)'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('cnt_cs1')
          IF (itype == 0) ncinfo = 'Number of cloudy columns (cs1)'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('w_cs1')
          IF (itype == 0) ncinfo = 'Conditional average of w over cs1'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('t_cs1 ')
          IF (itype == 0) ncinfo = 'Average of t over cs1'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('tv_cs1')
          IF (itype == 0) ncinfo = 'Conditional average of theta_v over cs1'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('rt_cs1')
          IF (itype == 0) ncinfo = 'Conditional average of rt over cs1'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('rc_cs1')
          IF (itype == 0) ncinfo = 'Conditional average of total condensate over cs1'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('wt_cs1')
          IF (itype == 0) ncinfo = 'Covariance of wtheta_l flux and cs1'
          IF (itype == 1) ncinfo = 'K*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('wtv_cs1')
          IF (itype == 0) ncinfo = 'Covariance of wtheta_v flux and cs1'
          IF (itype == 1) ncinfo = 'K*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('wrt_cs1')
          IF (itype == 0) ncinfo = 'Covariance of wr_t flux and cs1'
          IF (itype == 1) ncinfo = 'kg/kg*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('cs2')
          IF (itype == 0) ncinfo = 'Fraction of cloud core columns (cs2)'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('cnt_cs2')
          IF (itype == 0) ncinfo = 'Number of cloud core columns (cs2)'
          IF (itype == 1) ncinfo = '#'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('w_cs2')
          IF (itype == 0) ncinfo = 'Conditional average of w over cs2'
          IF (itype == 1) ncinfo = 'm/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('t_cs2')
          IF (itype == 0) ncinfo = 'Conditional average of theta_l over cs2'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('tv_cs2')
          IF (itype == 0) ncinfo = 'Conditional average of theta_v over cs2'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('rt_cs2')
          IF (itype == 0) ncinfo = 'Conditional average of rt over cs2'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('rc_cs2')
          IF (itype == 0) ncinfo = 'Conditional average of rl over cs2'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('wt_cs2')
          IF (itype == 0) ncinfo = 'Covariance of wtheta_l flux and cs2'
          IF (itype == 1) ncinfo = 'K*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('wtv_cs2')
          IF (itype == 0) ncinfo = 'Covariance of wtheta_v flux and cs2'
          IF (itype == 1) ncinfo = 'K*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('wrt_cs2')
          IF (itype == 0) ncinfo = 'Covariance of wr_t flux and cs2'
          IF (itype == 1) ncinfo = 'kg/kg*m/s'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('Nc')
          IF (itype == 0) ncinfo = 'Cloud droplet number concentration'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('Nr')
          IF (itype == 0) ncinfo = 'Rain drop number concentration'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
      case('Ni')
          if (itype == 0) ncinfo = 'Ice number concentration'
          if (itype == 1) ncinfo = 'kg^-1'
          if (itype == 2) ncinfo = 'tttt'
       case('Ns')
          if (itype == 0) ncinfo = 'Snow number concentration'
          if (itype == 1) ncinfo = 'kg^-1'
          if (itype == 2) ncinfo = 'tttt'
       CASE('rr')
          IF (itype == 0) ncinfo = 'Rain water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('rrate')
          IF (itype == 0) ncinfo = 'Precipitation Flux (positive downward)'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       case('srate')
          if (itype == 0) ncinfo = 'Snow deposition flux (positive downward)'
          if (itype == 1) ncinfo = 'W/m^2'
          if (itype == 2) ncinfo = 'ttmt'
       CASE('evap')
          IF (itype == 0) ncinfo = 'Net evap of rain-water'
          IF (itype == 1) ncinfo = 's^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('frc_prc')
          IF (itype == 0) ncinfo = 'Conditionally sampled rain fraction'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('prc_prc')
          IF (itype == 0) ncinfo = 'Conditionally sampled precipitation flux'
          IF (itype == 1) ncinfo = 'W/m^2'
          IF (itype == 2) ncinfo = 'ttmt'
       CASE('frc_ran')
          IF (itype == 0) ncinfo = 'Rain water fraction'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('hst_srf')
          IF (itype == 0) ncinfo = 'Histogram of surface rain rates'
          IF (itype == 1) ncinfo = '-'
          IF (itype == 2) ncinfo = 'tttt'
       !
       !
       ! SALSA analysis fields
       CASE('S_RH')
          IF (itype == 0) ncinfo = 'SALSA Relative humidity'
          IF (itype == 1) ncinfo = '1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_RHI')
          IF (itype == 0) ncinfo = 'SALSA Relative humidity over ice'
          IF (itype == 1) ncinfo = '1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Nc')
          IF (itype == 0) ncinfo = 'SALSA cdnc'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Ncba')
          IF (itype == 0) ncinfo = 'SALSA cloud droplet size distribution, regime a'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttcla'
       CASE('S_Ncbb')
          IF (itype == 0) ncinfo = 'SALSA cloud droplet size distribution, regime b'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttclb'
       CASE('S_Rwca')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of cloud droplets, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwcb')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of cloud droplets, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwcba')
          IF (itype == 0) ncinfo = 'SALSA bin cloud droplet radius, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttcla'
       CASE('S_Rwcbb')
          IF (itype == 0) ncinfo = 'SALSA bin cloud droplet radius, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttclb'
       CASE('S_Np')
          IF (itype == 0) ncinfo = 'SALSA rdnc'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Npba')
          IF (itype == 0) ncinfo = 'SALSA precipitation size distribution'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttprc'
       CASE('S_Rwpa')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of precipitation particles'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwpba')
          IF (itype == 0) ncinfo = 'SALSA bin precipitation particle radius'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttprc'

       CASE('S_Ni')
          IF (itype == 0) ncinfo = 'SALSA ice nuclei'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Niba')
          IF (itype == 0) ncinfo = 'SALSA ice size distribution, regime a'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttica'
       CASE('S_Nibb')
          IF (itype == 0) ncinfo = 'SALSA ice size distribution, regime b'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttica'
       CASE('S_Rwia')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of ice particles, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwib')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of ice particles, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwiba')
          IF (itype == 0) ncinfo = 'SALSA bin ice particle radius, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttica'
       CASE('S_Rwibb')
          IF (itype == 0) ncinfo = 'SALSA bin ice particle radius, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttticb'
       CASE('S_Ns')
          IF (itype == 0) ncinfo = 'SALSA sdnc'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Nsba')
          IF (itype == 0) ncinfo = 'SALSA snow size distribution'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttsnw'
       CASE('S_Rwsa')
          IF (itype == 0) ncinfo = 'SALSA number mean radius of snow'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwsba')
          IF (itype == 0) ncinfo = 'SALSA bin snow radius'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttsnw'


       CASE('S_Na')
          IF (itype == 0) ncinfo = 'SALSA total number of soluble aerosols, (regime A)'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Nb')
          IF (itype == 0) ncinfo = 'SALSA total number of insoluble aerosols, (regime B)'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Naba')
          IF (itype == 0) ncinfo = 'Aerosol size distribution, regime A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttaea'
       CASE('S_Nabb')
          IF (itype == 0) ncinfo = 'Aerosol size distribution, regime B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttttaeb'
       CASE('S_Rwaa')
          IF (itype == 0) ncinfo = 'SALSA number mean wet radius of aerosols, regime A'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwab')
          IF (itype == 0) ncinfo = 'SALSA number mean wet radius of aerosols, regime B'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_Rwaba')
          IF (itype == 0) ncinfo = 'SALSA bin aerosol wet radius, regime A'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttaea'
       CASE('S_Rwabb')
          IF (itype == 0) ncinfo = 'SALSA bin aerosol wet radius, regime B'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'ttttaeb'
       CASE('S_Nact')
          IF (itype == 0) ncinfo = 'SALSA Number of newly activated droplets'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aSO4a')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of SO4, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aSO4b')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of SO4, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aNHa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of NH3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aNHb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of NH3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aNOa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of NO3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aNOb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of NO3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aOCa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of OC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aOCb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of OC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aBCa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of BC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aBCb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of BC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aDUa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of DU, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aDUb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of DU, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aSSa')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of SS, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_aSSb')
          IF (itype == 0) ncinfo = 'SALSA aerosol mass concentration of SS, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cSO4a')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of SO4, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cSO4b')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of SO4, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cNHa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of NH3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cNHb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of NH3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cNOa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of NO3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cNOb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of NO3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cOCa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of OC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cOCb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of OC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cBCa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of BC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cBCb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of BC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cDUa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of DU, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cDUb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of DU, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cSSa')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of SS, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_cSSb')
          IF (itype == 0) ncinfo = 'SALSA CCN mass concentration of SS, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iSO4a')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of SO4, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iSO4b')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of SO4, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iNHa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of NH3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iNHb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of NH3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iNOa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of NO3, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iNOb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of NO3, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iOCa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of OC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iOCb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of OC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iBCa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of BC, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iBCb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of BC, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iDUa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of DU, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iDUb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of DU, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iSSa')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of SS, regime A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('S_iSSb')
          IF (itype == 0) ncinfo = 'SALSA IN mass concentration of SS, regime B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'

       !
       !
       ! SALSA profile statistics
       CASE('P_Naa')
          IF (itype == 0) ncinfo = 'SALSA aerosol number concentration in regime A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Nab')
          IF (itype == 0) ncinfo = 'SALSA aerosol number concentration in regime B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Nca')
          IF (itype == 0) ncinfo = 'SALSA CDNC in regime A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Ncb')
          IF (itype == 0) ncinfo = 'SALSA CDNC in regime B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Np')
          IF (itype == 0) ncinfo = 'SALSA rdnc'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Nia')
          IF (itype == 0) ncinfo = 'SALSA ice number concentration in regime A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Nib')
          IF (itype == 0) ncinfo = 'SALSA ice number concentration in regime B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Ns')
          IF (itype == 0) ncinfo = 'SALSA snow number concentration'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwaa')
          IF (itype == 0) ncinfo = 'SALSA mean aerosol wet radius, regime A'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwab')
          IF (itype == 0) ncinfo = 'SALSA mean aerosol wet radius, regime B'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwca')
          IF (itype == 0) ncinfo = 'SALSA mean cloud droplet radius, regime A'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwcb')
          IF (itype == 0) ncinfo = 'SALSA mean cloud droplet radius, regime B'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwp')
          IF (itype == 0) ncinfo = 'SALSA mean drizzle drop radius'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwia')
          IF (itype == 0) ncinfo = 'SALSA mean ice radius, regime a'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Rwib')
          IF (itype == 0) ncinfo = 'SALSA mean ice radius, regime b'
          IF (itype == 1) ncinfo = 'm'
          IF (itype == 2) ncinfo = 'tttt'
       case('P_Rws')
          if (itype == 0) ncinfo = 'SALSA mean snow radius'
          if (itype == 1) ncinfo = 'm'
          if (itype == 2) ncinfo = 'tttt'
       CASE('P_cSO4a')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SO4 in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSO4c')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SO4 in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSO4p')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SO4 in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSO4i')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SO4 in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSO4s')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SO4 in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cOCa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of OC in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cOCc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of OC in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cOCp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of OC in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cOCi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of OC in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cOCs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of OC in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cBCa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of BC in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cBCc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of BC in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cBCp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of BC in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cBCi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of BC in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cBCs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of BC in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cDUa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of DU in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cDUc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of DU in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cDUp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of DU in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cDUi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of DU in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cDUs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of DU in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSSa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SS in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSSc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SS in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSSp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SS in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSSi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SS in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cSSs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of SS in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNHa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NH3 in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNHc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NH3 in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNHp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NH3 in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNHi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NH3 in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNHs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NH3 in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNOa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NO3 in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNOc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NO3 in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNOp')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NO3 in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNOi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NO3 in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cNOs')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of NO3 in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cH2Oa')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of H2O in aerosols'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cH2Oc')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of H2O in cloud droplets'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cH2Op')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of H2O in drizzle drops'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cH2Oi')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of H2O in ice'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cH2Os')
          IF (itype == 0) ncinfo = 'SALSA total mass mixing ratio of H2O in snow'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_rl')
          IF (itype == 0) ncinfo = 'Cloud water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_rr')
          IF (itype == 0) ncinfo = 'Precipitation mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_ri')
          IF (itype == 0) ncinfo = 'Ice water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_rs')
          IF (itype == 0) ncinfo = 'Snow water mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_rv')
          IF (itype == 0) ncinfo = 'Water vapor mixing ratio'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_RH')
          IF (itype == 0) ncinfo = 'Relative humidity'
          IF (itype == 1) ncinfo = '%'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_RHi')
          IF (itype == 0) ncinfo = 'Relative humidity over ice'
          IF (itype == 1) ncinfo = '%'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Na_c')
          IF (itype == 0) ncinfo = 'Aerosol number concentration in cloudy columns'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Nc_c')
          IF (itype == 0) ncinfo = 'Cloud droplet number concentration in cloudy columns'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Np_c')
          IF (itype == 0) ncinfo = 'Rain drop number concentration in cloudy columns'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_cfrac')
          IF (itype == 0) ncinfo = 'Fraction of cloudy columns'
          IF (itype == 1) ncinfo = ''
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_clw_c')
          IF (itype == 0) ncinfo = 'Cloud liquid water in cloudy columns'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_thl_c')
          IF (itype == 0) ncinfo = 'Liquid water potential temperature in cloudy columns'
          IF (itype == 1) ncinfo = 'K'
          IF (itype == 2) ncinfo = 'tttt'
       CASE('P_Naba')
          IF (itype == 0) ncinfo = 'Aerosol number concentration in size bins A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_SO4aa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_OCaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_BCaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_DUaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_SSaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_NHaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_NOaa')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in aerosol bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaea'
       CASE('P_Nabb')
          IF (itype == 0) ncinfo = 'Aerosol number concentration in size bins B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_SO4ab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_OCab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_BCab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_DUab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_SSab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_NHab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_NOab')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in aerosol bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztaeb'
       CASE('P_Ncba')
          IF (itype == 0) ncinfo = 'Cloud droplet number concentration in size bins A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_SO4ca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_OCca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of Oc in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_BCca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_DUca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_SSca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_NHca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_NOca')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in cloud bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztcla'
       CASE('P_Ncbb')
          IF (itype == 0) ncinfo = 'Cloud droplet number concentration in size bins B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_SO4cb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_OCcb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_BCcb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_DUcb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_SScb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_NHcb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_NOcb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in cloud bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztclb'
       CASE('P_Npb')
          IF (itype == 0) ncinfo = 'Number concentration of drizzle'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_SO4pb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_OCpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_BCpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_DUpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_SSpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_NHpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_NOpb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in drizzle bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztprc'
       CASE('P_Niba')
          IF (itype == 0) ncinfo = 'Ice number concentration in size bins A'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_SO4ia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_OCia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of Oc in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_BCia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_DUia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_SSia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_NHia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_NOia')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in ice bins A'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztica'
       CASE('P_Nibb')
          IF (itype == 0) ncinfo = 'Ice number concentration in size bins B'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_SO4ib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_OCib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_BCib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_DUib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_SSib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_NHib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_NOib')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in ice bins B'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttzticb'
       CASE('P_Nsb')
          IF (itype == 0) ncinfo = 'Number concentration of snow'
          IF (itype == 1) ncinfo = 'kg^-1'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_SO4sb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SO4 in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_OCsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of OC in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_BCsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of BC in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_DUsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of DU in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_SSsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of SS in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_NHsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NH3 in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       CASE('P_NOsb')
          IF (itype == 0) ncinfo = 'Mass mixing ratio of NO3 in snow bins'
          IF (itype == 1) ncinfo = 'kg/kg'
          IF (itype == 2) ncinfo = 'ttztsnw'
       ! -----
       CASE DEFAULT
          IF (myid == 0) PRINT *, 'ABORTING: ncinfo: variable not found ',trim(short_name)
          CALL appl_abort(0)
    END SELECT

 END FUNCTION ncinfo

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
 SUBROUTINE close_aero_nc(ncid)
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ncid

    INTEGER :: iret

    iret = nf90_close(ncid)

 END SUBROUTINE close_aero_nc

 END MODULE ncio
