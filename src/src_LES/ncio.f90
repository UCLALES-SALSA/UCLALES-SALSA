module ncio

  use netcdf
  use mpi_interface, only : appl_abort, myid, pecount, wrxid, wryid

  implicit none
  private

  public :: open_nc, define_nc, define_nc_cs, &
            open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_aero_nc

contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Open_NC: Opens a NetCDF File and identifies starting record
  !
  subroutine open_nc (fname, ename, time, npts, ncid, nrec, version, author, info)

    integer, intent(in)             :: npts
    integer, intent(out)            :: ncid
    integer, intent(out)            :: nrec
    real, intent (in)               :: time
    character (len=80), intent (in) :: fname, ename
    CHARACTER(LEN=80) :: version, author, info

    real, allocatable :: xtimes(:)

    character (len=8)  :: date
    character (len=88) :: lfname
    integer :: iret, ncall, VarID, RecordDimID
    logical :: exans

    if (pecount > 1) then
       write(lfname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',wrxid,wryid,'.nc'
    else
       write(lfname,'(a,a3)') trim(fname),'.nc'
    end if

    inquire(file=trim(lfname),exist=exans)

    ncall = 0
    if (.not.exans) then
       call date_and_time(date)
       iret = nf90_create(lfname,NF90_SHARE,ncid)

       iret = nf90_put_att(ncid,NF90_GLOBAL,'title',ename)
       iret = nf90_put_att(ncid,NF90_GLOBAL,'history','Created on '//date)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'Source','UCLALES-SALSA '//trim(version))
       if (len(author)>0) iret = nf90_put_att(ncid, NF90_GLOBAL, 'Author',trim(author)) ! Optional
       if (len(info)>0) iret = nf90_put_att(ncid, NF90_GLOBAL, 'Info',trim(info)) ! Optional
       iret = nf90_put_att(ncid, NF90_GLOBAL, '_FillValue',-999.)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPTS',npts)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'NPROCS',pecount)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'PROCID',myid)
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'IO_version',1.1)
    else
       iret = nf90_open (trim(lfname), NF90_WRITE, ncid)
       iret = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
       iret = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
       ncall=1
       iret = nf90_inq_varid(ncid,'time',VarID)
       allocate (xtimes(nrec+1))
       iret = nf90_get_var(ncid, VarId, xtimes(1:nrec))
       ncall = 1
       do while(ncall <= nrec .and. xtimes(ncall) < time - spacing(1.))
          ncall=ncall+1
       end do
       deallocate(xtimes)
    end if
    nrec = ncall
    iret = nf90_sync(ncid)

  end subroutine open_nc
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Define_NC: Defines the structure of the nc file (if not
  ! already open)
  !
  ! Juha: Added more dimensions to represent bins for aerosol, cloud and
  ! precipitation particles.
  !
  subroutine define_nc(ncID, nRec, nVar, sx, n1, n2, n3, &
                       inae_a,incld_a,inprc,           &
                       inae_b,incld_b,inice_a,inice_b,insnw         )

    integer, intent (in)           :: nVar, ncID
    integer, optional, intent (in) :: n1, n2, n3
    ! Juha: Added
    INTEGER, OPTIONAL, INTENT(in)  :: inae_a,incld_a,inprc, &
                                      inae_b,incld_b,       &
                                      inice_a,inice_b,insnw
    ! --
    integer, intent (inout)        :: nRec
    character (len=7), intent (in) :: sx(nVar)

    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,&
         dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0  ,&
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

    character (len=7) :: xnm
    integer :: iret, n, VarID


    if (nRec == 0) then
       iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
       if (present(n1)) then
          iret = nf90_def_dim(ncID, 'zt', n1, ztID)
          iret = nf90_def_dim(ncID, 'zm', n1, zmID)
       end if
       if (present(n2)) then
          iret = nf90_def_dim(ncID, 'xt', n2, xtID)
          iret = nf90_def_dim(ncID, 'xm', n2, xmID)
       end if
       if (present(n3)) then
          iret = nf90_def_dim(ncID, 'yt', n3, ytID)
          iret = nf90_def_dim(ncID, 'ym', n3, ymID)
       end if
       ! If this is analysis file, dont write binned output by default!
       ! --------------------------------------------------------------
       IF (PRESENT(inae_a)) THEN
          iret = nf90_def_dim(ncID, 'aea', inae_a, aeaID)
       END IF
       IF (PRESENT(inae_b)) THEN
          iret = nf90_def_dim(ncID, 'aeb', inae_b, aebID)
       END IF
       IF (PRESENT(incld_a)) THEN
          iret = nf90_def_dim(ncID, 'cla', incld_a, claID)
       END IF
       IF (PRESENT(incld_b)) THEN
          iret = nf90_def_dim(ncID, 'clb', incld_b, clbID)
       END IF
       IF (PRESENT(inprc)) THEN
          iret = nf90_def_dim(ncID, 'prc', inprc, prcID)
       END IF
       IF (PRESENT(inice_a)) THEN
          iret = nf90_def_dim(ncID, 'ica', inice_a, icaID)
       END IF
       IF (PRESENT(inice_b)) THEN
          iret = nf90_def_dim(ncID, 'icb', inice_b, icbID)
       END IF
       IF (PRESENT(insnw)) THEN
          iret = nf90_def_dim(ncID, 'snow', insnw, snowID)
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

       do n=1,nVar
          select case(trim(ncinfo(2,sx(n))))
          case ('time')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,timeID  ,VarID)
          case ('zt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ztID    ,VarID)
          case ('zm')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,zmID    ,VarID)
          case ('xt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,xtID    ,VarID)
          case ('xm')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,xmID    ,VarID)
          case ('yt')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ytID    ,VarID)
          case ('ym')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,ymID    ,VarID)
          ! Juha: added for size distributions
          case ('aea')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,aeaID   ,VarID)
          case ('aeb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,aebID   ,VarID)
          case ('cla')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,claID   ,VarID)
          case ('clb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,clbID   ,VarID)
          case ('prc')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,prcID   ,VarID)
          !Jaakko added for ice and snow
          case ('ica')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,icaID   ,VarID)
          case ('icb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,icbID   ,VarID)
          case ('snow')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,snowID   ,VarID)
          !Juha added
          case ('ttttaea')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaea,VarID)
          case ('ttttaeb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaeb,VarID)
          case ('ttttcla')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttcla,VarID)
          case ('ttttclb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttclb,VarID)
          case ('ttttprc')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttprc,VarID)
          !Jaakko added
          case ('ttttica')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttica,VarID)
          case ('tttticb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttticb,VarID)
          case ('ttttsnw')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttsnw,VarID)
          ! ---
          case ('tttt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('mttt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('tmtt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tmtt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('ttmt')
             if (present(n2) .and. present(n3)) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttmt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mt,VarID)
             end if
          case ('ttaea')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttaea,VarID)
          case ('ttaeb')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttaeb,VarID)
          case ('ttcla')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttcla,VarID)
          case ('ttclb')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttclb,VarID)
          case ('ttprc')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttprc,VarID)
          case ('ttztaea')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaea,VarID)
          case ('ttztaeb')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaeb,VarID)
          case ('ttztcla')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztcla,VarID)
          case ('ttztclb')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztclb,VarID)
          case ('ttztprc')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztprc,VarID)
          case default
             if (myid == 0) print *, '  ABORTING: NCIO: Bad dimensional information ',trim(ncinfo(2,sx(n)))
             call appl_abort(0)
          end select
          iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,sx(n)))
          iret=nf90_put_att(ncID,VarID,'units'   ,ncinfo(1,sx(n)))
       end do
       iret  = nf90_enddef(ncID)
       iret  = nf90_sync(ncID)
       nRec = 1
    else
       iret = nf90_inquire(ncID, nVariables=n)
       if (n /= nVar) then
          iret = nf90_close(ncID)
          if (myid == 0) print *, '  ABORTING: Incompatible Netcdf File',n,nVar
          call appl_abort(0)
       else
          do n=1,nVar
             xnm=sx(n)
             iret = nf90_inquire_variable(ncID, n, name=xnm)
          end do
          iret = nf90_sync(ncID)
       end if
    end if

  end subroutine define_nc
  !
  ! ----------------------------------------------------------------------
  ! Subroutine define_nc_cs: Defines the structure of a new column statistics nc file
  !
  subroutine define_nc_cs(ncID, nRec, n2, n3, level, rad_level, spec_list, nspec )
    integer, intent (in) :: ncID, n2, n3, level, rad_level, nspec
    integer, intent (inout) :: nRec ! nRec=0 means new files
    CHARACTER(LEN=3), intent (in) :: spec_list(nspec) ! SALSA species (e.g. SO4, Org,...)

    integer, save :: timeID=0, xtID=0, ytID=0
    integer, save :: dim_ttt(3)
    CHARACTER(LEN=7) nam
    integer :: iret, n, VarID, ss, si

    if (nRec == 0) then
       iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
       iret = nf90_def_dim(ncID, 'xt', n2, xtID)
       iret = nf90_def_dim(ncID, 'yt', n3, ytID)

       dim_ttt= (/xtID,ytID,timeID/)

       iret=nf90_def_var(ncID,'time',NF90_FLOAT,timeID  ,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'time'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'time'))

       iret=nf90_def_var(ncID,'xt',NF90_FLOAT,xtID    ,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'xt'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'xt'))

       iret=nf90_def_var(ncID,'yt',NF90_FLOAT,ytID    ,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'yt'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'yt'))

       iret=nf90_def_var(ncID,'lwp',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'lwp'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'lwp'))

       iret=nf90_def_var(ncID,'rwp',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'rwp'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'rwp'))

       iret=nf90_def_var(ncID,'Nc',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Nc'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'Nc'))

       iret=nf90_def_var(ncID,'Nr',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'Nr'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'Nr'))

       iret=nf90_def_var(ncID,'nccnt',NF90_INT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nccnt'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'nccnt'))

       iret=nf90_def_var(ncID,'nrcnt',NF90_INT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'nrcnt'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'nrcnt'))

       iret=nf90_def_var(ncID,'zb',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zb'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'zb'))

       iret=nf90_def_var(ncID,'zc',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zc'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'zc'))

       iret=nf90_def_var(ncID,'zi1',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zi1_bar'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'zi1_bar'))

       iret=nf90_def_var(ncID,'lmax',NF90_FLOAT,dim_ttt,VarID)
       iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'lmax'))
       iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'lmax'))

       ! Can add: maximum/minimum vertical velocities and their variances,
       ! surface heat and humidity fluxes, buoyancy statistics,...

       IF (rad_level==3) THEN
          iret=nf90_def_var(ncID,'albedo',NF90_FLOAT,dim_ttt,VarID)
          iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'albedo'))
          iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'albedo'))
       ENDIF

       IF (level>=4) THEN
          ! Aerosol and water removal
           iret=nf90_def_var(ncID,'rmH2Odr',NF90_FLOAT,dim_ttt,VarID)
           iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'rmH2Odr'))
           iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'rmH2Odr'))

           iret=nf90_def_var(ncID,'rmH2Ocl',NF90_FLOAT,dim_ttt,VarID)
           iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'rmH2Ocl'))
           iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'rmH2Ocl'))

           iret=nf90_def_var(ncID,'rmH2Opr',NF90_FLOAT,dim_ttt,VarID)
           iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'rmH2Opr'))
           iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'rmH2Opr'))

           DO ss = 1,nspec
              nam='rm'//spec_list(ss)//'dr'
              iret=nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
              iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
              iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))

              nam='rm'//spec_list(ss)//'cl'
              iret=nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
              iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
              iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))

              nam='rm'//spec_list(ss)//'pr'
              iret=nf90_def_var(ncID,nam,NF90_FLOAT,dim_ttt,VarID)
              iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,nam))
              iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,nam))
           END DO
       ENDIF

       iret  = nf90_enddef(ncID)
       iret  = nf90_sync(ncID)
       nRec = 1
    end if

  end subroutine define_nc_cs
  !
  ! ----------------------------------------------------------------------
  ! Subroutine nc_info: Gets long_name, units and dimension info given a
  ! short name.
  !
  character (len=80) function ncinfo(itype,short_name)

    character (len=40) :: v_lnm ='scalar xx mixing ratio                  '

    integer, intent (in) :: itype
    character (len=*), intent (in) :: short_name

    integer :: scalar_number

    select case (trim(short_name))
    case ('sxx')
       read (short_name(2:3),'(i2.2)') scalar_number
       write(v_lnm(8:9),'(i2.2)') scalar_number
       if (itype==0) ncinfo = v_lnm
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('time')
       if (itype==0) ncinfo = 'Time'
       if (itype==1) ncinfo = 's'
       if (itype==2) ncinfo = 'time'
    case('zt')
       if (itype==0) ncinfo = 'Vertical displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'zt'
    case('zm')
       if (itype==0) ncinfo = 'Vertical displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'zm'
    case('xt')
       if (itype==0) ncinfo = 'East-west displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'xt'
    case('xm')
       if (itype==0) ncinfo = 'East-west displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'xm'
    case('yt')
       if (itype==0) ncinfo = 'North-south displacement of cell centers'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'yt'
    case('ym')
       if (itype==0) ncinfo = 'North-south displacement of cell edges'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ym'
    ! Juha: added for SALSA
    case('aea')
       if (itype==0) ncinfo = 'Aerosol size bins, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'aea'
    case('aeb')
       if (itype==0) ncinfo = 'Aerosol size bins, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'aeb'
    case('cla')
       if (itype==0) ncinfo = 'Cloud droplet size bins, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'cla'
    case('clb')
       if (itype==0) ncinfo = 'Cloud droplet size bins, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'clb'
    case('prc')
       if (itype==0) ncinfo = 'Precipitation size bins'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'prc'
    case('ica')
       if (itype==0) ncinfo = 'Ice cloud droplet size bins, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ica'
    case('icb')
       if (itype==0) ncinfo = 'Ice cloud droplet size bins, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'icb'
    case('snw')
       if (itype==0) ncinfo = 'Snow size bins'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'snow'
    !----
    case('u0')
       if (itype==0) ncinfo = 'Geostrophic zonal wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'zt'
    case('v0')
       if (itype==0) ncinfo = 'Geostrophic meridional wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'zt'
    case('dn0')
       if (itype==0) ncinfo = 'Base-state density'
       if (itype==1) ncinfo = 'kg/m^3'
       if (itype==2) ncinfo = 'zt'
    case('u')
       if (itype==0) ncinfo = 'Zonal wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'mttt'
    case('v')
       if (itype==0) ncinfo = 'Meridional wind'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'tmtt'
    case('w')
       if (itype==0) ncinfo = 'Vertical velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('t')
       if (itype==0) ncinfo = 'Potential temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('p')
       if (itype==0) ncinfo = 'Pressure'
       if (itype==1) ncinfo = 'Pa'
       if (itype==2) ncinfo = 'tttt'
    case('q')
       if (itype==0) ncinfo = 'Total water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('l')
       if (itype==0) ncinfo = 'Liquid water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('r')
       if (itype==0) ncinfo = 'Rain-water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('f')
       if (itype==0) ncinfo = 'Total ice mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('i')
       if (itype==0) ncinfo = 'Ice mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('s')
       if (itype==0) ncinfo = 'Snow mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('n')
       if (itype==0) ncinfo = 'Rain-drop number mixing ratio'
       if (itype==1) ncinfo = '#/kg'
       if (itype==2) ncinfo = 'tttt'
    case('stke')
       if (itype==0) ncinfo = 'Sub-filter scale TKE'
       if (itype==1) ncinfo = 'J/kg'
       if (itype==2) ncinfo = 'mttt'
    case('cfl')
       if (itype==0) ncinfo = 'Courant number'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('maxdiv')
       if (itype==0) ncinfo = 'Maximum divergence'
       if (itype==1) ncinfo = '1/s'
       if (itype==2) ncinfo = 'time'
    case('zi1_bar')
       if (itype==0) ncinfo = 'Height of maximum theta gradient'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zi2_bar')
       if (itype==0) ncinfo = 'Height of maximum theta variance'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zi3_bar')
       if (itype==0) ncinfo = 'Height of minimum buoyancy flux'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('vtke')
       if (itype==0) ncinfo = 'Vertical integral of total TKE'
       if (itype==1) ncinfo = 'kg/s'
       if (itype==2) ncinfo = 'time'
    case('sfcbflx')
       if (itype==0) ncinfo = 'Surface Buoyancy Flux'
       if (itype==1) ncinfo = 'm/s^2'
       if (itype==2) ncinfo = 'time'
    case('wmax')
       if (itype==0) ncinfo = 'Maximum vertical velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('tsrf')
       if (itype==0) ncinfo = 'Surface temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'time'
    case('ustar')
       if (itype==0) ncinfo = 'Surface friction velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('shf_bar')
       if (itype==0) ncinfo = 'Sensible heat flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('lhf_bar')
       if (itype==0) ncinfo = 'Latent heat flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('zi_bar')
       if (itype==0) ncinfo = 'Height of maximum scalar gradient'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('lwp_bar','lwp')
       if (itype==0) ncinfo = 'Liquid-water path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('lwp_var')
       if (itype==0) ncinfo = 'Liquid-water path variance'
       if (itype==1) ncinfo = 'kg^2/m^4'
       if (itype==2) ncinfo = 'time'
    case('zc')
       if (itype==0) ncinfo = 'Cloud-top height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zb')
       if (itype==0) ncinfo = 'Cloud-base height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('cfrac')
       if (itype==0) ncinfo = 'Cloud fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('lmax')
       if (itype==0) ncinfo = 'Maximum liquid water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('albedo')
       if (itype==0) ncinfo = 'Reflected (TOA) shortwave radiation'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('rwp_bar','rwp')
       if (itype==0) ncinfo = 'Rain-water path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('prcp')
       if (itype==0) ncinfo = 'Surface precipitation rate'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('pfrac')
       if (itype==0) ncinfo = 'Surface precipitation fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'time'
    case('CCN')
       if (itype==0) ncinfo = 'Cloud condensation nuclei'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nrain')
       if (itype==0) ncinfo = 'Conditionally sampled rain number mixing ratio'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nrcnt')
       if (itype==0) ncinfo = 'Rain cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('nccnt')
       if (itype==0) ncinfo = 'Cloud cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    !
    !
    ! SALSA temporal statistics
    case('Nc_ic')
       if (itype==0) ncinfo = 'In-cloud CDNC'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('Na_oc')
       if (itype==0) ncinfo = 'Aerosol number concentration outside clouds'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('Na_int')
       if (itype==0) ncinfo = 'In-cloud intersitial aerosol number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('SO4_ic')
       if (itype==0) ncinfo = 'Cloud droplet SO4 mass mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('SO4_oc')
       if (itype==0) ncinfo = 'Aerosol SO4 mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('SO4_int')
       if (itype==0) ncinfo = 'SO4 mass mixing ratio in intersitial aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('OC_ic')
       if (itype==0) ncinfo = 'Cloud droplet OC mass mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('OC_oc')
       if (itype==0) ncinfo = 'Aerosol OC mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('OC_int')
       if (itype==0) ncinfo = 'OC mass mixing ratio in interstitial aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('BC_ic')
       if (itype==0) ncinfo = 'Cloud droplet BC mass mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('BC_oc')
       if (itype==0) ncinfo = 'Aerosol BC mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('BC_int')
       if (itype==0) ncinfo = 'BC mass mixing ratio in interstitial aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('DU_ic')
       if (itype==0) ncinfo = 'Cloud droplet DU mass mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('DU_oc')
       if (itype==0) ncinfo = 'Aerosol DU mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('DU_int')
       if (itype==0) ncinfo = 'DU mass mixing ration in interstitial aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('SS_ic')
       if (itype==0) ncinfo = 'Cloud droplet mass mixing ratio of SS'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('SS_oc')
       if (itype==0) ncinfo = 'Aerosol SS mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('SS_int')
       if (itype==0) ncinfo = 'SS mass mixing ratio in interstitial particles'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NH3_ic')
       if (itype==0) ncinfo = 'Cloud droplet mass mixing ratio of NH3'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NH3_oc')
       if (itype==0) ncinfo = 'Aerosol NH3 mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NH3_int')
       if (itype==0) ncinfo = 'NH3 mass mixing ratio in interstitial particles'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NO3_ic')
       if (itype==0) ncinfo = 'Cloud droplet mass mixing ratio of NO3'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NO3_oc')
       if (itype==0) ncinfo = 'Aerosol NO3 mass mixing ratio outside clouds'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('NO3_int')
       if (itype==0) ncinfo = 'NO3 mass mixing ratio in interstitial particles'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('rmH2Oae','rmH2Odr')
       if (itype==0) ncinfo = 'Deposition of H2O with aerosols'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmH2Ocl')
       if (itype==0) ncinfo = 'Deposition of H2O with cloud droplets'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmH2Opr')
       if (itype==0) ncinfo = 'Deposition of water with rain'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSO4dr')
       if (itype==0) ncinfo = 'Aerosol deposition of SO4'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSO4cl')
       if (itype==0) ncinfo = 'Cloud deposition of SO4'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSO4pr')
       if (itype==0) ncinfo = 'Precipitation deposition of SO4'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSO4wt')
       if (itype==0) ncinfo = 'Total wet deposition of SO4'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSO4tt')
       if (itype==0) ncinfo = 'Total deposition of SO4'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmOCdr')
       if (itype==0) ncinfo = 'Aerosol deposition of OC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmOCcl')
       if (itype==0) ncinfo = 'Cloud deposition of OC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmOCpr')
       if (itype==0) ncinfo = 'Precipitation deposition of OC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmOCwt')
       if (itype==0) ncinfo = 'Total wet deposition of OC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmOCtt')
       if (itype==0) ncinfo = 'Total deposition of OC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmBCdr')
       if (itype==0) ncinfo = 'Aerosol deposition of BC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmBCcl')
       if (itype==0) ncinfo = 'Cloud deposition of BC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmBCpr')
       if (itype==0) ncinfo = 'Precipitation deposition of BC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmBCwt')
       if (itype==0) ncinfo = 'Total wet deposition of BC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmBCtt')
       if (itype==0) ncinfo = 'Total deposition of BC'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmDUdr')
       if (itype==0) ncinfo = 'Aerosol deposition of DU'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmDUcl')
       if (itype==0) ncinfo = 'Cloud deposition of DU'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmDUpr')
       if (itype==0) ncinfo = 'Precipitation deposition of DU'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmDUwt')
       if (itype==0) ncinfo = 'Total wet deposition of DU'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmDUtt')
       if (itype==0) ncinfo = 'Total deposition of DU'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSSdr')
       if (itype==0) ncinfo = 'Aerosol deposition of SS'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSScl')
       if (itype==0) ncinfo = 'Cloud deposition of SS'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSSpr')
       if (itype==0) ncinfo = 'Precipitation deposition of SS'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSSwt')
       if (itype==0) ncinfo = 'Total wet deposition of SS'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmSStt')
       if (itype==0) ncinfo = 'Total deposition of SS'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNH3dr')
       if (itype==0) ncinfo = 'Aerosol deposition of NH3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNH3cl')
       if (itype==0) ncinfo = 'Cloud deposition of NH3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNH3pr')
       if (itype==0) ncinfo = 'Precipitation deposition of NH3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNH3wt')
       if (itype==0) ncinfo = 'Total wet deposition of NH3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNH3tt')
       if (itype==0) ncinfo = 'Total deposition of NH3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNO3dr')
       if (itype==0) ncinfo = 'Aerosol deposition of NO3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNO3cl')
       if (itype==0) ncinfo = 'Cloud deposition of NO3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNO3pr')
       if (itype==0) ncinfo = 'Precipitation deposition of NO3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNO3wt')
       if (itype==0) ncinfo = 'Total wet deposition of NO3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('rmNO3tt')
       if (itype==0) ncinfo = 'Total deposition of NO3'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    ! // SALSA temporal
    case('fsttm')
       if (itype==0) ncinfo = 'First sample time'
       if (itype==1) ncinfo = 's'
       if (itype==2) ncinfo = 'time'
    case('lsttm')
       if (itype==0) ncinfo = 'Last sample time'
       if (itype==1) ncinfo = 's'
       if (itype==2) ncinfo = 'time'
    case('nsmp')
       if (itype==0) ncinfo = 'Number of samples'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('u_2')
       if (itype==0) ncinfo = 'Variance of u wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'tttt'
    case('v_2')
       if (itype==0) ncinfo = 'Variance of v wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'tttt'
    case('w_2')
       if (itype==0) ncinfo = 'Second raw moment of w wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('t_2')
       if (itype==0) ncinfo = 'Variance of theta'
       if (itype==1) ncinfo = 'K^2'
       if (itype==2) ncinfo = 'tttt'
    case('w_3')
       if (itype==0) ncinfo = 'Third raw moment of w wind'
       if (itype==1) ncinfo = 'm^3/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('t_3')
       if (itype==0) ncinfo = 'Third moment of theta'
       if (itype==1) ncinfo = 'K^3'
       if (itype==2) ncinfo = 'tttt'
    case('tot_tw')
       if (itype==0) ncinfo = 'Total vertical flux of theta'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_tw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of theta'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_uw')
       if (itype==0) ncinfo = 'Total vertical flux of u-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_uw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of u-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_vw')
       if (itype==0) ncinfo = 'Total vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_vw')
       if (itype==0) ncinfo = 'SGS vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_ww')
       if (itype==0) ncinfo = 'Total vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_ww')
       if (itype==0) ncinfo = 'SGS vertical flux of w-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('km')
       if (itype==0) ncinfo = 'Eddy viscosity'
       if (itype==1) ncinfo = 'm^2/s'
       if (itype==2) ncinfo = 'ttmt'
    case('kh')
       if (itype==0) ncinfo = 'Eddy diffusivity'
       if (itype==1) ncinfo = 'm^2/s'
       if (itype==2) ncinfo = 'ttmt'
    case('lmbd')
       if (itype==0) ncinfo = 'Mixing lengthscale'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('lmbde')
       if (itype==0) ncinfo = 'Dissipation lengthscale'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_tke')
       if (itype==0) ncinfo = 'Sub-filter scale TKE'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_boy')
       if (itype==0) ncinfo = 'Subfilter Buoyancy production of TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_shr')
       if (itype==0) ncinfo = 'Shear production of SGS TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('boy_prd')
       if (itype==0) ncinfo = 'Buoyancy production of resolved TKE '
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('shr_prd')
       if (itype==0) ncinfo = 'Shear production of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('trans')
       if (itype==0) ncinfo = 'Net transport of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('diss')
       if (itype==0) ncinfo = 'Dissipation rate of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('dff_u')
       if (itype==0) ncinfo = 'u(du/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('dff_v')
       if (itype==0) ncinfo = 'v(dv/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('dff_w')
       if (itype==0) ncinfo = 'w(dw/dt) from diffusion'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('adv_u')
       if (itype==0) ncinfo = 'u(du/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('adv_v')
       if (itype==0) ncinfo = 'v(dv/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('adv_w')
       if (itype==0) ncinfo = 'w(dw/dt) from advection'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('prs_u')
       if (itype==0) ncinfo = 'u(du/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('prs_v')
       if (itype==0) ncinfo = 'v(dv/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('prs_w')
       if (itype==0) ncinfo = 'w(dw/dt) from pressure'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('prd_uw')
       if (itype==0) ncinfo = 'uw shear production'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('storage')
       if (itype==0) ncinfo = 'Rate of increase of resolved TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'tttt'
    case('q_2')
       if (itype==0) ncinfo = 'Variance of total water'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('q_3')
       if (itype==0) ncinfo = 'Third moment of total water'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('tot_qw')
       if (itype==0) ncinfo = 'Total vertical flux of q'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_qw')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of q'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('rflx')
       if (itype==0) ncinfo =  'Total Radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('rflx2')
       if (itype==0) ncinfo = 'Variance of total radiative flux'
       if (itype==1) ncinfo = 'W^2/m^4'
       if (itype==2) ncinfo = 'ttmt'
    case('sflx')
       if (itype==0) ncinfo = 'Shortwave radiative flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sflx2')
       if (itype==0) ncinfo = 'Variance of shortwave radiative flux'
       if (itype==1) ncinfo = 'W^2/m^4'
       if (itype==2) ncinfo = 'ttmt'
    case('l_2')
       if (itype==0) ncinfo = 'Variance of liquid water mixing ratio'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('l_3')
       if (itype==0) ncinfo = 'Third moment of liquid water mixing ratio'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('tot_lw')
       if (itype==0) ncinfo = 'Resolved turbulent flux of liquid water mixing ratio'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sed_lw')
       if (itype==0) ncinfo = 'Sedimentation flux of r_l'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('cs1')
       if (itype==0) ncinfo = 'Fraction of cloudy columns (cs1)'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('cnt_cs1')
       if (itype==0) ncinfo = 'Number of cloudy columns (cs1)'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'tttt'
    case('w_cs1')
       if (itype==0) ncinfo = 'Conditional average of w over cs1'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs1')
       if (itype==0) ncinfo = 'Conditional average of theta_l over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'ttmt'
    case('tv_cs1')
       if (itype==0) ncinfo = 'Conditional average of theta_v over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs1')
       if (itype==0) ncinfo = 'Conditional average of rt over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rl_cs1')
       if (itype==0) ncinfo = 'Conditional average of rl over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wt_cs1')
       if (itype==0) ncinfo = 'Covariance of wtheta_l flux and cs1'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wv_cs1')
       if (itype==0) ncinfo = 'Covariance of wtheta_v flux and cs1'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wr_cs1')
       if (itype==0) ncinfo = 'Covariance of wr_t flux and cs1'
       if (itype==1) ncinfo = 'kg/kg*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('cs2')
       if (itype==0) ncinfo = 'Fraction of cloud core columns (cs2)'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('cnt_cs2')
       if (itype==0) ncinfo = 'Number of cloud core columns (cs2)'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'tttt'
    case('w_cs2')
       if (itype==0) ncinfo = 'Conditional average of w over cs2'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs2')
       if (itype==0) ncinfo = 'Conditional average of theta_l over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('tv_cs2')
       if (itype==0) ncinfo = 'Conditional average of theta_v over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs2')
       if (itype==0) ncinfo = 'Conditional average of rt over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rl_cs2')
       if (itype==0) ncinfo = 'Conditional average of rl over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wt_cs2')
       if (itype==0) ncinfo = 'Covariance of wtheta_l flux and cs2'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wv_cs2')
       if (itype==0) ncinfo = 'Covariance of wtheta_v flux and cs2'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wr_cs2')
       if (itype==0) ncinfo = 'Covariance of wr_t flux and cs2'
       if (itype==1) ncinfo = 'kg/kg*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('Nc')
       if (itype==0) ncinfo = 'Cloud droplet number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Nr')
       if (itype==0) ncinfo = 'Rain drop number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('rr')
       if (itype==0) ncinfo = 'Rain water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('precip')
       if (itype==0) ncinfo = 'Precipitation Flux (positive downward'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('evap')
       if (itype==0) ncinfo = 'Net evap  of rain-water'
       if (itype==1) ncinfo = 's^-1'
       if (itype==2) ncinfo = 'tttt'
    case('frc_prc')
       if (itype==0) ncinfo = 'Conditionally sampled rain fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'ttmt'
    case('prc_prc')
       if (itype==0) ncinfo = 'Conditionally sampled rain rate'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('frc_ran')
       if (itype==0) ncinfo = 'Rain water fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('hst_srf')
       if (itype==0) ncinfo = 'Histogram of surface rain rates'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    !
    !
    ! SALSA analysis fields
    case('S_RH')
       if (itype==0) ncinfo = 'SALSA Relative humidity'
       if (itype==1) ncinfo = '1'
       if (itype==2) ncinfo = 'tttt'
    case('S_RHI')
       if (itype==0) ncinfo = 'SALSA Relative humidity over ice'
       if (itype==1) ncinfo = '1'
       if (itype==2) ncinfo = 'tttt'
    case('S_Nc')
       if (itype==0) ncinfo = 'SALSA cdnc'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('S_Ncba')
       if (itype==0) ncinfo = 'SALSA cloud droplet size distribution, regime a'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttttcla'
    case('S_Ncbb')
       if (itype==0) ncinfo = 'SALSA cloud droplet size distribution, regime b'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttttclb'
    case('S_Rwca')
       if (itype==0) ncinfo = 'SALSA number mean radius of cloud droplets, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwcb')
       if (itype==0) ncinfo = 'SALSA number mean radius of cloud droplets, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwcba')
       if (itype==0) ncinfo = 'SALSA bin cloud droplet radius, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttcla'
    case('S_Rwcbb')
       if (itype==0) ncinfo = 'SALSA bin cloud droplet radius, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttclb'
    case('S_Np')
       if (itype==0) ncinfo = 'SALSA rdnc'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('S_Npba')
       if (itype==0) ncinfo = 'SALSA precipitation size distribution'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttttprc'
    case('S_Rwpa')
       if (itype==0) ncinfo = 'SALSA number mean radius of precipitation particles'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwpba')
       if (itype==0) ncinfo = 'SALSA bin precipitation particle radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttprc'

    case('S_Nic')
       if (itype==0) ncinfo = 'SALSA ice nuclei'
       if (itype==1) ncinfo = 'm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('S_Niba')
       if (itype==0) ncinfo = 'SALSA ice particle size distribution, regime a'
       if (itype==1) ncinfo = 'm^-3'
       if (itype==2) ncinfo = 'ttttica'
    case('S_Nibb')
       if (itype==0) ncinfo = 'SALSA ice particle size distribution, regime b'
       if (itype==1) ncinfo = 'm^-3'
       if (itype==2) ncinfo = 'ttttica'
    case('S_Rwia')
       if (itype==0) ncinfo = 'SALSA number mean radius of ice particles, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwib')
       if (itype==0) ncinfo = 'SALSA number mean radius of ice particles, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwiba')
       if (itype==0) ncinfo = 'SALSA bin ice particle radius, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttcla'
    case('S_Rwibb')
       if (itype==0) ncinfo = 'SALSA bin ice particle radius, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttcla'
    case('S_Ns')
       if (itype==0) ncinfo = 'SALSA sdnc'
       if (itype==1) ncinfo = 'm^-3'
       if (itype==2) ncinfo = 'tttt'
    case('S_Nsba')
       if (itype==0) ncinfo = 'SALSA snow size distribution'
       if (itype==1) ncinfo = 'm^-3'
       if (itype==2) ncinfo = 'ttttsnw'
    case('S_Rwsa')
       if (itype==0) ncinfo = 'SALSA number mean radius of snow particles'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwsba')
       if (itype==0) ncinfo = 'SALSA bin snow particle radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttsnw'


    case('S_Na')
       if (itype==0) ncinfo = 'SALSA total number of soluble aerosols, (regime A)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('S_Nb')
       if (itype==0) ncinfo = 'SALSA total number of insoluble aerosols, (regime B)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('S_Naba')
       if (itype==0) ncinfo = 'Aerosol size distribution, regime A'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttttaea'
    case('S_Nabb')
       if (itype==0) ncinfo = 'Aerosol size distribution, regime B'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttttaeb'
    case('S_Rwaa')
       if (itype==0) ncinfo = 'SALSA number mean wet radius of aerosols, regime A'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwab')
       if (itype==0) ncinfo = 'SALSA number mean wet radius of aerosols, regime B'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('S_Rwaba')
       if (itype==0) ncinfo = 'SALSA bin aerosol wet radius, regime A'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttaea'
    case('S_Rwabb')
       if (itype==0) ncinfo = 'SALSA bin aerosol wet radius, regime B'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'ttttaeb'
    case('S_Nact')
       if (itype==0) ncinfo = 'SALSA Number of newly activated droplets'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('S_aSO4a')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of SO4, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aSO4b')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of SO4, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aNH3a')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of NH3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aNH3b')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of NH3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aNO3a')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of NO3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aNO3b')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of NO3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aOCa')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of OC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aOCb')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of OC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aBCa')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of BC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aBCb')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of BC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aDUa')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of DU, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aDUb')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of DU, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aSSa')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of SS, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_aSSb')
       if (itype==0) ncinfo = 'SALSA aerosol mass concentration of SS, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cSO4a')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of SO4, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cSO4b')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of SO4, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cNH3a')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of NH3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cNH3b')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of NH3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cNO3a')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of NO3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cNO3b')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of NO3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cOCa')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of OC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cOCb')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of OC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cBCa')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of BC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cBCb')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of BC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cDUa')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of DU, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cDUb')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of DU, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cSSa')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of SS, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_cSSb')
       if (itype==0) ncinfo = 'SALSA CCN mass concentration of SS, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iSO4a')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of SO4, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iSO4b')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of SO4, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iNH3a')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of NH3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iNH3b')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of NH3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iNO3a')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of NO3, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iNO3b')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of NO3, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iOCa')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of OC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iOCb')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of OC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iBCa')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of BC, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iBCb')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of BC, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iDUa')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of DU, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iDUb')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of DU, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iSSa')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of SS, regime A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('S_iSSb')
       if (itype==0) ncinfo = 'SALSA IN mass concentration of SS, regime B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'

    !
    !
    ! SALSA profile statistics
    case('P_Naa')
       if (itype==0) ncinfo = 'SALSA aerosol number concentration in regime A'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('P_Nab')
       if (itype==0) ncinfo = 'SALSA aerosol number concentration in regime B'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('P_Nca')
       if (itype==0) ncinfo = 'SALSA CDNC in regime A'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('P_Ncb')
       if (itype==0) ncinfo = 'SALSA CDNC in regime B'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('P_Np')
       if (itype==0) ncinfo = 'SALSA rdnc'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwaa')
       if (itype==0) ncinfo = 'SALSA mean aerosol wet radius, regime A'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwab')
       if (itype==0) ncinfo = 'SALSA mean aerosol wet radius, regime B'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwca')
       if (itype==0) ncinfo = 'SALSA mean cloud droplet radius, regime A'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwcb')
       if (itype==0) ncinfo = 'SALSA mean cloud droplet radius, regime B'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwp')
       if (itype==0) ncinfo = 'SALSA mean drizzle drop radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwia')
       if (itype==0) ncinfo = 'SALSA mean ice particle radius, regime a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_Rwib')
       if (itype==0) ncinfo = 'SALSA mean ice particle radius, regime b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSO4a')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SO4 in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSO4c')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SO4 in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSO4p')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SO4 in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cOCa')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of OC in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cOCc')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of OC in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cOCp')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of OC in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cBCa')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of BC in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cBCc')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of BC in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cBCp')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of BC in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cDUa')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of DU in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cDUc')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of DU in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cDUp')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of DU in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSSa')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SS in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSSc')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SS in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cSSp')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of SS in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNH3a')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NH3 in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNH3c')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NH3 in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNH3p')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NH3 in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNO3a')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NO3 in aerosols'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNO3c')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NO3 in drizzle drops'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_cNO3p')
       if (itype==0) ncinfo = 'SALSA total mass mixing ratio of NO3 in cloud droplets'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_rl')
       if (itype==0) ncinfo = 'Level 4 cloud water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_rr')
       if (itype==0) ncinfo = 'Level 4 precipitation mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_rv')
       if (itype==0) ncinfo = 'Level 4 water vapor mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('P_RH')
       if (itype==0) ncinfo = 'Level 4 relative humidity'
       if (itype==1) ncinfo = '%'
       if (itype==2) ncinfo = 'tttt'
    case('P_Naba')
       if (itype==0) ncinfo = 'Aerosol number concentration in size bins A'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_SO4aa')
       if (itype==0) ncinfo = 'Mass mixing ratio of SO4 in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_OCaa')
       if (itype==0) ncinfo = 'Mass mixing ratio of OC in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_BCaa')
       if (itype==0) ncinfo = 'Mass mixing ratio of BC in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_DUaa')
       if (itype==0) ncinfo = 'Mass mixing ratio of DU in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_SSaa')
       if (itype==0) ncinfo = 'Mass mixing ratio of SS in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_NH3aa')
       if (itype==0) ncinfo = 'Mass mixing ratio of NH3 in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_NO3aa')
       if (itype==0) ncinfo = 'Mass mixing ratio of NO3 in aerosol bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaea'
    case('P_Nabb')
       if (itype==0) ncinfo = 'Aerosol number concentration in size bins B'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_SO4ab')
       if (itype==0) ncinfo = 'Mass mixing ratio of SO4 in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_OCab')
       if (itype==0) ncinfo = 'Mass mixing ratio of OC in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_BCab')
       if (itype==0) ncinfo = 'Mass mixing ratio of BC in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_DUab')
       if (itype==0) ncinfo = 'Mass mixing ratio of DU in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_SSab')
       if (itype==0) ncinfo = 'Mass mixing ratio of SS in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_NH3ab')
       if (itype==0) ncinfo = 'Mass mixing ratio of NH3 in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_NO3ab')
       if (itype==0) ncinfo = 'Mass mixing ratio of NO3 in aerosol bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztaeb'
    case('P_Ncba')
       if (itype==0) ncinfo = 'Cloud droplet number concentration in size bins A'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_SO4ca')
       if (itype==0) ncinfo = 'Mass mixing ratio of SO4 in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_OCca')
       if (itype==0) ncinfo = 'MAss mixing ratio of Oc in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_BCca')
       if (itype==0) ncinfo = 'Mass mixing ratio of BC in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_DUca')
       if (itype==0) ncinfo = 'Mass mixing ratio of DU in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_SSca')
       if (itype==0) ncinfo = 'Mass mixing ratio of SS in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_NH3ca')
       if (itype==0) ncinfo = 'Mass mixing ratio of NH3 in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_NO3ca')
       if (itype==0) ncinfo = 'Mass mixing ratio of NO3 in cloud bins A'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztcla'
    case('P_Ncbb')
       if (itype==0) ncinfo = 'Cloud droplet number concentration in size bins B'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_SO4cb')
       if (itype==0) ncinfo = 'Mass mixing ratio of SO4 in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_OCcb')
       if (itype==0) ncinfo = 'Mass mixing ratio of OC in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_BCcb')
       if (itype==0) ncinfo = 'Mass mixing ratio of BC in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_DUcb')
       if (itype==0) ncinfo = 'Mass mixing ratio of DU in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_SScb')
       if (itype==0) ncinfo = 'Mass mixing ratio of SS in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_NH3cb')
       if (itype==0) ncinfo = 'Mass mixing ratio of NH3 in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_NO3cb')
       if (itype==0) ncinfo = 'Mass mixing ratio of NO3 in cloud bins B'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztclb'
    case('P_Npb')
       if (itype==0) ncinfo = 'Number concentration of drizzle particles'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_SO4pb')
       if (itype==0) ncinfo = 'Mass mixing ratio of SO4 in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_OCpb')
       if (itype==0) ncinfo = 'Mass mixing ratio of OC in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_BCpb')
       if (itype==0) ncinfo = 'Mass mixing ratio of BC in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_DUpb')
       if (itype==0) ncinfo = 'Mass mixing ratio of DU in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_SSpb')
       if (itype==0) ncinfo = 'Mass mixing ratio of SS in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_NH3pb')
       if (itype==0) ncinfo = 'Mass mixing ratio of NH3 in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    case('P_NO3pb')
       if (itype==0) ncinfo = 'Mass mixing ratio of NO3 in drizzle bins'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'ttztprc'
    ! -----
    case default
       if (myid==0) print *, 'ABORTING: ncinfo: variable not found ',trim(short_name)
       call appl_abort(0)
    end select

  end function ncinfo

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

end module ncio
