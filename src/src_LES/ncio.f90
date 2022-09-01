module ncio

  use netcdf
  use mpi_interface, only : appl_abort, myid, pecount, wrxid, wryid

  implicit none
  private

  public :: open_nc, define_nc, ncinfo, &
            open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_aero_nc

contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine Open_NC: Opens a NetCDF File and identifies starting record
  !
  subroutine open_nc (fname, ename, time, npts, ncid, nrec, version, author, info, par)

    integer, intent(in)             :: npts
    integer, intent(out)            :: ncid
    integer, intent(out)            :: nrec
    real, intent (in)               :: time
    character (len=80), intent (in) :: fname, ename
    CHARACTER(LEN=80) :: version, author, info
    logical, optional, intent(in) :: par

    real, allocatable :: xtimes(:)

    character (len=8)  :: date
    character (len=88) :: lfname
    integer :: iret, ncall, VarID, RecordDimID
    logical :: exans, parallel

    ! Produce output from each PU (column statistics and analysis files when parallel jobs)
    parallel = pecount > 1
    IF (PRESENT(par)) parallel = par .AND. pecount > 1

    if (parallel) then
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
       iret = nf90_put_att(ncid, NF90_GLOBAL, 'IO_version',1.2)
    else
       iret = nf90_open (trim(lfname), NF90_WRITE, ncid)
       iret = nf90_inquire(ncid, unlimitedDimId = RecordDimID)
       iret = nf90_inquire_dimension(ncid, RecordDimID, len=nrec)
       iret = nf90_inq_varid(ncid,'time',VarID)
       allocate (xtimes(nrec+1))
       iret = nf90_get_var(ncid, VarId, xtimes(1:nrec))
       ncall = 1
       do while(ncall <= nrec .and. xtimes(ncall) < time + 0.01)
          ncall=ncall+1
       end do
       deallocate(xtimes)
       IF (time<0.01) ncall=1 ! If time about 0 s, then start from record 1
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
                       n1a, n2a, n2b, nprc, nsnw, nchist, nihist)

    integer, intent (in)           :: nVar, ncID
    integer, intent (inout)        :: nRec
    character (len=7), intent (in) :: sx(nVar)
    integer, optional, intent (in) :: n1, n2, n3,            &
                                      n1a,n2a,n2b,nprc,nsnw, &
                                      nchist,nihist

    integer, save :: timeID=0, ztID=0, zmID=0, xtID=0, xmID=0, ytID=0, ymID=0,&
         dim_mttt(4) = 0, dim_tmtt(4) = 0, dim_ttmt(4) = 0, dim_tttt(4) = 0  ,&
         dim_tt(2)  = 0, dim_mt(2)  = 0

    ! Juha: added
    INTEGER, SAVE :: aeaID=0, aebID=0, prcID=0, snowID=0, hcID=0, hiID=0,                &
         dim_ttttaea(5) = 0, dim_ttttaeb(5) = 0, dim_ttttprc(5) = 0, dim_ttttsnw(5) = 0, & ! z, x, y, bin, time
         dim_ttztaea(3) = 0, dim_ttztaeb(3) = 0, dim_ttztprc(3) = 0, dim_ttztsnw(3) = 0, & ! z, bin, time
         dim_tttzhct(3) = 0, dim_tttzhit(3) = 0, & ! z, histogram_bins, time
         dim_ttt(3) = 0    ! x, y, time

    character (len=7) :: xnm
    integer :: iret, n, VarID, dims

    if (nRec == 0) then
       iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
       dims = 0
       if (present(n1)) then
          iret = nf90_def_dim(ncID, 'zt', n1, ztID)
          iret = nf90_def_dim(ncID, 'zm', n1, zmID)
          dims = 1
       end if
       if (present(n2)) then
          iret = nf90_def_dim(ncID, 'xt', n2, xtID)
          iret = nf90_def_dim(ncID, 'xm', n2, xmID)
          dims = 3
       end if
       if (present(n3)) then
          iret = nf90_def_dim(ncID, 'yt', n3, ytID)
          iret = nf90_def_dim(ncID, 'ym', n3, ymID)
          IF (.not.present(n1)) dims = 2
       end if
       IF (PRESENT(n1a) .AND. PRESENT(n2a) .AND. PRESENT(n2b)) THEN
          iret = nf90_def_dim(ncID, 'B_Rd12a', n1a+n2a, aeaID) ! 1a+2a (a-aerosol only)
          iret = nf90_def_dim(ncID, 'B_Rd2ab', n2b, aebID) ! 2a and 2b (all other species)
       END IF
       IF (PRESENT(nprc)) THEN
          iret = nf90_def_dim(ncID, 'B_Rwprc', nprc, prcID)
       END IF
       IF (PRESENT(nsnw)) THEN
          iret = nf90_def_dim(ncID, 'B_Rwsnw', nsnw, snowID)
       END IF
       IF (PRESENT(nchist)) THEN
          IF (nchist>0) iret = nf90_def_dim(ncID, 'P_hRc', nchist, hcID)
       END IF
       IF (PRESENT(nihist)) THEN
          IF (nihist>0) iret = nf90_def_dim(ncID, 'P_hRi', nihist, hiID)
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
       dim_ttttprc = (/ztID,xtID,ytID,prcID,timeID/)
       dim_ttttsnw = (/ztID,xtID,ytID,snowID,timeID/)
       ! Zubair: dimension environments for avegare size distribution variables per bin - ps files
       dim_ttztaea = (/ztID,aeaID,timeID/)
       dim_ttztaeb = (/ztID,aebID,timeID/)
       dim_ttztprc = (/ztID,prcID,timeID/)
       dim_ttztsnw = (/ztID,snowID,timeID/)
       ! Histograms
       dim_tttzhct = (/ztID,hcID,timeID/)
       dim_tttzhit = (/ztID,hiID,timeID/)
       ! Column statistics
       dim_ttt = (/xtID,ytID,timeID/)

       do n=1,nVar
          select case(trim(ncinfo(2,sx(n),dimensions=dims)))
          case ('time')
             if (sx(n)/='time' .and. dims == 2) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttt,VarID)
             ELSE
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,timeID,VarID)
             ENDIF
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
          case ('B_Rd12a')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,aeaID   ,VarID)
          case ('B_Rd2ab')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,aebID   ,VarID)
          case ('B_Rwprc')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,prcID   ,VarID)
          case ('B_Rwsnw')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,snowID   ,VarID)
          CASE ('hcr')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,hcID    ,VarID)
          CASE ('hir')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,hiID    ,VarID)
          case ('ttttaea')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaea,VarID)
          case ('ttttaeb','ttttcla','ttttclb','ttttica','tttticb')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttaeb,VarID)
          case ('ttttprc')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttprc,VarID)
          case ('ttttsnw')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttttsnw,VarID)
          case ('tttzhct')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttzhct,VarID)
          case ('tttzhit')
             iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttzhit,VarID)
          ! ---
          case ('tttt')
             if (dims == 3) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tttt,VarID)
             ELSEIF (dims == 2) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('mttt')
             if (dims == 3) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('tmtt')
             if (dims == 3) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tmtt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_tt,VarID)
             end if
          case ('ttmt')
             if (dims == 3) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttmt,VarID)
             ELSEif (dims == 2) then
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttt,VarID)
             else
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_mt,VarID)
             end if
          case ('ttztaea')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaea,VarID)
          case ('ttztaeb','ttztcla','ttztclb','ttztica','ttzticb')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztaeb,VarID)
          case ('ttztprc')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztprc,VarID)
          case ('ttztsnw')
                iret=nf90_def_var(ncID,sx(n),NF90_FLOAT,dim_ttztsnw,VarID)
          case default
             if (myid == 0) print *, '  ABORTING: NCIO: Bad dimensional information ',trim(ncinfo(2,sx(n),dimensions=dims))
             call appl_abort(0)
          end select
          iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,sx(n),dimensions=dims))
          iret=nf90_put_att(ncID,VarID,'units'   ,ncinfo(1,sx(n),dimensions=dims))
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
             iret = nf90_inquire_variable(ncID, n, name=xnm)
             IF (xnm /= sx(n)) THEN
                if (myid == 0) print *, '  ABORTING: Incompatible Netcdf File',n,xnm,sx(n)
                call appl_abort(0)
             END IF
          end do
          iret = nf90_sync(ncID)
       end if
    end if

  end subroutine define_nc
  !
  ! ----------------------------------------------------------------------
  ! Subroutine nc_info: Gets long_name, units and dimension info given a
  ! short name.
  !
  character (len=80) function ncinfo(itype,short_name,dimensions)

    character (len=40) :: v_lnm ='scalar xx mixing ratio                  '

    integer, intent (in) :: itype
    character (len=*), intent (in) :: short_name
    integer, optional, intent (in) :: dimensions

    integer :: scalar_number
    INTEGER :: dims

    ! Number of dimensions in addition to time (0=*.ts.nc, 1=*.ps.nc, 3=*.nc)
    dims = 1
    IF (PRESENT(dimensions)) dims = dimensions

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
    case('B_Rd12a')
       if (itype==0) ncinfo = 'Dry size bins, regimes 1a and 2a'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'B_Rd12a'
    case('B_Rd2ab')
       if (itype==0) ncinfo = 'Dry size bins, regime 2a or 2b'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'B_Rd2ab'
    case('B_Rwprc')
       if (itype==0) ncinfo = 'Precipitation size bins'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'B_Rwprc'
    case('B_Rwsnw')
       if (itype==0) ncinfo = 'Snow size bins'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'B_Rwsnw'
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
    case('thi')
       if (itype==0) ncinfo = 'Ice-liquid water potential temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('thl')
       if (itype==0) ncinfo = 'Liquid water potential temperature'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('theta')
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
    case('rv')
       if (itype==0) ncinfo = 'Water vapor mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
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
    case('zi1_bar','zi1')
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
    case('tkeint')
       if (itype==0) ncinfo = 'Vertical integral of total TKE non-weighted'
       if (itype==1) ncinfo = 'm3/s2'
       if (itype==2) ncinfo = 'time'
    case('sfcbflx')
       if (itype==0) ncinfo = 'Surface Buoyancy Flux'
       if (itype==1) ncinfo = 'm/s^2'
       if (itype==2) ncinfo = 'time'
    case('wmax')
       if (itype==0) ncinfo = 'Maximum vertical velocity'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('wbase')
       if (itype==0) ncinfo = 'Cloud base vertical velocity'
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
    case('flx_aer')
       if (itype==0) ncinfo = 'Aerosol sea surface flux'
       if (itype==1) ncinfo = '#/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('flx_iso')
       if (itype==0) ncinfo = 'Isoprene sea surface flux'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('flx_mt')
       if (itype==0) ncinfo = 'Monoterpene sea surface flux'
       if (itype==1) ncinfo = 'kg/m^2/s'
       if (itype==2) ncinfo = 'time'
    case('u10')
       if (itype==0) ncinfo = 'Average 10 m wind speed'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'time'
    case('zi_bar')
       if (itype==0) ncinfo = 'Height of maximum total water mixing ratio gradient'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('wvp_bar')
       if (itype==0) ncinfo = 'Water vapor path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('lwp_bar','lwp')
       if (itype==0) ncinfo = 'Liquid-water path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('lwp_var')
       if (itype==0) ncinfo = 'Liquid-water path variance'
       if (itype==1) ncinfo = 'kg^2/m^4'
       if (itype==2) ncinfo = 'time'
    case('iwp_bar','iwp')
       if (itype==0) ncinfo = 'Ice-water path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('zc')
       if (itype==0) ncinfo = 'Cloud-top height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zb')
       if (itype==0) ncinfo = 'Cloud-base height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zcmn')
       if (itype==0) ncinfo = 'Mean cloud-top height'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'time'
    case('zbmn')
       if (itype==0) ncinfo = 'Mean cloud-base height'
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
    case('imax')
       if (itype==0) ncinfo = 'Maximum ice water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'time'
    case('smax')
       if (itype==0) ncinfo = 'Maximum snow water mixing ratio'
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
    case('swp_bar','swp')
       if (itype==0) ncinfo = 'Snow-water path'
       if (itype==1) ncinfo = 'kg/m^2'
       if (itype==2) ncinfo = 'time'
    case('prcp')
       if (itype==0) ncinfo = 'Surface precipitation rate'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('iprcp')
       if (itype==0) ncinfo = 'Surface ice precipitation rate'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('sprcp')
       if (itype==0) ncinfo = 'Surface snow precipitation rate'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('prcp_bc')
       if (itype==0) ncinfo = 'Below cloud precipitation rate'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'time'
    case('CCN')
       if (itype==0) ncinfo = 'Cloud condensation nuclei'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('ncloud')
       if (itype==0) ncinfo = 'Conditionally sampled cloud number mixing ratio'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nrain')
       if (itype==0) ncinfo = 'Conditionally sampled rain number mixing ratio'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nice')
       if (itype==0) ncinfo = 'Conditionally sampled ice number mixing ratio'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nsnow')
       if (itype==0) ncinfo = 'Conditionally sampled snow number mixing ratio'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'time'
    case('nrcnt')
       if (itype==0) ncinfo = 'Rain cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('nscnt')
       if (itype==0) ncinfo = 'Snow cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('nccnt')
       if (itype==0) ncinfo = 'Cloud cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('nicnt')
       if (itype==0) ncinfo = 'Ice cell counts'
       if (itype==1) ncinfo = '#'
       if (itype==2) ncinfo = 'time'
    case('SS_max')
       if (itype==0) ncinfo = 'Maximum supersaturation'
       if (itype==1) ncinfo = '%'
       if (itype==2) ncinfo = 'time'
    case('SSi_max')
       if (itype==0) ncinfo = 'Maximum supersaturation over ice'
       if (itype==1) ncinfo = '%'
       if (itype==2) ncinfo = 'time'
    case('thl_int')
       if (itype==0) ncinfo = 'Integrated liquid water potential temperature change'
       if (itype==1) ncinfo = 'Km'
       if (itype==2) ncinfo = 'time'
    case('thi_int')
       if (itype==0) ncinfo = 'Integrated ice-liquid water potential temperature'
       if (itype==1) ncinfo = 'Km'
       if (itype==2) ncinfo = 'time'

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
    case('theta_2')
       if (itype==0) ncinfo = 'Variance of theta'
       if (itype==1) ncinfo = 'K^2'
       if (itype==2) ncinfo = 'tttt'
    case('thl_2')
       if (itype==0) ncinfo = 'Variance of liquid water potential temperature'
       if (itype==1) ncinfo = 'K^2'
       if (itype==2) ncinfo = 'tttt'
    case('thi_2')
       if (itype==0) ncinfo = 'Variance of ice-liquid water potential temperature'
       if (itype==1) ncinfo = 'K^2'
       if (itype==2) ncinfo = 'tttt'
    case('w_3')
       if (itype==0) ncinfo = 'Third raw moment of w wind'
       if (itype==1) ncinfo = 'm^3/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('theta_3')
       if (itype==0) ncinfo = 'Third moment of theta'
       if (itype==1) ncinfo = 'K^3'
       if (itype==2) ncinfo = 'tttt'
    case('thl_3')
       if (itype==0) ncinfo = 'Third moment of liquid water potential temperature'
       if (itype==1) ncinfo = 'K^3'
       if (itype==2) ncinfo = 'tttt'
    case('thi_3')
       if (itype==0) ncinfo = 'Third moment of ice-liquid water potential temperature'
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
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of v-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('tot_ww')
       if (itype==0) ncinfo = 'Total vertical flux of w-wind'
       if (itype==1) ncinfo = 'm^2/s^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_ww')
       if (itype==0) ncinfo = 'Sub-filter scale vertical flux of w-wind'
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
       if (itype==0) ncinfo = 'Sub-filter scale buoyancy production of TKE'
       if (itype==1) ncinfo = 'm^2/s^3'
       if (itype==2) ncinfo = 'ttmt'
    case('sfs_shr')
       if (itype==0) ncinfo = 'Sub-filter scale shear production of TKE'
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
       if (itype==0) ncinfo = 'Total vertical flux of water'
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
    case('sw_up')
       if (itype==0) ncinfo = 'Upwelling shortwave radiation'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('sw_down')
       if (itype==0) ncinfo = 'Downwelling shortwave radiation'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
       if (itype==2) ncinfo = 'ttmt'
    case('lw_up')
       if (itype==0) ncinfo = 'Upwelling longwave radiation'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('lw_down')
       if (itype==0) ncinfo = 'Downwelling longwave radiation'
       if (itype==1) ncinfo = 'W/m^2'
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
       if (itype==0) ncinfo = 'Total vertical flux of liquid water'
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
       if (itype==0) ncinfo = 'Average of w over cs1'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs1')
       if (itype==0) ncinfo = 'Average of theta_l over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'ttmt'
    case('tv_cs1')
       if (itype==0) ncinfo = 'Average of theta_v over cs1'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs1')
       if (itype==0) ncinfo = 'Average of total water over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rc_cs1')
       if (itype==0) ncinfo = 'Average of total condensate over cs1'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wtl_cs1')
       if (itype==0) ncinfo = 'Average of vertical theta_l flux over cs1'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wtv_cs1')
       if (itype==0) ncinfo = 'Average of vertical theta_v flux over cs1'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wrt_cs1')
       if (itype==0) ncinfo = 'Average of vertical total water flux over cs1'
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
       if (itype==0) ncinfo = 'Average of w over cs2'
       if (itype==1) ncinfo = 'm/s'
       if (itype==2) ncinfo = 'ttmt'
    case('tl_cs2')
       if (itype==0) ncinfo = 'Average of theta_l over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('tv_cs2')
       if (itype==0) ncinfo = 'Average of theta_v over cs2'
       if (itype==1) ncinfo = 'K'
       if (itype==2) ncinfo = 'tttt'
    case('rt_cs2')
       if (itype==0) ncinfo = 'Average of total water over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rc_cs2')
       if (itype==0) ncinfo = 'Average of total condensate over cs2'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('wtl_cs2')
       if (itype==0) ncinfo = 'Average of vertical theta_l flux over cs2'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wtv_cs2')
       if (itype==0) ncinfo = 'Average of vertical theta_v flux over cs2'
       if (itype==1) ncinfo = 'K*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('wrt_cs2')
       if (itype==0) ncinfo = 'Average of vertical total water flux over cs2'
       if (itype==1) ncinfo = 'kg/kg*m/s'
       if (itype==2) ncinfo = 'ttmt'
    case('Nc_ic')
       if (itype==0) ncinfo = 'Conditionally sampled cloud droplet number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Nr_ir')
       if (itype==0) ncinfo = 'Conditionally sampled rain drop number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Ni_ii')
       if (itype==0) ncinfo = 'Conditionally sampled ice number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Ns_is')
       if (itype==0) ncinfo = 'Conditionally sampled snow number concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('frac_ic')
       if (itype==0) ncinfo = 'Cloud fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('frac_ir')
       if (itype==0) ncinfo = 'Rain fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('frac_ii')
       if (itype==0) ncinfo = 'Ice fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('frac_is')
       if (itype==0) ncinfo = 'Snow fraction'
       if (itype==1) ncinfo = '-'
       if (itype==2) ncinfo = 'tttt'
    case('rr')
       if (itype==0) ncinfo = 'Rain water mixing ratio'
       if (itype==1) ncinfo = 'kg/kg'
       if (itype==2) ncinfo = 'tttt'
    case('rrate')
       if (itype==0) ncinfo = 'Rain water deposition flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('crate')
       if (itype==0) ncinfo = 'Cloud water deposition flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('irate')
       if (itype==0) ncinfo = 'Ice water deposition flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('srate')
       if (itype==0) ncinfo = 'Snow water deposition flux'
       if (itype==1) ncinfo = 'W/m^2'
       if (itype==2) ncinfo = 'ttmt'
    case('evap')
       if (itype==0) ncinfo = 'Net evap of rain-water'
       if (itype==1) ncinfo = 's^-1'
       if (itype==2) ncinfo = 'tttt'
    case('rho_air')
       if (itype==0) ncinfo = 'Air density'
       if (itype==1) ncinfo = 'kg/m^3'
       if (itype==2) ncinfo = 'tttt'

    ! SALSA column statistics
    case('Nc')
       if (itype==0) ncinfo = 'Number-mean cloud droplet concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Nr')
       if (itype==0) ncinfo = 'Number-mean rain drop concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Ni')
       if (itype==0) ncinfo = 'Number-mean ice concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Ns')
       if (itype==0) ncinfo = 'Number-mean snow concentration'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttt'
    case('Rwc')
       if (itype==0) ncinfo = 'Number-mean cloud droplet radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('Rwr')
       if (itype==0) ncinfo = 'Number-mean rain drop radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('Rwi')
       if (itype==0) ncinfo = 'Number-mean ice radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'
    case('Rws')
       if (itype==0) ncinfo = 'Number-mean snow radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'tttt'

    ! SALSA profile statistics
    case('P_hNca')
       if (itype==0) ncinfo = 'Cloud droplets per radius bin (a bins)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttzhct'
    case('P_hNcb')
       if (itype==0) ncinfo = 'Cloud droplets per radius bin (b bins)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttzhct'
    case('P_hRc')
       if (itype==0) ncinfo = 'Cloud histogram bin mean radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'hcr'
    case('P_hNia')
       if (itype==0) ncinfo = 'Ice particles per radius bin (a bins)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttzhit'
    case('P_hNib')
       if (itype==0) ncinfo = 'Ice particles per radius bin (b bins)'
       if (itype==1) ncinfo = 'kg^-1'
       if (itype==2) ncinfo = 'tttzhit'
    case('P_hRi')
       if (itype==0) ncinfo = 'Ice histogram bin mean radius'
       if (itype==1) ncinfo = 'm'
       if (itype==2) ncinfo = 'hir'
    ! -----
    case default
       ! Automatically generated microphysical process rate statistics
       ncinfo=TRIM( get_rate_info(itype,trim(short_name),dims) )
       ! ... and some other species and bin-dependent SALSA variables
       IF (LEN(TRIM(ncinfo))<1) ncinfo=TRIM( get_salsa_info(itype,trim(short_name),dims) )
       ! Additional scalars (3D only)
       IF (LEN(TRIM(ncinfo))<1 .AND. short_name(1:1)=='s' .AND. dims==3) THEN
          read (UNIT=short_name(2:3),FMT='(i2.2)',iostat=dims) scalar_number
          IF (dims==0) THEN ! If valid name
             write(v_lnm(8:9),'(i2.2)') scalar_number
             if (itype==0) ncinfo = v_lnm
             if (itype==1) ncinfo = 'kg/kg'
             if (itype==2) ncinfo = 'tttt'
          ENDIF
       ENDIF
       IF (LEN(TRIM(ncinfo))<1) THEN
          if (myid==0) print *, 'ABORTING: ncinfo: variable not found ',trim(short_name)
          call appl_abort(0)
       END IF
    end select

  end function ncinfo

  !
  ! ----------------------------------------------------------------------
  ! Function that determines information related to microphysics process rate statistics.
  ! These have the same variable name for all outputs, but the units can be different.
  !
  ! 1) Units
  !     *.ts.nc: average column integrated rate of change
  !     *.ps.nc: average rate of change per volume
  !     *.nc: rate of change
  ! 2) Examples of names
  !     diag_ri     Change in ice water mixing ratio due to diagnostics
  !     diag_ni     Change in ice water number concentration due to diagnostics
  !     diag_tt     Change in [total] temperature ...
  !     r=mixing ratio, n=number concentration, t=temperature,..
  !     i=ice, t=total,..
  !
  character (len=80) function get_rate_info(itype,short_name,dims)
    implicit none
    integer, intent (in) :: itype, dims
    character (len=*), intent (in) :: short_name

    character (len=20) :: pros, spec
    integer :: i
    logical :: numc, mixr, temp
    CHARACTER :: idc

    IF (LEN(TRIM(short_name))<6 .OR. INDEX(short_name,'_')==0) THEN
        get_rate_info=''
        RETURN
    ENDIF

    ! Which process
    IF ('coag_'==short_name(1:5)) THEN
        pros='coagulation'
    ELSEIF ('oxid_'==short_name(1:5)) THEN
        pros='oxidation'
    ELSEIF ('ocon_'==short_name(1:5)) THEN
        pros='condensation'
    ELSEIF ('cond_'==short_name(1:5)) THEN
        pros='condensation'
    ELSEIF ('auto_'==short_name(1:5)) THEN
        pros='autoconversion'
    ELSEIF ('act_'==short_name(1:4) .OR. 'cact_'==short_name(1:5)) THEN
        pros='activation'
    ELSEIF ('nucl_'==short_name(1:5)) THEN
        pros='nucleation'
    ELSEIF ('melt_'==short_name(1:5)) THEN
        pros='melting'
    ELSEIF ('dist_'==short_name(1:5)) THEN
        pros='distribution update'
    ELSEIF ('sedi_'==short_name(1:5)) THEN
        pros='sedimentation'
    ELSEIF ('diag_'==short_name(1:5)) THEN
        pros='diagnostics'
    ELSEIF ('advf_'==short_name(1:5)) THEN
        pros='advection'
    ELSEIF ('srfc_'==short_name(1:5)) THEN
        pros='surface'
    ELSEIF ('diff_'==short_name(1:5)) THEN
        pros='diffusion'
    ELSEIF ('forc_'==short_name(1:5)) THEN
        pros='forcings'
    ELSEIF ('mcrp_'==short_name(1:5)) THEN
        pros='microphysics'
    ELSEIF ('nudg_'==short_name(1:5)) THEN
        pros='nudging'
    ELSE
        get_rate_info=''
        RETURN
    ENDIF

    ! Number concentration (n), water mixing ratio (r) or mixing ratio of component x (x, where x=1,2,...,9)
    ! temperature (t)
    numc=.FALSE.
    mixr=.FALSE.
    temp=.FALSE.
    i=LEN(TRIM(short_name))-1
    select case (short_name(i:i))
    CASE('n')
        numc=.TRUE.
    CASE('r')
        mixr=.TRUE.
    CASE('t')
        temp=.TRUE.
    CASE('1','2','3','4','5','6','7','8','9')
        idc=short_name(i:i)
    case default
        get_rate_info=''
        RETURN
    end select

    ! Species (a, c, r, i or s)
    i=i+1
    select case (short_name(i:i))
    CASE('a')
        spec='aerosol'
    CASE('c')
        spec='cloud'
    CASE('r')
        spec='rain'
    CASE('i')
        spec='ice'
    CASE('s')
        spec='snow'
    CASE('g')
        spec='gas'
    CASE('t')
        spec='total'
    case default
        get_rate_info=''
        RETURN
    end select

    ! Valid microphysical process identified, formulate the output
    if (itype==0) THEN
        ! Long name
        IF (dims==0 .OR. dims==2) THEN
            ! Vertical integral of number concentration or mixing ratio (#/m^2/s or kg/m^2/s)
            IF (numc) THEN
                get_rate_info='Change in column '//TRIM(spec)//' number due to '//TRIM(pros)
            ELSEIF (mixr) THEN
                get_rate_info='Change in column '//TRIM(spec)//' water due to '//TRIM(pros)
            ELSEIF (temp) THEN
                get_rate_info='Change in column '//TRIM(spec)//' temperature due to '//TRIM(pros)
            ELSE
                get_rate_info='Change in column '//TRIM(spec)//' component #'//idc//' mixing ratio due to '//TRIM(pros)
            ENDIF
        ELSEIF (dims==1 .OR. dims==3) THEN
            ! Profiles give the average rate per volume (#/m^3/s or kg/m^3/s), so just different unit
            ! 3D data in original units (#/kg/s or kg/kg/s)
            IF (numc) THEN
                get_rate_info='Change in '//TRIM(spec)//' number concentration due to '//TRIM(pros)
            ELSEIF (mixr) THEN
                get_rate_info='Change in '//TRIM(spec)//' water mixing ratio due to '//TRIM(pros)
            ELSEIF (temp) THEN
                get_rate_info='Change in '//TRIM(spec)//' temperature due to '//TRIM(pros)
            ELSE
                get_rate_info='Change in '//TRIM(spec)//' component #'//idc//' mixing ratio due to '//TRIM(pros)
            ENDIF
        ELSE
            get_rate_info = ''
        ENDIF
    ELSEIF (itype==1) THEN
        ! Unit
        IF (dims==0 .OR. dims==2) THEN
            ! Vertical integral of number concentration or mixing ratio (#/m^2/s or kg/m^2/s)
            IF (numc) THEN
                get_rate_info = '#/m^2/s'
            ELSEIF (temp) THEN
                get_rate_info = 'W/m^2'
            ELSE
                get_rate_info = 'kg/m^2/s'
            ENDIF
        ELSEIF (dims==1) THEN
            ! Profiles give the average rate per volume (#/m^3/s or kg/m^3/s)
            IF (numc) THEN
                get_rate_info = '#/m^3/s'
            ELSEIF (temp) THEN
                get_rate_info = 'W/m^3'
            ELSE
                get_rate_info = 'kg/m^3/s'
            ENDIF
        ELSEIF (dims==3) THEN
            ! 3D data in original units (#/kg/s or kg/kg/s)
            IF (numc) THEN
                get_rate_info = '#/kg/s'
            ELSEIF (temp) THEN
                get_rate_info = 'W/kg'
            ELSE
                get_rate_info = 'kg/kg/s'
            ENDIF
        ELSE
            get_rate_info = ''
        ENDIF
    ELSEIF (itype==2) THEN
        ! NetCDF dimensions
        IF (dims==0) THEN
            ! Time series
            get_rate_info = 'time'
        ELSEIF (dims==1) THEN
            ! Profiles
            get_rate_info = 'tttt'
        ELSEIF (dims==3) THEN
            ! 3D outputs
            get_rate_info = 'tttt'
        ELSEIF (dims==2) THEN
            ! 2D outputs
            get_rate_info = 'tttt'
        ELSE
            get_rate_info = ''
        ENDIF
    ELSE
        get_rate_info = ''
    ENDIF

  END function get_rate_info
  !
  ! ----------------------------------------------------------------------
  ! Function that determines information related to SALSA variables.
  ! 1) Time series and column stats
  !      Removal rates (e.g. rmNOcl)
  !      Vertically integrated or averaged (conditional) mass/number concentrations and sizes (e.g. tcNOaa, icNOaa)
  ! 2) Profiles
  !     Horizontally averaged (conditional) concentrations (e.g. tcNOaa, icNOat, icNcb)
  !     Horizontally averaged (conditional) concentration in each bin (e.g. B_OCaa)
  ! 3) Analysis files
  !     Grid cell total mass and number concentrations and mean radii (e.g. C_H2Ost, C_Nab, C_Rwca)
  !     Grid cell bin mass concentrations and mean radii (e.g. B_H2Oaa, B_Naa, B_Rwca)
  character (len=80) function get_salsa_info(itype,short_name,dims)
    implicit none
    integer, intent (in) :: itype
    character (len=*), intent (in) :: short_name
    INTEGER, INTENT(IN) :: dims

    character (len=20) :: phase, tmp
    character (len=3) :: spec
    integer :: i

    get_salsa_info = ''
    i=LEN(TRIM(short_name))
    spec='   '

    ! All SALSA species: 'SO4','OC ','BC ','DU ','SS ','NO ','NH ','H2O'

    IF (i<5) THEN
        RETURN
    ELSEIF (INDEX(short_name,'tc')==1 .OR. INDEX(short_name,'ic')==1 .OR. &
            INDEX(short_name,'ir')==1 .OR. INDEX(short_name,'ii')==1 .OR. &
            INDEX(short_name,'is')==1 .OR. INDEX(short_name,'C_')==1) THEN
        ! Averaged/integrated conditional mass concentrations and particle mean wet radii and number concentrations
        !   name=<condition>//<species>//<bin>
        !       condition: tc=total, ic=in cloud, ir=in rain, ii=in ice, is=in snow, C_=grid-cell mean 3D analysis data
        !       species: aerosol (SO4, OC, BC, ..) or gas name for mass concentration,
        !                R for mean radius and N for number concentration
        !       bin: aa,ab,at,ca,cb,ct,rt,ia,ib,it,st,gt (aerosol/cloud/rain/ice/snow/gas a-bins/b-bins/total)
        ! Mass concentrations are integrated vertically for 1D (ts) and 3D (cs) outputs.
        ! Otherwise values (mass, N or R) are just averaged.
        !
        ! Species name (or N or Rw)
        spec=short_name(3:i-2)
        ! .. .and bin
        SELECT CASE (short_name(i:i))
        CASE('t')
            tmp='total'
        CASE('a')
            tmp='a-bin'
        CASE('b')
            tmp='b-bin'
        case default
            RETURN
        END SELECT
        SELECT CASE (short_name(i-1:i-1))
        CASE('a')
            phase=TRIM(tmp)//' aerosol'
        CASE('c')
            phase=TRIM(tmp)//' cloud'
        CASE('r','p')
            phase=TRIM(tmp)//' rain'
        CASE('i')
            phase=TRIM(tmp)//' ice'
        CASE('s')
            phase=TRIM(tmp)//' snow'
        CASE('g')
            ! Gas: <bin> should be just 'gt' meaning total gas concentration
            phase='vapor phase'
        case default
            RETURN
        END SELECT
        !
        ! Condition
        SELECT CASE (short_name(1:2))
        CASE('tc')
            tmp=' '
        CASE('ic')
            tmp=' in cloud'
        CASE('ir')
            tmp=' in rain'
        CASE('ii')
            tmp=' in ice'
        CASE('is')
            tmp=' in snow'
        CASE('C_') ! 3D analysis data
            tmp='  '
        case default
            RETURN
        END SELECT
        !
        ! Generate output
        if (itype==0) THEN
            ! Long name
            IF (spec=='N  ') THEN
                get_salsa_info = 'Mean number concentration of '//TRIM(phase)//tmp
            ELSEIF (spec=='R  ' .OR. spec=='Rw ') THEN
                get_salsa_info = 'Mean radius of '//TRIM(phase)//tmp
            ELSE
                IF (dims==0 .OR. dims==2) THEN
                    get_salsa_info = 'Vertically integrated '//TRIM(phase)//' '//TRIM(spec)//' mass'//tmp
                ELSE
                    get_salsa_info = 'Mean '//TRIM(phase)//' '//TRIM(spec)//' mass'//tmp
                ENDIF
            ENDIF
        ELSEIF (itype==1) THEN
            ! Unit
            IF (spec=='N  ') THEN
                get_salsa_info = '#/kg'
            ELSEIF (spec=='R  ' .OR. spec=='Rw ') THEN
                get_salsa_info = 'm'
            ELSE
                IF (dims==0 .OR. dims==2) THEN
                    get_salsa_info = 'kg/m^2'
                ELSE
                    get_salsa_info = 'kg/kg'
                ENDIF
            ENDIF
        ELSEIF (itype==2) THEN
            ! NetCDF dimensions
            IF (dims==0) THEN
                get_salsa_info = 'time'
            ELSE
                get_salsa_info = 'tttt'
            ENDIF
        ENDIF
    ELSEIF (INDEX(short_name,'B_')==1 .AND. (dims==1 .OR. dims==3)) THEN
        ! Bin-dependent data (profiles or analysis)
        !   name='B_'//<species>//<bin>
        !       species: aerosol (SO4, OC, BC, ..) name for mass concentration,
        !                R for mean radius and N for number concentration
        !       bin: aa,ab,ca,cb,rt,ia,ib,st (aerosol/cloud/rain/ice/snow a-bins/b-bins/total)
        !
        ! Species name (or N or R)
        spec=short_name(3:i-2)
        ! .. .and bin
        SELECT CASE (short_name(i-1:i))
        CASE('aa')
            phase='Aerosol a-bin'
            tmp='aea' ! Analysis data starts with 'tttt' and profile data starts with 'ttzt'
        CASE('ab')
            phase='Aerosol b-bin'
            tmp='aeb'
        CASE('ca')
            phase='Cloud a-bin'
            tmp='cla'
        CASE('cb')
            phase='Cloud b-bin'
            tmp='clb'
        CASE('rt','pt')
            phase='Rain bin'
            tmp='prc'
        CASE('ia')
            phase='Ice a-bin'
            tmp='ica'
        CASE('ib')
            phase='Ice a-bin'
            tmp='icb'
        CASE('st')
            phase='Snow bin'
            tmp='snw'
        case default
            RETURN
        END SELECT
        !
        ! Generate output
        if (itype==0) THEN
            ! Long name
            IF (spec=='N  ') THEN
                get_salsa_info = TRIM(phase)//' number concentration'
            ELSEIF (spec=='R  ' .OR. spec=='Rw ') THEN
                get_salsa_info = TRIM(phase)//' wet radius'
            ELSE
                get_salsa_info = TRIM(phase)//' '//TRIM(spec)//' mass mixing ratio'
            ENDIF
        ELSEIF (itype==1) THEN
            ! Unit
            IF (spec=='N  ') THEN
                get_salsa_info = '#/kg'
            ELSEIF (spec=='R  ' .OR. spec=='Rw ') THEN
                get_salsa_info = 'm'
            ELSE
                get_salsa_info = 'kg/kg'
            ENDIF
        ELSEIF (itype==2) THEN
            ! NetCDF dimensions
            IF (dims==1) THEN
                get_salsa_info = 'ttzt'//TRIM(tmp)
            ELSE
                get_salsa_info = 'tttt'//TRIM(tmp)
            ENDIF
        ENDIF
    ELSEIF (INDEX(short_name,'rm')==1) THEN
        ! Removal rate time series
        !   name='rm'//<spec>//<x> where x='dr', 'cl', 'pr', 'ic' or 'sn'
        ! Species name
        spec(1:i-4)=short_name(3:i-2)
        ! Phase
        select case (short_name(i-1:i))
        CASE('dr')
            phase='aerosol'
        CASE('cl')
            phase='cloud'
        CASE('pr')
            phase='rain'
        CASE('ic')
            phase='ice'
        CASE('sn')
            phase='snow'
        case default
            RETURN
        end select
        ! Generate output
        if (itype==0) THEN
            ! Long name
            get_salsa_info = 'Removal of '//TRIM(spec)//' with '//TRIM(phase)
        ELSEIF (itype==1) THEN
            ! Unit
            get_salsa_info = 'kg/m^2/s'
        ELSEIF (itype==2) THEN
            ! NetCDF dimensions
            get_salsa_info = 'time'
        ENDIF
    ENDIF
    END function get_salsa_info
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
