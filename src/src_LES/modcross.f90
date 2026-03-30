!===============================================================================!
!
! 2D cross section outputs for UCLALES-SALSA
!
! January 15, 2026: modified the original UCLALES codes obtained from
! https://github.com/uclales/uclales/blob/master/src/ice_sb.F90
! https://github.com/uclales/uclales/blob/master/src/modcross.f90
! (last access: January 15, 2026) to work with UCLALES-SALSA.
!
! Tomi Raatikainen, Finnish Meteorological Institute, Helsinki, Finland
! (tomi.raatikainen@fmi.fi)
!
!===============================================================================!
!
!  This file is part of MicroHH.
!
! MicroHH is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! MicroHH is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!  Copyright 2010-2011 Chiel van Heerwaarden and Thijs Heus
!> All kinds of auxilarily functions and variables for
module modcross
  use grid, only : maxn_list
  implicit none
  ! NAMELIST inputs
  logical, public :: lcross = .false., lxy = .false., lxz = .false., lyz = .false.
  real, public    :: xcross = 0., ycross = 0., zcross(10) = 0., frqcross = 3600.
  character (len=7), dimension(maxn_list), public :: crossvars='       '

  PRIVATE
  integer, SAVE :: icross, jcross, kcross(10), nkcross, ncrossvars
  integer, SAVE :: nccrossxzid, nccrossyzid, nccrossxyid, nccrossrec

  PUBLIC  triggercross, initcross, close_cross
contains

  subroutine initcross(rtimee, expname)
    USE mpi_interface, ONLY : myid
    use ncio, only : open_nc
    use grid, only : nzp, nxp, nyp, zt, xt, yt, xm, ym
    real, intent(in) :: rtimee
    character(len=*), intent(in) :: expname
    integer :: n, i, j, k
    character(len=80) :: fname, title

    ! Update variables
    ncrossvars=0 ! Total number
    DO i=1,SIZE(crossvars)
        IF (LEN_TRIM(crossvars(i))>0) THEN
            ncrossvars=ncrossvars+1
            IF (ncrossvars<i) crossvars(ncrossvars)=crossvars(i)
        ENDIF
    ENDDO
    IF (ncrossvars==0) THEN
        lcross = .false.
        RETURN
    ENDIF

    if(myid==0) print "(//' ',49('-')/,' ')"

    if (lxy) then
      nkcross = 0
      do n = 1, 10
        IF (zcross(n)>0.0) THEN
            ! Find the closest thermo point
            k=1
            do i=2,nzp
                if (abs(zt(i)-zcross(n))<abs(zt(k)-zcross(n))) k=i
            end do
            IF (k>1 .AND. k<nzp) THEN
                ! Valid
                nkcross = nkcross + 1
                kcross(nkcross) = k
                zcross(nkcross) = zt(k) ! The actual height level
            ELSEIF (myid==0) THEN
                write(*,'("     Module modcross: ignoring z=",i5," m")') nint(zcross(n))
            ENDIF
        ENDIF
      end do
      !
      if (nkcross>0) THEN
        IF (nkcross>1) THEN
            title = 'xy cross sections'
        ELSE
            WRITE(title,'("xy cross sections at z=",i5," m")') nint(zt(kcross(1)))
        ENDIF
        fname = trim(expname)//'.xy'
        if(myid==0) print "('   Initializing: ',A20,'  N=',I3)", fname, ncrossvars
        call open_nc(fname, title, rtimee, nccrossxyid, nccrossrec)
        ! Dimensions: xt, yt (and zt=zcross)
        CALL define_nc(nccrossxyid, nccrossrec, nkcross, nxp-4, nyp-4)
      ELSE
        lxy = .FALSE.
      end if
    end if

    if (lxz) then
      if (ycross < ym(2) .or. ycross >= ym(nyp - 2)) then
        lxz = .false.
      else
        j=3
        do i=4,nyp-2
          if (abs(yt(i)-ycross)<abs(yt(j)-ycross)) j=i
        end do
        jcross = j
        fname = trim(expname)//'.xz'
        WRITE(title,'("xz cross sections at y=",i5," m")') nint(yt(j))
        if(myid==0) print "('   Initializing: ',A20,'  N=',I3)", fname, ncrossvars
        call open_nc(fname, title, rtimee, nccrossxzid, nccrossrec)
        ! Dimensions: zt, xt
        CALL define_nc(nccrossxzid, nccrossrec, nzp, nxp-4, 0)
      end if
    end if

    if (lyz) then
      if (xcross < xm(2) .or. xcross >= xm(nxp - 2)) then
        lyz = .false.
      else
        i=3
        do j=4,nxp-2
          if (abs(xt(j)-xcross)<abs(xt(i)-xcross)) i=j
        end do
        icross = i
        fname = trim(expname)//'.yz'
        WRITE(title,'("yz cross sections at x=",i5," m")') nint(xt(i))
        if(myid==0) print "('   Initializing: ',A20,'  N=',I3)", fname, ncrossvars
        call open_nc(fname, title, rtimee, nccrossyzid, nccrossrec)
        ! Dimensions: zt, yt
        CALL define_nc(nccrossyzid, nccrossrec, nzp, 0, nyp-4)
      end if
    end if

    ! Ready to write the first record (the same for all output files)
    IF (nccrossrec==0) nccrossrec=1

    if (myid == 0) print *,'   ...starting record: ', nccrossrec

  end subroutine initcross


  SUBROUTINE define_nc(ncID, nRec, n1, n2, n3)
    USE mpi_interface, ONLY : myid, appl_abort
    USE ncio, ONLY : ncinfo
    use grid, only : nzp, nxp, nyp, zt, xt, yt
    USE netcdf, ONLY : nf90_def_dim, nf90_def_var, nf90_put_att, nf90_enddef, &
        nf90_sync, nf90_inquire_variable, nf90_inquire, nf90_float, nf90_unlimited, &
        nf90_inq_varid, nf90_put_var
    integer, intent(in) :: ncID, nRec, n1, n2, n3
    INTEGER :: timeID, ztID, xtID, ytID, ids(4), VarID, iret
    INTEGER :: i, n
    character(len=7) :: name
    character (len=80) :: longname, units
    !
    IF (nRec==0) THEN
        ! Dimensions (time and x-y, x-z or y-z)
        iret = nf90_def_dim(ncID, 'time', NF90_UNLIMITED, timeID)
        IF (n1>1) iret = nf90_def_dim(ncID, 'zt', n1, ztID) ! n1=nkcross or n1=nzp
        IF (n2>0) iret = nf90_def_dim(ncID, 'xt', n2, xtID)
        IF (n3>0) iret = nf90_def_dim(ncID, 'yt', n3, ytID)
        !
        ! Dimension variables
        n=0
        IF (n1>1) THEN
            iret=nf90_def_var(ncID,'zt',NF90_FLOAT,ztID,VarID)
            iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'zt'))
            iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'zt'))
            n=n+1
            ids(n)=ztID
        ENDIF
        IF (n2>0) THEN
            iret=nf90_def_var(ncID,'xt',NF90_FLOAT,xtID,VarID)
            iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'xt'))
            iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'xt'))
            n=n+1
            ids(n)=xtID
        ENDIF
        IF (n3>0) THEN
            iret=nf90_def_var(ncID,'yt',NF90_FLOAT,ytID,VarID)
            iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'yt'))
            iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'yt'))
            n=n+1
            ids(n)=ytID
        ENDIF
        iret=nf90_def_var(ncID,'time',NF90_FLOAT,timeID,VarID)
        iret=nf90_put_att(ncID,VarID,'longname',ncinfo(0,'time'))
        iret=nf90_put_att(ncID,VarID,'units',ncinfo(1,'time'))
        n=n+1
        ids(n)=timeID
        !
        ! Active variables
        DO i=1,ncrossvars
          name=crossvars(i)
          longname=ncinfo(0,name,dimensions=3)
          units=ncinfo(1,name,dimensions=3)
          iret=nf90_def_var(ncID,name,NF90_FLOAT,ids(1:n),VarID)
          iret=nf90_put_att(ncID,VarID,'longname',longname)
          iret=nf90_put_att(ncID,VarID,'units',units)
        ENDDO
        !
        ! Done defining
        iret  = nf90_enddef(ncID)
        !
        ! Put dimensions
        IF (n1==nzp) THEN
            iret = nf90_inq_varid(ncID, 'zt',VarID)
            iret = nf90_put_var(ncID, VarID, zt)
        ELSEIF (n1==nkcross) THEN
            iret = nf90_inq_varid(ncID, 'zt',VarID)
            iret = nf90_put_var(ncID, VarID, zcross(1:nkcross))
        ENDIF
        IF (n2>0) THEN
            iret = nf90_inq_varid(ncID, 'xt',VarID)
            iret = nf90_put_var(ncID, VarID, xt(3:nxp-2))
        ENDIF
        IF (n3>0) THEN
            iret = nf90_inq_varid(ncID, 'yt',VarID)
            iret = nf90_put_var(ncID, VarID, yt(3:nyp-2))
        ENDIF
        !
        iret  = nf90_sync(ncID)
    ELSE
        ! Dimensions should be there, but confim the other variables
        i=COUNT((/n1>1,n2>0,n3>0/))+1 ! The number of dimensions
        iret = nf90_inquire(ncID, nVariables=n)
        if (n /= i+ncrossvars) then
            if (myid == 0) print *, '  ABORTING: Incompatible Netcdf File',n,i+ncrossvars
            call appl_abort(0)
        else
            DO n=1,ncrossvars
                i=i+1
                iret = nf90_inquire_variable(ncID, i, name=name)
                IF (name /= crossvars(n)) THEN
                    if (myid == 0) print *, '  ABORTING: Incompatible Netcdf File', &
                                        n,crossvars(n),i,NAME
                    call appl_abort(0)
                END IF
            ENDDO
        ENDIF
        iret = nf90_sync(ncID)
    end if
    !
  END SUBROUTINE define_nc


  subroutine triggercross(rtimee)
    USE mpi_interface, ONLY : myid
    use grid, only : level, nxp, nyp, nzp, dzm, dzt, a_up, a_vp, a_wp, umean, vmean, & 
        a_press, a_qp, a_theta, a_temp, a_tp, a_rflx, a_sflx, a_fus, a_fds, a_fuir, a_fdir, &
        a_rv, a_rsl, a_rsi, a_rp, a_rc, a_ri, a_ap, &
        ccn, a_rpp, a_npp, a_rip, a_nip, a_rsp, a_nsp, a_rgp, a_ngp, a_rhp, a_nhp, & ! SB microphysics
        a_ncloudp, a_mcloudp, a_nprecpp, a_mprecpp, a_nicep, a_micep, a_nsnowp, a_msnowp, & ! SALSA
        ncld, nprc, nice, nsnw
    USE defs, ONLY : cp, alvi
    USE stat, ONLY : calc_user_data, sflg
    USE netcdf, only : nf90_inq_varid, nf90_put_var
    real, intent(in) :: rtimee
    real, dimension(nzp,nxp,nyp) :: interp
    integer :: n, i, j, k, iret, VarID
    LOGICAL :: fail, mask(nzp,nxp,nyp)

    ! Time
    if (lxy) then
        iret = nf90_inq_varid(nccrossxyid, 'time',VarID)
        iret = nf90_put_var(nccrossxyid, VarID, rtimee, start=(/nccrossrec/))
    end if
    if (lxz) then
        iret = nf90_inq_varid(nccrossxzid, 'time',VarID)
        iret = nf90_put_var(nccrossxzid, VarID, rtimee, start=(/nccrossrec/))
    end if
    if (lyz) then
        iret = nf90_inq_varid(nccrossyzid, 'time',VarID)
        iret = nf90_put_var(nccrossyzid, VarID, rtimee, start=(/nccrossrec/))
    end if

    do n = 1, ncrossvars
      select case(trim(crossvars(n)))
      case('u')
        do j=3,nyp-2
          do i=3,nxp-2
              do k=1,nzp
                interp(k,i,j) = 0.5*(a_up(k,i-1,j) + a_up(k,i,j)) + umean
              end do
          end do
        end do
        call writecross_3D(crossvars(n), interp)
      case('v')
        do j=3,nyp-2
          do i=3,nxp-2
              do k=1,nzp
                interp(k,i,j) = 0.5*(a_vp(k,i,j-1) + a_vp(k,i,j)) + vmean
              end do
          end do
        end do
        call writecross_3D(crossvars(n), interp)
      case('w')
        do j=3,nyp-2
          do i=3,nxp-2
              interp(1,i,j) = 0.0
              do k=2,nzp
                interp(k,i,j) = 0.5*dzt(k) * (a_wp(k-1,i,j) / dzm(k) + a_wp(k,i,j) / dzm(k-1))
              end do
          end do
        end do
        call writecross_3D(crossvars(n), interp)
      case('theta') ! Potential temperature
        call writecross_3D(crossvars(n), a_theta)
      case('thl') ! Liquid water potential temperature
        IF (level==0 .OR. level==5) THEN
            WHERE(a_temp>0.) interp = a_tp + (a_theta/a_temp)*alvi/cp*a_ri
            call writecross_3D(crossvars(n), interp)
        ELSE
            call writecross_3D(crossvars(n), a_tp)
        ENDIF
      case('thi') ! Ice-liquid water potential temperature
        call writecross_3D(crossvars(n), a_tp)
      case('temp') ! Absolute temperature
        call writecross_3D(crossvars(n), a_temp)
      case ('SS') ! Supersaturation (%) over liquid water
            interp=0.
            IF (level<4) THEN
                WHERE(a_rsl>1e-10) interp=(a_rv/a_rsl-1.0)*100.
            ELSE
                WHERE(a_rsl>1e-10) interp=(a_rp/a_rsl-1.0)*100.
            ENDIF
            call writecross_3D(crossvars(n), interp)
      case ('SSi') ! Supersaturation (%) over ice
            interp=0.
            IF (level<4) THEN
                WHERE(a_rsi>1e-10) interp=(a_rv/a_rsi-1.0)*100.
            ELSE
                WHERE(a_rsi>1e-10) interp=(a_rp/a_rsi-1.0)*100.
            ENDIF
            call writecross_3D(crossvars(n), interp)
      case('p')
        call writecross_3D(crossvars(n), a_press)
      case('rflx')
        call writecross_3D(crossvars(n), a_rflx)
      case('sflx')
        call writecross_3D(crossvars(n), a_sflx)
      case('sw_up')
        call writecross_3D(crossvars(n), a_fus(1:nzp,:,:))
      case('sw_down')
        call writecross_3D(crossvars(n), a_fds(1:nzp,:,:))
      case('lw_up')
        call writecross_3D(crossvars(n), a_fuir(1:nzp,:,:))
      case('lw_down')
        call writecross_3D(crossvars(n), a_fdir(1:nzp,:,:))
      case('stke')
        call writecross_3D(crossvars(n), a_qp)
      case('q') ! Total water
        IF (level<4) THEN
            call writecross_3D(crossvars(n), a_rp)
        ELSEIF (level==4) THEN
            interp = a_rp + a_rc
            call writecross_3D(crossvars(n), interp)
        ELSEIF (level==5) THEN
            interp = a_rp + a_rc + a_ri
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('rc') ! cloud water
        IF (level<4) THEN
            call writecross_3D(crossvars(n), a_rc)
        ELSE
            interp = SUM(a_mcloudp(:,:,:,1:ncld),DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('rv') ! Water vapor
        IF (level<4) THEN
            call writecross_3D(crossvars(n), a_rv)
        ELSE
            call writecross_3D(crossvars(n), a_rp)
        ENDIF
      case('nc') ! CDNC
        IF (level<4) THEN
            interp = CCN
        ELSE
            interp = SUM(a_ncloudp,DIM=4)
        ENDIF
        call writecross_3D(crossvars(n), interp)
      case('nr')
        IF (level<4) THEN
            call writecross_3D(crossvars(n), a_npp)
        ELSE
            interp = SUM(a_nprecpp,DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('rr','r')
        IF (level<4) THEN
            call writecross_3D(crossvars(n), a_rpp)
        ELSE
            interp = SUM(a_mprecpp(:,:,:,1:nprc),DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('ni')
        IF (level==0) THEN
            call writecross_3D(crossvars(n), a_nip)
        ELSEIF (level==5) THEN
            interp = SUM(a_nicep,DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('ri')
        IF (level==0) THEN
            call writecross_3D(crossvars(n), a_rip)
        ELSEIF (level==5) THEN
            interp = SUM(a_micep(:,:,:,1:nice),DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('ns')
        IF (level==0) THEN
            call writecross_3D(crossvars(n), a_nsp)
        ELSEIF (level==5) THEN
            interp = SUM(a_nsnowp,DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('rs')
        IF (level==0) THEN
            call writecross_3D(crossvars(n), a_rsp)
        ELSEIF (level==5) THEN
            interp = SUM(a_msnowp(:,:,:,1:nsnw),DIM=4)
            call writecross_3D(crossvars(n), interp)
        ENDIF
      case('ng')
        IF (level==0) call writecross_3D(crossvars(n), a_ngp)
      case('rg')
        IF (level==0) call writecross_3D(crossvars(n), a_rgp)
      case('nh')
        IF (level==0) call writecross_3D(crossvars(n), a_nhp)
      case('rh')
        IF (level==0) call writecross_3D(crossvars(n), a_rhp)
      CASE('s01','s02','s03','s04','s05','s06','s07','s08','s09','s10') ! etc.
        read (crossvars(n)(2:3),'(i2.2)') i
        interp = a_ap(:,:,:,i)
        call writecross_3D(crossvars(n), interp)
      CASE DEFAULT
        ! Outputs like icNca
        fail = calc_user_data(crossvars(n),interp,mask)
        ! Rates like coag_Na - available for statistics time steps
        IF (fail .AND. sflg) fail = test_if_rate(crossvars(n),interp)
        !
        IF (fail) THEN
            ! Not found
            IF (myid==0) WRITE(*,*) "Variable "//trim(crossvars(n))//" not found!"
        ELSE
            WHERE(.NOT.mask) interp=-999. ! Set masked to -999
            call writecross_3D(crossvars(n), interp)
        ENDIF
      end select
    end do
    !
    if (myid==0) print "(/' ',12('-'),'   Record ',I4,' to cross sections')",nccrossrec
    !
    ! Next
    nccrossrec = nccrossrec +1
  end subroutine triggercross

  LOGICAL FUNCTION test_if_rate(short_name,res)
    USE stat, only : out_mcrp_nout, out_mcrp_list, out_mcrp_data
    use grid, only : nzp, nxp, nyp, maxn_list, out_an_list, out_an_data
    CHARACTER(LEN=7), INTENT(IN) :: short_name ! Variable name
    REAL, INTENT(INOUT) :: res(nzp,nxp,nyp)
    INTEGER :: i
    ! 4D array out_mcrp_data contains 3D data arrays whose names are
    ! specified in the out_mcrp_list containing out_mcrp_nout items
    DO i=1,out_mcrp_nout
        ! Calculate different outputs
        IF (short_name==out_mcrp_list(i)) THEN
            res(:,:,:) = out_mcrp_data(:,:,:,i)
            test_if_rate = .FALSE.
            RETURN
        ENDIF
    ENDDO
    !
    ! 4D analysis data can be available (e.g., forc_xx)
    DO i=1,maxn_list
        IF (short_name==out_an_list(i)) THEN
            ! Analysis data as is
            res(:,:,:) = out_an_data(:,:,:,i)
            test_if_rate = .FALSE.
            RETURN
        ENDIF
    ENDDO
    !
    ! Not found (fail=.TRUE.)
    test_if_rate = .TRUE.
    !
  END FUNCTION test_if_rate

  subroutine writecross_3D(crossname, am)
    use grid, only : nxp, nyp, nzp
    use netcdf, only : nf90_inq_varid, nf90_put_var, nf90_sync
    character(*), intent(in) :: crossname
    real, intent(in) :: am(nzp,nxp,nyp)
    integer :: n, VarID, iret

    ! XZ crosssection
    if (lxz) then
        iret = nf90_inq_varid(nccrossxzid, trim(crossname),VarID)
        iret = nf90_put_var(nccrossxzid, VarID, am(:,3:nxp-2,jcross), start=(/1,1,nccrossrec/))
        iret = nf90_sync(nccrossxzid)
    end if

    ! YZ crosssection
    if (lyz) then
        iret = nf90_inq_varid(nccrossyzid, trim(crossname),VarID)
        iret = nf90_put_var(nccrossyzid, VarID, am(:,icross,3:nyp-2), start=(/1,1,nccrossrec/))
        iret = nf90_sync(nccrossyzid)
    end if

    ! XY crosssections
    if (lxy) then
        iret = nf90_inq_varid(nccrossxyid, trim(crossname),VarID)
        if (nkcross>1) then
            do n=1,nkcross
                iret = nf90_put_var(nccrossxyid, VarID, am(kcross(n),3:nxp-2,3:nyp-2), &
                                    start=(/n,1,1,nccrossrec/), count=(/1,nxp-4,nyp-4,1/))
            end do
        else
            iret = nf90_put_var(nccrossxyid, VarID, am(kcross(1),3:nxp-2,3:nyp-2), &
                                start=(/1,1,nccrossrec/))
        end if
        iret = nf90_sync(nccrossxyid)
    end if

  end subroutine writecross_3D


  subroutine close_cross
    use netcdf, ONLY : nf90_close
    INTEGER :: iret
    if (lxy) iret = nf90_close(nccrossxyid)
    if (lxz) iret = nf90_close(nccrossxzid)
    if (lyz) iret = nf90_close(nccrossyzid)
  end subroutine close_cross

end module modcross
