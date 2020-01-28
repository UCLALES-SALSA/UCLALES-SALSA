!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999, 2001, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module init

  use grid

  integer, parameter    :: nns = 500
  integer               :: ns
  integer               :: iseed = 0
  integer               :: ipsflg = 1
  integer               :: itsflg = 1
       ! itsflg = 0 :potential temperature in kelvin
       !          1 :liquid water potential temperature in kelvin
       !          2 :temperature
  real, dimension(nns)  :: us,vs,ts,thds,ps,hs,rts,rss,tks,xs
  real                  :: zrand = 200.
  real                  :: zrndamp = 0.2 ! the amplitude of random temperature fluctuations
  real                  :: zrndampq = 5.0e-5 ! the amplitude of random humidity fluctuations
  logical               :: zrandnorm = .FALSE. ! normalize the data after inserting random fluctuations
  character  (len=80)   :: hfilin = 'test.'

contains
  !
  ! ----------------------------------------------------------------------
  ! INITLZ:  this is the main driver for the model's initializ-
  ! ation routines.  it initializes the model according to runtype
  !
  subroutine initialize

    use step, only : time, outflg, salsa_diag_update, anl_start, nudging
    use stat, only : init_stat
    use sgsm, only : tkeinit
    use mpi_interface, only : appl_abort, myid
    use thrm, only : thermo
    USE mo_salsa_driver, ONLY : run_SALSA
    USE mo_submctl, ONLY : in2b, fn2b, nlim, stat_b_bins

    implicit none

    ! Local variables for SALSA basic state
    INTEGER :: n4

    if (runtype == 'INITIAL') then
       time=0.
       call arrsnd
       call basic_state
       call fldinit

       ! If SALSA is used, call SALSA with full configuration once before beginning
       ! spin-up period to set up aerosol and cloud fields.
       IF (level >= 4) THEN

          n4 = nspec+1 ! Aerosol components + water

          ! Initialize SALSA species
          CALL aerosol_init
          CALL init_gas_tracers

          ! Update diagnostic SALSA tracers
          CALL SALSA_diag_update
          CALL thermo(level)
          CALL run_SALSA(nxp,nyp,nzp,n4,nbins,ncld,nprc,nice,nsnw,&
                  a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_dn,   &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_gaerop,  a_gaerot,  1, dtl, time, level,   &
                  coag_ra, coag_na, coag_rc, coag_nc, coag_rr, coag_nr, &
                  coag_ri, coag_ni, coag_rs, coag_ns, &
                  cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, &
                  auto_rr, auto_nr, auto_rs, auto_ns, &
                  cact_rc, cact_nc, nucl_ri, nucl_ni, &
                  melt_ri, melt_ni, melt_rs, melt_ns)

          CALL SALSAInit
          CALL SALSA_diag_update
          CALL thermo(level)

       END IF !level >= 4

       ! Initialize nudging
       CALL nudging(time)

    else if (runtype == 'HISTORY') then
       if (isgstyp == 2) call tkeinit(nxyzp,a_qp)
       call hstart
       ! Update diagnostic SALSA tracers
       IF (level >= 4) THEN
          CALL SALSAInit
          CALL SALSA_diag_update
          CALL thermo(level)
       END IF
    else
       if (myid == 0) print *,'  ABORTING:  Invalid Runtype'
       call appl_abort(0)
    end if ! runtype


    ! When SALSA b-bin outputs are needed?
    !   -level >= 4
    !   -outputs are forced (stat_b_bins=.true.) or b-bins are initialized with non-zero concentration
    IF (level >= 4) THEN
        IF (stat_b_bins .AND. no_b_bins) STOP 'Error: b-bins not prognostic, but outputs requested!'
        IF (.not. stat_b_bins .AND. .not.no_b_bins) stat_b_bins=any( a_naerop(:,:,:,in2b:fn2b)>nlim )
    ENDIF

    call sponge_init
    call init_stat(time,filprf,expnme,nzp)
    !
    ! write analysis and history files from restart if appropriate
    !
    if (outflg) then
       if (runtype == 'INITIAL') then
          !call write_hist(1, time)
          call init_anal(time)
          call thermo(level)
          IF (time >= anl_start) call write_anal(time)
       else
          call init_anal(time+dtl)
          call write_hist(0, time)
       end if
    end if !outflg

    return
  end subroutine initialize
  !
  !----------------------------------------------------------------------
  ! FLDINIT: Initializeds 3D fields, mostly from 1D basic state
  !
  !          Modified for level 4.
  !          Juha Tonttila, FMI, 20140828
  !
  subroutine fldinit

    use defs, only : alvl, cpr, cp, p00
    use sgsm, only : tkeinit
    use thrm, only : thermo, rslf

    implicit none

    integer :: i,j,k
    real    :: exner, pres, tk, rc, xran(nzp)

    call htint(ns,ts,hs,nzp,th0,zt)

    do j=1,nyp
       do i=1,nxp
          a_ustar(i,j) = 0.
          do k=1,nzp
             a_up(k,i,j)    = u0(k)
             a_vp(k,i,j)    = v0(k)
             a_tp(k,i,j)    = (th0(k)-th00)
             if (associated (a_rp)) a_rp(k,i,j)   = rt0(k)
             a_theta(k,i,j) = th0(k)
             a_pexnr(k,i,j) = 0.
          end do
       end do
    end do

    ! Juha: Added select-case for level 4
    SELECT CASE(level)
       CASE(1,2,3)
          if ( allocated (a_rv)) a_rv = a_rp

          if ( allocated (a_rc)) then
             do j=1,nyp
                do i=1,nxp
                   do k=1,nzp
                      exner = (pi0(k)+pi1(k))/cp
                      pres  = p00 * (exner)**cpr
                      if (itsflg == 0) then
                         tk    = th0(k)*exner
                         rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                         a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                         a_rv(k,i,j) = a_rp(k,i,j)-rc
                      end if
                      if (itsflg == 2) then
                         tk    = th0(k)
                         a_theta(k,i,j) = tk/exner
                         rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                         a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                         a_rv(k,i,j) = a_rp(k,i,j)-rc
                      end if
                   end do
                end do
             end do
          end if

       CASE(4,5)
          ! Condensation will be calculated by the initial call of SALSA, so use the
          ! saturation adjustment method to estimate the amount of liquid water,
          ! which is needed for theta_l
          DO j = 1,nyp
             DO i = 1,nxp
                DO k = 1,nzp
                   exner = (pi0(k)+pi1(k))/cp
                   pres  = p00 * (exner)**cpr
                   IF (itsflg == 0) THEN
                      tk = th0(k)*exner
                      rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                      a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                   END IF
                   IF (itsflg == 2) THEN
                      tk = th0(k)
                      a_theta(k,i,j) = tk/exner
                      rc  = max(0.,a_rp(k,i,j)-rslf(pres,tk))
                      a_tp(k,i,j) = a_theta(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                   END IF
                END DO !k
             END DO !i
          END DO !j

    END SELECT

    k=1
    do while( zt(k+1) <= zrand .and. k+1 < nzp)
       k=k+1
       xran(k) = zrndamp*(zrand - zt(k))/zrand
    end do
    call random_pert(nzp,nxp,nyp,zt,a_tp,xran,k,'temperature')

    if (associated(a_rp)) then
       k=1
       do while( zt(k+1) <= zrand .and. k+1 < nzp)
          k=k+1
          xran(k) = zrndampq*(zrand - zt(k))/zrand
       end do
       call random_pert(nzp,nxp,nyp,zt,a_rp,xran,k,'humidity')
    end if

    a_wp=0.
    if(isgstyp == 2) call tkeinit(nxyzp,a_qp)
    !
    ! initialize thermodynamic fields
    !
    call thermo (level)

    a_uc=a_up
    a_vc=a_vp
    a_wc=a_wp

    return
  end subroutine fldinit
  !----------------------------------------------------------------------
  ! SPONGE_INIT: Initializes variables for sponge layer
  !
  subroutine sponge_init

    use mpi_interface, only: myid

    implicit none

    integer :: k,kk

    if (nfpt > 0) then
       allocate (spng_tfct(max(1,nfpt)), spng_wfct(max(1,nfpt)))

       do k=nzp-nfpt,nzp-1
          kk = k + 1 - (nzp-nfpt)
          spng_tfct(kk)=max(0.,(zm(nzp)-zt(k))/((zm(nzp)-zm(nzp-nfpt))*distim))
          spng_wfct(kk)=max(0.,(zm(nzp)-zm(k))/((zm(nzp)-zm(nzp-nfpt))*distim))
          spng_tfct(kk) = max(0.,(1./distim - spng_tfct(kk)))
          spng_wfct(kk) = max(0.,(1./distim - spng_wfct(kk)))
       end do

       if(myid == 0) then
          print "(//' ',49('-')/)"
          print '(2X,A17)', 'Sponge Layer Init '
          print '(3X,A12,F6.1,A1)', 'Starting at ', zt(nzp-nfpt), 'm'
          print '(3X,A18,F6.1,A1)', 'Minimum timescale ', 1/spng_wfct(nfpt),'s'
       end if
    end if

    return
  end subroutine sponge_init

  !
  !
  ! ----------------------------------------------------------------------
  ! ARRSND: Arranges the sounding input into proper arrays
  !
  subroutine arrsnd

    use defs, only          : p00,p00i,cp,cpr,rcp,r,g,ep2,alvl,Rm,ep
    use thrm, only          : rslf
    use mpi_interface, only : appl_abort, myid

    implicit none

    integer :: k, iterate
    real    :: tavg, zold2, zold1, x1, xx, yy, zz, til
    character (len=245) :: fm0 = &
         "(/,' -------------------------------------------------',/,"       //&
         "'  Sounding Input: ',//,7x,'ps',9x,'hs',7x,'ts',6x ,'thds',6x," // &
         "'us',7x,'vs',7x,'rts',5x,'rel hum',/,6x,'(Pa)',7X,'(m)',6X,'(K)'"// &
         ",6X,'(K)',6X,'(m/s)',4X,'(m/s)',3X,'(kg/kg)',5X,'(%)',/,1x/)"
    character (len=36) :: fm1 = "(f11.1,f10.1,2f9.2,2f9.2,f10.5,f9.1)"
    !
    ! arrange the input sounding
    !
    if (ps(1) == 0.) then
       open (1,file='sound_in',status='old',form='formatted')
       do ns=1,nns
          read (1,*,end=100) ps(ns),ts(ns),rts(ns),us(ns),vs(ns)
       end do
       close (1)
    end if
100 continue

    ns=1
    do while (ps(ns) /= 0. .and. ns <= nns)
       !
       ! filling relative humidity array only accepts sounding in mixing
       ! ratio (g/kg) converts to (kg/kg)
       !
       rts(ns)=rts(ns)*1.e-3
       !
       ! filling pressure array:
       ! ipsflg = 0 :pressure in millibars
       ! 1 :pressure array is height in meters (ps(1) is surface pressure)
       !
       select case (ipsflg)
       case (0)
          ps(ns)=ps(ns)*100.
       case default
          xs(ns)=(1.+ep2*rts(ns))
          if (ns == 1)then
             ps(ns)=ps(ns)*100.
             zold2=0.
             hs(1) = 0.
          else
             hs(ns) = ps(ns)
             zold1=zold2
             zold2=ps(ns)
             IF ( itsflg==0 .OR. itsflg==1) THEN
                ! ts=potential or liquid water potential temperature (condensation not included here)
                tavg=0.5*(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1))*((p00/ps(ns-1))**rcp)
             ELSE
                ! ts=T [K]
                tavg=0.5*(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1))
             ENDIF
             ps(ns)=(ps(ns-1)**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
          end if
       end select
       !
       ! filling temperature array:
       ! itsflg = 0 :potential temperature in kelvin
       !          1 :liquid water potential temperature in kelvin
       !          2 :temperature
       !
       select case (itsflg)
       case (0)
          tks(ns)=ts(ns)*(ps(ns)*p00i)**rcp
       case (1)
          til=ts(ns)*(ps(ns)*p00i)**rcp
          xx=til
          yy=rslf(ps(ns),xx)
          zz=max(rts(ns)-yy,0.)
          if (zz > 0.) then
             do iterate=1,3
                x1=alvl/(cp*xx)
                xx=xx - (xx - til*(1.+x1*zz))/(1. + x1*til                &
                     *(zz/xx+(1.+yy*ep)*yy*alvl/(Rm*xx*xx)))
                yy=rslf(ps(ns),xx)
                zz=max(rts(ns)-yy,0.)
             enddo
          endif
          tks(ns)=xx
       case (2)
          tks(ns) = ts(ns) ! a long way of saying do nothing
       case default
          if (myid == 0) print *, '  ABORTING: itsflg not supported'
          call appl_abort(0)
       end select
       ns = ns+1
    end do
    ns=ns-1
    !
    ! compute height levels of input sounding.
    !
    if (ipsflg == 0) then
       do k=2,ns
          hs(k)=hs(k-1)-r*.5 *(tks(k)*(1.+ep2*rts(k))                      &
               +tks(k-1)*(1.+ep2*rts(k-1)))*(log(ps(k))-log(ps(k-1)))/g
       end do
    end if

    if (hs(ns) < zt(nzp)) then
       if (myid == 0) print *, '  ABORTING: Model top above sounding top'
       if (myid == 0) print '(2F12.2)', hs(ns), zt(nzp)
       call appl_abort(0)
    end if

    do k=1,ns
       thds(k)=tks(k)*(p00/ps(k))**rcp
    end do

    do k=1,ns
       xs(k)=100.*rts(k)/rslf(ps(k),tks(k))
    end do

    if(myid == 0) then
       write(6,fm0)
       write(6,fm1)(ps(k),hs(k),tks(k),thds(k),us(k),vs(k),rts(k),xs(k),k=1,ns)
    endif

    return
  end subroutine arrsnd
  !
  !----------------------------------------------------------------------
  ! BASIC_STATE: This routine computes the basic state values
  ! of pressure, density, moisture and temperature. The basic state
  ! temperature is assumed to be a the volume weighted average value of
  ! the sounding
  !
  subroutine basic_state

    use defs, only : cp, rcp, cpr, r, g, p00, p00i, ep2
    use mpi_interface, only : myid

    implicit none

    integer k
    real :: v1da(nzp), v1db(nzp), v1dc(nzp), exner

    character (len=305) :: fmt =  &
         "(/,' -------------------------------------------------',/,"     //&
         "'  Basic State: ',//,4X,'Z',6X,'U0',6X,'V0',6X,'DN0',6X,' P0'"   //&
         ",6X,'PRESS',4X,'TH0',6X,'THV',5X,'RT0',/,3X,'(m)',5X,'(m/s)'"     //&
         ",3X,'(m/s)',2X,'(kg/m3)',2X,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X"      //&
         ",'(K)',4X,'(g/kg)',//,(1X,F7.1,2F8.2,F8.3,2F10.2,2F8.2,F7.2))"

    !

    call htint(ns,thds,hs,nzp,th0,zt)
    call htint(ns,us,hs,nzp,u0,zt)
    call htint(ns,vs,hs,nzp,v0,zt)

    if (level >= 1) then
       call htint(ns,rts,hs,nzp,rt0,zt)
       rt0(1)=rt0(2)
    else
       do k=1,nzp
          rt0(k)=0.
       end do
    end if
    !
    ! calculate theta_v for an unsaturated layer, neglecting condensate here is
    ! okay as this is only used for the first estimate of pi1, which will be
    ! updated in a consistent manner on the first dynamic timestep
    !
    do k=1,nzp
       v1dc(k)=th0(k) * (1.+ep2*rt0(k)) ! theta_v assuming unsaturated
    end do
    !
    ! calculate pressure for actual initial state
    !
    pi1(1)=cp*(ps(1)*p00i)**rcp+g*(hs(1)-zt(1))/v1dc(1)
    do k=2,nzp
       pi1(k) = pi1(k-1)-g/(dzm(k-1)*0.5*(v1dc(k)+v1dc(k-1)))
    end do
    !
    ! calculate hydrostatic exner function associated with th00 constant along
    ! with associated basic state density
    !
    pi0(1)=cp*(ps(1)*p00i)**rcp + g*(hs(1)-zt(1))/th00
    dn0(1)=((cp**(1.-cpr))*p00)/(r*th00*pi0(1)**(1.-cpr))
    do k=2,nzp
       pi0(k)=pi0(1) + g*(zt(1) - zt(k))/th00
       dn0(k)=((cp**(1.-cpr))*p00)/(r*th00*pi0(k)**(1.-cpr))
       u0(k)=u0(k)-umean
       v0(k)=v0(k)-vmean
    end do
    !
    ! define pi1 as the difference between pi associated with th0 and pi
    ! associated with th00, thus satisfying pi1+pi0 = pi = cp*(p/p00)**(R/cp)
    !
    do k=1,nzp
       pi1(k) = pi1(k)-pi0(k)
    end do
    !
    do k=1,nzp
       exner = (pi0(k) + pi1(k))/cp
       v1db(k)=p00*(exner)**cpr      ! pressure
       v1da(k)=p00*(pi0(k)/cp)**cpr  ! pressure associated with pi0
    end do

    u0(1) = u0(2)
    v0(1) = v0(2)
    psrf  = ps(1)

    if(myid == 0) write (*,fmt) (zt(k),u0(k),v0(k),dn0(k),v1da(k),v1db(k), &
         th0(k),v1dc(k),rt0(k)*1000.,k=1,nzp)

    return
  end subroutine basic_state
  !
  !---------------------------------------------------------------------
  ! HTINT: Height interpolation of field on one grid, to field on another
  !
  subroutine htint(na,xa,za,nb,xb,zb)

    implicit none
    integer, intent (in) :: na, nb
    real, intent (in)    :: xa(na),za(na),zb(nb)
    real, intent (out)   :: xb(nb)

    integer :: l, k
    real    :: wt

    l = 1
    do k=1,nb
       if (zb(k) <= za(na)) then
          do while ( zb(k) > za(l+1) .and. l < na)
             l=l+1
          end do
          wt=(zb(k)-za(l))/(za(l+1)-za(l))
          xb(k)=xa(l)+(xa(l+1)-xa(l))*wt
       else
          wt=(zb(k)-za(na))/(za(na-1)-za(na))
          xb(k)=xa(na)+(xa(na-1)-xa(na))*wt
       end if
    end do

    return
  end subroutine htint
  !
  ! -----------------------------------------------------------------------
  ! HTINT2d: Same as HTINT but for 2d variables
  !
  SUBROUTINE htint2d(na,xa,za,nb,xb,zb,nx)
    IMPLICIT NONE

    integer, intent (in) :: na, nb, nx
    real, intent (in)    :: xa(na,nx),za(na),zb(nb)
    real, intent (out)   :: xb(nb,nx)

    integer :: l, k, i
    real    :: wt

    do i=1,nx
       l = 1
       do k = 1,nb
          if (zb(k) <= za(na)) then
             do while ( zb(k) > za(l+1) .and. l < na)
                l=l+1
             end do
             wt=(zb(k)-za(l))/(za(l+1)-za(l))
             xb(k,i)=xa(l,i)+(xa(l+1,i)-xa(l,i))*wt
          else
             wt=(zb(k)-za(na))/(za(na-1)-za(na))
             xb(k,i)=xa(na,i)+(xa(na-1,i)-xa(na,i))*wt
          end if
       end do
    end do

  END SUBROUTINE htint2d


  !
  !----------------------------------------------------------------------
  ! HSTART:  This subroutine reads a history file and does
  ! a history start
  !
  subroutine hstart

    use step, only : time
    use mpi_interface, only : myid

    implicit none

    call read_hist(time, hfilin)

    if(myid == 0) &
         print "(//' ',49('-')/,' ',/,' History read from: ',A60)",hfilin

    return
  end subroutine hstart
  !
  !----------------------------------------------------------------------
  ! RANDOM_PERT: initialize field between k=2 and kmx with a
  ! random perturbation of specified magnitude
  !
  subroutine random_pert(n1,n2,n3,zt,fld,xmag,kmx,target_name)

    use mpi_interface, only :  nypg,nxpg,myid,wrxid,wryid,xoffset,yoffset, &
         double_scalar_par_sum

    use util, only : sclrset
    implicit none

    integer, intent(in) :: n1,n2,n3,kmx
    real, intent(inout) :: fld(n1,n2,n3)
    real, intent(in)    :: zt(n1),xmag(n1)
    character(len=*), intent(in) :: target_name

    real (kind=8) :: rand(3:n2-2,3:n3-2),  xx, xxl, tot
    real (kind=8), allocatable :: rand_temp(:,:)

    integer, dimension (:), allocatable :: seed

    integer :: i,j,k,n2g,n3g,isize

    rand=0.0
    ! seed must be a double precision odd whole number greater than
    ! or equal to 1.0 and less than 2**48.

    call random_seed(size=isize)
    allocate (seed(isize))
    seed(:) = iseed
    call random_seed(put=seed)
    deallocate (seed)
    n2g = nxpg
    n3g = nypg

    tot =0.
    do k=2,kmx
       allocate (rand_temp(3:n2g-2,3:n3g-2))
       call random_number(rand_temp)
       rand(3:n2-2, 3:n3-2)=rand_temp(3+xoffset(wrxid):n2+xoffset(wrxid)-2, &
            3+yoffset(wryid):n3+yoffset(wryid)-2)
       deallocate (rand_temp)

       xx = 0.
       do j=3,n3-2
          do i=3,n2-2
             fld(k,i,j) = fld(k,i,j) + rand(i,j)*xmag(k)
             xx = xx + rand(i,j)*xmag(k)
          end do
       end do

       xxl = xx
       call double_scalar_par_sum(xxl,xx)
       xx = xx/real((n2g-4)*(n3g-4))
       IF (zrandnorm) fld(k,:,:)= fld(k,:,:) - xx

       tot = tot + xx/(kmx-1) ! Average perturbation
    end do

    if(myid == 0) then
       print *
       print *,'-------------------------------------------------'
       print *,' Inserting random '//target_name//' perturbations'
       print 600,zt(kmx),rand(3,3),tot
       print *,'-------------------------------------------------'
    endif

    call sclrset('cnst',n1,n2,n3,fld)

    return

600 format( &
         /3x,'Below: ',F7.2,' meters;',                        &
         /3x,'with test value of: ',E12.5,                     &
         /3x,'and a magnitude of: ',E12.5)
  end subroutine random_pert


  !
  !--------------------------------------------------------------------
  ! CLDINIT: Apply the tendencies from the initialization call of SALSA
  !          instantaneously to account for the basic state thermodynamics
  !          and microphysics.
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE SALSAInit
    IMPLICIT NONE
    INTEGER :: k,i,j

    DO j=1,nyp
       DO i=1,nxp
          DO k=1,nzp ! Apply tendencies
             a_naerop(k,i,j,:) = MAX( a_naerop(k,i,j,:)   + dtl*a_naerot(k,i,j,:), 0. )
             a_ncloudp(k,i,j,:) = MAX( a_ncloudp(k,i,j,:) + dtl*a_ncloudt(k,i,j,:), 0. )
             a_nprecpp(k,i,j,:) = MAX( a_nprecpp(k,i,j,:) + dtl*a_nprecpt(k,i,j,:), 0. )
             a_maerop(k,i,j,:)  = MAX( a_maerop(k,i,j,:)  + dtl*a_maerot(k,i,j,:), 0. )
             a_mcloudp(k,i,j,:) = MAX( a_mcloudp(k,i,j,:) + dtl*a_mcloudt(k,i,j,:), 0. )
             a_mprecpp(k,i,j,:) = MAX( a_mprecpp(k,i,j,:) + dtl*a_mprecpt(k,i,j,:), 0. )
             a_gaerop(k,i,j,:)  = MAX( a_gaerop(k,i,j,:)  + dtl*a_gaerot(k,i,j,:), 0. )
             a_rp(k,i,j) = a_rp(k,i,j) + dtl*a_rt(k,i,j)

             IF(level < 5) cycle

             a_nicep(k,i,j,:)   = MAX( a_nicep(k,i,j,:)   + dtl*a_nicet(k,i,j,:), 0. )
             a_nsnowp(k,i,j,:)  = MAX( a_nsnowp(k,i,j,:)  + dtl*a_nsnowt(k,i,j,:), 0. )
             a_micep(k,i,j,:)   = MAX( a_micep(k,i,j,:)   + dtl*a_micet(k,i,j,:), 0. )
             a_msnowp(k,i,j,:)  = MAX( a_msnowp(k,i,j,:)  + dtl*a_msnowt(k,i,j,:), 0. )
          END DO
       END DO
    END DO

  END SUBROUTINE SALSAInit

  ! --------------------------------------------------------------------------------------------------
  ! Replacement for SUBROUTINE init_aero_sizedist (init.f90): initilize altitude-dependent aerosol
  ! size distributions and compositions.
  !
  ! Tomi Raatikainen, FMI, 29.2.2016
  !
  SUBROUTINE aerosol_init

    USE mo_submctl, ONLY : pi6,in1a,in2a,in2b,fn1a,fn2a,fn2b,aerobins, &
                           nmod, sigmag, dpg, n, volDistA, volDistB, nf2a, isdtyp, &
                           iso, rhosu, ioc, rhooc, nspec, dens, zspec, nlim, salsa1a_SO4_OC
    USE mpi_interface, ONLY : myid

    IMPLICIT NONE
    REAL :: core(fn2a), nsect(fn2a)                   ! Size of the bin mid aerosol particle, local aerosol size dist
    REAL :: pndist(nzp,fn2a)                          ! Aerosol size dist as a function of height
    REAL :: pvf2a(nzp,nspec), pvf2b(nzp,nspec)        ! Mass distributions of aerosol species for a and b-bins
    REAL :: pnf2a(nzp)                                ! Number fraction for bins 2a
    REAL :: mass(2*nspec)
    INTEGER :: ss,ee,i,j,k,nc
    CHARACTER(len=600) :: fmt
    LOGICAL :: bbins

    ! Bin mean aerosol particle volume
    core(1:fn2a) = 4.*pi6*(aerobins(1:fn2a)**3+aerobins(2:fn2a+1)**3) ! = 4/3*pi*(rmin**3+rmax**3)/2

    ! Set concentrations to zero
    pndist = 0.
    pvf2a = 0.; pvf2b = 0.
    pnf2a = 0.

    a_maerop(:,:,:,:)=0.0
    a_naerop(:,:,:,:)=0.0

    !
    ! Read  size distributions and compositions
    ! ----------------------------------------------------------
    IF (isdtyp == 1) THEN
       ! Altitude dependent profiles from text or NetCDF files

       CALL READ_AERO_INPUT(pndist,pvf2a,pvf2b,pnf2a)

    ELSE IF (isdtyp == 0) THEN
       ! Uniform profiles based on NAMELIST parameters

       ! Mass fractions for species in a and b-bins
       DO ss = 1,nspec
          pvf2a(:,ss) = volDistA(ss)
          pvf2b(:,ss) = volDistB(ss)
       END DO

       ! Number fraction for 2a
       pnf2a(:) = nf2a
       !
       ! Uniform aerosol size distribution with height.
       ! Using distribution parameters (n, dpg and sigmag) from the SALSA namelist
       !
       ! Convert to SI
       n = n*1.e6
       dpg = dpg*1.e-6
       CALL size_distribution(nmod, n, dpg, sigmag, nsect)
       DO ss = 1,fn2a
          pndist(:,ss) = nsect(ss)
       END DO

    END IF

    !
    ! Initialize concentrations
    ! ----------------------------------------------------------
    DO k = 2,nzp

       DO j = 1,nyp
          DO i = 1,nxp

             ! a) Number concentrations
             ! Region 1
             a_naerop(k,i,j,in1a:fn1a) = pndist(k,in1a:fn1a)

             ! Region 2
             a_naerop(k,i,j,in2a:fn2a) = max(0.0,pnf2a(k))*pndist(k,in2a:fn2a)
             IF (pnf2a(k)<1.) a_naerop(k,i,j,in2b:fn2b) = max(0.0,1.0-pnf2a(k))*pndist(k,in2a:fn2a)

             !
             ! b) Aerosol mass concentrations
             DO nc=1,nspec
                ! 1a and 2a
                ss = nc*nbins + in1a; ee = nc*nbins + fn2a
                a_maerop(k,i,j,ss:ee) = max( 0.0,pvf2a(k,nc) )*pnf2a(k) * &
                     pndist(k,in1a:fn2a)*core(in1a:fn2a)*dens(nc+1)
                ! 2b
                IF (pnf2a(k)<1.) THEN
                   ss = nc*nbins + in2b; ee = nc*nbins + fn2b
                   a_maerop(k,i,j,ss:ee) = max( 0.0,pvf2b(k,nc) )*(1.0-pnf2a(k)) * &
                        pndist(k,in2a:fn2a)*core(in2a:fn2a)*dens(nc+1)
                ENDIF
             END DO

             ! Modify 1a so that there can be sulfate and/or OC?
             IF (salsa1a_SO4_OC) THEN
                ! Set to zero
                a_maerop(k,i,j,in1a:nspec*nbins + fn1a:nbins) = 0.
                IF (iso > 0 .AND. ioc > 0) THEN
                  ! Sulfate and OC mixture
                  ss = (iso-1)*nbins + in1a; ee = (iso-1)*nbins + fn1a
                  a_maerop(k,i,j,ss:ee) = max(0.0,pvf2a(k,iso-1)/(pvf2a(k,ioc-1)+pvf2a(k,iso-1)) )* &
                      pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhosu
                  ss = (ioc-1)*nbins + in1a; ee = (ioc-1)*nbins + fn1a
                  a_maerop(k,i,j,ss:ee) = max(0.0,pvf2a(k,ioc-1)/(pvf2a(k,ioc-1)+pvf2a(k,iso-1)) )* &
                      pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhooc
                ELSEIF (iso > 0) THEN
                  ! Pure sulfate
                  ss = (iso-1)*nbins + in1a; ee = (iso-1)*nbins + fn1a
                  a_maerop(k,i,j,ss:ee) = pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhosu
                ELSEIF (ioc > 0) THEN
                  ! Pure OC
                  ss = (ioc-1)*nbins + in1a; ee = (ioc-1)*nbins + fn1a
                  a_maerop(k,i,j,ss:ee) = pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhooc
                ELSE
                  STOP 'Either OC or SO4 must be active for aerosol region 1a!'
                ENDIF
             ENDIF

             ! Apply concentration threshold
             DO nc = 1,nbins
                IF (a_naerop(k,i,j,nc)*a_dn(k,i,j) < nlim) THEN
                   a_naerop(k,i,j,nc) = 0.
                   a_maerop(k,i,j,nc:nspec*nbins+nc:nbins) = 0.
                END IF
             END DO
          END DO ! i
       END DO ! j
    END DO ! k
    a_naerop(1,:,:,:) = a_naerop(2,:,:,:)
    a_maerop(1,:,:,:) = a_maerop(2,:,:,:)

    ! Put out some info about the initial state
    ! ----------------------------------------------------------
    IF (myid == 0 .AND. isdtyp == 0) THEN
        WRITE(*,'(/,A)') 'Aerosol properties from the NAMELIST (constant size distribution)'

        ! Are b-bins used?
        IF (no_b_bins) THEN
            bbins = .FALSE.
        ELSE
            bbins = SUM(a_naerop(2,3,3,in2b:fn2b))>0.
        ENDIF

        WRITE(*,'(/,A)') ' Initial aerosol number [1e6/kg] and species mass mixing ratios [ug/kg]:'
        ! Header
        WRITE(fmt,"(A12,I2,A8,I2,A4)") "(A7,A12,A11,",nspec,"A12,A11,",nspec,"A12)"
        IF (bbins) THEN
            WRITE(*,fmt) 'Bin','Dmin [m]','Na',zspec(2:nspec+1),'Nb',zspec(2:nspec+1)
        ELSE
            WRITE(*,fmt) 'Bin','Dmin [m]','Na',zspec(2:nspec+1)
        ENDIF
        ! Lines
        WRITE(fmt,"(A17,I2,A13,I2,A7)") "(I7,ES12.2,F11.2,",nspec,"ES12.3,F11.2,",nspec,"ES12.3)"
        DO i=1,fn2a
            IF (i<in2a .OR. .NOT.bbins) THEN
                WRITE(*,FMT) i, aerobins(i)*2., a_naerop(2,3,3,i)*1e-6, a_maerop(2,3,3,i+nbins:i+nspec*nbins:nbins)*1e9
            ELSE
                j = in2b + i - in2a
                WRITE(*,FMT) i, aerobins(i)*2., a_naerop(2,3,3,i)*1e-6, a_maerop(2,3,3,i+nbins:i+nspec*nbins:nbins)*1e9, &
                                                a_naerop(2,3,3,j)*1e-6, a_maerop(2,3,3,j+nbins:j+nspec*nbins:nbins)*1e9
            ENDIF
        ENDDO
        ! Total number and mass
        DO i=1,nspec
            mass(i)=SUM( a_maerop(2,3,3,i*nbins+in1a:i*nbins+fn2a) )*1e9
            IF (bbins) mass(nspec+i)=SUM( a_maerop(2,3,3,i*nbins+in2b:i*nbins+fn2b) )*1e9
        ENDDO
        WRITE(*,'(A)')'    --------------------------'
        WRITE(fmt,"(A11,I2,A13,I2,A7)") "(A19,F11.2,",nspec,"ES12.3,F11.2,",nspec,"ES12.3)"
        IF (bbins) THEN
            WRITE(*,fmt) 'total', SUM(a_naerop(2,3,3,in1a:fn2a))*1.e-6, mass(1:nspec), &
                                  SUM(a_naerop(2,3,3,in2b:fn2b))*1.e-6, mass(nspec+1:2*nspec)
        ELSE
            WRITE(*,fmt) 'total', SUM(a_naerop(2,3,3,in1a:fn2a))*1.e-6, mass(1:nspec)
        ENDIF
    ELSEIF (myid == 0) THEN
        WRITE(*,'(/,A)') 'Aerosol properties from file aerosol_in'

        ! Are b-bins used?
        IF (no_b_bins) THEN
            bbins = .FALSE.
        ELSE
            bbins = SUM(a_naerop(2:nzp-1,3,3,in2b:fn2b))>0.
        ENDIF

        IF (bbins) THEN
            WRITE(*,'(/,A)') ' Initial aerosol profile (total number [1e6/kg] and mass [ug/kg] for a and b bins):'
            ! Header
            WRITE(fmt,"(A9,I2,A8,I2,A4)") "(A12,A10,",nspec,"A12,A10,",nspec,"A12)"
            WRITE(*,fmt) 'Height (m)','Na',zspec(2:nspec+1),'Nb',zspec(2:nspec+1)
            !
            ! Number and mass concentrations for each active species in a and b bins
            WRITE(fmt,"(A13,I2,A13,I2,A7)") "(F12.1,F10.1,",nspec,"ES12.3,F10.1,",nspec,"ES12.3)"
        ELSE
            WRITE(*,'(/,A)') ' Initial aerosol profile (total number [1e6/kg] and mass [ug/kg] for a bins):'
            ! Header
            WRITE(fmt,"(A9,I2,A4)") "(A12,A10,",nspec,"A12)"
            WRITE(*,fmt) 'Height (m)','Na',zspec(2:nspec+1)
            !
            ! Number and mass concentrations for each active species in a bins
            WRITE(fmt,"(A13,I2,A7)") "(F12.1,F10.1,",nspec,"ES12.3)"
        ENDIF
        !
        DO k=1,nzp
            ! Calculate mass (1e-6g/kg)
            DO i=1,nspec
                mass(i)=SUM( a_maerop(k,3,3,i*nbins+in1a:i*nbins+fn2a) )*1e9
                IF (bbins) mass(nspec+i)=SUM( a_maerop(k,3,3,i*nbins+in2b:i*nbins+fn2b) )*1e9
            ENDDO
            ! Print
            IF (bbins)  THEN
                WRITE(*,fmt) zt(k), SUM(a_naerop(k,3,3,in1a:fn2a))*1.e-6, mass(1:nspec), &
                                SUM(a_naerop(k,3,3,in2b:fn2b))*1.e-6, mass(nspec+1:2*nspec)
            ELSE
                WRITE(*,fmt) zt(k), SUM(a_naerop(k,3,3,in1a:fn2a))*1.e-6, mass(1:nspec)
            ENDIF
        ENDDO
    ENDIF
  END SUBROUTINE aerosol_init
  !
  ! -------------------------------------------------------------------------
  ! Reads vertical profiles of aerosol size distribution parameters, aerosol species volume fractions and
  ! number concentration fractions between a and b bins
  !
  SUBROUTINE READ_AERO_INPUT(ppndist,ppvf2a,ppvf2b,ppnf2a)
    USE ncio, ONLY : open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_aero_nc
    USE mo_submctl, ONLY : fn2a, nspec, maxspec, nmod
    USE mpi_interface, ONLY : appl_abort, myid
    IMPLICIT NONE

    REAL, INTENT(out) :: ppndist(nzp,fn2a)                    ! Aerosol size dist as a function of height
    REAL, INTENT(out) :: ppvf2a(nzp,nspec), ppvf2b(nzp,nspec) ! Volume distributions of aerosol species for a and b-bins
    REAL, INTENT(out) :: ppnf2a(nzp)                          ! Number fraction for bins 2a

    REAL :: nsect(fn2a)

    INTEGER :: ncid, k, i
    INTEGER :: nc_levs=500, nc_nspec, nc_nmod

    ! Stuff that will be read from the file
    REAL, ALLOCATABLE :: zlevs(:),        &  ! Levels in meters
                         zvolDistA(:,:),  &  ! Volume distribution of aerosol species in a
                         zvoldistB(:,:),  &  ! ... and b bins
                         znf2a(:),        &  ! Number fraction for bins 2a
                         zn(:,:),         &  ! Aerosol mode number concentrations
                         zsigmag(:,:),    &  ! Geometric standard deviations
                         zdpg(:,:),       &  ! Mode mean diameters
                         znsect(:,:)         ! Helper for binned number concentrations
    LOGICAL :: READ_NC

    ! Read the NetCDF input when it is available
    INQUIRE(FILE='aerosol_in.nc',EXIST=READ_NC)

    ! Allocate input variables
    ALLOCATE( zlevs(nc_levs),              &
              zvolDistA(nc_levs,maxspec),  &
              zvolDistB(nc_levs,maxspec),  &
              znf2a(nc_levs),              &
              zn(nc_levs,nmod),            &
              zsigmag(nc_levs,nmod),       &
              zdpg(nc_levs,nmod)           )
    zlevs = 0.; zvolDistA = 0.; zvolDistB = 0.; znf2a = 0.; zn = 0.; zsigmag = 0.; zdpg = 0.

    IF (READ_NC) THEN
       ! Open the input file
       CALL open_aero_nc(ncid, nc_levs, nc_nspec, nc_nmod)

       ! Read the aerosol profile data
       CALL read_aero_nc_1d(ncid,'levs',nc_levs,zlevs)
       CALL read_aero_nc_2d(ncid,'volDistA',nc_levs,maxspec,zvolDistA)
       CALL read_aero_nc_2d(ncid,'volDistB',nc_levs,maxspec,zvolDistB)
       CALL read_aero_nc_1d(ncid,'nf2a',nc_levs,znf2a)
       CALL read_aero_nc_2d(ncid,'n',nc_levs,nmod,zn)
       CALL read_aero_nc_2d(ncid,'dpg',nc_levs,nmod,zdpg)
       CALL read_aero_nc_2d(ncid,'sigmag',nc_levs,nmod,zsigmag)

       CALL close_aero_nc(ncid)
    ELSE
       ! Read the profile data from a text file
       open (11,file='aerosol_in',status='old',form='formatted')
       do i=1,nc_levs
          read (11,*,end=100) zlevs(i)
          read (11,*,end=100) (zvolDistA(i,k),k=1,nspec) ! Note: reads just "nspec" values from the current line
          read (11,*,end=100) (zvolDistB(i,k),k=1,nspec) ! -||-
          read (11,*,end=100) (zn(i,k),k=1,nmod)
          read (11,*,end=100) (zdpg(i,k),k=1,nmod)
          read (11,*,end=100) (zsigmag(i,k),k=1,nmod)
          read (11,*,end=100) znf2a(i)
       end do
100    continue
       close (11)
       !
       ! The true number of altitude levels
       nc_levs=i-1
    END IF
    !
    IF (zlevs(nc_levs)<zt(nzp)) then
       if (myid == 0) print *, '  ABORTING: Model top above aerosol sounding top'
       if (myid == 0) print '(2F12.2)', zlevs(nc_levs), zt(nzp)
       call appl_abort(0)
    END IF

    ! Convert to SI
    zn = zn*1.e6
    zdpg = zdpg*1.e-6

    ! Get the binned size distribution
    ALLOCATE( znsect(nc_levs,fn2a) )
    znsect = 0.
    DO k = 1,nc_levs
       CALL size_distribution(nmod,zn(k,:),zdpg(k,:),zsigmag(k,:),nsect)
       znsect(k,:) = nsect(:)
    END DO

    ! Interpolate the input variables to model levels
    ! ------------------------------------------------
    CALL htint2d(nc_levs,zvolDistA(1:nc_levs,1:nspec),zlevs(1:nc_levs),nzp,ppvf2a,zt,nspec)
    CALL htint2d(nc_levs,zvolDistB(1:nc_levs,1:nspec),zlevs(1:nc_levs),nzp,ppvf2b,zt,nspec)
    CALL htint2d(nc_levs,znsect(1:nc_levs,:),zlevs(1:nc_levs),nzp,ppndist,zt,fn2a)
    CALL htint(nc_levs,znf2a(1:nc_levs),zlevs(1:nc_levs),nzp,ppnf2a,zt)

    DEALLOCATE( zlevs, zvolDistA, zvolDistB, znf2a, zn, zsigmag, zdpg, znsect )

  END SUBROUTINE READ_AERO_INPUT

  !
  !------------------------------------------------------------------
  ! INIT_GAS_TRACERS: Set initial values for gas compound tracers
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE init_gas_tracers
    USE mpi_interface, ONLY : myid
    USE mo_submctl, ONLY : avog, mws_gas, zgas, &
        part_h2so4, conc_h2so4, part_ocnv, conc_ocnv, &
        ox_prescribed, conc_oh, conc_o3, conc_no3, mair, &
        nvocs, nvbs, naqsoa, conc_voc, conc_vbsg, conc_aqsoag, &
        ngases_diag, zgas_diag, model_lat, start_doy
    USE step, ONLY : cntlat, strtim
    IMPLICIT NONE

    ! Local variables
    INTEGER :: i, j, k, iout
    CHARACTER(LEN=300) :: fmt
    REAL :: array(nzp,10)

    ! Gases initialized in SALSA
    IF (ngases==0) RETURN

    ! Note: gas phase concentration units for prognostic and diagnostic variables
    ! are kg/kg and mol/kg, respectively. Inputs can have different units!

    ! Sulfate and non-volatile organics
    !  - Input as molecules/kg
    i=0
    IF (part_h2so4) THEN
        ! Sulfate or H2SO4
        i=i+1
        a_gaerop(:,:,:,i) = conc_h2so4/avog*mws_gas(i)
    ENDIF
    IF (part_ocnv) THEN
        ! Non-volatile organic vapor
        i=i+1
        a_gaerop(:,:,:,i) = conc_ocnv/avog*mws_gas(i)
    ENDIF


    ! VBS: oxidants (OH, O3 and NOx) + VOCs(g) + VBS(g) [+ aqSOA(g)]
    !   a) Oxidants: initial concentration given as a number mixing ratio
    j=0
    IF (conc_oh>=0.) THEN
        IF (ox_prescribed) THEN
            ! Diagnostic (constant)
            j=j+1
            zgas_diag(j) = conc_oh/mair ! Note: mol/kg
        ELSE
            ! Prognostic (variable)
            i=i+1
            a_gaerop(:,:,:,i) = conc_oh*mws_gas(i)/mair ! kg/kg
        ENDIF
    ENDIF
    IF (conc_o3>0.) THEN
        IF (ox_prescribed) THEN
            j=j+1
            zgas_diag(j) = conc_o3/mair
        ELSE
            i=i+1
            a_gaerop(:,:,:,i) = conc_o3*mws_gas(i)/mair
        ENDIF
    ENDIF
    IF (conc_no3>0.) THEN
        IF (ox_prescribed) THEN
            j=j+1
            zgas_diag(j) = conc_o3/mair
        ELSE
            i=i+1
            a_gaerop(:,:,:,i) = conc_no3*mws_gas(i)/mair
        ENDIF
    ENDIF
    !   b) VOCs: initial concentration given as a mass mixing ratio (the same for VBS and aqSOA)
    !   See if the input profile data file vocg_in exists
    IF (nvocs>0) CALL read_input_array('vocg_in',nzp,zt,nvocs,array(:,1:nvocs),iout)
    DO j=1,nvocs
        i=i+1
        IF (iout==0) THEN
            a_gaerop(:,:,:,i) = conc_voc(j)
        ELSE
            DO k=1,nzp
                a_gaerop(k,:,:,i) = array(k,j)
            ENDDO
        ENDIF
    ENDDO
    !   c) VBS bins
    IF (nvbs>0) CALL read_input_array('vbsg_in',nzp,zt,nvbs,array(:,1:nvbs),iout)
    DO j=1,nvbs
        i=i+1
        IF (iout==0) THEN
            a_gaerop(:,:,:,i) = conc_vbsg(j)
        ELSE
            DO k=1,nzp
                a_gaerop(k,:,:,i) = array(k,j)
            ENDDO
        ENDIF
    ENDDO
    !   d) aqSOA
    IF (naqsoa>0) CALL read_input_array('aqsoag_in',nzp,zt,naqsoa,array(:,1:naqsoa),iout)
    DO j=1,naqsoa
        i=i+1
        IF (iout==0) THEN
            a_gaerop(:,:,:,i) = conc_aqsoag(j)
        ELSE
            DO k=1,nzp
                a_gaerop(k,:,:,i) = array(k,j)
            ENDDO
        ENDIF
    ENDDO

    ! Additional VBS parameters
    start_doy=strtim ! Start time as decimal day of year
    model_lat=cntlat ! Center latitude (degrees)

    ! Info
    IF (myid == 0 .AND. ngases>0) THEN
        WRITE(*,*) ''
        WRITE(*,'(/,A)') ' Initial gas/vapor concentration profile [ug/kg]:'
        ! Header
        WRITE(fmt,"(A9,I2,A8)") "(A12,",ngases,"A12)"
        WRITE(*,fmt) 'Height (m)',zgas(1:ngases)
        !
        ! Number and mass concentrations for each active species in a and b bins
        WRITE(fmt,"(A7,I2,A7)") "(F12.1,",ngases,"ES12.3)"
        !
        DO i=1,nzp
            ! Print
            WRITE(*,fmt) zt(i), a_gaerop(i,3,3,:)*1.e9 ! kg => ug
            IF (i==5 .AND. ALL( ABS(a_gaerop(4,3,3,:)-a_gaerop(5,3,3,:))<1e-13 )) THEN
                WRITE(*,'(A14)')'...'
                EXIT
            ENDIF
        ENDDO
        !
        ! Diagnostic
        IF (ngases_diag>0) THEN
            WRITE(*,*) ''
            WRITE(*,*) 'Diagnostic gas phase species (id, name and concentration [mol/kg]):'
            DO i=1,ngases_diag
                WRITE(*,"(4X,I2,2X,A3,2X,ES12.3)")i,zgas(ngases+i),zgas_diag(i)
            ENDDO
        ENDIF
    ENDIF

  END SUBROUTINE init_gas_tracers

  !
  !------------------------------------------------------------------
  ! read_input_array: read an input text file where the first column is altitude (m) and
  ! the remaining columns contain the data. The data is interpolated for the current grid.
  !
  ! Tomi Raatikainen, FMI, 2020
  !
  SUBROUTINE read_input_array(fname,nzp,zt,ncols,output,istat)
    ! Read standard input text files and interpolate these to the current vertical grid
    USE mpi_interface, ONLY : appl_abort, myid
    IMPLICIT NONE
    ! Inputs and outputs
    CHARACTER(LEN=*), INTENT(IN) :: fname ! Data file name
    INTEGER, INTENT(IN) :: nzp,ncols ! Known output dimensions
    REAL, INTENT(IN) :: zt(nzp) ! Target or output grid
    REAL, INTENT(OUT) :: output(nzp,ncols) ! Output data
    INTEGER, INTENT(OUT) :: istat ! Flag: 1 = OK, 0 = file not found, which can be OK
    ! Local variables
    INTEGER, PARAMETER :: max_levs=1000 ! Maximum number or rows
    REAL :: zlevs(max_levs), ztmp(max_levs,ncols) ! z and data
    INTEGER :: k, i, nc_levs
    LOGICAL :: read_file
    !
    ! File may not exists
    istat = 0
    INQUIRE(FILE=fname,EXIST=read_file)
    IF (.NOT. read_file) RETURN
    !
    ! File exists, so read the data
    open (11,file=fname,status='old',form='formatted')
    do i=1,max_levs
        ! Each row contains altitude and then data (ncols columns)
        read (11,*,end=100) zlevs(i), (ztmp(i,k),k=1,ncols)
    end do
100   continue
    close (11)
    ! The true number of altitude levels
    nc_levs=i-1
    !
    IF (nc_levs<1) THEN
       if (myid == 0) print *, '  ABORTING: empty input file '//TRIM(fname)
       call appl_abort(0)
    ELSEIF (zlevs(nc_levs)<zt(nzp)) then
       if (myid == 0) print *, '  ABORTING: Model top above input top in file '//TRIM(fname)
       if (myid == 0) print '(2F12.2)', zlevs(nc_levs), zt(nzp)
       call appl_abort(0)
    END IF
    !
    ! Interpolate input variables to model levels
    CALL htint2d(nc_levs,ztmp(1:nc_levs,:),zlevs(1:nc_levs),nzp,output,zt,ncols)
    !
    ! All OK
    istat = 1
  END SUBROUTINE read_input_array


  !
  !------------------------------------------------------------------
  ! SIZE_DISTRIBUTION: Converts multimodal log-normal size distribution into size bins
  !
  ! From module mo_salsa_sizedist
  !
  SUBROUTINE size_distribution(nmod, n, dpg, sigmag, naero)
    USE mo_submctl, ONLY : pi, fn2a, aerobins
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: nmod ! number of modes
    REAL, INTENT(IN) :: n(nmod), dpg(nmod), sigmag(nmod) ! Mode number concentration, mean diameter and width
    REAL, INTENT(OUT) :: naero(fn2a) ! number concentration

    !-- local variables
    REAL :: deltadp,d1,d2,delta_d,dmid
    INTEGER :: kk,ib

    naero(:)=0.
    DO kk = 1, fn2a ! Bin
        d1 = aerobins(kk)*2.
        d2 = aerobins(kk+1)*2.
        delta_d=(d2-d1)/10.
        DO ib = 1,10
            d1=aerobins(kk)*2.+(ib-1.)*delta_d
            d2=d1+delta_d
            dmid=(d1+d2)/2.
            deltadp = log(d2/d1)
            !-- size distribution
            !   ntot = total number, total area, or total volume concentration
            !   dpg = geometric-mean number, area, or volume diameter
            !   n(kk) = number, area, or volume concentration in a bin
            naero(kk) = naero(kk)+sum(n*deltadp/        &
                 (sqrt(2.*pi)*log(sigmag))*                   &
                 exp(-log(dmid/dpg)**2/(2.*log(sigmag)**2)))
        END DO
    END DO
  END SUBROUTINE size_distribution

end module init
