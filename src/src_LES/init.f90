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
  character  (len=80)   :: hfilin = 'test.'

contains
  !
  ! ----------------------------------------------------------------------
  ! INITLZ:  this is the main driver for the model's initializ-
  ! ation routines.  it initializes the model according to runtype
  !
  subroutine initialize

    use step, only : time, outflg
    use stat, only : init_stat, mcflg, acc_massbudged
    use sgsm, only : tkeinit
    use mpi_interface, only : appl_abort, myid
    use thrm, only : thermo
    USE mo_salsa_driver, ONLY : run_SALSA
    USE mo_submctl, ONLY : nbins ! Olis parempi jos ei tarttis
    USE util, ONLY : maskactiv
    USE class_ComponentIndex, ONLY : GetNcomp

    implicit none

    ! Local variables for SALSA basic state
    REAL :: zwp(nzp,nxp,nyp), ztkt(nzp,nxp,nyp)
    LOGICAL :: zactmask(nzp,nxp,nyp)
    LOGICAL :: TMP
    INTEGER :: n4
    
    ztkt = 0.

    TMP = .false.

    ! Set vertical velocity as 0.5 m/s to intialize cloud microphysical properties with
    ! SALSA
    zwp(:,:,:) = 0.5

    if (runtype == 'INITIAL') then
       time=0.
       call arrsnd
       call basic_state
       call fldinit ! Juha: aerosol size distributions are initialized here.
                    !       Also thermodynamics!

       ! If SALSA is used, call SALSA with full configuration once before beginning
       ! spin-up period to set up aerosol and cloud fields.
       IF (level >= 4) THEN

          ! Not needed when using interst. acivation?
          CALL maskactiv(zactmask,nxp,nyp,nzp,nbins,1,prtcl,a_rh)

          n4 = GetNcomp(prtcl) + 1 ! Aerosol compoenents + water

          IF ( nxp == 5 .and. nyp == 5 ) THEN
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_scr1,ztkt,a_rp,a_rt,a_scr2,a_rsi,zwp,a_dn, &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  a_Radry,   a_Rcdry,   a_Rpdry,               &
                  a_Ridry,   a_Rsdry,                          &
                  a_Rawet,   a_Rcwet,   a_Rpwet,               &
                  a_Riwet,   a_Rswet,                          &
                  1, prtcl, dtlt, .false., 0., level   )
          ELSE
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_scr1,ztkt,a_rp,a_rt,a_scr2,a_rsi,a_wp,a_dn, &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  a_Radry,   a_Rcdry,   a_Rpdry,               &
                  a_Ridry,   a_Rsdry,                          &
                  a_Rawet,   a_Rcwet,   a_Rpwet,               &
                  a_Riwet,   a_Rswet,                          &
                  1, prtcl, dtlt, .false., 0., level   )

          END IF
          CALL SALSAInit

          
       END IF !level >= 4

    else if (runtype == 'HISTORY') then
       if (isgstyp == 2) call tkeinit(nxyzp,a_qp)
       call hstart
    else
       if (myid == 0) print *,'  ABORTING:  Invalid Runtype'
       call appl_abort(0)
    end if ! runtype

    call sponge_init
    call init_stat(time+dtl,filprf,expnme,nzp)
    !
    IF (mcflg) THEN
       ! Juha:
       ! Calculate some numbers for mass concervation experiments
       mc_Vdom = deltax*deltay*deltaz*(nxp-4)*(nyp-4)*(nzp-1)
       mc_Adom = deltax*deltay*(nxp-4)*(nyp-4)
       mc_ApVdom = mc_Adom/mc_Vdom
       ! Get the initial mass of atmospheric water
       CALL acc_massbudged(nzp,nxp,nyp,0,dtlt,dzt,a_dn,     &
            rv=a_rp,rc=a_rc,prc=a_srp)
    END IF ! mcflg
    !
    ! write analysis and history files from restart if appropriate
    !
    if (outflg) then
       if (runtype == 'INITIAL') then
          call write_hist(1, time)
          call init_anal(time)
          call thermo(level)
          call write_anal(time)
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

    use defs, only : alvl, cpr, cp, p00, R
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
          ! Condensation will be calculated by the initial call of SALSA
          ! Assume no condensation at this point
          DO j = 1,nyp
             DO i = 1,nxp
                DO k = 1,nzp
                   exner = (pi0(k)+pi1(k))/cp
                   pres  = p00 * (exner)**cpr
                   IF (itsflg == 0) THEN
                      tk = th0(k)*exner
                      a_tp(k,i,j) = a_theta(k,i,j) - th00
                   END IF
                   IF (itsflg == 2) THEN
                      tk = th0(k)
                      a_theta(k,i,j) = tk/exner
                      a_tp(k,i,j) = a_theta(k,i,j) - th00
                   END IF
                END DO !k
             END DO !i
          END DO !j

    END SELECT

    k=1
    do while( zt(k+1) <= zrand .and. k < nzp)
       k=k+1
       xran(k) = 0.2*(zrand - zt(k))/zrand
    end do
    call random_pert(nzp,nxp,nyp,zt,a_tp,xran,k)

    if (associated(a_rp)) then
       k=1
       do while( zt(k+1) <= zrand .and. k < nzp)
          k=k+1
          xran(k) = 5.0e-5*(zrand - zt(k))/zrand
       end do
       call random_pert(nzp,nxp,nyp,zt,a_rp,xran,k)
    end if

    a_wp=0.
    if(isgstyp == 2) call tkeinit(nxyzp,a_qp)
    !
    ! initialize thermodynamic fields
    !
    call thermo (level)

    !
    ! Initialize aerosol size distributions
    !
    IF (level >= 4) THEN
       CALL aerosol_init
       CALL liq_ice_init      ! This should be replaced by physical processing!
       CALL init_gas_tracers
    END IF

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

    zold1 = 0.
    zold2 = 0.

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
             tavg=(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1)*(p00**rcp)             &
                  /ps(ns-1)**rcp)*.5
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
  ! of pressure, density, moisture and temperature.  The basi!state
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

    dtlv=2.*dtl
    dtlt=dtl

    if(myid == 0) &
         print "(//' ',49('-')/,' ',/,' History read from: ',A60)",hfilin

    return
  end subroutine hstart
  !
  !----------------------------------------------------------------------
  ! RANDOM_PERT: initialize field between k=2 and kmx with a
  ! random perturbation of specified magnitude
  !
  subroutine random_pert(n1,n2,n3,zt,fld,xmag,kmx)

    use mpi_interface, only :  nypg,nxpg,myid,wrxid,wryid,xoffset,yoffset, &
         double_scalar_par_sum

    use util, only : sclrset
    implicit none

    integer, intent(in) :: n1,n2,n3,kmx
    real, intent(inout) :: fld(n1,n2,n3)
    real, intent(in)    :: zt(n1),xmag(n1)

    real (kind=8) :: rand(3:n2-2,3:n3-2),  xx, xxl
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
          end do
       end do

       xxl = xx
       call double_scalar_par_sum(xxl,xx)
       xx = xx/real((n2g-4)*(n3g-4))
       fld(k,:,:)= fld(k,:,:) - xx
    end do

    if(myid == 0) then
       print *
       print *,'-------------------------------------------------'
       print 600,zt(kmx),rand(3,3),xx
       print *,'-------------------------------------------------'
    endif

    call sclrset('cnst',n1,n2,n3,fld)

    return

600 format(2x,'Inserting random temperature perturbations',    &
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
    USE mo_submctl, ONLY : ncld,nbins,nice
    IMPLICIT NONE

    INTEGER :: k,i,j,bb

    DO j=1,nyp
       DO i=1,nxp
          DO k=1,nzp ! Apply tendencies
             a_naerop(k,i,j,:) = MAX( a_naerop(k,i,j,:) + dtlt*a_naerot(k,i,j,:), 0. )
             a_ncloudp(k,i,j,:) = MAX( a_ncloudp(k,i,j,:) + dtlt*a_ncloudt(k,i,j,:), 0. )
             a_nprecpp(k,i,j,:) = MAX( a_nprecpp(k,i,j,:) + dtlt*a_nprecpt(k,i,j,:), 0. )
             a_maerop(k,i,j,:)  = MAX( a_maerop(k,i,j,:)  + dtlt*a_maerot(k,i,j,:), 0. )
             a_mcloudp(k,i,j,:) = MAX( a_mcloudp(k,i,j,:) + dtlt*a_mcloudt(k,i,j,:), 0. )
             a_mprecpp(k,i,j,:) = MAX( a_mprecpp(k,i,j,:) + dtlt*a_mprecpt(k,i,j,:), 0. )
             a_gaerop(k,i,j,:)  = MAX( a_gaerop(k,i,j,:)  + dtlt*a_gaerot(k,i,j,:), 0. )
             a_rp(k,i,j) = a_rp(k,i,j) + dtlt*a_rt(k,i,j)

            IF(level < 5) cycle

             a_nicep(k,i,j,:)   = MAX( a_nicep(k,i,j,:)   + dtlt*a_nicet(k,i,j,:), 0. )
             a_nsnowp(k,i,j,:)  = MAX( a_nsnowp(k,i,j,:)  + dtlt*a_nsnowt(k,i,j,:), 0. )
             a_micep(k,i,j,:)   = MAX( a_micep(k,i,j,:)   + dtlt*a_micet(k,i,j,:), 0. )
             a_msnowp(k,i,j,:)  = MAX( a_msnowp(k,i,j,:)  + dtlt*a_msnowt(k,i,j,:), 0. )
          END DO
       END DO
    END DO

    ! Activation + diagnostic array initialization
    ! Clouds and aerosols
    a_rc(:,:,:) = 0.
    DO bb = 1, ncld
       CALL DiagInitCloud(bb)
    END DO
    DO bb = 1,nbins
       CALL DiagInitAero(bb)
    END DO

    ! Ice
    a_ri(:,:,:) = 0.
    do bb = 1,nice
        call DiagInitIce(bb)
    end do

  END SUBROUTINE SALSAInit
  !-------------------------------------------
  SUBROUTINE ActInit(b,bpar,pactmask)
    USE mo_submctl, ONLY : ncld,nbins
    USE class_componentIndex, ONLY : GetNcomp, GetIndex
    IMPLICIT NONE

    INTEGER, INTENT(in) :: b, bpar
    LOGICAL, INTENT(in) :: pactmask(nzp,nxp,nyp)
    REAL :: frac(nzp,nxp,nyp)
    INTEGER :: m,mpar,nc,s

    ! --------------------
    ! Initial activation
    ! --------------------
    ! Get the number of activated droplets
    a_ncloudp(:,:,:,b) = MERGE(a_nactd(:,:,:,b), 0., pactmask(:,:,:))

    DO s = 1,GetNcomp(prtcl) ! Dry aerosol mass
       m = (s-1)*ncld + b
       a_mcloudp(:,:,:,m) = MERGE(a_vactd(:,:,:,m), 0., pactmask(:,:,:))
    END DO

    nc = GetIndex(prtcl,'H2O')
    m = (nc-1)*ncld + b
    mpar = (nc-1)*nbins + bpar
    ! Water mass
    a_mcloudp(:,:,:,m) = MERGE(a_vactd(:,:,:,m), 0., pactmask(:,:,:))

    frac(:,:,:) = a_nactd(:,:,:,b)/MAX(a_ncloudp(:,:,:,b),1.)

    ! Remove water from the vapor phase
    a_rp(:,:,:) = a_rp(:,:,:) -    &
         MERGE( a_vactd(:,:,:,m)-frac(:,:,:)*a_maerop(:,:,:,mpar), 0.,  &
         pactmask(:,:,:))

    ! Remove aerosol particles due to activation
    a_naerop(:,:,:,bpar) = a_naerop(:,:,:,bpar) -   &
         MERGE(a_nactd(:,:,:,b), 0., pactmask(:,:,:))

    ! Water from aerosol phase
    a_maerop(:,:,:,mpar) = a_maerop(:,:,:,mpar) -   &
         MERGE(frac(:,:,:)*a_maerop(:,:,:,mpar), 0., pactmask(:,:,:))

    ! Aerosol mass
    DO s = 1,GetNcomp(prtcl)
       m = (s-1)*ncld + b
       mpar = (s-1)*nbins + bpar
       a_maerop(:,:,:,mpar) = a_maerop(:,:,:,mpar) -  &
            MERGE(a_vactd(:,:,:,m), 0., pactmask(:,:,:))
    END DO

  END SUBROUTINE ActInit
  !------------------------------------------
  SUBROUTINE DiagInitCloud(b)
    USE mo_submctl, ONLY : ncld, pi6, rhowa, rhosu, nlim
    USE class_componentIndex, ONLY : GetIndex
    IMPLICIT NONE

    INTEGER, INTENT(in) :: b

    INTEGER :: str,k,i,j,nc
    REAL :: zvol

    nc = GetIndex(prtcl,'H2O')
    str = (nc-1)*ncld+b

       DO j = 1,nyp
          DO i = 1,nxp
             DO k = 1,nzp

                IF (a_ncloudp(k,i,j,b)  > nlim) THEN
                   CALL binMixrat('cloud','dry',b,i,j,k,zvol)
                   zvol = zvol/rhosu
                   a_Rcdry(k,i,j,b) = 0.5*( zvol/(pi6*a_ncloudp(k,i,j,b)) )**(1./3.)
                   CALL binMixrat('cloud','wet',b,i,j,k,zvol)
                   zvol = zvol/rhowa
                   a_Rcwet(k,i,j,b) = 0.5*( zvol/(pi6*a_ncloudp(k,i,j,b)) )**(1./3.)
                ELSE
                   a_Rcdry(k,i,j,b) = 1.e-10
                   a_Rcwet(k,i,j,b) = 1.e-10
                END IF

                ! Cloud water
                a_rc(k,i,j) = a_rc(k,i,j) + a_mcloudp(k,i,j,str)

             END DO ! k
          END DO ! i
       END DO ! j
  END SUBROUTINE DiagInitCloud
    !------------------------------------------
  SUBROUTINE DiagInitIce(b)
    USE mo_submctl, ONLY : nice, pi6, rhosu,rhoic, prlim
    USE class_componentIndex, ONLY : IsUsed, GetIndex
    IMPLICIT NONE

    INTEGER, INTENT(in) :: b

    INTEGER :: str,k,i,j,nc
    REAL :: zvol

    nc = GetIndex(prtcl,'H2O')
    str = (nc-1)*nice + b

       DO j = 1,nyp
          DO i = 1,nxp
             DO k = 1,nzp

                IF (a_nicep(k,i,j,b)  > prlim) THEN
                   CALL binMixrat('ice','dry',b,i,j,k,zvol)
                    zvol = zvol/rhosu
                   a_Ridry(k,i,j,b) = 0.5*( zvol/(pi6*a_nicep(k,i,j,b)) )**(1./3.)
                   CALL binMixrat('ice','wet',b,i,j,k,zvol)
                    zvol = zvol/rhoic
                   a_Riwet(k,i,j,b) = 0.5*( zvol/(pi6*a_nicep(k,i,j,b)) )**(1./3.)
                ELSE
                   a_Ridry(k,i,j,b) = 1.e-10
                   a_Riwet(k,i,j,b) = 1.e-10
                END IF

                ! Cloud ice
                a_ri(k,i,j) = a_ri(k,i,j) + a_micep(k,i,j,str)

             END DO ! k
          END DO ! i
       END DO ! j
  END SUBROUTINE DiagInitIce
    !------------------------------------------
  SUBROUTINE DiagInitAero(b)
    USE mo_submctl, ONLY : nbins, pi6, rhowa, nlim
    USE class_componentIndex, ONLY : IsUsed, GetIndex
    IMPLICIT NONE

    INTEGER, INTENT(in) :: b

    INTEGER :: k,i,j,nc, str
    REAL :: zvol

    nc = GetIndex(prtcl,'H2O')
    str = (nc-1)*nbins+b

       DO j = 1,nyp
          DO i = 1,nxp
             DO k = 1,nzp

                IF (a_naerop(k,i,j,b) > nlim) THEN
                   CALL binMixrat('aerosol','dry',b,i,j,k,zvol)
                   a_Radry(k,i,j,b) = 0.5*( zvol/(pi6*a_naerop(k,i,j,b)) )**(1./3.)
                   CALL binMixrat('aerosol','wet',b,i,j,k,zvol)
                   a_Rawet(k,i,j,b) = 0.5*( zvol/(pi6*a_naerop(k,i,j,b)) )**(1./3.)
                ELSE
                   a_Radry(k,i,j,b) = 1.e-10
                   a_Rawet(k,i,j,b) = 1.e-10
                END IF

                ! To cloud water
                a_rc(k,i,j) = a_rc(k,i,j) + a_maerop(k,i,j,str)

             END DO ! k
          END DO ! i
       END DO ! j
  END SUBROUTINE DiagInitAero

    !
  ! --------------------------------------------------------------------------------------------------
  ! Replacement for SUBROUTINE init_aero_sizedist (init.f90): initilize altitude-dependent aerosol
  ! size distributions and compositions.
  !
  ! Tomi Raatikainen, FMI, 29.2.2016
  !
  SUBROUTINE aerosol_init

    USE class_componentIndex, ONLY : getIndex,IsUsed
    USE mo_salsa_sizedist, ONLY : size_distribution
    USE mo_salsa_driver, ONLY : aero
    USE mo_submctl, ONLY : pi6, nbins, in1a,in2a,in2b,fn1a,fn2a,fn2b,  &
                               sigmag, dpg, n, volDistA, volDistB, nf2a, nreg,isdtyp,nspec,maxspec, &
                               rhosu, rhooc, rhobc, rhodu, rhoss, rhono, rhonh
    USE mpi_interface, ONLY : myid

    IMPLICIT NONE
    REAL :: core(nbins), nsect(1,1,nbins)             ! Size of the bin mid aerosol particle, local aerosol size dist
    REAL :: pndist(nzp,nbins)                         ! Aerosol size dist as a function of height
    REAL :: pvf2a(nzp,nspec), pvf2b(nzp,nspec)        ! Mass distributions of aerosol species for a and b-bins
    REAL :: pnf2a(nzp)                                ! Number fraction for bins 2a
    REAL :: pvfOC1a(nzp)                              ! Mass distribution between SO4 and OC in 1a
    INTEGER :: ss,ee,i,j,k
    INTEGER :: iso4=-1, ioc=-1, ibc=-1, idu=-1, &
               iss=-1, inh=-1, ino=-1

    CHARACTER(len=600) :: fmt = &
         "(/,' -------------------------------------------------',/," // &
         "' Initial aerosol profile: ',//, 4X, 'Height', 6X, 'Na', 9X, 'Nb'," // &
         "7X, 'SO4a', 8X, 'OCa', 9X, 'BCa', 9X, 'DUa', 9X, 'SSa', 9X, 'NH3a', 8X, 'HNO3a'," // &
         "7X, 'SO4b', 8X, 'OCb', 9X, 'BCb', 9X, 'DUb', 9X, 'SSb', 9X, 'NH3b', 8X, 'HNO3b'//," // &
         "(3F10.2,14ES12.3))"
    !
    ! Bin mean aerosol particle volume
    core(1:nbins) = pi6 * aero(1,1,1:nbins)%dmid**3

    ! Set concentrations to zero
    pndist = 0.
    pvf2a = 0.; pvf2b = 0.
    pnf2a = 0.; pvfOC1a = 0.

    a_maerop(:,:,:,:)=0.0
    a_naerop(:,:,:,:)=0.0

    ! Indices (-1 = not used)
    i=0
    IF (IsUsed(prtcl,'SO4')) THEN
       iso4 = GetIndex(prtcl,'SO4')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'OC')) THEN
       ioc = GetIndex(prtcl, 'OC')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'BC')) THEN
       ibc = GetIndex(prtcl,'BC')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'DU')) THEN
       idu = GetIndex(prtcl,'DU')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'SS')) THEN
       iss = GetIndex(prtcl,'SS')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'NO')) THEN
       ino = GetIndex(prtcl,'NO')
       i=i+1
    END IF
    IF (IsUsed(prtcl,'NH')) THEN
       inh = GetIndex(prtcl,'NH')
       i=i+1
    END IF

    ! All species must be known
    IF (i /= nspec) STOP 'Unknown aerosol species given in the initialization!'

    !
    ! Altitude dependent size distributions and compositions.
    ! Read and interpolate/extrapolate size distribution and composition for altitude level k: z=(zt(k)
    ! ---------------------------------------------------------------------------------------------------
    IF (isdtyp == 1) THEN

       CALL READ_AERO_INPUT(iso4,ioc,pndist,pvfOC1a,pvf2a,pvf2b,pnf2a)

    !
    ! Uniform profiles based on NAMELIST parameters
    ! ---------------------------------------------------------------------------------------------------
    ELSE IF (isdtyp == 0) THEN

       IF (ioc>0 .AND. iso4>0) THEN
          ! Both are there, so use the given "massDistrA"
          pvfOC1a(:) = volDistA(ioc)/(volDistA(ioc)+volDistA(iso4)) ! Normalize
       ELSE IF (ioc>0) THEN
          ! Pure OC
          pvfOC1a(:) = 1.0
       ELSE IF (iso4>0) THEN
          ! Pure SO4
          pvfOC1a(:) = 0.0
       ELSE
          STOP 'Either OC or SO4 must be active for aerosol region 1a!'
       ENDIF

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
       CALL size_distribution(1,1,1, n, dpg, sigmag, nsect)
       DO ss = 1,nbins
          pndist(:,ss) = nsect(1,1,ss)
       END DO

    END IF

    ! ----------------------------------------------------------

    !
    ! Initialize concentrations
    ! ----------------------------------------------------------
    DO k = 2,nzp  ! DONT PUT STUFF INSIDE THE GROUND
       DO j = 1,nyp
          DO i = 1,nxp

             ! a) Number concentrations
             ! Region 1
             a_naerop(k,i,j,in1a:fn1a) = pndist(k,in1a:fn1a)

             ! Region 2
             IF (nreg>1) THEN
                a_naerop(k,i,j,in2a:fn2a) = max(0.0,pnf2a(k))*pndist(k,in2a:fn2a)
                a_naerop(k,i,j,in2b:fn2b) = max(0.0,1.0-pnf2a(k))*pndist(k,in2a:fn2a)
             END IF

             !
             ! b) Aerosol mass concentrations
             ! bin regime 1, done here separately because of the SO4/OC convention
             ! SO4
             IF (iso4 > 0) THEN
                ss = (iso4-1)*nbins + in1a; ee = (iso4-1)*nbins + fn1a
                a_maerop(k,i,j,ss:ee) = max(0.0,1.0-pvfOC1a(k))*pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhosu
             END IF
             ! OC
             IF (ioc > 0) THEN
                ss = (ioc-1)*nbins + in1a; ee = (ioc-1)*nbins + fn1a
                a_maerop(k,i,j,ss:ee) = max(0.0,pvfOC1a(k))*pndist(k,in1a:fn1a)*core(in1a:fn1a)*rhooc
             END IF

          END DO ! i
       END DO ! j
    END DO ! k

    !
    ! c) Aerosol mass concentrations
    ! bin regime 2
    IF (nreg>1) THEN

       ! 1: SO4
       IF (iso4>0) THEN
          CALL setAeroMass(iso4,pvf2a,pvf2b,pnf2a,pndist,core,rhosu)
       END IF

       ! 2: OC
       IF (ioc>0) THEN
          CALL setAeroMass(ioc,pvf2a,pvf2b,pnf2a,pndist,core,rhooc)
       END IF

       ! 3: BC
       IF (ibc>0) THEN
          CALL setAeroMass(ibc,pvf2a,pvf2b,pnf2a,pndist,core,rhobc)
       END IF

       ! 4: DU
       IF (idu>0) THEN
          CALL setAeroMass(idu,pvf2a,pvf2b,pnf2a,pndist,core,rhodu)
       END IF

       ! 5: SS
       IF (iss>0) THEN
          CALL setAeroMass(iss,pvf2a,pvf2b,pnf2a,pndist,core,rhoss)
       END IF

       ! 6: NO3
       IF (ino>0) THEN
          CALL setAeroMass(ino,pvf2a,pvf2b,pnf2a,pndist,core,rhono)
       END IF

       ! 7: NH3
       IF (inh>0) THEN
          CALL setAeroMass(inh,pvf2a,pvf2b,pnf2a,pndist,core,rhonh)
       END IF

    END IF

    ! Put out some info about the initial state
    ! ---------------------------------------------------------------------------------------------------------------------
    IF (myid == 0)                   WRITE(*,*) ''
    IF (myid == 0 .AND. isdtyp == 0) WRITE(*,*) 'AEROSOL PROPERTIES TAKEN FROM A NAMELIST'
    IF (myid == 0 .AND. isdtyp == 1) WRITE(*,*) 'AEROSOL PROPERTIES READ FROM aerosol_in.nc'

    IF (myid == 0) WRITE(*,fmt) &
         ( zt(k), SUM(a_naerop(k,3,3,in1a:fn2a))*1.e-6, SUM(a_naerop(k,3,3,in2b:fn2b))*1.e-6,                 &

         MERGE( SUM( a_maerop(k,3,3,MAX(iso4-1,0)*nbins+in1a:MAX(iso4-1,0)*nbins+fn2a) ), -999., iso4>0 ),    &
         MERGE( SUM( a_maerop(k,3,3,MAX(ioc-1,0)*nbins+in1a:MAX(ioc-1,0)*nbins+fn2a) ), -999., ioc>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(ibc-1,0)*nbins+in1a:MAX(ibc-1,0)*nbins+fn2a) ), -999., ibc>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(idu-1,0)*nbins+in1a:MAX(idu-1,0)*nbins+fn2a) ), -999., idu>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(iss-1,0)*nbins+in1a:MAX(iss-1,0)*nbins+fn2a) ), -999., iss>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(ino-1,0)*nbins+in1a:MAX(ino-1,0)*nbins+fn2a) ), -999., ino>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(inh-1,0)*nbins+in1a:MAX(inh-1,0)*nbins+fn2a) ), -999., inh>0 ),       &

         MERGE( SUM( a_maerop(k,3,3,MAX(iso4-1,0)*nbins+in2b:MAX(iso4-1,0)*nbins+fn2b) ), -999., iso4>0 ),    &
         MERGE( SUM( a_maerop(k,3,3,MAX(ioc-1,0)*nbins+in2b:MAX(ioc-1,0)*nbins+fn2b) ), -999., ioc>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(ibc-1,0)*nbins+in2b:MAX(ibc-1,0)*nbins+fn2b) ), -999., ibc>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(idu-1,0)*nbins+in2b:MAX(idu-1,0)*nbins+fn2b) ), -999., idu>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(iss-1,0)*nbins+in2b:MAX(iss-1,0)*nbins+fn2b) ), -999., iss>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(ino-1,0)*nbins+in2b:MAX(ino-1,0)*nbins+fn2b) ), -999., ino>0 ),       &
         MERGE( SUM( a_maerop(k,3,3,MAX(inh-1,0)*nbins+in2b:MAX(inh-1,0)*nbins+fn2b) ), -999., inh>0 ),       &
         k=1,nzp )

  END SUBROUTINE aerosol_init


!!!
!!! initialize liquid and ice cloud particles ie. move particle fractions from aerosol bins to liquid cloud and ice particle bins
!!!

 ! ---------- Juha: This should be replaced ASAP with a physical treatment. Do NOT use for liquid clouds.
  SUBROUTINE liq_ice_init

    !USE mo_salsa_driver, ONLY : aero
    USE mo_submctl, ONLY : nbins, in2a,in2b,fn2a,fn2b,  &
                               nspec, ncld, nice, initliqice, &
                               liqFracA,iceFracA,liqFracB,iceFracB

    IMPLICIT NONE

    INTEGER :: i,j,k,bb,m
    REAL :: zumA, zumB, zumCumIce, zumCumLiq, &
            excessIce, excessLiq,excessFracIce,excessFracLiq

    ! initialize liquid and ice only if it is determinded so in the namelist.salsa
    IF(initliqice) THEN
        IF(level==4) THEN
            iceFracA = 0.0; iceFracB = 0.0;
        END IF
        !#cloudinit
        DO k = 2,nzp  ! DONT PUT STUFF INSIDE THE GROUND
           DO j = 1,nyp
              DO i = 1,nxp
                 IF (a_rh(k,i,j)<1.0  .or. a_scr1(k,i,j) > 273.15) CYCLE
                 zumA=sum(a_naerop(k,i,j,in2a:fn2a))
                 zumB=sum(a_naerop(k,i,j,in2b:fn2b))

                 zumCumIce = 0.0
                 zumCumLiq = 0.0
                 excessIce = 0.0
                 excessFracIce = 1.0
                 excessFracLiq = 1.0

                 DO bb=fn2a,in2a,-1

                    IF(a_scr1(k,i,j) < 273.15 .and. zumCumIce<iceFracA*zumA .and. a_naerop(k,i,j,bb)>10e-10) THEN !initialize ice if it is cold enough

                       excessIce =min(abs(zumCumIce-iceFracA*zumA),a_naerop(k,i,j,bb))
                       a_nicep(k,i,j,bb-3) = a_nicep(k,i,j,bb-3) + excessIce
                       zumCumIce = zumCumIce + excessIce

                       excessFracIce = excessIce/a_naerop(k,i,j,bb)
                       excessFracIce = MAX(0.0,MIN(1.0,excessFracIce))
                       a_naerop(k,i,j,bb)=(1.0-excessFracIce)*a_naerop(k,i,j,bb)

                       DO m=1,nspec
                          a_micep(k,i,j,(m-1)*nice+bb-3) = &
                             excessFracIce*a_maerop(k,i,j,(m-1)*nbins+bb)

                          a_maerop(k,i,j,(m-1)*nbins+bb) = &
                            (1.0-excessFracIce)*a_maerop(k,i,j,(m-1)*nbins+bb)
                       END DO

                    END IF

                    IF (a_rh(k,i,j)>1.0 .and. zumCumLiq<liqFracA*zumA .and. a_naerop(k,i,j,bb)>10e-10) THEN

                       excessLiq =min(abs(zumCumLiq-liqFracA*zumA),a_naerop(k,i,j,bb))
                       a_ncloudp(k,i,j,bb-3) = a_ncloudp(k,i,j,bb-3) + excessLiq
                       zumCumLiq = zumCumLiq + excessLiq

                       excessFracLiq = excessLiq/a_naerop(k,i,j,bb)
                       excessFracLiq = MAX(0.0,MIN(1.0,excessFracLiq))

                       a_naerop(k,i,j,bb)=(1.0-excessFracLiq)*a_naerop(k,i,j,bb)

                       DO m=1,nspec
                          a_mcloudp(k,i,j,(m-1)*ncld+bb-3) = &
                              excessFracLiq*a_maerop(k,i,j,(m-1)*nbins+bb)

                          a_maerop(k,i,j,(m-1)*nbins+bb) = &
                             (1.0-excessFracLiq)*a_maerop(k,i,j,(m-1)*nbins+bb)
                       END DO

                    END IF

                 END DO ! fn2a

                zumCumIce = 0.0
                zumCumLiq = 0.0
                excessFracIce = 1.0
                excessFracLiq = 1.0
                DO bb=fn2b,in2b,-1
                   IF(a_scr1(k,i,j) < 273.15 .and. zumCumIce<iceFracB*zumB  .and. a_naerop(k,i,j,bb)>10e-10) THEN !initialize ice if it is cold enough

                      excessIce =min(abs(zumCumIce-iceFracB*zumB),a_naerop(k,i,j,bb))
                      a_nicep(k,i,j,bb-3) = a_nicep(k,i,j,bb-3) + excessIce
                      zumCumIce = zumCumIce + excessIce

                      excessFracIce = excessIce/a_naerop(k,i,j,bb)
                      excessFracIce = MAX(0.0,MIN(1.0,excessFracIce))

                      a_naerop(k,i,j,bb)=(1.0-excessFracIce)*a_naerop(k,i,j,bb)

                      DO m=1,nspec
                         a_micep(k,i,j,(m-1)*nice+bb-3) = &
                             excessFracIce*a_maerop(k,i,j,(m-1)*nbins+bb)
 
                         a_maerop(k,i,j,(m-1)*nbins+bb) = &
                            (1.0-excessFracIce)*a_maerop(k,i,j,(m-1)*nbins+bb)
                      END DO
                   END IF

                   IF (a_rh(k,i,j)>1.0 .and. zumCumLiq<liqFracB*zumB .and. a_naerop(k,i,j,bb)>10e-10) THEN

                      excessLiq =min(abs(zumCumLiq-liqFracB*zumB),a_naerop(k,i,j,bb))
                      a_ncloudp(k,i,j,bb-3) = a_ncloudp(k,i,j,bb-3) + excessLiq
                      zumCumLiq = zumCumLiq + excessLiq

                      excessFracLiq = excessLiq/a_naerop(k,i,j,bb)
                      excessFracLiq = MAX(0.0,MIN(1.0,excessFracLiq))

                      a_naerop(k,i,j,bb)=(1.0-excessFracLiq)*a_naerop(k,i,j,bb)

                      DO m=1,nspec
                         a_mcloudp(k,i,j,(m-1)*ncld+bb-3) = &
                             excessFracLiq*a_maerop(k,i,j,(m-1)*nbins+bb)

                         a_maerop(k,i,j,(m-1)*nbins+bb) = &
                             (1.0-excessFracLiq)*a_maerop(k,i,j,(m-1)*nbins+bb)
                      END DO

                   END IF
                END DO ! fn2b

             END DO ! i
          END DO ! j
       END DO ! 
    END IF

END SUBROUTINE liq_ice_init

  !
  ! ----------------------------------------------------------
  ! Sets the mass concentrations to aerosol arrays in 2a and 2b
  !
  !
  SUBROUTINE setAeroMass(ispec,ppvf2a,ppvf2b,ppnf2a,ppndist,pcore,prho)
    USE mo_submctl, ONLY : nbins, in2a,fn2a,in2b,fn2b, nspec
    IMPLICIT NONE

    INTEGER, INTENT(in) :: ispec                             ! Aerosol species index
    REAL, INTENT(in) :: ppvf2a(nzp,nspec), ppvf2b(nzp,nspec) ! Mass distributions for a and b bins
    REAL, INTENT(in) :: ppnf2a(nzp)                          ! Number fraction for 2a
    REAL, INTENT(in) :: ppndist(nzp,nbins)                   ! Aerosol size distribution
    REAL, INTENT(in) :: pcore(nbins)                         ! Aerosol bin mid core volume
    REAL, INTENT(in) :: prho                                 ! Aerosol density

    INTEGER :: ss,ee
    INTEGER :: i,j,k

    DO k = 2,nzp ! DONT PUT STUFF INSIDE THE GROUND
       DO j = 1,nyp
          DO i = 1,nxp
             ! 2a
             ss = (ispec-1)*nbins + in2a; ee = (ispec-1)*nbins + fn2a
             a_maerop(k,i,j,ss:ee) =      &
                  max( 0.0,ppvf2a(k,ispec) )*ppnf2a(k) * &
                  ppndist(k,in2a:fn2a)*pcore(in2a:fn2a)*prho
             ! 2b
             ss = (ispec-1)*nbins + in2b; ee = (ispec-1)*nbins + fn2b
             a_maerop(k,i,j,ss:ee) =      &
                  max( 0.0,ppvf2b(k,ispec) )*(1.0-ppnf2a(k)) * &
                  ppndist(k,in2a:fn2a)*pcore(in2a:fn2a)*prho
          END DO
       END DO
    END DO

  END SUBROUTINE setAeroMass
  !
  ! -------------------------------------------------------------------------
  ! Reads vertical profiles of aerosol size distribution parameters, aerosol species volume fractions and
  ! number concentration fractions between a and b bins
  !
  SUBROUTINE READ_AERO_INPUT(piso4,pioc,ppndist,ppvfOC1a,ppvf2a,ppvf2b,ppnf2a)
    USE ncio, ONLY : open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_aero_nc
    USE mo_submctl, ONLY : nbins,  &
                               nspec, maxspec, nmod
    USE mo_salsa_sizedist, ONLY : size_distribution
    USE mpi_interface, ONLY : appl_abort, myid
    IMPLICIT NONE

    INTEGER, INTENT(in) :: piso4,pioc
    REAL, INTENT(out) :: ppndist(nzp,nbins)                   ! Aerosol size dist as a function of height
    REAL, INTENT(out) :: ppvf2a(nzp,nspec), ppvf2b(nzp,nspec) ! Volume distributions of aerosol species for a and b-bins
    REAL, INTENT(out) :: ppnf2a(nzp)                          ! Number fraction for bins 2a
    REAL, INTENT(out) :: ppvfOC1a(nzp)                        ! Volume distribution between SO4 and OC in 1a

    REAL :: nsect(1,1,nbins)

    INTEGER :: ncid, k, i
    INTEGER :: nc_levs=500, nc_nspec, nc_nmod

    ! Stuff that will be read from the file
    REAL, ALLOCATABLE :: zlevs(:),         &  ! Levels in meters
                         zvolDistA(:,:),  &  ! Volume distribution of aerosol species in a and b bins
                         zvoldistB(:,:),  &  ! (Don't mess these with the ones in namelist.salsa -
                                             !  they are not used here!)
                         znf2a(:),        &  ! Number fraction for bins 2a
                         zn(:,:),         &  ! Aerosol mode number concentrations
                         zsigmag(:,:),    &  ! Geometric standard deviations
                         zdpg(:,:),       &  ! Mode mean diameters
                         znsect(:,:),     &  ! Helper for binned number concentrations
                         helper(:,:)         ! nspec helper
    LOGICAL :: READ_NC

    ! Read the NetCDF input when it is available
    INQUIRE(FILE='aerosol_in.nc',EXIST=READ_NC)

    ! Open the input file
    IF (READ_NC) CALL open_aero_nc(ncid, nc_levs, nc_nspec, nc_nmod)

    ! Check that the input dimensions are compatible with SALSA initialization
    ! ....

    ! Allocate input variables
    ALLOCATE( zlevs(nc_levs),               &
              zvolDistA(nc_levs,maxspec),  &
              zvolDistB(nc_levs,maxspec),  &
              znf2a(nc_levs),              &
              zn(nc_levs,nmod),            &
              zsigmag(nc_levs,nmod),       &
              zdpg(nc_levs,nmod),          &
              ! Couple of helper arrays
              znsect(nc_levs,nbins),        &
              helper(nc_levs,nspec)        )

    zlevs = 0.; zvolDistA = 0.; zvolDistB = 0.; znf2a = 0.; zn = 0.; zsigmag = 0.
    zdpg = 0.; znsect = 0.; helper = 0.

    IF (READ_NC) THEN
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
    znsect = 0.
    DO k = 1,nc_levs
       CALL size_distribution(1,1,1,zn(k,:),zdpg(k,:),zsigmag(k,:),nsect)
       znsect(k,:) = nsect(1,1,:)
    END DO

    ! Interpolate the input variables to model levels
    ! ------------------------------------------------
    CALL htint2d(nc_levs,zvolDistA(1:nc_levs,1:nspec),zlevs(1:nc_levs),nzp,ppvf2a,zt,nspec)
    CALL htint2d(nc_levs,zvolDistB(1:nc_levs,1:nspec),zlevs(1:nc_levs),nzp,ppvf2b,zt,nspec)
    CALL htint2d(nc_levs,znsect(1:nc_levs,:),zlevs(1:nc_levs),nzp,ppndist,zt,nbins)
    CALL htint(nc_levs,znf2a(1:nc_levs),zlevs(1:nc_levs),nzp,ppnf2a,zt)

    ! Since 1a bins by SALSA convention can only contain SO4 or OC,
    ! get renormalized mass fractions.
    ! --------------------------------------------------------------
    IF (pioc>0 .AND. piso4>0) THEN
       ! Both are there, so use the given "massDistrA"
       ppvfOC1a(:) = ppvf2a(:,pioc)/(ppvf2a(:,pioc)+ppvf2a(:,piso4)) ! Normalize
    ELSE IF (pioc>0) THEN
       ! Pure OC
       ppvfOC1a(:) = 1.0
    ELSE IF (piso4>0) THEN
       ! Pure SO4
       ppvfOC1a(:) = 0.0
    ELSE
       STOP 'Either OC or SO4 must be active for aerosol region 1a!'
    ENDIF

    DEALLOCATE( zlevs, zvolDistA, zvolDistB, znf2a, zn, zsigmag, zdpg, znsect, helper )

  END SUBROUTINE READ_AERO_INPUT

  !
  !------------------------------------------------------------------
  ! INIT_GAS_TRACERS: Set initial values for gas compound tracers
  !
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE init_gas_tracers
    IMPLICIT NONE

    INTEGER :: j,i,k

    ! Nääkin voisi lukea tiedostosta
    ! Taken as molecules/kg
    DO j = 1,nyp
       DO i = 1,nxp
          DO k = 1,nzp
             a_gaerop(k,i,j,1) = 5.E14/dn0(k) !SO4
             a_gaerop(k,i,j,2) = 0./dn0(k)    !NO3
             a_gaerop(k,i,j,3) = 0./dn0(k)    !NH4
             a_gaerop(k,i,j,4) = 5.E14/dn0(k) !OCNV
             a_gaerop(k,i,j,5) = 1.E14/dn0(k) !OCSV
          END DO
       END DO
    END DO


  END SUBROUTINE init_gas_tracers



end module init
