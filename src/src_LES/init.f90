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
MODULE init

   USE grid
   USE mo_progn_state
   USE mo_vector_state
   USE mo_diag_state
   USE mo_aux_state
   USE mo_history, ONLY : read_hist, write_hist
   
   INTEGER, PARAMETER    :: nns = 1500
   INTEGER               :: ns
   INTEGER               :: iseed = 0
   INTEGER               :: ipsflg = 1
   INTEGER               :: itsflg = 1
        ! itsflg = 0 :potential temperature in kelvin
        !          1 :liquid water potential temperature in kelvin
        !          2 :temperature
   REAL, DIMENSION(nns)  :: us,vs,ts,thds,ps,hs,rts,rss,tks,xs
   REAL                  :: zrand = 200.
   REAL                  :: zrndamp = 0.2 ! the amplitude of pseudorandom fluctuations
   CHARACTER  (len=80)   :: hfilin = 'test.'

   INTEGER               :: init_type = 1 ! Switch for how to initialize mixing: 1) random perturbations, 2) warm bubble

CONTAINS
   !
   ! ----------------------------------------------------------------------
   ! INITLZ:  this is the main driver for the model's initializ-
   ! ation routines.  it initializes the model according to runtype
   !
   SUBROUTINE initialize

      USE step, ONLY : time, outflg
      !USE stat, ONLY : init_stat, mcflg, acc_massbudged, salsa_b_bins
      USE sgsm, ONLY : tkeinit
      USE mpi_interface, ONLY : appl_abort, myid
      USE thrm, ONLY : thermo
      USE mo_salsa_driver, ONLY : run_SALSA
      USE mo_submctl, ONLY : in2b, fn2b, nlim, prlim, spec
      USE util, ONLY : maskactiv
      USE nudg, ONLY : init_nudg
      USE emission_init, ONLY : init_emission
      USE constrain_SALSA, ONLY : SALSA_diagnostics
      USE mo_structured_datatypes
      USE mo_output, ONLY : init_main, write_main, init_ps, write_ps
      USE mo_field_types, ONLY : Diag, Prog
      
      IMPLICIT NONE

      ! Local variables for SALSA basic state
      REAL    :: zwp(nzp,nxp,nyp)
      INTEGER :: n4
   
      ! Set vertical velocity as 0.5 m/s to intialize cloud microphysical properties with
      ! SALSA
      zwp(:,:,:) = 0.5

      IF (runtype == 'INITIAL') THEN
         time = 0.
         CALL arrsnd
         CALL basic_state
         CALL fldinit ! Juha: aerosol size distributions are initialized here.
                      !       Also thermodynamics!
         
         ! If SALSA is used, call SALSA with full configuration once before beginning
         ! spin-up period to set up aerosol and cloud fields.
         IF (level >= 4) THEN

            n4 = spec%getNSpec(type="wet")

            IF ( nxp == 5 .AND. nyp == 5 ) THEN
               CALL run_SALSA(Diag,Prog,nzp,nxp,nyp,n4,   &
                              zwp,a_nactd,a_vactd,dtlt,   &
                              time,level,.TRUE.           )
            ELSE
               CALL run_SALSA(Diag,Prog,nzp,nxp,nyp,n4,   &
                              a_wp%d,a_nactd,a_vactd,dtlt,  &
                              time, level,.TRUE.          )
            END IF
            CALL SALSAInit

         END IF !level >= 4
         
      ELSE IF (runtype == 'HISTORY') THEN
         IF (isgstyp == 2) CALL tkeinit(nxyzp,a_qp%d)
         CALL hstart
      ELSE
         IF (myid == 0) PRINT *,'  ABORTING:  Invalid Runtype'
         CALL appl_abort(0)
      END IF ! runtype

     ! When SALSA b-bin outputs are needed?
     !   -level >= 4
     !   -outputs are forced (lsalsabbins=.true.)
     !   -b-bins initialized with non-zero concentration
     !   -nucleation set to produce particles to b bins (currently only a bins)
     IF (level >= 4 .AND. (.NOT. lsalsabbins)) &
        lsalsabbins = ANY( a_naerop%d(:,:,:,in2b:fn2b) > nlim ) 

     CALL sponge_init
     !CALL init_stat(time+dtl,filprf,expnme,nzp)
     !
     ! Initialize nudging profiles
     ! ----------------------------
     !Ali, for history run parmeters of init_nudg are read from history file
     IF ((runtype == 'INITIAL') .AND. (lnudging)) CALL init_nudg()
     !
     ! Initialize aerosol emissions
     ! -----------------------------
     IF (lemission .AND. level >= 4) CALL init_emission()

     !
     !IF (mcflg) THEN
     !   ! Juha:
     !   ! Calculate some numbers for mass concervation experiments
     !   mc_Vdom = deltax*deltay*deltaz*(nxp-4)*(nyp-4)*(nzp-1)
     !   mc_Adom = deltax*deltay*(nxp-4)*(nyp-4)
     !   mc_ApVdom = mc_Adom/mc_Vdom
     !   ! Get the initial mass of atmospheric water
     !   CALL acc_massbudged(nzp,nxp,nyp,0,dtlt,dzt,a_dn,     &
     !                       rv=a_rp,rc=a_rc,prc=a_srp)
     !END IF ! mcflg

     ! Diagnostic calculations that should take place (with SALSA) both for INITIAL and HISTORY
     IF ( (level >= 4) ) THEN
        !CALL thermo(level)
        CALL SALSA_diagnostics(onlyDiag=.TRUE.)
        CALL thermo(level)
     END IF

     !
      ! write analysis and history files from restart if appropriate
      !
      IF (outflg) THEN
         IF (runtype == 'INITIAL') THEN
            CALL write_hist(1, time)
            CALL init_main(time)
            CALL init_ps(time)
            CALL thermo(level)
            CALL write_main(time)
         ELSE
            CALL init_main(time+dtl)
            CALL write_hist(0, time)
         END IF
      END IF !outflg

      RETURN
   END SUBROUTINE initialize
   !
   !----------------------------------------------------------------------
   ! FLDINIT: Initializeds 3D fields, mostly from 1D basic state
   !
   !          Modified for level 4.
   !          Juha Tonttila, FMI, 20140828
   !
   SUBROUTINE fldinit

      USE defs, ONLY : alvl, cpr, cp, p00
      USE sgsm, ONLY : tkeinit
      USE thrm, ONLY : thermo, rslf
      USE init_warm_bubble, ONLY : warm_bubble

      IMPLICIT NONE

      INTEGER :: i,j,k
      REAL    :: exner, pres, tk, rc, xran(nzp)
      
      CALL htint(ns,ts,hs,nzp,th0%d,zt%d)
      
      DO j = 1, nyp
         DO i = 1, nxp
            a_ustar%d(i,j) = 0.
            DO k = 1, nzp
               a_up%d(k,i,j)    = u0%d(k)
               a_vp%d(k,i,j)    = v0%d(k)
               a_tp%d(k,i,j)    = (th0%d(k)-th00)
               IF (associated (a_rp%d)) a_rp%d(k,i,j)   = rt0%d(k)
               a_theta%d(k,i,j) = th0%d(k)
               a_pexnr%d(k,i,j) = 0.
            END DO
         END DO
      END DO

      ! Juha: Added SELECT-CASE for level 4
      SELECT CASE(level)
      CASE(1,2,3)
         IF ( ASSOCIATED (a_rv%d)) a_rv%d = a_rp%d
         
         IF ( ASSOCIATED (a_rc%d)) THEN
            DO j = 1, nyp
               DO i = 1, nxp
                  DO k = 1, nzp
                     exner = (pi0%d(k)+pi1%d(k))/cp
                     pres  = p00 * (exner)**cpr
                     IF (itsflg == 0) THEN
                        tk = th0%d(k)*exner
                        rc = max(0.,a_rp%d(k,i,j)-rslf(pres,tk))
                        a_tp%d(k,i,j) = a_theta%d(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                        a_rv%d(k,i,j) = a_rp%d(k,i,j)-rc
                     END IF
                     IF (itsflg == 2) THEN
                        tk = th0%d(k)
                        a_theta%d(k,i,j) = tk/exner
                        rc = max(0.,a_rp%d(k,i,j)-rslf(pres,tk))
                        a_tp%d(k,i,j) = a_theta%d(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                        a_rv%d(k,i,j) = a_rp%d(k,i,j)-rc
                     END IF
                  END DO
               END DO
            END DO
         END IF
         
      CASE(4,5)
         ! Condensation will be calculated by the initial call of SALSA, so use the
         ! saturation adjustment method to estimate the amount of liquid water,
         ! which is needed for theta_l
         DO j = 1, nyp
            DO i = 1, nxp
               DO k = 1, nzp
                  exner = (pi0%d(k)+pi1%d(k))/cp
                  pres  = p00 * (exner)**cpr
                  IF (itsflg == 0) THEN
                     tk = th0%d(k)*exner
                     rc = max(0.,a_rp%d(k,i,j)-rslf(pres,tk))
                     a_tp%d(k,i,j) = a_theta%d(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                  END IF
                  IF (itsflg == 2) THEN
                     tk = th0%d(k)
                     a_theta%d(k,i,j) = tk/exner
                     rc = max(0.,a_rp%d(k,i,j)-rslf(pres,tk))
                     a_tp%d(k,i,j) = a_theta%d(k,i,j)*exp(-(alvl/cp)*rc/tk) - th00
                  END IF
               END DO !k
            END DO !i
         END DO !j
          
      END SELECT
      
      IF (init_type == 1) THEN
         ! Initialize with random perturbations
         k = 1
         DO WHILE( zt%d(k+1) <= zrand .AND. k+1 < nzp)
            k = k+1
            xran(k) = zrndamp*(zrand - zt%d(k))/zrand
         END DO
         CALL random_pert(nzp,nxp,nyp,zt,a_tp,xran,k)
         
         IF (associated(a_rp%d)) THEN
            k = 1
            DO WHILE( zt%d(k+1) <= zrand .AND. k+1 < nzp)
               k = k+1
               xran(k) = 5.0e-5*(zrand - zt%d(k))/zrand
            END DO
            CALL random_pert(nzp,nxp,nyp,zt,a_rp,xran,k)
         END IF
      ELSE IF (init_type == 2) THEN
         ! Initialize with warm bubble (for convection)
         CALL warm_bubble()
      END IF

      a_wp%d = 0.
      IF(isgstyp == 2) CALL tkeinit(nxyzp,a_qp%d)
      !
      ! initialize thermodynamic fields
      !
      CALL thermo (level)
      
      !
      ! Initialize aerosol size distributions
      !
      IF (level >= 4) THEN
         CALL aerosol_init()
         CALL init_gas_tracers
      END IF

      a_uc%d = a_up%d
      a_vc%d = a_vp%d
      a_wc%d = a_wp%d

      RETURN
   END SUBROUTINE fldinit
   !----------------------------------------------------------------------
   ! SPONGE_INIT: Initializes variables for sponge layer
   !
   SUBROUTINE sponge_init

      USE mpi_interface, ONLY : myid

      IMPLICIT NONE

      INTEGER :: k,kk

      IF (nfpt > 0) THEN
         ALLOCATE (spng_tfct(max(1,nfpt)), spng_wfct(max(1,nfpt)))

         DO k = nzp-nfpt, nzp-1
            kk = k + 1 - (nzp-nfpt)
            spng_tfct(kk) = max(0.,(zm%d(nzp)-zt%d(k))/((zm%d(nzp)-zm%d(nzp-nfpt))*distim))
            spng_wfct(kk) = max(0.,(zm%d(nzp)-zm%d(k))/((zm%d(nzp)-zm%d(nzp-nfpt))*distim))
            spng_tfct(kk) = max(0.,(1./distim - spng_tfct(kk)))
            spng_wfct(kk) = max(0.,(1./distim - spng_wfct(kk)))
         END DO

         IF(myid == 0) THEN
            PRINT "(//' ',49('-')/)"
            PRINT '(2X,A17)', 'Sponge Layer Init '
            PRINT '(3X,A12,F8.1,A1)', 'Starting at ', zt%d(nzp-nfpt), 'm'
            PRINT '(3X,A18,F6.1,A1)', 'Minimum timescale ', 1/spng_wfct(nfpt),'s'
         END IF
      END IF

      RETURN
   END SUBROUTINE sponge_init

   !
   !
   ! ----------------------------------------------------------------------
   ! ARRSND: Arranges the sounding input into proper arrays
   !
   SUBROUTINE arrsnd

      USE defs, ONLY          : p00,p00i,cp,cpr,rcp,r,g,ep2,alvl,Rm,ep
      USE thrm, ONLY          : rslf
      USE mpi_interface, ONLY : appl_abort, myid

      IMPLICIT NONE

      INTEGER :: k, iterate
      REAL    :: tavg, zold2, zold1, x1, xx, yy, zz, til
      CHARACTER (len=245) :: fm0 = &
         "(/,' -------------------------------------------------',/,"       //&
         "'  Sounding Input: ',//,7x,'ps',9x,'hs',7x,'ts',6x ,'thds',6x," // &
         "'us',7x,'vs',7x,'rts',5x,'rel hum',/,6x,'(Pa)',7X,'(m)',6X,'(K)'"// &
         ",6X,'(K)',6X,'(m/s)',4X,'(m/s)',3X,'(kg/kg)',5X,'(%)',/,1x/)"
      CHARACTER (len=36) :: fm1 = "(f11.1,f10.1,2f9.2,2f9.2,f10.5,f9.1)"
      !
      ! arrange the input sounding
      !
      IF (ps(1) == 0.) THEN
         OPEN(1,file='sound_in',status='old',form='formatted')

         DO ns = 1, nns
            READ(1,*,end=100) ps(ns),ts(ns),rts(ns),us(ns),vs(ns)
         END DO

         CLOSE(1)
      END IF
100 CONTINUE

    zold1 = 0.
    zold2 = 0.

    ns = 1
    DO WHILE (ps(ns) /= 0. .AND. ns <= nns)
       !
       ! filling relative humidity array only accepts sounding in mixing
       ! ratio (g/kg) converts to (kg/kg)
       !
       rts(ns) = rts(ns)*1.e-3
       !
       ! filling pressure array:
       ! ipsflg = 0 :pressure in millibars
       ! 1 :pressure array is height in meters (ps(1) is surface pressure)
       !
       SELECT CASE (ipsflg)

          CASE (0)
             ps(ns) = ps(ns)*100.

          CASE DEFAULT
             xs(ns) = (1.+ep2*rts(ns))
             IF (ns == 1)then
                ps(ns) = ps(ns)*100.
                zold2 = 0.
                hs(1) = 0.
             ELSE
                hs(ns) = ps(ns)
                zold1 = zold2
                zold2 = ps(ns)
                IF ( itsflg == 0 .OR. itsflg == 1) THEN
                   ! ts=potential or liquid water potential temperature (condensation not included here)
                   tavg=0.5*(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1))*((p00/ps(ns-1))**rcp)
                ELSE
                   ! ts=T [K]
                   tavg = 0.5*(ts(ns)*xs(ns)+ts(ns-1)*xs(ns-1))
                END IF
                ps(ns) = (ps(ns-1)**rcp-g*(zold2-zold1)*(p00**rcp)/(cp*tavg))**cpr
             END IF

       END SELECT
       !
       ! filling temperature array:
       ! itsflg = 0 :potential temperature in kelvin
       !          1 :liquid water potential temperature in kelvin
       !          2 :temperature
       !
       SELECT CASE (itsflg)

          CASE (0)
             tks(ns) = ts(ns)*(ps(ns)*p00i)**rcp

          CASE (1)
             til = ts(ns)*(ps(ns)*p00i)**rcp
             xx = til
             yy = rslf(ps(ns),xx)
             zz = max(rts(ns)-yy,0.)
             IF (zz > 0.) THEN
                DO iterate = 1, 3
                   x1 = alvl/(cp*xx)
                   xx = xx - (xx - til*(1.+x1*zz))/(1. + x1*til                &
                        *(zz/xx+(1.+yy*ep)*yy*alvl/(Rm*xx*xx)))
                   yy = rslf(ps(ns),xx)
                   zz = max(rts(ns)-yy,0.)
                END DO
             END IF
             tks(ns) = xx

          CASE (2)
             tks(ns) = ts(ns) ! a long way of saying do nothing

          CASE DEFAULT
             IF (myid == 0) PRINT *, '  ABORTING: itsflg not supported'
             CALL appl_abort(0)

       END SELECT
       ns = ns+1
    END DO
    ns = ns-1
    !
    ! compute height levels of input sounding.
    !
    IF (ipsflg == 0) THEN
       DO k = 2, ns
          hs(k) = hs(k-1)-r*.5 *(tks(k)*(1.+ep2*rts(k))                      &
                  +tks(k-1)*(1.+ep2*rts(k-1)))*(log(ps(k))-log(ps(k-1)))/g
       END DO
    END IF

    IF (hs(ns) < zt%d(nzp)) THEN
       IF (myid == 0) PRINT *, '  ABORTING: Model top above sounding top'
       IF (myid == 0) PRINT '(2F12.2)', hs(ns), zt%d(nzp)
       CALL appl_abort(0)
    END IF

    DO k = 1, ns
       thds(k) = tks(k)*(p00/ps(k))**rcp
    END DO

    DO k = 1, ns
       xs(k) = 100.*rts(k)/rslf(ps(k),tks(k))
    END DO

    IF(myid == 0) THEN
       WRITE(6,fm0)
       WRITE(6,fm1)(ps(k),hs(k),tks(k),thds(k),us(k),vs(k),rts(k),xs(k),k=1,ns)
    END IF

    RETURN
 END SUBROUTINE arrsnd
 !
 !----------------------------------------------------------------------
 ! BASIC_STATE: This routine computes the basic state values
 ! of pressure, density, moisture and temperature.  The basic state
 ! temperature is assumed to be a the volume weighted average value of
 ! the sounding
 !
 SUBROUTINE basic_state

    USE defs, ONLY : cp, rcp, cpr, r, g, p00, p00i, ep2
    USE mpi_interface, ONLY : myid
    USE thrm, ONLY : rslf

    IMPLICIT NONE

    INTEGER :: k
    REAL    :: v1da(nzp), v1db(nzp), v1dc(nzp), exner, zrh(nzp),ztk(nzp)

    CHARACTER (len=328) :: fmt =  &
       "(/,' -------------------------------------------------',/,"     //&
       "'  Basic State: ',//,4X,'Z',6X,'U0',6X,'V0',6X,'DN0',6X,' P0'"   //&
       ",6X,'PRESS',4X,'TH0',6X,'THV',5X,'RT0','RH',/,3X,'(m)',5X,'(m/s)'"     //&
       ",3X,'(m/s)',2X,'(kg/m3)',2X,'(J/kgK)',4X,'(Pa)',5X,'(K)',5X"      //&
       ",'(K)',4X,'(g/kg)',4X,'1',4X,'K',//,(1X,F7.1,2F8.2,F8.3,2F10.2,2F8.2,F7.2,F7.2,F7.2))"

    CALL htint(ns,thds,hs,nzp,th0%d,zt%d)
    CALL htint(ns,us,hs,nzp,u0%d,zt%d)
    CALL htint(ns,vs,hs,nzp,v0%d,zt%d)

    IF (level >= 1) THEN
       CALL htint(ns,rts,hs,nzp,rt0%d,zt%d)
       rt0%d(1) = rt0%d(2)
    ELSE
       DO k = 1, nzp
          rt0%d(k) = 0.
       END DO
    END IF
    !
    ! calculate theta_v for an unsaturated layer, neglecting condensate here is
    ! okay as this is only used for the first estimate of pi1, which will be
    ! updated in a consistent manner on the first dynamic timestep
    !
    DO k = 1, nzp
       v1dc(k) = th0%d(k) * (1.+ep2*rt0%d(k)) ! theta_v assuming unsaturated
    END DO
    !
    ! calculate pressure for actual initial state
    !
    pi1%d(1) = cp*(ps(1)*p00i)**rcp+g*(hs(1)-zt%d(1))/v1dc(1)
    DO k = 2, nzp
       pi1%d(k) = pi1%d(k-1)-g/(dzm%d(k-1)*0.5*(v1dc(k)+v1dc(k-1)))
    END DO
    !
    ! calculate hydrostatic exner function associated with th00 constant along
    ! with associated basic state density
    !
    pi0%d(1) = cp*(ps(1)*p00i)**rcp + g*(hs(1)-zt%d(1))/th00
    dn0%d(1) = ((cp**(1.-cpr))*p00)/(r*th00*pi0%d(1)**(1.-cpr))
    DO k = 2, nzp
       pi0%d(k) = pi0%d(1) + g*(zt%d(1) - zt%d(k))/th00
       dn0%d(k) = ((cp**(1.-cpr))*p00)/(r*th00*pi0%d(k)**(1.-cpr))
       u0%d(k) = u0%d(k)-umean
       v0%d(k) = v0%d(k)-vmean
    END DO
    !
    ! define pi1 as the difference between pi associated with th0 and pi
    ! associated with th00, thus satisfying pi1+pi0 = pi = cp*(p/p00)**(R/cp)
    !
    DO k = 1, nzp
       pi1%d(k) = pi1%d(k)-pi0%d(k)
    END DO
    !
    DO k = 1, nzp
       exner = (pi0%d(k) + pi1%d(k))/cp
       v1db(k) = p00*(exner)**cpr      ! pressure
       v1da(k) = p00*(pi0%d(k)/cp)**cpr  ! pressure associated with pi0
    END DO

    u0%d(1) = u0%d(2)
    v0%d(1) = v0%d(2)
    psrf  = ps(1)

    ! Juha: For debugging
    zrh = 0.
    ztk = 0.
    DO k = 1,nzp
       exner = (pi0%d(k) + pi1%d(k))/cp
       ztk(k) = th0%d(k)*(v1db(k)*p00i)**rcp
       zrh(k) = 100.*rt0%d(k)/rslf(v1db(k),ztk(k))
    END DO
       
    IF(myid == 0) WRITE(*,fmt) (zt%d(k),u0%d(k),v0%d(k),dn0%d(k),v1da(k),v1db(k), &
                                th0%d(k),v1dc(k),rt0%d(k)*1000.,zrh(k),ztk(k),k=1,nzp)  ! rt0(k)*1000.

    RETURN
 END SUBROUTINE basic_state
 !
 !---------------------------------------------------------------------
 ! HTINT: Height interpolation of field on one grid, to field on another
 !
 SUBROUTINE htint(na,xa,za,nb,xb,zb)

    IMPLICIT NONE
    INTEGER, INTENT (in) :: na, nb
    REAL, INTENT (in)    :: xa(na),za(na),zb(nb)
    REAL, INTENT (out)   :: xb(nb)

    INTEGER :: l, k
    REAL    :: wt

    l = 1
    DO k = 1, nb
       IF (zb(k) <= za(na)) THEN
          DO WHILE ( zb(k) > za(l+1) .AND. l < na)
             l = l+1
          END DO
          wt = (zb(k)-za(l))/(za(l+1)-za(l))
          xb(k) = xa(l)+(xa(l+1)-xa(l))*wt
       ELSE
          wt = (zb(k)-za(na))/(za(na-1)-za(na))
          xb(k) = xa(na)+(xa(na-1)-xa(na))*wt
       END IF
    END DO

    RETURN
 END SUBROUTINE htint
 !
 ! -----------------------------------------------------------------------
 ! HTINT2d: Same as HTINT but for 2d variables
 !
 SUBROUTINE htint2d(na,xa,za,nb,xb,zb,nx)
    IMPLICIT NONE

    INTEGER, INTENT (in) :: na, nb, nx
    REAL, INTENT (in)    :: xa(na,nx),za(na),zb(nb)
    REAL, INTENT (out)   :: xb(nb,nx)

    INTEGER :: l, k, i
    REAL    :: wt

    DO i = 1, nx
       l = 1
       DO k = 1, nb
          IF (zb(k) <= za(na)) THEN
             DO WHILE ( zb(k) > za(l+1) .AND. l < na)
                l = l+1
             END DO
             wt = (zb(k)-za(l))/(za(l+1)-za(l))
             xb(k,i) = xa(l,i)+(xa(l+1,i)-xa(l,i))*wt
          ELSE
             wt = (zb(k)-za(na))/(za(na-1)-za(na))
             xb(k,i) = xa(na,i)+(xa(na-1,i)-xa(na,i))*wt
          END IF
       END DO
    END DO

 END SUBROUTINE htint2d

 !
 !----------------------------------------------------------------------
 ! HSTART:  This subroutine reads a history file and does
 ! a history start
 !
 SUBROUTINE hstart

    USE step, ONLY : time
    USE mpi_interface, ONLY : myid

    IMPLICIT NONE

    CALL read_hist(time, hfilin)

    dtlv = 2.*dtl
    dtlt = dtl

    IF(myid == 0) &
       PRINT "(//' ',49('-')/,' ',/,' History read from: ',A60)",hfilin

    RETURN
 END SUBROUTINE hstart
 !
 !----------------------------------------------------------------------
 ! RANDOM_PERT: initialize field between k=2 and kmx with a
 ! random perturbation of specified magnitude
 !
 SUBROUTINE random_pert(n1,n2,n3,zt,fld,xmag,kmx)

    USE mpi_interface, ONLY : nypg,nxpg,myid,wrxid,wryid,xoffset,yoffset, &
                              double_scalar_par_sum

    USE util, ONLY : sclrset
    IMPLICIT NONE

    INTEGER, INTENT(in) :: n1,n2,n3,kmx
    TYPE(FloatArray3d), INTENT(inout) :: fld
    TYPE(FloatArray1d), INTENT(in) :: zt
    REAL, INTENT(in)    :: xmag(n1)

    REAL (kind=8) :: rand(3:n2-2,3:n3-2),  xx, xxl
    REAL (kind=8), ALLOCATABLE :: rand_temp(:,:)

    INTEGER, DIMENSION (:), ALLOCATABLE :: seed

    INTEGER :: i,j,k,n2g,n3g,isize

    rand = 0.0
    ! seed must be a double precision odd whole number greater than
    ! or equal to 1.0 and less than 2**48.

    CALL random_seed(size=isize)
    ALLOCATE (seed(isize))
    seed(:) = iseed
    CALL random_seed(put=seed)
    DEALLOCATE (seed)
    n2g = nxpg
    n3g = nypg

    DO k = 2, kmx
       ALLOCATE (rand_temp(3:n2g-2,3:n3g-2))
       CALL random_number(rand_temp)
       rand(3:n2-2, 3:n3-2) = rand_temp(3+xoffset(wrxid):n2+xoffset(wrxid)-2, &
                              3+yoffset(wryid):n3+yoffset(wryid)-2)
       DEALLOCATE (rand_temp)

       xx = 0.
       DO j = 3, n3-2
          DO i = 3, n2-2
             fld%d(k,i,j) = fld%d(k,i,j) + rand(i,j)*xmag(k)
             xx = xx + rand(i,j)*xmag(k)
          END DO
       END DO

       xxl = xx
       CALL double_scalar_par_sum(xxl,xx)
       xx = xx/REAL((n2g-4)*(n3g-4))
       fld%d(k,:,:)= fld%d(k,:,:) - xx
    END DO

    IF(myid == 0) THEN
       PRINT *
       PRINT *,'-------------------------------------------------'
       PRINT 600,zt%d(kmx),rand(3,3),xx
       PRINT *,'-------------------------------------------------'
    END IF

    CALL sclrset('cnst',n1,n2,n3,fld%d)

    RETURN

600 FORMAT(2x,'Inserting random temperature perturbations',      &
           /3x,'Below: ',F7.2,' meters;',                        &
           /3x,'with test value of: ',E12.5,                     &
           /3x,'and a magnitude of: ',E12.5)
 END SUBROUTINE random_pert


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
    USE util, ONLY : getMassIndex
    IMPLICIT NONE
    INTEGER :: k,i,j,bb,nc

    DO j = 3, nyp-2
       DO i = 3, nxp-2
          DO k = 2, nzp ! Apply tendencies
             a_naerop%d(k,i,j,:)  = MAX( a_naerop%d(k,i,j,:)  + dtlt*a_naerot%d(k,i,j,:), 0. )
             a_ncloudp%d(k,i,j,:) = MAX( a_ncloudp%d(k,i,j,:) + dtlt*a_ncloudt%d(k,i,j,:), 0. )
             a_nprecpp%d(k,i,j,:) = MAX( a_nprecpp%d(k,i,j,:) + dtlt*a_nprecpt%d(k,i,j,:), 0. )
             a_maerop%d(k,i,j,:)  = MAX( a_maerop%d(k,i,j,:)  + dtlt*a_maerot%d(k,i,j,:), 0. )
             a_mcloudp%d(k,i,j,:) = MAX( a_mcloudp%d(k,i,j,:) + dtlt*a_mcloudt%d(k,i,j,:), 0. )
             a_mprecpp%d(k,i,j,:) = MAX( a_mprecpp%d(k,i,j,:) + dtlt*a_mprecpt%d(k,i,j,:), 0. )
             a_gaerop%d(k,i,j,:)  = MAX( a_gaerop%d(k,i,j,:)  + dtlt*a_gaerot%d(k,i,j,:), 0. )
             a_rp%d(k,i,j) = a_rp%d(k,i,j) + dtlt*a_rt%d(k,i,j)

             IF(level == 5) THEN 
                a_nicep%d(k,i,j,:)   = MAX( a_nicep%d(k,i,j,:)   + dtlt*a_nicet%d(k,i,j,:), 0. )
                a_micep%d(k,i,j,:)   = MAX( a_micep%d(k,i,j,:)   + dtlt*a_micet%d(k,i,j,:), 0. )
             END IF

          END DO
       END DO
    END DO

    nc = spec%getIndex('H2O')
    ! Activation + diagnostic array initialization
    ! Clouds and aerosols
    a_rc%d(:,:,:) = 0.
    DO bb = 1, ncld
       a_rc%d(:,:,:) = a_rc%d(:,:,:) + a_mcloudp%d(:,:,:,getMassIndex(ncld,bb,nc))
    END DO
    DO bb = 1, nbins
       a_rc%d(:,:,:) = a_rc%d(:,:,:) + a_maerop%d(:,:,:,getMassIndex(nbins,bb,nc))
    END DO
    
    ! Ice
    IF ( level == 5 ) THEN
       a_ri%d(:,:,:) = 0.
       a_riri%d(:,:,:) = 0.
       DO bb = 1, nice
          a_ri%d(:,:,:) = a_ri%d(:,:,:) + a_micep%d(:,:,:,getMassIndex(nice,bb,nc))
          a_riri%d(:,:,:) = a_riri%d(:,:,:) + a_micep%d(:,:,:,getMassIndex(nice,bb,nc+1))
       END DO
    END IF
    
 END SUBROUTINE SALSAInit

 ! --------------------------------------------------------------------------------------------------
 ! Replacement for subroutine init_aero_sizedist (init.f90): initilize altitude-dependent aerosol
 ! size distributions and compositions.
 !
 ! Tomi Raatikainen, FMI, 29.2.2016
 !
 SUBROUTINE aerosol_init()

    USE mo_salsa_sizedist, ONLY : size_distribution
    USE mo_salsa_types, ONLY : aero
    USE mo_submctl, ONLY : pi6, nmod, nbins, nspec_dry, in1a,in2a,in2b,fn1a,fn2a,fn2b,  &
                           sigmagA, dpgA, nA, sigmagB, dpgB, nB, volDistA, volDistB, nreg,isdtyp
    USE mpi_interface, ONLY : myid
    USE util, ONLY : getMassIndex

    IMPLICIT NONE
    REAL :: core(nbins), nsectA(1,1,nbins), nsectB(1,1,nbins)   ! Size of the bin mid aerosol particle, local aerosol size dist
    REAL :: pndist(nzp,nbins)                         ! Aerosol size dist as a function of height
    REAL :: pvf2a(nzp,nspec_dry), pvf2b(nzp,nspec_dry)        ! Mass distributions of aerosol species for a and b-bins
    REAL :: pvfOC1a(nzp)                              ! Mass distribution between SO4 and OC in 1a
    INTEGER :: ss,ee,i,j,k
    INTEGER :: iso4 = -1, ioc = -1, ibc = -1, idu = -1, &
               iss = -1, inh = -1, ino = -1

    CHARACTER(len=600) :: fmt = &
       "(/,' -------------------------------------------------',/," // &
       "' Initial aerosol profile: ',//, 4X, 'Height', 6X, 'Na', 9X, 'Nb'," // &
       "7X, 'SO4a', 8X, 'OCa', 9X, 'BCa', 9X, 'DUa', 9X, 'SSa', 9X, 'NH3a', 8X, 'HNO3a'," // &
       "7X, 'SO4b', 8X, 'OCb', 9X, 'BCb', 9X, 'DUb', 9X, 'SSb', 9X, 'NH3b', 8X, 'HNO3b'//," // &
       "(3F10.2,14ES12.3))"
    !
    ! Bin mean aerosol particle volume
    core = 0.
    core(1:nbins) = pi6 * aero(1,1,1:nbins)%dmid**3

    ! Set concentrations to zero
    pndist = 0.
    pvf2a = 0.; pvf2b = 0.
    pvfOC1a = 0.

    a_maerop%d(:,:,:,:) = 0.0
    a_naerop%d(:,:,:,:) = 0.0

    ! Indices (-1 = not used)
    i = 0
    IF (spec%isUsed('SO4')) THEN
       iso4 = spec%getIndex('SO4')
       i = i+1
    END IF
    IF (spec%isUsed('OC')) THEN
       ioc = spec%getIndex('OC')
       i = i+1
    END IF
    IF (spec%isUsed('BC')) THEN
       ibc = spec%getIndex('BC')
       i = i+1
    END IF
    IF (spec%isUsed('DU')) THEN
       idu = spec%getIndex('DU')
       i = i+1
    END IF
    IF (spec%isUsed('SS')) THEN
       iss = spec%getIndex('SS')
       i = i+1
    END IF
    IF (spec%isUsed('NO')) THEN
       ino = spec%getIndex('NO')
       i = i+1
    END IF
    IF (spec%isUsed('NH')) THEN
       inh = spec%getIndex('NH')
       i = i+1
    END IF

    ! All species must be known
    IF (i /= nspec_dry) THEN
       WRITE(*,*) i,nspec_dry
       STOP 'Unknown aerosol species given in the initialization!'
    END IF

    !
    ! Altitude dependent size distributions and compositions.
    ! Read and interpolate/extrapolate size distribution and composition for altitude level k: z=(zt(k)
    ! ---------------------------------------------------------------------------------------------------
    IF (isdtyp == 1) THEN

       CALL READ_AERO_INPUT(pndist,pvfOC1a,pvf2a,pvf2b)

    !
    ! Uniform profiles based on namelist parameters
    ! ---------------------------------------------------------------------------------------------------
    ELSE IF (isdtyp == 0) THEN

       IF (spec%isUsed("OC") .AND. spec%isUsed("SO4")) THEN
          ! Both are there, so use the given "massDistrA"
          pvfOC1a(:) = volDistA(spec%getIndex("OC"))/(volDistA(spec%getIndex("OC"))+volDistA(spec%getIndex("SO4"))) ! Normalize
       ELSE IF (spec%isUsed("OC")) THEN
          ! Pure OC
          pvfOC1a(:) = 1.0
       ELSE IF (spec%isUSed("SO4")) THEN
          ! Pure SO4
          pvfOC1a(:) = 0.0
       ELSE
          STOP 'Either OC or SO4 must be active for aerosol region 1a!'
       END IF

       ! Mass fractions for species in a and b-bins
       DO ss = 1,nspec_dry
          pvf2a(:,ss) = volDistA(ss)
          pvf2b(:,ss) = volDistB(ss)
       END DO

       !
       ! Uniform aerosol size distribution with height.
       ! Using distribution parameters (n, dpg and sigmag) from the SALSA namelist
       !
       ! Convert to SI
       nsectA = 0.
       nsectB = 0.
       nA = nA*1.e6
       nB = nB*1.e6
       dpgA = dpgA*1.e-6
       dpgB = dpgB*1.e-6 
       CALL size_distribution(1,1,1, nmod, in1a, fn2a, nA, dpgA, sigmagA, nsectA)
       CALL size_distribution(1,1,1, nmod, in2b, fn2b, nB, dpgB, sigmagB, nsectB)
       DO ss = 1, nbins
          pndist(:,ss) = nsectA(1,1,ss) + nsectB(1,1,ss)
       END DO

    END IF

    ! ----------------------------------------------------------

    !
    ! Initialize concentrations
    ! ----------------------------------------------------------
    DO k = 2, nzp  ! DONT PUT STUFF INSIDE THE GROUND
       DO j = 1, nyp
          DO i = 1, nxp

             ! a) Number concentrations
             ! Region 1
             a_naerop%d(k,i,j,in1a:fn1a) = pndist(k,in1a:fn1a)

             ! Region 2
             IF (nreg > 1) THEN
                a_naerop%d(k,i,j,in2a:fn2a) = pndist(k,in2a:fn2a)
                a_naerop%d(k,i,j,in2b:fn2b) = pndist(k,in2b:fn2b)
             END IF

             !
             ! b) Aerosol mass concentrations
             ! bin regime 1, done here separately because of the SO4/OC convention
             ! SO4
             IF (spec%isUsed("SO4")) THEN
                ss = getMassIndex(nbins,in1a,spec%getIndex("SO4")); ee = getMassIndex(nbins,fn1a,spec%getIndex("SO4"))
                a_maerop%d(k,i,j,ss:ee) = max(0.0,1.0-pvfOC1a(k))*pndist(k,in1a:fn1a)*core(in1a:fn1a)*spec%rhosu
             END IF
             ! OC
             IF (spec%isUsed("OC")) THEN
                ss = getMassIndex(nbins,in1a,spec%getIndex("OC")); ee = getMassIndex(nbins,fn1a,spec%getIndex("OC"))
                a_maerop%d(k,i,j,ss:ee) = max(0.0,pvfOC1a(k))*pndist(k,in1a:fn1a)*core(in1a:fn1a)*spec%rhooc
             END IF

          END DO ! i
       END DO ! j
    END DO ! k

    !
    ! c) Aerosol mass concentrations
    ! bin regime 2

    IF (nreg > 1) THEN
       DO ss = 1,nspec_dry
          CALL setAeroMass(spec%ind(ss),pvf2a,pvf2b,pndist,core,spec%rholiq(ss))
       END DO
    END IF
       
    ! Put out some info about the initial state
    ! ---------------------------------------------------------------------------------------------------------------------
    IF (myid == 0)                   WRITE(*,*) ''
    IF (myid == 0 .AND. isdtyp == 0) WRITE(*,*) 'AEROSOL PROPERTIES TAKEN FROM A NAMELIST'
    IF (myid == 0 .AND. isdtyp == 1) WRITE(*,*) 'AEROSOL PROPERTIES READ FROM aerosol_in.nc'

    IF (myid == 0) WRITE(*,fmt) &
       ( zt%d(k), SUM(a_naerop%d(k,3,3,in1a:fn2a))*1.e-6, SUM(a_naerop%d(k,3,3,in2b:fn2b))*1.e-6,                 &

       MERGE( SUM( a_maerop%d(k,3,3,MAX(iso4-1,0)*nbins+in1a:MAX(iso4-1,0)*nbins+fn2a) ), -999., iso4>0 ),    &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ioc-1,0)*nbins+in1a:MAX(ioc-1,0)*nbins+fn2a) ), -999., ioc>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ibc-1,0)*nbins+in1a:MAX(ibc-1,0)*nbins+fn2a) ), -999., ibc>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(idu-1,0)*nbins+in1a:MAX(idu-1,0)*nbins+fn2a) ), -999., idu>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(iss-1,0)*nbins+in1a:MAX(iss-1,0)*nbins+fn2a) ), -999., iss>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ino-1,0)*nbins+in1a:MAX(ino-1,0)*nbins+fn2a) ), -999., ino>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(inh-1,0)*nbins+in1a:MAX(inh-1,0)*nbins+fn2a) ), -999., inh>0 ),       &

       MERGE( SUM( a_maerop%d(k,3,3,MAX(iso4-1,0)*nbins+in2b:MAX(iso4-1,0)*nbins+fn2b) ), -999., iso4>0 ),    &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ioc-1,0)*nbins+in2b:MAX(ioc-1,0)*nbins+fn2b) ), -999., ioc>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ibc-1,0)*nbins+in2b:MAX(ibc-1,0)*nbins+fn2b) ), -999., ibc>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(idu-1,0)*nbins+in2b:MAX(idu-1,0)*nbins+fn2b) ), -999., idu>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(iss-1,0)*nbins+in2b:MAX(iss-1,0)*nbins+fn2b) ), -999., iss>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(ino-1,0)*nbins+in2b:MAX(ino-1,0)*nbins+fn2b) ), -999., ino>0 ),       &
       MERGE( SUM( a_maerop%d(k,3,3,MAX(inh-1,0)*nbins+in2b:MAX(inh-1,0)*nbins+fn2b) ), -999., inh>0 ),       &
       k=1,nzp )

 END SUBROUTINE aerosol_init

 !
 ! ----------------------------------------------------------
 ! Sets the mass concentrations to aerosol arrays in 2a and 2b
 !
 !
 SUBROUTINE setAeroMass(ispec,ppvf2a,ppvf2b,ppndist,pcore,prho)
    USE mo_submctl, ONLY : nbins, in2a,fn2a,in2b,fn2b,nspec_dry
    USE util, ONLY : getMassIndex
    
    IMPLICIT NONE
    
    INTEGER, INTENT(in) :: ispec                             ! Aerosol species index
    REAL, INTENT(in) :: ppvf2a(nzp,nspec_dry), ppvf2b(nzp,nspec_dry) ! Mass distributions for a and b bins
    REAL, INTENT(in) :: ppndist(nzp,nbins)                   ! Aerosol size distribution
    REAL, INTENT(in) :: pcore(nbins)                         ! Aerosol bin mid core volume
    REAL, INTENT(in) :: prho                                 ! Aerosol density

    INTEGER :: ss,ee
    INTEGER :: i,j,k

    DO k = 2, nzp ! DONT PUT STUFF INSIDE THE GROUND
       DO j = 1, nyp
          DO i = 1, nxp
             ! 2a
             ss = getMassIndex(nbins,in2a,ispec); ee = getMassIndex(nbins,fn2a,ispec)
             a_maerop%d(k,i,j,ss:ee) =      &
                max( 0.0,ppvf2a(k,ispec) ) * &
                ppndist(k,in2a:fn2a)*pcore(in2a:fn2a)*prho
             ! 2b
             ss = getMassIndex(nbins,in2b,ispec); ee = getMassIndex(nbins,fn2b,ispec)
             a_maerop%d(k,i,j,ss:ee) =      &
                max( 0.0,ppvf2b(k,ispec) ) * &
                ppndist(k,in2b:fn2b)*pcore(in2b:fn2b)*prho
          END DO
       END DO
    END DO

 END SUBROUTINE setAeroMass
 !
 ! -------------------------------------------------------------------------
 ! Reads vertical profiles of aerosol size distribution parameters, aerosol species volume fractions and
 ! number concentration fractions between a and b bins
 !
 SUBROUTINE READ_AERO_INPUT(ppndist,ppvfOC1a,ppvf2a,ppvf2b)
    USE ncio, ONLY : open_aero_nc, read_aero_nc_1d, read_aero_nc_2d, close_nc
    USE mo_submctl, ONLY : nbins, in1a, fn2a, in2b, fn2b,  &
                           nspec_dry, maxspec, nmod
    USE mo_salsa_sizedist, ONLY : size_distribution
    USE mpi_interface, ONLY : appl_abort, myid
    IMPLICIT NONE

    REAL, INTENT(out) :: ppndist(nzp,nbins)                   ! Aerosol size dist as a function of height
    REAL, INTENT(out) :: ppvf2a(nzp,nspec_dry), ppvf2b(nzp,nspec_dry) ! Volume distributions of aerosol species for a and b-bins
    REAL, INTENT(out) :: ppvfOC1a(nzp)                        ! Volume distribution between SO4 and OC in 1a

    REAL :: nsectA(1,1,nbins), nsectB(1,1,nbins)

    INTEGER :: ncid, k, i
    INTEGER :: nc_levs=500, nc_nspec, nc_nmod

    REAL :: pndistA(nzp,nbins), pndistB(nzp,nbins)

    ! Stuff that will be read from the file
    REAL, ALLOCATABLE :: zlevs(:),        &  ! Levels in meters
                         zvolDistA(:,:),  &  ! Volume distribution of aerosol species in a and b bins
                         zvoldistB(:,:),  &  ! (Don't mess these with the ones in namelist.salsa -
                                             !  they are not used here!)
                         znA(:,:),         &  ! Aerosol mode number concentrations, regime A
                         zsigmagA(:,:),    &  ! Geometric standard deviations, regime A
                         zdpgA(:,:),       &  ! Mode mean diameters, regime A
                         znB(:,:),         &  ! number concentration for regime B
                         zsigmagB(:,:),    &  ! getometric std, regime B
                         zdpgB(:,:),       &  ! Mode mean diameter, regime B
                         znsectA(:,:),     &  ! Helper for binned number concentrations regime A
                         znsectB(:,:),     &  ! - '' - regime B
                         helper(:,:)         ! nspec helper
    LOGICAL :: READ_NC

    ! Read the NetCDF input when it is available
    INQUIRE(FILE='aerosol_in.nc',EXIST=READ_NC)

    ! Open the input file
    IF (READ_NC) CALL open_aero_nc(ncid, nc_levs, nc_nspec, nc_nmod)

    ! Check that the input dimensions are compatible with SALSA initialization
    ! ....

    ! Allocate input variables
    ALLOCATE( zlevs(nc_levs),              &
              zvolDistA(nc_levs,maxspec),  &
              zvolDistB(nc_levs,maxspec),  &
              znA(nc_levs,nmod),           &
              zsigmagA(nc_levs,nmod),      &
              zdpgA(nc_levs,nmod),         &
              znB(nc_levs,nmod),           &
              zsigmagB(nc_levs,nmod),      &
              zdpgB(nc_levs,nmod),         &
        ! Couple of helper arrays
              znsectA(nc_levs,nbins),      &
              znsectB(nc_levs,nbins),      &
              helper(nc_levs,nspec_dry)    )

    zlevs = 0.; zvolDistA = 0.; zvolDistB = 0.; znA = 0.; zsigmagA = 0.
    zdpgA = 0.; znB = 0.; zsigmagB = 0.; zdpgB = 0.; znsectA = 0.; znsectB = 0.; helper = 0.

    IF (READ_NC) THEN
       ! Read the aerosol profile data
       CALL read_aero_nc_1d(ncid,'levs',nc_levs,zlevs)
       CALL read_aero_nc_2d(ncid,'volDistA',nc_levs,maxspec,zvolDistA)
       CALL read_aero_nc_2d(ncid,'volDistB',nc_levs,maxspec,zvolDistB)
       CALL read_aero_nc_2d(ncid,'nA',nc_levs,nmod,znA)
       CALL read_aero_nc_2d(ncid,'nB',nc_levs,nmod,znB)
       CALL read_aero_nc_2d(ncid,'dpgA',nc_levs,nmod,zdpgA)
       CALL read_aero_nc_2d(ncid,'dpgB',nc_levs,nmod,zdpgB)
       CALL read_aero_nc_2d(ncid,'sigmagA',nc_levs,nmod,zsigmagA)
       CALL read_aero_nc_2d(ncid,'sigmagB',nc_levs,nmod,zsigmagB)

       CALL close_nc(ncid)
    ELSE
       ! Read the profile data from a text file
       OPEN(11,file='aerosol_in',status='old',form='formatted')
       DO i = 1, nc_levs
          READ(11,*,end=100) zlevs(i)
          READ(11,*,end=100) (zvolDistA(i,k),k=1,nspec_dry) ! Note: reads just "nspec_dry" values from the current line
          READ(11,*,end=100) (zvolDistB(i,k),k=1,nspec_dry) ! -||-
          READ(11,*,end=100) (znA(i,k),k=1,nmod)
          READ(11,*,end=100) (znB(i,k),k=1,nmod)
          READ(11,*,end=100) (zdpgA(i,k),k=1,nmod)
          READ(11,*,end=100) (zdpgB(i,k),k=1,nmod)
          READ(11,*,end=100) (zsigmagA(i,k),k=1,nmod)
          READ(11,*,end=100) (zsigmagB(i,k),k=1,nmod)
       END DO
100 CONTINUE
    CLOSE(11)
    !
    ! The true number of altitude levels
    nc_levs = i-1
 END IF
 !
 IF (zlevs(nc_levs) < zt%d(nzp)) THEN
    IF (myid == 0) PRINT *, '  ABORTING: Model top above aerosol sounding top'
    IF (myid == 0) PRINT '(2F12.2)', zlevs(nc_levs), zt%d(nzp)
    CALL appl_abort(0)
 END IF

 ! Convert to SI
 znA = znA*1.e6
 znB = znB*1.e6
 zdpgA = zdpgA*1.e-6
 zdpgB = zdpgB*1.e-6

 ! Get the binned size distribution
 znsectA = 0.
 znsectB = 0.
 DO k = 1, nc_levs
    CALL size_distribution(1,1,1,nmod,in1a,fn2a,znA(k,:),zdpgA(k,:),zsigmagA(k,:),nsectA)
    CALL size_distribution(1,1,1,nmod,in2b,fn2b,znB(k,:),zdpgB(k,:),zsigmagB(k,:),nsectB)
    znsectA(k,:) = nsectA(1,1,:)
    znsectB(k,:) = nsectB(1,1,:)
 END DO

 ! Interpolate the input variables to model levels
 ! ------------------------------------------------
 CALL htint2d(nc_levs,zvolDistA(1:nc_levs,1:nspec_dry),zlevs(1:nc_levs),nzp,ppvf2a,zt%d,nspec_dry)
 CALL htint2d(nc_levs,zvolDistB(1:nc_levs,1:nspec_dry),zlevs(1:nc_levs),nzp,ppvf2b,zt%d,nspec_dry)
 CALL htint2d(nc_levs,znsectA(1:nc_levs,:),zlevs(1:nc_levs),nzp,pndistA,zt%d,nbins)
 CALL htint2d(nc_levs,znsectB(1:nc_levs,:),zlevs(1:nc_levs),nzp,pndistB,zt%d,nbins)
 ppndist = pndistA + pndistB
 
 ! Since 1a bins by SALSA convention can only contain SO4 or OC,
 ! get renormalized mass fractions.
 ! --------------------------------------------------------------
 IF (spec%isUsed("OC") .AND. spec%isUsed("SO4")) THEN
    ! Both are there, so use the given "massDistrA"
    ppvfOC1a(:) = ppvf2a(:,spec%getIndex("OC"))/(ppvf2a(:,spec%getIndex("OC"))+ppvf2a(:,spec%getIndex("SO4"))) ! Normalize
 ELSE IF (spec%isUsed("OC")) THEN
    ! Pure OC
    ppvfOC1a(:) = 1.0
 ELSE IF (spec%isUsed("SO4")) THEN
    ! Pure SO4
    ppvfOC1a(:) = 0.0
 ELSE
    STOP 'Either OC or SO4 must be active for aerosol region 1a!'
 END IF

 DEALLOCATE( zlevs, zvolDistA, zvolDistB, znA, znB, zsigmagA, zsigmagB, zdpgA, zdpgB, znsectA, znsectB, helper )

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

    ! These could be read from a file
    ! Taken as molecules/kg
    DO j = 1, nyp
       DO i = 1, nxp
          DO k = 1, nzp
             a_gaerop%d(k,i,j,1) = 5.E14/dn0%d(k) !SO4
             a_gaerop%d(k,i,j,2) = 0./dn0%d(k)    !NO3
             a_gaerop%d(k,i,j,3) = 0./dn0%d(k)    !NH4
             a_gaerop%d(k,i,j,4) = 5.E14/dn0%d(k) !OCNV
             a_gaerop%d(k,i,j,5) = 1.E14/dn0%d(k) !OCSV
          END DO
       END DO
    END DO


 END SUBROUTINE init_gas_tracers



 END MODULE init
