!****************************************************************
!*                                                            *
!*   module MO_SALSA_PROPERTIES                                 *
!*                                                              *
!*   Contains subroutines and functions that are used           *
!*   to calculate particle properties during simulation         *
!*                                                              *
!****************************************************************
MODULE mo_salsa_properties


CONTAINS

  ! ****************************************************************
  ! Kohler-theory based function for calculating equilibrium droplet water content
  ! for SALSA aerosol bins. Saturation ratio limited between 0.1 and 0.98.
  !
  ! Tomi Raatikainen (FMI) 2018
  SUBROUTINE equilibration(kbdim, klev, Sw, temp, paero, init)

    USE mo_submctl, ONLY : &
         t_section,    &
         in1a, fn1a, fn2b, &
         pi6, rg ,surfw0, nlim, &
         diss, mws, dens, rhowa, mwa

    IMPLICIT NONE

    !-- input variables -------------
    INTEGER, INTENT(in) ::    &
         kbdim,               & ! dimension for arrays
         klev                   ! number of vertical levels

    REAL, INTENT(in) ::       &
         Sw(kbdim,klev),      & ! equilibrium saturation ratio [-]
         temp(kbdim,klev)       ! temperature [K]

    LOGICAL, INTENT(in) :: init  ! TRUE: update water content for 2a and 2b bins in addition to 1a

    !-- inpu/output variables -------------
    TYPE(t_section), INTENT(inout) :: paero(kbdim,klev,fn2b) ! SALSA aerosol


    !-- local variables --------------
    INTEGER :: ii, jj, kk, last, n

    REAL ::      &
         zke,    &   ! Kelvin term [-]
         ns,     &   ! Solute molar concentration [mol]
         zcore,  &   ! Solute volume [m^3]
         zw,     &   ! Water volume [m^3]
         zdwet,  &   ! Droplet diameter [m]
         zdold,  &   ! Old droplet diameter [m]
         zrh         ! Current saturation ratio [-]

    ! 1a always, but 2a and 2b bins only when init = .TRUE.
    last = fn1a
    IF (init) last = fn2b

    ! Loop over all aerosol bins
    DO kk = in1a,last
      DO jj = 1,klev      ! vertical grid
         DO ii = 1,kbdim ! horizontal grid

            IF (paero(ii,jj,kk)%numc > nlim) THEN

               ! -- Saturation ratio (limited between 0.1 and 0.98)
               zrh = MAX(0.1,MIN(Sw(ii,jj),0.98))

               !-- total volume of solutes and water in one particle (m^3)
               zcore = sum(paero(ii,jj,kk)%volc(2:))/paero(ii,jj,kk)%numc
               zw = paero(ii,jj,kk)%volc(1)/paero(ii,jj,kk)%numc

               !-- total dissolved solute molar concentration in one particle (mol)
               ns = SUM( paero(ii,jj,kk)%volc(2:)*diss(2:)*dens(2:)/mws(2:) )/paero(ii,jj,kk)%numc

               zdold = 1.
               DO n=1,1000
                  !-- particle wet radius [m]
                  zdwet = ( (zw+zcore)/pi6 )**(1./3.)

                  !-- stop iteration when droplet diameter is not changing
                  IF (abs(zdwet/zdold-1.) < 1.e-12) EXIT
                  zdold = max(zdwet,1.e-20)

                  !-- Kelvin effect
                  zke = exp(4.*surfw0*mwa/(rhowa*rg*temp(ii,jj)*zdwet))

                  !-- calculate liquid water volume in one particle (m^3)
                  ! zrh = xw*zke = nw/(nw+ns)*zke <=> nw = zrh/(zke-zrh)*ns
                  zw = mwa/rhowa*zrh/(zke-zrh)*ns

               END DO
			   if(zw < 0.)THEN
                  WRITE(*,*)'  bin ',kk,', parameters:',zcore, ns, zw, zdwet, zke
                  STOP 'SALSA equilibration: negative water in aerosol!!'
               ENDIF
               IF (n > 1000) THEN
                  WRITE(*,*)'  bin ',kk,', parameters:',zcore, ns, zw, zdwet, zke
                  STOP 'SALSA equilibration: no convergence!!'
               ENDIF

               ! Update water volume concentration
               paero(ii,jj,kk)%volc(1) = zw*paero(ii,jj,kk)%numc

            END IF

         END DO
      END DO
    END DO

  END SUBROUTINE equilibration


END MODULE mo_salsa_properties
