MODULE emission_main
  use emission_types, ONLY : EmitConfig, EmitSizeDist, EmitType3Config,       &
                             emitModes, emitData, emitType3, nEmissionModes,  &
                             emitPristineIN
  
  USE mo_seasalt_emission

  USE mo_submctl, ONLY : pi6, in1a, fn2a, in2b, fn2b, nbins, spec, pi6, prlim, &
                         ice_theta_dist 

  USE mo_salsa_types, ONLY : aero

  USE mo_salsa_sizedist, ONLY : size_distribution  ! Could this be packaged somehow differently?

  USE mo_aux_state, ONLY : dzt,zt,xt,yt
  USE mo_progn_state, ONLY : a_maerot, a_naerot, a_naerop, a_indefp, a_indeft
  USE mo_diag_state, ONLY : a_dn
  !USE mo_vector_state, ONLY : a_up, a_vp ! needed for the seasalt thing
  USE grid, ONLY: deltax, deltay, deltaz, dtlt, &                  
                  nxp,nyp,nzp, level
    
  USE util, ONLY: smaller, closest, getMassIndex
  USE exceptionHandling, ONLY: errorMessage
  USE mpi_interface, ONLY : myid
  
  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "emission_main"
   
  CONTAINS

  !
  ! -------------------------------------------------------------------
  ! subroutine aerosol_emission:  calls methods to calculate emitted
  !                               aerosols from ground/sea
  !  
  ! Adapted from the original code by Antti Kukkurainen
  ! Juha Tonttila, FMI, 2017
  !
  SUBROUTINE aerosol_emission(time_in)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "aerosol_emission"

    REAL, INTENT(in) :: time_in   ! time in seconds
    LOGICAL :: condition
    INTEGER :: pr
    !Ali, addition of emitType 3  
    TYPE(EmitType3Config), POINTER :: emdT3
    INTEGER :: conditionT3 = 0

    emdT3 => NULL()
    
    ! Loop over all specified emission profiles
    DO pr = 1,nEmissionModes
       ASSOCIATE(emd => emitModes(pr), edt => emitData(pr))
         IF (emd%emitType == 2) THEN
            condition = getCondition(emd,time_in)
            IF (condition) CALL custom_emission(edt,emd)
         END IF
         
         !Ali, addition of emitType 3 
         IF (emd%emitType == 3) THEN
            emdT3 => emitType3(pr)
            conditionT3 = getConditionT3(emdT3,time_in)
            IF (conditionT3 > 0) THEN
               CALL custom_emission_typ3(edt,emd,emdT3,time_in,conditionT3)
            END IF
         END IF
         
       END ASSOCIATE
       ! Ali, addition of emission type
       emdT3 => NULL()
       
    END DO
    
  END SUBROUTINE aerosol_emission

  !
  ! ---------------------------------------------------------------
  ! Simulates the emission of seasalt particles from
  ! an ocean surface as a function of the 10-m wind
  ! speed.
  !
  !SUBROUTINE surface_emission()
  !  IMPLICIT NONE

 !   CHARACTER(len=50), PARAMETER :: name = "surface_emission"
 !   
 !   REAL :: mass_flux(1,nbins) !mass flux at given radius
 !   REAL :: numb_flux(1,nbins)        !number flux at given radius
 !   REAL :: pseaice(1) = 0            !sea ice fraction
 !   REAL :: velo10m_salsa(1,1)        !wind speed

!    INTEGER :: nc, st, en, ii, jj
!    INTEGER :: in, fn
    
!    ! Surface seasalt emissions possible only if sea salt aerosol is used
!    IF (spec%isUsed(spec%nss)) THEN
!       nc=spec%getIndex(spec%nss)
       
!       IF (esrfc%regime == 1) THEN
!          in = in1a
!          fn = fn2a
!       ELSE IF (esrfc%regime == 2) THEN
!          in = in2b
!          fn = fn2b
!       END IF
!       st=(nc-1)*nbins+in
!       en=(nc-1)*nbins+fn
!       
!       DO jj=3,nyp-2
!          DO ii=3,nxp-2
!             
!             velo10m_salsa = SQRT(a_up(2,ii,jj)**2 + a_vp(2,ii,jj)**2)
!             
!             CALL seasalt_emissions_lsce_salsa(1, 1, 1, pseaice, velo10m_salsa, mass_flux, numb_flux)
!             
!             !number of particles + more particles per unit of time * scaling factor [#/kg]
!             a_naerot(2,ii,jj,in:fn) = a_naerot(2,ii,jj,in:fn) + numb_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             !mass + more mass per unit of time * scaling factor [kg/kg]
!             a_maerot(2,ii,jj,st:en) = a_maerot(2,ii,jj,st:en) + mass_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             
!          END DO
!       END DO
!       
!    END IF
!    
!  END SUBROUTINE surface_emission

  !
  ! ----------------------------------------------------------------
  ! Subroutine custom_emissio: "Customized" emission routine, mainly
  !                            for simulating atnhropogenic emissions,
  !                            such as ship or aircraft emissions etc.
  !                            Support for point sources will be included
  !                            soon. Now only does domain-wide emissions
  !                            at specified altitude and time.
  !
  SUBROUTINE custom_emission(edt,emd)
    IMPLICIT NONE
    
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding"
    
    TYPE(EmitSizeDist), INTENT(in) :: edt  ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd   ! Emission configuration instance
    REAL :: hlp1, hlp2  ! helper variables
    
    INTEGER :: i,j,k,bb,ss,mm
    
    IF (myid == 0) THEN
      WRITE(*,*) '========================'
      WRITE(*,*) 'CALCULATING EMISSIONS'
      WRITE(*,*) '========================'
    END IF

    ASSOCIATE( k1 => emd%emitLevMin, k2 => emd%emitLevMax,   &
               x1 => emd%emitXmin, x2 => emd%emitXmax,       &
               y1 => emd%emitYmin, y2 => emd%emitYmax )
    
      DO bb = 1,nbins
         DO j = y1,y2
            DO i = x1,x2
               DO k = k1,k2

                  ! With level 5 using the contact angle distribution for ice nucletion, update the
                  ! IN nucleated fraction assuming emitted particles contain IN and comprise pirstine INP.
                  ! this contributuion can also be switched off using
                  ! --------------------------------------------------------------------------------
                  IF (level == 5 .AND. ice_theta_dist .AND. emitPristineIN) THEN
                     hlp1 = 0.; hlp2 = 0.
                     IF (a_naerop%d(k,i,j,bb) < prlim) THEN
                        ! Check for empty bins with a rather small limit. This should minimize the contact angle for current bin
                        a_indeft%d(k,i,j,bb) = a_indeft%d(k,i,j,bb) - a_indefp%d(k,i,j,bb)/dtlt
                     ELSE
                        ! hlp1 and hlp2 gets the new IN nucleated fraction after emission, i.e. emission should decrease it
                        hlp1 = (a_naerop%d(k,i,j,bb)*a_indefp%d(k,i,j,bb) + edt%numc(bb)*0.*dtlt) ! latter term obv. symbolic...
                        hlp2 = (a_naerop%d(k,i,j,bb)+edt%numc(bb)*dtlt)
                        ! Conver the emission contribution into tendency
                        a_indeft%d(k,i,j,bb) = a_indeft%d(k,i,j,bb) +  &
                             ( (hlp1 / hlp2) - a_indefp%d(k,i,j,bb) ) / dtlt
                     END IF
                  END IF
                  ! --------------------------------------------------------------------------------

                  ! Emission contribution to number concentration
                  a_naerot%d(k,i,j,bb) = a_naerot%d(k,i,j,bb) + edt%numc(bb)

                  ! Contribution to particle composition
                  DO ss = 1,spec%getNSpec(type="wet")
                     mm = getMassIndex(nbins,bb,ss)
                     a_maerot%d(k,i,j,mm) = a_maerot%d(k,i,j,mm) + edt%mass(mm)
                  END DO

               END DO
            END DO
         END DO
      END DO
      
    END ASSOCIATE

  END SUBROUTINE custom_emission
 
 
! -----------------------------------------------------------------------------
! Subroutine custom_emission_typ3:
!
  SUBROUTINE custom_emission_typ3(edt,emd,emdT3,time,conditionT3)
    IMPLICIT NONE
  
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding_typ3"
    INTEGER, INTENT(in) :: conditionT3
    TYPE(EmitSizeDist), INTENT(in) :: edt       ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd         ! Emission configuration instance
    TYPE(EmitType3Config), INTENT(inout) :: emdT3  ! Emission type 3 configuration instance

    REAL :: hlp1, hlp2
    
    INTEGER :: j,bb,ss,mm, xx,yy,zz1,zz2,k
    REAL :: dt, t_str,t_end
    INTEGER :: ind, i_str,i_end, di
    
    IF (myid == 0) THEN
       WRITE(*,*) '========================'
       WRITE(*,*) 'CALCULATING EMISSIONS TYPE 3'
       WRITE(*,*) '========================'
    END IF
    
    ASSOCIATE( ix => emdT3%ix, iy => emdT3%iy, iz => emdT3%iz, t => emdT3%t, np => emdT3%np, &
               t_trac => emdT3%t_trac,t_in => emdT3%t_in, t_out => emdT3%t_out, &
               z_expan_up => emd%z_expan_up, z_expan_dw => emd%z_expan_dw)
      
      t_str = MAX(t_in(conditionT3), t_trac)
      t_end = MIN(t_out(conditionT3), (time + dtlt) )  
      i_str = MAXLOC(t, 1, mask = t <= t_str)
      i_end = MAXLOC(t, 1, mask = t < t_end)
      
      di = i_end - i_str + 1
      
      DO bb = 1,nbins   
         DO j = 1, di
            
            dt  = ( MIN(t_end, t(i_str+j)) - MAX(t_str, t(i_str+j-1)) )/dtlt
            ind = i_str+j-1
            xx = ix(ind); yy = iy(ind)
            zz1 = iz(ind)-z_expan_dw; zz2 = iz(ind)+z_expan_up

            DO k = zz1,zz2
               a_naerot%d(k,xx,yy,bb) = &
                    a_naerot%d(k,xx,yy,bb) + edt%numc(bb) * dt
               ! If level=5 and ice_theta_dist=TRUE, relax the ice nucleation IN "deficit"
               ! ratio used for evolving contact angle by taking number wghted avg and
               ! assuming zero for the emitted population. This is executed regardless of
               ! the emitted species, but the information is only used for ice nucleating
               ! species. Aerosol comes first in the indeft array.
               IF (level == 5 .AND. ice_theta_dist) THEN
                  hlp1 = 0; hlp2 = 0.
                  IF (a_naerop%d(k,xx,yy,bb) < prlim) THEN
                     ! Check for empty bins with a rather small limit. This should minimize the contact angle for current bin
                     a_indeft%d(k,xx,yy,bb) = a_indeft%d(k,xx,yy,bb) - a_indefp%d(k,xx,yy,bb)/dtlt
                  ELSE                  
                     hlp1 = (a_naerop%d(k,xx,yy,bb)* &
                          a_indefp%d(k,xx,yy,bb) + edt%numc(bb)*0.*dt)
                     hlp2 = (a_naerop%d(k,xx,yy,bb) + edt%numc(bb)*dt)
                     
                     a_indeft%d(k,xx,yy,bb) = &
                          a_indeft%d(k,xx,yy,bb) +  &                 
                          ( (hlp1 / hlp2) - a_indefp%d(k,xx,yy,bb) ) / dt
                  END IF
               END IF
                  
               DO ss = 1,spec%getNSpec(type="wet")
                  mm = getMassIndex(nbins,bb,ss)
                  a_maerot%d(k,xx,yy,mm) = &
                       a_maerot%d(k,xx,yy,mm) + edt%mass(mm) * dt
               END DO
            END DO
         END DO
      END DO
      
      t_trac = time + dtlt
      
    END ASSOCIATE
  END SUBROUTINE custom_emission_typ3
  
  ! ----------------------------------------------------------
  
  FUNCTION getCondition(emd,time)
    IMPLICIT NONE
    LOGICAL :: getCondition
    CLASS(EmitConfig), INTENT(in) :: emd
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "getCondition"
    
    getCondition = (                                &
                     emd%start_time <= time    .AND. &
                     emd%end_time > time            &
                   )   
    
  END FUNCTION getCondition
  
  ! -----------------------------------------------------------------------------
  FUNCTION getConditionT3(emdT3,time)
    IMPLICIT NONE
    INTEGER :: getConditionT3
    INTEGER :: N, i
    CLASS(EmitType3Config), INTENT(in) :: emdT3
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "getConditiontyp3"
    
    ASSOCIATE(t_in => emdT3%t_in, t_out => emdT3%t_out)
      
      getConditionT3 = 0
      N = SIZE(t_in)

      DO i =1,N
        IF (t_in(i) <= time .AND. t_out(i) > time) getConditionT3 = i
      END DO
       
    END ASSOCIATE
    
  END FUNCTION getConditionT3
! -----------------------------------------------------------------------------
!
        
END MODULE emission_main
