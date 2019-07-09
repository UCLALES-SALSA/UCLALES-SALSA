MODULE mo_ice_shape
  USE mo_submctl, ONLY : spec,pi6,pi
  IMPLICIT NONE
 
  ! -----------------------------------------------
  ! Container for functions describing the
  ! mass-diameter relationships of particles
  ! -----------------------------------------------

  ! Suggested values for shape parameters for non-spherical ice from Khvorostyanov and Curry 2002
  ! for "Crystal with sector-like branches". Can be specified from SALSA namelist. Mass coeffs changed
  ! to SI units. Also make sure to use SI units when specifying these in the namelist. Pretty much all the
  ! articles report these in non-SI units!
  
  ! For Mass:
  REAL, SAVE :: iceShapeAlpha = 15.56999e-3
  REAL, SAVE :: iceShapeBeta = 2.02

  ! For cross sectional area
  REAL, SAVE :: iceShapeGamma = 0.55
  REAL, SAVE :: iceShapeSigma = 1.97

  TYPE t_shape_coeffs
     REAL :: alpha, beta, gamma, sigma     
  END TYPE t_shape_coeffs

  
  CONTAINS
    !
    !
    !
    !
    FUNCTION getDiameter(mpri,mrim,numc)
      REAL, INTENT(in) :: mpri, mrim
      REAL, INTENT(in) :: numc
      REAL :: getDiameter
      
      REAL :: Mtot,     &
              Dth,Mth,      &
              Dgr,Mgr,      &
              Dcr,Mcr,      &
              rho_b,    &
              Fr

      getDiameter = 1.e-9 
      
      IF (numc < 1.e-6) RETURN
                 
      Mtot = ( mpri + mrim )/numc
      
      Fr = MAX( MIN(mrim/(Mtot*numc),0.99),0.01 )
      rho_b = getBulkRho(Fr)  ! "almost" like graupel density...
      Dth = (pi6*spec%rhoic/iceShapeAlpha)**(1./(iceShapeBeta-3.))
      Dgr = (iceShapeAlpha/(pi6*rho_b))**(1./(3.-iceShapeBeta))
      Dcr = ( (1./(1.-Fr))*iceShapeAlpha/(pi6*rho_b) )**(1./(3.-iceShapeBeta))
      Mth = spec%rhoic*pi6*Dth**3
      Mgr = rho_b*pi6*Dgr**3
      Mcr = rho_b*pi6*Dcr**3
         
      IF (Mgr < Mth .OR. Mgr > Mcr) WRITE(*,*) 'ICE SHAPE VAARIN', Dth, Dgr, Dcr
      
      IF ( Mtot < Mth ) THEN
         ! Small spherical particles
         getDiameter = D_spherical(spec%rhoic,Mtot)         
      ELSE IF (Mtot >= Mth .AND. Mtot < Mgr) THEN
         ! Dense non-spherical
         getDiameter = D_nonsphericalIce(Mtot)
      ELSE IF (Mtot >= Mgr .AND. Mtot < Mcr) THEN
         ! Graupel
         getDiameter = D_spherical(rho_b,Mtot)
      ELSE IF (Mtot > Mcr) THEN
         ! Partially rimed crystals
         getDiameter = D_nonsphericalRimed(Fr,Mtot)
      ELSE
         WRITE(*,*) 'ICE SHAPE VAARIN 2'
      END IF
      
    END FUNCTION getDiameter      

    !
    !-------------------------------------------------------------
    ! Function getDcr
    ! Gets the critical size for complete infilling by rime, i.e.
    ! formation of graupel
    
    !
    !-------------------------------------------------------------
    ! Function getPartiallyRimedAlpha
    ! Makes a linear interpolation for the alpha parameter
    ! between alpha1 and alpha2 based on the rime fraction
    ! DONT USE
    !!REAL FUNCTION getPartiallyRimedAlpha(Fr)
    !  REAL, INTENT(in) :: Fr      
    !  REAL :: hlp
    !  hlp = alpha1 + (alpha2-alpha1)*Fr
    !  getPartiallyRimedAlpha = hlp           
    !END FUNCTION getPartiallyRimedAlpha
    !
    !-------------------------------------------------------------
    ! Function getBulkRho
    ! Gets the bulk ice density based on the particle's rime fraction
    ! Note that this simplifies Morrison and Milbrandt 2015 by neglecting
    ! the Eq 17 and just using pristine ice bulk density for "rho_d"
    !
    REAL FUNCTION getBulkRho(Fr)
      REAL, INTENT(in) :: Fr
      getBulkRho = Fr*spec%rhori + (1.-Fr)*spec%rhoic
    END FUNCTION getBulkRho
    !
    !-------------------------------------------------------------
    ! Function D_nonsphericalIce
    ! Get the diameter (~maximum particle dimension) of non-spherical
    ! ice particles using the coefficients for an m-D relationship.
    !
    REAL FUNCTION D_nonsphericalIce(mass)
      REAL, INTENT(in) :: mass
      D_nonsphericalIce = ( mass/iceShapeAlpha )**(1./iceShapeBeta)
    END FUNCTION D_nonsphericalIce
    !
    !-------------------------------------------------------------
    ! Function D_nonsphericalRimed
    ! Same as above but for large rimed particles
    !
    REAL FUNCTION D_nonsphericalRimed(Fr,mass)
      REAL, INTENT(in) :: Fr
      REAL, INTENT(in) :: mass
      REAL :: hlp
      hlp = 1./(1.-Fr)
      D_nonsphericalRimed = ( mass/(hlp*iceShapeAlpha) )**(1./iceShapeBeta)
    END FUNCTION D_nonsphericalRimed
      
    !
    !-------------------------------------------------------------
    ! Function D_spherical
    ! Get the diameter of a regular spherical particle with
    ! given density
    !
    REAL FUNCTION D_spherical(rho,mass)
      REAL, INTENT(in) :: rho
      REAL, INTENT(in) :: mass
      D_spherical = ( mass/(pi6*rho) )**(1./3.)      
    END FUNCTION D_spherical
      
    ! --------------------------

    SUBROUTINE getShapeCoefficients(ishape,mpri,mrim,numc)
      TYPE(t_shape_coeffs), INTENT(out) :: ishape
      REAL, INTENT(in) :: mpri, mrim
      REAL, INTENT(in) :: numc
      
      REAL :: Mtot,     &
              Dth,Mth,      &
              Dgr,Mgr,      &
              Dcr,Mcr,      &
              rho_b,    &
              Fr

      ! Initially, set values for spherical
      ishape%alpha = pi6*spec%rhoic
      ishape%beta = 3.
      ishape%gamma = pi/4.
      ishape%sigma = 2.      
      
      IF (numc < 1.e-6) RETURN
                 
      Mtot = ( mpri + mrim )/numc
      
      Fr = MAX( MIN(mrim/(Mtot*numc),0.99),0.01 )
      rho_b = getBulkRho(Fr)  ! "almost" like graupel density...
      Dth = (pi6*spec%rhoic/iceShapeAlpha)**(1./(iceShapeBeta-3.))
      Dgr = (iceShapeAlpha/(pi6*rho_b))**(1./(3.-iceShapeBeta))
      Dcr = ( (1./(1.-Fr))*iceShapeAlpha/(pi6*rho_b) )**(1./(3.-iceShapeBeta))
      Mth = spec%rhoic*pi6*Dth**3
      Mgr = rho_b*pi6*Dgr**3
      Mcr = rho_b*pi6*Dcr**3

      ! Modify values for non-spherical particles
      IF ( (Mtot >= Mth .AND. Mtot < Mgr) .OR. (Mtot > Mcr) ) THEN
         ishape%alpha = iceShapeAlpha; ishape%beta = iceShapeBeta
         ishape%gamma = iceShapeGamma; ishape%sigma = iceShapeSigma         
      END IF
      
    END SUBROUTINE getShapeCoefficients

    
END MODULE mo_ice_shape
