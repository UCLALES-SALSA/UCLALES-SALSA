MODULE mo_ice_shape
  USE mo_submctl, ONLY : spec
  IMPLICIT NONE
 
  ! -----------------------------------------------
  ! Container for functions describing the
  ! mass-diameter relationships of particles
  ! -----------------------------------------------

  ! Shape parameters for non-spherical ice (from Morrison and Milbrandt 2015, based on Mithcell and Erfani 2014)
  REAL, PARAMETER :: alpha1 = 2.63e-3  ! For unrimed
  REAL, PARAMETER :: alpha2 = 9.88e-3  ! for heavily rimed
  REAL, PARAMETER :: beta = 1.78
  
  CONTAINS

    !
    !
    !
    !
    SUBROUTINE getDiameter(nb,mpri,mrim,numc)
      INTEGER, INTENT(in) :: nb
      REAL, INTENT(in) :: mpri(nb), mrim(nb)
      REAL, INTENT(in) :: numc(nb)
      REAL, INTENT(out) :: diam(nb)
      
      INTEGER :: bb      
      REAL :: Mtot,     &
              Mth,      &
              rho_b,    &
              Fr
      
      diam(:) = 0.
      
      ! Loop over bins
      DO bb = 1,nb
         Mtot = ( mpri(nb) + mrim(nb) )/numc(nb)
         Fr = mrim(nb)/(Mtot*numc(nb))
         rho_b = getBulkRho(Fr)
         Mth = rho_b*pi6*getDth(Fr,rho_b)**3
         
         IF ( Mtot < Mth ) THEN
            ! Small spherical particles
            diam(nb) = D_spherical(rho_b,Mtot)

         ELSE IF (Mtot >= Mth) THEN
            IF ( Fr < 0.8 ) THEN
               diam(nb) = D_nonsphericalIce(Fr,mass)
            ELSE IF (Fr >= 0.8) THEN
               diam(nb) = D_spherical(rho_b,Mtot)
            END IF
         END IF         
      END DO
      
    END SUBROUTINE getDiameter

      
    !
    !-------------------------------------------------------------
    ! Function getDth
    ! Get the maximum threshold diameter for "small spherical"
    ! ice
    ! 
    REAL FUNCTION getDth(Fr,prho)
      REAL, INTENT(in) :: Fr
      REAL, INTENT(in) :: prho
      REAL :: alpha
      alpha = getPartiallyRimedAlpha(Fr)
      getDth = ( (pi6*prho/alpha)*(1.-Fr) )**(1./(beta-3.))           
    END FUNCTION getDth
    !
    !-------------------------------------------------------------
    ! Function getPartiallyRimedAlpha
    ! Makes a linear interpolation for the alpha parameter
    ! between alpha1 and alpha2 based on the rime fraction
    !
    REAL FUNCTION getPartiallyRimedAlpha(Fr)
      REAL, INTENT(in) :: Fr      
      REAL :: hlp
      hlp = alpha1 + (alpha2-alpha1)*Fr
      getPartiallyRimedAlpha = hlp           
    END FUNCTION getPartiallyRimedAlpha
    !
    !-------------------------------------------------------------
    ! Function getBulkRho
    ! Gets the bulk ice density based on the particle's rime fraction
    !
    REAL FUNCTION getBulkRho(Fr)
      REAL, INTENT(in) :: Fr
      bulkRho = Fr*spec%rhori + (1.-Fr)*spec%rhoic
    END FUNCTION getBulkRho
    !
    !-------------------------------------------------------------
    ! Function D_nonsphericalIce
    ! Get the diameter (~maximum particle dimension) of non-spherical
    ! ice particles using the coefficients for an m-D relationship.
    !
    REAL FUNCTION D_nonsphericalIce(Fr,mass)
      REAL, INTENT(in) :: Fr
      REAL, INTENT(in) :: mass
      REAL :: alpha,hlp
      alpha = getPartiallyRimedAlpha(Fr)
      hlp = 1./(1.-Fr)
      D_nonsphericalIce = ( mass/(hlp*alpha) )**(1./beta)
    END FUNCTION D_nonsphericalIce
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
      
      
END MODULE mo_ice_shape
