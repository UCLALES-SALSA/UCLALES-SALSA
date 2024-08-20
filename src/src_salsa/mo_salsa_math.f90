MODULE mo_salsa_math
  USE mo_submctl, ONLY : pi,pi6
  IMPLICIT NONE

  ! -------------------------------------------------------
  ! Contains a selection of math functions needed in SALSA.
  ! Would be good to see if some of these could be replaced
  ! by existing library functions.
  ! 

  CONTAINS

    ! ----------------------------------------------
    ! Get volume from diameter for spherical particle
    !
    REAL FUNCTION V2D(vol)
      REAL, INTENT(in) :: vol
      V2D = (vol/pi6)**(1./3.)
    END FUNCTION V2D

    ! ------------------------------------------------
    ! Get diameter from volume for spherical particle
    !       
    REAL FUNCTION D2V(diam)
      REAL, INTENT(in) :: diam
      D2V = pi6*diam**3
    END FUNCTION D2V

END MODULE mo_salsa_math
