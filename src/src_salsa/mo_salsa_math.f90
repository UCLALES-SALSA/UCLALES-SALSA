MODULE mo_salsa_math
  USE mo_submctl, ONLY : pi,pi6
  IMPLICIT NONE

  ! -------------------------------------------------------
  ! Contains a selection of math functions needed in SALSA.
  ! Would be good to see if some of these could be replaced
  ! by existing library functions.
  ! 

  CONTAINS

    ! --------------------------------------------------------
    ! Inverse error function - this is the approximation
    ! based on elementrary functions taken from wikipedia...
    !
    REAL FUNCTION erfm1(x)
      REAL, INTENT(in) :: x
      REAL, PARAMETER :: a = 0.140012 
      REAL :: sgn, t1, t2, t3, t4, t5,xx
      
      sgn = 1.
      IF (x < 0.) sgn = -1.
      xx = MAX( MIN( x, 1.-1.e-6  ), -1.+1.e-6 )
      t1 = LOG(1.-xx**2)
      t2 = 2./(pi*a)
      
      t3 = t2 + 0.5*t1
      t4 = t1/a
      t5 = t2 + 0.5*t1
      erfm1 = sgn * SQRT( SQRT(t3**2 - t4) - t5 )
      
    END FUNCTION erfm1

    ! ------------------------------------------------
    ! Normalized frequency from Gaussian distribution
    !
    REAL FUNCTION f_gauss(x,sig,mean)
      REAL, INTENT(in) :: x, sig, mean
      REAL :: a, b
      
      a = SQRT(2.*pi)*sig
      b = (x - mean)**2
      b = b/(2.*sig**2)
      
      f_gauss = (1./a) * EXP(-b)          
    END FUNCTION f_gauss
    

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

    ! -----------------------------------------------
    ! Cumulative lognormal function
    !
    REAL FUNCTION cumlognorm(dg,sigmag,dpart)
      IMPLICIT NONE
      REAL, INTENT(in) :: dg
      REAL, INTENT(in) :: sigmag
      REAL, INTENT(in) :: dpart
      
      REAL :: hlp1,hlp2
      
      hlp1 = ( LOG(dpart) - LOG(dg) )
      hlp2 = SQRT(2.)*LOG(sigmag)
      cumlognorm = 0.5 + 0.5*ERF( hlp1/hlp2 )
      
    END FUNCTION cumlognorm
    
END MODULE mo_salsa_math
