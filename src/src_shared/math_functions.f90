MODULE math_functions

  ! Contains a collection of general purpose math functions that can/should be used
  ! by both UCLALES and SALSA. Note that both UCLALES and SALSA have also their dedicated
  ! utility modules that should be used to hold model-specific functions.
  
  REAL, PARAMETER :: pi = 3.14159265358979323846264338327
  
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
      xx = MAX(1.e-20,1.-x**2)
      t1 = LOG(xx)
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
      b = b/sig**2
      f_gauss = (1./a) * EXP(-0.5*b)          

    END FUNCTION f_gauss

    !-------------------------------------------------
    ! Non-normalized frequency from Gaussian distribution
    !
    REAL FUNCTION f_gauss_NN(x,sig,mean)
      REAL, INTENT(in) :: x,sig,mean
      REAL :: b

      b = (x-mean)**2
      b = b/sig**2
      f_gauss_NN = EXP(-0.5*b)
      
    END FUNCTION f_gauss_NN

    
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
    


    
END MODULE math_functions
