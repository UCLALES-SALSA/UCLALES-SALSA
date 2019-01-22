MODULE mo_aux_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : in1a, fn2a, in2b, fn2b, fca, ica, fca, icb, fcb, nprc, nice, aerobins, cloudbins, precpbins, icebins
  IMPLICIT NONE

  ! LES grid displacements and spacings
  TYPE(FloatArray1d), TARGET :: xt, xm, yt, ym, zt, zm, dzt, dzm  
  TYPE(FloatArray1d), TARGET :: aea, aeb, aetot, cla, clb, cltot, prc, ice
  
  ! Initialization profiles
  TYPE(FloatArray1d), TARGET :: u0, v0, pi0, pi1, th0, dn0, rt0

  CONTAINS

    SUBROUTINE setInitialProfiles(nzp)
      INTEGER, INTENT(in) :: nzp

      REAL :: zeros1d(nzp)

      u0 = FloatArray1d(zeros1d,store=.TRUE.)
      v0 = FloatArray1d(zeros1d,store=.TRUE.)
      pi0 = FloatArray1d(zeros1d,store=.TRUE.)
      pi1 = FloatArray1d(zeros1d,store=.TRUE.)
      th0 = FloatArray1d(zeros1d,store=.TRUE.)
      dn0 = FloatArray1d(zeros1d,store=.TRUE.)
      rt0 = FloatArray1d(zeros1d,store=.TRUE.)
            
    END SUBROUTINE setInitialProfiles

    ! ------------------------------------------

    SUBROUTINE setGridSpacings(level,nzp,nxp,nyp)
      INTEGER, INTENT(in) :: level,nzp,nxp,nyp
      REAL :: zeros_z(nzp), zeros_x(nxp), zeros_y(nyp)

      zeros_z = 0.; zeros_x(:) = 0.; zeros_y(:) = 0.
      
      xt = FloatArray1d(zeros_x,store=.TRUE.)
      xm = FloatArray1d(zeros_x,store=.TRUE.)
      yt = FloatArray1d(zeros_y,store=.TRUE.)
      ym = FloatArray1d(zeros_y,store=.TRUE.)
      zt = FloatArray1d(zeros_z,store=.TRUE.)
      zm = FloatArray1d(zeros_z,store=.TRUE.)
      dzt = FloatArray1d(zeros_z,store=.TRUE.)
      dzm = FloatArray1d(zeros_z,store=.TRUE.)      

      aea = FloatArray1d(aerobins(in1a:fn2a),store=.TRUE.)
      aeb = FloatArray1d(aerobins(in2b:fn2b),store=.TRUE.)
      aetot = FloatArray1d(aerobins(in1a:fn2b),store=.TRUE.)
      cla = FloatArray1d(cloudbins(ica%cur:fca%cur),store=.TRUE.)
      clb = FloatArray1d(cloudbins(icb%cur:fcb%cur),store=.TRUE.)
      cltot = FloatArray1d(cloudbins(ica%cur:fcb%cur),store=.TRUE.)
      prc = FloatArray1d(precpbins(1:nprc),store=.TRUE.)
      IF (level == 5) ice = FloatArray1d(icebins(1:nice),store=.TRUE.)     
      
    END SUBROUTINE setGridSpacings
      
  
END MODULE mo_aux_state
