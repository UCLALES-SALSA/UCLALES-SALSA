MODULE mo_aux_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : in1a, fn2a, in2b, fn2b, fca, ica, fca, icb, fcb, nprc, nice, aerobins, cloudbins, precpbins, icebins
  IMPLICIT NONE

  SAVE
  
  ! LES grid displacements and spacings
  TYPE(FloatArray1d), TARGET :: xt, xm, yt, ym, zt, zm, dzt, dzm  
  TYPE(FloatArray1d), TARGET :: aea, aeb, aetot, cla, clb, cltot, prc, ice
  
  ! Initialization profiles
  TYPE(FloatArray1d), TARGET :: u0, v0, pi0, pi1, th0, dn0, rt0

  CONTAINS

    SUBROUTINE setInitialProfiles(BasicState,nzp)
      INTEGER, INTENT(in) :: nzp
      TYPE(FieldArray), INTENT(inout) :: BasicState
      
      REAL :: zeros1d(nzp)

      CLASS(*), POINTER :: pipeline
      
      u0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => u0
      CALL BasicState%newField('u0','basic state u wind','m/s','zt',.FALSE.,pipeline)
      
      v0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => v0
      CALL BasicState%newField('v0','basic state v wind','m/s','zt',.FALSE.,pipeline)
      
      pi0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => pi0
      CALL BasicState%newField('pi0','pi0','Pa','zt',.FALSE.,pipeline)
      
      pi1 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => pi1
      CALL BasicState%newField('pi1','pi1','Pa','zt',.FALSE.,pipeline)
      
      th0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => th0
      CALL BasicState%newField('th0','basic state potential temperature','K','zt',.FALSE.,pipeline)
      
      dn0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => dn0
      CALL BasicState%newField('dn0','basic state air density','kg/m3','zt',.FALSE.,pipeline)
      
      rt0 = FloatArray1d(zeros1d,store=.TRUE.)
      pipeline => rt0
      CALL BasicState%newField('u0','basic state specific humidity','kg/kg','zt',.FALSE.,pipeline)
      
            
    END SUBROUTINE setInitialProfiles

    ! ------------------------------------------

    SUBROUTINE setGridSpacings(Axes,lbinanl,level,nzp,nxp,nyp)
      INTEGER, INTENT(in) :: level,nzp,nxp,nyp
      LOGICAL, INTENT(in) :: lbinanl
      TYPE(FieldArray), INTENT(inout) :: Axes

      REAL :: zeros_z(nzp), zeros_x(nxp), zeros_y(nyp)

      ClASS(*), POINTER :: pipeline
      
      zeros_z = 0.; zeros_x(:) = 0.; zeros_y(:) = 0.
      
      xt = FloatArray1d(zeros_x,store=.TRUE.)
      pipeline => xt
      CALL Axes%newField('xt','xt','m','xt',.TRUE.,pipeline)

      xm = FloatArray1d(zeros_x,store=.TRUE.)
      pipeline => xm
      CALL Axes%newField('xm','xm','m','xm',.TRUE.,pipeline)

      yt = FloatArray1d(zeros_y,store=.TRUE.)
      pipeline => yt
      CALL Axes%newField('yt','yt','m','yt',.TRUE.,pipeline)

      ym = FloatArray1d(zeros_y,store=.TRUE.)
      pipeline => ym
      CALL Axes%newField('ym','ym','m','ym',.TRUE.,pipeline)

      zt = FloatArray1d(zeros_z,store=.TRUE.)
      pipeline => zt
      CALL Axes%newField('zt','zt','m','zt',.TRUE.,pipeline)

      zm = FloatArray1d(zeros_z,store=.TRUE.)
      pipeline => zm
      CALL Axes%newField('zm','zm','m','zm',.TRUE.,pipeline)

      dzt = FloatArray1d(zeros_z,store=.TRUE.)
      pipeline => dzt
      CALL Axes%newField('dzt','dzt','m','dzt',.FALSE.,pipeline)

      dzm = FloatArray1d(zeros_z,store=.TRUE.)
      pipeline => dzt
      CALL Axes%newField('dzt','dzt','m','dzt',.FALSE.,pipeline)

      IF (level >= 4) THEN
         aea = FloatArray1d(aerobins(in1a:fn2a),store=.TRUE.)
         pipeline => aea
         CALL Axes%newField('aea','Aerosol A lower limits','m','aea',lbinanl,pipeline)
         
         aeb = FloatArray1d(aerobins(in2b:fn2b),store=.TRUE.)
         pipeline => aeb
         CALL Axes%newField('aeb','Aerosol B lower limits','m','aeb',lbinanl,pipeline)
         
         aetot = FloatArray1d(aerobins(in1a:fn2b),store=.TRUE.)
         pipeline => aetot
         CALL Axes%newField('aetot','All aerosol lower limits','m','aetot',.FALSE.,pipeline)
         
         cla = FloatArray1d(cloudbins(ica%cur:fca%cur),store=.TRUE.)
         pipeline => cla
         CALL Axes%newField('cla','Cloud A lower limits','m','cla',lbinanl,pipeline)
         
         clb = FloatArray1d(cloudbins(icb%cur:fcb%cur),store=.TRUE.)
         pipeline => clb
         CALL Axes%newField('clb','Cloud B lower limits','m','clb',lbinanl,pipeline)
         
         cltot = FloatArray1d(cloudbins(ica%cur:fcb%cur),store=.TRUE.)
         pipeline => cltot
         CALL Axes%newField('cltot','All cloud lower limits','m','cltot',.FALSE.,pipeline)
         
         prc = FloatArray1d(precpbins(1:nprc),store=.TRUE.)
         pipeline => prc
         CALL Axes%newField('prc','Precip lower limits','m','prc',lbinanl,pipeline)
         
         IF (level == 5) THEN
            ice = FloatArray1d(icebins(1:nice),store=.TRUE.)     
            pipeline => ice
            CALL Axes%newField('ice','Ice lower limits','m','ice',lbinanl,pipeline)
         END IF
      END IF
         
    END SUBROUTINE setGridSpacings
      
  
END MODULE mo_aux_state
