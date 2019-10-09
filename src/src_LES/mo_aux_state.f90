MODULE mo_aux_state
  USE classFieldArray
  USE mo_structured_datatypes, ONLY : FloatArray1d, FloatArray2d, FloatArray3d, FloatArray4d
  USE mo_submctl, ONLY : in1a, fn2a, in2b, fn2b, fca, ica, fca, icb, fcb,  &
                         nbins, ncld, nprc, nice
  IMPLICIT NONE

  SAVE
  
  ! LES grid displacements and spacings
  TYPE(FloatArray1d), TARGET :: xt, xm, yt, ym, zt, zm, dzt, dzm  
  TYPE(FloatArray1d), TARGET :: aea, aeb, aetot, cla, clb, cltot, prc, ice
  REAL, ALLOCATABLE, TARGET :: a_grid(:)
  
  ! Initialization profiles
  TYPE(FloatArray1d), TARGET :: u0, v0, pi0, pi1, th0, dn0, rt0
  REAL, ALLOCATABLE, TARGET :: a_basicstate(:)

  
  CONTAINS

    SUBROUTINE setInitialProfiles(BasicState,nzp)
      INTEGER, INTENT(in) :: nzp
      TYPE(FieldArray), INTENT(inout) :: BasicState
      
      CLASS(*), POINTER :: pipeline => NULL()
      INTEGER :: nvar, n1,n2
      
      ! Allocate the source array. Remember to extend if adding new variables!
      nvar = 7*nzp
      ALLOCATE(a_basicstate(nvar))
      a_basicstate = 0.
      
      n1 = 1; n2 = nzp
      pipeline => NULL()
      u0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => u0
      CALL BasicState%newField('u0','basic state u wind','m/s','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      v0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => v0
      CALL BasicState%newField('v0','basic state v wind','m/s','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      pi0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => pi0
      CALL BasicState%newField('pi0','pi0','Pa','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      pi1 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => pi1
      CALL BasicState%newField('pi1','pi1','Pa','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      th0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => th0
      CALL BasicState%newField('th0','basic state potential temperature','K','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      dn0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => dn0
      CALL BasicState%newField('dn0','basic state air density','kg/m3','zt',.FALSE.,pipeline)

      n1 = n1+nzp; n2 = n2+nzp
      pipeline => NULL()
      rt0 = FloatArray1d(a_basicstate(n1:n2))
      pipeline => rt0
      CALL BasicState%newField('u0','basic state specific humidity','kg/kg','zt',.FALSE.,pipeline)

      pipeline => NULL()
            
    END SUBROUTINE setInitialProfiles

    ! ------------------------------------------

    SUBROUTINE setGridSpacings(Axes,lbinanl,lsalsabbins,level,nzp,nxp,nyp)
      INTEGER, INTENT(in) :: level,nzp,nxp,nyp
      LOGICAL, INTENT(in) :: lbinanl, lsalsabbins
      TYPE(FieldArray), INTENT(inout) :: Axes

      ClASS(*), POINTER :: pipeline => NULL()
      INTEGER :: nvar, n1,n2,n1aux,n2aux

      ! Allocate the source array, remember to extend if adding new variables
      nvar = 2*nxp + 2*nyp + 4*nzp
      IF (level >= 4) &
           nvar = nvar + nbins + ncld + nprc
      IF (level > 4) &
           nvar = nvar + nice
      ALLOCATE(a_grid(nvar))
      a_grid(:) = 0.
      
      n1 = 1; n2 = nxp
      pipeline => NULL()
      xt = FloatArray1d(a_grid(n1:n2))
      pipeline => xt
      CALL Axes%newField('xt','xt','m','xt',.TRUE.,pipeline)

      n1 = n2 + 1; n2 = n2 + nxp
      pipeline => NULL()
      xm = FloatArray1d(a_grid(n1:n2))
      pipeline => xm
      CALL Axes%newField('xm','xm','m','xm',.TRUE.,pipeline)

      n1 = n2 + 1; n2 = n2 + nyp
      pipeline => NULL()
      yt = FloatArray1d(a_grid(n1:n2))
      pipeline => yt
      CALL Axes%newField('yt','yt','m','yt',.TRUE.,pipeline)

      n1 = n2 + 1; n2 = n2 + nyp
      pipeline => NULL()
      ym = FloatArray1d(a_grid(n1:n2))
      pipeline => ym
      CALL Axes%newField('ym','ym','m','ym',.TRUE.,pipeline)

      n1 = n2 + 1; n2 = n2 + nzp
      pipeline => NULL()
      zt = FloatArray1d(a_grid(n1:n2))
      pipeline => zt
      CALL Axes%newField('zt','zt','m','zt',.TRUE.,pipeline,in_group=["ps"])

      n1 = n2 + 1; n2 = n2 + nzp
      pipeline => NULL()
      zm = FloatArray1d(a_grid(n1:n2))
      pipeline => zm
      CALL Axes%newField('zm','zm','m','zm',.TRUE.,pipeline,in_group=["ps"])

      n1 = n2 + 1; n2 = n2 + nzp
      pipeline => NULL()
      dzt = FloatArray1d(a_grid(n1:n2))
      pipeline => dzt
      CALL Axes%newField('dzt','dzt','m','dzt',.FALSE.,pipeline)

      n1 = n2 + 1; n2 = n2 + nzp
      pipeline => NULL()
      dzm = FloatArray1d(a_grid(n1:n2))
      pipeline => dzm
      CALL Axes%newField('dzm','dzm','m','dzm',.FALSE.,pipeline)

      IF (level >= 4) THEN
         n1 = n2 + 1; n2 = n2 + fn2a
         n1aux = n1
         pipeline => NULL()
         aea = FloatArray1d(a_grid(n1:n2))
         pipeline => aea
         CALL Axes%newField('aea','Aerosol A lower limits','m','aea',lbinanl,pipeline,in_group=["ps","ts"])

         n1 = n2 + 1; n2 = n2 + (fn2b-fn2a)
         n2aux = n2
         pipeline => NULL()
         aeb = FloatArray1d(a_grid(n1:n2))
         pipeline => aeb
         CALL Axes%newField('aeb','Aerosol B lower limits','m','aeb',(lbinanl .AND. lsalsabbins),pipeline,in_group=["ps","ts"])

         pipeline => NULL()
         aetot = FloatArray1d(a_grid(n1aux:n2aux))
         pipeline => aetot
         CALL Axes%newField('aetot','All aerosol lower limits','m','aetot',.FALSE.,pipeline)

         n1 = n2 + 1; n2 = n2 + fca%cur
         n1aux = n1
         pipeline => NULL()
         cla = FloatArray1d(a_grid(n1:n2))
         pipeline => cla
         CALL Axes%newField('cla','Cloud A lower limits','m','cla',lbinanl,pipeline,in_group=["ps","ts"])

         n1 = n2 + 1; n2 = n2 + (fcb%cur-fca%cur)
         n2aux = n2
         pipeline => NULL()
         clb = FloatArray1d(a_grid(n1:n2))
         pipeline => clb
         CALL Axes%newField('clb','Cloud B lower limits','m','clb',(lbinanl .AND. lsalsabbins),pipeline,in_group=["ps","ts"])

         pipeline => NULL()
         cltot = FloatArray1d(a_grid(n1aux:n2aux))
         pipeline => cltot
         CALL Axes%newField('cltot','All cloud lower limits','m','cltot',.FALSE.,pipeline)

         n1 = n2 + 1; n2 = n2 + nprc
         pipeline => NULL()
         prc = FloatArray1d(a_grid(n1:n2))
         pipeline => prc
         CALL Axes%newField('prc','Precip lower limits','m','prc',lbinanl,pipeline,in_group=["ps","ts"])
         
         IF (level == 5) THEN
            n1 = n2 + 1; n2 = n2 + nice
            pipeline => NULL()
            ice = FloatArray1d(a_grid(n1:n2))     
            pipeline => ice
            CALL Axes%newField('ice','Ice lower limits','m','ice',lbinanl,pipeline,in_group=["ps","ts"])
         END IF
      END IF

      pipeline => NULL()
         
    END SUBROUTINE setGridSpacings
      
  
END MODULE mo_aux_state
