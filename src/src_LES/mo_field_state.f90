MODULE mo_field_state
  USE grid, ONLY : memsize,varlist_main,varlist_ps,    &
                   level,iradtyp,isgstyp,lbinanl,      &
                   lsalsabbins,nzp,nxp,nyp,nscl,       &
                   a_sclrp,a_sclrt
  USE mo_diag_state, ONLY : setDiagnosticVariables
  USE mo_progn_state, ONLY : setPrognosticVariables
  USE mo_vector_state, ONLY : setVectorVariables
  USE mo_derived_state, ONLY : setDerivedVariables
  USE mo_ps_state, ONLY : setPSVariables
  USE mo_field_types, ONLY : Prog, Vector, Diag, Derived,    &
                             outProg, outVector, outDiag,    &
                             outDerived, outPS, SALSA_tracers_4d
  USE classFieldArray, ONLY : FieldArray
  IMPLICIT NONE
  
  CONTAINS
  
    SUBROUTINE initialize_FieldArrays()

      ! Instantiate the field arrays
      Prog = FieldArray()
      Vector = FieldArray()
      Diag = FieldArray()
      Derived = FieldArray()
      
      IF (level >= 4) &
           SALSA_tracers_4d = FieldArray()

      CALL setPrognosticVariables(a_sclrp,a_sclrt,Prog,varlist_main,memsize,  &
                                  level,isgstyp,lbinanl,nzp,nxp,nyp,nscl      )
      CALL setVectorVariables(Vector,varlist_main,nzp,nxp,nyp)
      CALL setDiagnosticVariables(Diag,varlist_main,memsize,level,iradtyp,nzp,nxp,nyp)      
      CALL setDerivedVariables(Derived,varlist_main,varlist_ps,lsalsabbins,level,nzp,nxp,nyp)
      
      !  Create some subsets of the full FieldArrays. note that these are pointers, not copies!
      IF (level >= 4) &
           CALL Prog%getByGroup("SALSA_4d",SALSA_tracers_4d)      
      CALL Prog%getByOutputstatus(outProg)
      CALL Vector%getByOutputstatus(outVector)
      CALL Diag%getByOutputstatus(outDiag)
      CALL Derived%getByOutputstatus(outDerived)
      
      CALL setPSVariables(outPS,varlist_ps,level,nzp)
      
    END SUBROUTINE initialize_FieldArrays
    
END MODULE mo_field_state
