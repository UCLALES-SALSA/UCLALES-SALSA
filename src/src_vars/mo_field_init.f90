MODULE mo_field_init
  USE mpi_interface, ONLY : appl_abort
  USE grid
  USE mo_progn_state
  USE mo_vector_state
  USE mo_diag_state
  USE mo_derived_state
  USE mo_aux_state
  USE mo_ps_state
  USE mo_ts_state
  USE mo_field_state
  USE mo_check_state
  IMPLICIT NONE

  CONTAINS
  
    ! INITIALIZE FIELDARRAY INSTANCES AND VARIABLES
    SUBROUTINE initialize_FieldArrays()
      LOGICAL, ALLOCATABLE :: main_exist(:), ps_exist(:), ts_exist(:)
      INTEGER :: N_main,N_ps,N_ts
      INTEGER :: F_main,F_ps,F_ts
      LOGICAL :: error
      INTEGER :: i
      
      ! Instantiate the field arrays
      Prog = FieldArray()
      Vector = FieldArray()
      Diag = FieldArray()
      Derived = FieldArray()
      PS = FieldArray()
      TS = FieldArray()
      
      IF (level >= 4) &
           SALSA_tracers_4d = FieldArray()
      
      CALL setPrognosticVariables(a_sclrp,a_sclrt,Prog,varlist_main,  &
           level,isgstyp,lpback,nzp,nxp,nyp,nscl)
      CALL setVectorVariables(Vector,varlist_main,nzp,nxp,nyp)
      CALL setDiagnosticVariables(Diag,varlist_main,memsize,level,iradtyp,lpback,nzp,nxp,nyp)      
      CALL setDerivedVariables(Derived,varlist_main,level)
      CALL setPSVariables(PS,varlist_ps,level,lpback)
      CALL setTSVariables(TS,varlist_ts,level,lpback)

      F_main = SIZE(varlist_main)
      F_ps = SIZE(varlist_ps)
      F_ts = SIZE(varlist_ts)
      
      ! Check for faulty output variable list
      ALLOCATE(main_exist(F_main), ps_exist(F_ps),  &
           ts_exist(F_ts))
      main_exist = .FALSE.
      ts_exist = .FALSE.
      ps_exist = .FALSE.
      N_main = COUNT( varlist_main /= '' )
      N_ps = COUNT( varlist_ps /= '' )
      N_ts = COUNT( varlist_ts /= '' )

      CALL checkOutputs(N_main,main_exist,varlist_main,Prog)
      CALL checkOutputs(N_main,main_exist,varlist_main,Diag)
      CALL checkOutputs(N_main,main_exist,varlist_main,Vector)
      CALL checkOutputs(N_main,main_exist,varlist_main,Derived)
      CALL checkOutputs(N_ps,ps_exist,varlist_ps,PS)
      CALL checkOutputs(N_ts,ts_exist,varlist_ts,TS)
      
      IF ( ANY(.NOT. main_exist(1:N_main)) ) THEN
         IF (breakUndefOutput) THEN
            WRITE(*,*) 'ERROR: the following output (main) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_main(1:N_main),MASK=(.NOT. main_exist(1:N_main)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            CALL appl_abort(11)
         ELSE
            WRITE(*,*) 'WARNING: the following output (main) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_main(1:N_main),MASK=(.NOT. main_exist(1:N_main)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            WRITE(*,*)
         END IF
      END IF
      
      IF ( ANY(.NOT. ts_exist(1:N_ts)) ) THEN
         IF (breakUndefOutput) THEN
            WRITE(*,*) 'ERROR: the following output (ts) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_ts(1:N_ts),MASK=(.NOT. ts_exist(1:N_ts)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            CALL appl_abort(11)
         ELSE
            WRITE(*,*) 'WARNING: the following output (ts) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_ts(1:N_ts),MASK=(.NOT. ts_exist(1:N_ts)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            WRITE(*,*)
         END IF
      END IF
      
      IF ( ANY(.NOT. ps_exist(1:N_ps)) ) THEN
         IF (breakUndefOutput) THEN
            WRITE(*,*) 'ERROR: the following output (ps) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_ps(1:N_ps),MASK=(.NOT. ps_exist(1:N_ps)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            CALL appl_abort(11)
         ELSE
            WRITE(*,*) 'WARNING: the following output (ps) variable names are not valid: '
            WRITE(*,'(A50)') PACK(varlist_ps(1:N_ps),MASK=(.NOT. ps_exist(1:N_ps)))
            WRITE(*,*) ''
            WRITE(*,*) 'Check for typos, asking ice variables with level < 5, and so on.'
            WRITE(*,*)
         END IF
      END IF
            
      DEALLOCATE(main_exist,ps_exist,ts_exist)
      
      !  Create some subsets of the full FieldArrays. note that these are pointers, not copies!
      IF (level >= 4) &
           CALL Prog%getByGroup("SALSA_4d",SALSA_tracers_4d)      
      CALL Prog%getByOutputstatus(outProg)
      CALL Vector%getByOutputstatus(outVector)
      CALL Diag%getByOutputstatus(outDiag)
      CALL Derived%getByOutputstatus(outDerived)
      CALL PS%getByOutputstatus(outPS)
      CALL TS%getByOutputstatus(outTS)
      
    END SUBROUTINE initialize_FieldArrays
         
END MODULE mo_field_init
