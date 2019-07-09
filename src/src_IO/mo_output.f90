MODULE mo_output
  USE netcdf
  USE mpi_interface, ONLY : myid, mpiroot, ver, author, info
  USE ncio
  USE grid, ONLY : outAxes,outAxesPS,outAxesTS,expnme,nzp,nxp,nyp,filprf,  &
                   lbinanl,level,lsalsabbins
  USE mo_field_state, ONLY : outProg, outVector, outDiag, outDerived, outPS, outTS
  USE classFieldArray, ONLY : FieldArray
  USE mo_structured_datatypes
  USE mo_submctl, ONLY : in1a,fn2a,in2b,fn2b,ica,fca,icb,fcb,nprc,nice
  
  IMPLICIT NONE

  ! NOTE ABOUT TS AND PS FILES: These are only written by the rank mpiroot. Initialization and closing calls should be
  ! done only for myid = mpiroot. Writing procedures however must be called by all processes, because the onDemand functions
  ! for the statistics use MPI's collective reduction procedures. However, only mpiroot will actually write the output, i.e.
  ! if the variable dimensions imply statistical output, only mpiroot will call the ncio.f90/write_nc subroutine.

  
  LOGICAL :: tsflg = .FALSE., psflg = .FALSE. ! Flags to accumulate data for statistical processes from within physics routines; should coincide with the statistics output timestep
  REAL :: ps_intvl = 120.
  REAL :: ts_intvl = 120.
  REAL :: main_intvl = 3600.
  INTEGER, PRIVATE :: ncid_main, ncid_ps, ncid_ts
  INTEGER, PRIVATE :: nrec_main, nvar_main, nrec_ps, nvar_ps, nrec_ts, nvar_ts
  CHARACTER(len=150), PRIVATE :: fname_main, fname_ps, fname_ts      
  
  CONTAINS

    !
    ! ----------------------------------------------------------------------
    ! Subroutine init_main:  Defines the netcdf Analysis file
    !
    ! Modified for level 4.
    ! Juha Tonttila, FMI, 2014
    !
    !
    SUBROUTINE init_main(time)
      REAL, INTENT(in) :: time

      INTEGER :: npoints
      
      npoints = (nxp-4)*(nyp-4)

      fname_main = trim(filprf)

      IF ( .NOT. ANY([outProg%Initialized,     &
                      outDiag%Initialized,     &
                      outDerived%Initialized]  &
                    ) ) RETURN ! If no output variables defined, do not open the files 
      
      IF(myid == 0) WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname_main)
      
      CALL open_nc(fname_main,expnme,time,npoints,ncid_main,nrec_main,ver,author,info)
      IF (level < 4 .OR. .NOT. lbinanl) THEN
         CALL define_nc(ncid_main,nrec_main,nvar_main,          &
                        outProg=outProg,outVector=outVector,    &
                        outDiag=outDiag,outDerived=outDerived,  &
                        outAxes=outAxes, n1=nzp,n2=nxp-4,n3=nyp-4)

      ELSE IF (level == 4 .AND. lbinanl) THEN
         IF (lsalsabbins) THEN           
            CALL define_nc(ncid_main,nrec_main,nvar_main,          &
                           outProg=outProg,outVector=outVector,    &
                           outDiag=outDiag,outDerived=outDerived,  &
                           outAxes=outAxes, n1=nzp,n2=nxp-4,         &
                           n3=nyp-4,inae_a=fn2a,inae_b=fn2b-fn2a,    &
                           incld_a=fca%cur,incld_b=fcb%cur-fca%cur,  &
                           inprc=nprc                              )
               
         ELSE              
            CALL define_nc(ncid_main,nrec_main,nvar_main,          &
                           outProg=outProg,outVector=outVector,    &
                           outDiag=outDiag,outDerived=outDerived,  &
                           outAxes=outAxes,n1=nzp,n2=nxp-4,        &
                           n3=nyp-4,inae_a=fn2a,incld_a=fca%cur,   &
                           inprc=nprc                              )                

         END IF
            
      ELSE IF (level == 5 .AND. lbinanl) THEN
         IF (lsalsabbins) THEN
            CALL define_nc(ncid_main,nrec_main,nvar_main,          &
                           outProg=outProg,outVector=outVector,    &
                           outDiag=outDiag,outDerived=outDerived,  &
                           outAxes=outAxes,n1=nzp,n2=nxp-4,        &
                           n3=nyp-4,inae_a=fn2a,inae_b=fn2b-fn2a,  &
                           incld_a=fca%cur,incld_b=fcb%cur-fca%cur,        &                              
                           inprc=nprc,inice=nice                   )
              
         ELSE
            CALL define_nc(ncid_main,nrec_main,nvar_main,          &
                           outProg=outProg,outVector=outVector,    &
                           outDiag=outDiag,outDerived=outDerived,  &
                           outAxes=outAxes,n1=nzp,n2=nxp-4,        &
                           n3=nyp-4,inae_a=fn2a,incld_a=fca%cur,   &                              
                           inprc=nprc,inice=nice                   )
         END IF
         
      END IF
      IF (myid == 0) WRITE(*,*) '   ...starting record: ', nrec_main

      
    END SUBROUTINE init_main

    ! ------------------------------------------------------------------------

    SUBROUTINE init_ps(time)
      REAL, INTENT(in) :: time
      INTEGER :: npoints

      npoints = (nxp-4)*(nyp-4)
      
      fname_ps = TRIM(filprf)//'.ps'

      IF (.NOT. outPS%Initialized) RETURN ! If no variables defined for output, do not create the file 
      
      IF(myid == 0) WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname_ps)

      CALL open_nc(fname_ps,expnme,time,npoints,ncid_ps,nrec_ps,ver,author,info)

      IF (level < 4 .OR. .NOT. lbinanl) THEN
         CALL define_nc(ncid_ps,nrec_ps,nvar_ps,         &
                        outPS=outPS,outAxes=outAxesPS,   &  ! Use the PS subset of the axis variables                               
                        n1=nzp                           )
      END IF

      IF (myid == 0) WRITE(*,*) '   ...starting record: ', nrec_ps
      
    END SUBROUTINE init_ps

    ! ------------------------------------------------------------------------

    SUBROUTINE init_ts(time)
      REAL, INTENT(in) :: time
      INTEGER :: npoints

      npoints = (nxp-4)*(nyp-4)
      
      fname_ts = TRIM(filprf)//'.ts'

      IF (.NOT. outTS%Initialized) RETURN ! If no variables defined for output, do not create the file 
      
      IF(myid == 0) WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname_ts)

      CALL open_nc(fname_ts,expnme,time,npoints,ncid_ts,nrec_ts,ver,author,info)

      IF (level < 4 .OR. .NOT. lbinanl) THEN
         CALL define_nc(ncid_ts,nrec_ts,nvar_ts,         &
                        outTS=outTS,outAxes=outAxesTS    ) ! Use the TS subset of the axis variables                               
      END IF

      IF (myid == 0) write(*,*) '   ...starting record: ', nrec_ts
      
    END SUBROUTINE init_ts
    
    ! -------------------------------------------------------------------------
    
    SUBROUTINE close_main()
      CALL close_nc(ncid_main)
    END SUBROUTINE close_main

    SUBROUTINE close_ps()
      CALL close_nc(ncid_ps)
    END SUBROUTINE close_ps

    SUBROUTINE close_ts()
      CALL close_nc(ncid_ts)
    END SUBROUTINE close_ts
    
    !
    ! ----------------------------------------------------------------------
    ! Subroutine Write_main:  Writes the netcdf Analysis file
    !
    ! Modified for levels 4 and 5
    ! Juha Tonttila, FMI, 2014
    !
    !
    SUBROUTINE write_main(time)
      REAL, INTENT(in) :: time
      INTEGER :: ibeg0(1)

      ibeg0 = [nrec_main]

      IF ( .NOT. ANY([outProg%Initialized,     &
                      outDiag%Initialized,     &
                      outDerived%Initialized]  &
                    ) ) RETURN ! If no output variables defined, do not try to write
      
      ! write time
      CALL write_nc(ncid_main,'time',time,ibeg0)
      
      IF (nrec_main == 1) THEN
         ! First entry -> write axis variables
         CALL write_output(outAxes,nrec_main,ncid_main)
      END IF
      
      IF (outProg%Initialized) &
           CALL write_output(outProg,nrec_main,ncid_main)
      IF (outVector%Initialized) &
           CALL write_output(outVector,nrec_main,ncid_main)
      IF (outDiag%Initialized) &
           CALL write_output(outDiag,nrec_main,ncid_main)
      IF (outDerived%Initialized) &
           CALL write_output(outDerived,nrec_main,ncid_main)

      IF (myid == 0) WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")    &
         nrec_main,fname_main

      CALL sync_nc(ncid_main)
      nrec_main = nrec_main+1      
      
    END SUBROUTINE write_main

    ! --------------------------------------------------------------------------

    SUBROUTINE write_ps(time)   
      REAL, INTENT(in) :: time
      INTEGER :: ibeg0(1)

      ibeg0 = [nrec_ps]

      IF ( .NOT. outPS%Initialized ) RETURN
      
      ! write time
      IF (myid == mpiroot) &
           CALL write_nc(ncid_ps,'time',time,ibeg0)

      IF (nrec_ps == 1) THEN
         ! First entry -> write axis variables
         CALL write_output(outAxesPS,nrec_ps,ncid_ps)
      END IF

      IF (outPS%Initialized) &
           CALL write_output(outPS,nrec_ps,ncid_ps)

      IF (myid == 0) WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")     &
         nrec_ps,fname_ps      

      IF (myid == mpiroot) &
           CALL sync_nc(ncid_ps)
      nrec_ps = nrec_ps+1
      
    END SUBROUTINE write_ps

    ! ----------------------------------------------------------------------------

    SUBROUTINE write_ts(time)
      REAL, INTENT(in) :: time
      INTEGER :: ibeg0(1)

      ibeg0 = [nrec_ts]

      IF ( .NOT. outTS%Initialized ) RETURN

      ! Write time
      IF (myid == mpiroot) &
           CALL write_nc(ncid_ts,'time',time,ibeg0)

      IF (nrec_ts == 1) THEN
         ! First entry -> write axis variables
         CALL write_output(outAxesTS,nrec_ts,ncid_ts)
      END IF

      IF (outTS%Initialized)  &
           CALL write_output(outTS,nrec_ts,ncid_ts)

      IF (myid == 0) WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")    &
         nrec_ts,fname_ts   

      IF (myid == mpiroot) &
           CALL sync_nc(ncid_ts)
      nrec_ts = nrec_ts+1
      
    END SUBROUTINE write_ts    
    !
    ! Subroutine WRITE_OUTPUT: A general wrap-around routine used to to the writing of all kinds of
    !                          output variables
    !
    !
    SUBROUTINE write_output(varArray,nrec0,ncid0)
      TYPE(FieldArray), INTENT(in) :: varArray
      INTEGER, INTENT(in) :: ncid0
      INTEGER, INTENT(inout) :: nrec0
      
      INTEGER :: icnt0dsd(2), icnt1d(2), icnt1dsd(3), icnt2d(3), icnt3d(4), icnt3dsd(5) ! count arrays for 3d variables and 4d size distribution variables
      INTEGER :: ibeg0d(1), ibeg1d(2), ibeg2d(3), ibeg3d(4), ibeg4d(5) ! same for beginning indices

      CHARACTER(len=50) :: vname

      INTEGER :: n, nvar

      TYPE(FloatArray0d), POINTER :: var0d => NULL()
      TYPE(FloatArray1d), POINTER :: var1d => NULL()
      TYPE(FloatArray2d), POINTER :: var2d => NULL()
      TYPE(FloatArray3d), POINTER :: var3d => NULL()
      TYPE(FloatArray4d), POINTER :: var4d => NULL()

      REAL :: out2d(nxp,nyp),out3d(nzp,nxp,nyp)
      REAL :: out1d(nzp)
      REAL :: out0d
      REAL, ALLOCATABLE :: out3dsd(:,:,:,:), outsd(:), out1dsd(:,:)
      
      INTEGER :: i1,i2,j1,j2,nstr,nend

      IF ( .NOT. varArray%Initialized) RETURN  ! No variables assigned
      
      i1 = 3; i2 = nxp-2
      j1 = 3; j2 = nyp-2
      
      ! for 3d outputs these are always the same
      icnt1d = [nzp,1]
      icnt2d = [nxp-4,nyp-4,1]
      icnt3d = [nzp,nxp-4,nyp-4,1]      
      !icntXsd has to be defined on a case by case basis!

      ibeg1d = [1,nrec0]
      ibeg2d = [1,1,nrec0]
      ibeg3d = [1,1,1,nrec0]
      ibeg0d = [nrec0]
      ibeg4d = [1,1,1,1,nrec0]
      
      nvar = varArray%count
      DO n = 1,nvar
         vname = varArray%list(n)%name
         SELECT CASE(varArray%list(n)%dimension)
         CASE('time')
            CALL varArray%getData(1,var0d,index=n)
            IF (ASSOCIATED(var0d%onDemand)) THEN
               CALL var0d%onDemand(vname,out0d)
            ELSE
               out0d = var0d%d
            END IF
            IF (myid == mpiroot) &
                 CALL write_nc(ncid0,vname,out0d,ibeg0d)
            
         CASE('zt','zm')
            CALL varArray%getData(1,var1d,index=n)
            CALL write_nc(ncid0,vname,var1d%d(:),ibeg0d)
            
         CASE('xt','xm')
            CALL varArray%getData(1,var1d,index=n)
            CALL write_nc(ncid0,vname,var1d%d(i1:i2),ibeg0d)
            
         CASE('yt','ym')
            CALL varArray%getData(1,var1d,index=n)
            CALL write_nc(ncid0,vname,var1d%d(j1:j2),ibeg0d)

         CASE('aea','aeb','cla','clb','prc','ice')
            CALL varArray%getData(1,var1d,index=n)
            CALL write_nc(ncid0,vname,var1d%d(:),ibeg0d)
            
         CASE('ztt','zmt')
            CALL varArray%getData(1,var1d,index=n)
            IF (ASSOCIATED(var1d%onDemand)) THEN
               CALL var1d%onDemand(vname,out1d)
            ELSE
               out1d = var1d%d
            END IF
            IF (myid == mpiroot) &
                 CALL write_nc(ncid0,vname,out1d(:),ibeg1d,icnt=icnt1d)

         CASE('xtytt')
            CALL varArray%getData(1,var2d,index=n)
            IF (ASSOCIATED(var2d%onDemand)) THEN
               CALL var2d%onDemand(vname,out2d)
            ELSE
               out2d = var2d%d
            END IF
            CALL write_nc(ncid0,vname,out2d(i1:i2,j1:j2),ibeg2d,icnt=icnt2d)
            
         CASE('zttaea','zttaeb','zttcla','zttclb','zttprc','zttice')
            CALL getSDdim(varArray%list(n)%dimension,nstr,nend)
            icnt1dsd = [nzp,nend-nstr+1,1]
            ALLOCATE(out1dsd(nzp,nend-nstr+1))
            CALL varArray%getData(1,var2d,index=n)
            IF (ASSOCIATED(var2d%onDemand)) THEN
               CALL var2d%onDemand(vname,out1dsd)
            ELSE
               out1dsd = var2d%d
            END IF
            IF (myid == mpiroot) &
                 CALL write_nc(ncid0,vname,out1dsd(:,nstr:nend),ibeg2d,icnt=icnt1dsd)
            DEALLOCATE(out1dsd)
            
         CASE('taea','taeb','tcla','tclb','tprc','tice')
            CALL getSDdim(varArray%list(n)%dimension,nstr,nend)
            icnt0dsd = [nend-nstr+1,1]
            ALLOCATE(outsd(nend-nstr+1)); outsd = 0.
            CALL varArray%getData(1,var1d,index=n)
            IF (ASSOCIATED(var1d%onDemand)) THEN
               CALL var1d%onDemand(vname,outsd)
            ELSE
               outsd = var1d%d
            END IF
            IF (myid == mpiroot) &
                 CALL write_nc(ncid0,vname,outsd(nstr:nend),ibeg1d,icnt=icnt0dsd)
            DEALLOCATE(outsd)
               
         CASE('tttt','mttt','tmtt','ttmt')
            CALL varArray%getData(1,var3d,index=n)
            IF (ASSOCIATED(var3d%onDemand)) THEN
               CALL var3d%onDemand(vname,out3d)
            ELSE
               out3d = var3d%d
            END IF            
            CALL write_nc(ncid0,vname,out3d(:,i1:i2,j1:j2),ibeg3d,icnt=icnt3d)
            
         CASE('ttttaea','ttttaeb','ttttcla','ttttclb','ttttprc','ttttice')
            CALL getSDdim(varArray%list(n)%dimension,nstr,nend)
            icnt3dsd = [nzp,nxp-4,nyp-4,nend-nstr+1,1]
            ALLOCATE(out3dsd(nzp,nxp,nyp,nend-nstr+1)); out3dsd = 0.
            CALL varArray%getData(1,var4d,index=n)
            IF (ASSOCIATED(var4d%onDemand)) THEN
               CALL var4d%onDemand(vname,out3dsd,nstr,nend)
            ELSE
               out3dsd = var4d%d
            END IF
            CALL write_nc(ncid0,vname,out3dsd(:,i1:i2,j1:j2,1:nend-nstr+1),ibeg4d,icnt=icnt3dsd)
            DEALLOCATE(out3dsd)
            
         END SELECT                  
      END DO
      
      var0d => NULL(); var1d => NULL(); var2d => NULL(); var3d => NULL(); var4d => NULL()
            
    END SUBROUTINE write_output

    SUBROUTINE getSDdim(dim,nstr,nend)
      CHARACTER(len=*), INTENT(in) :: dim
      INTEGER, INTENT(out) :: nstr,nend

      nstr = 0
      nend = 0
      
      SELECT CASE(dim)
      CASE('ttttaea','zttaea','taea')
         nstr = in1a; nend = fn2a
      CASE('ttttaeb','zttaeb','taeb')
         nstr = in2b; nend = fn2b
      CASE('ttttcla','zttcla','tcla')
         nstr = ica%cur; nend = fca%cur
      CASE('ttttclb','zttclb','tclb')
         nstr = icb%cur; nend = fcb%cur
      CASE('ttttprc','zttprc','tprc')
         nstr = 1; nend = nprc
      CASE('ttttice','zttice','tice')
         nstr = 1; nend = nice
      END SELECT
              
    END SUBROUTINE getSDdim
    
    
END MODULE mo_output
