MODULE mo_output
  USE netcdf
  USE mpi_interface, ONLY : myid, mpiroot, ver, author, info, wrxid, wryid, pecount
  USE ncio, ONLY : StreamDef, close_nc, sync_nc
  USE grid, ONLY : outAxes,outAxesPS,outAxesTS,expnme,nzp,nxp,nyp,filprf,  &
                   level
  USE mo_field_state, ONLY : outProg, outVector, outDiag, outDerived, outPS, outTS
  USE classFieldArray, ONLY : FieldArray
  USE mo_structured_datatypes
  USE mo_submctl, ONLY : in1a,fn2a,in2b,fn2b,ica,fca,icb,fcb,nprc,nice
  
  IMPLICIT NONE

  ! NOTE ABOUT TS AND PS FILES: These are only written by the rank mpiroot. Initialization and closing calls should be
  ! done only for myid = mpiroot. Writing procedures however must be called by all processes, because the onDemand functions
  ! for the statistics use MPI's collective reduction procedures. However, only mpiroot will actually write the output, i.e.
  ! if the variable dimensions imply statistical output, only mpiroot will call the ncio.f90/write_nc subroutine.

  
  LOGICAL :: tsflg = .FALSE., psflg = .FALSE. ! Flags to fetch data for statistical outputs from within physics routines; should coincide with the statistics output timestep
  REAL :: ps_intvl = 120.
  REAL :: ts_intvl = 120.
  REAL :: main_intvl = 3600.
  TYPE(StreamDef) :: StreamMain, StreamPS, StreamTS

  
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
      CHARACTER(len=150) :: fname
      
      npoints = (nxp-4)*(nyp-4)

      fname = trim(filprf)

      IF (pecount > 1 ) THEN
         WRITE(fname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',wrxid,wryid,'.nc'
      ELSE
         WRITE(fname,'(a,a3)') trim(fname),'.nc'
      END IF
      
      IF ( .NOT. ANY([outProg%Initialized,     &
                      outDiag%Initialized,     &
                      outDerived%Initialized,  &
                      outVector%Initialized]   &
                    ) ) RETURN ! If no output variables defined, do not open the files 
      
      IF(myid == mpiroot) WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname)
      
      CALL StreamMain%open_nc(fname,expnme,time,npoints,ver,author,info)

      IF (level < 4) THEN
         CALL StreamMain%define_nc_dims(n1=nzp,n2=nxp-4,n3=nyp-4   )
      ELSE IF (level == 4) THEN
         CALL StreamMain%define_nc_dims(n1=nzp,n2=nxp-4,n3=nyp-4,                  &
                                        inae_a=fn2a,inae_b=fn2b-fn2a,              &
                                        incld_a=fca%cur,incld_b=fcb%cur-fca%cur,   &
                                        inprc=nprc                                 )               
      ELSE IF (level == 5) THEN
         CALL StreamMain%define_nc_dims(n1=nzp,n2=nxp-4,n3=nyp-4,                  &
                                        inae_a=fn2a,inae_b=fn2b-fn2a,              &
                                        incld_a=fca%cur,incld_b=fcb%cur-fca%cur,   &                              
                                        inprc=nprc,inice=nice                      )              
      END IF

      CALL StreamMain%define_nc_vars([outProg,outDiag,outDerived,outVector,outAxes])
      
      IF (myid == mpiroot) WRITE(*,*) '   ...starting record: ', StreamMain%nrec
      
    END SUBROUTINE init_main

    ! ------------------------------------------------------------------------

    SUBROUTINE init_ps(time)
      REAL, INTENT(in) :: time
      INTEGER :: npoints
      CHARACTER(len=150) :: fname
      
      npoints = (nxp-4)*(nyp-4)
      
      fname = TRIM(filprf)//'.ps.nc'
      
      IF (.NOT. outPS%Initialized) RETURN ! If no variables defined for output, do not create the file 
      
      IF(myid == mpiroot) THEN
         
         WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname)
      
         CALL StreamPS%open_nc(fname,expnme,time,npoints,ver,author,info)
      
         IF (level < 4) THEN
            CALL StreamPS%define_nc_dims(n1=nzp)         
         ELSE IF (level == 4) THEN
            CALL StreamPS%define_nc_dims(n1=nzp,inae_a=fn2a,inae_b=fn2b-fn2a,     &
                                         incld_a=fca%cur,incld_b=fcb%cur-fca%cur, &
                                         inprc=nprc                               )         
         ELSE IF (level == 5) THEN
            CALL StreamPS%define_nc_dims(n1=nzp,inae_a=fn2a,inae_b=fn2b-fn2a,     &
                                         incld_a=fca%cur,incld_b=fcb%cur-fca%cur, &
                                         inprc=nprc,inice=nice                    )                 
         END IF

         CALL StreamPS%define_nc_vars([outPS,outAxesPS])
         
         WRITE(*,*) '   ...starting record: ', StreamPS%nrec

      END IF
         
    END SUBROUTINE init_ps

    ! ------------------------------------------------------------------------

    SUBROUTINE init_ts(time)
      REAL, INTENT(in) :: time
      INTEGER :: npoints
      CHARACTER(len=150) :: fname
      
      npoints = (nxp-4)*(nyp-4)
      
      fname = TRIM(filprf)//'.ts.nc'

      IF (.NOT. outTS%Initialized) RETURN ! If no variables defined for output, do not create the file 
      
      IF(myid == mpiroot) THEN
         
         WRITE(*,"(//' ',49('-')/,' ',/,'   Initializing: ',A20)") trim(fname)

         CALL StreamTS%open_nc(fname,expnme,time,npoints,ver,author,info)

         IF (level < 4) THEN
            CALL StreamTS%define_nc_dims()
         ELSE IF (level == 4) THEN
            CALL StreamTS%define_nc_dims(inae_a=fn2a,inae_b=fn2b-fn2a,            &
                                         incld_a=fca%cur,incld_b=fcb%cur-fca%cur, &
                                         inprc=nprc                               )            
         ELSE IF (level == 5) THEN
            CALL StreamTS%define_nc_dims(inae_a=fn2a,inae_b=fn2b-fn2a,            &
                                         incld_a=fca%cur,incld_b=fcb%cur-fca%cur, &
                                         inprc=nprc,inice=nice                    )            
         END IF

         CALL StreamTS%define_nc_vars([outTS,outAxesTS])
         
         write(*,*) '   ...starting record: ', StreamTS%nrec

      END IF
         
    END SUBROUTINE init_ts
    
    ! -------------------------------------------------------------------------
    
    SUBROUTINE close_main()
      CALL close_nc(StreamMain%ncid)
    END SUBROUTINE close_main

    SUBROUTINE close_ps()
      CALL close_nc(StreamPS%ncid)
    END SUBROUTINE close_ps

    SUBROUTINE close_ts()
      CALL close_nc(StreamTS%ncid)
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

      ibeg0 = [StreamMain%nrec]

      IF ( .NOT. ANY([outProg%Initialized,     &
                      outDiag%Initialized,     &
                      outDerived%Initialized]  &
                    ) ) RETURN ! If no output variables defined, do not try to write
      
      ! write time
      CALL StreamMain%write_nc('time',time,ibeg0)
      
      IF (StreamMain%nrec == 1) THEN
         ! First entry -> write axis variables
         CALL write_output(outAxes,StreamMain)
      END IF
      
      IF (outProg%Initialized) &
           CALL write_output(outProg,StreamMain)
      IF (outVector%Initialized) &
           CALL write_output(outVector,StreamMain)
      IF (outDiag%Initialized) &
           CALL write_output(outDiag,StreamMain)
      IF (outDerived%Initialized) &
           CALL write_output(outDerived,StreamMain)

      IF (myid == 0) WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")    &
         StreamMain%nrec,StreamMain%fname

      CALL sync_nc(StreamMain%ncid)
      StreamMain%nrec = StreamMain%nrec+1      
      
    END SUBROUTINE write_main

    ! --------------------------------------------------------------------------

    SUBROUTINE write_ps(time)   
      REAL, INTENT(in) :: time
      INTEGER :: ibeg0(1)

      ibeg0 = [StreamTS%nrec]

      IF ( .NOT. outPS%Initialized ) RETURN
      
      ! write time
      IF (myid == mpiroot) &         
           CALL StreamPS%write_nc('time',time,ibeg0)

      IF (StreamPS%nrec == 1) &
           ! First entry -> write axis variables
           CALL write_output(outAxesPS,StreamPS)

      IF (outPS%Initialized) &
           CALL write_output(outPS,StreamPS)

      IF (myid == mpiroot) &
           WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")     &
           StreamPS%nrec,StreamPS%fname      

      IF (myid == mpiroot) &
           CALL sync_nc(StreamPS%ncid)

      StreamPS%nrec = StreamPS%nrec+1
      
    END SUBROUTINE write_ps

    ! ----------------------------------------------------------------------------

    SUBROUTINE write_ts(time)
      REAL, INTENT(in) :: time
      INTEGER :: ibeg0(1)

      ibeg0 = [StreamTS%nrec]

      IF ( .NOT. outTS%Initialized ) RETURN

      ! Write time
      IF (myid == mpiroot) &
           CALL StreamTS%write_nc('time',time,ibeg0)

      IF (StreamTS%nrec == 1) &
           ! First entry -> write axis variables
           CALL write_output(outAxesTS,StreamTS)

      IF (outTS%Initialized)  &
           CALL write_output(outTS,StreamTS)

      IF (myid == mpiroot) &
           WRITE(*,"(//' ',12('-'),'   Record ',I3,' to: ',A60 //)")    &
           StreamTS%nrec,StreamTS%fname   

      IF (myid == mpiroot) &
           CALL sync_nc(StreamTS%ncid)
         
      StreamTS%nrec = StreamTS%nrec+1
     
    END SUBROUTINE write_ts    
    !
    ! Subroutine WRITE_OUTPUT: A general wrap-around routine used to to the writing of all kinds of
    !                          output variables
    !
    !
    SUBROUTINE write_output(varArray,stream)
      TYPE(FieldArray), INTENT(in) :: varArray
      TYPE(StreamDef), INTENT(in)  :: stream
      
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
      
      ! Count arrays corresponding to 1 record. For non-binned outputs these are always the same
      icnt1d = [nzp,1]
      icnt2d = [nxp-4,nyp-4,1]
      icnt3d = [nzp,nxp-4,nyp-4,1]      
      !icntXXsd for binned variables have to be allocated on a case by case basis in the loop below!

      ! Start indices for a record; i.e. ones except for the record dimension (time)
      ibeg1d = [1,stream%nrec]
      ibeg2d = [1,1,stream%nrec]
      ibeg3d = [1,1,1,stream%nrec]
      ibeg0d = [stream%nrec]
      ibeg4d = [1,1,1,1,stream%nrec]
      
      nvar = varArray%count
      DO n = 1,nvar
         vname = varArray%list(n)%name
         SELECT CASE(varArray%list(n)%dimension)

         ! Axis variables
         ! -------------------------------------------   
         CASE('time') ! Time is identical for all processes -> write only from root
            CALL varArray%getData(1,var0d,index=n)
            IF (ASSOCIATED(var0d%onDemand)) THEN  ! Check whether data is stored or obtained via on-demand function
               CALL var0d%onDemand(out0d)
            ELSE
               out0d = var0d%d
            END IF
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,out0d,ibeg0d)
            
         CASE('zt','zm') ! Vertical axis is indentical for all processes -> write only from root
            CALL varArray%getData(1,var1d,index=n)
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,var1d%d(:),ibeg0d)
            
         CASE('xt','xm') ! Lateral axes are divided between processes -> write from all with appr. indices
            CALL varArray%getData(1,var1d,index=n)
            CALL stream%write_nc(vname,var1d%d(i1:i2),ibeg0d)
            
         CASE('yt','ym') ! - '' -
            CALL varArray%getData(1,var1d,index=n)
            CALL stream%write_nc(vname,var1d%d(j1:j2),ibeg0d)

         CASE('aea','aeb','cla','clb','prc','ice') ! Bin axes identical -> only from root
            CALL varArray%getData(1,var1d,index=n)
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,var1d%d(:),ibeg0d)
            
         CASE('ztt','zmt') ! Time-height dimension is exclusively for ps-stream -> only from root
            CALL varArray%getData(1,var1d,index=n)
            IF (ASSOCIATED(var1d%onDemand)) THEN
               CALL var1d%onDemand(out1d)
            ELSE
               out1d = var1d%d
            END IF
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,out1d(:),ibeg1d,icnt=icnt1d)

            
         ! Physical variables
         ! -------------------------------------------------   
         CASE('xtytt')
            CALL varArray%getData(1,var2d,index=n)
            IF (ASSOCIATED(var2d%onDemand)) THEN
               CALL var2d%onDemand(vname,out2d)
            ELSE
               out2d = var2d%d
            END IF
            CALL stream%write_nc(vname,out2d(i1:i2,j1:j2),ibeg2d,icnt=icnt2d)
            
         CASE('zttaea','zttaeb','zttcla','zttclb','zttprc','zttice') ! Binned time-height for ps-stream -> from root
            CALL getSDdim(varArray%list(n)%dimension,nstr,nend) ! Check the start and end indices for the bins
            icnt1dsd = [nzp,nend-nstr+1,1]                      ! Set the count array based on the bin indices
            ALLOCATE(out1dsd(nzp,nend-nstr+1)); out1dsd = 0.    ! Allocate the output array based on the bin indices
            CALL varArray%getData(1,var2d,index=n)
            IF (ASSOCIATED(var2d%onDemand)) THEN                ! Check if data is stored or obtained via on-demand function
               CALL var2d%onDemand(vname,out1dsd,nstr,nend)
            ELSE
               out1dsd = var2d%d
            END IF
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,out1dsd(:,:),ibeg2d,icnt=icnt1dsd)
            DEALLOCATE(out1dsd)
            
         CASE('taea','taeb','tcla','tclb','tprc','tice') ! Binned temporal statistics for ts-stream -> from root
            CALL getSDdim(varArray%list(n)%dimension,nstr,nend)
            icnt0dsd = [nend-nstr+1,1]
            ALLOCATE(outsd(nend-nstr+1)); outsd = 0.
            CALL varArray%getData(1,var1d,index=n)
            IF (ASSOCIATED(var1d%onDemand)) THEN
               CALL var1d%onDemand(outsd)
            ELSE
               outsd = var1d%d
            END IF
            IF (myid == mpiroot) &
                 CALL stream%write_nc(vname,outsd(:),ibeg1d,icnt=icnt0dsd)
            DEALLOCATE(outsd)
               
         CASE('tttt','mttt','tmtt','ttmt')
            CALL varArray%getData(1,var3d,index=n)
            IF (ASSOCIATED(var3d%onDemand)) THEN
               CALL var3d%onDemand(vname,out3d)
            ELSE
               out3d = var3d%d
            END IF            
            CALL stream%write_nc(vname,out3d(:,i1:i2,j1:j2),ibeg3d,icnt=icnt3d)
            
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
            CALL stream%write_nc(vname,out3dsd(:,i1:i2,j1:j2,:),ibeg4d,icnt=icnt3dsd)
            DEALLOCATE(out3dsd)
            
         END SELECT                  
      END DO
      
      var0d => NULL(); var1d => NULL(); var2d => NULL(); var3d => NULL(); var4d => NULL()
            
    END SUBROUTINE write_output


    !
    !---------------------------------------------------------------------------
    ! Determines the start and end bin indices for size distribution variables
    !
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
