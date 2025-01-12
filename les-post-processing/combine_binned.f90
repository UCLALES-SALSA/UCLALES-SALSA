!************************************************************************
! Post-processing UCLALES-SALSA analysis and column outputs.
!    Tomi Raatikainen, March 2020, FMI
!
! Versions
!   20200316    The first version
!
! Compile and run
! ===============
! a) Cygwin
!  Compile
!    gfortran -O2 -I/usr/include -o ./pples combine.f90 -lnetcdff
!  Run
!    ./pples <file name prefix>
! b) Puhti
!  Compile
!    module load intel/19.0.4
!    module load netcdf-fortran/4.4.4
!    ifort -O2 -o ./pples combine.f90 -lnetcdff
!  Run
!    srun --ntasks=1 --time=00:0:10 --partition=<partition> --account=<project> pples <file name prefix> <key1=value1 key2=value2 ...>
!  Examples
!    srun --ntasks=1 --time=00:10:00 --mem=4000 --partition=fmi --account=project_2001823 pples dycoms_l4 deflate_level=4
!
!************************************************************************
PROGRAM combine
    USE netcdf
    IMPLICIT NONE

    ! Current version of this file
    CHARACTER(LEN=8), PARAMETER :: version='20200316'

    CHARACTER(len=100) :: fname, iname, oname, vname, aname, variable_name
    character (len=8)  :: date
    INTEGER :: ii, jj, kk, i, j, imax, jmax, n, nmax
    INTEGER :: ncid, ncid_src, iret, ityp, binned_id
    INTEGER :: ndims, nvars, natt
    INTEGER , DIMENSION(5):: iarray, dims, fst, lst ! Max 5 dimensions
    REAL, ALLOCATABLE :: vec_new(:), vec_old(:), mval(:,:,:,:,:)
    INTEGER, ALLOCATABLE :: indices(:,:,:), ind_old(:), dimsize_new(:), dimsize_old(:)
    LOGICAL :: file_exists
    ! Default settings
    !    Compression
    !       a) If shuffle is non-zero, turn on the shuffle filter
    !       b) If deflate is non-zero, turn on the deflate filter at the level specified by the deflate_level (0-9)
    INTEGER :: shuffle=0, deflate=0, deflate_level=0
    !
    ! Number of input arguments
    n = command_argument_count()
    IF (n==0) STOP 'At least file name prefix is needed as input!'
    !
    ! File name as the first  input argument
    CALL get_command_argument(1, fname)
    !
    ! Other optional arguments (key=value)
    DO i=2,n
        CALL get_command_argument(i,iname)
        j=INDEX(iname,'=')
        IF (j==0) THEN
            WRITE(*,*) 'Bad input argument: '//TRIM(iname)
            STOP
        ENDIF
        SELECT CASE (iname(1:j-1))
            CASE ('deflate_level')
                READ(iname(j+1:),*) deflate_level
                IF (deflate_level>0) deflate=1 ! Compression is on
            CASE ('shuffle')
                READ(iname(j+1:),*) shuffle
            CASE ('variable_name')
                READ(iname(j+1:),*) variable_name
                WRITE(*,*) 'Variable requested is ' 
                WRITE(*,*) variable_name
            CASE DEFAULT
                WRITE(*,*)'Bad argument',i,': '//TRIM(iname)
                STOP
        END SELECT
    ENDDO
    !
    ! ************************ Check the input data ************************
    !
    WRITE(*,*) "Checking files '"//TRIM(fname)//"'..."
    ! Count the number of files (0:imax and 0:jmax)
    ! The first file
    i=0;  j=0
    WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
    INQUIRE(FILE=iname, EXIST=file_exists)
    IF (.NOT. file_exists) THEN
        WRITE(*,*) 'Examining file '//iname
        STOP 'Data not found!'
    ENDIF
    ! Dimension i
    file_exists= .TRUE.
    Do WHILE (file_exists)
        i=i+1
        WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
        INQUIRE(FILE=iname, EXIST=file_exists)
    ENDDO
    imax=i-1
    ! Dimension j
    file_exists= .TRUE.
    Do WHILE (file_exists)
        i=0
        j=j+1
        WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
        INQUIRE(FILE=iname, EXIST=file_exists)
        IF (file_exists) THEN
            DO i=0,imax
                WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
                INQUIRE(FILE=iname, EXIST=file_exists)
                IF (.NOT. file_exists) THEn
                    WRITE(*,*) 'Examining file '//iname
                    STOP 'Data not found!'
                ENDIF
            ENDDO
        ENDIF
    ENDDO
    jmax=j-1
    IF (imax==0 .AND. jmax==0) THEN
        STOP 'Only one file found?'
    ENDIF
    WRITE(*,'(A8,I3,A6)') '  Found ',(imax+1)*(jmax+1),' files'
    WRITE(*,*) 'Done'
    !
    ! ************************ Define output file ************************
    !
    WRITE(*,*) ' '
    WRITE(*,*) 'Defining output netCDF...'
    ! Open the first input file for reading some basic information
    i=0; j=0
    WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
    iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
    IF (iret/=nf90_noerr) THEN
        WRITE(*,*) 'Error in opening the first input file ('//TRIM(iname)//')!'
        WRITE(*,*) nf90_strerror(iret)
        STOP
    ENDIF
    iret = nf90_inquire(ncid_src,nDimensions=ndims, nVariables=nvars,nAttributes=natt,formatNum=ii)
    IF (iret/=nf90_noerr) THEN
        WRITE(*,*) 'Error in reading the first input file ('//TRIM(iname)//')!'
        WRITE(*,*) nf90_strerror(iret)
        STOP
    ENDIF
    ! File size
    INQUIRE(FILE=iname, SIZE=jj)
    !
    ! Generate output
    WRITE(oname,'(a,a3)') trim(fname),'.nc'
    IF (shuffle>0 .OR. deflate>0) THEN
        ! Compression requires netCDF4
        ii=nf90_netcdf4
    ELSEIF ((imax+1)*(jmax+1)*jj>2e9) THEN
        WRITE(*,"('   Warning: file size exceeds 2 GB (total approx. ',F3.1, ' GB) - changing to netCDF4')") &
                REAL((imax+1)*(jmax+1)*jj/1073741824)
        ! Format supporting large data files : nf90_64bit_offset and nf90_netcdf4
        !ii=nf90_64bit_offset ! Not readable by Igor?
        ii=nf90_netcdf4 ! This works
    ENDIF
    oname = TRIM(fname) // TRIM('_') 
    oname = TRIM(oname) // TRIM(variable_name)
    oname = TRIM(oname) // TRIM('.nc') 
    WRITE(*,*)'  Creating file '//TRIM(oname)//' in mode',ii
    iret = nf90_create(oname,ii,ncid)
    IF (iret/=nf90_noerr) THEN
        WRITE(*,*) 'Error in creating the output file ('//TRIM(oname)//')!'
        WRITE(*,*) nf90_strerror(iret)
        STOP
    ENDIF
    ! New global attributes about post-processing
    call date_and_time(date)
    iret = nf90_put_att(ncid,NF90_GLOBAL,'PP_date','Post-processing date '//date)
    iret = nf90_put_att(ncid,NF90_GLOBAL,'PP_version','Post-processing code '//version)
    !
    ! Copy global attributes
    DO i=1,natt
        ! Name
        iret = nf90_inq_attname(ncid_src,NF90_GLOBAL,i,vname)
        ! Copy
        iret = nf90_copy_att(ncid_src,NF90_GLOBAL,vname,ncid,NF90_GLOBAL)
    ENDDO

    ! Read information about the requested variable
    iret = nf90_inq_varid(ncid_src, variable_name, binned_id)
    WRITE(*,*) binned_id
    iret = nf90_close(ncid_src)
    !
    !
    ! a) Dimensions
    ! Generate a map based on dimensions
    ALLOCATE( indices(1:ndims,0:imax,0:jmax), ind_old(ndims), dimsize_new(ndims), dimsize_old(ndims) )
    DO ii=1,ndims
        ! All files
        DO i=0,imax
            DO j=0,jmax
                ! Source
                WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
                iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
                IF (iret/=NF90_NOERR) THEN
                    WRITE(*,*) 'Examining file '//trim(iname)
                    WRITE(*,*) nf90_strerror(iret)
                    STOP 'Error in opening data file!'
                ENDIF
                ! Source data
                IF (j==0 .AND. i==0) THEN
                    ! Length of dimension
                    iret = nf90_inquire_dimension(ncid_src,ii,name=vname,len=n)
                    !
                    ! Variable index (ii:th dimension may not be the ii:th variable)
                    iret = nf90_inq_varid(ncid_src,vname,jj)
                    ind_old(ii) = jj
                    !
                    ! Variable type
                    iret = nf90_inquire_variable(ncid_src,jj,xtype=ityp)
                    !
                    ! Allocate arrays for current and new dimension (final size not known)
                    ALLOCATE( vec_old(n), vec_new((imax+1)*(jmax+1)*n) )
                    vec_new(:)=-999.
                    nmax=0
                ENDIF
                iret = nf90_get_var(ncid_src,jj,vec_old)
                !
                IF (nmax==0) THEN
                    ! Initialize output
                    vec_new(1:n)=vec_old(1:n)
                    indices(ii,i,j)=1
                    nmax=n
                ELSEIF (MINVAL(vec_old)>MAXVAL(vec_new(1:nmax))) THEN
                    ! Append after previous data
                    vec_new(nmax+1:nmax+n)=vec_old(1:n)
                    indices(ii,i,j)=nmax+1
                    nmax=nmax+n
                ELSEIF (ALL(ABS(vec_old(1:n)-vec_new(1:n))<1e-10)) THEN
                    ! Identical with the first n values, i.e. constant dimension
                    indices(ii,i,j)=1
                ELSEIF (MAXVAL(vec_old)<MINVAL(vec_new(1:nmax))) THEN
                    ! Append before => error!
                    WRITE(*,*) 'Examining file '//trim(iname)//' and variable '//trim(vname)
                    STOP 'Monotonic order expected (values before current range)!'
                ELSE
                    ! There should be a matching sequence
                    indices(ii,i,j)=find_sequence(vec_old,n,vec_new(1:nmax),nmax)
                    IF (indices(ii,i,j)<0) THEN
                        WRITE(*,*) 'Examining file '//trim(iname)//' and variable '//trim(vname)
                        stop 'Matching sequence of values not found!'
                    ENDIF
                ENDIF
                iret = nf90_close(ncid_src)
            ENDDO
        ENDDO
        dimsize_new(ii) = nmax ! Output dimensions
        dimsize_old(ii) = n ! Input dimensions
        !
        ! Create dimension (ii=jj=kk)
        iret = nf90_def_dim(ncid,TRIM(vname),len=nmax,dimid=jj)
        iret = nf90_def_var(ncid,TRIM(vname),xtype=ityp,dimids=jj,varid=kk)
        !
        ! Copy attributes
        i=0; j=0
        WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
        iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
        iret = nf90_inq_varid(ncid_src,TRIM(vname),jj)
        iret = nf90_inquire_variable(ncid_src,jj,nAtts=natt)
        DO i=1,natt
            ! Name
            iret = nf90_inq_attname(ncid_src,jj,i,aname)
            ! Copy
            iret = nf90_copy_att(ncid_src,jj,TRIM(aname),ncid,kk)
        ENDDO
        iret = nf90_close(ncid_src)
        !
        ! Clean
        DEALLOCATE(vec_old,vec_new)
    ENDDO
    !
    !
    ! b) Other variables
    i=0; j=0
    WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
    iret = nf90_open(iname,NF90_NOWRITE,ncid_src)

    DO ii=1,nvars
        ! Current dimensions for variable ii
        iret = nf90_inquire_variable(ncid_src,ii,name=vname,xtype=ityp,ndims=n,dimids=iarray)
        IF (n>5) STOP 'Max 5 dimensions!'
        IF (n<5) CYCLE
        !
        ! Skip dimension variables
        iret = nf90_inq_dimid(ncid_src,vname,jj)
        IF (iret == nf90_noerr) CYCLE
        IF (ii /= binned_id) CYCLE
        !
        ! Create variable (the same dimension variables as before)
        iret = nf90_def_var(ncid,TRIM(vname),xtype=ityp,dimids=iarray(1:n),varid=jj)
        ! Copy attributes
        iret = nf90_inquire_variable(ncid_src,ii,nAtts=natt)
        DO i=1,natt
            ! Name
            iret = nf90_inq_attname(ncid_src,ii,i,aname)
            ! Copy
            iret = nf90_copy_att(ncid_src,ii,TRIM(aname),ncid,jj)
        ENDDO
        !
        ! Compression (optional)
        IF (shuffle>0 .OR. deflate>0) THEN
            iret = nf90_def_var_deflate(ncid,jj,shuffle,deflate,deflate_level)
            IF (iret/=nf90_noerr) THEN
                WRITE(*,*) 'Error when setting compression for variable '//TRIM(vname)//'!'
                WRITE(*,*) nf90_strerror(iret)
                STOP
            ENDIF
        ENDIF
    ENDDO
    iret = nf90_close(ncid_src)
    WRITE(*,*) 'Done'
    !
    !
    ! ************************ Generate output data ************************
    iret = nf90_enddef(ncid)
    !
    !
    ! a) Dimensions
    WRITE(*,*) ' '
    WRITE(*,*) 'Creating dimensions...'
    DO ii=1,ndims
        ! Allocate arrays for current and new dimension
        ALLOCATE( vec_old(dimsize_old(ii)), vec_new(dimsize_new(ii)) )
        n=dimsize_old(ii)
        nmax=dimsize_new(ii)
        !
        ! Variable ID
        jj = ind_old(ii)
        !
        ! All files
        DO i=0,imax
            DO j=0,jmax
                ! Source
                WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
                iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
                !
                ! Source data
                iret = nf90_get_var(ncid_src,jj, vec_new(indices(ii,i,j):indices(ii,i,j)+n-1) )
                !
                ! Close input
                iret = nf90_close(ncid_src)
            ENDDO
        ENDDO
        !
        ! Save updated dimensions
        iret = nf90_put_var(ncid,ii,vec_new(1:nmax))
        !
        ! Clean
        DEALLOCATE(vec_old,vec_new)
        !
        ! Print info
        iret = nf90_inquire_dimension(ncid,ii,name=vname)
        WRITE(*,"('  ',I3,' ',A8,I4,' => ',I4)") ii,vname,n,nmax
    ENDDO
    iret = nf90_sync(ncid)
    WRITE(*,*) 'Done'
    !
    !
    ! b) Other variables
    WRITE(*,*) ' '
    WRITE(*,*) 'Creating variables...'
    DO ii=1,nvars
        ! The first source file
        i=0; j=0
        WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
        iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
        !
        ! Information about variable ii
        iret = nf90_inquire_variable(ncid_src,ii,name=vname,ndims=n,dimids=iarray)
        !
        ! Skip dimension variables
        iret = nf90_inq_dimid(ncid_src,TRIM(variable_name),jj)
        IF (iret == nf90_noerr) CYCLE
        IF (ii /= binned_id) CYCLE
        !
        ! Determine output dimensions
        dims(:) = 1
        dims(1:n) = dimsize_new(iarray(1:n))
        ! Allocate output data array
        ALLOCATE( mval(dims(1),dims(2),dims(3),dims(4),dims(5)) )
        mval = -999.
        !
        ! Close
        iret = nf90_close(ncid_src)
        !
        ! All files
        DO i=0,imax
            DO j=0,jmax
                ! Source
                WRITE(iname,'(a,a1,i4.4,i4.4,a3)') trim(fname),'.',i,j,'.nc'
                iret = nf90_open(iname,NF90_NOWRITE,ncid_src)
                !
                ! Indices for the current PU (first:last)
                fst(:) = 1
                fst(1:n) = indices(iarray(1:n),i,j)
                lst(:) = 1
                lst(1:n) = fst(1:n)+dimsize_old(iarray(1:n))-1
                !
                ! Read data directly to the output array
                iret = nf90_get_var(ncid_src,ii, mval(fst(1):lst(1),fst(2):lst(2),fst(3):lst(3),fst(4):lst(4),fst(5):lst(5)) )
                !
                ! Close input
                iret = nf90_close(ncid_src)
            ENDDO
        ENDDO
        ! Save the data
        iret = nf90_inq_varid(ncid,vname,jj)
        iret = nf90_put_var(ncid,jj,mval(:,:,:,:,:))
        iret = nf90_sync(ncid)
        !
        ! Clean
        DEALLOCATE(mval)
        !
        ! Print info
        WRITE(*,"('  ',I3,' ',A8,I3,'D:',5I3)") jj,vname,n,iarray(1:n)
    ENDDO
    !
	! Close file
    iret = nf90_close(ncid)
    WRITE(*,*) 'Done'

CONTAINS

    INTEGER FUNCTION find_sequence(seq,nseq,base,nbase)
        ! Find the starting index of vector seq in a longer vector base
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: nseq, nbase
        REAL, INTENT(IN) :: seq(nseq), base(nbase)
        INTEGER :: i
        !
        find_sequence = -1
        DO i=1,nbase-nseq+1
            IF ( ALL(abs(base(i:i+nseq-1)-seq(1:nseq))<1e-5 ) ) THEN
                find_sequence = i
                RETURN
            ENDIF
        ENDDO
        !
    END FUNCTION find_sequence

END PROGRAM combine
