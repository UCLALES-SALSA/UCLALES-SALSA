MODULE mo_mpi_io
  USE mpi_interface, ONLY : REAL_SIZE, INT_SIZE, MY_REAL, myid, pecount, mpiroot
  USE mpi
  IMPLICIT NONE
  
  INTERFACE write_hist_mpi
     MODULE PROCEDURE :: write_hist_parameters_int,     &
                         write_hist_parameters_int_1d,  &
                         write_hist_parameters_real,    &
                         write_hist_field_1d,           &
                         write_hist_field_2d,           &
                         write_hist_field_3d
  END INTERFACE write_hist_mpi

  INTERFACE read_hist_mpi
     MODULE PROCEDURE :: read_hist_parameters_int,    &
                         read_hist_parameters_int_1d, &
                         read_hist_parameters_real,   &
                         read_hist_field_1d,          &
                         read_hist_field_2d,          &
                         read_hist_field_3d
  END INTERFACE read_hist_mpi

  
  TYPE mpi_file_parameters
     INTEGER :: id
     INTEGER(kind=MPI_OFFSET_KIND) :: disp  ! Make sure the displacement gets updated uniquely for each mpi task
  END TYPE mpi_file_parameters

  CONTAINS 

    !-------------------------------------------------------------------------
    ! SUBROUTINE CREATE_MPI_HIST: Create a new file for writing restart files
    !
    SUBROUTINE create_mpi_hist(filename,fhist)
      CHARACTER(len=*), INTENT(in) :: filename
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr
      
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename,           &
                         MPI_MODE_CREATE + MPI_MODE_WRONLY,  &
                         MPI_INFO_NULL, fhist%id, ierr)

      ! Initialize displacement
      fhist%disp = 0
      
    END SUBROUTINE create_mpi_hist

    !------------------------------------------------------------------
    ! SUBROUTINE OPEN_MPI_HIST: Open existin restart file for reading
    !
    SUBROUTINE open_mpi_hist(filename,fhist)
      CHARACTER(len=*), INTENT(in) :: filename
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,    &
                         MPI_INFO_NULL, fhist%id, ierr)

      ! Initialize displacement
      fhist%disp = 0
      
    END SUBROUTINE open_mpi_hist
    
    !---------------------------------------------------
    ! SUBROUTINE CLOSE_MPI_HIST: Close the restart file
    !
    SUBROUTINE close_mpi_hist(fhist)
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr
      CALL MPI_FILE_CLOSE(fhist%id,ierr)      
    END SUBROUTINE close_mpi_hist
    
    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_INT: Write an individual parameter value, integer implementation
    !    
    SUBROUTINE write_hist_parameters_int(var,onlyroot,fhist)
      INTEGER, INTENT(in) :: var
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * INT_SIZE 
      
      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, 1,    &
                                  MPI_INTEGER, MPI_STATUS_IGNORE, ierr  )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, 1,   &
                                      MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
      
      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + INT_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount*INT_SIZE
      END IF
         
    END SUBROUTINE write_hist_parameters_int

    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_INT_1D: Write an individual parameter value, integer implementation
    !    
    SUBROUTINE write_hist_parameters_int_1d(nn,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: nn
      INTEGER, INTENT(in) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = nn
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * INT_SIZE
      
      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, nn,    &
                                  MPI_INTEGER, MPI_STATUS_IGNORE, ierr  )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, nn,   &
                                      MPI_INTEGER, MPI_STATUS_IGNORE, ierr)

      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * INT_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * INT_SIZE
      END IF
         
    END SUBROUTINE write_hist_parameters_int_1d
    
    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_REAL: Write an individual parameter value, float implementation
    !       
    SUBROUTINE write_hist_parameters_real(var,onlyroot,fhist)
      REAL, INTENT(in) :: var
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * REAL_SIZE  

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, 1,    &
                                  MY_REAL, MPI_STATUS_IGNORE, ierr )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, 1,   &
                                      MY_REAL, MPI_STATUS_IGNORE, ierr)

      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * REAL_SIZE
      END IF
         
    END SUBROUTINE write_hist_parameters_real    

    !------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_1D: Write 1d vectors to restart file 
    !    
    SUBROUTINE write_hist_field_1d(nn,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: nn
      REAL, INTENT(in) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp
      
      lsize = nn
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, lsize,    &
                                  MY_REAL, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, lsize,      &
                                      MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE
      END IF

    END SUBROUTINE write_hist_field_1d

    !-----------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_2D: Write 2d fields to restart file
    !    
    SUBROUTINE write_hist_field_2d(n2,n3,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(in) :: var(n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = n2*n3
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE  

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, lsize,    &
                                  MY_REAL, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, lsize,     &
                                      MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE
      END IF
         
    END SUBROUTINE write_hist_field_2d

    !------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_3D: Write 3d fields to restart file
    !
    SUBROUTINE write_hist_field_3d(n1,n2,n3,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in) :: var(n1,n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = n1*n2*n3
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(fhist%id, disp, var, lsize,    &
                                  MY_REAL, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &      
           CALL MPI_FILE_WRITE_AT_ALL(fhist%id, disp, var, lsize,     &
                                      MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE     
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE
      END IF
         
    END SUBROUTINE write_hist_field_3d
    
    ! -----------------------------------------
    ! IF ONLYROOT == TRUE, DISPLACEMENTS FOR ALL PROCESSORS ARE IDENTICAL
    SUBROUTINE read_hist_parameters_int(var,onlyroot,fhist)
      INTEGER, INTENT(out) :: var
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * INT_SIZE

      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, 1,    &
                                MPI_INTEGER, MPI_STATUS_IGNORE, ierr  )

      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + INT_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * INT_SIZE
      END IF

    END SUBROUTINE read_hist_parameters_int

    ! -----------------------------------------
    ! IF ONLYROOT == TRUE, DISPLACEMENTS FOR ALL PROCESSORS ARE IDENTICAL
    SUBROUTINE read_hist_parameters_int_1d(nn,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: nn
      INTEGER, INTENT(out) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = nn
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * INT_SIZE

      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, nn,    &
                                MPI_INTEGER, MPI_STATUS_IGNORE, ierr  )
      
      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * INT_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * INT_SIZE
      END IF

    END SUBROUTINE read_hist_parameters_int_1d

    
    ! ----------------------------------------------
    
    SUBROUTINE read_hist_parameters_real(var,onlyroot,fhist)
      REAL, INTENT(out) :: var
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * REAL_SIZE

      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, 1,    &
                                MY_REAL, MPI_STATUS_IGNORE, ierr  )

      ! Update displacement
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * REAL_SIZE
      END IF
    END SUBROUTINE read_hist_parameters_real

    ! -----------------------------------------------
    
    SUBROUTINE read_hist_field_1d(nn,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: nn
      REAL, INTENT(out) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = nn
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE
      
      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, lsize,      &
                                MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE
      END IF
    END SUBROUTINE read_hist_field_1d

    ! -------------------------------------------------
    
    SUBROUTINE read_hist_field_2d(n2,n3,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: n2,n3
      REAL, INTENT(out) :: var(n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = n2*n3
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE

      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, lsize,      &
                                MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE 
      END IF
    END SUBROUTINE read_hist_field_2d

    ! -------------------------------------------------
    
    SUBROUTINE read_hist_field_3d(n1,n2,n3,var,onlyroot,fhist)
      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(out) :: var(n1,n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      TYPE(mpi_file_parameters), INTENT(inout) :: fhist
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: lsize
      INTEGER(KIND=MPI_OFFSET_KIND) :: disp

      lsize = n1*n2*n3
      
      disp = fhist%disp
      IF (.NOT. onlyroot) &
           disp = fhist%disp + myid * lsize * REAL_SIZE
      
      CALL MPI_FILE_READ_AT_ALL(fhist%id, disp, var, lsize,      &
                                MY_REAL, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         fhist%disp = fhist%disp + lsize * REAL_SIZE
      ELSE
         fhist%disp = fhist%disp + pecount * lsize * REAL_SIZE 
      END IF
    END SUBROUTINE read_hist_field_3d

    ! --------------------------------------------------
        
    
END MODULE mo_mpi_io
