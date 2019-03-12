MODULE mo_mpi_io
  USE mpi_interface, ONLY : REAL_SIZE, INT_SIZE, myid, pecount, mpiroot
  USE mpi
  IMPLICIT NONE
  
  INTERFACE write_hist
     MODULE PROCEDURE :: write_hist_parameters_int,     &
                         write_hist_parameters_int_1d,  &
                         write_hist_parameters_real,    &
                         write_hist_parameters_real_1d, &  
                         write_hist_field_1d,           &
                         write_hist_field_2d,           &
                         write_hist_field_3d
  END INTERFACE write_hist

  INTERFACE read_hist
     MODULE PROCEDURE :: read_hist_parameters_int,    &
                         read_hist_parameters_real,   &
                         read_hist_field_1d,          &
                         read_hist_field_2d,          &
                         read_hist_field_3d
  END INTERFACE read_hist

  
  TYPE mpi_file_parameters
     INTEGER :: id
     INTEGER(kind=MPI_OFFSET_KIND) :: disp  ! Make sure the displacement gets updated uniquely for each mpi task
  END TYPE mpi_file_parameters
  
  TYPE(mpi_file_parameters), SAVE :: f_hist

  CONTAINS 

    !-------------------------------------------------------------------------
    ! SUBROUTINE CREATE_MPI_HIST: Create a new file for writing restart files
    !
    SUBROUTINE create_mpi_hist(filename)
      CHARACTER(len=*), INTENT(in) :: filename
      INTEGER :: ierr
      
      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename,           &
                         MPI_MODE_CREATE + MPI_MODE_WRONLY,  &
                         MPI_INFO_NULL, f_hist%id, ierr)

      ! Initialize displacement
      f_hist%disp = 0
      
    END SUBROUTINE create_mpi_hist

    !------------------------------------------------------------------
    ! SUBROUTINE OPEN_MPI_HIST: Open existin restart file for reading
    !
    SUBROUTINE open_mpi_hist(filename)
      CHARACTER(len=*), INTENT(in) :: filename
      INTEGER :: ierr

      CALL MPI_FILE_OPEN(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,    &
                         MPI_INFO_NULL, f_hist%id, ierr)

      ! Initialize displacement
      f_hist%disp = 0
      
    END SUBROUTINE open_mpi_hist
    
    !---------------------------------------------------
    ! SUBROUTINE CLOSE_MPI_HIST: Close the restart file
    !
    SUBROUTINE close_mpi_hist()
      INTEGER :: ierr
      CALL MPI_FILE_CLOSE(f_hist%id,ierr)      
    END SUBROUTINE close_mpi_hist
    
    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_INT: Write an individual parameter value, integer implementation
    !    
    SUBROUTINE write_hist_parameters_int(var,onlyroot)
      INTEGER, INTENT(in) :: var
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,1,INT_SIZE)
      
      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, 1,    &
                                  INT_SIZE, MPI_STATUS_IGNORE, ierr  )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, 1,   &
                                      INT_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! Update displacement
      IF (onlyroot) THEN
         f_hist%disp = f_hist%disp + INT_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,1,INT_SIZE)
      END IF
         
    END SUBROUTINE write_hist_parameters_int

    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_INT_1D: Write an individual parameter value, integer implementation
    !    
    SUBROUTINE write_hist_parameters_int_1d(nn,var,onlyroot)
      INTEGER, INTENT(in) :: nn
      INTEGER, INTENT(in) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,nn,INT_SIZE)
      
      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, nn,    &
                                  INT_SIZE, MPI_STATUS_IGNORE, ierr  )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, nn,   &
                                      INT_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! Update displacement
      IF (onlyroot) THEN
         f_hist%disp = f_hist%disp + nn*INT_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,nn,INT_SIZE)
      END IF
         
    END SUBROUTINE write_hist_parameters_int_1d
    
    !------------------------------------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_PARAMETERS_REAL: Write an individual parameter value, float implementation
    !       
    SUBROUTINE write_hist_parameters_real(var,onlyroot)
      REAL, INTENT(in) :: var
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,1,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, 1,    &
                                  REAL_SIZE, MPI_STATUS_IGNORE, ierr )

      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, 1,   &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)

      ! Update displacement
      IF (onlyroot) THEN
         f_hist%disp = f_hist%disp + REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,1,REAL_SIZE)
      END IF
         
    END SUBROUTINE write_hist_parameters_real    

    !------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_1D: Write 1d vectors to restart file 
    !    
    SUBROUTINE write_hist_field_1d(nn,var,onlyroot)
      INTEGER, INTENT(in) :: nn
      REAL, INTENT(in) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = nn
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                  REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, lsize,      &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
            
    END SUBROUTINE write_hist_field_1d

    !-----------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_2D: Write 2d fields to restart file
    !    
    SUBROUTINE write_hist_field_2d(n2,n3,var,onlyroot)
      INTEGER, INTENT(in) :: n2,n2
      REAL, INTENT(in) :: var(n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = n1*n2
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                  REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, lsize,     &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
         
    END SUBROUTINE write_hist_field_2d

    !------------------------------------------------------------------
    ! SUBROUTINE WRITE_HIST_FIELD_3D: Write 3d fields to restart file
    !
    SUBROUTINE write_hist_field_3d(n1,n2,n3,var,onlyroot)
      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(in) :: var(n1,n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = n1*n2*n3
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_WRITE_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                  REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &      
           CALL MPI_FILE_WRITE_AT_ALL(f_hist%id, f_hist%disp, var, lsize,     &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE     
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
         
    END SUBROUTINE write_hist_field_3d
    
    ! -----------------------------------------

    SUBROUTINE read_hist_parameters_int(var,onlyroot)
      INTEGER, INTENT(out) :: var
      LOGICAL, INTENT(in) :: onlyroot
      INTEGER :: ierr

      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,1,INT_SIZE)

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_READ_AT(f_hist%id, f_hist%disp, var, 1,    &
                                 INT_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_READ_AT_ALL(f_hist%id, f_hist%disp, var, 1,    &
                                     INT_SIZE, MPI_STATUS_IGNORE, ierr  )

      ! Update displacement
      IF (onlyroot) THEN
         f_hist%disp = f_hist%disp + INT_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,1,INT_SIZE)
      END IF

    END SUBROUTINE read_hist_parameters_int

    ! ----------------------------------------------
    
    SUBROUTINE read_hist_parameters_real(var,onlyroot)
      REAL, INTENT(out) :: var
      LOGICAL, INTENT(in) :: onlyroot
      INTEGER :: ierr
      LOGICAL :: onlyroot
      INTEGER :: ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,1,REAL_SIZE)

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_READ_AT(f_hist%id, f_hist%disp, var, 1,    &
                                 REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_READ_AT_ALL(f_hist%id, f_hist%disp, var, 1,    &
                                     REAL_SIZE, MPI_STATUS_IGNORE, ierr  )

      ! Update displacement
      IF (onlyroot) THEN
         f_hist%disp = f_hist%disp + REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,1,REAL_SIZE)
      END IF
    END SUBROUTINE read_hist_parameters_real

    ! -----------------------------------------------
    
    SUBROUTINE read_hist_field_1d(nn,var,onlyroot)
      INTEGER, INTENT(in) :: nn
      REAL, INTENT(out) :: var(nn)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = nn
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_READ_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                 REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_READ_AT_ALL(f_hist%id, f_hist%disp, var, lsize,      &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
    END SUBROUTINE read_hist_field_1d

    ! -------------------------------------------------
    
    SUBROUTINE read_hist_field_2d(n1,n2,var,onlyroot)
      INTEGER, INTENT(in) :: n1,n2
      REAL, INTENT(out) :: var(n1,n2)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = n1*n2
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_READ_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                 REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_READ_AT_ALL(f_hist%id, f_hist%disp, var, lsize,      &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
    END SUBROUTINE read_hist_field_2d

    ! -------------------------------------------------
    
    SUBROUTINE read_hist_field_3d(n1,n2,n3,var,onlyroot)
      INTEGER, INTENT(in) :: n1,n2,n3
      REAL, INTENT(out) :: var(n1,n2,n3)
      LOGICAL, INTENT(in) :: onlyroot
      LOGICAL :: l_root
      INTEGER :: lsize, ierr
      INTEGER(KIND=MPI_OFFSET_KIND) :: init_disp

      lsize = n1*n2*n3
      
      IF (onlyroot) THEN
         l_root = (myid == mpiroot)
      ELSE
         l_root = .TRUE.
      END IF

      init_disp = f_hist%disp
      IF (.NOT. onlyroot) &
           f_hist%disp = getRankDisplacement(f_hist%disp,lsize,REAL_SIZE)      

      IF (l_root .AND. onlyroot) &
           CALL MPI_FILE_READ_AT(f_hist%id, f_hist%disp, var, lsize,    &
                                 REAL_SIZE, MPI_STATUS_IGNORE, ierr )
      
      IF (l_root .AND. .NOT. onlyroot) &
           CALL MPI_FILE_READ_AT_ALL(f_hist%id, f_hist%disp, var, lsize,      &
                                      REAL_SIZE, MPI_STATUS_IGNORE, ierr)
      
      ! update the displacement to a next common value (same for all processes!) depending on the total process count,
      ! and the field size
      IF (onlyroot) THEN
         f_hist%disp = init_disp + lsize * REAL_SIZE
      ELSE
         f_hist%disp = getGlobalDisplacement(init_disp,lsize,REAL_SIZE)
      END IF
    END SUBROUTINE read_hist_field_3d

    ! --------------------------------------------------
    
    FUNCTION getRankDisplacement(idisp,isize,ikind)
      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(in) :: idisp
      INTEGER, INTENT(in) :: isize,ikind
      INTEGER(KIND=MPI_OFFSET_KIND) :: getRankDisplacement
      getRankDisplacement = idisp + myid * isize * ikind
    END FUNCTION getRankDisplacement

    FUNCTION getGlobalDisplacement(idisp,isize,ikind)
      INTEGER(KIND=MPI_OFFSET_KIND), INTENT(in) :: idisp
      INTEGER, INTENT(in) :: isize,ikind
      INTEGER(KIND=MPI_OFFSET_KIND) :: getGlobalDisplacement
      getGlobalDisplacement = idisp + pecount * lsize * ikind
    END FUNCTION
    
    
END MODULE mo_mpi_io
