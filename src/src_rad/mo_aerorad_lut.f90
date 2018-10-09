MODULE mo_aerorad_lut

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: init_aerorad_lut

CONTAINS



  SUBROUTINE check(istatus)
    !!Checks the return status of a nf90 function. If an error
    !!occured, the according error message is written to the
    !!screen and the program exits

    USE netcdf
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: istatus
    IF (istatus /= nf90_noerr) THEN
       WRITE(*,*) TRIM(ADJUSTL(nf90_strerror(istatus)))
       WRITE(*,*) 'exiting'
       STOP 
    END IF
  END SUBROUTINE check


  SUBROUTINE init_aerorad_lut(filename,dalpha, dn_re, dn_im,         & 
                              n_re, n_im, alpha,            &
                              lut_sigma,lut_asym,lut_omega  )

    USE netcdf
    IMPLICIT NONE

    CHARACTER (len = *), INTENT(in) :: filename
    REAL, INTENT(out), ALLOCATABLE :: n_im(:),n_re(:),alpha(:),lut_sigma(:,:,:),lut_asym(:,:,:),lut_omega(:,:,:)
    INTEGER, INTENT(out) :: dalpha,dn_im,dn_re

    !dimension names
    CHARACTER (len = *), PARAMETER :: n_im_name = 'n_im'
    CHARACTER (len = *), PARAMETER :: n_re_name = 'n_re'
    CHARACTER (len = *), PARAMETER :: alpha_name = 'alpha'

    !file id
    INTEGER :: ncid

    !dimension ids
    INTEGER :: n_im_did, n_re_did, alpha_did

    !variable names
    CHARACTER (len = *), PARAMETER :: sigma_name = 'sigma'
    CHARACTER (len = *), PARAMETER :: asym_name = 'asym'
    CHARACTER (len = *), PARAMETER :: omega_name = 'omega'

    !variable ids
    INTEGER :: &
         n_im_vid, n_re_vid, alpha_vid, &
         sigma_vid, asym_vid, omega_vid



    !opening the file
    CALL check(nf90_open(filename, nf90_nowrite, ncid))
    !!nf90_open(filename, nf90_nowrite, ncid)

    !getting the dimension ids and then the dimension lengths
    !n_im
    CALL check(nf90_inq_dimid(ncid, n_im_name, n_im_did))
    CALL check(nf90_inquire_dimension(ncid, n_im_did, len=dn_im))

    !n_re
    CALL check(nf90_inq_dimid(ncid, n_re_name, n_re_did))
    CALL check(nf90_inquire_dimension(ncid, n_re_did, len=dn_re))

    !alpha
    CALL check(nf90_inq_dimid(ncid, alpha_name, alpha_did))
    CALL check(nf90_inquire_dimension(ncid, alpha_did, len=dalpha))


    !WRITE(*,*) n_im_name, n_im_did, dn_im
    !WRITE(*,*) n_re_name, n_re_did, dn_re
    !WRITE(*,*) alpha_name, alpha_did, dalpha

    !allocating memory:
    !dimensions
    ALLOCATE(n_im(dn_im))
    ALLOCATE(n_re(dn_re))
    ALLOCATE(alpha(dalpha))

    !variables
    ALLOCATE(lut_sigma(dn_re,dn_im,dalpha))
    ALLOCATE(lut_asym(dn_re,dn_im,dalpha))
    ALLOCATE(lut_omega(dn_re,dn_im,dalpha))

    !getting the variable ids
    CALL check(nf90_inq_varid(ncid, n_im_name, n_im_vid)) 
    CALL check(nf90_inq_varid(ncid, n_re_name, n_re_vid))
    CALL check(nf90_inq_varid(ncid, alpha_name, alpha_vid))
    CALL check(nf90_inq_varid(ncid, sigma_name, sigma_vid))
    CALL check(nf90_inq_varid(ncid, asym_name, asym_vid))
    CALL check(nf90_inq_varid(ncid, omega_name, omega_vid))

    !reading the data
    CALL check(nf90_get_var(ncid, n_im_vid, n_im))
    CALL check(nf90_get_var(ncid, n_re_vid, n_re))
    CALL check(nf90_get_var(ncid, alpha_vid, alpha))
    CALL check(nf90_get_var(ncid, sigma_vid, lut_sigma))
    CALL check(nf90_get_var(ncid, asym_vid, lut_asym))
    CALL check(nf90_get_var(ncid, omega_vid, lut_omega))

    !closing file
    CALL check(nf90_close(ncid))


  END SUBROUTINE init_aerorad_lut

END MODULE mo_aerorad_lut
