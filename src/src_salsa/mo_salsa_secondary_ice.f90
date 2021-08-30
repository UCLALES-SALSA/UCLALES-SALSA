MODULE mo_salsa_secondary_ice
  USE classSection, ONLY : Section
  USE mo_submctl, ONLY : nice, pi6, spec
  IMPLICIT NONE

  REAL, PARAMETER :: c1 = 350.   ! Splinters generated per milligram of new rime
  REAL, PARAMETER :: tmax = 270.16, tmid = 268.16, tmin = 265.16
  REAL, PARAMETER :: Dsplint = 1.e-6  ! Assumed splinter diameter

  
  CONTAINS

    SUBROUTINE halletmossop(kbdim,kproma,klev,any_sizerange,ptemp,drimdt,ice)
      TYPE(Section), POINTER, INTENT(inout)  :: ice(:,:,:)
      INTEGER, INTENT(in) :: kbdim,kproma,klev
      LOGICAL, INTENT(in) :: any_sizerange(kbdim,klev) ! Logical whether liquid hydrometeor exist in the sizeranges required by the H-M process
      REAL, INTENT(in) :: ptemp(kbdim,klev)
      REAL, INTENT(in) :: drimdt(kbdim,klev,nice)    ! Volume change in rime due to droplets larger than 20 um

      INTEGER :: ii,jj,bb
      REAL :: dN  ! Number of splintered rime
      REAL :: dm  ! Mass of splintered rime
      INTEGER :: iwa, iri

      iwa = spec%getIndex("H2O")
      iri = spec%getIndex("rime")
      
      DO bb = 1,nice
         DO jj = 1,klev
            DO ii = 1,kproma

               IF ( drimdt(ii,jj,bb) > 0. .AND. any_sizerange(ii,jj) ) THEN
                  ! Number of generated splinters. The maximum rate is at 268.16 K
                  IF ( ptemp(ii,jj) < tmax .AND. ptemp(ii,jj) >= tmid ) THEN
                     
                     dN = c1 * (spec%rhori*drimdt(ii,jj,bb)*1.e-6) * (tmax - ptemp(ii,jj))/2.
                     
                  ELSE IF ( ptemp(ii,jj) <= tmid .AND. ptemp(ii,jj) > tmin ) THEN
                     
                     dN = c1 * (spec%rhori*drimdt(ii,jj,bb)*1.e-6) * (ptemp(ii,jj) - tmin)/3.
                     
                  END IF

                  ! This will assume that the splinters consist of frozen spheres 1 um in diameter. The mass for these
                  ! is taken from the "new" rime mass, which should be more than enough, since only the rime generated
                  ! by > 20 um droplets is considered here. Change the type of ice to pristine for the splinters.
                  dm = spec%rhori*dN*pi6*Dsplint**3
                  
                  ice(ii,jj,1)%numc = ice(ii,jj,1)%numc + dN
                  ice(ii,jj,1)%volc(iwa) = ice(ii,jj,1)%volc(iwa) + dm/spec%rhoic
                  ice(ii,jj,bb)%volc(iri) = ice(ii,jj,bb)%volc(iri) - dm/spec%rhori
               END IF
                  
            END DO
         END DO
      END DO
   
      
    END SUBROUTINE halletmossop
  

END MODULE mo_salsa_secondary_ice
