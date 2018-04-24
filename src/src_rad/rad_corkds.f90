!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 1999-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE ckd

   USE defs, ONLY : mair, nv, nv1, mb, totalpower
   IMPLICIT NONE
   PRIVATE

   CHARACTER (len=20) :: gasfile = 'datafiles/ckd.dat'
   LOGICAL, SAVE      :: Initialized = .FALSE.
   REAL, SAVE         :: bllmx, brlmn
   INTEGER, SAVE      :: ngases

   TYPE ckd_properties
      CHARACTER (len=5) :: name
      INTEGER           :: ng, nt, np, noverlap, iband
      REAL              :: mweight, default_conc, tbase
      REAL, ALLOCATABLE :: hk(:),sp(:),xk(:,:,:,:)
   END TYPE ckd_properties

   TYPE band_properties
      PRIVATE
      INTEGER              :: kg, ngases
      REAL                 :: llimit, rlimit, center, power = 0.
      INTEGER, ALLOCATABLE :: gas_id(:)
      REAL, ALLOCATABLE    :: hk(:)
   END TYPE band_properties

   TYPE (ckd_properties),  ALLOCATABLE :: gas(:)
   TYPE (band_properties), ALLOCATABLE :: band(:), solar_bands(:), ir_bands(:)

   PUBLIC :: band_properties, bllmx, brlmn, &
             band, solar_bands, ir_bands,   &
             init_ckd, gases,               &
             kg, llimit, rlimit, center, power, gPointWeight, isSolar
CONTAINS
   !
   ! ---------------------------------------------------------------------------
   ! Subroutine ckd_init:  Reads the correlated K distribution data an assures
   ! that it conforms to expected properties
   !
   SUBROUTINE init_ckd
      USE mpi_interface, ONLY : myid
      IMPLICIT NONE

      INTEGER :: i, j, k, l, n, ib, ii, mbs, mbir
      LOGICAL :: check
      REAL    :: bllmx, brlmn
      REAL, DIMENSION(2000) :: REALVars

      INTEGER, ALLOCATABLE :: gasesinband(:)

      OPEN( unit = 66, file = trim(gasfile), status = 'old' )
      READ(66, '(2I5)') mb, ngases
      ALLOCATE (band(mb), gas(ngases))
      READ(66, '(300(6E12.4,/))') REALVars(1:mb+1)
      DO ib = 1, mb
         band(ib)%llimit = REALVars(ib)
         band(ib)%rlimit = REALVars(ib+1)
      END DO
      READ(66, '(300(6E12.4,/))') REALVars(1:mb)
      DO ib = 1, mb
         band(ib)%power = REALVars(ib)
      END DO

      mbs = 0
      DO ib = 1, mb-1
         band(ib)%center =(band(ib)%rlimit+band(ib)%llimit)*0.5
         IF (band(ib)%power > 0.) mbs = mbs + 1
      END DO
      band(mb)%center = (band(mb)%rlimit+band(mb)%llimit)*0.5

      IF (band(mb)%power > 0.) mbs = mbs + 1
      mbir = mb - mbs
      IF (myid == 0) PRINT 600, trim(gasfile), ngases, mb, mbs, mbir, sum(band%power)

      DO n = 1, ngases
         READ(66,'(A5,I4)') gas(n)%name,gas(n)%iband
         READ(66,'(4I4)') gas(n)%noverlap,gas(n)%ng,gas(n)%np,gas(n)%nt
         READ(66,'(3E13.5)') gas(n)%mweight, gas(n)%default_conc, gas(n)%tbase
         ALLOCATE (gas(n)%hk(gas(n)%ng))
         ALLOCATE (gas(n)%sp(gas(n)%np))
         ALLOCATE (gas(n)%xk(gas(n)%nt,gas(n)%np,gas(n)%ng,gas(n)%noverlap))
         READ(66, '(300(6E12.4,/))') gas(n)%hk(1:gas(n)%ng)
         READ(66,'(300(6E12.4,/))') gas(n)%sp(1:gas(n)%np)
         READ(66,'(300(6E12.4,/))')                                            &
            REALVars(:(gas(n)%ng * gas(n)%np * gas(n)%nt * gas(n)%noverlap))
         ib = 1
         DO l = 1, gas(n)%noverlap
            DO k = 1, gas(n)%nt
               DO j = 1, gas(n)%np
                  DO i = 1, gas(n)%ng
                     gas(n)%xk(k,j,i,l) = REALVars(ib)
                     ib = ib + 1
                  END DO
               END DO
            END DO
         END DO 

         IF (abs(sum(gas(n)%hk) - 1.) <= 1.1 * spacing(1.) ) THEN
            IF (myid == 0) PRINT 601, gas(n)%name, gas(n)%iband, gas(n)%noverlap, &
               gas(n)%ng, gas(n)%np, gas(n)%nt
         ELSE
            PRINT *, gas(n)%hk, sum(gas(n)%hk(:))
            STOP 'TERMINATING: gas did not occur with probability one in band'
         END IF
      END DO

      bllmx = tiny(1.)
      brlmn = huge(1.)
      ALLOCATE (gasesinband(ngases))
      totalpower = 0.
      DO ib = 1, mb
         check = .FALSE.
         ii = 0
         DO n = 1, ngases
            IF (gas(n)%iband == ib) THEN
               check = .TRUE.
               ii = ii+1
               gasesinband(ii) = n
               IF (gas(n)%ng > 1 .AND. .NOT. allocated(band(ib)%hk)) THEN
                  ALLOCATE(band(ib)%hk(gas(n)%ng))
                  band(ib)%hk = gas(n)%hk
                  band(ib)%kg = gas(n)%ng
               END IF
            END IF
         END DO

         IF (.NOT. check) STOP 'TERMINATING: Gases do not span bands'
         band(ib)%ngases = ii
         ALLOCATE(band(ib)%gas_id(ii))
         band(ib)%gas_id(:) = gasesinband(1:ii)
         IF (band(ib)%power < epsilon(1.)) THEN
            bllmx = max(band(ib)%llimit, bllmx)
            brlmn = min(band(ib)%rlimit, brlmn)
         END IF
         totalpower = totalpower+band(ib)%power
      END DO
      DEALLOCATE (gasesinband)

      !
      ! Make separate solar and IR spectra
      !
      ALLOCATE(solar_bands(mbs), ir_bands(mbir))
      i = 1; j = 1
      DO ib = 1, mb
         IF(isSolar(band(ib))) THEN
            IF(i > size(solar_bands)) STOP 'TERMINATING: mismatch in solar bands'
            solar_bands(i) = band(ib)
            i = i + 1
         ELSE
            IF(j > size(ir_bands))    STOP 'TERMINATING: mismatch in IR bands'
            ir_bands(j) = band(ib)
            j = j + 1
         END IF
      END DO
    
      DO ib = 1, mb
         IF (myid == 0) PRINT 602, ib, band(ib)%power, band(ib)%llimit, band(ib)%rlimit, &
            band(ib)%ngases, band(ib)%kg
      END DO
      IF (myid == 0) PRINT 604

      Initialized = .TRUE.

600   FORMAT ('-----------------------------------------------------------', &
              /3x,'Computing with file: ',A20,' containing ',I3,' gases',   &
              /3x,'Number of Bands (Total/Solar/IR): ',I3,'/',I3,'/',I3     &
              /3x,'Total Power in Solar Bands: ',F8.2)
601   FORMAT ('                                                -----------', &
              /3x,'Reading Gas: ',A5,', Band ',I3,', Overlap Type ',I1,     &
              /3x,'# of g-points/ P Levels/ T Coeffs :',3I3)
602   FORMAT ('-----------                                                ', &
              /3x,'Band: ',I3,': ',F8.2,' Wm^-2',', between ',F6.0,' and ',F6.0,&
              ' cm^-1',/3x,I3,' gase(s): and ',I3,' g-points')
604   FORMAT ('---------------------------------------- Finished band init ')

   END SUBROUTINE init_ckd
   !
   ! ---------------------------------------------------------------------------
   ! Subroutine gases:  given an atmospheric state, a band and g point number
   ! this routroutine calculates the optical depth at that g-point within that
   ! band for all the gases (nongray gaseous absorbers) in the band.
   !
   SUBROUTINE gases ( this_band, ig, pp, pt, ph, po, tg )

      IMPLICIT NONE

      TYPE(band_properties), &
         INTENT (in) :: this_band
      INTEGER, INTENT (in) :: ig
      REAL, INTENT(in)     :: pp(nv1), pt(nv), ph(nv), po(nv)
      REAL, INTENT (out)   :: tg(nv)

      REAL, DIMENSION (nv) :: fkg, fkga, fkgb, pq

      INTEGER :: k, n, igg, nn
      REAL    :: xfct

      IF (.NOT. Initialized) STOP 'TERMINATING:  ckd_gases not initialized'
      DO k = 1, nv
         tg(k) = 0.
      END DO
      !
      ! loop through gases finding the cumulative probability for the band and g
      ! point, weight this by the power if in the solar bands, for which power>0.
      !
      DO nn = 1, this_band%ngases

         n = this_band%gas_id(nn)
         igg = min(ig,gas(n)%ng)
         SELECT CASE(gas(n)%noverlap)
            CASE (1)
               igg = min(ig,gas(n)%ng)
               CALL qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
                        gas(n)%xk(1,1,igg,1), pp, pt, fkg )
               CALL select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
                               pp, ph, po, pq)
               xfct = (2.24e4/gas(n)%mweight) * 10./9.81
               DO k = 1, nv
                  tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
               END DO

            CASE (2)
               CALL qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
                        gas(n)%xk(1,1,igg,1), pp, pt, fkga )
               CALL qk (gas(n)%nt, gas(n)%np, gas(n)%sp, gas(n)%tbase,            &
                        gas(n)%xk(1,1,igg,2), pp, pt, fkgb )
               CALL select_gas(gas(n)%name, gas(n)%default_conc, gas(n)%mweight,  &
                               pp, ph, po, pq)
               DO k = 1, nv
                  fkg(k) = fkga(k) + pq(k) * fkgb(k)
               END DO
               CALL select_gas('  CO2', gas(n)%default_conc, gas(n)%mweight,       &
                               pp, ph, po, pq)
               xfct = (2.24e4/gas(n)%mweight) * 10./9.81
               DO k = 1, nv
                  tg(k) = tg(k) + fkg(k)*pq(k)*(pp(k+1)-pp(k))*xfct
               END DO

            CASE DEFAULT
               STOP 'TERMINATING: overlap type not supported'
         END SELECT
      END DO

   END SUBROUTINE gases
   !
   ! ---------------------------------------------------------------------------
   ! Subroutine select_gas determintes the mixing ratio to use in the optical
   ! depth calculation.  For simple overlaps, this amounts to selecting the
   ! correct input array, ph for water vapor, po for ozone. For fixed gases
   ! this converts a concentration to a mixing ratio given the molecular weight
   ! of the gas and some specified background concentration.  For C02 the
   ! mixing ratio is chosen so that if conc = 330.e-6, the pq*xfct in subroutine
   ! po = 0.5
   !
   SUBROUTINE select_gas (name, conc_x, mx, pp, ph, po, pq)
    
      CHARACTER (len=5), INTENT (in) :: name
      REAL, INTENT (in)  :: conc_x,mx
      REAL, INTENT (in)  :: pp(nv1), ph(nv), po(nv)
      REAL, INTENT (out) :: pq(nv)

      INTEGER :: k
      REAL    :: xx

      SELECT CASE(name)
         CASE ('  H2O')
            pq(:) = ph(:)
         CASE ('   O3')
            pq(:) = po(:)
         CASE ('  CO2')
            xx = conc_x/330.e-6 / ((2.24e4/mx) * 10./9.81)
            pq(:) = xx
         CASE ('OVRLP')
            DO k = 1, nv
               IF ( pp(k) >= 63.1 ) THEN
                  pq(k) = ph(k)
               ELSE
                  pq(k) = 0.0
               END IF
            END DO
         CASE DEFAULT
            xx = conc_x*(mx/mair)
            pq(:) = xx
      END SELECT

      RETURN
   END SUBROUTINE select_gas
   !
   ! ---------------------------------------------------------------------------
   ! Subroutine qk: interpolates the gasesous absorption coefficients in units
   ! of (cm-atm)**-1 to the given temperature and pressure in each layer
   ! following: ln k = a + b * ( t - tbase ) + c * ( t - tbase ) ** 2 in
   ! temperature and  linear interpolation in pressure.
   !
   SUBROUTINE qk (nt, np, stanp, tbase, coefki, pp, pt, fkg )

      IMPLICIT NONE

      INTEGER, INTENT (in) :: np, nt
      REAL, INTENT(in)     :: pp(nv1), pt(nv), coefki(nt,np), stanp(np), tbase
      REAL, INTENT(out)    :: fkg(nv)

      INTEGER :: i1, k
      REAL    :: x1, x2, y1, y2, xp, pmid

      IF (nt < 3) THEN
         DO k = 1, nv
            fkg(k) = coefki(1,1)
         END DO
      ELSE
         i1 = 1
         x1 = 0.
         xp = 0.
         DO k = 1, nv
            pmid = 0.5*(pp(k)+pp(k+1))
            DO WHILE ( pmid >= stanp(i1) .AND. i1 < np)
               i1 = i1 + 1
            END DO
            y1 = (pt(k)-tbase)
            y2 = y1 * y1
            x2 = exp (coefki(1,i1) + coefki(2,i1)*y1 + coefki(3,i1)*y2)
            IF (i1 > 1) THEN
               x1 = exp (coefki(1,i1-1) + coefki(2,i1-1)*y1 + coefki(3,i1-1)*y2)
               xp = stanp(i1-1)
            END IF
            fkg(k) = x1 + (x2-x1)*(pmid-xp)/(stanp(i1)-xp)
         END DO
      END IF
      RETURN
   END SUBROUTINE qk
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL INTEGER FUNCTION kg(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      kg = thisBand%kg
   END FUNCTION kg
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL REAL FUNCTION power(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      power = thisBand%power
   END FUNCTION power
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL REAL FUNCTION llimit(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      llimit = thisBand%llimit
   END FUNCTION llimit
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL REAL FUNCTION rlimit(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      rlimit = thisBand%rlimit
   END FUNCTION rlimit
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL REAL FUNCTION center(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      center = (thisBand%llimit + thisBand%rlimit) / 2.
   END FUNCTION center
   !
   ! ----------------------------------------------------------------------
   !
   ELEMENTAL REAL FUNCTION gPointWeight(thisBand, ig)
      TYPE(band_properties), INTENT(in) :: thisBand
      INTEGER,               INTENT(in) :: ig
    
      gPointWeight = thisBand%hk(ig)
   END FUNCTION gPointWeight
   !
   ! ----------------------------------------------------------------------
   ! Function isSolar:  Returns True if the band is in the Solar region
   !
   ELEMENTAL LOGICAL FUNCTION isSolar(thisBand)
      TYPE(band_properties), INTENT(in) :: thisBand
    
      isSolar = thisBand%power > 0.
   END FUNCTION isSolar

END MODULE ckd
