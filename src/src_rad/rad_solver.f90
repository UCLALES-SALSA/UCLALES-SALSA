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
! Copyright 1999-2008, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
MODULE solver

   USE defs, ONLY : nv, nv1, pi
   IMPLICIT NONE
   !
   ! u denotes the Double-Gauss quadratures and weights (Sykes, 1951).
   ! the p0-p3 denote the legendre polynomials for x=u(1) to u(4)
   ! respectively.  p11d(4,4), p22d(4,4), and p33d(4,4) are defined as
   ! 0.5*p1d(i)*p1d(j), 0.5*p2d(i)*p2d(j), and 0.5*p3d(i)*p3d(j) resp.
   !
   REAL, PARAMETER :: u(4) = (/-0.7886752,-0.2113247,0.2113247,0.7886752/)
   REAL, PARAMETER :: p0d(4) = (/ 1., 1., 1., 1./)
   REAL, PARAMETER :: p1d(4) = (/-0.788675, -.211325,  .211325, .788675/)
   REAL, PARAMETER :: p2d(4) = (/ 0.433013, -.433013, -.433013, .433013/)
   REAL, PARAMETER :: p3d(4) = (/-.0433940,  .293394, -.293394, .043394/)
   REAL, PARAMETER, DIMENSION(4, 4) ::                       &
      p11d = reshape(                                         &
      (/ .311004E+00, .833334E-01,-.833334E-01,-.311004E+00,  &
      0.833334E-01, .223291E-01,-.223291E-01,-.833334E-01, &
      -.833334E-01,-.223291E-01, .223291E-01, .833334E-01, &
      -.311004E+00,-.833334E-01, .833334E-01, .311004E+00  /), (/ 4, 4 /))
   REAL, PARAMETER, DIMENSION(4, 4) ::                       &
      p22d = reshape(                                          &
      (/ .937501E-01,-.937501E-01,-.937501E-01, .937501E-01,  &
      -.937501E-01, .937501E-01, .937501E-01,-.937501E-01, &
      -.937501E-01, .937501E-01, .937501E-01,-.937501E-01, &
      0.937501E-01,-.937501E-01,-.937501E-01, .937501E-01  /), (/ 4, 4 /))
   REAL, PARAMETER, DIMENSION(4, 4) ::                       &
      p33d = reshape(                                         &
      (/ .941520E-03,-.636577E-02, .636577E-02,-.941520E-03,  &
      -.636577E-02, .430400E-01,-.430400E-01, .636577E-02, &
      0.636577E-02,-.430400E-01, .430400E-01,-.636577E-02, &
      -.941520E-03, .636577E-02,-.636577E-02, .941520E-03  /), (/ 4, 4 /))

CONTAINS

   ! **********************************************************************
   ! coefficient calculations for four first-order differential equations.
   !
   ! See the paper by Liou, Fu and Ackerman (1988) for the formulation of
   ! the delta-four-stream approximation in a homogeneous layer.
   ! **********************************************************************
   SUBROUTINE coefft(solar,w,w1,w2,w3,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2)

      LOGICAL, INTENT(in) :: solar
      REAL, INTENT(in)    :: w, w1, w2, w3, t0, t1, u0, f0
      REAL, INTENT(out)   :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

      INTEGER :: i, j
      REAL :: a(2,2,2), b(4,3), c(4,5), d(4), z(4), b1, c1, fw, fw1, fw2, &
              fw3, fw4, w0w, w1w, w2w, w3w, q1, q2, q3, zx, x, y, a2, b2,   &
              fq0, fq1, fq, dt

      x = 0.5 * w
      w0w = x
      w1w = x * w1
      w2w = x * w2
      w3w = x * w3
      IF ( solar ) THEN
         fw = u0 * u0
         q1 = - w1w * u0
         q2 =   w2w * ( 1.5 * fw - 0.5 )
         q3 = - w3w * ( 2.5 * fw - 1.5 ) * u0
         DO i = 1, 4
            c(i,5) = (w0w + q1*p1d(i) + q2*p2d(i) + q3*p3d(i))/u(i)
         END DO
      ELSE
         DO i = 1, 4
            c(i,5) = 1.0/ u(i)
         END DO
      END IF

      fq = 0.5 * w0w
      DO i = 3, 4
         DO j = 1, 4
            c(i,j) = fq + w1w*p11d(i,j) + w2w*p22d(i,j) + w3w*p33d(i,j)
            IF ( i == j ) THEN
               c(i,j) = ( c(i,j) - 1.0 ) / u(i)
            ELSE
               c(i,j) = c(i,j) / u(i)
            END IF
         END DO
      END DO

      b(1,1) = c(4,4) - c(4,1)
      b(1,2) = c(4,4) + c(4,1)
      b(2,1) = c(4,3) - c(4,2)
      b(2,2) = c(4,3) + c(4,2)
      b(3,1) = c(3,4) - c(3,1)
      b(3,2) = c(3,4) + c(3,1)
      b(4,1) = c(3,3) - c(3,2)
      b(4,2) = c(3,3) + c(3,2)
      b(1,3) = c(4,5) - c(1,5)
      b(2,3) = c(3,5) - c(2,5)
      b(3,3) = c(3,5) + c(2,5)
      b(4,3) = c(4,5) + c(1,5)
      ! **********************************************************************
      ! coefficient calculations for second order differential equations.
      ! **********************************************************************
      fw1 = b(1,1) * b(1,2)
      fw2 = b(2,1) * b(3,2)
      fw3 = b(3,1) * b(2,2)
      fw4 = b(4,1) * b(4,2)
      a(2,2,1) = fw1 + fw2
      a(2,1,1) = b(1,1) * b(2,2) + b(2,1) * b(4,2)
      a(1,2,1) = b(3,1) * b(1,2) + b(4,1) * b(3,2)
      a(1,1,1) = fw3 + fw4
      a(2,2,2) = fw1 + fw3
      a(2,1,2) = b(1,2) * b(2,1) + b(2,2) * b(4,1)
      a(1,2,2) = b(3,2) * b(1,1) + b(4,2) * b(3,1)
      a(1,1,2) = fw2 + fw4
      d(1) = b(3,2) * b(4,3) + b(4,2) * b(3,3) + b(2,3) / u0
      d(2) = b(1,2) * b(4,3) + b(2,2) * b(3,3) + b(1,3) / u0
      d(3) = b(3,1) * b(1,3) + b(4,1) * b(2,3) + b(3,3) / u0
      d(4) = b(1,1) * b(1,3) + b(2,1) * b(2,3) + b(4,3) / u0

      ! **********************************************************************
      ! coefficient calculations for fourth-order differential equations.
      ! **********************************************************************
      x = u0 * u0
      b1 = a(2,2,1) + a(1,1,1)
      c1 = a(2,1,1) * a(1,2,1) - a(1,1,1) * a(2,2,1)
      z(1) = a(2,1,1) * d(3) + d(4) / x - a(1,1,1) * d(4)
      z(2) = a(1,2,1) * d(4) - a(2,2,1) *d(3) + d(3) / x
      z(3) = a(2,1,2) * d(1) + d(2) / x - a(1,1,2) * d(2)
      z(4) = a(1,2,2) * d(2) - a(2,2,2) * d(1) + d(1) / x

      ! **********************************************************************
      ! fk1 and fk2 are the eigenvalues.
      ! **********************************************************************
      dt = t1 - t0
      x = sqrt ( b1 * b1 + 4.0 * c1 )
      fk1 = sqrt ( ( b1 + x ) * 0.5 )
      fk2 = sqrt ( ( b1 - x ) * 0.5 )
      fw = u0 * u0
      x = 1.0 / ( fw * fw ) - b1 / fw - c1
      fw = 0.5 * f0 / x
      z(1) = fw * z(1)
      z(2) = fw * z(2)
      z(3) = fw * z(3)
      z(4) = fw * z(4)
      z1(1) = 0.5 * ( z(1) + z(3) )
      z1(2) = 0.5 * ( z(2) + z(4) )
      z1(3) = 0.5 * ( z(2) - z(4) )
      z1(4) = 0.5 * ( z(1) - z(3) )
      a2 = ( fk1 * fk1 - a(2,2,1) ) / a(2,1,1)
      b2 = ( fk2 * fk2 - a(2,2,1) ) / a(2,1,1)
      x = b(1,1) * b(4,1) - b(3,1) * b(2,1)
      fw1 = fk1 / x
      fw2 = fk2 / x
      y = fw2 * ( b2 * b(2,1) - b(4,1) )
      zx = fw1 * ( a2 * b(2,1) - b(4,1) )
      a1(1,1) = 0.5 * ( 1 - y )
      a1(1,2) = 0.5 * ( 1 - zx )
      a1(1,3) = 0.5 * ( 1 + zx )
      a1(1,4) = 0.5 * ( 1 + y )
      y = fw2 * ( b(3,1) - b2 * b(1,1) )
      zx = fw1 * ( b(3,1) - a2 * b(1,1) )
      a1(2,1) = 0.5 * ( b2 - y )
      a1(2,2) = 0.5 * ( a2 - zx )
      a1(2,3) = 0.5 * ( a2 + zx )
      a1(2,4) = 0.5 * ( b2 + y )
      a1(3,1) = a1(2,4)
      a1(3,2) = a1(2,3)
      a1(3,3) = a1(2,2)
      a1(3,4) = a1(2,1)
      a1(4,1) = a1(1,4)
      a1(4,2) = a1(1,3)
      a1(4,3) = a1(1,2)
      a1(4,4) = a1(1,1)
	! Should we dig deeper into the commented segments? -Juha
      IF ( solar ) THEN
         !fq0 = 0.
         !fq1 = 0.
         !IF(t0/u0 < 300) fq0 = exp ( - t0 / u0 ) !Stopping denormal values AZ
         !IF(t0/u0 < 300) fq1 = exp ( - t1 / u0 )
         fq0 = exp ( - t0 / u0 )
         fq1 = exp ( - t1 / u0 )
      ELSE
         fq0 = 1.0
         !fq1 = 0.
         !IF(dt/u0 < 300) fq1 = exp ( - dt / u0 )
         fq1 = exp ( - dt / u0 )
      END IF
      !x = 0.
      !y = 0.
      !IF (fk1*dt < 300) x = exp ( - fk1 * dt )
      !IF (fk2*dt < 300) y = exp ( - fk2 * dt )
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )
      DO i = 1, 4
         zz(i,1) = z1(i) * fq0
         zz(i,2) = z1(i) * fq1
         aa(i,1,1) = a1(i,1)
         aa(i,2,1) = a1(i,2)
         aa(i,3,1) = a1(i,3) * x
         aa(i,4,1) = a1(i,4) * y
         aa(i,3,2) = a1(i,3)
         aa(i,4,2) = a1(i,4)
         aa(i,1,2) = a1(i,1) * y
         aa(i,2,2) = a1(i,2) * x
      END DO

   END SUBROUTINE coefft

   ! **********************************************************************
   ! In the limits of no scattering ( Fu, 1991 ), fk1 = 1.0 / u(3) and
   ! fk2 = 1.0 / u(4).
   ! **********************************************************************
   SUBROUTINE coefft0( solar,t0,t1,u0,f0,aa,zz,a1,z1,fk1,fk2)

      LOGICAL, INTENT(in) :: solar
      REAL, INTENT(in)    :: t0, t1, u0, f0
      REAL, INTENT(out)   :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fk1, fk2

      INTEGER :: i, j, k, jj
      REAL    :: x, y, dt, fw

      fk1 = 4.7320545
      fk2 = 1.2679491
	! This too? -Juha
      !y = 0.
      !IF( (t1-t0)/u0 < 300)  y = exp ( -   ( t1 - t0 ) / u0 ) !Stop denormal AZ
      y = exp ( - ( t1 - t0 ) / u0 )
      fw = 0.5 * f0
      DO i = 1, 4
         IF ( solar ) THEN
            z1(i) = 0.0
            zz(i,1) = 0.0
            zz(i,2) = 0.0
         ELSE
            jj = 5 - i
            z1(i) = fw / ( 1.0 + u(jj) / u0 )
            zz(i,1) = z1(i)
            zz(i,2) = z1(i) * y
         END IF
         DO j = 1, 4
            a1(i,j) = 0.0
            DO k = 1, 2
               aa(i,j,k) = 0.0
            END DO
         END DO
      END DO
      DO i = 1, 4
         j = 5 - i
         a1(i,j) = 1.0
      END DO
      dt = t1 - t0
	! Again? -Juha
      !x = 0.
      !y = 0.
      !IF (fk1*dt < 300)  x = exp ( - fk1 * dt ) !stop denormal AZ
      !IF (fk2*dt < 300) y = exp ( - fk2 * dt )
      x = exp ( - fk1 * dt )
      y = exp ( - fk2 * dt )
      aa(1,4,1) = y
      aa(2,3,1) = x
      aa(3,2,1) = 1.0
      aa(4,1,1) = 1.0
      aa(1,4,2) = 1.0
      aa(2,3,2) = 1.0
      aa(3,2,2) = x
      aa(4,1,2) = y

   END SUBROUTINE coefft0

   ! **********************************************************************
   ! In the solar band  asbs is the surface albedo, while in the infrared
   ! band asbs is  blackbody intensity emitted at the surface temperature
   ! times surface emissivity.  In this subroutine, the delta-four-stream
   ! is applied to nonhomogeneous atmospheres. See comments in subroutine
   ! 'qcfel' for array AB(13,4*n).
   ! **********************************************************************
   SUBROUTINE qccfe (solar,asbs,ee,t,w,w1,w2,w3,u0,f0,fk1,fk2,a4,g4,z4 )

      LOGICAL, INTENT (in)              :: solar
      REAL, INTENT (in)                 :: asbs, ee
      REAL, DIMENSION (nv), INTENT(in)  :: t, w, w1, w2, w3, u0, f0
      REAL, DIMENSION (nv), INTENT(out) :: fk1, fk2
      REAL, INTENT(out)                 :: a4(4,4,nv), g4(4,nv), z4(4,nv)
      INTEGER :: i, j, k, m18, m28, n4, kf, j1, j2, j3, i1
      INTEGER :: i2, i3, i8, m2, m1
      REAL    :: fu(4,4), wu(4), ab(13,4*nv), bx(4*nv), xx(4*nv)
      REAL    :: aa(4,4,2), zz(4,2), a1(4,4), z1(4), fw1, fw2, wn, v1, v2, v3
      n4 = 4*nv
      DO i = 1, n4
         DO j = 1, 13
            ab(j,i) = 0.0
         END DO
      END DO
      wn = w(1)
      IF ( wn <= 1.0e-4 ) THEN
         CALL coefft0(solar,0.0,t(1),u0(1),f0(1),aa,zz,a1,z1,fk1(1),fk2(1))
      ELSE
         IF ( wn >= 0.999999 ) wn = 0.999999
         CALL coefft(solar,wn,w1(1),w2(1),w3(1),0.0,t(1),u0(1),f0(1), &
                     aa,zz,a1,z1,fk1(1),fk2(1))
      END IF
      DO  i = 1, 4
         z4(i,1) = z1(i)
         DO  j = 1, 4
            a4(i,j,1) = a1(i,j)
         END DO
      END DO
      DO  i = 1, 2
         bx(i) = - zz(i+2,1)
         i8 = i + 8
         DO j = 1, 4
            ab(i8-j,j) = aa(i+2,j,1)
         END DO
      END DO
      DO i = 1, 4
         wu(i) = zz(i,2)
         DO j = 1, 4
            fu(i,j) = aa(i,j,2)
         END DO
      END DO
      DO k = 2, nv
         wn = w(k)
         !WRITE(*,*) k
         IF ( w(k) <= 1.0e-4 ) THEN
            CALL coefft0(solar,t(k-1),t(k),u0(k),f0(k),aa,zz,a1,z1,fk1(k),fk2(k))
         ELSE
            IF ( wn >= 0.999999 ) wn = 0.999999
            CALL coefft(solar,wn,w1(k),w2(k),w3(k),t(k-1),t(k),u0(k),f0(k),  &
                        aa,zz,a1,z1,fk1(k),fk2(k))
         END IF
         DO  i = 1, 4
            z4(i,k) = z1(i)
            DO  j = 1, 4
               a4(i,j,k) = a1(i,j)
            END DO
         END DO
         kf = k + k + k + k
         i1 = kf - 5
         i2 = i1 + 3
         j1 = kf - 7
         j2 = j1 + 3
         i3 = 0
         DO  i = i1, i2
            i3 = i3 + 1
            bx(i) = - wu(i3) + zz(i3,1)
            j3 = 0
            i8 = i + 8
            DO  j = j1, j2
               j3 = j3 + 1
               ab(i8-j,j) = fu(i3,j3)
            END DO
            j3 = 0
            DO j = j2 + 1, j2 + 4
               j3 = j3 + 1
               ab(i8-j,j) = - aa(i3,j3,1)
            END DO
         END DO
         DO  i = 1, 4
            wu(i) = zz(i,2)
            DO j = 1, 4
               fu(i,j) = aa(i,j,2)
            END DO
         END DO
      END DO
      IF ( solar ) THEN
         v1 = 0.2113247 * asbs
         v2 = 0.7886753 * asbs
         !v3 = 0.
         !IF (t(nv)/u0(1) < 300) v3 = asbs * u0(1) * f0(1) * exp ( - t(nv) / u0(1) ) !Stop denormal AZ
         v3 = asbs * u0(1) * f0(1) * exp ( - t(nv) / u0(1) )
      ELSE
         v1 = 0.2113247 * ( 1.0 - ee )
         v2 = 0.7886753 * ( 1.0 - ee )
         v3 = asbs
      END IF
      m1 = n4 - 1
      m2 = n4
      m18 = m1 + 8
      m28 = m2 + 8
      fw1 = v1 * wu(3)
      fw2 = v2 * wu(4)
      bx(m1) = - ( wu(1) - fw1 - fw2 - v3 )
      bx(m2) = - ( wu(2) - fw1 - fw2 - v3 )
      DO j = 1, 4
         j1 = n4 - 4 + j
         fw1 = v1 * fu(3,j)
         fw2 = v2 * fu(4,j)
         ab(m18-j1,j1) = fu(1,j) - fw1 - fw2
         ab(m28-j1,j1) = fu(2,j) - fw1 - fw2
      END DO
      CALL qcfel (ab, bx, xx)
      DO k = 1, nv
         j = k + k + k + k - 4
         DO i = 1, 4
            j = j + 1
            g4(i,k) = xx(j)
         END DO
      END DO

   END SUBROUTINE qccfe

   ! **********************************************************************
   ! 1. `qcfel' is the abbreviation of ` qiu constants for each layer'.
   ! 2. The inhomogeneous atmosphere is divided into n adjacent homogeneous
   !    layers where the  single scattering properties are constant in each
   !    layer and allowed to vary from one to another. Delta-four-stream is
   !    employed for each homogeneous layer. The boundary conditions at the
   !    top and bottom of the atmosphere,  together with  continuity condi-
   !    tions  at  layer interfaces lead to a system of algebraic equations
   !    from which 4*n unknown constants in the problom can be solved.
   ! 3. This subroutine is used for solving the 4*n unknowns of A *X = B by
   !    considering the fact that the coefficient matrix is a sparse matrix
   !    with the precise pattern in this special problom.
   ! 4. The method is not different in principle from the general scheme of
   !    Gaussian elimination with backsubstitution, but carefully optimized
   !    so as to minimize arithmetic operations.  Partial  pivoting is used
   !    to quarantee  method's numerical stability,  which will  not change
   !    the basic pattern of sparsity of the matrix.
   ! 5. Scaling special problems so as to make  its nonzero matrix elements
   !    have comparable magnitudes, which will ameliorate the stability.
   ! 6. a, b and x present A, B and X in A*X=B, respectively. and n4=4*n.
   ! 7. AB(13,4*n) is the matrix A in band storage, in rows 3 to 13; rows 1
   !    and 2 and other unset elements should be set to zero on entry.
   ! 8. The jth column of A is stored in the jth column of the array AB  as
   !    follows:
   !            AB(8+i-j,j) = A(i,j) for max(1,j-5) <= i <= min(4*n,j+5).
   !    Reversedly, we have
   !            A(ii+jj-8,jj) = AB(ii,jj).
   ! **********************************************************************
   SUBROUTINE qcfel(ab, b, x)

      REAL, INTENT(inout) :: ab(13,4*nv),b(4*nv)
      REAL, INTENT(out)   :: x(4*nv)

      INTEGER :: i, j, k, l, m, n1, n2, n3, n4, m1, m2, m3, m4
      INTEGER :: k44, n44, m18, m28, m38, m48, m1f, im1, i0m1, i0, i0f, ifq
      REAL    :: xx, yy, t, p

      n4 = 4*nv
      DO  k = 1, nv - 1
         k44 = 4 * k - 4
         DO  l= 1, 4

            m1 = k44 + l
            p = 0.0
            DO  i = 8, 14 - l
               IF ( abs ( ab(i,m1) ) > abs ( p ) ) THEN
                  p = ab(i,m1)
                  i0 = i
               END IF
            END DO

            i0m1 = i0 + m1
            m18 = m1 + 8
            IF ( i0 /= 8 ) THEN
               DO  j = m1, m1 + 8 - l
                  i0f = i0m1 - j
                  m1f = m18 - j
                  t = ab(i0f,j)
                  ab(i0f,j) = ab(m1f,j)
                  ab(m1f,j) = t
               END DO
               i0f = i0m1 - 8
               t = b(i0f)
               b(i0f) = b(m1)
               b(m1) = t
            END IF

            yy = ab(8,m1)
            ab(8,m1) = 1.0
            DO  j = m1 + 1, m1 + 8 - l
               m1f = m18 - j
               ab(m1f,j) = ab(m1f,j) / yy
            END DO

            b(m1) = b(m1) / yy
            DO i = 9, 14 - l
               xx = ab(i,m1)
               ab(i,m1) = 0.0
               im1 = i + m1
               DO  j = m1 + 1, m1 + 8 - l
                  ifq = im1 - j
                  m1f = m18 - j
                  ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
               END DO
               ifq = im1 - 8
               b(ifq) = b(ifq) - b(m1) * xx
            END DO

         END DO
      END DO

      n44 = n4 - 4
      DO l = 1, 3

         m1 = n44 + l
         p = 0.0
         DO i = 8, 12 - l
            IF ( abs ( ab(i,m1) ) > abs ( p ) ) THEN
               p = ab(i,m1)
               i0 = i
            END IF
         END DO

         i0m1 = i0 + m1
         m18 = m1 + 8
         IF( i0 /= 8 ) THEN
            DO  j = m1, m1 + 4 - l
               i0f = i0m1 - j
               m1f = m18 - j
               t = ab(i0f,j)
               ab(i0f,j) = ab(m1f,j)
               ab(m1f,j) = t
            END DO
            i0f = i0m1 - 8
            t = b(i0f)
            b(i0f) = b(m1)
            b(m1) = t
         END IF

         yy = ab(8,m1)
         ab(8,m1) = 1.0
         DO  j = m1 + 1, m1 + 4 - l
            m1f = m18 - j
            ab(m1f,j) = ab(m1f,j) / yy
         END DO

         b(m1) = b(m1) / yy
         DO  i = 9, 12 - l
            xx = ab(i,m1)
            ab(i,m1) = 0.0
            im1 = i + m1
            DO  j = m1 + 1, m1 + 4 - l
               ifq = im1 - j
               m1f = m18 - j
               ab(ifq,j) = ab(ifq,j) - ab(m1f,j) * xx
            END DO
            ifq = im1 - 8
            b(ifq) = b(ifq) - b(m1) * xx
         END DO
      END DO

      yy = ab(8,n4)
      ab(8,n4) = 1.0
      b(n4) = b(n4) / yy
      n3 = n4 - 1
      n2 = n3 - 1
      n1 = n2 - 1
      x(n4) = b(n4)
      x(n3) = b(n3) - ab(7,n4) * x(n4)
      x(n2) = b(n2) - ab(7,n3) * x(n3) - ab(6,n4) * x(n4)
      x(n1) = b(n1) - ab(7,n2) * x(n2) - ab(6,n3) * x(n3) - ab(5,n4) * x(n4)
      DO  k = 1, nv - 1
         m4 = 4 * ( nv - k )
         m3 = m4 - 1
         m2 = m3 - 1
         m1 = m2 - 1
         m48 = m4 + 8
         m38 = m3 + 8
         m28 = m2 + 8
         m18 = m1 + 8
         x(m4) = b(m4)
         DO  m = m4 + 1, m4 + 4
            x(m4) = x(m4) - ab(m48-m,m) * x(m)
         END DO
         x(m3) = b(m3)
         DO  m = m3 + 1, m3 + 5
            x(m3) = x(m3) - ab(m38-m,m) * x(m)
         END DO
         x(m2) = b(m2)
         DO  m = m2 + 1, m2 + 6
            x(m2) = x(m2) - ab(m28-m,m) * x(m)
         END DO
         x(m1) = b(m1)
         DO  m = m1 + 1, m1 + 7
            x(m1) = x(m1) - ab(m18-m,m) * x(m)
         END DO
      END DO

   END SUBROUTINE qcfel

   ! **********************************************************************
   ! In this subroutine, we incorporate a delta-function adjustment to
   ! account for the  forward  diffraction  peak in the context of the
   ! four-stream approximation ( Liou, Fu and Ackerman, 1988 ). w1(n),
   ! w2(n), w3(n), w(n), and t(n) are the adjusted PARAMETERs.
   ! **********************************************************************
   SUBROUTINE adjust (tt,ww,ww1,ww2,ww3,ww4,t,w,w1,w2,w3)

      REAL, DIMENSION (nv), INTENT (in) :: tt,ww,ww1,ww2,ww3,ww4
      REAL, DIMENSION (nv), INTENT (out):: t,w,w1,w2,w3

      INTEGER :: k
      REAL    :: tt0, f, fw, dt(nv)

      tt0 = 0.0;f = 0.0; fw = 0.0; dt = 0.0; k = 0
      DO  k = 1, nv
         f = ww4(k) / 9.0
         fw = 1.0 - f * ww(k)
         w1(k) = ( ww1(k) - 3.0 * f ) / ( 1.0 - f )
         w2(k) = ( ww2(k) - 5.0 * f ) / ( 1.0 - f )
         w3(k) = ( ww3(k) - 7.0 * f ) / ( 1.0 - f )
         w(k) = ( 1.0 - f ) * ww(k) / fw
         dt(k) = (tt(k) - tt0) * fw
         tt0 = tt(k)
      END DO

      t(1) = dt(1)
      DO k = 2, nv
         t(k) = dt(k) + t(k-1)
      END DO

   END SUBROUTINE adjust
   ! --------------------------------------------------------------------------
   ! Subroutine qft: Delta 4-stream solver for fluxes
   !
   SUBROUTINE qft (solar, ee, as, u0, bf, tt, ww, ww1, ww2, ww3, ww4, ffu, ffd)

      LOGICAL, INTENT (in) :: solar
      REAL, INTENT (in)    :: ee, as, u0
      REAL, DIMENSION (nv), INTENT (in)   :: tt,ww,ww1,ww2,ww3,ww4
      REAL, DIMENSION (nv1), INTENT (in)  :: bf
      REAL, DIMENSION (nv1), INTENT (out) :: ffu, ffd

      REAL, DIMENSION (nv) :: t,w,w1,w2,w3,u0a,f0a,fk1,fk2
      INTEGER :: k, kk, ii, jj
      REAL    :: x(4), fi(4), a4(4,4,nv), z4(4,nv), g4(4,nv)
      REAL    :: tkm1, fw3, fw4, y1, xy, xas, xee
      REAL, PARAMETER :: fw1 = 0.6638960, fw2 = 2.4776962
      t = 0.0; w = 0.0; w1 = 0.0; w2 = 0.0; w3 = 0.0; u0a = 0.0; f0a = 0.0; fk1 = 0.0; fk2 = 0.0;

      CALL adjust(tt,ww,ww1,ww2,ww3,ww4,t,w,w1,w2,w3)

      IF (solar) THEN
         fw3 = u0
         xee = 0.0
         xas = as
         DO k = 1, nv
            u0a(k) = u0
            f0a(k) = 1./pi
         END DO
      ELSE
         fw3 = 0.
         xas = bf(nv1) * ee
         xee = ee
         tkm1 = 0.0
         DO k = 1, nv
            f0a(k) = 2.0 * ( 1.0 - w(k) ) * bf(k)
            u0a(k) = -(t(k)-tkm1) / ( alog( bf(k+1)/bf(k) ) + epsilon(1.))
            tkm1 = t(k)
         END DO
      END IF
      CALL qccfe (solar,xas,xee,t,w,w1,w2,w3,u0a,f0a,fk1,fk2,a4,g4,z4 )

      tkm1 = 0.
      DO k = 1, nv1

         IF ( k == 1 ) THEN
            x(1) = 1.0
            x(2) = 1.0
            !x(3) = 0.0
            !x(4) = 0.0
            !IF (fk1(1)*t(1) < 300) x(3) = exp ( - fk1(1) * t(1) ) !stop denormal values AZ
            !IF (fk2(1)*t(1) < 300) x(4) = exp ( - fk2(1) * t(1) )
            x(3) = exp ( - fk1(1) * t(1) )
            x(4) = exp ( - fk2(1) * t(1) )
            kk = 1
            xy = 1.0
         ELSE
            kk = k - 1
            y1 = t(kk) - tkm1
            !x(1) = 0.0
            !x(2) = 0.0
            !IF (fk2(kk)*y1 < 300) x(1) = exp ( - fk2(kk) * y1 ) !stop denormal values
            !IF (fk1(kk)*y1 < 300) x(2) = exp ( - fk1(kk) * y1 )
            x(1) = exp ( - fk2(kk) * y1 )
            x(2) = exp ( - fk1(kk) * y1 )
            x(3) = 1.0
            x(4) = 1.0
            IF (solar) y1 = t(kk)
            !xy = 0.0
            !IF (y1/u0a(kk) < 300) xy =  exp ( - y1 / u0a(kk) ) !AZ
            xy =  exp ( - y1 / u0a(kk) )
         END IF
         IF (kk > 1) tkm1 = t(kk)

         DO  jj = 1, 4
            fi(jj) = z4(jj,kk) * xy
         END DO

         DO ii = 1, 4
            fw4 = g4(ii,kk) * x(ii)
            DO jj = 1, 4
               fi(jj) = fi(jj) + a4(jj,ii,kk) * fw4
            END DO
         END DO

         ffu(k)= fw1 * fi(2) + fw2 * fi(1)
         ffd(k)= fw1 * fi(3) + fw2 * fi(4)  + fw3 * xy
      END DO

   END SUBROUTINE qft

END MODULE solver
