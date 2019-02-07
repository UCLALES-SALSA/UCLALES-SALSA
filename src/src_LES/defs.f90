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
!
MODULE defs

  INTEGER :: nv = 1, nv1 = 2, mb = 1
  REAL    :: totalpower
  REAL, PARAMETER :: R      = 287.04
  REAL, PARAMETER :: Rm     = 461.5
  REAL, PARAMETER :: ep     = R/Rm
  REAL, PARAMETER :: ep2    = Rm/R - 1.
  REAL, PARAMETER :: cp     = 1005.
  REAL, PARAMETER :: cv     = cp-R
  REAL, PARAMETER :: rcp    = R/cp
  REAL, PARAMETER :: cpr    = cp/R
  REAL, PARAMETER :: g      = 9.81
  REAL, PARAMETER :: p00    = 1.e+05
  REAL, PARAMETER :: p00i   = 1./p00
  REAL, PARAMETER :: omega  = 7.292e-05
  REAL, PARAMETER :: alvl   = 2.5e+06  ! latent heat of vaporization
  REAL, PARAMETER :: alvi   = 2.834e+06 ! latent heat of sublimation
  REAL, PARAMETER :: rowt   = 1.e+3
  REAL, PARAMETER :: roice  = 0.9e+3
  REAL, PARAMETER :: vonk   = 0.40
  REAL, PARAMETER :: stefan = 5.6696e-8
  REAL, PARAMETER :: SolarConstant  = 1.365d+3
  REAL, PARAMETER :: mair   = 28.967 ! molar mass of air
  REAL, PARAMETER :: pi     = 3.14159265358979323846264338327
  REAL, PARAMETER :: kb     = 1.3807e-23 !Boltzman constant

END MODULE defs
