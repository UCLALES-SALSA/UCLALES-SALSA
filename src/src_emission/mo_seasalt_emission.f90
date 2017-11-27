MODULE mo_seasalt_emission
IMPLICIT NONE

CONTAINS
SUBROUTINE seasalt_emissions_lsce_salsa(kproma, kbdim, krow, pseaice, velo10m_salsa, mass_flux, numb_flux)  !(number of horiz. grid kproma, dimension for arrays, local latitude index, sea ice fraction, wind speed, mass flux at given radius, number flux at given radius)
    !     -----------------------------------------------------------------------
    !     
    !   Author:
    !   -------
    !   Michael Schulz 
    !   Laboratoire des Sciences du Climat et de l'Environnement / Saclay
    !   10.1.2002
    !
    !   Modifications:
    !   --------------
    !   Philip Stier, MPI-MET  (Adaption to the ECHAM/HAM structure)       2002
    !   Michael Schulz, LSCE   (Modified source coefficients)        02/08/2002
    !   Harri Kokkola, MPI-MET (Conversion to sectional approach)          2007
    !   Tommi Bergman, CSC     (Updated sectional seasalt scheme)    14/09/2009
    !   Antti Kukkurainen, FMI (Modified for LES integration)        28/06/2016
    !
    !   Purpose:
    !   --------
    !   Describe source flux of sea salt aerosol mass and number flux 
    !   as a function of wind speed
    !   for salsa approach (salsa subranges 2a and 3a: 50nm-10um)
    !
    !   Interface: TB ok
    !   ----------
    !
    !      Input
    !       kproma         number of horizontal grid points
    !       kbdim          dimension for arrays
    !       krow           local latitude index
    !       pseaice        sea ice fraction         [ % 0.0 (sea) - 1. (ice)]
    !       velo10m_salsa  wind speed at 10 m       [m/s]
    !
    !      Output
    !       mass_flux      mass flux at given radius     [kg m-2 s-1]
    !       numb_flux      number flux at given radius   [# m-2 s-1]
    !
    !
    !   Method:
    !   -------
    !
    !TB Seasalt emissions calculated on-line by Monahan(1986) parameterization with improvements by Gong(2003).
    !   
    !     

    USE mo_submctl,        ONLY: aerobins, pi6, &
                                 rhoss, nbins,  &
                                 in2a,fn2a, fn2b          



    !TB


    IMPLICIT NONE

    !--- Parameters:

    INTEGER,     INTENT(in)    :: kproma, kbdim, krow
    REAL,    INTENT(in)    :: pseaice(kbdim)
    REAL,    INTENT(in)    :: velo10m_salsa(kbdim,krow)
    REAL,    INTENT(out)   :: mass_flux(kbdim,fn2b)
    REAL,    INTENT(out)   :: numb_flux(kbdim,fn2b)

    !--- Local Variables:
    REAL, DIMENSION(kbdim,krow)  :: slf    ! sea land fraction        [ % 0.0 (sea) - 1. (land)]
    REAL, DIMENSION(kbdim,krow)  :: alake  ! lake fraction            [%]
    REAL    :: zseafrac(kbdim)             ! fraction of the gridcell covered by non-iced sea water [0.-1.]

    INTEGER     :: ii, jj                      ! Loop indices

    REAL    :: B,                        & ! Auxilliary variable
                   B2                          ! Auxilliary variable

    REAL    :: dfdr,dfdr_a,dfdr_b          ! Numerflux per dr 
    REAL    :: ra,                       & ! Lowerlimit of bin aerosol radius
                   rb,                       & ! Upperlimit of bin aerosol radius
                   deltadp,                  & ! Width of a bin
                   A                           ! Auxilliary variable
    REAL    :: theta                       ! Parameterization parameter
    real    :: GF                          ! Hygroscopic growthfactor(2.0 in Guelle et al.

!!!! refined bin loop parameters for integration of dfdr

    REAL    :: velo10m           ! Width of a bin
    REAL    :: vlolim(nbins)     !low volume limit for bins in regimes 1 and 2 [fxm]
    REAL    :: vhilim(nbins)     !high volume limit for bins in regimes 1 and 2 [fxm]
    REAL    :: dpmid(nbins)      !center point of bin

    ! initilization

    slf=0.0 
    alake=0.0
    theta=30.0
    GF=2.0
    numb_flux(:,:)=0.
    mass_flux(:,:)=0.

    !Define boundaries:

    vlolim(in2a:fn2a) = pi6*aerobins(in2a:fn2a)**3
    vhilim(in2a:fn2a-1) = pi6*aerobins(in2a+1:fn2a)**3
    vhilim(fn2a) = pi6*(aerobins(fn2a) + 2.*aerobins(fn2a-1))**3
    dpmid(in2a:fn2a) = (  (vhilim(in2a:fn2a) + vlolim(in2a:fn2a)) / (2.0*pi6) )**(1.0/3.0) 


    !--- Calculate fraction of the gridcell of non ice-covered water:

    zseafrac(1:kproma)=(1.-slf(1:kproma,krow)-alake(1:kproma,krow))*(1.-pseaice(1:kproma))
    zseafrac(1:kproma)=MAX(0.,MIN(zseafrac(1:kproma),1.))

    !corinna: avoid emission from land
    WHERE (slf(1:kproma,krow)+alake(1:kproma,krow) .GT. 0.95) zseafrac(1:kproma)=0.
    !velo=min(velo10m_salsa(ii,krow),20.)
    ! TB: Monahan
    ! subregion 2a first

   
!!!!-------------
!!!! subregion 2
!!!!-------------

    do jj=in2a,fn2a
       do ii=1,kproma
          velo10m=min(velo10m_salsa(ii,krow),32.)
!          velo10m=velo10m_salsa(ii,krow)
          !calculate the width for current bin

          deltadp = (vhilim(jj)**(1./3.)-vlolim(jj)**(1./3.))&
               /pi6**(1./3.)*0.5e6 

          ! all of the particles will have radius of the average particle
          ! in the parameterization the radius is in micrometers

          ! We calculate the integral over the range by using simple 
          ! trapezoidal rule
          !    dr80dr0=0.506r0^-0.024
          ! lower limit
          ra=GF*((vlolim(jj)/pi6)**(1./3.))*0.5e6
          ! upper limit
          rb=GF*((vhilim(jj)/pi6)**(1./3.))*0.5e6


          ! Calculate the mass flux [kg m-2 s-1]

          !r=ra, r := mean radius of bin
          if (ra>0.4)then

             !Monahan

             ! temporary variables
             B=(0.38-log10(ra))/0.65
             B2=B**2.0

             ! massflux/dr: Guelle er al 2003(Monahan et al.) eq. A1 
             dfdr_a=1.373*velo10m**3.41*ra**(-3.0)*              &
                  (1.0+0.057*ra**1.05)*10**(1.19*exp(-B2)) * &
                  vlolim(jj)*rhoss
          else ! ra<0.4 -> Gong
             !Gong
             A=4.7*(1+theta*ra)**(-0.017*ra**(-1.44))
             B=(0.433-log10(ra))/0.433
             B2=B**2.0
          
             dfdr_a=1.373*velo10m**3.41*ra**(-A)*&
                  (1.0+0.057*ra**3.45)*10**(1.607*exp(-B2)) * &
                  vlolim(jj)*rhoss
          end if
          if (rb>0.4)then

             B=(0.38-log10(rb))/0.65
             B2=B**2.0
          
          ! Calculate the mass flux [kg m-2 s-1]
          
          ! massflux/dr: Guelle er al 2003(Monahan et al.) eq. A1 
          dfdr_b=1.373*velo10m**3.41*rb**(-3.0)*              &
               (1.0+0.057*rb**1.05)*10**(1.19*exp(-B2)) *  &
               vhilim(jj)*rhoss
          else ! ra<0.4 -> Gong
             ! flux within the section: 
             A=4.7*(1+theta*rb)**(-0.017*rb**(-1.44))
             B=(0.433-log10(rb))/0.433
             B2=B**2.0
             dfdr_b=1.373*velo10m**3.41*rb**(-A)*&
                  (1.0+0.057*rb**3.45)*10**(1.607*exp(-B2)) * &
                  vlolim(jj)*rhoss
          end if

          dfdr=0.5*(dfdr_a+dfdr_b)*deltadp

          ! seasalt production only within sea gridpoints
          dfdr=dfdr*zseafrac(ii)

          !mass flux = number flux time mass of one particle
          mass_flux(ii,jj)=dfdr

          ! Calculate the number flux[# m-2 s-1]
          numb_flux(ii,jj)=dfdr/((dpmid(jj))**3.0*pi6*rhoss)
          
       end do
    end do



END SUBROUTINE seasalt_emissions_lsce_salsa

END MODULE mo_seasalt_emission
