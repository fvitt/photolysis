      module module_params

      use phot_kind_mod, only: rk => kind_phot

      implicit none

! BROADLY USED PARAMETERS:
!_________________________________________________
! i/o file unit numbers
! output
!      INTEGER, PARAMETER :: kout=53
! input
      INTEGER, PARAMETER :: kin=12
!_________________________________________________
! altitude, wavelength, time (or solar zenith angle) grids
! altitude
      integer, PARAMETER :: kz=125
! wavelength
      integer, PARAMETER :: kw=1000
! time/sza
!      integer, PARAMETER :: kt=100
!_________________________________________________
! number of weighting functions
!  wavelength dependent
      integer, PARAMETER :: ks=60
!  wavelength and altitude dependent
      integer, PARAMETER :: kj=150
!  wavelength dependent DOM (dissolved organic matter) spectra
      integer, PARAMETER :: kdom=200
! delta for adding points at beginning or end of data grids
      real(rk), PARAMETER :: deltax = 1.E-5_rk

! some constants...

! pi:
      real(rk), PARAMETER :: pi=3.1415926535898_rk

! radius of the earth, km:
      real(rk), PARAMETER :: radius=6.371E+3_rk

! Planck constant x speed of light, J m
      real(rk), PARAMETER :: hc = 6.626068E-34_rk * 2.99792458E8_rk

! largest number of the machine:
      real(rk), PARAMETER :: largest=1.E+36_rk

! small numbers (positive and negative)
      real(rk), PARAMETER :: pzero = +10._rk/largest
      real(rk), PARAMETER :: nzero = -10._rk/largest

! machine precision
      real(rk), PARAMETER :: precis = 1.e-7_rk

      end module module_params
