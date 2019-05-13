      module PARAMS_MOD

      use phot_kind_mod, only: rk => kind_phot

      implicit none

      real(rk), parameter :: m2km = .001_rk                          ! meters to km
      real(rk), parameter :: ppm2vmr = 1.e-6_rk                      ! ppm to vmr
      real(rk), parameter :: o2vmr = .2095_rk                        ! o2 vmr
      real(rk), parameter :: km2cm = 1.e5_rk                         ! km to centimeters
      real(rk), parameter :: m2s   = 60._rk                          ! minutes to seconds

      REAL(rk), PARAMETER :: pi = 3.1415926535898_rk
      REAL(rk), PARAMETER :: radius = 6.371E+3_rk                    ! km
      REAL(rk), PARAMETER :: hc = 6.626068E-34_rk * 2.99792458E8_rk
      REAL(rk), PARAMETER :: largest=1.E+36_rk
      real(rk), parameter :: kboltz= 1.38064852e-16_rk ! boltzmann constant (erg/K)
      real(rk), parameter :: R=2.8704e6_rk       ! gas constant (erg/g/K)
      real(rk), parameter :: g=980.616_rk        ! grav acceleration (cm/sec2)


      REAL(rk), PARAMETER :: pzero = +10._rk/largest
      REAL(rk), PARAMETER :: nzero = -10._rk/largest

      REAL(rk), PARAMETER :: precis = 1.e-7_rk

      real(rk) :: lambda_cutoff                        ! nm

      character(len=256) :: input_data_root = 'NOT_SET'

      end module PARAMS_MOD
