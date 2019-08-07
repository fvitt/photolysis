module jno_mod
  use phot_kind_mod, only: rk => kind_phot

  implicit none

  private
  public :: jno_nbins
  public :: jno_we
  public :: jno_init
  public :: jno_run
  public :: jno_timestep_init

  integer, parameter :: jno_nbins = 4       ! wavelength bins for MS, 93

  real(rk), parameter :: jno_we(jno_nbins+1) = (/ 181.6_rk, 183.1_rk, 184.6_rk, 190.2_rk, 192.5_rk /)

  real(rk) :: wtno50(6,2)
  real(rk) :: wtno90(6,2)
  real(rk) :: wtno100(6,2)
  real(rk) :: csno50(6,2)
  real(rk) :: csno90(6,2)
  real(rk) :: csno100(6,2)
  real(rk) :: etfphot_ms93(jno_nbins)  ! photons cm-2 sec-1 nm-1
  
contains

  subroutine jno_init
    !-------------------------------------------------------------
    !     ... Local variables
    !-------------------------------------------------------------
    real(rk), dimension(24) :: a, b, c

    !------------------------------------------------------------------------------
    !   	... 6 sub-intervals for O2 5-0 at 265K,
    !	    2 sub-sub-intervals for NO 0-0 at 250K
    !------------------------------------------------------------------------------
    a(:) = (/    0._rk,       0._rk,       0._rk,       0._rk, &
           5.12e-02_rk, 5.68e-03_rk, 1.32e-18_rk, 4.41e-17_rk, &
           1.36e-01_rk, 1.52e-02_rk, 6.35e-19_rk, 4.45e-17_rk, &
           1.65e-01_rk, 1.83e-02_rk, 7.09e-19_rk, 4.50e-17_rk, &
           1.41e-01_rk, 1.57e-02_rk, 2.18e-19_rk, 2.94e-17_rk, &
           4.50e-02_rk, 5.00e-03_rk, 4.67e-19_rk, 4.35e-17_rk /)

    !------------------------------------------------------------------------------
    !   	... sub-intervals for o2 9-0 band,
    !	    2 sub-sub-intervals for no 1-0 at 250 k
    !------------------------------------------------------------------------------
    b(:) = (/    0._rk,       0._rk,       0._rk,       0._rk, &
                 0._rk,       0._rk,       0._rk,       0._rk, &
           1.93e-03_rk, 2.14e-04_rk, 3.05e-21_rk, 3.20e-21_rk, &
           9.73e-02_rk, 1.08e-02_rk, 5.76e-19_rk, 5.71e-17_rk, &
           9.75e-02_rk, 1.08e-02_rk, 2.29e-18_rk, 9.09e-17_rk, &
           3.48e-02_rk, 3.86e-03_rk, 2.21e-18_rk, 6.00e-17_rk /)

    !------------------------------------------------------------------------------
    ! 	... sub-intervals for o2 10-0 band,
    !	    2 sub-sub-intervals for no 1-0 at 250 k
    !------------------------------------------------------------------------------
    c(:) = (/  4.50e-02_rk, 5.00e-03_rk, 1.80e-18_rk, 1.40e-16_rk, &
               1.80e-01_rk, 2.00e-02_rk, 1.50e-18_rk, 1.52e-16_rk, &
               2.25e-01_rk, 2.50e-02_rk, 5.01e-19_rk, 7.00e-17_rk, &
               2.25e-01_rk, 2.50e-02_rk, 7.20e-20_rk, 2.83e-17_rk, &
               1.80e-01_rk, 2.00e-02_rk, 6.72e-20_rk, 2.73e-17_rk, &
               4.50e-02_rk, 5.00e-03_rk, 1.49e-21_rk, 6.57e-18_rk /)

    wtno50 (1:6,1) = a(1:24:4)
    wtno50 (1:6,2) = a(2:24:4)
    csno50 (1:6,1) = a(3:24:4)
    csno50 (1:6,2) = a(4:24:4)
    wtno90 (1:6,1) = b(1:24:4)
    wtno90 (1:6,2) = b(2:24:4)
    csno90 (1:6,1) = b(3:24:4)
    csno90 (1:6,2) = b(4:24:4)
    wtno100(1:6,1) = c(1:24:4)
    wtno100(1:6,2) = c(2:24:4)
    csno100(1:6,1) = c(3:24:4)
    csno100(1:6,2) = c(4:24:4)

  end subroutine jno_init

  subroutine jno_timestep_init( etfphot_in )

    real(rk), intent(in) :: etfphot_in(:) ! photons cm-2 sec-1 nm-1

    ! (externally?) rebin etf to wavelength bins for MS, 93 (jno_we)
    etfphot_ms93(:) = etfphot_in(:)
    
  end subroutine jno_timestep_init

  subroutine jno_run( nlev, zen, n2vmr, o2vmr, o3vmr, novmr, press, temp, zkm, jno_out )
    use params_mod, only : kboltz, R, g
    use phot_util_mod, only: sphers
    
    integer,  intent(in) :: nlev
    real(rk), intent(in) :: zen
    real(rk), intent(in) :: n2vmr(nlev)
    real(rk), intent(in) :: o2vmr(nlev)
    real(rk), intent(in) :: o3vmr(nlev)
    real(rk), intent(in) :: novmr(nlev)
    real(rk), intent(in) :: press(nlev)
    real(rk), intent(in) :: temp(nlev)
    real(rk), intent(in) :: zkm(nlev)
    real(rk), intent(out) :: jno_out(nlev)

    real(rk) :: delz(nlev)               ! layer thickness (cm)
    real(rk), parameter    :: km2cm = 1.e5_rk
    
    real(rk) :: n2cc(nlev)
    real(rk) :: o2cc(nlev)
    real(rk) :: o3cc(nlev)
    real(rk) :: nocc(nlev)
    real(rk) :: dens(nlev)

    real(rk) :: o2scol(nlev)
    real(rk) :: o3scol(nlev)
    real(rk) :: noscol(nlev)
    real(rk) :: jno(nlev)

    
    integer  :: nid(0:nlev-1)
    real(rk) :: dsdh(0:nlev-1,nlev-1)
    real(rk) :: scaleh
    integer :: k, nlyr

    delz = 0._rk
    
    dens(:) = 10._rk * press(:) / ( kboltz * temp(:) ) ! molecules / cm3
    n2cc(:) = n2vmr(:) * dens(:)
    o2cc(:) = o2vmr(:) * dens(:)
    o3cc(:) = o3vmr(:) * dens(:)
    nocc(:) = novmr(:) * dens(:)

    !------------------------------------------------------------------------------
    !     ... Derive Slant Path for Spherical Atmosphere
    !------------------------------------------------------------------------------
    nlyr = nlev - 1
    call sphers( nlyr, zkm(nlev:1:-1), zen, dsdh, nid )

    !------------------------------------------------------------------------------
    !     ... Derive O2, O3, and NO Slant Column
    !------------------------------------------------------------------------------
    delz(1:nlev-1) = km2cm*(zkm(1:nlev-1) - zkm(2:nlev))
    scaleh = 10.e5_rk ! cm
    call slant_col( nlev, delz, dsdh, nid, o2cc, scaleh, o2scol )
    call slant_col( nlev, delz, dsdh, nid, o3cc, scaleh, o3scol )
    call slant_col( nlev, delz, dsdh, nid, nocc, scaleh, noscol )

    do k = 1,nlev
      write(10,fmt='(a,i2,a,e22.16)') ' o2scol(', k, '): ', o2scol(k)
    end do
 
    ! calc J rates for NO + hv -> N + O
    call jno_calc( nlev, etfphot_ms93, n2cc, o2scol, o3scol,  noscol, jno )

    jno_out(1:nlev) = jno(nlev:1:-1)

  end subroutine jno_run

  subroutine jno_calc( nlev, etfphot, n2cc, o2scol, o3scol,  noscol, jno )
    !-----------------------------------------------------------------------------!
    !   PURPOSE:                                                                  !
    !   Compute the total photolytic rate constant for NO in the SR bands         !
    !     - following the approach of Minshwanner and Siskind, JGR,               !
    !       98, D11, 20401-20412, 1993.                                           !
    !                                                                             !
    !-----------------------------------------------------------------------------!
    !   PARAMETERS:                                                               !
    !   NZ           - INTEGER, number of specified altitude levels               !
    !                                                                             !
    !   etfphot      - Extraterrestrial Flux, within the MS 1993 Grid             !
    !                  units of photons cm-2 sec-1 nm-1                           !
    !   n2cc         - N2 conc (molecules cm-3)                                   !
    !   o3scol       - Ozone Slant Column (molecules cm-2)                        !
    !   o2scol       - Oxygen Slant Column (molecules cm-2)                       !
    !   noscol       - Nitric Oxide Slant Column(molecules cm-2)                  !
    !                                                                             !
    !   LOCAL VARIABLES:                                                          !
    !   tauo3        - Transmission factor in the Hartley Band of O3              !
    !   etfphot_ms93 - Solar Irr. on Minschwaner and Siskind 1993 (MS93) Grid     !
    !   xs_o3ms93    - O3 cross section on the MS93 Grid                          !
    !                                                                             !
    !   OUTPUT VARIABLES:                                                         !
    !   jno          - photolytic rate constant                                   !
    !                  each specified altitude                                    !
    !                                                                             !
    !-----------------------------------------------------------------------------!
    !   EDIT HISTORY:                                                             !
    !   08/01  Created, Doug Kinnison, NCAR, ACD                                  !
    !-----------------------------------------------------------------------------!

    !------------------------------------------------------------------------------
    !       ... Dummy arguments
    !------------------------------------------------------------------------------
    integer, intent(in) :: nlev
    real(rk), intent(in)    :: etfphot(jno_nbins)
    real(rk), intent(in)    :: n2cc(nlev)
    real(rk), intent(in)    :: o3scol(nlev)
    real(rk), intent(in)    :: o2scol(nlev)
    real(rk), intent(in)    :: noscol(nlev)
    real(rk), intent(out)   :: jno(nlev)

    !------------------------------------------------------------------------------
    !	... Local variables
    !------------------------------------------------------------------------------
    integer     :: i, iw, lev
    real(rk)    :: jno50
    real(rk)    :: jno90
    real(rk)    :: jno100
    real(rk)    :: tauo3(nlev,jno_nbins)

    !------------------------------------------------------------------------------
    !   	... O3 SRB Cross Sections from WMO 1985, interpolated onto MS, 1993 grid
    !------------------------------------------------------------------------------
    real(rk), parameter :: xso3_ms93(jno_nbins) = (/ 7.3307600e-19_rk, 6.9660105E-19_rk, 5.9257699E-19_rk, 4.8372219E-19_rk /)

    !------------------------------------------------------------------------------
    !   	... delta wavelength of the MS, 1993 grid
    !------------------------------------------------------------------------------
    real(rk), parameter :: wlintv_ms93(jno_nbins) = (/ 1.50_rk, 1.50_rk, 5.6_rk, 2.3_rk /)

    !------------------------------------------------------------------------------
    !   	... O2 SRB Cross Sections for the six ODF regions, MS, 1993
    !------------------------------------------------------------------------------
    real(rk), parameter :: cs250(6)  = (/ 1.117e-23_rk, 2.447e-23_rk, 7.188e-23_rk, 3.042e-22_rk, 1.748e-21_rk, 1.112e-20_rk /)
    real(rk), parameter :: cs290(6)  = (/ 1.350e-22_rk, 2.991e-22_rk, 7.334e-22_rk, 3.074e-21_rk, 1.689e-20_rk, 1.658e-19_rk /)
    real(rk), parameter :: cs2100(6) = (/ 2.968e-22_rk, 5.831e-22_rk, 2.053e-21_rk, 8.192e-21_rk, 4.802e-20_rk, 2.655e-19_rk /)

    !------------------------------------------------------------------------------
    !     ... derive tauo3 for the three o2 srb
    !     ... iw = 1,2, and 4 are used below for jno
    !------------------------------------------------------------------------------
    do iw = 1,jno_nbins
       tauo3(:,iw) = exp( -xso3_ms93(iw)*o3scol(:) )
    end do

    !------------------------------------------------------------------------------
    !   	... Call PJNO Function to derive SR Band JNO contributions
    !         Called in order of wavelength interval (shortest firs)
    !------------------------------------------------------------------------------
    do lev = 1,nlev
       jno100   = pjno( 1, cs2100, wtno100, csno100 )
       jno90    = pjno( 2, cs290,  wtno90,  csno90 )
       jno50    = pjno( 4, cs250,  wtno50,  csno50 )
       jno(lev) = jno50 + jno90 + jno100
    end do

  contains

    function pjno( w, cso2, wtno, csno )
      !------------------------------------------------------------------------------
      !   	... uses xsec at center of g subinterval for o2
      !           uses mean values for no
      !------------------------------------------------------------------------------
      implicit none

      !------------------------------------------------------------------------------
      !	... parameters
      !------------------------------------------------------------------------------
      integer, parameter :: ngint = 6
      integer, parameter :: nno = 2

      !----------------------------------------------------------------
      !	... Dummy arguments
      !----------------------------------------------------------------
      integer, intent(in)     :: w
      real(rk),    intent(in) :: cso2(ngint)
      real(rk),    intent(in) :: csno(ngint,nno)
      real(rk),    intent(in) :: wtno(ngint,nno)

      !----------------------------------------------------------------
      !	... Function declarations
      !----------------------------------------------------------------
      real(rk) :: pjno

      !----------------------------------------------------------------
      !	... Local variables
      !----------------------------------------------------------------
      integer  ::  jj, i, k
      real(rk) :: tauno
      real(rk) :: transno
      real(rk) :: transo2
      real(rk) :: tauo2
      real(rk) :: jno
      real(rk) :: jno1

      !----------------------------------------------------------------
      !	... derive the photolysis frequency for no within a given
      !         srb (i.e., 5-0, 9-0, 10-0)
      !----------------------------------------------------------------
      jno = 0._rk
      do k = 1,ngint
	 tauo2 = o2scol(lev) * cso2(k)
	 if( tauo2 < 50._rk ) then
	    transo2 = exp( -tauo2 )
	 else
	    transo2 = 0._rk
	 end if
         jno1 = 0._rk
         do jj = 1,nno
            tauno = noscol(lev)*csno(k,jj)
	    if( tauno < 50._rk ) then
	       transno = exp( -tauno )
	    else
	       transno = 0._rk
	    end if
            jno1 = jno1 + csno(k,jj) * wtno(k,jj) * transno
         end do
         jno = jno + jno1*transo2
      end do

      pjno = wlintv_ms93(w)*etfphot(w)*tauo3(lev,w)*jno

      !----------------------------------------------------------------
      !	... correct for the predissociation of the deltq 1-0
      !         transition in the srb (5-0)
      !----------------------------------------------------------------
      if( w == 4 ) then
         pjno = 1.65e9_rk/(5.1e7_rk + 1.65e9_rk + (1.5e-9_rk*n2cc(nlev-lev+1)))*pjno
      end if

    end function pjno

  end subroutine jno_calc
  
      subroutine slant_col( nlev, delz, dsdh, nid, absden, hscale, scol )
!=============================================================================!
!   PURPOSE:                                                                  !
!   Derive Column
!=============================================================================!
!   PARAMETERS:                                                               !
!   NLEV   - INTEGER, number of specified altitude levels in the working  (I) !
!            grid                                                             !
!   DELZ   - REAL, specified altitude working grid (km)                   (I) !
!   DSDH   - REAL, slant path of direct beam through each layer crossed  (O)  !
!             when travelling from the top of the atmosphere to layer i;      !
!             DSDH(i,j), i = 0..NZ-1, j = 1..NZ-1                             !
!   NID    - INTEGER, number of layers crossed by the direct beam when   (O)  !
!             travelling from the top of the atmosphere to layer i;           !
!             NID(i), i = 0..NZ-1                                             !
!            specified altitude at each specified wavelength                  !
!   absden - REAL, absorber concentration, molecules cm-3                     !
!   SCOL   - REAL, absorber Slant Column, molecules cm-2                      !
!=============================================================================!
!   EDIT HISTORY:                                                             !
!   09/01  Read in profile from an input file, DEK                            !
!   01/02  Taken from Sasha Madronich's TUV code                              !
!=============================================================================!
       use phot_util_mod, only: airmas
!------------------------------------------------------------------------------
!       ... Dummy arguments
!------------------------------------------------------------------------------
      integer, intent(in) :: nlev
      integer, intent(in) :: nid(0:nlev-1)                ! see above
      real(rk), intent(in)    :: delz(nlev)	        ! layer thickness (cm)
      real(rk), intent(in)    :: dsdh(0:nlev-1,nlev-1)	! see above
      real(rk), intent(in)    :: absden(nlev)           ! absorber concentration (molec. cm-3)
      real(rk), intent(in)    :: hscale
      real(rk), intent(out)   :: scol(nlev)		! absorber Slant Column (molec. cm-2)

!------------------------------------------------------------------------------
!       ... Local variables
!------------------------------------------------------------------------------
      real(rk), parameter :: largest = 1.e+36_rk

      real(rk) :: sum
      real(rk) :: numer, denom
      real(rk) :: cz(nlev)
      real(rk) :: vcol(nlev)

      integer :: id
      integer :: j
      integer :: k
      integer :: nlyr

!------------------------------------------------------------------------------
!     ... compute column increments (logarithmic integrals)
!------------------------------------------------------------------------------
      do k = 1,nlev-1
	if( absden(k) /= 0._rk .and. absden(k+1) /= 0._rk ) then
           cz(nlev-k) = (absden(k) - absden(k+1))/log( absden(k)/absden(k+1) ) * delz(k)
	else
           cz(nlev-k) = .5_rk*(absden(k) + absden(k+1)) * delz(k)
	end if
      end do

!------------------------------------------------------------------------------
!     ... Include exponential tail integral from infinity to model top
!         specify scale height near top of model
!------------------------------------------------------------------------------
      cz(nlev-1) = cz(nlev-1) + hscale * absden(1)

      nlyr = nlev-1
      vcol = 0._rk
      scol = 0._rk
      call airmas( nlyr, dsdh, nid, cz(:nlyr), vcol(:nlyr), scol(:nlyr))

      scol(nlev) = .95_rk*scol(nlev-1)

      end subroutine slant_col


end module jno_mod
