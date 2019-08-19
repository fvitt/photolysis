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

  subroutine jno_run( nlev, zen, n2vmr, o2vmr, o3vmr, novmr, press, temp, alt, jno_out )
    use phot_util_mod, only : sphers, airmas
    use params_mod,    only : kboltz, R, g
    
    integer,  intent(in) :: nlev
    real(rk), intent(in) :: zen
    real(rk), intent(in) :: n2vmr(nlev)
    real(rk), intent(in) :: o2vmr(nlev)
    real(rk), intent(in) :: o3vmr(nlev)
    real(rk), intent(in) :: novmr(nlev)
    real(rk), intent(in) :: press(nlev)
    real(rk), intent(in) :: temp(nlev)
    real(rk), intent(in) :: alt(nlev)
    real(rk), intent(out) :: jno_out(nlev)
    
    real(rk) :: dpress(nlev)
    real(rk) :: n2cc(nlev)
    real(rk) :: aircol(nlev)
    real(rk) :: o2col(nlev)
    real(rk) :: o3col(nlev)
    real(rk) :: nocol(nlev)
    real(rk) :: o2scol(nlev)
    real(rk) :: o3scol(nlev)
    real(rk) :: noscol(nlev)
    real(rk) :: jno(nlev)
    real(rk) :: zlev(nlev)

    integer  :: nid(0:nlev)
    real(rk) :: dsdh(0:nlev,nlev)
    real(rk) :: vcol(nlev)

    integer :: k

    n2cc(nlev:1:-1) = n2vmr(1:nlev) * 10._rk * press(1:nlev) / ( kboltz * temp(1:nlev) )
    
    ! compute O2, O3 slant columns

    dpress(1) = press(1)
    dpress(2:nlev) = press(2:nlev) - press(1:nlev)

    aircol(:) = 10._rk*dpress(:)*R/(kboltz*g)

    o2col(1) = o2vmr(1)*aircol(1)
    o3col(1) = o3vmr(1)*aircol(1)
    nocol(1) = novmr(1)*aircol(1)

    do k = 2,nlev
       o2col(k) = 0.5_rk*(o2vmr(k)+o2vmr(k-1))*aircol(k)
       o3col(k) = 0.5_rk*(o3vmr(k)+o3vmr(k-1))*aircol(k)
       nocol(k) = 0.5_rk*(novmr(k)+novmr(k-1))*aircol(k)
    end do
    
    aircol(1:nlev) = aircol(nlev:1:-1)
    o2col(1:nlev)  = o2col(nlev:1:-1)
    o3col(1:nlev)  = o3col(nlev:1:-1)
    nocol(1:nlev)  = nocol(nlev:1:-1)
    zlev(1:nlev)  = alt(nlev:1:-1)*1.e-3_rk ! m -> km

    nid = 0
    dsdh = 0._rk
    call sphers( nlev, zlev, zen, dsdh, nid )

    ! calc slant columns above each level
    call airmas( nlev, dsdh, nid, o2col, vcol, o2scol )
    call airmas( nlev, dsdh, nid, o3col, vcol, o3scol )
    call airmas( nlev, dsdh, nid, nocol, vcol, noscol )

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
         pjno = 1.65e9_rk/(5.1e7_rk + 1.65e9_rk + (1.5e-9_rk*n2cc(lev)))*pjno
      end if

    end function pjno

  end subroutine jno_calc

end module jno_mod
