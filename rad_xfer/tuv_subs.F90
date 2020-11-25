      MODULE tuv_subs

      use phot_kind_mod, only: dp
      use phot_kind_mod, only: rk => kind_phot

      IMPLICIT none

      private
      public :: tuv_radfld, fery, futr, addpnt, inter2

      CONTAINS

      SUBROUTINE tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlev, nwave, &
                             zenith, z, albedo, &
                             aircol, o3col, so2col, no2col, &
                             dtaer, omaer, gaer, dtcld, omcld, gcld, &
                             qll, dobsi, o3_xs, no2_xs, o2_xs, &
                             so2_xs, wmin, wc, tlev, dto2, radfld, efld, &
                             e_dir, e_dn, e_up, &
                             dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )
!-----------------------------------------------------------------------------
!     ... calculate the radiation field
!-----------------------------------------------------------------------------
  
      use molec_ox_xsect,only : molec_ox_xsect_run
      use rad_trans,     only : rtlink
      use phot_util_mod, only : sphers, airmas
      
!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in)  :: nlambda_start
      integer, intent(in)  :: nlev
      integer, intent(in)  :: nwave
      integer, intent(in)  :: cld_od_opt
      real(rk), intent(in)  :: zenith
      real(rk), intent(in)  :: dobsi
      real(rk), intent(in)  :: wmin
      real(rk), intent(in)  :: z(:)
      real(rk), intent(in)  :: albedo(:)
      real(rk), intent(in)  :: aircol(:)
      real(rk), intent(in)  :: o3col(:)
      real(rk), intent(in)  :: so2col(:)
      real(rk), intent(in)  :: no2col(:)
      real(rk), intent(in)  :: qll(:)
      real(rk), intent(in)  :: wc(:)
      real(rk), intent(in)  :: tlev(:)
      real(rk), intent(in)  :: cldfrac(:)
      real(rk), intent(in)  :: o2_xs(:)
      real(rk), intent(in)  :: so2_xs(:)
      real(rk), intent(in)  :: o3_xs(:,:)
      real(rk), intent(in)  :: no2_xs(:,:)
      real(rk), intent(in)  :: dto2(:,:)
      real(rk), intent(out) :: radfld(:,:)
      real(rk), intent(out) :: efld(:,:)
      real(rk), intent(inout)  :: dir_fld(:,:), dwn_fld(:,:), up_fld(:,:)
      real(rk), intent(inout)  :: e_dir(:,:), e_dn(:,:), e_up(:,:)
      real(rk), intent(inout)  :: dt_cld(:)
      real(rk), intent(in)  :: dtaer(:,:), omaer(:,:), gaer(:,:)
      real(rk), intent(inout)  :: dtcld(:,:), omcld(:,:), gcld(:,:)

      character(len=*), intent(out)   :: errmsg
      integer,          intent(out)   :: errflg

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      integer :: wn
      integer :: n_radlev, n_radlevp1
      integer :: nid(0:nlev)
      real(rk) :: dtrl(nlev,nwave)
      real(rk) :: dto3(nlev,nwave)
      real(rk) :: dtso2(nlev,nwave)
      real(rk) :: dtno2(nlev,nwave)
!     real :: dtcld(nlev,nwave)
!     real :: dtaer(nlev,nwave)
      real(rk) :: dtsnw(nlev,nwave)

!     real :: omcld(nlev,nwave)
!     real :: gcld(nlev,nwave)
!     real :: omaer(nlev,nwave)
!     real :: gaer(nlev,nwave)
      real(rk) :: omsnw(nlev,nwave)
      real(rk) :: gsnw(nlev,nwave)

      real(rk) :: edir(nlev+1)
      real(rk) :: edn(nlev+1)
      real(rk) :: eup(nlev+1)
      real(rk) :: fdir(nlev+1)
      real(rk) :: fdn(nlev+1)
      real(rk) :: fup(nlev+1)
      real(rk) :: dsdh(0:nlev,nlev)

      errmsg = ' '
      errflg = 0
      
      n_radlev = size( radfld,dim=2 )
      n_radlevp1 = n_radlev + 1

      do wn = 1,nwave
        omcld(:,wn) = 0._rk
!        omaer(:,wn) = 0._rk
        omsnw(:,wn) = 0._rk
        gcld(:,wn)  = 0._rk
!        gaer(:,wn)  = 0._rk
        gsnw(:,wn)  = 0._rk
        dtcld(:,wn) = 0._rk
!        dtaer(:,wn) = 0._rk
        dtsnw(:,wn) = 0._rk
      end do

      call odrl( wc, aircol, dtrl )
      call odo3( o3col, o3_xs, dto3, dobsi )
      call setso2( so2col, so2_xs, dtso2 )
      call setno2( no2col, no2_xs, dtno2 )
!-------------------------------------------------------------
! aerosol optical depths
!-------------------------------------------------------------
!      if( has_aer_ra_feedback ) then
!        call setaer( nlambda_start, wc, tauaer300, tauaer400, &
!                     tauaer600, tauaer999, waer300, &
!                     waer400, waer600, waer999,     &
!                     gaer300, gaer400, gaer600,     &
!                     gaer999, dtaer, omaer, gaer )
!      endif
!-------------------------------------------------------------
! cloud optical depths (cloud water units = g/m3)
!
      call setcld( nlambda_start, cld_od_opt, z, qll, cldfrac, &
                   dtcld, omcld, gcld, errmsg, errflg )
      if (errflg .ne. 0) return
      
!      dt_cld(:n_radlev) = dtcld(2:n_radlevp1,1)

      call sphers( nlev, z, zenith, dsdh, nid )

      do wn = nlambda_start,nwave
        call rtlink( &
           nlev+1, nlev, nwave, &
           wn, albedo(wn), zenith, &
           dsdh, nid, &
           dtrl,  &
           dto3,  &
           dto2, &
           dtso2, &
           dtno2,  &
           dtcld, omcld, gcld, &
           dtaer, omaer, gaer, &
           dtsnw, omsnw, gsnw, &
           edir, edn, eup, fdir, fdn, fup, errmsg, errflg )
        
        radfld(wn,1:n_radlev) = fdir(1:n_radlev) + fdn(1:n_radlev) + fup(1:n_radlev)
        efld(1:n_radlev,wn)    = edir(1:n_radlev) + edn(1:n_radlev) + eup(1:n_radlev)
        dir_fld(1:n_radlev,wn) = fdir(:n_radlev)
        dwn_fld(1:n_radlev,wn) = fdn(1:n_radlev)
        up_fld(1:n_radlev,wn)  = fup(1:n_radlev)
        e_dir(1:n_radlev,wn)   = edir(1:n_radlev)
        e_dn(1:n_radlev,wn)    = edn(1:n_radlev)
        e_up(1:n_radlev,wn)    = eup(1:n_radlev)
        
      end do

      END SUBROUTINE tuv_radfld

      SUBROUTINE odrl( wc, aircol, dtrl )
!-----------------------------------------------------------------------------*
!=  PURPOSE:                                                                 =*
!=  Compute Rayleigh optical depths as a function of altitude and wavelength =*
!-----------------------------------------------------------------------------*
!=  PARAMETERS:                                                              =*
!=  C       - REAL, number of air molecules per cm^2 at each specified    (O)=*
!=            altitude layer                                                 =*
!=  DTRL    - REAL, Rayleigh optical depth at each specified altitude     (O)=*
!=            and each specified wavelength                                  =*
!-----------------------------------------------------------------------------*

!-----------------------------------------------------------------------------*
!     ...dummy arguments
!-----------------------------------------------------------------------------*
      REAL(rk),    intent(in)  :: aircol(:)
      REAL(rk),    intent(in)  :: wc(:)
      REAL(rk),    intent(out) :: dtrl(:,:)

!-----------------------------------------------------------------------------*
!     ...local variables
!-----------------------------------------------------------------------------*
      INTEGER :: nwave, nlyr
      INTEGER :: wn
      REAL(rk)    :: srayl, wmicrn, xx 
      
      nwave = size( wc )
      nlyr  = size( aircol )
!-----------------------------------------------------------------------------*
! compute Rayleigh cross sections and depths:
!-----------------------------------------------------------------------------*
      DO wn = 1,nwave
!-----------------------------------------------------------------------------*
! Rayleigh scattering cross section from WMO 1985 (originally from
! Nicolet, M., On the molecular scattering in the terrestrial atmosphere:
! An empirical formula for its calculation in the homoshpere, Planet.
! Space Sci., 32, 1467-1468, 1984.
!-----------------------------------------------------------------------------*
        wmicrn =  wc(wn)*1.E-3_rk
        IF( wmicrn <= 0.55_rk ) THEN
          xx = 3.6772_rk + 0.389_rk*wmicrn + 0.09426_rk/wmicrn
        ELSE
          xx = 4.04_rk
        ENDIF
        srayl = 4.02e-28_rk/(wmicrn)**xx
!-----------------------------------------------------------------------------*
! alternate (older) expression from
! Frohlich and Shaw, Appl.Opt. v.11, p.1773 (1980).
!-----------------------------------------------------------------------------*
        dtrl(:nlyr,wn) = aircol(:nlyr)*srayl
      END DO

      END SUBROUTINE odrl

      SUBROUTINE odo3( o3col, o3xs, dto3, dobsi )
!-----------------------------------------------------------------------------
!=  NAME:  Optical Depths of O3
!=  PURPOSE:
!=  Compute ozone optical depths as a function of altitude and wavelength
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  O3XS   - REAL, molecular absoprtion cross section (cm^2) of O3 at     (I)
!=           each specified wavelength and altitude
!=  C      - REAL, ozone vertical column increments, molec cm-2, for each (I)
!=           layer
!=  DTO3   - REAL, optical depth due to ozone absorption at each          (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)    :: dobsi
      REAL(rk), intent(in)    :: o3col(:)
      REAL(rk), intent(in)    :: o3xs(:,:)
      REAL(rk), intent(inout) :: dto3(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn
      REAL(rk)    :: dob_at_grnd, scale_fac

      nwave = size(o3xs,dim=1)
      nlyr  = size(o3col)

      if( dobsi == 0._rk ) then
!-----------------------------------------------------------------------------
!  no scaling
!-----------------------------------------------------------------------------
        DO wn = 1,nwave
          dto3(:nlyr,wn) = o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      else
!-----------------------------------------------------------------------------
!  scale model o3 column to dobsi
!-----------------------------------------------------------------------------
        dob_at_grnd = sum( o3col(:nlyr) )/2.687e16_rk
        scale_fac   = dobsi/dob_at_grnd
        DO wn = 1,nwave
          dto3(:nlyr,wn) = scale_fac * o3col(:nlyr) * o3xs(wn,:nlyr)
        END DO
      endif

      END SUBROUTINE odo3

     SUBROUTINE setso2( colso2, so2_xs, dtso2 )
!-----------------------------------------------------------------------------
!=  PURPOSE:
!=  Set up an altitude profile of SO2 molecules, and corresponding absorption
!=  optical depths.  Subroutine includes a shape-conserving scaling method
!=  that allows scaling of the entire profile to a given overhead SO2
!=  column amount.
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  SO2_XS - REAL, molecular absoprtion cross section (cm^2) of O2 at     (I)
!=           each specified wavelength
!=  DTSO2  - REAL, optical depth due to SO2 absorption at each            (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk),    intent(in)  :: colso2(:)
      REAL(rk),    intent(in)  :: so2_xs(:)
      REAL(rk),    intent(out) :: dtso2(:,:)

!-----------------------------------------------------------------------------
!     ... local variables
!-----------------------------------------------------------------------------
      integer :: nwave, nlyr
      integer :: wn

      nwave = size( so2_xs )
      nlyr  = size( colso2 )

      DO wn = 1,nwave
        dtso2(:nlyr,wn) = colso2(:nlyr)*so2_xs(wn)
      END DO

      END SUBROUTINE setso2

      SUBROUTINE setno2( colno2, no2_xs, dtno2 )
!-----------------------------------------------------------------------------
!=  NAME:  Optical Depths of no2
!=  PURPOSE:
!=  Compute no2 optical depths as a function of altitude and wavelength
!-----------------------------------------------------------------------------
!=  PARAMETERS:
!=  NO2_XS - REAL, molecular absoprtion cross section (cm^2) of no2 at    (I)
!=           each specified wavelength and altitude
!=  COLNO2 - REAL, no2 vertical column increments, molec cm-2, for each   (I)
!=           layer
!=  DTNO2  - REAL, optical depth due to no2 absorption at each            (O)
!=           specified altitude at each specified wavelength
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
!     ... dummy arguments
!-----------------------------------------------------------------------------
      REAL(rk), intent(in)    :: colno2(:)
      REAL(rk), intent(in)    :: no2_xs(:,:)
      REAL(rk), intent(inout) :: dtno2(:,:)

      INTEGER :: nlyr, nwave
      INTEGER :: wn

      nwave = size(no2_xs,dim=1)
      nlyr  = size(colno2)

      DO wn = 1,nwave
        dtno2(:nlyr,wn) = colno2(:nlyr) * no2_xs(wn,:nlyr)
      END DO

      END SUBROUTINE setno2

      subroutine setaer( nlambda_start, wc, tauaer300, tauaer400, &
                         tauaer600, tauaer999,               &
                         waer300, waer400, waer600, waer999, &
                         gaer300, gaer400, gaer600, gaer999, &
                         dtaer, omaer, gaer )
!----------------------------------------------------------------------
! The routine is based on aerosol treatment in module_ra_rrtmg_sw.F
! INPUT: 
! nzlev: number of specified altitude levels in the working grid
! z: specified altitude working grid   
! Aerosol optical properties at 300, 400, 600 and 999 nm. 
!   tauaer300, tauaer400, tauaer600, tauaer999: Layer AODs
!   waer300, waer400, waer600, waer999: Layer SSAs
!   gaer300, gaer400, gaer600, gaer999: Layer asymmetry parameters

! OUTPUT:
! dtaer: Layer AOD at FTUV wavelengths
! omaer: Layer SSA at FTUV wavelengths
! gaer : Layer asymmetry parameters at FTUV wavelengths
!------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! Dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nlambda_start
      real(rk), intent(in)  :: wc(:)
      real(rk), intent(in)  :: tauaer300(:), tauaer400(:),    &
                           tauaer600(:), tauaer999(:)
      real(rk), intent(in)  :: waer300(:), waer400(:),        &
                           waer600(:), waer999(:)
      real(rk), intent(in)  :: gaer300(:), gaer400(:),        &
                           gaer600(:), gaer999(:)
      real(rk), intent(out) :: dtaer(:,:), omaer(:,:), gaer(:,:)

!-----------------------------------------------------------------------------
! Local Variables
!-----------------------------------------------------------------------------
      real(rk), parameter :: thresh = 1.e-9_rk
      integer     :: k, wn, nlyr, nwave
      real(rk)        :: ang, slope, wfac

      nlyr =  size(dtaer,dim=1)
      nwave = size(dtaer,dim=2)

wave_loop: &
      do wn = nlambda_start,nwave
        wfac = wc(wn)*1.e-3_rk - .6_rk
        do k = 1,nlyr-1
!-----------------------------------------------------------------------------
! use angstrom exponent to calculate aerosol optical depth; wc is in nm.  
!-----------------------------------------------------------------------------
          if( tauaer300(k) > thresh .and. tauaer999(k) > thresh ) then
            ang = log(tauaer300(k)/tauaer999(k))/log(0.999_rk/0.3_rk)
            dtaer(k,wn) = tauaer400(k)*(0.4_rk/(wc(wn)*1.e-3_rk))**ang
!-----------------------------------------------------------------------------
! ssa - use linear interpolation/extrapolation
!-----------------------------------------------------------------------------
            slope = 5._rk*(waer600(k) - waer400(k))
            omaer(k,wn) = slope*wfac + waer600(k)
            omaer(k,wn) = max( .4_rk,min( 1._rk,omaer(k,wn) ) )
!-----------------------------------------------------------------------------
! asymmetry parameter - use linear interpolation/extrapolation
!-----------------------------------------------------------------------------
            slope = 5._rk*(gaer600(k) - gaer400(k))
            gaer(k,wn) = slope*wfac + gaer600(k)
            gaer(k,wn) = max( .5_rk,min( 1._rk,gaer(k,wn) ) )
          endif
        end do
      end do wave_loop
      
! CGB - For debug, set everything to the 400 nm values.
do k = 1,nlyr-1
  dtaer(k,:) = tauaer400(k)
  omaer(k,:) = waer400(k)
  gaer(k,:) = gaer400(k)
end do

      end subroutine setaer

      subroutine setcld( nlambda_start, cld_od_opt, z, xlwc, cldfrac, &
                         dtcld, omcld, gcld, errmsg, errflg )
!-----------------------------------------------------------------------------
!= PURPOSE:
!= Set up cloud optical depth, single albedo and g
!-----------------------------------------------------------------------------
!= PARAMETERS:
!= PARAMETERS:
!= NZ - INTEGER, number of specified altitude levels in the working (I)
!= grid
!= Z - real(dp), specified altitude working grid (km) (I)
!= XLWC Cloud water content g/M3 (I)
!=
!= dtcld - cloud optical depth
!= omcld - cloud droplet single albedo
!= gcld  - g
!-----------------------------------------------------------------------------
!
! VERTICAL DOMAIN is from bottom(1) to TOP (TOP=nz)
! CCM from top(1) to bottom(nz)
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------------
! ... Dummy arguments
!-----------------------------------------------------------------------------
      integer, intent(in) :: nlambda_start
      integer, intent(in) :: cld_od_opt
      real(rk), intent(in)  :: z(:)
      real(rk), intent(in)  :: xlwc(:)
      real(rk), intent(in)  :: cldfrac(:)
      real(rk), intent(inout) :: dtcld(:,:)
      real(rk), intent(inout) :: omcld(:,:)
      real(rk), intent(inout) :: gcld(:,:)

      character(len=*), intent(out)   :: errmsg
      integer,            intent(out)   :: errflg

!-----------------------------------------------------------------------------
! ... Local variables
!-----------------------------------------------------------------------------
      real(rk), parameter :: km2m = 1.e3_rk        ! kilometer to meter
      real(rk), parameter :: wden = 1.e6_rk        ! g/m3 (1 m3 water = 1e6 g water)
      real(rk), parameter :: re = 10.0_rk * 1.e-6_rk  ! assuming cloud drop radius = 10 um to M
      real(rk), parameter :: fac = 1._rk/(wden*re)

      integer  :: astat
      integer  :: wn
      integer  :: nlyr, nwave
      real(rk), allocatable :: wrk(:), layer_cldfrac(:)

      errmsg = ' '
      errflg = 0

      nlyr  = size(dtcld,dim=1)
      nwave = size(dtcld,dim=2)

      allocate( wrk(nlyr),layer_cldfrac(nlyr),stat=astat )
      if( astat /= 0 ) then
         errmsg = 'setcld: failed to allocate wrk'
         errflg = 1
         return
      endif

!-----------------------------------------------------------------------------
! ... calculate optical depth
!-----------------------------------------------------------------------------     
      wrk(1:nlyr-1) = (z(2:nlyr) - z(1:nlyr-1))*km2m   !  (km -> m)
      wrk(1:nlyr-1) = 1.5_rk * .5_rk*(xlwc(1:nlyr-1) + xlwc(2:nlyr))*wrk(1:nlyr-1)*fac
      wrk(1:nlyr-1) = max( wrk(1:nlyr-1),0._rk )
      if( cld_od_opt == 2 ) then
        layer_cldfrac(1:nlyr-1) = .5_rk*(cldfrac(1:nlyr-1) + cldfrac(2:nlyr))
        wrk(1:nlyr-1) = wrk(1:nlyr-1)*layer_cldfrac(1:nlyr-1)*sqrt( layer_cldfrac(1:nlyr-1) )
      endif
!----------------------------------------------------
! ....calculate cloud optical depth T
! following Liao et al. JGR, 104, 23697, 1999
!----------------------------------------------------
      if( any( wrk(1:nlyr-1) > 0._rk ) ) then
        do wn = nlambda_start,nwave
          dtcld(1:nlyr-1,wn) = wrk(1:nlyr-1)
          omcld(1:nlyr-1,wn) = .9999_rk
          gcld (1:nlyr-1,wn) = .85_rk
        end do
      endif

      if( allocated( wrk ) ) then
        deallocate( wrk )
      endif
      if( allocated( layer_cldfrac ) ) then
        deallocate( layer_cldfrac )
      endif

      end subroutine setcld
      
      
!-----------------------------------------------------------------------------*
!  PURPOSE:                                                                 =*
!  Calculate the action spectrum value for erythema at a given wavelength   =*
!  according to: McKinlay, A.F and B.L.Diffey, A reference action spectrum  =*
!  for ultraviolet induced erythema in human skin, CIE Journal, vol 6,      =*
!  pp 17-22, 1987.                                                          =*
!  Value at 300 nm = 0.6486                                                 =*
!----------------------------------------------------------------------------*
!  PARAMETERS:                                                              =*
!  W - REAL, wavelength (nm)                                             (I)=*
!----------------------------------------------------------------------------*

      function fery(w)


      implicit none

! input:
      real(rk) w 

! function value:
      real(rk) fery

      if (w .lt. 250._rk) then
          fery = 1._rk
! outside the ery spectrum range
      elseif ((w .ge. 250._rk) .and. (w .lt. 298._rk)) then
          fery = 1._rk
      elseif ((w .ge. 298._rk) .and. (w .lt. 328._rk)) then
          fery = 10._rk**( 0.094_rk*(298._rk - w) )
      elseif ((w .ge. 328._rk) .and. (w .lt. 400._rk)) then
          fery = 10._rk**( 0.015_rk*(139._rk - w) )
      else
         fery = 1.e-36_rk
! outside the ery spectrum range
      endif

      return
      end
      

!-----------------------------------------------------------------------------*
!  PURPOSE:                                                                 =*
!  Calculate the action spectrum value for skin cancer of albino hairless   =*
!  mice at a given wavelength according to:  deGRuijl, F.R., H.J.C.M.Steren-=*
!  borg, P.D.Forbes, R.E.Davies, C.Colse, G.Kelfkens, H.vanWeelden,         =*
!  and J.C.van der Leun, Wavelength dependence of skin cancer induction by  =*
!  ultraviolet irradiation of albino hairless mice, Cancer Research, vol 53,=*
!  pp. 53-60, 1993                                                          =*
!  (Action spectrum for carcinomas)                                         =*
!-----------------------------------------------------------------------------*
!  PARAMETERS:                                                              =*
!  W  - REAL, wavelength (nm)                                            (I)=*
!-----------------------------------------------------------------------------*

      function futr(w)

      implicit none

! input:
      real(rk) w

! function value:
      real(rk) futr

! local:
      real(rk) :: a1, a2, a3, a4, a5
      real(rk) :: x1, x2, x3, x4, x5
      real(rk) :: t1, t2, t3, t4, t5
      real(rk) :: b1, b2, b3, b4, b5
      real(rk) :: p

      a1 = -10.91_rk
      a2 = - 0.86_rk
      a3 = - 8.60_rk
      a4 = - 9.36_rk
      a5 = -13.15_rk

      x1 = 270._rk
      x2 = 302._rk
      x3 = 334._rk
      x4 = 367._rk
      x5 = 400._rk

      t1 = (w-x2)*(w-x3)*(w-x4)*(w-x5)
      t2 = (w-x1)*(w-x3)*(w-x4)*(w-x5)
      t3 = (w-x1)*(w-x2)*(w-x4)*(w-x5)
      t4 = (w-x1)*(w-x2)*(w-x3)*(w-x5)
      t5 = (w-x1)*(w-x2)*(w-x3)*(w-x4)

      b1 = (x1-x2)*(x1-x3)*(x1-x4)*(x1-x5)
      b2 = (x2-x1)*(x2-x3)*(x2-x4)*(x2-x5)
      b3 = (x3-x1)*(x3-x2)*(x3-x4)*(x3-x5)
      b4 = (x4-x1)*(x4-x2)*(x4-x3)*(x4-x5)
      b5 = (x5-x1)*(x5-x2)*(x5-x3)*(x5-x4)

      p = a1*t1/b1 + a2*t2/b2 + a3*t3/b3 + a4*t4/b4 + a5*t5/b5

      futr  = exp(p)

      return
      end
      
      

!-----------------------------------------------------------------------------*
!  PURPOSE:                                                                 =*
!  Map input data given on single, discrete points onto a set of target     =*
!  bins.                                                                    =*
!  The original input data are given on single, discrete points of an       =*
!  arbitrary grid and are being linearly interpolated onto a specified set  =*
!  of target bins.  In general, this is the case for most of the weighting  =*
!  functions (action spectra, molecular cross section, and quantum yield    =*
!  data), which have to be matched onto the specified wavelength intervals. =*
!  The average value in each target bin is found by averaging the trapezoi- =*
!  dal area underneath the input data curve (constructed by linearly connec-=*
!  ting the discrete input values).                                         =*
!  Some caution should be used near the endpoints of the grids.  If the     =*
!  input data set does not span the range of the target grid, an error      =*
!  message is printed and the execution is stopped, as extrapolation of the =*
!  data is not permitted.                                                   =*
!  If the input data does not encompass the target grid, use ADDPNT to      =*
!  expand the input array.                                                  =*
!-----------------------------------------------------------------------------*
!  PARAMETERS:                                                              =*
!  NG  - INTEGER, number of bins + 1 in the target grid                  (I)=*
!  XG  - REAL, target grid (e.g., wavelength grid);  bin i is defined    (I)=*
!        as [XG(i),XG(i+1)] (i = 1..NG-1)                                   =*
!  YG  - REAL, y-data re-gridded onto XG, YG(i) specifies the value for  (O)=*
!        bin i (i = 1..NG-1)                                                =*
!  N   - INTEGER, number of points in input grid                         (I)=*
!  X   - REAL, grid on which input data are defined                      (I)=*
!  Y   - REAL, input y-data                                              (I)=*
!-----------------------------------------------------------------------------*

      subroutine inter2(ng,xg,yg,n,x,y,ierr)

      implicit none

      ! input:
      integer ng, n
      real(rk) x(n), y(n), xg(ng)

      ! output:
      real(rk) yg(ng)

      ! local:
      real(rk) area, xgl, xgu
      real(rk) darea, slope
      real(rk) a1, a2, b1, b2
      integer ngintv
      integer i, k, jstart
      integer ierr

      ierr = 0

      ! test for correct ordering of data, by increasing value of x
      do 10, i = 2, n
         if (x(i) .le. x(i-1)) then
            ierr = 1
            write(*,*) 'TUV::inter2 - ERROR data not sorted'
            return
         end if
   10 continue     

      do i = 2, ng
        if (xg(i) .le. xg(i-1)) then
           ierr = 2
           write(*,*) 'TUV::inter2 - ERROR  xg-grid not sorted!'
           return
        end if
      end do

      ! check for xg-values outside the x-range
      if ( (x(1) .gt. xg(1)) .or. (x(n) .lt. xg(ng)) ) then
          ierr = 3
          write(*,*) ('TUV::inter2 - ERROR data do not span grid. Use addpnt to expand data and re-run.')
          return
      end if

      ! find the integral of each grid interval and use this to 
      ! calculate the average y value for the interval      
      ! xgl and xgu are the lower and upper limits of the grid interval
      jstart = 1
      ngintv = ng - 1
      do 50, i = 1,ngintv

         ! initalize:
         area = 0.0_rk
         xgl = xg(i)
         xgu = xg(i+1)

         ! discard data before the first grid interval and after the 
         ! last grid interval
         ! for internal grid intervals, start calculating area by interpolating
         ! between the last point which lies in the previous interval and the
         ! first point inside the current interval
         k = jstart
         if (k .le. n-1) then

         ! if both points are before the first grid, go to the next point
   30       continue
            if (x(k+1) .le. xgl) then
               jstart = k - 1
               k = k+1
               if (k .le. n-1) go to 30
            end if


            ! if the last point is beyond the end of the grid, complete and go to the next
            ! grid
   40       continue
            if ((k .le. n-1) .and. (x(k) .lt. xgu)) then          

               jstart = k-1

               ! compute x-coordinates of increment
               a1 = max(x(k),xgl)
               a2 = min(x(k+1),xgu)

               ! if points coincide, contribution is zero
               if (x(k+1).eq.x(k)) then
                  darea = 0.e0_rk
               else
                  slope = (y(k+1) - y(k))/(x(k+1) - x(k))
                  b1 = y(k) + slope*(a1 - x(k))
                  b2 = y(k) + slope*(a2 - x(k))
                  darea = (a2 - a1)*(b2 + b1)/2._rk
               end if

               ! find the area under the trapezoid from a1 to a2
               area = area + darea

               ! go to next point
               k = k+1
               go to 40

             end if
          end if

          ! Calculate the average y after summing the areas in the interval
          yg(i) = area/(xgu - xgl)
   50 continue

      return
      end


!-----------------------------------------------------------------------------*
!  PURPOSE:                                                                 =*
!  Add a point <xnew,ynew> to a set of data pairs <x,y>.  x must be in      =*
!  ascending order                                                          =*
!-----------------------------------------------------------------------------*
!  PARAMETERS:                                                              =*
!  X    - REAL vector of length LD, x-coordinates                       (IO)=*
!  Y    - REAL vector of length LD, y-values                            (IO)=*
!  LD   - INTEGER, dimension of X, Y exactly as declared in the calling  (I)=*
!         program                                                           =*
!  N    - INTEGER, number of elements in X, Y.  On entry, it must be:   (IO)=*
!         N < LD.  On exit, N is incremented by 1.                          =*
!  XNEW - REAL, x-coordinate at which point is to be added               (I)=*
!  YNEW - REAL, y-value of point to be added                             (I)=*
!-----------------------------------------------------------------------------*

    subroutine addpnt( x, y, ld, n, xnew, ynew, ierr)
      implicit none

      ! calling parameters
      integer ld, n, ierr
      real(rk) x(ld), y(ld)
      real(rk) xnew, ynew

      ! local variables
      integer insert
      integer i
      
      ierr = 0

      ! check n<ld to make sure x will hold another point

      if (n .ge. ld) then
         ierr = 1
         write(*,*) 'TUV::addpnt - ERROR Cannot expand array all elements used.'
      endif

      insert = 1
      i = 2

      ! check, whether x is already sorted.
      ! also, use this loop to find the point at which xnew needs to be inserted
      ! into vector x, if x is sorted.

 10   continue
      if (i .lt. n) then
        if (x(i) .lt. x(i-1)) then
           ierr = 2
           write(*,*) 'TUV_addpnt - ERROR x-data must be in ascending order!'
        else
           if (xnew .gt. x(i)) insert = i + 1
        endif
        i = i+1
        goto 10
      endif

      ! if <xnew,ynew> needs to be appended at the end, just do so,
      !otherwise, insert <xnew,ynew> at position insert
      if ( xnew .gt. x(n) ) then
         x(n+1) = xnew
         y(n+1) = ynew
      else

         ! shift all existing points one index up
         do i = n, insert, -1
           x(i+1) = x(i)
           y(i+1) = y(i)
         enddo

         ! insert new point
         x(insert) = xnew
         y(insert) = ynew
      endif

      ! increase total number of elements in x, y
      n = n+1
    end


  END MODULE tuv_subs
