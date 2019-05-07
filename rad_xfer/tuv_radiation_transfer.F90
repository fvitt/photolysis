module tuv_radiation_transfer
  use machine, only: rk => kind_phys

  implicit none

  private
  public :: tuv_radiation_transfer_init
  public :: tuv_radiation_transfer_run
  
  integer :: nlev, nlyr
  
contains
  
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  subroutine tuv_radiation_transfer_init( nlev_in, errmsg, errflg )
    use params_mod, only: input_data_root
    use rad_abs_xsect, only: rad_abs_xsect_init
    use module_xsections, only: rdxs_init
    use rad_abs_xsect,    only: nwave, wl, wc

    integer,          intent(in)  :: nlev_in
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    character(len=256) :: filepath
    
    errmsg = ''
    errflg = 0

    nlev = nlev_in
    nlyr = nlev-1

    filepath = trim(input_data_root)//'/wrf_tuv_xsqy.nc'

    call rad_abs_xsect_init( filepath, errmsg, errflg )
    if (errflg.ne.0) return
    
    call rdxs_init( nwave, wl, errmsg, errflg )
    if (errflg.ne.0) return
   
  end subroutine tuv_radiation_transfer_init
  
  !-------------------------------------------------------------------------------------------------
  !-------------------------------------------------------------------------------------------------
  subroutine tuv_radiation_transfer_run( &
       zenith, albedo, press_mid, alt, temp, o2vmr, o3vmr, so2vmr, no2vmr, radfld, srb_o2_xs, errmsg, errflg )

    use tuv_subs,         only: tuv_radfld
    use rad_abs_xsect,    only: o2_xs, so2_xs, nwave, wl, wc
    use module_xsections, only: o3xs, no2xs_jpl06a

    real(rk),         intent(in)  :: zenith
    real(rk),         intent(in)  :: albedo
    real(rk),         intent(in)  :: press_mid(:)
    real(rk),         intent(in)  :: alt(:)  ! km
    real(rk),         intent(in)  :: temp(:) ! K
    real(rk),         intent(in)  :: o2vmr(:)
    real(rk),         intent(in)  :: o3vmr(:)
    real(rk),         intent(in)  :: so2vmr(:)
    real(rk),         intent(in)  :: no2vmr(:)
    real,             intent(out) :: radfld(:,:) ! /sec
    real,             intent(out) :: srb_o2_xs(:,:)
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer, parameter :: nlambda_start=1
    integer, parameter :: cld_od_opt=1
    logical, parameter :: has_aer_ra_feedback = .false.
    real,    parameter :: dobsi = 0.
    
    real, parameter :: kboltz= 1.38064852e-16 ! boltzmann constant (erg/K)
    real, parameter :: R=2.8704e6       ! gas constant (erg/g/K)
    real, parameter :: g=980.616        ! grav acceleration (cm/sec2)

    real :: zen
    real :: alb(nwave)
    real :: zlev(nlev) ! km 
    real :: tlev(nlev)
    real :: aircol(nlyr)  ! # molecules / cm2 in each layer
    real :: o2col(nlyr)  
    real :: o3col(nlyr) 
    real :: so2col(nlyr)
    real :: no2col(nlyr)
    real :: dpress(nlyr)

    real :: tauaer300(nlev) ! aerosol properties
    real :: tauaer400(nlev)
    real :: tauaer600(nlev)
    real :: tauaer999(nlev)
    real :: waer300(nlev)
    real :: waer400(nlev)
    real :: waer600(nlev)
    real :: waer999(nlev)
    real :: gaer300(nlev)
    real :: gaer400(nlev)
    real :: gaer600(nlev)
    real :: gaer999(nlev)
    
    real :: dtaer(nlyr,nwave), omaer(nlyr,nwave), gaer(nlyr,nwave)
    real :: dtcld(nlyr,nwave), omcld(nlyr,nwave), gcld(nlyr,nwave)
    real :: dt_cld(nlyr)
    
    real :: qll(nlev) ! cld water content (g/m3)
    real :: cldfrac(nlev)
    real :: efld(nlev,nwave)
    real :: e_dir(nlev,nwave)
    real :: e_dn(nlev,nwave)
    real :: e_up(nlev,nwave)
    real :: dir_fld(nlev,nwave)
    real :: dwn_fld(nlev,nwave)
    real :: up_fld(nlev,nwave)

    integer :: k, kk

    real :: o3_xs(nwave,nlev)
    real :: no2_xs(nwave,nlev)
    real :: o3_xs_tpose(nlev,nwave)
    real :: no2_xs_tpose(nlev,nwave)

    errmsg = ''
    errflg = 0

    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       kk=nlyr-k+1
       aircol(k) = 10.*dpress(k)*R/(kboltz*g)
       o3col(kk)  = 0.5*(o3vmr(k)+o3vmr(k+1))*aircol(k)
       o2col(kk)  = 0.5*(o2vmr(k)+o2vmr(k+1))*aircol(k)
       so2col(kk) = 0.5*(so2vmr(k)+so2vmr(k+1))*aircol(k)
       no2col(kk) = 0.5*(no2vmr(k)+no2vmr(k+1))*aircol(k)
    end do

    ! inputs need to be bottom up vert coord
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev) ! km

    qll=0.0
    cldfrac=0.0
    tauaer300=0.0
    tauaer400=0.0
    tauaer600=0.0
    tauaer999=0.0
    waer300=1.0
    waer400=1.0
    waer600=1.0
    waer999=1.0
    gaer300=0.0
    gaer400=0.0
    gaer600=0.0
    gaer999=0.0

    zen = zenith
    alb(:) = albedo

    call o3xs( nlev,tlev,nwave,wl,o3_xs_tpose )
    call no2xs_jpl06a( nlev,tlev,nwave,wl,no2_xs_tpose )
    o3_xs  = transpose( o3_xs_tpose )
    no2_xs = transpose( no2_xs_tpose )

    call tuv_radfld( nlambda_start, cld_od_opt, cldfrac, nlyr, nwave, &
         zen, zlev, alb, &
         aircol, o2col, o3col, so2col, no2col, &
         tauaer300, tauaer400, tauaer600, tauaer999, &
         waer300, waer400, waer600, waer999, &
         gaer300, gaer400, gaer600, gaer999, &
         dtaer, omaer, gaer, dtcld, omcld, gcld, &
         has_aer_ra_feedback, &
         qll, dobsi, o3_xs, no2_xs, o2_xs, &
         so2_xs, wl(1), wc, tlev, srb_o2_xs, radfld, efld, &
         e_dir, e_dn, e_up, &
         dir_fld, dwn_fld, up_fld, dt_cld, errmsg, errflg )

  end subroutine tuv_radiation_transfer_run
  
end module tuv_radiation_transfer
