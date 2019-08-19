module molec_ox_xsect
  use phot_kind_mod, only: rk => kind_phot
  use wavelength_grid, only: nwave, wl

  implicit none

contains

!> \section arg_table_molec_ox_xsect_init Argument Table
!! | local_name | standard_name             | long_name                 | units   | rank | type      | kind      | intent | optional |
!! |------------|---------------------------|---------------------------|---------|------|-----------|-----------|--------|----------|
!! | errmsg     | ccpp_error_message        | CCPP error message        | none    |    0 | character | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag           | flag    |    0 | integer   |           | out    | F        |
!!
  subroutine molec_ox_xsect_init( errmsg, errflg )
    use la_srb_mod, only : la_srb_init

    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call la_srb_init( errmsg, errflg )

  end subroutine molec_ox_xsect_init

!> \section arg_table_molec_ox_xsect_run Argument Table
!! | local_name | standard_name             | long_name                          | units     | rank | type        | kind      | intent | optional |
!! |------------|---------------------------|------------------------------------|-----------|------|-------------|-----------|--------|----------|
!! | nlev       | num_levels_for_photolysis | number of column layers            | count     |    0 | integer     |           | in     | F        |
!! | zen        | solar_zenith              | solar zenith angle                 | degrees   |    0 | real        | kind_phys | in     | F        |
!! | alt        | layer_altitude            | mid-point layer altitude           | m         |    1 | real        | kind_phys | in     | F        |
!! | temp       | layer_temperature         | mid-point layer temperature        | K         |    1 | real        | kind_phys | in     | F        |
!! | press_mid  | layer_pressure            | mid-point layer pressure           | Pa        |    1 | real        | kind_phys | in     | F        |
!! | o2vmr      | O2_vmr_col                | O2 volume mixing ratio column      | mole/mole |    1 | real        | kind_phys | in     | F        |
!! | dto2       | O2_optical_depth          | optical depth due to O2 absorption | cm        |    2 | real        | kind_phys | in     | F        |
!! | srb_o2_xs  | O2_xsect                  | O2 effective cross section         | cm2       |    2 | real        | kind_phys | in     | F        |
!! | errmsg     | ccpp_error_message        | CCPP error message                 | none      |    0 | character   | len=*     | out    | F        |
!! | errflg     | ccpp_error_flag           | CCPP error flag                    | flag      |    0 | integer     |           | out    | F        |
!!
  subroutine molec_ox_xsect_run( nlev, zen, alt, temp, press_mid, press_top, o2vmr, dto2, srb_o2_xs, errmsg, errflg )
    use module_xsections, only: o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod,    only : R, g, kboltz

    integer,          intent(in)    :: nlev
    real(rk),         intent(in)    :: zen
    real(rk),         intent(in)    :: alt(:)  ! m
    real(rk),         intent(in)    :: temp(:) ! K
    real(rk),         intent(in)    :: press_mid(:)
    real(rk),         intent(in)    :: press_top
    real(rk),         intent(in)    :: o2vmr(:)
    real(rk),         intent(inout) :: dto2(:,:)
    real(rk),         intent(inout) :: srb_o2_xs(:,:)
    character(len=*), intent(out)   :: errmsg
    integer,          intent(out)   :: errflg

    integer  :: nlyr, k
    integer  :: wn
    integer  :: nid(0:nlev)
    real(rk) :: dsdh(0:nlev,nlev)
    real(rk) :: vcol(nlev)
    real(rk) :: scol(nlev)
    real(rk) :: tlev(nlev)
    real(rk) :: zlev(nlev+1)
    real(rk) :: o2lev(nlev+1)

    real(rk) :: aircol(nlev)  ! # molecules / cm2 in each layer
    real(rk) :: o2col(nlev)  
    real(rk) :: dpress(nlev-1)
    real(rk) :: delz_cm, delz_km

    errmsg = ''
    errflg = 0

    dto2(:,:) = 0._rk
    srb_o2_xs(:,:) = 0._rk

    nlyr = nlev-1

    dpress(nlyr:1:-1) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    zlev(nlev:1:-1) = alt(1:nlev) *1.e-3_rk ! m -> km
    o2lev(nlev:1:-1) = o2vmr(1:nlev)

    delz_km = zlev(nlev) - zlev(nlev-1)
    delz_cm = delz_km*1.e5_rk ! cm

    zlev(nlev+1) = zlev(nlev) + delz_km ! km

    do k=1,nlyr
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_rk*(o2lev(k)+o2lev(k+1))*aircol(k)
    end do

    aircol(nlev) = delz_cm * 10._rk * press_top / ( kboltz * temp(1) ) ! molecules / cm2
    o2col(nlev)  = o2lev(nlev) * aircol(nlev)

    tlev(nlev:1:-1) = temp(1:nlev)

    do wn = 1,nwave
       dto2(:nlev,wn) = o2col(:nlev) * o2_xs(wn)
    end do

    call sphers( nlev, zlev, zen, dsdh, nid )
    call airmas( nlev, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlev, wl(1), tlev, vcol, scol, o2lev, o2_xs, dto2, srb_o2_xs )

  end subroutine molec_ox_xsect_run

!> \section arg_table_molec_ox_xsect_finalize Argument Table
!! | local_name | standard_name      | long_name          | units     | rank | type      | kind  | intent | optional |
!! |------------|--------------------|--------------------|-----------|------|-----------|-------|--------|----------|
!! | errmsg     | ccpp_error_message | CCPP error message | none      |    0 | character | len=* | out    | F        |
!! | errflg     | ccpp_error_flag    | CCPP error flag    | flag      |    0 | integer   |       | out    | F        |
!!
  subroutine molec_ox_xsect_finalize( errmsg, errflg )

    !--- arguments
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    !--- initialize CCPP error handling variables
    errmsg = ''
    errflg = 0

  end subroutine molec_ox_xsect_finalize

end module molec_ox_xsect
