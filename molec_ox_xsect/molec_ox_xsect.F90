module molec_ox_xsect
  use phot_kind_mod, only: rk => kind_phot
  use rad_abs_xsect, only: o2_xs, nwave, wl

  implicit none

contains
  subroutine molec_ox_xsect_init( errmsg, errflg )
    use la_srb_mod, only : la_srb_init
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    call la_srb_init( errmsg, errflg )

  end subroutine molec_ox_xsect_init

  subroutine molec_ox_xsect_run( nlev, zen, alt, temp, press_mid, o2vmr, dto2, srb_o2_xs )
    use rad_abs_xsect, only : o2_xs
    use phot_util_mod, only : sphers, airmas
    use la_srb_mod,    only : la_srb_comp
    use params_mod, only: R, g, kboltz

    integer,  intent(in)  :: nlev
    real(rk), intent(in)  :: zen
    real(rk), intent(in)  :: alt(:)  ! km
    real(rk), intent(in)  :: temp(:) ! K
    real(rk), intent(in)  :: press_mid(:)
    real(rk), intent(in)  :: o2vmr(:)
    real(rk), intent(inout) :: dto2(:,:)
    real(rk), intent(inout) :: srb_o2_xs(:,:)

    integer :: nlyr, k
    integer :: wn
    integer  :: nid(0:nlev-1)
    real(rk) :: dsdh(0:nlev-1,nlev-1)
    real(rk) :: vcol(nlev-1)
    real(rk) :: scol(nlev-1)
    real(rk) :: tlev(nlev)
    real(rk) :: zlev(nlev)

    real(rk) :: aircol(nlev-1)  ! # molecules / cm2 in each layer
    real(rk) :: o2col(nlev)  
    real(rk) :: dpress(nlev)

    dto2(:,:) = 0._rk
    srb_o2_xs(:,:) = 0._rk

    nlyr = nlev-1
    dpress(1:nlyr) = press_mid(2:nlyr+1) - press_mid(1:nlyr)
    do k=1,nlyr
       aircol(k) = 10._rk*dpress(k)*R/(kboltz*g)
       o2col(k)  = 0.5_rk*(o2vmr(k)+o2vmr(k+1))*aircol(k)
    end do
    aircol(1:nlyr) = aircol(nlyr:1:-1)
    o2col(1:nlyr)  = o2col(nlyr:1:-1)
    tlev(nlev:1:-1) = temp(1:nlev)
    zlev(nlev:1:-1) = alt(1:nlev) ! km

    do wn = 1,nwave
       dto2(:nlyr,wn) = o2col(:nlyr) * o2_xs(wn)
    end do

    call sphers( nlyr, zlev, zen, dsdh, nid )
    call airmas( nlyr, dsdh, nid, aircol, vcol, scol )

    call la_srb_comp( nlyr, wl(1), tlev, vcol, scol, o2_xs, dto2, srb_o2_xs )

  end subroutine molec_ox_xsect_run

end module molec_ox_xsect
