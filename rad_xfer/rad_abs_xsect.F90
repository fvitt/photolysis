module rad_abs_xsect

  implicit none

  public

  real, protected, allocatable :: o2_xs(:)
  real, protected, allocatable :: so2_xs(:)
  real, protected, allocatable :: wl(:)
  real, protected, allocatable :: wc(:)

  integer, protected :: nwave

contains

  subroutine rad_abs_xsect_init( filepath, errmsg, errflg )

    use netcdf
    use la_srb_mod, only : ila, isrb
    use la_srb_mod, only : nchebev_term, nchebev_wave
    use la_srb_mod, only : chebev_ac, chebev_bc
    use la_srb_mod, only : init_srb

    character(len=*), intent(in ) :: filepath
    character(len=*), intent(out) :: errmsg
    integer,          intent(out) :: errflg

    integer :: ncid, dimid, varid
    integer :: astat, ret

    ! open file
    ret = nf90_open( trim(filepath), nf90_noclobber, ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to open '//trim(filepath)
       return
    end if

    ! get dimensions
    ret = nf90_inq_dimid( ncid, 'nwave', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get nwave id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nwave )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get nwave'
       return
    end if

    ! allocate memory
    allocate( wl(nwave+1), wc(nwave), o2_xs(nwave), so2_xs(nwave), stat=astat )
    if( astat /= 0 ) then
       errflg = astat
       errmsg = 'rad_abs_xsect_init: failed to allocate memory'
       return
    end if

    ! read xsect data
    ret = nf90_inq_varid( ncid, 'o2_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get o2_xs variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, o2_xs )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read o2_xs variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'so2_xs', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get so2_xs variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, so2_xs )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read so2_xs variable'
       return
    end if

    ret = nf90_inq_varid( ncid, 'wl', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get wl variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wl )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read wl variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'wc', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to get wc variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, wc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to read wc variable'
       return
    end if

    ! for la_srb
    ret = nf90_inq_dimid( ncid, 'nchebev_term', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nchebev_term id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nchebev_term )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nchebev'
       return
    end if
    ret = nf90_inq_dimid( ncid, 'nchebev_wave', dimid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nchebev_wave id'
       return
    end if
    ret = nf90_inquire_dimension( ncid, dimid, len=nchebev_wave )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get nchebev'
       return
    end if
    allocate( chebev_ac(nchebev_term,nchebev_wave), chebev_bc(nchebev_term,nchebev_wave), stat=astat )
    if( astat /= 0 ) then
       errflg = astat
       errmsg = 'rad_abs_xsect_init: failed to allocate la_srb memory'
       return
    end if
    ret = nf90_inq_varid( ncid, 'chebev_ac', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get chebev_ac variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, chebev_ac )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read chebev_ac variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'chebev_bc', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get chebev_bc variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, chebev_bc )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read chebev_bc variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'ila', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get ila variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, ila )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read ila variable'
       return
    end if
    ret = nf90_inq_varid( ncid, 'isrb', varid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to get isrb variable id'
       return
    end if
    ret = nf90_get_var( ncid, varid, isrb )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'get_xsqy_tab: failed to read isrb variable'
       return
    end if

    ! close the file
    ret = nf90_close( ncid )
    if( ret /= nf90_noerr ) then
       errflg = 1
       errmsg = 'rad_abs_xsect_init: failed to close '//trim(filepath)
       return
    end if

    call init_srb()

  end subroutine rad_abs_xsect_init

end module rad_abs_xsect
