program driver

  use machine, only: r8 => kind_phys
  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use params_mod, only: input_data_root
  use environ_conditions_mod, only: environ_conditions_create, environ_conditions
  use rad_abs_xsect,    only: nwave
  use tuv_photolysis,   only: tuv_photolysis_readnl, tuv_photolysis_init, tuv_photolysis_run
  use tuv_photolysis,   only: tuv_n_phot
  use tuv_photolysis,   only: tuv_jnames
  use module_prates_tuv,only: rxn_ndx
  use module_rxn, only : xsqy_table => xsqy_tab
  use molec_ox_xsect, only: molec_ox_xsect_init
  use molec_ox_xsect, only: molec_ox_xsect_run

  implicit none

  character(len=50) :: rxn_string
  integer :: i,jndx
  integer :: errflg
  character(len=444) :: errmsg

  real(r8) :: zenith
  real(r8) :: albedo
  real(r8), allocatable :: alt(:)
  real(r8), allocatable :: press_mid(:)
  real(r8), allocatable :: press_int(:)
  real(r8), allocatable :: temp(:)
  real(r8), allocatable :: o2vmrcol(:)
  real(r8), allocatable :: o3vmrcol(:)
  real(r8), allocatable :: so2vmrcol(:)
  real(r8), allocatable :: no2vmrcol(:)
  real(r8), allocatable :: srb_o2_xs(:,:)
  real(r8), allocatable :: dto2(:,:)
  real(r8), allocatable :: radfld(:,:)
  real(r8), allocatable :: tuv_prates(:,:)
  
  character(len=*), parameter :: nml_file = '/terminator-data1/home/fvitt/tuvdev/driver/test_nml'
  character(len=*), parameter :: env_conds_file = &
       '/terminator-data1/fvitt/micm_inputs/FW2000climo.f09_f09_mg17.cam6_0_030.n01.cam.h2.0001-01-01-00000.nc'
  type(environ_conditions),pointer :: colEnvConds => null()

  integer :: nlevels,k

  write(*,*) 'BEGIN TEST'

  colEnvConds => environ_conditions_create( env_conds_file, lat=45.0, lon=180. )

  nlevels = colEnvConds%nlevels()

  call tuv_photolysis_readnl(nml_file)
  
  errflg=0
  errmsg=' '

  call molec_ox_xsect_init( errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_photolysis_init( r8, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if
  
  call tuv_radiation_transfer_init( r8, nlevels, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  allocate(srb_o2_xs(nwave,nlevels), dto2(nlevels-1,nwave))
  allocate(radfld(nwave,nlevels))
  allocate(tuv_prates(nlevels, tuv_n_phot ) )
 
  allocate(alt(nlevels))
  allocate(press_mid(nlevels))
  allocate(press_int(nlevels))
  allocate(temp(nlevels))
  allocate(o2vmrcol(nlevels))
  allocate(o3vmrcol(nlevels))
  allocate(so2vmrcol(nlevels))
  allocate(no2vmrcol(nlevels))

  zenith = colEnvConds%getsrf('SZA')
  albedo = colEnvConds%getsrf('ASDIR')
  press_mid(:nlevels) = colEnvConds%press_mid(nlevels)
  alt(:nlevels) = colEnvConds%getcol('Z3',nlevels)*1.e-3 ! m --> km
  temp(:nlevels) = colEnvConds%getcol('T',nlevels)
  o2vmrcol(:nlevels) = colEnvConds%getcol('O2',nlevels)
  o3vmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  so2vmrcol(:nlevels) = colEnvConds%getcol('SO2',nlevels)
  no2vmrcol(:nlevels) = colEnvConds%getcol('NO2',nlevels)

  call molec_ox_xsect_run( nlevels, zenith, alt, temp, press_mid, o2vmrcol, dto2, srb_o2_xs )

  call  tuv_radiation_transfer_run( &
       zenith, albedo, press_mid, alt, temp, o3vmrcol, so2vmrcol, no2vmrcol, dto2, radfld, errmsg, errflg )

  call tuv_photolysis_run( nlevels, temp, press_mid, radfld, srb_o2_xs, tuv_prates, errmsg, errflg )
  
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  do i=1,tuv_n_phot

     jndx = rxn_ndx(i)
     if (jndx>0) then
        rxn_string = trim(xsqy_table(jndx)%label)
     else
        rxn_string = 'O2 -> O + O'
     end if
     write( *,*) trim(tuv_jnames(i))//'   '//trim(rxn_string)
     write(10,*) trim(tuv_jnames(i))//'   '//trim(rxn_string)
     write( *,'("  rate = ",e12.4," /sec")' ) tuv_prates(nlevels,i)
     do k=1,nlevels
         write(10,'("  rate = ",e24.16," /sec")' ) tuv_prates(k,i)
     end do

     write(10,*) ' '
     print*,' '
  end do

  write(*,*) 'END TEST'
  
end program driver
