program driver

  use machine, only: r8 => kind_phys
  use tuv_radiation_transfer, only: tuv_radiation_transfer_init
  use tuv_radiation_transfer, only: tuv_radiation_transfer_run
  use params_mod, only: input_data_root
  use environ_conditions_mod, only: environ_conditions_create, environ_conditions
  use tuv_photolysis,   only: tuv_photolysis_readnl, tuv_photolysis_init, tuv_photolysis_run
  use tuv_photolysis,   only: tuv_n_phot, tuv_n_wavelen
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
  real(r8), allocatable :: cldfrac(:)
  real(r8), allocatable :: cldwat(:)
  real(r8), allocatable :: srb_o2_xs(:,:)
  real(r8), allocatable :: dto2(:,:)
  real(r8), allocatable :: radfld(:,:)
  real(r8), allocatable :: tuv_prates(:,:)

  character(len=*), parameter :: nml_file = 'test_nml'
  character(len=*), parameter :: env_conds_file = &
       '/terminator-data1/fvitt/micm_inputs/MusicBox_env_cond_1col_c190109.nc'
  type(environ_conditions),pointer :: colEnvConds => null()

  integer :: nlevels,k

  character(len=16), parameter :: my_jnames(119) = &
       (/'j_o2            ' &
       , 'jo3_a           ' &
       , 'jo3_b           ' &
       , 'jno2            ' &
       , 'jno3_a          ' &
       , 'jno3_b          ' &
       , 'jn2o5_a         ' &
       , 'jn2o5_b         ' &
       , 'jhno2           ' &
       , 'jhno3           ' &
       , 'jhno4           ' &
       , 'jh2o2           ' &
       , 'jchbr3          ' &
       , 'jch3cho         ' &
       , 'j_ch3cho_a      ' &
       , 'j_ch3cho_b      ' &
       , 'j_ch3cho_c      ' &
       , 'j_c2h5cho       ' &
       , 'j_gly_a         ' &
       , 'j_gly_b         ' &
       , 'j_gly_c         ' &
       , 'jmgly           ' &
       , 'j_ch3coch3      ' &
       , 'jacet           ' &
       , 'j_ch3ooh        ' &
       , 'j_ch3ono2       ' &
       , 'jch3ooh         ' &
       , 'jpan            ' &
       , 'j_pan_a         ' &
       , 'j_pan_b         ' &
       , 'j_ccl2o         ' &
       , 'jccl4           ' &
       , 'j_cclfo         ' &
       , 'j_cf2o          ' &
       , 'jcfc113         ' &
       , 'jcfc114         ' &
       , 'jcfc115         ' &
       , 'jcfcl3          ' &
       , 'jcf2cl2         ' &
       , 'jch3br          ' &
       , 'jch3ccl3        ' &
       , 'jch3cl          ' &
       , 'j_cloo          ' &
       , 'j_cf3chcl2      ' &
       , 'j_cf3chfcl      ' &
       , 'jhcfc141b       ' &
       , 'jhcfc142b       ' &
       , 'j_cf3cf2chcl2   ' &
       , 'j_cf2clcf2chfcl ' &
       , 'jhcfc22         ' &
       , 'j_ho2           ' &
       , 'j_cf2bf2        ' &
       , 'j_cf2brcl       ' &
       , 'jcf3br          ' &
       , 'jh2402          ' &
       , 'jn2o            ' &
       , 'jclono2_a       ' &
       , 'jclono2_b       ' &
       , 'jbrono2_b       ' &
       , 'jbrono2_a       ' &
       , 'jcl2            ' &
       , 'j_glyald_a      ' &
       , 'j_glyald_b      ' &
       , 'j_glyald_c      ' &
       , 'j_biacetyl      ' &
       , 'jmvk            ' &
       , 'j_macr          ' &
       , 'j_ch3cocooh     ' &
       , 'j_ch3ch2ono2    ' &
       , 'j_ch3chono2ch3  ' &
       , 'j_ch2ohch2ono2  ' &
       , 'j_ch3coch2ono2  ' &
       , 'j_bnit1         ' &
       , 'j_bnit2         ' &
       , 'j_cloocl        ' &
       , 'j_hyac_a        ' &
       , 'j_hyac_b        ' &
       , 'jhobr           ' &
       , 'jbro            ' &
       , 'j_br2           ' &
       , 'j_no3_aq_a      ' &
       , 'j_no3_aq_b      ' &
       , 'j_no3_aq_c      ' &
       , 'j_mek           ' &
       , 'j_ppn_a         ' &
       , 'j_ppn_b         ' &
       , 'j_hoch2ooh      ' &
       , 'j_acrol         ' &
       , 'j_ch3coooh      ' &
       , 'j_amine         ' &
       , 'j_clo_a         ' &
       , 'jclo            ' &
       , 'j_clno2         ' &
       , 'j_brno          ' &
       , 'j_brno2         ' &
       , 'jbrono_b        ' &
       , 'jbrono_a        ' &
       , 'j_nocl          ' &
       , 'joclo           ' &
       , 'jbrcl           ' &
       , 'j_ch3oono2      ' &
       , 'j_clono         ' &
       , 'j_ch2o_r        ' &
       , 'j_ch2o_m        ' &
       , 'j_ch3cooh       ' &
       , 'j_ch3ocl        ' &
       , 'j_chcl3         ' &
       , 'j_c2h5ono2      ' &
       , 'j_nc3h7ono2     ' &
       , 'j_1c4h9ono2     ' &
       , 'j_2c4h9ono2     ' &
       , 'j_perfluoro     ' &
       , 'j_i2            ' &
       , 'j_io            ' &
       , 'j_ioh           ' &
       , 'jhocl           ' &
       , 'jh2o_a          ' &
       , 'jh2o_b          ' &
       , 'jh2o_c          ' &
       /)
  
  write(*,*) 'BEGIN TEST'

  colEnvConds => environ_conditions_create( env_conds_file, lat=45.0, lon=180. )

  nlevels = colEnvConds%nlevels()
  
  errflg=0
  errmsg=' '

  call tuv_photolysis_readnl(nml_file, errmsg, errflg)
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call molec_ox_xsect_init( errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_radiation_transfer_init( r8, nlevels, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  call tuv_photolysis_init( r8, nlevels, my_jnames, errmsg, errflg )
  if (errflg/=0) then
      write(*,*) 'FAILURE tuv_photolysis_init: '//trim(errmsg)
     call abort()
  end if

  allocate(srb_o2_xs(tuv_n_wavelen,nlevels), dto2(nlevels-1,tuv_n_wavelen))
  allocate(radfld(tuv_n_wavelen,nlevels))
  allocate(tuv_prates(nlevels, tuv_n_phot ) )
 
  allocate(alt(nlevels))
  allocate(press_mid(nlevels))
  allocate(press_int(nlevels))
  allocate(temp(nlevels))
  allocate(o2vmrcol(nlevels))
  allocate(o3vmrcol(nlevels))
  allocate(so2vmrcol(nlevels))
  allocate(no2vmrcol(nlevels))
  allocate(cldfrac(nlevels))
  allocate(cldwat(nlevels))

  zenith = colEnvConds%getsrf('SZA')
  albedo = colEnvConds%getsrf('ASDIR')
  press_mid(:nlevels) = colEnvConds%press_mid(nlevels)
  alt(:nlevels) = colEnvConds%getcol('Z3',nlevels) ! meters
  temp(:nlevels) = colEnvConds%getcol('T',nlevels)
  o2vmrcol(:nlevels) = colEnvConds%getcol('O2',nlevels)
  o3vmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  so2vmrcol(:nlevels) = colEnvConds%getcol('SO2',nlevels)
  no2vmrcol(:nlevels) = colEnvConds%getcol('NO2',nlevels)
  cldfrac = 0._r8
  cldwat = 0._r8

  call molec_ox_xsect_run( nlevels, zenith, alt, temp, press_mid, o2vmrcol, dto2, srb_o2_xs, errmsg, errflg )

  call  tuv_radiation_transfer_run( &
       zenith, albedo, press_mid, alt, temp, o3vmrcol, so2vmrcol, no2vmrcol, cldfrac, cldwat, dto2, radfld, errmsg, errflg )

  call tuv_photolysis_run( nlevels, temp, press_mid, radfld, srb_o2_xs, tuv_prates, errmsg, errflg )
  
  if (errflg/=0) then
      write(*,*) 'FAILURE: '//trim(errmsg)
     call abort()
  end if

  do i=1,tuv_n_phot

     jndx = rxn_ndx(i)
     if (jndx>0) then
        rxn_string = trim(xsqy_table(jndx)%equation)
     else
        rxn_string = 'O2 -> O + O'
     end if
     write( *,*) trim(my_jnames(i))//'   '//trim(rxn_string)
     write(10,*) trim(my_jnames(i))//'   '//trim(rxn_string)
     write( *,'("  rate = ",e12.4," /sec")' ) tuv_prates(nlevels,i)
     do k=1,nlevels
         write(10,'("  rate = ",e24.16," /sec")' ) tuv_prates(k,i)
     end do

     write(10,*) ' '
     print*,' '
  end do

  write(*,*) 'END TEST'
  
end program driver
