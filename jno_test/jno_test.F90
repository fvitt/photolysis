program driver

  use machine, only: r8 => kind_phys
  use environ_conditions_mod, only: environ_conditions_create, environ_conditions
  use jno_mod

  implicit none

  integer :: i,jndx
  integer :: errflg
  character(len=444) :: errmsg

  real(r8) :: zenith
  real(r8), allocatable :: alt(:)
  real(r8), allocatable :: press(:)
  real(r8), allocatable :: temp(:)
  real(r8), allocatable :: n2vmrcol(:)
  real(r8), allocatable :: o2vmrcol(:)
  real(r8), allocatable :: o3vmrcol(:)
  real(r8), allocatable :: novmrcol(:)
  real(r8), allocatable :: jno(:)

  character(len=*), parameter :: env_conds_file = &
       '/terminator-data1/fvitt/micm_inputs/MusicBox_env_cond_1col_c190109.nc'
  type(environ_conditions),pointer :: colEnvConds => null()

  integer :: nlevels,k

  real(r8) :: fluxes(jno_nbins) = (/ 2.1471364E+11,  2.1680430E+11,  2.8858957E+11,  3.9553672E+11 /)

  
  write(*,*) 'BEGIN TEST'

  colEnvConds => environ_conditions_create( env_conds_file, lat=45.0, lon=180. )

  nlevels = colEnvConds%nlevels()
  
  errflg=0
  errmsg=' '
 
  allocate(alt(nlevels))
  allocate(press(nlevels))
  allocate(temp(nlevels))
  allocate(o2vmrcol(nlevels))
  allocate(o3vmrcol(nlevels))
  allocate(n2vmrcol(nlevels))
  allocate(novmrcol(nlevels))
  allocate(jno(nlevels))

  zenith = colEnvConds%getsrf('SZA')
  press(:nlevels) = colEnvConds%press_mid(nlevels)
  alt(:nlevels) = colEnvConds%getcol('Z3',nlevels) *1.e-3 ! kmeters
  temp(:nlevels) = colEnvConds%getcol('T',nlevels)
  n2vmrcol(:nlevels) = colEnvConds%getcol('N2',nlevels)
  o2vmrcol(:nlevels) = colEnvConds%getcol('O2',nlevels)
  o3vmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  novmrcol(:nlevels) = colEnvConds%getcol('O3',nlevels)
  
  call jno_init()

  call jno_timestep_init( fluxes )

  call jno_run( nlevels, zenith, n2vmrcol, o2vmrcol, o3vmrcol, novmrcol, press, temp, alt, jno )

  print*,'jno: ',jno

  do k = 1, nlevels
     write(10,fmt='(a,i2,a,e22.16)') ' jno(', k, '): ',jno(k)
  end do
  
  write(*,*) 'END TEST'
  
end program driver
