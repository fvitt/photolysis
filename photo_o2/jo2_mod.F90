module jo2_mod
!==============================================================================!
!   JO2 parameterizations based on:                          !
!        Lyman alpha... Chabrillat and Kockarts, GRL, 25, 2659, 1998           !
!        SRC .......... Brasseur and Solomon, 1986 (from TUV)                  !
!        SRB .......... Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996    !
!                        (supplied by Dan Marsh, NCAR ACD                      !
!------------------------------------------------------------------------------
!     ... Derive the O2 rate constant and apply branching ratio (QY)
!------------------------------------------------------------------------------
!     ... SRC and SRB QY
!         Longward  of 174.65 the product is O2 + hv => O(3P) + O(3P)
!         Shortward of 174.65 the product is O2 + hv => O(3P) + O(1D)
!         The QY is assumed to be unity in both wavelength ranges.
!
!     ... Lyman Alpha QY
!         O2 + hv -> O(3P) + O(3P) at Lyman Alpha has a QY = 0.47
!         O2 + hv -> O(3P) + O(1D) at Lyman Alpha has a QY = 0.53
!         Lacoursiere et al., J. Chem. Phys. 110., 1949-1958, 1999.
!------------------------------------------------------------------------------
!     ... Molecular Oxygen, SRB
!------------------------------------------------------------------------------
!     ... Koppers Grid (see Koppers and Murtagh, Ann. Geophys., 14, 68-79, 1996)
!        #    wl(i)       wl(i+1)
!        1   174.4        177.0
!        2   177.0        178.6
!        3   178.6        180.2
!        4   180.2        181.8
!        5   181.8        183.5
!        6   183.5        185.2
!        7   185.2        186.9
!        8   186.9        188.7
!        9   188.7        190.5
!        10  190.5        192.3
!        11  192.3        194.2
!        12  194.2        196.1
!        13  196.1        198.0
!        14  198.0        200.0  <<last wl bin for <200nm
!        ----------------------
!        15  200.0        202.0
!        16  202.0        204.1
!        17  204.1        205.8
!------------------------------------------------------------------------------

  implicit none

  integer, parameter :: r8 = 8

  real(r8), allocatable :: etfphot(:) ! /cm2/sec
  integer :: la_num

  integer :: n_src, n_srb

  integer :: nw

contains

  subroutine jo2_init( wc, wl )

    real(r8), intent(in) :: wc(:)
    real(r8), intent(in) :: wl(:)

    integer :: n

    nw = size( wc )
    allocate( etfphot(nw) )

    n_src = -1

    find_src: do n = 1,nw
       if (wc(n) < 174.65) then
          n_src = n
          exit find_src
       endif
    end do find_src

    n_srb = nw - n_src ! Schumann-Runge 
 
    find_la: do  n = 1,nw
       if (wc(n) > 121.4 .and. wc(n) < 121.9 ) then
          la_num = n
          exit find_la
       endif
    end do find_la

  end subroutine jo2_init

  subroutine jo2_set_etf( etf )
    real(r8), intent(in) :: etf(:)

    etfphot(:nw) = etf(:nw)
    
  end subroutine jo2_set_etf
  
  subroutine jo2_calc( debug, nlev, srb_o2_xs, fnorm, ro2la, jo2 )

    logical,  intent(in) :: debug
    integer,  intent(in) :: nlev
    real(r8), intent(in) :: srb_o2_xs(nw,nlev)
    real(r8), intent(in) :: fnorm(nw,nlev)
    real(r8), intent(in) :: ro2la(nlev)
    real(r8), intent(out) :: jo2(nlev,2) ! 2 branches for O2 dissociation
    
    integer :: k

    real(r8) :: wrk(nw)	     ! wrk array
    real(r8) :: jo2_lya(nlev)
    real(r8) :: jo2_src(nlev)
    real(r8) :: jo2_srb(nlev)
    
    !------------------------------------------------------------------------------
    !     ... Lyman Alpha
    !------------------------------------------------------------------------------
    jo2_lya(:nlev) = etfphot(la_num)*ro2la(:nlev)

    if (debug) then
!!$       write(*,*) 'jo2_mod  etfphot: ',etfphot
!!$       write(*,*) 'jo2_mod  la_num : ',la_num
!!$       
       write(*,*) 'jo2_mod  etfphot(la_num): ',etfphot(la_num)
       write(*,*) 'jo2_mod  ro2la: ',ro2la
       write(*,*) 'jo2_mod  jo2_lya: ',jo2_lya
    end if
    
    !------------------------------------------------------------------------------
    !     ... o2 src and srb photolysis
    !------------------------------------------------------------------------------
    do k = 1,nlev
       wrk(1:nw) = srb_o2_xs(1:nw,k)*etfphot(1:nw)
       wrk(la_num)= 0._r8
       
       jo2_src(k) = dot_product( fnorm(1:n_src   ,k),wrk(1:n_src   ) )
       jo2_srb(k) = dot_product( fnorm(n_src+1:nw,k),wrk(n_src+1:nw) )
    end do

    !------------------------------------------------------------------------------
    !     ... total o2 photolysis
    !------------------------------------------------------------------------------
    !------------------------------------------------------------------------------
    !     ... Branch 1, O2 + hv => O(3P) + O(3P); wavelengths >175nm
    !------------------------------------------------------------------------------
    jo2(nlev:1:-1,1) = jo2_lya(:)*.47_r8 !+ jo2_srb(:)
    !------------------------------------------------------------------------------
    !     ... Branch 2, O2 + hv => O(3P) + O(1D);  wavelengths <175nm
    !------------------------------------------------------------------------------
    jo2(nlev:1:-1,2) = jo2_lya(:)*.53_r8 !+ jo2_src(:)

  end subroutine jo2_calc
    
end module jo2_mod
