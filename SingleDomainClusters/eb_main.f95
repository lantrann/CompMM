!----------------------------------------------------------------------------------------
! @brief : main program for Exchang Bias simulation in magnetic core/shell nanoparticles
! @author: Le Bin Ho
! @author: Tran Nguyen Lan
! @date  : ...
!----------------------------------------------------------------------------------------

program EB_main 
  
  use eb_module, only: site_core, site_shell
  use eb_module, only: pi
  use input

  implicit none
  
  real(kind=8) :: Rtot, Rsh ! radius of whole particle and thick of shell
  integer :: ispin, trial_kind, config_kind
  real(kind=8) :: tmp

  !open(unit=10,file='coor_test.dat',status='new')

  Rtot = 11*a ! radius of particle
  Rsh  =  2*a ! thick of shell
  
  ! generate information of atomic spins in core/shell particle
  call set_info_site ( L, a, Rtot, Rsh )                   

  trial_kind  = 2 ! 1: eta_max; 2: three steps
  config_kind = 2 ! 1: zfc or ; 2: fc
  call  set_init_config( trial_kind, config_kind )         

  call hysteresis ( trial_kind, tmes )

end program EB_main

