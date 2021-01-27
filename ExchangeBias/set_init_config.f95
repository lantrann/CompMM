!----------------------------------------------------------------------------------------
! @brief : set initial configuration (zfc or fc) before running hysteresis
! @author: Tran Nguyen Lan
! @date  : ...
!----------------------------------------------------------------------------------------

subroutine set_init_config( trial_kind, config_kind )
  
  use eb_module, only : site_core, site_shell
  use eb_module, only : pi
  use input

  implicit none
  
  integer, intent(in)    :: trial_kind  ! 1 : eta_max, 2 : three steps    
  integer, intent(in)    :: config_kind ! 1 : zfc, 2 : fc

  integer      :: imcs, ispin, itmp, ntmp, ntmes
  real(kind=8) :: h, tmp, accpt_prob
  real(kind=8) :: E, Ea, E_trial, dE
  real(kind=8) :: mcore, mshell, mcore_interf, mshell_interf

  if (trial_kind < 1 .or. trial_kind > 2) then
     stop "trial_kind should be 1 or 2"
  end if

  if (config_kind < 1 .or. config_kind > 2) then
     stop "config_kind should be 1 or 2"
  end if

  if (config_kind == 1 ) then
     h = 0 ! zfc
  else if (config_kind == 2) then
     h = h_fc ! fc
  end if

  ntmp  = int( tmax / dtmp )
  ntmes = int( tmes / dtmp )
  do itmp = ntmp, ntmes, -1
     tmp = itmp * dtmp
     mcore         = 0.00d+00
     mcore_interf  = 0.00d+00
     mshell        = 0.00d+00
     mshell_interf = 0.00d+00
     do imcs = 1, nequ
        call sweep_core (imcs, trial_kind, tmp, h, mcore , mcore_interf , accpt_prob)
        call sweep_shell(imcs, trial_kind, tmp, h, mshell, mshell_interf, accpt_prob)
     end do ! do imcs

     mcore  = mcore  / ( nequ * size(site_core ) )
     mshell = mshell / ( nequ * size(site_shell) )
     !mcore_interf = mcore_interf / ( nmcs * size(site_core) )
     print*, tmp, mcore, mshell
  end do 
  
end subroutine set_init_config

