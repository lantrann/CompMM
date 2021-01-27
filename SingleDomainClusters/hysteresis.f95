!----------------------------------------------------------------------------------------
! @brief : hysteresis
! @author: Tran Nguyen Lan
! @date  : ...
!----------------------------------------------------------------------------------------

subroutine hysteresis( trial_kind, tmp )
  
  use eb_module, only : site_core, site_shell
  use eb_module, only : pi
  use input

  implicit none
  
  integer      , intent(in)    :: trial_kind  ! 1 : eta_max, 2 : three steps    
  real(kind=8) , intent(in)    :: tmp

  integer      :: imcs, ispin, ihext, nhext
  real(kind=8) :: hext, accpt_prob
  real(kind=8) :: E, Ea, E_trial, dE
  real(kind=8) :: mtot, mcore, mshell, mcore_interf, mshell_interf

  if (trial_kind < 1 .or. trial_kind > 2) then
     stop "trial_kind should be 1 or 2"
  end if

  nhext = int( hmax / dhext )
  ! a haft first loop
  do ihext = nhext, -nhext, -1
     hext = ihext * dhext

     do imcs = 1, nmcs
        call sweep_core (imcs, trial_kind, tmp, hext, mcore , mcore_interf , accpt_prob)
        call sweep_shell(imcs, trial_kind, tmp, hext, mshell, mshell_interf, accpt_prob)
     end do ! do imcs

     mtot          = 0.00d+00
     mcore         = 0.00d+00
     mcore_interf  = 0.00d+00
     mshell        = 0.00d+00
     mshell_interf = 0.00d+00
     do imcs = 1, nmcs
        call sweep_core (imcs, trial_kind, tmp, hext, mcore , mcore_interf , accpt_prob)
        call sweep_shell(imcs, trial_kind, tmp, hext, mshell, mshell_interf, accpt_prob)
     end do ! do imcs

     mtot  = mcore + mshell 
     mtot  = mtot  / ( nmcs * ( size(site_core) + size(site_shell) ))
     mcore = mcore / ( nmcs * size(site_core) )
     !mcore_interf = mcore_interf / ( nmcs * size(site_core) )
     print*, hext, mtot, mcore
  end do 

  ! second later loop
  do ihext = -nhext, nhext, 1
     hext = ihext * dhext

     do imcs = 1, nmcs
        call sweep_core (imcs, trial_kind, tmp, hext, mcore , mcore_interf , accpt_prob)
        call sweep_shell(imcs, trial_kind, tmp, hext, mshell, mshell_interf, accpt_prob)
     end do ! do imcs

     mtot          = 0.00d+00
     mcore         = 0.00d+00
     mcore_interf  = 0.00d+00
     mshell        = 0.00d+00
     mshell_interf = 0.00d+00
     do imcs = 1, nmcs
        call sweep_core (imcs, trial_kind, tmp, hext, mcore , mcore_interf , accpt_prob)
        call sweep_shell(imcs, trial_kind, tmp, hext, mshell, mshell_interf, accpt_prob)
     end do ! do imcs

     mtot  = mcore + mshell 
     mtot  = mtot  / ( nmcs * ( size(site_core) + size(site_shell) ))
     mcore = mcore / ( nmcs * size(site_core) )
     !mcore_interf = mcore_interf / ( nmcs * size(site_core) )
     print*, hext, mtot, mcore
  end do 
  
end subroutine hysteresis

