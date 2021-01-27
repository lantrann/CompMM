! @brief : uniform trial steps (UTS) procedure for spin update
! @author: Le Bin Ho
! @date  : ...

SUBROUTINE uniform_trial_step (theta_trial, phi_trial)

  use eb_module, only : pi
  implicit none
  
  real(kind=8), intent(inout) :: theta_trial, phi_trial
  real(kind=8) :: rdn1, rdn2
  
  CALL Random_number(rdn1)
  CALL Random_number(rdn2)
  theta_trial = 1.0d+00 * pi * rdn1
  phi_trial   = 2.0d+00 * pi * rdn2
  
end SUBROUTINE uniform_trial_step
