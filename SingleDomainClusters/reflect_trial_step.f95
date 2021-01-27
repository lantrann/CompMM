! @brief : reflect trial steps (RTS) procedure for spin update
! @author: Le Bin Ho
! @date  : ...

subroutine reflect_trial_step (theta, theta_trial, phi, phi_trial)

  use eb_module, only : pi
  implicit none

   real(kind=8), intent(in)  :: theta, phi
   real(kind=8), intent(inout) :: theta_trial, phi_trial

   theta_trial = pi - theta
   phi_trial   = phi

 end SUBROUTINE reflect_trial_step
