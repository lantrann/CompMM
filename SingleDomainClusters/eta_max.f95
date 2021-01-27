! @brief :  eta_max procedure for spin update
! @author: Tran Nguyen Lan
! @date  : ...

subroutine eta_max (eta, theta, theta_trial, phi, phi_trial)

  use eb_module, only : pi
  implicit none
  
  real(kind=8), intent(in)    :: eta
  real(kind=8), intent(in  )  :: theta, phi
  real(kind=8), intent(inout) :: theta_trial, phi_trial
  
  real(kind=8) :: rdn1, rdn2
  real(kind=8) :: dtheta
  real(kind=8) :: dphi
  
  ! radomly chosen the trial orientation from uniform distribution [-eta, eta]
  call random_number(rdn1)
  call random_number(rdn2)
  dtheta = eta * (2.0d+00 * rdn1 - 1)
  dphi   = 2.0d+00 * pi * rdn2
  
  theta_trial = theta + dtheta
  phi_trial   = phi   + dphi
  
end subroutine eta_max
