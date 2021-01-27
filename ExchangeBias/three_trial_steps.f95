! @brief :  eta_max procedure for spin update
! @author: Tran Nguyen Lan
! @author: Le Bin Ho
! @date  : ...

subroutine three_trial_steps (imcs, open_angle, theta, theta_trial, phi, phi_trial)

  implicit none
  
  integer     , intent(in)    :: imcs
  real(kind=8), intent(in)    :: theta, phi, open_angle
  real(kind=8), intent(inout) :: theta_trial, phi_trial
  
  integer :: i, trial_step

  do i = 0 ,4
     if (mod(imcs,5) .eq. i)then
        trial_step = i
     endif
  enddo
  
  select case( trial_step )
     
  case(1)
     !Small trial step (STS)                                                                                                                                                                                                      
     call small_trial_step (open_angle, theta, theta_trial, phi, phi_trial)
     
  case(3)
     !Reflection trial step (RTS)                                                                                                                                                                                                 
     call reflect_trial_step (theta, theta_trial, phi, phi_trial)
     
  case default 
     !Uniform trial step (UTS)                                                                                                                                                                                                    
     call uniform_trial_step (theta_trial, phi_trial)
  end select
  
end subroutine three_trial_steps
