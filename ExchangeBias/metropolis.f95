!----------------------------------------------------
! @brief : Metropolis algorithm 
! @author: Tran Nguyen Lan
! @date  : ...
!----------------------------------------------------

subroutine metropolis(dE, T, naccept,         & ! deviation of energy, temperature, # of accepted spins
                      theta, phi ,            & ! old and trial orientation
                      theta_trial, phi_trial)   ! old and trial orientation

! In Metropolis algorithm, the new trial will be generate, if the energy converges (Etrial < E or dE < 0),
! it means that the system tends to thermal equilibrium => new trial is accepted. 
! Else, if dE > 0, the new trial will be accepted with probability p = min[1,exp(- dE / (kB * T) )]

real(kind=8), intent(inout) :: theta, phi
real(kind=8), intent(in)    :: theta_trial, phi_trial
real(kind=8), intent(in)    :: dE, T
integer,      intent(inout) :: naccept

real(kind=8) :: rdn, w

! we can use only one "IF loop" but should avoid very huge probability (no physical meaning) when dE < 0  
if (dE < 0) then
   theta = theta_trial
   phi   = phi_trial
   naccept = naccept + 1
else
   call random_number(rdn)
   ! calculate Boltzmann thermal probability w
   w = exp( - dE / T )
   !print*, w, rdn
   if (w >= rdn) then
      theta = theta_trial
      phi   = phi_trial
      naccept = naccept + 1
   end if
end if

end subroutine metropolis
