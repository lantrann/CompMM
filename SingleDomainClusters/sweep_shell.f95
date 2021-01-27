!-------------------------------------------------------------------------
! @brief  : sweeping shell spins
! @author: Tran Nguyen Lan
! @date   : ...
!-------------------------------------------------------------------------

subroutine sweep_shell(imcs, trial_kind, tmp, h, mshell, mshell_interf, accpt_prob)

  use eb_module, only : site_shell
  use input    , only : Jc, Jsh, Jint, Kc, Ksh, eta, open_angle

  implicit none
  
  integer     , intent(in)    :: trial_kind ! 1: eta_max, 2: three steps
  integer     , intent(in)    :: imcs ! current MC step
  real(kind=8), intent(in)    :: tmp, h ! temperature and external field
  real(kind=8), intent(inout) :: mshell, mshell_interf
  real(kind=8), intent(out)   :: accpt_prob

  real(kind=8) :: Ea, Eh, E, Etrial, dE
  real(kind=8) :: Exc_shell_shell, Exc_shell_core
  integer      :: ispin, jspin, k, naccept
  
  !print*, 'we are here'
  naccept = 0
  do ispin = 1, size(site_shell)
     
     ! first, we calculate the energy with given orientation
     !------------------------------------------------------
     
     ! zeeman energy
     call zeemann_z (site_shell(ispin)%s, site_shell(ispin)%theta, h, Eh) 
     
     ! uniaxial anisotropy energy
     call uniaxial_anisotropy_z (site_shell(ispin)%s, site_shell(ispin)%theta, Ksh, Ea)
     
     ! exchange coupling in shell
     call exchange_coupling (Jsh, Exc_shell_shell, ispin, 1, 2) ! itype = 2, itrial = 1
     
     ! exchange coupling at interface
     if ( size( site_shell(ispin)%list_nrsc ) /= 0 ) then
        call exchange_coupling (Jint, Exc_shell_core, ispin, 1, 4)  ! itype = 4, itrial = 1
     else
        Exc_shell_core = 0
     end if
     
     ! total energy
     E = Exc_shell_shell + Exc_shell_core + Ea + Eh
     !E = Exc_shell_shell + Ea + Eh
     !E = Ea + Eh

     !updating spin orientation
     !------------------------------------------------------

     ! call trial procedure: eta_max or three_steps
     if (trial_kind == 1) then
        call eta_max ( eta, site_shell(ispin)%theta, site_shell(ispin)%theta_trial, & 
                            site_shell(ispin)%phi  , site_shell(ispin)%phi_trial )
     else if (trial_kind == 2) then
        call three_trial_steps ( imcs, open_angle, &
                                 site_shell(ispin)%theta, site_shell(ispin)%theta_trial, & 
                                 site_shell(ispin)%phi  , site_shell(ispin)%phi_trial )
     end if

     !print*, site_shell(ispin)%theta, site_shell(ispin)%theta_trial

     ! now, turn to calculate energy with trial orientation 
     !------------------------------------------------------
     
     ! zeeman energy 
     call zeemann_z (site_shell(ispin)%s, site_shell(ispin)%theta_trial, h, Eh) ! in case of zfc, h = 0
     
     ! uniaxial anisotropy energy
     call uniaxial_anisotropy_z (site_shell(ispin)%s, site_shell(ispin)%theta_trial, Ksh, Ea)
     
     ! exchange coupling in shell
     call exchange_coupling (Jsh, Exc_shell_shell, ispin, 2, 2) ! itype = 2, itrial = 2
     
     ! exchange coupling at interface
     if ( size( site_shell(ispin)%list_nrsc ) /= 0 ) then
        call exchange_coupling (Jint, Exc_shell_core, ispin, 2, 4)  ! itype = 4, itrial = 2
     else 
        Exc_shell_core = 0
     end if
     
     ! total trial energy
     Etrial = Exc_shell_shell + Exc_shell_core + Ea + Eh
     !Etrial = Exc_shell_shell + Ea + Eh
     !Etrial = Ea + Eh

     ! calculate difference of energy and using Metropolis algorithm for chosing the trial
     dE = Etrial - E
     !print*, Etrial, E
     !print*, eta, site_shell(ispin)%theta, site_shell(ispin)%theta_trial, dE
     call metropolis( dE, tmp, naccept, &                                
                      site_shell(ispin)%theta      , site_shell(ispin)%phi , & 
                      site_shell(ispin)%theta_trial, site_shell(ispin)%phi_trial )
     
     !magnetization of shell spins
     mshell = mshell + cos( site_shell(ispin)%theta )

     !magnetization of shell spins at interface
     if ( size( site_shell(ispin)%list_nrsc ) /= 0 ) then
        mshell_interf = mshell_interf + cos( site_shell(ispin)%theta )
     end if

  end do ! do ispin

  accpt_prob =  dfloat( naccept ) / dfloat( size(site_shell) ) ! accepted probability
  !print*, accpt_prob
  !print*, tmp, mshell
  
end subroutine sweep_shell
