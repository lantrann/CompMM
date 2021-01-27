!-------------------------------------------------------------------------
! @brief  : sweeping core spins
! @author: Tran Nguyen Lan
! @date   : ...
!-------------------------------------------------------------------------

subroutine sweep_core(imcs, trial_kind, tmp, h, mcore, mcore_interf, accpt_prob)

  use eb_module, only : site_core
  use input    , only : Jc, Jsh, Jint, Kc, Ksh, eta, open_angle

  implicit none
  
  integer     , intent(in)    :: trial_kind ! 1: eta_max, 2: three steps
  integer     , intent(in)    :: imcs ! current MC step
  real(kind=8), intent(in)    :: tmp, h ! temperature and external field
  real(kind=8), intent(inout) :: mcore, mcore_interf
  real(kind=8), intent(out)   :: accpt_prob

  real(kind=8) :: Ea, Eh, E, Etrial, dE
  real(kind=8) :: Exc_core_core, Exc_core_shell
  integer      :: ispin, jspin, k, naccept
  
  !print*, 'we are here'
  naccept = 0
  do ispin = 1, size(site_core)
     
     ! first, we calculate the energy with given orientation
     !------------------------------------------------------
     
     ! zeeman energy
     call zeemann_z (site_core(ispin)%s, site_core(ispin)%theta, h, Eh) 
     
     ! uniaxial anisotropy energy
     call uniaxial_anisotropy_z (site_core(ispin)%s, site_core(ispin)%theta, Kc, Ea)
     
     ! exchange coupling in core
     call exchange_coupling (Jc, Exc_core_core, ispin, 1, 1) ! itype = 1, itrial = 1
     
     ! exchange coupling at interface
     if ( size( site_core(ispin)%list_nrcs ) /= 0 ) then
        call exchange_coupling (Jint, Exc_core_shell, ispin, 1, 3)  ! itype = 3, itrial = 1
     else
        Exc_core_shell = 0
     end if
     
     ! total energy
     E = Exc_core_core + Exc_core_shell + Ea + Eh
     !E = Exc_core_core + Ea + Eh
     !E = Ea + Eh

     !updating spin orientation
     !------------------------------------------------------

     ! call trial procedure: eta_max or three_steps
     if (trial_kind == 1) then
        call eta_max ( eta, site_core(ispin)%theta, site_core(ispin)%theta_trial, & 
                            site_core(ispin)%phi  , site_core(ispin)%phi_trial )
     else if (trial_kind == 2) then
        call three_trial_steps ( imcs, open_angle, &
                                 site_core(ispin)%theta, site_core(ispin)%theta_trial, & 
                                 site_core(ispin)%phi  , site_core(ispin)%phi_trial )
     end if

     !print*, site_core(ispin)%theta, site_core(ispin)%theta_trial

     ! now, turn to calculate energy with trial orientation 
     !------------------------------------------------------
     
     ! zeeman energy 
     call zeemann_z (site_core(ispin)%s, site_core(ispin)%theta_trial, h, Eh) ! in case of zfc, h = 0
     
     ! uniaxial anisotropy energy
     call uniaxial_anisotropy_z (site_core(ispin)%s, site_core(ispin)%theta_trial, Kc, Ea)
     
     ! exchange coupling in core
     call exchange_coupling (Jc, Exc_core_core, ispin, 2, 1) ! itype = 1, itrial = 2
     
     ! exchange coupling at interface
     if ( size( site_core(ispin)%list_nrcs ) /= 0 ) then
        call exchange_coupling (Jint, Exc_core_shell, ispin, 2, 3)  ! itype = 3, itrial = 2
     else 
        Exc_core_shell = 0
     end if
     
     ! total trial energy
     Etrial = Exc_core_core + Exc_core_shell + Ea + Eh
     !Etrial = Exc_core_core + Ea + Eh
     !Etrial = Ea + Eh

     ! calculate difference of energy and using Metropolis algorithm for chosing the trial
     !------------------------------------------------------------------------------------
     dE = Etrial - E
     !print*, Etrial, E
     !print*, eta, site_core(ispin)%theta, site_core(ispin)%theta_trial, dE
     call metropolis( dE, tmp, naccept, &                                
                      site_core(ispin)%theta      , site_core(ispin)%phi , & 
                      site_core(ispin)%theta_trial, site_core(ispin)%phi_trial )
     
     !magnetization of core spins
     mcore = mcore + cos( site_core(ispin)%theta )

     !magnetization of core spins at interface
     if ( size( site_core(ispin)%list_nrcs ) /= 0 ) then
        mcore_interf = mcore_interf + cos( site_core(ispin)%theta )
     end if

  end do ! do ispin

  accpt_prob =  dfloat( naccept ) / dfloat( size(site_core) ) ! accepted probability
  !print*, accpt_prob
  !print*, tmp, mcore
  
end subroutine sweep_core
