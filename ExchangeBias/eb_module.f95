!----------------------------------------------------------------------------------------
! @brief : information about atomic spin in  magnetic core/shell nanoparticles
! @author: Tran Nguyen Lan
! @date  : ...
!----------------------------------------------------------------------------------------

module eb_module

  real(kind=8) :: pi = 3.1415926535897932384626433832795d+00

  ! information about the 'site-th' spin at core
  type info_site_core
     real(kind=8)     :: r(1:3)       ! coordinate of atomic spins
     real(kind=8)     :: s            ! module of vector magnetic moment
     real(kind=8)     :: theta        ! orientation with respect to z axis
     real(kind=8)     :: phi          ! orientation in xy plane
     real(kind=8)     :: theta_trial  ! trial orientation with respect to z axis
     real(kind=8)     :: phi_trial    ! trial orientation in xy plane
     integer, allocatable :: list_nrcc(:) ! list of neareast spins core-core
     integer, allocatable :: list_nrcs(:) ! list of neareast spins core-shell (interface)
  end type info_site_core

  ! information about the 'site-th' spin at shell 
  type info_site_shell
     real(kind=8)     :: r(1:3)       ! coordinate of atomic spins
     real(kind=8)     :: s            ! module of vector magnetic moment (atomic spin)
     real(kind=8)     :: theta        ! orientation with respect to z axis
     real(kind=8)     :: phi          ! orientation in xy plane
     real(kind=8)     :: theta_trial  ! trial orientation with respect to z axis
     real(kind=8)     :: phi_trial    ! trial orientation in xy plane
     integer, allocatable :: list_nrss(:) ! list of neareast spins shell-shell
     integer, allocatable :: list_nrsc(:) ! list of neareast spins shell-core (interface)
  end type info_site_shell
  
  type(info_site_core ), allocatable :: site_core (:) ! spin core
  type(info_site_shell), allocatable :: site_shell(:) ! spin shell

end module eb_module
