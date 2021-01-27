!------------------------------------------------------------------
! @brief : calculate the uniaxial anisotropy energy
!         ( with the easy axis parrallel to z axis )
! @author: Le Bin Ho
! @date  : ...
!------------------------------------------------------------------

subroutine uniaxial_anisotropy_z(s, theta, k, Ea)
  
  real(kind=8), intent(in)  :: s, theta, k
  real(kind=8), intent(out) :: Ea
  
  Ea = - k * ( s * cos(theta) ) ** 2

end subroutine uniaxial_anisotropy_z
