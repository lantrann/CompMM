!------------------------------------------------------------------
! @brief : calculate the Zeemann energy
!        ( with the external field  parrallel to z axis ) 
! @author: Le Bin Ho
! @date  : ...
!------------------------------------------------------------------

subroutine zeemann_z(s, theta, h, Eh)
  
  real(kind=8), intent(in)  :: s, theta, h
  real(kind=8), intent(out) :: Eh
  
  Eh = - h * s * cos(theta) 

end subroutine zeemann_z
