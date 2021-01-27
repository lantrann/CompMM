!------------------------------------------------------------------
! @brief : small trial steps (STS) procedure for spin update
! @author: Le Bin Ho
! @date  : ...
!------------------------------------------------------------------

SUBROUTINE small_trial_step (open_angle, theta, theta_trial, phi, phi_trial)

   use eb_module, only : pi
   implicit none

   real(kind=8), intent(in)    :: theta, phi
   real(kind=8), intent(inout) :: theta_trial, phi_trial
   real(kind=8), intent(in)    :: open_angle

   real(kind=8)	:: x, y, z, u, v, w, tcos ! tcos = 1- cos(d_phi)
   real(kind=8)	:: R(1:3,1:3)	          ! Rotation random matrix
   real(kind=8)	:: A(1:3)	          ! New vector, Ax, Ay, Az
   real(kind=8)	:: re		          ! renomalization
   real(kind=8)	:: rdn1, rdn2, rdn3 
   real(kind=8) :: d_theta, d_phi, d_phi_o

   CALL Random_number(rdn1)
   d_theta = open_angle * ( 2 * rdn1 - 1 )
   CALL Random_number(rdn2)
   CALL Random_number(rdn3)
   d_phi_o = 2.0d+00 * pi * rdn2 
   d_phi   = 2.0d+00 * pi * rdn3 
   tcos = 1.0 - cos(d_phi)
   
   x = sin(theta) * cos(phi)
   y = sin(theta) * sin(phi)
   z = cos(theta)

   !Create a rotation random matrix 3 x 3!
   R(1,1) = cos(d_phi) + tcos*(x**2)
   R(1,2) = x * y * tcos - z * sin(d_phi)
   R(1,3) = x * z * tcos + y * sin(d_phi)
   R(2,1) = y * x * tcos + z * sin(d_phi)
   R(2,2) = cos(d_phi) + tcos*(y**2)
   R(2,3) = y * z * tcos - x * sin(d_phi)
   R(3,1) = z * x * tcos - y * sin(d_phi)
   R(3,2) = z * y * tcos + x * sin(d_phi)
   R(3,3) = cos(d_phi) + tcos * (z**2)
   !print*, R(1,1),R(1,2),R(1,3),"hang 1"
   !print*, R(2,1),R(2,2),R(2,3)
   !print*, R(3,1),R(3,2),R(3,3)
   !print*, R(1,1)*R(2,2)*R(3,3) + R(1,2)*R(2,3)*R(3,1) +R(1,3)*R(2,1)*R(3,2) -&
   !	R(1,3)*R(2,2)*R(3,1) - R(1,2)*R(2,1)*R(3,3)-R(1,1)*R(2,3)*R(3,2)  

   !Rotate initial vector with d_theta   !
   A(1) = sin(theta + d_theta) * cos(phi + d_phi_o)
   A(2) = sin(theta + d_theta) * sin(phi + d_phi_o)
   A(3) = cos(theta + d_theta)

   u = R(1,1) * A(1) + R(1,2) * A(2) + R(1,3)*A(3)
   v = R(2,1) * A(1) + R(2,2) * A(2) + R(2,3)*A(3)
   w = R(3,1) * A(1) + R(3,2) * A(2) + R(3,3)*A(3)
   re = u * u + v * v + w * w

   ! Renomalization
   u = u / sqrt(re)
   v = v / sqrt(re)
   w = w / sqrt(re)
   
   ! trial orientation of vector magnetic moment 
   theta_trial = acos(w)
   phi_trial   = atan2(v,u)
   !print*, d_theta,d_phi, spins_array(2,i), spins_array(3,i),theta_new, phi_new ,i
   
 end SUBROUTINE small_trial_step
 
