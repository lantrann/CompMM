
module input

  implicit none

  real(kind=8), parameter :: L     = 12.00d+00  
  real(kind=8), parameter :: a     =  1.00d+00  
  real(kind=8), parameter :: Jc    = 10.00d+00   
  real(kind=8), parameter :: Jint  =  5.00d+00 ! + or -
  real(kind=8), parameter :: Jsh   = -5.00d+00   
  real(kind=8), parameter :: hmax  =  4.00d+00  
  real(kind=8), parameter :: Kc    =  1.00d+00   
  real(kind=8), parameter :: Ksh   = 10.00d+00  
  real(kind=8), parameter :: tmax  =  5.00d+00  
  real(kind=8), parameter :: tmes  =  0.1d+00  
  real(kind=8), parameter :: h_fc  = -4.00d+00  
  real(kind=8), parameter :: dtmp  =  0.10d+00  
  real(kind=8), parameter :: dhext =  0.10d+00  
  real(kind=8), parameter :: eta   =  1.00d+00
  real(kind=8), parameter :: open_angle   =  0.1 !3.14159 / 24.00d+00 

  integer     , parameter :: nequ  = 100  ! first MC steps for thermal equilibrium
  integer     , parameter :: nmcs  = 200   ! MC steps for measurment
  
end module input

