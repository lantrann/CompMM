!-------------------------------------------------------------------------
! @brief  : set up information of atomic spin in core/shell particles
! @author: Tran Nguyen Lan
! @date   : ...
!-------------------------------------------------------------------------

subroutine set_info_site ( L, a, R, Rsh)

  ! in this subroutine, we will generate essential information for atomic spins in core-shell nanoparticles: 
  !    + coordinate of atomic spins (r)
  !    + module of atomic spins (m)     
  !    + orientation of atomic spin (theta, phi)
  !    + nearest neighbor spins in core, shell and core-shell interface

  ! Algorithm:
  !    + select and store to array
  !    + set information (reduce size of seleted array at the same time))

  use eb_module, only : site_core, site_shell
  use eb_module, only : pi
  implicit none

  real(kind=8), intent(in) :: L, a    
  real(kind=8), intent(in) :: R, Rsh  
  
  ! coordinate of atomic spins
  real(kind=8), allocatable  :: coord_core (:,:)
  real(kind=8), allocatable  :: coord_shell(:,:)

  ! vector magnetic moments
  real(kind=8) :: m
  real(kind=8), allocatable :: theta_core (:), phi_core (:)
  real(kind=8), allocatable :: theta_shell(:), phi_shell(:)

  !intermediate variables
  integer, allocatable :: tmp_coord_core (:,:)
  integer, allocatable :: tmp_coord_shell(:,:)
  integer, target :: tmp_list_nearest_spin(1:6)

  real(kind=8) :: rdn
  integer :: N, Nn, Nc, Nsh
  integer :: Rc, rsqt
  integer :: i, j, counter
  integer :: xrdn, yrdn, zrdn
  integer, allocatable :: x (:), y (:), z (:)


  N = (2*(L-1)+1)**3	! total number of spins in cubic of size L
  allocate( tmp_coord_core (N, 3) ) 
  allocate( tmp_coord_shell(N, 3) ) 
  allocate( x(N), y(N), z(N) )      

  ! create an array which stores all positons in a simle cubic lattice of size L
  do i = 1 , N
15   CALL Random_number(rdn)
     xrdn = (2*rdn - 1)*L
     xrdn = xrdn*a
     CALL Random_number(rdn)
     yrdn = (2*rdn - 1)*L
     yrdn = yrdn*a
     CALL Random_number(rdn)
     zrdn = (2*rdn - 1)*L
     zrdn = zrdn*a
     
     if ( i > 1 ) then ! begin avoiding the overlap from second loop
        do j = 1 , i-1
           if ((xrdn .EQ. x(j)).and.(yrdn .EQ. y(j)).and.(zrdn.eq.z(j))) then
              goto 15
           endif
        enddo ! do j
     end if ! if i > 1

     x(i) = xrdn
     y(i) = yrdn
     z(i) = zrdn
  enddo ! do i

  !-----------------------------
  !   coordinate of spins
  !----------------------------
  
  ! create core-shell nanoparticle
  Rc  = R - Rsh
  Nc  = 0
  Nsh = 0
  do i = 1, N
     rsqt = x(i)**2 + y(i)**2 + z(i)**2
     if (rsqt .le. Rc**2) then 
        Nc = Nc + 1
        tmp_coord_core(Nc,1:3)  = (/x(i), y(i), z(i)/)  ! coordinate of core spins
     else if ((rsqt > Rc**2) .and. (rsqt <= R**2)) then
        Nsh = Nsh + 1
        tmp_coord_shell(Nsh,1:3) = (/x(i), y(i), z(i)/) ! coordinate of shell spins
     endif
     Nn = Nc + Nsh	
  enddo
  
  ! allocation for core spin and shell spin arrays
  allocate( site_core (Nc ) )
  allocate( site_shell(Nsh) )

  allocate(coord_core (Nc,3) ) ! intermediate allocation for core spins coordinate
  coord_core(1:Nc,1:3) = tmp_coord_core(1:Nc,1:3) ! reduce the size from N to Nc
  do i = 1, Nc
     site_core(i)%r = coord_core(i,1:3)  ! set coordinate information for core spins
  end do

  allocate(coord_shell(Nsh,3)) ! intermediate allocation for shell spins coordinate
  coord_shell(1:Nc,1:3) = tmp_coord_shell(1:Nsh,1:3) ! reduce the size from N to Nsh
  do i = 1, Nsh
     site_shell(i)%r = coord_shell(i,1:3) ! set coordinate information for shell spins
  end do

  !print*, N, Nc, Nsh, Nn, size(site_core), size(site_shell)

  !-----------------------------
  !   nearest neigbor spins 
  !-----------------------------

  ! core
  do i = 1, Nc
  counter = 0
  tmp_list_nearest_spin = 0 ! initializing the list
  do j = 1, Nc
     rsqt  = ( ( coord_core(i,1) - coord_core(j,1) )**2 + &
               ( coord_core(i,2) - coord_core(j,2) )**2 + &
               ( coord_core(i,3) - coord_core(j,3) )**2 )
     if ((j .ne. i).and.(rsqt .eq. a*a)) then
        counter = counter + 1
        tmp_list_nearest_spin(counter) = j
     end if
  enddo ! do j
  allocate(site_core(i)%list_nrcc(counter))
  site_core(i)%list_nrcc = tmp_list_nearest_spin(1:counter) ! set nearest spins list for core spins
  !print*, i, size(site_core(i)%list_nrcc), site_core(i)%list_nrcc(1:size(site_core(i)%list_nrcc))
  enddo ! do i

  ! shell
  do i = 1, Nsh
  counter = 0
  tmp_list_nearest_spin = 0 ! initializing the list
  do j = 1, Nsh
     rsqt  = ( ( coord_shell(i,1) - coord_shell(j,1) )**2 + &
               ( coord_shell(i,2) - coord_shell(j,2) )**2 + &
               ( coord_shell(i,3) - coord_shell(j,3) )**2 )
     if ((j .ne. i).and.(rsqt .eq. a*a)) then
        counter = counter + 1
        tmp_list_nearest_spin(counter) = j
     endif
  enddo ! do j
  allocate(site_shell(i)%list_nrss(counter))
  site_shell(i)%list_nrss = tmp_list_nearest_spin(1:counter) ! set nearest spins list for shell spins
  ! print*, i, counter, size(site_shell(i)%list_nrss), ":", site_shell(i)%list_nrss(1:size(site_shell(i)%list_nrss))
  enddo ! do i

  ! core-shell
  do i = 1, Nc
  counter = 0
  tmp_list_nearest_spin = 0
  do j = 1, Nsh
     rsqt  = ( ( coord_core(i,1) - coord_shell(j,1) )**2 + &
               ( coord_core(i,2) - coord_shell(j,2) )**2 + &
               ( coord_core(i,3) - coord_shell(j,3) )**2 )
     if ((j .ne. i).and.(rsqt .eq. a*a)) then
        counter = counter + 1
        tmp_list_nearest_spin(counter) = j
     endif
  enddo ! do j
  allocate(site_core(i)%list_nrcs(counter))
  site_core(i)%list_nrcs = tmp_list_nearest_spin(1:counter) ! set nearest shell spins list for core spin (at interface)
enddo ! do i

  ! shell-core
  do i = 1, Nsh
  counter = 0
  tmp_list_nearest_spin = 0
  do j = 1, Nc
     rsqt  = ( ( coord_shell(i,1) - coord_core(j,1) )**2 + &
               ( coord_shell(i,2) - coord_core(j,2) )**2 + &
               ( coord_shell(i,3) - coord_core(j,3) )**2 )
     if ((j .ne. i).and.(rsqt .eq. a*a)) then
        counter = counter + 1
        tmp_list_nearest_spin(counter) = j
     endif
  enddo ! do j
  allocate(site_shell(i)%list_nrsc(counter))
  site_shell(i)%list_nrsc = tmp_list_nearest_spin(1:counter) ! set nearest core spins list for shell spin (at interface)
  !print*, i, size(site_shell(i)%list_nrsc), site_shell(i)%list_nrsc( 1 : size(site_shell(i)%list_nrsc) )
  
  enddo ! do i

  !-----------------------------
  !   vector magnetic moment
  !----------------------------

  allocate(theta_core (Nc ), phi_core (Nc ))
  allocate(theta_shell(Nsh), phi_shell(Nsh))

  ! module of vector
  m = 0.5d+00

  ! core
  do i = 1, Nc
     CALL RANDOM_NUMBER(rdn)
     theta_core(i) = 1.00d+00 * rdn * pi
     phi_core  (i) = 2.00d+00 * rdn * pi
     site_core(i)%s           = m
     site_core(i)%theta       = theta_core(i)  ! randomly initialize
     site_core(i)%phi         = phi_core  (i)  ! randomly initialize
     site_core(i)%theta_trial = theta_core(i)  ! randomly initialize
     site_core(i)%phi_trial   = phi_core  (i)  ! randomly initialize
  end do

  ! shell
  do i = 1, Nsh
     CALL RANDOM_NUMBER(rdn)
     theta_shell(i) = 1.00d+00 * rdn * pi
     phi_shell  (i) = 2.00d+00 * rdn * pi
     site_shell(i)%s           = m
     site_shell(i)%theta       = theta_shell(i) ! randomly initialize
     site_shell(i)%phi         = phi_shell  (i) ! randomly initialize
     site_shell(i)%theta_trial = theta_shell(i) ! randomly initialize
     site_shell(i)%phi_trial   = phi_shell  (i) ! randomly initialize
  end do

  ! deallocate intermediate dynamic arrays
  deallocate( x, y, z,  &
              coord_core, coord_shell, & 
              tmp_coord_core, tmp_coord_shell, & 
              theta_core , phi_core, & 
              theta_shell, phi_shell )

end subroutine set_info_site
