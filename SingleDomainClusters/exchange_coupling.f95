!----------------------------------------------------------------------
! @brief : calculate the exchange coupling in core/shell nanoparticles
! @author: Tran Nguyen Lan
! @date  : ...
!---------------------------------------------------------------------

subroutine exchange_coupling (Jxc, Exc, i, itrial, itype)

  ! Exc = Jxc ( mx1 * mx2 + my1 * my2 + mz1 * mz2)
  ! The memeber list_nrcc, list_nrss, list_nrcs, and list_nrsc of type site_core and site_shell 
  !     store nearest spins of given spin i

  use eb_module, only : site_core, site_shell
  implicit none
  
  integer,      intent(in)  :: i ! given spin
  integer,      intent(in)  :: itype  ! 1 : core-core, 2 : shell-shell, 3 : core-shell, 4 : shell-core 
  integer,      intent(in)  :: itrial ! 1 : no trial , 2 : trial 
  real(kind=8), intent(in)  :: Jxc
  real(kind=8), intent(out) :: Exc

  integer      :: j, k
  real(kind=8) :: mx1, my1, mz1 ! components of vector magnetic moment
  real(kind=8) :: mx2, my2, mz2 ! components of vector magnetic moment

 if(itrial < 1 .or. itrial > 2) then
     stop "itrial should be 1 or 2"
  end if
   
  if(itype < 1 .or. itype > 4) then
     stop " itype should be 1, 2, 3, or 4"
  end if

  ! core-core
  if (itype == 1) then 

     ! components of vector magnetic moment of given spin
     if (itrial == 1) then
        mx1 = site_core(i)%s * sin( site_core(i)%theta ) * cos( site_core(i)%phi )
        my1 = site_core(i)%s * sin( site_core(i)%theta ) * sin( site_core(i)%phi )
        mz1 = site_core(i)%s * cos( site_core(i)%theta )
     else if (itrial == 2) then
        mx1 = site_core(i)%s * sin( site_core(i)%theta_trial ) * cos( site_core(i)%phi_trial )
        my1 = site_core(i)%s * sin( site_core(i)%theta_trial ) * sin( site_core(i)%phi_trial )
        mz1 = site_core(i)%s * cos( site_core(i)%theta_trial )
     end if
        
     Exc = 0
     do j = 1, size(site_core(i)%list_nrcc)
        k   = site_core(i)%list_nrcc(j)
        mx2 = site_core(k)%s * sin( site_core(k)%theta ) * cos( site_core(k)%phi )
        my2 = site_core(k)%s * sin( site_core(k)%theta ) * sin( site_core(k)%phi )
        mz2 = site_core(k)%s * cos( site_core(k)%theta )
        Exc = Exc - Jxc * ( mx1 * mx2 + my1 * my2 + mz1 * mz2 )
     end do
     
  ! shell-shell
  else if (itype == 2) then 

     ! components of vector magnetic moment of given spin
     if (itrial == 1) then
        mx1 = site_shell(i)%s * sin( site_shell(i)%theta ) * cos( site_shell(i)%phi )
        my1 = site_shell(i)%s * sin( site_shell(i)%theta ) * sin( site_shell(i)%phi )
        mz1 = site_shell(i)%s * cos( site_shell(i)%theta )
     else if (itrial == 2) then
        mx1 = site_shell(i)%s * sin( site_shell(i)%theta_trial ) * cos( site_shell(i)%phi_trial )
        my1 = site_shell(i)%s * sin( site_shell(i)%theta_trial ) * sin( site_shell(i)%phi_trial )
        mz1 = site_shell(i)%s * cos( site_shell(i)%theta_trial )
     end if
     
     Exc = 0
     do j = 1, size(site_shell(i)%list_nrss)
        k   = site_shell(i)%list_nrss(j)
        mx2 = site_shell(k)%s * sin( site_shell(k)%theta ) * cos( site_shell(k)%phi )
        my2 = site_shell(k)%s * sin( site_shell(k)%theta ) * sin( site_shell(k)%phi )
        mz2 = site_shell(k)%s * cos( site_shell(k)%theta )
        Exc = Exc - Jxc * ( mx1 * mx2 + my1 * my2 + mz1 * mz2 )
     end do

  ! core-shell
  else if (itype == 3) then 

     ! components of vector magnetic moment of given spin
     if (itrial == 1) then
        mx1 = site_core(i)%s * sin( site_core(i)%theta ) * cos( site_core(i)%phi )
        my1 = site_core(i)%s * sin( site_core(i)%theta ) * sin( site_core(i)%phi )
        mz1 = site_core(i)%s * cos( site_core(i)%theta )
     else if (itrial == 2) then
        mx1 = site_core(i)%s * sin( site_core(i)%theta_trial ) * cos( site_core(i)%phi_trial )
        my1 = site_core(i)%s * sin( site_core(i)%theta_trial ) * sin( site_core(i)%phi_trial )
        mz1 = site_core(i)%s * cos( site_core(i)%theta_trial )
     end if

     Exc = 0
     do j = 1, size(site_core(i)%list_nrcs)
        k   = site_core(i)%list_nrcs(j)
        mx2 = site_shell(k)%s * sin( site_shell(k)%theta ) * cos( site_shell(k)%phi )
        my2 = site_shell(k)%s * sin( site_shell(k)%theta ) * sin( site_shell(k)%phi )
        mz2 = site_shell(k)%s * cos( site_shell(k)%theta )
        Exc = Exc - Jxc * ( mx1 * mx2 + my1 * my2 + mz1 * mz2 )
     end do

  ! shell-core
  else if (itype == 4) then 

     ! components of vector magnetic moment of given spin
     if (itrial == 1) then
        mx1 = site_shell(i)%s * sin( site_shell(i)%theta ) * cos( site_shell(i)%phi )
        my1 = site_shell(i)%s * sin( site_shell(i)%theta ) * sin( site_shell(i)%phi )
        mz1 = site_shell(i)%s * cos( site_shell(i)%theta )
     else if (itrial == 2) then
        mx1 = site_shell(i)%s * sin( site_shell(i)%theta_trial ) * cos( site_shell(i)%phi_trial )
        my1 = site_shell(i)%s * sin( site_shell(i)%theta_trial ) * sin( site_shell(i)%phi_trial )
        mz1 = site_shell(i)%s * cos( site_shell(i)%theta_trial )
     end if

     Exc = 0
     do j = 1, size(site_shell(i)%list_nrsc)
        k   = site_shell(i)%list_nrsc(j)
        mx2 = site_core(k)%s * sin( site_core(k)%theta ) * cos( site_core(k)%phi )
        my2 = site_core(k)%s * sin( site_core(k)%theta ) * sin( site_core(k)%phi )
        mz2 = site_core(k)%s * cos( site_core(k)%theta )
        Exc = Exc - Jxc * ( mx1 * mx2 + my1 * my2 + mz1 * mz2 )
     end do

  end if ! itype
  
end subroutine exchange_coupling
