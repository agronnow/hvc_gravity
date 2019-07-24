!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine rho_ana(x,d,dx,ncell)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Number of cells
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector)       ::d ! Density
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates analytical Poisson source term.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! d(i) is the density field in user units.
  !================================================================
  integer::i
  real(dp)::r,rx,ry,rho0,R_s,xmass,ymass,zmass,dnfw,dmax

  !emass=2.*boxlen*0.5d0**nlevelmax
  !dmax=1.0/emass/(1.0+emass)**2
  dmax=rho0/((0.001/R_s)*(1+0.001/R_s)**2)
  rho0 =gravity_params(1)
  R_s  =gravity_params(2)
  xmass=gravity_params(3)*boxlen
  ymass=gravity_params(4)*boxlen
!  zmass=gravity_params(5)

  do i=1,ncell
     rx=x(i,1)-xmass
     ry=x(i,2)-ymass
     !rz=x(i,3)-zmass
     r=sqrt(rx**2+ry**2)!+rz**2)
     dnfw=rho0/((r/R_s)*(1+r/R_s)**2)
     d(i)=min(dnfw,dmax)
  end do

end subroutine rho_ana
