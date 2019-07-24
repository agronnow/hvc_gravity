!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine gravana(x,f,dx,ncell)
  use amr_parameters
  use poisson_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp),dimension(1:nvector,1:ndim)::f ! Gravitational acceleration
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the acceleration using analytical models.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! f(i,1:ndim) is the gravitational acceleration in user units.
  !================================================================
  real(dp)::UNIT_KPCPERMYR2,UNIT_CMPERS2_TO_KPCPERTSIM2
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  integer::i

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_T2)

  UNIT_KPCPERMYR2 = 1.0!(scale_t/3.154e13)**2
  UNIT_CMPERS2_TO_KPCPERTSIM2 = scale_t**2/scale_l

   do i=1,ncell
     f(i,1)=0
     f(i,2)=-gravity_params(1)*UNIT_CMPERS2_TO_KPCPERTSIM2!*UNIT_KPCPERMYR2
#if NDIM == 3
     f(i,3)=0
#endif
  end do

end subroutine gravana
!#########################################################
!#########################################################
!#########################################################
!#########################################################
subroutine phi_ana(rr,pp,ngrid)
  use amr_commons
  use poisson_commons
  implicit none
  integer::ngrid
  real(dp),dimension(1:nvector)::rr,pp
  ! -------------------------------------------------------------------
  ! This routine set up boundary conditions for fine levels.
  ! -------------------------------------------------------------------

  integer :: i
  real(dp):: fourpi

  fourpi=4.D0*ACOS(-1.0D0)

#if NDIM==1
  do i=1,ngrid
     pp(i)=multipole(1)*fourpi/2d0*rr(i)
  end do
#endif
#if NDIM==2
  do i=1,ngrid
     pp(i)=multipole(1)*2d0*log(rr(i))
  end do
#endif
#if NDIM==3
  do i=1,ngrid
     pp(i)=-multipole(1)/rr(i)
  end do
#endif
end subroutine phi_ana
