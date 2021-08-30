!================================================================
!================================================================
!================================================================
!================================================================
#define P3_EQUILIBRIUM
subroutine condinit(x,u,dx,nn)
  use amr_parameters
  use hydro_parameters
  use poisson_parameters
  use cooling_module, ONLY: twopi, kb, mh, mu_hot
  implicit none
  integer ::nn                            ! Number of cells
  real(dp)::dx                            ! Cell size
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+ndim)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:nvar+ndim),save::q   ! Primitive variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
  real(dp),dimension(1:nvector,1:nvar),save::q   ! Primitive variables
#endif
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine generates initial conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! Q is the primitive variable vector. Conventions are here:
  ! Q(i,1): d, Q(i,2:ndim+1):u,v,w and Q(i,ndim+2): P.
  ! If nvar >= ndim+3, remaining variables are treated as passive
  ! scalars in the hydro solver.
  ! U(:,:) and Q(:,:) are in user units.
  !================================================================
#if NENER>0 || NVAR>NDIM+2+NENER
  integer::ivar
#endif

  integer::i
  real(dp)::currad,width,ZRad,dist,zprime_l,zprime_r
  real(dp)::xc,yc,zc,P_wind,ndens,ndens_wind
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_b,scale_prs
  real(dp)::mu,temp,nH,U_grav,E_mag

  ! Call built-in initial condition generator
  !call region_condinit(x,q,dx,nn)

  ! Add here, if you wish, some user-defined initial conditions

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_T2)
  scale_prs = scale_v*scale_v * scale_d
  scale_b = 1d6*dsqrt(2.0*twopi*scale_d)*scale_v !code magnetic unit to microGauss

  xc = x1_c*boxlen
  yc = x2_c*boxlen
  zc = x3_c*boxlen

  width = rad_cloud/density_steepness

  ! Calculate mean molecular weight in the high temperature (fully ionised) limit
  nH = 1d-3*scale_nH
  call GetMuFromTemperature(T_wind,nH,mu)
  mu_hot = mu

  do i=1,nn
#if NDIM==3
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2 + (x(i,3)-zc)**2)
#else
    currad = dsqrt((x(i,1)-xc)**2 + (x(i,2)-yc)**2)
#endif

    dist = x(i,2) + dist_init

#ifdef SOLVERmhd
    !Calculate heights on staggered grid, B_j-1/2, B_j+1/2
    zprime_l = (dist - 0.5*dx - 1.5)/4.0
    zprime_r = (dist + 0.5*dx - 1.5)/4.0
    !Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
    q(i,ndim+3) = 2.0*(4.0/bfield_r0)*dexp(-(bfield_r0-4.0)/4.0)/(1.0+zprime_l*zprime_l)/scale_b
    q(i,nvar+1) = 2.0*(4.0/bfield_r0)*dexp(-(bfield_r0-4.0)/4.0)/(1.0+zprime_r*zprime_r)/scale_b
    q(i,ndim+4) = 0.0
    q(i,nvar+2) = 0.0
#if NDIM>2
    q(i,ndim+5) = 0.0
    q(i,nvar+3) = 0.0
#endif
#endif

#ifdef P3_EQUILIBRIUM
#ifdef SOLVERmhd
    E_mag = 0.125d0*(q(i,ndim+3)+q(i,nvar+1))**2  !Magnetic field energy density
#else
    E_mag = 0.0
#endif
    ndens_wind = nwind0*dexp(-gravity_params(1)*mu_hot*scale_d*dist*scale_l/(kb*T_wind))!nwind0*dexp(-dist/y0)
    U_grav = -ndens_wind*T_wind/scale_T2!-mu_hot*ndens_wind*gravity_params(1)*y0*scale_l*scale_d/scale_prs !Gravitational potential per unit volume at scale height
    P_wind = -(U_grav+E_mag) !Pressure from magnetohydrostatic equilibrium: P=-rho*y0*g-B^2/(8*pi)
!    write(*,*)"n", ndens_wind, "P", P_wind, "T", P_wind/ndens_wind*scale_T2, "E_mag", E_mag, "B", 0.5*(q(i,ndim+3)+q(i,nvar+1)), "U_grav", U_grav
#else
    ndens_wind = nwind0*dexp(-gravity_params(1)*mu_hot*scale_d*dist*scale_l/(kb*T_wind))
    P_wind = ndens_wind*T_wind/scale_T2
#endif

    ndens = ndens_wind + (ndens_cloud - ndens_wind)*0.5*(1.0-tanh(1.0+(currad-(rad_cloud+width))/width))
    ZRad = (rad_cloud/density_steepness)*atanh((ndens_cloud-3*ndens_wind)/(ndens_cloud-ndens_wind)) + rad_cloud

    temp = P_wind/ndens*scale_T2
    nH = ndens*scale_nH

    call GetMuFromTemperature(temp,nH,mu)

!write(*,*) 'T: ', temp, " mu: ", mu, " n: ", ndens
    q(i,1) = ndens*mu
    q(i,2) = 0.0      !x-velocity
#if NDIM==3
    q(i,4) = 0.0      !z-velocity
#endif
    q(i,ndim+2) = P_wind
    !metallicity converted from relative to solar to absolute assuming Z_sol=0.02 as hardcoded in other parts of RAMSES
    if (currad < ZRad) then
       q(i,imetal) = Z_cloud*0.02
    else
       q(i,imetal) = Z_wind*0.02
    endif
    if (currad < rad_cloud*(1.0 + 3.0/density_steepness)) then
       q(i,3) = 0.0      !y-velocity
       if (currad < rad_cloud) then
          q(i,imetal+1) = 1.0 !tracer
       else
          q(i,imetal+1) = 0.0 !tracer
       endif
    else
       q(i,3) = vel_wind*1.e5/scale_v !y-velocity (given in km/s)
       q(i,imetal+1) = 0.0 !tracer
    endif
    !write(*,*) "i ", i, " r ", currad, " rho: ", q(i,1), " Z: ", q(i,imetal), " vy: ", q(i,3), " P: ", q(i,ndim+2)*scale_prs/(kb*ndens), " T ", T, " mu: ", mu

  enddo

  ! Convert primitive to conservative variables
  ! density -> density
  u(1:nn,1)=q(1:nn,1)
  ! velocity -> momentum
  u(1:nn,2)=q(1:nn,1)*q(1:nn,2)
#if NDIM>1
  u(1:nn,3)=q(1:nn,1)*q(1:nn,3)
#endif
#if NDIM>2
  u(1:nn,4)=q(1:nn,1)*q(1:nn,4)
#endif
  ! kinetic energy
  u(1:nn,ndim+2)=0.0d0
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,2)**2
#if NDIM>1
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,3)**2
#endif
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.5*q(1:nn,1)*q(1:nn,4)**2
#endif
  ! thermal pressure -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+q(1:nn,ndim+2)/(gamma-1.0d0)
#ifdef SOLVERmhd
  ! magnetic energy -> total fluid energy
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.125d0*(q(1:nn,ndim+3)+q(1:nn,nvar+1))**2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.125d0*(q(1:nn,ndim+4)+q(1:nn,nvar+2))**2
#if NDIM>2
  u(1:nn,ndim+2)=u(1:nn,ndim+2)+0.125d0*(q(1:nn,ndim+5)+q(1:nn,nvar+3))**2
#endif
  ! Magnetic field components
  u(1:nn,ndim+3:2*ndim+2)=q(1:nn,ndim+3:2*ndim+2)
  u(1:nn,nvar+1:nvar+ndim)=q(1:nn,nvar+1:nvar+ndim)
#endif
#if NENER>0
  ! radiative pressure -> radiative energy
  ! radiative energy -> total fluid energy
  do irad=1,nener
     u(1:nn,2*ndim+2+irad)=q(1:nn,2*ndim+2+irad)/(gamma_rad(irad)-1.0d0)
     u(1:nn,ndim+2)=u(1:nn,ndim+2)+u(1:nn,2*ndim+2+irad)
  enddo
#endif
#if NVAR>2*NDIM+2+NENER
  ! passive scalars
  do ivar=2*ndim+3+nener,nvar
     u(1:nn,ivar)=q(1:nn,1)*q(1:nn,ivar)
  end do
#endif

end subroutine condinit

!================================================================
!================================================================
!================================================================
!================================================================
subroutine velana(x,v,dx,t,ncell)
  use amr_parameters
  use hydro_parameters
  implicit none
  integer ::ncell                         ! Size of input arrays
  real(dp)::dx                            ! Cell size
  real(dp)::t                             ! Current time
  real(dp),dimension(1:nvector,1:3)::v    ! Velocity field
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !================================================================
  ! This routine computes the user defined velocity fields.
  ! x(i,1:ndim) are cell center position in [0,boxlen] (user units).
  ! v(i,1:3) is the imposed 3-velocity in user units.
  !================================================================
  integer::i
  real(dp)::xx,yy=0.,zz=0.,vx,vy,vz,aa,twopi
!!$  real(dp)::rr,tt,omega
  
  ! Add here, if you wish, some user-defined initial conditions
  aa=1.0+0.*t
  twopi=2d0*ACOS(-1d0)
  do i=1,ncell

     xx=x(i,1)+0.*dx
#if NDIM > 1
     yy=x(i,2)
#endif
#if NDIM > 2
     zz=x(i,3)
#endif
     ! ABC
     vx=aa*(cos(twopi*yy)+sin(twopi*zz))
     vy=aa*(sin(twopi*xx)+cos(twopi*zz))
     vz=aa*(cos(twopi*xx)+sin(twopi*yy))

!!$     ! 1D advection test
!!$     vx=1.0_dp
!!$     vy=0.0_dp
!!$     vz=0.0_dp

!!$     ! Ponomarenko
!!$     xx=xx-boxlen/2.0
!!$     yy=yy-boxlen/2.0
!!$     rr=sqrt(xx**2+yy**2)
!!$     if(yy>0)then
!!$        tt=acos(xx/rr)
!!$     else
!!$        tt=-acos(xx/rr)+twopi
!!$     endif
!!$     if(rr<1.0)then
!!$        omega=0.609711
!!$        vz=0.792624
!!$     else
!!$        omega=0.0
!!$        vz=0.0
!!$     endif
!!$     vx=-sin(tt)*rr*omega
!!$     vy=+cos(tt)*rr*omega

     v(i,1)=vx
#if NDIM > 1
     v(i,2)=vy
#endif
#if NDIM > 2
     v(i,3)=vz
#endif
  end do


end subroutine velana

subroutine GetMuFromTemperature(T,nH,mu)
!Note: T is assumed here to be actual temperature in Kelvin, NOT T/mu
  use amr_parameters, ONLY: dp, aexp
  use cooling_module, ONLY: set_rates, cmp_chem_eq
  implicit none
  real(dp)::T,nH,mu
  real(dp)::mu_old,err_mu,mu_left,mu_right,n_TOT
  real(dp),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(dp),dimension(1:6) :: n_spec
  integer::niter

    call set_rates(t_rad_spec,h_rad_spec,aexp)

    ! Iteration to find mu
    err_mu=1.
    mu_left=0.5
    mu_right=1.3
    niter=0
    do while (err_mu > 1.d-4 .and. niter <= 50)
       mu_old=0.5*(mu_left+mu_right)
       !T = T2*mu_old
       call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)
       err_mu = (mu-mu_old)/mu_old
       if(err_mu>0.)then
          mu_left =0.5*(mu_left+mu_right)
          mu_right=mu_right
       else
          mu_left =mu_left
          mu_right=0.5*(mu_left+mu_right)
       end if
       err_mu=ABS(err_mu)
       niter=niter+1
    end do
    if (niter > 50) then
       write(*,*) 'ERROR in calculation of mu : too many iterations.'
       STOP
    endif
end subroutine GetMuFromTemperature

subroutine GetMuAndTemperature(T2,nH,T,mu)
!Note: T2 is T/mu NOT temperature
  use amr_parameters, ONLY: dp, aexp
  use cooling_module, ONLY: set_rates, cmp_chem_eq
  implicit none
  real(dp)::T2,mu_old,err_mu,mu_left,mu_right,n_TOT,nH,T,mu
  real(dp),dimension(1:3) :: t_rad_spec,h_rad_spec
  real(dp),dimension(1:6) :: n_spec
  integer::niter

    call set_rates(t_rad_spec,h_rad_spec,aexp)

    ! Iteration to find mu
    err_mu=1.
    mu_left=0.5
    mu_right=1.3
    niter=0
    do while (err_mu > 1.d-4 .and. niter <= 50)
       mu_old=0.5*(mu_left+mu_right)
       T = T2*mu_old
       call cmp_chem_eq(T,nH,t_rad_spec,n_spec,n_TOT,mu)
       err_mu = (mu-mu_old)/mu_old
       if(err_mu>0.)then
          mu_left =0.5*(mu_left+mu_right)
          mu_right=mu_right
       else
          mu_left =mu_left
          mu_right=0.5*(mu_left+mu_right)
       end if
       err_mu=ABS(err_mu)
       niter=niter+1
    end do
    if (niter > 50) then
       write(*,*) 'ERROR in calculation of mu : too many iterations.'
       STOP
    endif
end subroutine GetMuAndTemperature

