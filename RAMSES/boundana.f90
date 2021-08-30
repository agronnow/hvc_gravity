!############################################################
!############################################################
!############################################################
!############################################################
#define P3_EQUILIBRIUM
subroutine boundana(x,xb,u,dx,ibound,ncell)
  use amr_commons, ONLY: t, myid
  use amr_parameters, ONLY: dp,ndim,nvector,ndens_cloud,T_wind,vel_wind,Z_wind,output_dir,prob_debug,dist_init,frame_dist,y0,nwind0
  use hydro_parameters, ONLY: nvar,nener,boundary_var,gamma,imetal
  use cooling_module, ONLY: twopi, kb, mu_hot
  use poisson_parameters, ONLY: gravity_params
  implicit none
  integer ::ibound                        ! Index of boundary region
  integer ::ncell                         ! Number of active cells
  real(dp)::dx                            ! Cell size
#ifdef SOLVERmhd
  real(dp),dimension(1:nvector,1:nvar+ndim)::u ! Conservative variables
#else
  real(dp),dimension(1:nvector,1:nvar)::u ! Conservative variables
#endif
  real(dp),dimension(1:nvector,1:ndim)::x ! Cell center position.
  !real(dp),dimension(1:nvector)::xb ! Cell center y-position at last non-ghost cell.
  real(dp)::xb ! Cell center y-position at last non-ghost cell.
  !================================================================
  ! This routine generates boundary conditions for RAMSES.
  ! Positions are in user units:
  ! x(i,1:3) are in [0,boxlen]**ndim.
  ! U is the conservative variable vector. Conventions are here:
  ! U(i,1): d, U(i,2:ndim+1): d.u,d.v,d.w and U(i,ndim+2): E.
  ! If MHD, then:
  ! U(i,1): d, U(i,2:4): d.u,d.v,d.w, U(i,5): E,
  ! U(i,6:8): Bleft, U(i,nvar+1:nvar+3): Bright
  ! U is in user units.
  ! ibound is the index of the boundary region defined in the namelist.
  !================================================================
  integer::ivar,i,idim,ilun
!  real(dp),dimension(1:nvector,1:nvar+ndim),save::q ! Primitive variables
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_b,scale_prs
  real(dp)::prs,temp,nWind,ekk,emag,erad,rho,e,delta_dist,T2,nH,mu,zprime_l,zprime_r,ydist,U_grav,E_mag
  real(dp)::ctime=0.0,cvel=0.0,ccurvel=0.0,cdist=0.0,ccom=0.0,cmass=0.0,pdist
  logical::file_exist
  character(LEN=256)::fileloc,dummyline
  logical,save::firstcall = .true.

!#ifdef SOLVERmhd
!  do ivar=1,nvar+3
!#else
!  do ivar=1,nvar
!#endif
!     do i=1,ncell
!        u(i,ivar)=boundary_var(ibound,ivar)
!     end do
!  end do

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_T2)
  scale_b = 1d6*dsqrt(2.0*twopi*scale_d)*scale_v !Code magnetic unit to microGauss
  scale_prs = scale_d * scale_v*scale_v

  if ((firstcall) .and. (frame_dist == 0.0)) then
     fileloc=trim(output_dir)//'framevel.dat'
     ilun=150
     inquire(file=fileloc,exist=file_exist)
     if(file_exist) then
        open(ilun, file=fileloc)
        read(ilun,*)dummyline
        do
           pdist = cdist
           read(150,*,end=100)ctime,cvel,ccurvel,cdist,ccom,cmass
           if (ctime .ge. t) exit
        end do
100     close(ilun)
        frame_dist = pdist
        if (myid==1)write(*,*)"Restart boundaries with frame_dist ",frame_dist
     endif
     nH = 1d-3*scale_nH
     call GetMuFromTemperature(T_wind,nH,mu)
     mu_hot = mu
  endif
  firstcall = .false.

  do i=1,ncell
     rho = u(i,1)
     e=u(i,ndim+2)
    ! Compute thermal pressure
     ekk=0.0d0
     do idim=1,ndim
        ekk=ekk+0.5*u(i,idim+1)**2/rho
     end do
     erad=0.0d0
#if NENER>0
     do irad=0,nener-1
        erad=erad+u(i,inener+irad)
     end do
#endif
     emag=0.0d0
     ydist = frame_dist + x(i,2) + dist_init

#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(u(i,idim+ndim+2)+u(i,idim+nvar))**2
     end do

     !Calculate height on staggered grid, B_j-1/2, B_j+1/2
     zprime_l = (ydist - 0.5*dx - 1.5)/4.0
     zprime_r = (ydist + 0.5*dx - 1.5)/4.0
     !Azimuthal field from Sun & Reich (2010) in cylindrical phi coordinate
     u(i,ndim+3) = 2.0/(1.0+zprime_l*zprime_l)/scale_b
     u(i,nvar+1) = 2.0/(1.0+zprime_r*zprime_r)/scale_b
     u(i,ndim+4) = 0.0
     u(i,nvar+2) = 0.0
#if NDIM>2
     u(i,ndim+5) = 0.0
     u(i,nvar+3) = 0.0
#endif

#endif

#ifdef P3_EQUILIBRIUM
#ifdef SOLVERmhd
     E_mag = 0.125d0*(u(i,ndim+3)+u(i,nvar+1))**2 !Magnetic field energy density
#else
     E_mag = 0.0
#endif
     nWind = nwind0*dexp(-gravity_params(1)*mu_hot*scale_d*ydist*scale_l/(kb*T_wind))!nwind0*dexp(-ydist/y0)
     U_grav = -nWind*T_wind/scale_T2!-mu_hot*ndens_wind*gravity_params(1)*y0*scale_l*scale_d/scale_prs !Gravitational potential per unit volume at scale height
!     nWind = nwind0*dexp(-ydist/y0)
!     U_grav = -mu_hot*nWind*gravity_params(1)*y0*scale_l*scale_d/scale_prs !Gravitational potential per unit volume at scale height
     prs = -(U_grav+E_mag) !Pressure from magnetohydrostatic equilibrium: P=-rho*y0*g-B^2/(8*pi)
     mu = mu_hot
#else
     ! Compute T in Kelvin
     T2=(gamma-1.0)*(e-ekk-erad-emag)/rho*scale_T2
     nH = rho*scale_nH

     !call GetMuAndTemperature(T2,nH,temp,mu)
     mu = mu_hot
     temp = T2*mu

     delta_dist = xb
     nWind = (rho/mu)*dexp(gravity_params(1)*mu*scale_d*delta_dist*scale_l/(kb*temp))
     prs = nWind*temp/scale_T2
#endif
     u(i,1) = nWind*mu
     ! Update momentum and energy to match new density and magnetic field at boundary
     ekk=0.0d0
     do idim=1,ndim
        u(i,idim+1) = u(i,idim+1)*u(i,1)/rho
        ekk=ekk+0.5*u(i,idim+1)**2/u(i,1)
     end do
     emag=0.0d0
#ifdef SOLVERmhd
     do idim=1,ndim
        emag=emag+0.125d0*(u(i,idim+ndim+2)+u(i,idim+nvar))**2
     end do
#endif
     u(i,ndim+2) = ekk+erad+emag+prs/(gamma-1.0)

     if ((i==1) .and. (prob_debug)) then
        ! Output to file
        fileloc=trim(output_dir)//'boundary.dat'
        ilun=140
        inquire(file=fileloc,exist=file_exist)
        if(.not.file_exist) then
           open(ilun, file=fileloc, form='formatted')
        else
           open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
        endif
        write(ilun,'(13E26.16)') t, x(i,2), ydist, zprime_l, zprime_r, u(i,1), u(i, 2), u(i, 3), u(i, 4), u(i, 5), rho, prs, mu
        close(ilun)
     endif
     !write(*,*) "y: ", x(i,2), " yb: ", xb, " delta_dist: ", delta_dist, " rho: ", u(i,1), " rhob: ", rho, " prs: ", prs, " prsb: ", (kb*rho*T/mu_hot)/scale_prs, " T: ", T

     u(i,imetal) = Z_wind*0.02*u(i,1)	!metallicity
     u(i,imetal+1) = 0.0	        !tracer
  end do

!  nWind = 1.d-4
!  Twind = 1.8d6
!  mu = 1.34
!  Pwind = ndens_cloud*kb*T_cloud/scale_prs

!  q(1:,1) = ndens_wind*mu_hot	!density
!  q(1:ncell,2) = 0.0	        !x-velocity
!  q(1:ncell,3) = cur_velocity !vel_wind*1.e5/scale_v	!y-velocity (given in km/s)
!#if NDIM == 3
!  q(1:ncell,4) = 0.0    !z-velocity
!#endif
!  q(1:ncell,ndim+2) = Pwind	!pressure
!  q(1:ncell,imetal) = Z_wind	!metallicity
!  q(1:ncell,imetal+1) = 0.0	!tracer

!  ! Convert primitive to conservative variables
!  ! density -> density
!  u(1:ncell,1)=q(1:ncell,1)
!  ! velocity -> momentum
!  u(1:ncell,2)=q(1:ncell,1)*q(1:ncell,2)
!#if NDIM>1
!  u(1:ncell,3)=q(1:ncell,1)*q(1:ncell,3)
!#endif
!#if NDIM>2
!  u(1:ncell,4)=q(1:ncell,1)*q(1:ncell,4)
!#endif
  ! kinetic energy
!  u(1:ncell,ndim+2)=0.0d0
!  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,2)**2
!#if NDIM>1
!  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,3)**2
!#endif
!#if NDIM>2
!  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+0.5*q(1:ncell,1)*q(1:ncell,4)**2
!#endif
!  ! thermal pressure -> total fluid energy
!  u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+q(1:ncell,ndim+2)/(gamma-1.0d0)
!#if NENER>0
!  ! radiative pressure -> radiative energy
!  ! radiative energy -> total fluid energy
!  do ivar=1,nener
!     u(1:ncell,ndim+2+ivar)=q(1:ncell,ndim+2+ivar)/(gamma_rad(ivar)-1.0d0)
!     u(1:ncell,ndim+2)=u(1:ncell,ndim+2)+u(1:ncell,ndim+2+ivar)
!  enddo
!#endif
!#if NVAR>NDIM+2+NENER
!  ! passive scalars
!  do ivar=ndim+3+nener,nvar
!     u(1:ncell,ivar)=q(1:ncell,1)*q(1:ncell,ivar)
!  end do
!#endif

end subroutine boundana
