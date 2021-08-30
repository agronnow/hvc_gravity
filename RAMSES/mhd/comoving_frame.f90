#define CALC_COMVEL

subroutine calc_center_of_mass
  use amr_commons
  use hydro_commons
  use cooling_module
  use mpi_mod
  implicit none
  integer::ilevel
  integer::ind,idim,ivar,ix,iy,iz,nx_loc,nleaf
  integer::i,icell,igrid,iskip,ngrid,ncache,ierr,ilev,ilun=150
  integer,save::iter = 0
  real(dp)::dx,dx_loc,scale
  real(dp),dimension(1:3)::skip_loc
  real(dp),dimension(1:twotondim,1:3)::xc
  real(dp),dimension(1:nvector,1:ndim),save::xx
  integer ,dimension(1:nvector),save::ind_grid,ind_cell,ind_leaf
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v,scale_m
  real(dp)::localcenterofmassy,localmass,globalcenterofmassy_level,globalmass_level,vol_loc
#ifdef CALC_COMVEL
  real(dp)::localcenterofmassvely,globalcenterofmassvely_level
  real(dp),dimension(1:MAXLEVEL),save::globalcenterofmassvely = 0.0
#endif
  real(dp),dimension(1:MAXLEVEL),save::globalcenterofmassy = 0.0,globalmass = 0.0
  logical,dimension(1:MAXLEVEL),save::UpdateVel = .false.
  real(dp)::ctime=0.0,cvel=0.0,ccurvel=0.0,cdist=0.0,ccom=0.0,cmass=0.0,ptime,pvel,pcom
  real(dp),save::t_prev = 0.0
  real(dp),save::com_prev = 0.0
  real(dp),save::vel_prev = 1d-12
  real(dp),save::frame_velocity = 0.0
  real(dp)::frame_com = 0.0,cur_velocity = 0.0, prev_momy = 0.0
  logical::file_exist
  logical,save::firstcall = .true.
  character(LEN=256)::fileloc,dummyline
#ifndef WITHOUTMPI
  integer::info
#endif

  call units(scale_l,scale_t,scale_d,scale_v,scale_nh,scale_T2)
  scale_m = scale_l**ndim*scale_d


  do ilevel=levelmin,nlevelmax
     ! Mesh size at level ilevel in coarse cell units
     dx=0.5D0**ilevel

     ! Set position of cell centers relative to grid center
     do ind=1,twotondim
        iz=(ind-1)/4
        iy=(ind-1-4*iz)/2
        ix=(ind-1-2*iy-4*iz)
        if(ndim>0)xc(ind,1)=(dble(ix)-0.5D0)*dx
        if(ndim>1)xc(ind,2)=(dble(iy)-0.5D0)*dx
        if(ndim>2)xc(ind,3)=(dble(iz)-0.5D0)*dx
     end do

     ! Local constants
     nx_loc=(icoarse_max-icoarse_min+1)
     skip_loc=(/0.0d0,0.0d0,0.0d0/)
     if(ndim>0)skip_loc(1)=dble(icoarse_min)
     if(ndim>1)skip_loc(2)=dble(jcoarse_min)
     if(ndim>2)skip_loc(3)=dble(kcoarse_min)
     scale=boxlen/dble(nx_loc)
     dx_loc=dx*scale

     localcenterofmassy = 0
     localmass = 0
     globalcenterofmassy_level = 0
     globalmass_level = 0
#ifdef CALC_COMVEL
     localcenterofmassvely = 0
     globalcenterofmassvely_level = 0
#endif
     ncache=active(ilevel)%ngrid
     do igrid=1,ncache,nvector
        ngrid=MIN(nvector,ncache-igrid+1)
        do i=1,ngrid
           ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
        end do
        ! Loop over cells
        do ind=1,twotondim
           iskip=ncoarse+(ind-1)*ngridmax
           do i=1,ngrid
              ind_cell(i)=iskip+ind_grid(i)
           end do

           ! Gather cell centre positions
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=xg(ind_grid(i),idim)+xc(ind,idim)
              end do
           end do
           ! Rescale position from code units to user units
           do idim=1,ndim
              do i=1,ngrid
                 xx(i,idim)=(xx(i,idim)-skip_loc(idim))*scale
              end do
           end do

           ! Gather leaf cells
           nleaf=0
           do i=1,ngrid
              if(son(ind_cell(i))==0)then
                 nleaf=nleaf+1
                 ind_leaf(nleaf)=ind_cell(i)
              end if
           end do
           if(nleaf.eq.0)cycle

           ! Compute mass
           do i=1,nleaf
              !           write(*,*)"time:",t,"imetal:",imetal
              localmass = localmass + uold(ind_leaf(i),imetal+1)
#ifdef CALC_COMVEL
              localcenterofmassvely=localcenterofmassvely + uold(ind_leaf(i),imetal+1)*uold(ind_leaf(i),3)/uold(ind_leaf(i),1)
#endif
              localcenterofmassy=localcenterofmassy + uold(ind_leaf(i),imetal+1)*xx(i,2)
              !  write(*,*) "density: ", uold(ind_leaf(i),1), "trc density: ", uold(ind_leaf(i),imetal+1), " y: ", xx(i,2)!, " Local CoM: ", localcenterofmassy, " Local mass: ", localmass, " Global CoM: ", globalcenterofmassy, " Global mass: ", globalmass
           end do
        end do
     end do

#ifndef WITHOUTMPI
     call MPI_ALLREDUCE(localcenterofmassy, globalcenterofmassy_level, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     call MPI_ALLREDUCE(localmass, globalmass_level, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)

     globalmass(ilevel-levelmin+1) = globalmass_level*dx_loc**3
     globalcenterofmassy(ilevel-levelmin+1) = globalcenterofmassy_level*dx_loc**3

#ifdef CALC_COMVEL
     call MPI_ALLREDUCE(localcenterofmassvely, globalcenterofmassvely_level, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, info)
     globalcenterofmassvely(ilevel-levelmin+1) = globalcenterofmassvely_level*dx_loc**3
#endif
#else
     globalcenterofmassy(ilevel-levelmin+1) = localcenterofmassy_level*dx_loc**3
     globalmass(ilevel-levelmin+1) = localmass_level*dx_loc**3

#ifdef CALC_COMVEL
     globalcenterofmassvely(ilevel-levelmin+1) = localcenterofmassvely_level*dx_loc**3
#endif
#endif
  enddo


  if (firstcall) then
     fileloc=trim(output_dir)//'framevel.dat'
     ilun=150
     inquire(file=fileloc,exist=file_exist)
     if(file_exist) then
        open(ilun, file=fileloc)
        read(ilun,*)dummyline
        do
           ptime=ctime
           pvel=cvel
           pcom=ccom
           read(150,*,end=100)ctime,cvel,ccurvel,cdist,ccom,cmass
           if (ctime .ge. t) exit
        end do
100     close(ilun)
        t_prev = ptime
        frame_velocity = pvel*1.e5/scale_v
        if(frame_dist == 0.0)frame_dist = cdist
        com_prev = pcom
        if (myid==1)write(*,*)"Restart comoving frame with t_prev ",t_prev," frame_velocity ",frame_velocity," com_prev ", com_prev, " frame_dist ",frame_dist
     endif
  endif
  firstcall = .false.

  if (myid==1)write(*,*)"CoM", sum(globalcenterofmassy)/sum(globalmass), " mass", sum(globalmass), " cmass", cmass

!  if (ilevel==levelmin) then
  iter = iter+1
  if (iter==com_iter) then
     frame_com = sum(globalcenterofmassy)/sum(globalmass)
     if (t_prev > 0) then
        ! Calculate CoM velocity
#ifdef CALC_COMVEL              
        cur_velocity = sum(globalcenterofmassvely)/sum(globalmass)
#else
        cur_velocity = (frame_com - com_prev)/(t - t_prev)
#endif
        !if ((vel_prev < 1.e-10) .or. (abs(cur_velocity/vel_prev) < 5.0)) then !Protect against spurious velocity changes
        frame_dist = frame_dist + frame_velocity*(t - t_prev)
        frame_velocity = frame_velocity + cur_velocity

        if (myid==1) then
           write(*,*) "current com: ", frame_com, " prev com: ", com_prev, " current time: ", t, " prev time: ", t_prev, " current vel: ", cur_velocity, " total vel: ", frame_velocity
           ! Output to file
           fileloc=trim(output_dir)//'framevel.dat'
           ilun=150
           inquire(file=fileloc,exist=file_exist)
           if(.not.file_exist) then
              open(ilun, file=fileloc, form='formatted')
              write(ilun, '(A53)') 'time frame_vel cur_vel frame_dist frame_CoM M_c,total (M_sun)'
           else
              open(ilun, file=fileloc, status="old", position="append", action="write", form='formatted')
           endif
           write(ilun,'(6E26.16)') t, frame_velocity*scale_v/1.e5, cur_velocity*scale_v/1.e5, frame_dist, frame_com, sum(globalmass)*scale_m/2e33
           close(ilun)
        endif
        !              do ilev=levelmin,nlevelmax
        !                 UpdateVel(ilev-levelmin+1) = .true.
        !              end do
        t_prev = t
        com_prev = frame_com
        vel_prev = cur_velocity
        !else if (abs(cur_velocity/vel_prev) > 5.0) then
        !   write(*,*) "Warning: new center of mass velocity too high, velocity not updated! curv/oldv=", cur_velocity/vel_prev
        !endif

        do ilevel=levelmin,nlevelmax
           ncache=active(ilevel)%ngrid
           ! Loop over cells to subtract velocity
           do igrid=1,ncache,nvector
              ngrid=MIN(nvector,ncache-igrid+1)
              do i=1,ngrid
                 ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
              end do
              do ind=1,twotondim
                 iskip=ncoarse+(ind-1)*ngridmax
                 do i=1,ngrid
                    ind_cell(i)=iskip+ind_grid(i)
                 end do
                 ! Gather leaf cells
                 nleaf=0
                 do i=1,ngrid
                    if(son(ind_cell(i))==0)then
                       nleaf=nleaf+1
                       ind_leaf(nleaf)=ind_cell(i)
                    end if
                 end do
                 if(nleaf.eq.0)cycle

                 ! Update velocity to CoM frame
                 do i=1,nleaf
                    prev_momy = uold(ind_leaf(i),3)
                    ! Update y-momentum
                    uold(ind_leaf(i),3) = prev_momy - uold(ind_leaf(i),1)*cur_velocity
                    ! Update total energy
                    uold(ind_leaf(i),ndim+2) = uold(ind_leaf(i),ndim+2) - 0.5*prev_momy**2/uold(ind_leaf(i),1) + 0.5*uold(ind_leaf(i),3)**2/uold(ind_leaf(i),1)!  + 0.5*uold(ind_leaf(i),1)*cur_velocity**2 - prev_momy*cur_velocity
                 end do
              end do
           end do
        enddo

     else
        t_prev = t
        com_prev = frame_com
     endif
     globalcenterofmassy = 0.0
#ifdef CALC_COMVEL
     globalcenterofmassvely = 0.0
#endif
     globalmass = 0.0
     iter = 0
  endif

  !write(*,*) "Level: ", ilevel, " cell vol: ", dx_loc**ndim, " Local CoM: ", localcenterofmassy, " Local mass (M_sol): ", localmass*dx_loc**ndim*scale_m/2e33, " Global CoM: ", globalcenterofmassy, " Global mass (M_sol): ", globalmass*dx_loc**ndim*scale_m/2e33
end subroutine calc_center_of_mass


