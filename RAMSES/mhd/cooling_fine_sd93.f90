subroutine cooling_fine(ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
  use mpi_mod
  implicit none
  integer::ilevel
  !-------------------------------------------------------------------
  ! Compute cooling for fine levels
  !-------------------------------------------------------------------
  integer::ncache,i,igrid,ngrid
  integer,dimension(1:nvector),save::ind_grid

  if(numbtot(1,ilevel)==0)return
  if(verbose)write(*,111)ilevel

  ! Operator splitting step for cooling source term
  ! by vector sweeps
  ncache=active(ilevel)%ngrid
  do igrid=1,ncache,nvector
     ngrid=MIN(nvector,ncache-igrid+1)
     do i=1,ngrid
        ind_grid(i)=active(ilevel)%igrid(igrid+i-1)
     end do
     call coolfine1(ind_grid,ngrid,ilevel)
  end do

  if((cooling.and..not.neq_chem).and.ilevel==levelmin.and.cosmo)then
#ifdef grackle
     if(use_grackle==0)then
        if(myid==1)write(*,*)'Computing new cooling table'
        call set_table(dble(aexp))
     endif
#else
     if(myid==1)write(*,*)'Computing new cooling table'
     call set_table(dble(aexp))
#endif
  endif

111 format('   Entering cooling_fine for level',i2)

end subroutine cooling_fine
!###########################################################
!###########################################################
!###########################################################
!###########################################################
subroutine coolfine1(ind_grid,ngrid,ilevel)
  use amr_commons
  use hydro_commons
  use cooling_module
#ifdef grackle
  use grackle_parameters
#endif
#ifdef ATON
  use radiation_commons, ONLY: Erad
#endif
#ifdef RT
  use rt_parameters, only: nGroups, iGroups
  use rt_hydro_commons
  use rt_cooling_module, only: rt_solve_cooling,iIR,rt_isIRtrap &
       ,rt_pressBoost,iIRtrapVar,kappaSc,kappaAbs,a_r,is_kIR_T,rt_vc
#endif
  use mpi_mod
  implicit none
#if defined(grackle) && !defined(WITHOUTMPI)
  integer::info
#endif
  integer::ilevel,ngrid
  integer,dimension(1:nvector)::ind_grid
  !-------------------------------------------------------------------
  !-------------------------------------------------------------------
  integer::i,ind,iskip,idim,nleaf,nx_loc
  real(dp)::scale_nH,scale_T2,scale_l,scale_d,scale_t,scale_v
  real(kind=8)::dtcool,nISM,nCOM,damp_factor,cooling_switch,t_blast
  real(dp)::polytropic_constant=1.
  integer,dimension(1:nvector),save::ind_cell,ind_leaf
  real(kind=8),dimension(1:nvector),save::T2,delta_T2,ekk,erad,emag,temperature,rho,ndens
  real(kind=8),dimension(1:nvector),save::T2min,Zsolar,boost
  real(dp),dimension(1:3)::skip_loc
  real(kind=8)::dx,dx_loc,scale,vol_loc
  !Temperatures below/above which mu=const. Should match the table used for mu(T).
  real(dp),parameter::T_neutral = 10000
  real(dp),parameter::T_ionized = 15000
#ifdef RT
  integer::ii,ig,iNp,il
  real(kind=8),dimension(1:nvector),save:: ekk_new,T2_new
  logical,dimension(1:nvector),save::cooling_on=.true.
  real(dp)::scale_Np,scale_Fp,work,Npc,Npnew, kIR, E_rad, TR
  real(dp),dimension(1:ndim)::Fpnew
  real(dp),dimension(nIons, 1:nvector),save:: xion
  real(dp),dimension(nGroups, 1:nvector),save:: Np, Np_boost=0d0, dNpdt=0d0
  real(dp),dimension(ndim, nGroups, 1:nvector),save:: Fp, Fp_boost=0d0, dFpdt
  real(dp),dimension(ndim, 1:nvector),save:: p_gas, u_gas
  real(kind=8)::f_trap, NIRtot, EIR_trapped, unit_tau, tau, Np2Ep
  real(kind=8)::aexp_loc, f_dust, xHII
  real(dp),dimension(nDim, nDim):: tEdd ! Eddington tensor
  real(dp),dimension(nDim):: flux
#endif
#ifdef grackle
  real(kind=8),dimension(1:nvector),save:: T2_new
#endif
#if NENER>0
  integer::irad
#endif

  ! Mesh spacing in that level
  dx=0.5D0**ilevel
  nx_loc=(icoarse_max-icoarse_min+1)
  skip_loc=(/0.0d0,0.0d0,0.0d0/)
  if(ndim>0)skip_loc(1)=dble(icoarse_min)
  if(ndim>1)skip_loc(2)=dble(jcoarse_min)
  if(ndim>2)skip_loc(3)=dble(kcoarse_min)
  scale=boxlen/dble(nx_loc)
  dx_loc=dx*scale
  vol_loc=dx_loc**ndim

  ! Conversion factor from user units to cgs units
  call units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)

!Adapted from PLUTO C code
!  -------------------------------------------
!        Read tabulated cooling function
!  -------------------------------------------
!Read multiple tables with different metallcities
!IMPORTANT: All tables must have the exact same number of lines
!and the same temperatures

  if (firstcall) then
     write(*,*) "Reading cooling tables from disk..."
     allocate(cooltable_T(7000))
     allocate(cooltable_L_0(7000))

     ntab = 0
     dummy = 0
     open(151, 'cooltable-z0.dat', status='old')
     ReadLoop0: do iline=1, 7000
        read (151, *, iostat=readerr) cooltable_T(iline), cooltable_L_0(iline)
        if (readerr /= 0) then
           if (readerr == iostat_end) then
              exit ReadLoop0
           else
              write (*, '( / "Error reading cooling table: ", I0 )') readerr
              stop
           end if
        end if
        ntab = ntab + 1
     end do ReadLoop0
     close(151)

     allocate(cooltable_L_05(7000))
     open(151, 'cooltable-z05.dat', status='old')
     ReadLoop05: do iline=1, 7000
        read (151, *, iostat=readerr) dummy, cooltable_L_05(iline)
        if (readerr /= 0) then
           if (readerr == iostat_end) then
              exit ReadLoop05
           else
              write (*, '( / "Error reading cooling table: ", I0 )') readerr
              stop
           end if
        end if
     end do ReadLoop05
     close(151)

     allocate(cooltable_L_1(7000))
     open(151, 'cooltable-z05.dat', status='old')
     ReadLoop1: do iline=1, 7000
        read (151, *, iostat=readerr) dummy, cooltable_L_1(iline)
        if (readerr /= 0) then
           if (readerr == iostat_end) then
              exit ReadLoop1
           else
              write (*, '( / "Error reading cooling table: ", I0 )') readerr
              stop
           end if
        end if
     end do ReadLoop1
     close(151)
     firstcall = .false.
  endif


  ! Loop over cells
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

     do i=1,nleaf
        ! ---------------------------------------------
        !   Get mean molecular weight and temperature 
        ! ---------------------------------------------
        rho(i)=MAX(uold(ind_leaf(i),1),smallr)
        ! Compute metallicity in solar units
        if(metal)then
           Zsolar(i)=uold(ind_leaf(i),imetal)/rho/0.02
        else
           Zsolar(i)=z_ave
        endif
        ! Compute thermal pressure
        e=uold(ind_leaf(i),ndim+2)
        ekk(i)=0.0d0
        do idim=1,ndim
           ekk(i)=ekk(i)+0.5*uold(ind_leaf(i),idim+1)**2/rho(i)
        end do
        erad=0.0d0
#if NENER>0
        do irad=0,nener-1
           erad(i)=erad(i)+uold(ind_leaf(i),inener+irad)
        end do
#endif
        emag(i)=0.0d0
#ifdef SOLVERmhd
        do idim=1,ndim
           emag(i)=emag(i)+0.125d0*(uold(ind_leaf(i),idim+ndim+2)+uold(ind_leaf(i),idim+nvar))**2
        end do
#endif
        T2(i)=(gamma-1.0)*(e-ekk(i)-erad(i)-emag(i))*scale_T2/rho(i)
        call GetMuAndTemperature(T2(i),temp,mu)
        ndens(i) = rho(i)/mu
	temperature(i) = temp
     end do

     ! Compute cooling time step in second
     dtcool = dtnew(ilevel)*scale_t

     if(cooling.and..not.neq_chem)then
        call solve_cooling_sd93(ndens(i),temperature(i),Zsolar(i),dtcool,delta_T2,nleaf)
     endif
   end do

     ! Deal with cooling
     if(cooling.or.neq_chem)then
        ! Compute net energy sink
        do i=1,nleaf
           delta_T2(i) = delta_T2(i)*ndens(i)/scale_T2/(gamma-1.0)
        end do
        ! Compute initial fluid internal energy
        do i=1,nleaf
           temperature(i) = temperature(i)*ndens(i)/scale_T2/(gamma-1.0)
        end do
        ! Turn off cooling in blast wave regions
        if(delayed_cooling)then
           do i=1,nleaf
              cooling_switch = uold(ind_leaf(i),idelay)/max(uold(ind_leaf(i),1),smallr)
              if(cooling_switch > 1d-3)then
                 delta_T2(i) = MAX(delta_T2(i),real(0,kind=dp))
              endif
           end do
        endif
     endif

     ! Update fluid internal energy
     if(cooling.or.neq_chem)then
        do i=1,nleaf
           temperature(i) = temperature(i) + delta_T2(i)
        end do
     endif

     ! Update total fluid energy
     if(cooling .or. neq_chem)then
        do i=1,nleaf
           uold(ind_leaf(i),ndim+2) = temperature(i) + ekk(i) + erad(i) + emag(i)
        end do
     endif

     ! Update delayed cooling switch
     if(delayed_cooling)then
        t_blast=t_diss*1d6*(365.*24.*3600.)
        damp_factor=exp(-dtcool/t_blast)
        do i=1,nleaf
           uold(ind_leaf(i),idelay)=max(uold(ind_leaf(i),idelay)*damp_factor,0d0)
        end do
     endif
  end do
  ! End loop over cells

end subroutine coolfine1

#ifdef RT
!************************************************************************
subroutine cmp_Eddington_tensor(Npc,Fp,T_Edd)

! Compute Eddington tensor for given radiation variables
! Npc     => Photon number density times light speed
! Fp     => Photon number flux
! T_Edd  <= Returned Eddington tensor
!------------------------------------------------------------------------
  use amr_commons
  implicit none
  real(dp)::Npc
  real(dp),dimension(1:ndim)::Fp ,u
  real(dp),dimension(1:ndim,1:ndim)::T_Edd
  real(dp)::iterm,oterm,Np_c_sq,Fp_sq,fred_sq,chi
  integer::p,q
!------------------------------------------------------------------------
  if(Npc .le. 0.d0) then
     write(*,*)'negative photon density in cmp_Eddington_tensor. -EXITING-'
     call clean_stop
  endif
  T_Edd(:,:) = 0.d0
  Np_c_sq = Npc**2
  Fp_sq = sum(Fp**2)              !  Sq. photon flux magnitude
  u(:) = 0.d0                           !           Flux unit vector
  if(Fp_sq .gt. 0.d0) u(:) = Fp/sqrt(Fp_sq)
  fred_sq = Fp_sq/Np_c_sq           !      Reduced flux, squared
  chi = max(4.d0-3.d0*fred_sq, 0.d0)   !           Eddington factor
  chi = (3.d0+ 4.d0*fred_sq)/(5.d0 + 2.d0*sqrt(chi))
  iterm = (1.d0-chi)/2.d0               !    Identity term in tensor
  oterm = (3.d0*chi-1.d0)/2.d0          !         Outer product term
  do p = 1, ndim
     do q = 1, ndim
        T_Edd(p,q) = oterm * u(p) * u(q)
     enddo
     T_Edd(p,p) = T_Edd(p,p) + iterm
  enddo

end subroutine cmp_Eddington_tensor
#endif
