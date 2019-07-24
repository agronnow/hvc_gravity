

!=======================================================================
subroutine solve_cooling_sd93(nH,T2,zsolar,dt,deltaT2,ncell)
!=======================================================================
  implicit none
  integer::ncell
  real(kind=8)::dt
  real(kind=8),dimension(1:ncell)::nH,T2,deltaT2,zsolar

  real(kind=8)::precoeff
  real(kind=8)::cool,heat
  real(kind=8)::cool_prime,heat_prime,wcool
  real(kind=8)::lambda,lambda_prime
  real(kind=8),dimension(1:ncell)::tau,tau_old
  real(kind=8),dimension(1:ncell)::time,time_old,tau_ini
  real(kind=8),dimension(1:ncell)::wmax,time_max
  real(kind=8)::varmax=4d0
  integer::i,i_T2,iter,n,n_active
  integer,dimension(1:ncell)::ind,i_nH
  logical::tau_negative

  ! Initializations
  precoeff=2d0*X/(3d0*kB)
  do i=1,ncell
     tau(i)=T2(i)
     tau_ini(i)=T2(i)
     time_max(i)=dt*precoeff*nH(i)
     time(i)=0d0
     wmax(i)=1d0/time_max(i)
     ind(i)=i
  end do

  ! Check positivity
  tau_negative=.false.
  do i=1,ncell
     if(tau(i)<=0.)tau_negative=.true.
  end do
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling_sd93 :'
     write(*,*)'Initial temperature is negative'
     STOP
  endif

  ! Loop over active cells
  iter=0
  n=ncell
  do while(n>0)

     iter=iter+1
     if (iter > 500) then
        write(*,*) 'Too many iterations in solve_cooling_sd93',iter,n
        do i=1,n
           write(*,*)i,tau(ind(i)),T2(ind(i)),nH(ind(i)),i_nH(ind(i))
        end do
        STOP
     endif

     n_active=0
     do i=1,n
        if (temp >= minCoolingTemp) then
           Zstatus = 0
           Zhi = 0.0
           Zlo = 0.0
           if (Zsolar > 1.0) then
              Zstatus = 2 !Metallicity greater than max table
           else if (Zsolar > 0.3162) then
              Zhi = 1.0
              Zlo = 0.3162
              L_tab_loZ => cooltable_L_05
              L_tab_hiZ => cooltable_L_0
           else if (Zsolar > 0.1) then
              Zhi = 0.3162
              Zlo = 0.1
              L_tab_loZ => cooltable_L_1
              L_tab_hiZ => cooltable_L_05
           else
              Zstatus = 1 !Metallicity lower than min table
           endif
           dZ = Zhi - Zlo

           ! ----------------------------------------------
           !          Table lookup by binary search  
           ! ----------------------------------------------
           klo = 0
           khi = ntab - 1

           if (temp > cooltable_T(khi)) then
              OutOfBounds_hi = OutOfBOunds_hi+1
              !    print (" ! T out of range   %12.6e\n",T);
              !    QUIT_PLUTO(1);
              if (Zstatus == 1) then
                 cool = cooltable_L_1(khi)
              else if (Zstatus == 2) then
                 cool = cooltable_L_0(khi)
              else
                 cool = L_tab_loZ(khi)*(Zhi - Zsolar)/dZ + L_tab_hiZ(khi)*(Zsolar - Zlo)/dZ
              endif
           else if (temp < cooltable_T(klo)) then
              OutOfBounds_low = OutOfBounds_low+1
              !    print (" ! T out of range   %12.6e\n",T);
              !    QUIT_PLUTO(1);
              if (Zstatus == 1) then
                 cool = cooltable_L_1(klo)
              else if (Zstatus == 2) then
                 cool = cooltable_L_0(klo)
              else
                 cool = L_tab_loZ(klo)*(Zhi - Zsolar)/dZ + L_tab_hiZ(klo)*(Zsolar - Zlo)/dZ
              endif
           else
              do while (klo /= (khi - 1))
                 kmid = (klo + khi)/2
                 Tmid = T_tab(kmid)
                 if (temp <= Tmid) then
                    khi = kmid
                 else if (temp > Tmid) then
                    klo = kmid
                 endif
              end do
              dT       = cooltable_T(khi) - cooltable_T(klo)
              if (Zstatus == 1) then
                 cool = L_tab_z1(klo)*(cooltable_T(khi) - temp)/dT + L_tab_z1(khi)*(temp - cooltable_T(klo))/dT
              else if (Zstatus == 2) then
                 cool = L_tab_z0(klo)*(cooltable_T(khi) - temp)/dT + L_tab_z0(khi)*(temp - cooltable_T(klo))/dT
              else
                 L_loZ = L_tab_loZ(klo)*(cooltable_T(khi) - temp)/dT + L_tab_loZ(khi)*(temp - cooltable_T(klo))/dT
                 L_hiZ = L_tab_hiZ(klo)*(cooltable_T(khi) - temp)/dT + L_tab_hiZ(khi)*(temp - cooltable_T(klo))/dT
                 cool  = L_loZ*(Zhi - Zsolar)/dZ + L_hiZ*(Zsolar - Zlo)/dZ;
              endif
              cool_prime = cool/tau(ind(i))
           endif

           ! Total net cooling
           lambda=cool-heat
           lambda_prime=cool_prime-heat_prime

        else
           lambda=0d0
           lambda_prime=0d0
        endif

        wcool=MAX(abs(lambda)/tau(ind(i))*varmax,wmax(ind(i)),-lambda_prime*varmax)

        tau_old(ind(i))=tau(ind(i))
        tau(ind(i))=tau(ind(i))*(1d0+lambda_prime/wcool-lambda/tau(ind(i))/wcool)/(1d0+lambda_prime/wcool)
        time_old(ind(i))=time(ind(i))
        time(ind(i))=time(ind(i))+1d0/wcool

!!$        if(i==1)then
!!$           write(10,'(I5,10(1PE10.3,1X))')iter,tau_old(ind(i)),cool+zzz(ind(i))*metal,heat,lambda
!!$        endif

        if(time(ind(i))<time_max(ind(i)))then
           n_active=n_active+1
           ind(n_active)=ind(i)
        end if

     end do
     n=n_active
  end do
  ! End loop over active cells

  ! Compute exact time solution
  do i=1,ncell
     tau(i)=tau(i)*(time_max(i)-time_old(i))/(time(i)-time_old(i))+tau_old(i)*(time(i)-time_max(i))/(time(i)-time_old(i))
  end do

  ! Check positivity
  tau_negative=.false.
  do i=1,ncell
     if (tau(i)<=0.)tau_negative=.true.
  end do
  if (tau_negative) then
     write(*,*)'ERROR in solve_cooling_sd93 :'
     write(*,*)'Final temperature is negative'
     STOP
  endif

  ! Compute delta T
  do i=1,ncell
     deltaT2(i)=tau(i)-tau_ini(i)
  end do

end subroutine solve_cooling


