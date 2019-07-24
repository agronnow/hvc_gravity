subroutine GetMuFromTemperature(T,mu)
!Note: T is assumed here to be actual temperature in Kelvin, NOT T/mu
  use amr_parameters, ONLY: dp
  implicit none
  real(dp)::T,mu
  integer::klo, khi, kmid, iline, readerr
  logical,save::firstcall=.true.
  real(dp)::Tmid, dT

  if (firstcall) then
     write(*,*) "Reading mu table from disk..."
     allocate(mutable_T(7000))
     allocate(mutable_mu(7000))

     ntabmu = 0
     open(151, 'mutable.dat', status='old')
     ReadLoop: do iline=1, 7000
        read (151, *, iostat=readerr) mutable_T(iline), mutable_mu(iline)
        if (readerr /= 0) then
           if (readerr == iostat_end) then
              exit ReadLoop
           else
              write (*, '( / "Error reading mu table: ", I0 )') readerr
              stop
           end if
        end if
        ntabmu = ntabmu + 1
     end do ReadLoop
     close(151)

  if (T > 2.e5) then
     mu = mutable_mu(ntabmu) !Above this temperature mu is approximately constant
  else if (T < 1.e4) then
     mu = mutable_mu(1) !Below this temperature mu is approximately constant
  else
  klo = 1
  khi = ntab

  if (T > mutable_T(khi)) then
    !g_oof++!//print (" ! T out of range   %12.6e\n",T);
    !//    QUIT_PLUTO(1);
    mu = mutable_mu(khi)
  else if (T < mutable_mu(klo)) then
    !g_oof++;//print (" ! T out of range   %12.6e\n",T);
    !//    QUIT_PLUTO(1);
    mu = mutable_mu(klo)
  else
!    ----------------------------------------------
!              Table lookup by binary search
!    ---------------------------------------------- 
    do while (klo /= (khi - 1))
      kmid = (klo + khi)/2
      Tmid = mutable_T(kmid)
      if (T <= Tmid) then
        khi = kmid
      else if (T > Tmid) then
        klo = kmid
      endif
    end do
    dT = mutable_T(khi) - mutable_T(klo)
    mu = mutable_mu(klo)*(mutable_T(khi) - T)/dT + mutable_mu(khi)*(T - mutable_T(klo))/dT
end subroutine GetMuFromTemperature

subroutine GetMuAndTemperature(T2,T,mu)
!Note: T2 is T/mu NOT temperature
  use amr_parameters, ONLY: aexp, dp
  use cooling_module, ONLY: set_rates, cmp_chem_eq
  implicit none
  real(dp)::T2,err_T,T_left,T_right,mu_left,mu_right,T,mu
  integer::niter


 
T_left = T2*0.5 ; T_right = T2*1.3 ; err_T = 1d-4
 
DO 
  IF (niter > 100) THEN
    WRITE(*,*) "T not converging!!!"
    EXIT
  END IF
  GetMuFromTemperature(T_right,mu_right)
  f_right = T2*scale_T2*mu_right - T_right
  GetMuFromTemperature(T_left,mu_left)
  f_left = T2*scale_T2*mu_left - T_left
  d = (T_right - T_left) / (T_right - T_left) * T_right
  IF (ABS(d) < err_T) THEN
    EXIT    
  END IF
  T_left = T_right
  T_right = T_right - d
  niter = niter + 1
END DO
mu = mu_right
T = T2*mu_right

end subroutine GetMuAndTemperature
