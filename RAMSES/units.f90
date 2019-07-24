subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  implicit none
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
!   Added by TTG, Feb 2016
!  real(dp)::scale_prs
  logical,save::first_call=.true.
! #ifndef WITHOUTMPI
!   include 'mpif.h'
! #endif
!-----------------------------------------------------------------------
! Conversion factors from user units into cgs units
! For gravity runs, make sure that G=1 in user units.
!-----------------------------------------------------------------------

! scale_d converts mass density from user units into g/cc
  scale_d = units_density
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

! scale_t converts time from user units into seconds
  scale_t = units_time
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)
  if(poisson) scale_t = 1.0/sqrt(6.6726d-8 * scale_d)

! scale_l converts distance from user units into cm
  scale_l = units_length
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  scale_v = scale_l / scale_t
! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_t2 = mh/kb * scale_v**2

! scale_nH converts rho in user units into nH in H/cc
  scale_nh = x/mh * scale_d

! scale_prs converts pressure in user units into dyne/cm^2
!  scale_prs = scale_t2 * scale_d

! output unit information; only main process
!NOTE: add pressure unit
  if ((first_call).and.(myid==1)) then
  first_call = .false.
    write(*,*)
    write(*,*)
    write(*,*) "Physical Units:"
    write(*,*) "---------------"
    write(*,'(a32,1pe10.2)') "Unit density [g/cm^3]: ", scale_d
    write(*,'(a32,1pe10.2)') "Unit time [s]: ", scale_t
    write(*,'(a32,1pe10.2)') "Unit length [cm]: ", scale_l
    write(*,'(a32,1pe10.2)') "Unit velocity [cm/s]: ", scale_v
    write(*,'(a32,1pe10.2)') "Unit temperature/mu [K]: ", scale_t2
    write(*,'(a32,1pe10.2)') "Unit hydrogen density [cm^-3]: ", scale_nh
!    write(*,'(a32,1pe10.2)') "Unit pressure [dyne/cm^2]: ", scale_prs
    write(*,*)
    write(*,*)
  end if

end subroutine units
