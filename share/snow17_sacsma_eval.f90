subroutine snow17_sacsma_eval(a, obj_val)
  use nrtype
  use constants, only: sec_day, sec_hour
  use gauge_calib, only: calc_pet_pt, sfc_pressure, &
                         calc_rmse, calc_mse, calc_nse, calc_kge,&
                         spin_up_first_year
  use snow17_sac

  implicit none
!input variables
  real(sp), dimension(30), intent(inout) 	:: a

!output variables
  real(dp), intent(out)            		:: obj_val


!local variables
  integer :: i,h,k,m,ntau,end_pt,cnt,ll

  real    :: dtuh

  real(dp) :: obj_val_tmp

!sac-sma state variables
  real(dp)              :: uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc


!single precision state variables
  real(sp)		:: uztwc_sp,uzfwc_sp,lztwc_sp
  real(sp)		:: lzfsc_sp,lzfpc_sp,adimc_sp

!single precision forcing variables
  real(sp)		:: tair_sp
  real(sp)		:: precip_sp
  real(sp)		:: pet_sp

!sac-sma output variables
!need to get back to allocatable
  real, dimension(36500)    :: qs,qg,eta,tci,route_tci

  real(dp), dimension(36500)    :: tci_dp,route_tci_dp,streamflow_dp

!snow-17 surface pressure
  real(dp) :: pa

!snow-17 output variables
!need to get back to allocatable
  real(sp), dimension(36500)    :: snowh, sneqv, snow, raim_snow17 	!output variables

!snow-17 carry over variables
  real(sp) :: tprev				!carry over variable
  real(sp), dimension(19)  :: cs		!carry over variable

!unit hydrograph
  real(sp),dimension(1000)		:: unit_hydro  !array for unit hydrograph call (duamel)

  real(dp)			:: spinup_crit	!spin up criteria

!!!!!!!!!!!!
!
!    code
!
!!!!!!!!!!!!

!spin sac-sma up using first year repeated
  spinup_crit = 0.1_dp

!  print *,'spin up'
  call spin_up_first_year(a,spinup_crit,uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
!  print *,'done spinup'
!  print *,uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc

!set single precision state variables
  uztwc_sp = real(uztwc,kind(sp))
  uzfwc_sp = real(uzfwc,kind(sp))
  lztwc_sp = real(lztwc,kind(sp))
  lzfsc_sp = real(lzfsc,kind(sp))
  lzfpc_sp = real(lzfpc,kind(sp))
  adimc_sp = real(adimc,kind(sp))

!reset snow-17 carryover variables
  tprev = 0.0
  cs = 0.0

!set sac output variables to zero
  qs = 0.0
  qg = 0.0
  tci = 0.0
  eta = 0.0
  unit_hydro = 0.0
  route_tci = 0.0

!get snow17 surface pressure
  call sfc_pressure(elev,pa)

! print *,'sfc_pres',pa,elev
!need to set non-optimized parameters here again due to way sce passes things after first iteration...
  a(21) = nmf(1)
  a(22) = tipm(1)
  a(23) = mbase(1)
  a(24) = plwhc(1)
  a(25) = daygm(1)
  a(26) = adimp(1)
  a(27) = pctim(1)
  a(28) = riva
  a(29) = side
  a(30) = rserv
  

  end_pt = sim_length

  call calc_pet_pt(a)

!print *,'pet',pet(100),end_pt

!print *,uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp

  !now run model
  do i = 1,end_pt

  !set single precision inputs
    tair_sp   = real(tair(i),kind(sp))
    precip_sp = real(precip(i),kind(sp))
    pet_sp    = real(pet(i),kind(sp))
      !precip multiplier
!	precip_sp = precip_sp*1.9
    CALL EXSNOW19(int(dt),int(dt/sec_hour),day(i),month(i),year(i),&
      !SNOW17 INPUT AND OUTPUT VARIABLES
			precip_sp,tair_sp,raim_snow17(i),sneqv(i),snow(i),snowh(i),&
      !SNOW17 PARAMETERS
!ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!			    alat,a(1),a(2),a(3),a(4),a(5),a(7),a(8),a(9),&
!			    a(6),a(10),a(11),elev,pa,adc(1),&
!			    alat,a(1),a(2),a(3),a(4),a(5),a(20),a(21),a(22),&
!			    a(6),a(23),a(24),elev,pa,adc(1),&
			real(lat,kind(sp)),a(1),a(2),a(3),a(4),a(5),a(21),a(22),a(23),&
			a(6),a(24),a(25),real(elev,kind(sp)),real(pa,kind(sp)),adc,&
      !SNOW17 CARRYOVER VARIABLES
			cs,tprev) 


 !print *,'here 4',i
    call exsac(1,real(dt),raim_snow17(i),tair_sp,pet_sp,&
      !SAC PARAMETERS
!UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
!REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
!SIDE,RSERV, &
!a(12),a(13),a(18),a(23),a(17),a(25),a(21),a(22),a(14),a(16),a(15),a(20),a(19),a(24),a(26),a(27)
!			a(12),a(13),a(18),a(23),a(17),a(25),a(21),&
!			a(22),a(14),a(16),a(15),a(20),a(19),a(24),&
!			a(26),a(27),&
!			a(7),a(8),a(12),a(26),a(25),a(27),a(15), &
!			a(16),a(9),a(11),a(10),a(14),a(13),a(17),&
!			a(28),a(29), &
		    a(7),a(8),a(12),a(27),a(26),a(28),a(15), &
		    a(16),a(9),a(11),a(10),a(14),a(13),a(17),&
		    a(29),a(30), &
      !SAC State variables
		    uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp,&
      !SAC OUTPUTS
		    qs(i),qg(i),tci(i),eta(i))


  enddo

!print *,'model done'

  dtuh = real(dt/sec_day)

  if (a(18) .le. 0.0 .and. a(19) .le. 0.0) THEN
    k = 0
    m = 1
  else
    k = 1
    m = 1000
  end if
  ntau = 0
!call unit hydrograph routine
  if(a(18) .gt. 0.0) then
    call DUAMEL(tci,1,unit_hydro,a(18),a(19),dtuh,end_pt-1,m,route_tci,k,ntau)
			         !shape,scale
  endif

!print *,'unit hydrograph',a(18)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!      calculate objective function for daily streamflow
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !need to pass kind(dp) to objective function subroutines
  route_tci_dp = real(route_tci, kind(dp))
  streamflow_dp = real(streamflow, kind(dp))
  tci_dp = real(tci, kind(dp))

  if(a(18) .gt. 0.0) then
    if(trim(metric) .eq. "rmse" .or. trim(metric) .eq. "RMSE") then
      call calc_rmse(route_tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "mse" .or. trim(metric) .eq. "MSE") then
      call calc_mse(route_tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "nse" .or. trim(metric) .eq. "NSE") then
      call calc_nse(route_tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "kge" .or. trim(metric) .eq. "KGE") then
      call calc_kge(route_tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "fms" .or. trim(metric) .eq. "FMS") then
      call fdc_fms(route_tci_dp,sort_streamflow,order,end_pt,valid,obj_val_tmp)
      obj_val_tmp = abs(obj_val_tmp)
    elseif(trim(metric) .eq. "fhv" .or. trim(metric) .eq. "FHV") then
      call fdc_fhvbias(route_tci_dp,sort_streamflow,order,end_pt,valid,obj_val_tmp)
      obj_val_tmp = abs(obj_val_tmp)
    endif
  else
    if(trim(metric) .eq. "rmse" .or. trim(metric) .eq. "RMSE") then
      call calc_rmse(tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "mse" .or. trim(metric) .eq. "MSE") then
      call calc_mse(tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "nse" .or. trim(metric) .eq. "NSE") then
      call calc_nse(tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "kge" .or. trim(metric) .eq. "KGE") then
      call calc_kge(tci_dp,streamflow_dp,end_pt,valid,obj_val_tmp)
    elseif(trim(metric) .eq. "fms" .or. trim(metric) .eq. "FMS") then
      call fdc_fms(tci_dp,sort_streamflow,order,end_pt,valid,obj_val_tmp)
      obj_val_tmp = abs(obj_val_tmp)
    elseif(trim(metric) .eq. "fhv" .or. trim(metric) .eq. "FHV") then
      call fdc_fhvbias(tci_dp,sort_streamflow,order,end_pt,valid,obj_val_tmp)
      obj_val_tmp = abs(obj_val_tmp)
    endif
  endif

  obj_val = obj_val_tmp
!print *,obj_val

  return
end subroutine