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
  real(sp)		:: fuse_sp

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

!total flow variables
  real(sp),dimension(36500)	:: total_tci,total_qg

!hru related variables
  real(dp)			:: tot_area

!hru parameter array
  real(sp),dimension(41)	:: params_in

!!!!!!!!!!!!
!
!    code
!
!!!!!!!!!!!!

!spin sac-sma up using first year repeated
  spinup_crit = 0.1_dp

!initialize total basin modeled streamflow variables
  total_tci = 0.0
  total_qg  = 0.0

  tot_area = sum(hru_area)
!print *,tot_area
!print *,hru_area
  do h = 1,num_hru
  !  print *,'spin up'
  !  print *,hru_tair(h,12),hru_precip(h,12),hru_pet(h,12),fuse_raim(h,12)

  !setup second parameter array to pass to spin_up 

  !setup parameters using init_params and multipliers in a vector

    !Snow-17 parameters
    params_in(1) = a(1)*init_params(h,1)
!    if(params_in(1) < blp(1)) params_in(1) = blp(1)
!    if(params_in(1) > bup(1)) params_in(1) = bup(1)

    params_in(2) = a(2)*init_params(h,2)
!    if(params_in(2) < blp(2)) params_in(2) = blp(2)
!    if(params_in(2) > bup(2)) params_in(2) = bup(2)

    params_in(3) = a(3)*init_params(h,3)
!    if(params_in(3) < blp(3)) params_in(3) = blp(3)
!    if(params_in(3) > bup(3)) params_in(3) = bup(3)

    params_in(4)  = a(4)*init_params(h,4)
!    if(params_in(4) < blp(4)) params_in(4) = blp(4)
!    if(params_in(4) > bup(4)) params_in(4) = bup(4)

    params_in(5)  = a(5)*init_params(h,5)
!    if(params_in(5) < blp(5)) params_in(5) = blp(5)
!    if(params_in(5) > bup(5)) params_in(5) = bup(5)

    params_in(6)  = a(6)*init_params(h,6)
!    if(params_in(6) < blp(6)) params_in(6) = blp(6)
!    if(params_in(6) > bup(6)) params_in(6) = bup(6)


    !SAC-SMA parameters
    params_in(7) = a(7)*init_params(h,23)
!    if(params_in(7) < blp(7)) params_in(7) = blp(7)
!    if(params_in(7) > bup(7)) params_in(7) = bup(7)

    params_in(8) = a(8)*init_params(h,24)
!    if(params_in(8) < blp(8)) params_in(8) = blp(8)
!    if(params_in(8) > bup(8)) params_in(8) = bup(8)

    params_in(9) = a(9)*init_params(h,25)
!    if(params_in(9) < blp(9)) params_in(9) = blp(9)
!    if(params_in(9) > bup(9)) params_in(9) = bup(9)

    params_in(10) = a(10)*init_params(h,26)
!    if(params_in(10) < blp(10)) params_in(10) = blp(10)
!    if(params_in(10) > bup(10)) params_in(10) = bup(10)

    params_in(11) = a(11)*init_params(h,27)
!    if(params_in(11) < blp(11)) params_in(11) = blp(11)
!    if(params_in(11) > bup(11)) params_in(11) = bup(11)

    params_in(12) = a(12)*init_params(h,29)
!    if(params_in(12) < blp(12)) params_in(12) = blp(12)
!    if(params_in(12) > bup(12)) params_in(12) = bup(12)

    params_in(13) = a(13)*init_params(h,30)
!    if(params_in(13) < blp(13)) params_in(13) = blp(13)
!    if(params_in(13) > bup(13)) params_in(13) = bup(13)

    params_in(14) = a(14)*init_params(h,31)
!    if(params_in(14) < blp(14)) params_in(14) = blp(14)
!    if(params_in(14) > bup(14)) params_in(14) = bup(14)

    params_in(15) = a(15)*init_params(h,32)
!    if(params_in(15) < blp(15)) params_in(15) = blp(15)
!    if(params_in(15) > bup(15)) params_in(15) = bup(15)

    params_in(16) = a(16)*init_params(h,33)
!    if(params_in(16) < blp(16)) params_in(16) = blp(16)
!    if(params_in(16) > bup(16)) params_in(16) = bup(16)

    params_in(17) = a(17)*init_params(h,35)
!    if(params_in(17) < blp(17)) params_in(17) = blp(17)
!    if(params_in(17) > bup(17)) params_in(17) = bup(17)

    params_in(21) = init_params(h,7)
    params_in(22) = init_params(h,8)
    params_in(23) = init_params(h,9)
    params_in(24) = init_params(h,10)
    params_in(25) = init_params(h,11)
    params_in(26) = init_params(h,28)
    params_in(27) = init_params(h,34)
    params_in(28) = init_params(h,36)
    params_in(29) = init_params(h,37)
    params_in(30) = init_params(h,38)

!setup ADC
    do i = 1,11
      adc(i) = real(init_params(h,i+11),kind(sp))
    enddo

!    print *,'ADC: ',adc

!    do i = 1,30
!      print *,params_in(i),init_params(h,i)
!    enddo

  !need to set non-optimized parameters here again due to way sce passes things after first iteration...
!    a(21) = nmf(1)
!    a(22) = tipm(1)
!    a(23) = mbase(1)
!    a(24) = plwhc(1)
!    a(25) = daygm(1)
!    a(26) = adimp(1)
!    a(27) = pctim(1)
!    a(28) = riva
!    a(29) = side
!    a(30) = rserv

    call spin_up_first_year(a,params_in,spinup_crit,hru_tair(h,:),hru_precip(h,:),hru_pet(h,:), &
                            fuse_raim(h,:),uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
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
    call sfc_pressure(hru_elev(h),pa)

  ! print *,'sfc_pres',pa,elev    

    end_pt = sim_length

!    call calc_pet_pt(a)

  !print *,'pet',pet(100),end_pt

  !print *,uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp

    !now run model
    do i = 1,end_pt

    !set single precision inputs
      tair_sp   = real(hru_tair(h,i),kind(sp))
      precip_sp = real(hru_precip(h,i),kind(sp))
      pet_sp    = real(hru_pet(h,i),kind(sp))
      fuse_sp   = real(fuse_raim(h,i),kind(sp))

	!precip multiplier
  !	precip_sp = precip_sp*1.9
!AJN fuse raim switch
      CALL EXSNOW19(int(dt),int(dt/sec_hour),day(i),month(i),year(i),&
	!SNOW17 INPUT AND OUTPUT VARIABLES
			  precip_sp,tair_sp,raim_snow17(i),sneqv(i),snow(i),snowh(i),&
	!SNOW17 PARAMETERS
!ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!			    alat,a(1),a(2),a(3),a(4),a(5),a(7),a(8),a(9),&
!			    a(6),a(10),a(11),elev,pa,adc(1),&
			  real(lat,kind(sp)),params_in(1),params_in(2),params_in(3),&
                          params_in(4),params_in(5),params_in(21),params_in(22),params_in(23),&
			  params_in(6),params_in(24),params_in(25),real(elev,kind(sp)),&
                          real(pa,kind(sp)),adc,&
	!SNOW17 CARRYOVER VARIABLES
			  cs,tprev) 

! print *,'here 4',raim_snow17,i,pet_sp
! print *,uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp
!print *,a(7),a(8),a(12),a(27),a(26),a(28),a(15), &
!		      a(16),a(9),a(11),a(10),a(14),a(13),a(17),&
!		      a(29),a(30)
!print *,real(dt)
!AJN fuse raim switch
      call exsac(1,real(dt),raim_snow17(i),tair_sp,pet_sp,&
!      call exsac(1,real(dt),fuse_sp,tair_sp,pet_sp,&
	!SAC PARAMETERS
!UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
!REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
!SIDE,RSERV, &
!a(12),a(13),a(18),a(23),a(17),a(25),a(21),a(22),a(14),a(16),a(15),a(20),a(19),a(24),a(26),a(27)
!			a(12),a(13),a(18),a(23),a(17),a(25),a(21),&
!			a(22),a(14),a(16),a(15),a(20),a(19),a(24),&
!			a(26),a(27),&
		      params_in(7),params_in(8),params_in(12),params_in(27),&
                      params_in(26),params_in(28),params_in(15), &
		      params_in(16),params_in(9),params_in(11),params_in(10),&
                      params_in(14),params_in(13),params_in(17),&
		      params_in(29),params_in(30), &
	!SAC State variables
		      uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp,&
	!SAC OUTPUTS
		      qs(i),qg(i),tci(i),eta(i))

    !get total streamflow for basin
    !add hru output to running total
    total_tci(i) = total_tci(i) + tci(i)*hru_area(h)/tot_area
    total_qg(i)  = total_qg(i)  +  qg(i)*hru_area(h)/tot_area

!      if(i .ge. 20 .and. i .le. 30) then
!	print *,day(i),month(i),year(i),raim_snow17(i),tci(i),hru_area(h),tot_area,streamflow(i)
!      endif
    enddo  !end model loop


  enddo  !end hru loop

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
    call DUAMEL(total_tci,1,unit_hydro,a(18),a(19),dtuh,end_pt-1,m,route_tci,k,ntau)
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

!print *,route_tci_dp(1000),streamflow(1000)

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