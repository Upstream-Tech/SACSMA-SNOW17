program lump_driver
  use nrtype
  use snow17_sac
  use constants, only: sec_hour,sec_day
  use gauge_calib, only: read_namelist, calc_pet_pt, sfc_pressure, &
                         read_areal_forcing,read_streamflow, get_sim_length, &
			 julian_day,write_snow17_state,write_sac_state, &
			 read_sac_state,read_snow17_state,read_uh_state, &
			 write_uh_state

  implicit none


!local variables

  character(len=2000)		:: arg  !command line arg for namelist file
  character(len=2000)		:: namelist_name  !command line arg for namelist file
  integer(I4B) :: i,ntau,k,m,j,cnt
  integer(I4B) :: end_pt		!end point for old uh state addition
  integer(I4B) :: error
  integer(I4B) :: sim_length	!length of simulation
  real(sp) :: dtuh	!for unit hydrograph

!single precision sac-sma state variables
  real(sp)		:: uztwc_sp
  real(sp)		:: uzfwc_sp
  real(sp)		:: lztwc_sp
  real(sp)		:: lzfsc_sp
  real(sp)		:: lzfpc_sp
  real(sp)		:: adimc_sp

  real(sp)		:: pet_sp
  real(sp)		:: tair_sp
  real(sp)		:: precip_sp


!snow-17 carry over variables
  real(sp) :: tprev				!carry over variable
  real(sp), dimension(19)  :: cs			!carry over variable

!snow-17 surface pressure
  real(dp) :: pa

!unit hydrograph
  real(sp),dimension(1000)       :: unit_hydro
  integer(I4B)			 :: uh_length

!restart array for unit hydrograph
  real(sp),dimension(:),allocatable		:: uh_state

!ALLOCATABLE VARIABLES

!sac-sma state variables
  real(dp),dimension(:),allocatable              :: uztwc_dp
  real(dp),dimension(:),allocatable              :: uzfwc_dp
  real(dp),dimension(:),allocatable              :: lztwc_dp
  real(dp),dimension(:),allocatable              :: lzfsc_dp
  real(dp),dimension(:),allocatable              :: lzfpc_dp
  real(dp),dimension(:),allocatable              :: adimc_dp


!sac-sma output variables and routed flow
  real(sp), dimension(:),allocatable	:: qs,qg,eta,tci,route_tci

!snow-17 output variables  single precision
  real(sp), dimension(:),allocatable    :: snowh, sneqv, snow !output variables for snow-17

!date variables
  integer, dimension(:),allocatable :: year,month,day,hour,jday
!  integer(I4B), dimension(36500) :: year,month,day,hour

!observed streamflow
  real(dp),dimension(:),allocatable :: streamflow
!  real(dp), dimension(36500)	:: streamflow



  real(sp),dimension(:),allocatable :: raim
!  real(dp), dimension(36500) :: raim

!Atmospheric forcing variables
  real(dp), dimension(:),allocatable    :: tmin,tmax,vpd,dayl,swdown,precip
!derived forcing variables
  real(dp), dimension(:),allocatable :: pet,tair


!
!   code starts below
!

  uh_length = size(unit_hydro,1)

!get namelist filename
  i = 0
  do
    call get_command_argument(i,arg)
    if(i .eq. 1) namelist_name=arg
    if(LEN_TRIM(arg) == 0) EXIT
    i = i + 1
  end do

!read namelist file
  call read_namelist(namelist_name)

!read parameter files
  call read_sac_params(sac_param_file)
  call read_snow17_params(snow17_param_file)
  call read_uhp_params(uhp_param_file)

  call get_sim_length(sim_length)

!allocate variables now
!forcing variables
  allocate(year(sim_length))
  allocate(month(sim_length))
  allocate(day(sim_length))
  allocate(jday(sim_length))
  allocate(hour(sim_length))
  allocate(tmax(sim_length))
  allocate(tmin(sim_length))
  allocate(vpd(sim_length))
  allocate(dayl(sim_length))
  allocate(swdown(sim_length))
  allocate(precip(sim_length))
  allocate(pet(sim_length))
  allocate(tair(sim_length))

!observations
  allocate(streamflow(sim_length))

!sac-sma state variables
  allocate(uztwc_dp(sim_length))
  allocate(uzfwc_dp(sim_length))
  allocate(lztwc_dp(sim_length))
  allocate(lzfsc_dp(sim_length))
  allocate(lzfpc_dp(sim_length))
  allocate(adimc_dp(sim_length))

!sac-sma output variables
  allocate(qg(sim_length))
  allocate(qs(sim_length))
  allocate(eta(sim_length))
  allocate(tci(sim_length+uh_length))
  allocate(route_tci(sim_length+uh_length))

!snow-17 output variables
  allocate(snowh(sim_length))
  allocate(sneqv(sim_length))
  allocate(snow(sim_length))
  allocate(raim(sim_length))


!print *,'here',sim_length

!read forcing data
  call read_areal_forcing(year,month,day,hour,tmin,tmax,&
			  vpd,dayl,swdown,precip)
!print *,'here'
!read streamflow
  if(stream_name == '-9999' .or. stream_name == '-9999.0' .or. stream_name == '-9' &
     .or. stream_name == '-99' .or. stream_name == '-999' .or. stream_name == '-9.0' &
     .or. stream_name == '-99.0' .or. stream_name == '-999.0') then
    streamflow = -9999.0
  else
    call read_streamflow(streamflow)
  endif

!print *,'here'
!get day of year array setup
  call julian_day(year,month,day,jday)
!print *,'here'
!compute pet
  tair = (tmax+tmin)/2.0_dp
!print *,'here'
  call calc_pet_pt(jday,tmax,tmin,tair,vpd,swdown,dayl,pet)
!print *,'here'
!run model

!get sfc_pressure
  call sfc_pressure(elev,pa)

  allocate(uh_state(uh_length))
  uh_state = 0.0
!set single precision sac state variables to initial values
  if(restart_run .eq. 0) then
    uztwc_sp = in_uztwc
    uzfwc_sp = in_uzfwc
    lztwc_sp = in_lztwc
    lzfsc_sp = in_lzfsc
    lzfpc_sp = in_lzfpc
    adimc_sp = in_adimc
    cs(1) = in_swe
  elseif(restart_run .gt. 0) then
    !if restart run: read in output state files
    !this overwrites namelist sac state variables and initial swe

    !get snow17 states
    call read_snow17_state(cs,tprev)

    !get sac-sma states
    call read_sac_state(uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp)

    !also need to get in water in channel
    call read_uh_state(uh_state,uh_length,sim_length)

  endif


  do i = 1,sim_length,1

  !set single precision inputs
    tair_sp   = real(tair(i),kind(sp))
    precip_sp = real(precip(i),kind(sp))
    pet_sp    = real(pet(i),kind(sp))
!print *,'Pressure',pa
!print *,day(i),month(i),year(i),precip_sp,tair_sp,lat,elev,raim(i)
!print *,tprev,cs
!print *,'melt!',mfmax,mfmin
!tair_sp = 10.0
    CALL EXSNOW19(int(dt),int(dt/sec_hour),day(i),month(i),year(i),&
	!SNOW17 INPUT AND OUTPUT VARIABLES
			  precip_sp,tair_sp,raim(i),sneqv(i),snow(i),snowh(i),&
	!SNOW17 PARAMETERS
!ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!			    alat,a(1),a(2),a(3),a(4),a(5),a(7),a(8),a(9),&
!			    a(6),a(10),a(11),elev,pa,adc(1),&
			  real(lat,kind(sp)),scf,mfmax,mfmin,uadj,si,nmf,tipm,mbase,&
			  pxtemp,plwhc,daygm,real(elev,kind(sp)),real(pa,kind(sp)),adc,&
	!SNOW17 CARRYOVER VARIABLES
			  cs,tprev) 

! print *,'here 4'

!print *,'snow',snow(i),snowh(i)

    call exsac(1,real(dt),raim(i),tair_sp,pet_sp,&
	!SAC PARAMETERS
!UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
!REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
!SIDE,RSERV, &
!a(12),a(13),a(18),a(23),a(17),a(25),a(21),a(22),a(14),a(16),a(15),a(20),a(19),a(24),a(26),a(27)
!			a(12),a(13),a(18),a(23),a(17),a(25),a(21),&
!			a(22),a(14),a(16),a(15),a(20),a(19),a(24),&
!			a(26),a(27),&
		      uztwm,uzfwm,uzk,pctim,adimp,riva,zperc, &
		      rexp,lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,&
		      side,rserv, &
	!SAC State variables
			uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp,&
	!SAC OUTPUTS
			qs(i),qg(i),tci(i),eta(i))


    !place state variables in output arrays
    uztwc_dp(i) = real(uztwc_sp,kind(dp))
    uzfwc_dp(i) = real(uzfwc_sp,kind(dp))
    lztwc_dp(i) = real(lztwc_sp,kind(dp))
    lzfsc_dp(i) = real(lzfsc_sp,kind(dp))
    lzfpc_dp(i) = real(lzfpc_sp,kind(dp))
    adimc_dp(i) = real(adimc_sp,kind(dp))
  enddo  !end simulation loop


  dtuh = real(dt/sec_day)


!route flow using UH

  if (unit_shape .le. 0.0 .and. unit_scale .le. 0.0) THEN
    k = 0
    m = 1
  else
    k = 1
    m = 1000
  end if
  ntau = 0

!call unit hydrograph routine
  if(unit_shape .gt. 0.0) then
    call DUAMEL(tci,1,unit_hydro,unit_shape,unit_scale,dtuh,sim_length+uh_length,m,route_tci,k,ntau)
  endif

!do i=-10,20
!do i = 1,sim_length+40
!!  print *, unit_hydro(i), route_tci(sim_length+i)
!  print *, route_tci(i)
!enddo
!print *,uh_state(1:10)

  if(restart_run .gt. 0) then
  !!!!!
  !add prior UH state to routed flow array...
      
    if(sim_length .lt. uh_length) then
      end_pt = sim_length
    else
      end_pt = uh_length
    endif

    do i = 1,end_pt
!        route_tci(i) = route_tci(i) + sum(uh_state(i:uh_length-1))
      route_tci(i) = route_tci(i) + uh_state(i)
    enddo



  endif


  31 FORMAT(I4.4, 3(1x,I2.2),16(F15.7))

  !output file open
  open(unit=45,FILE=trim(model_out),FORM='formatted')


  !write out
  !print header first
  write(45,*) 'year month day hour swe pcp pcp*scf raim et pet tair uztwc uzfwc lztwc lzfsc lzfpc adimc ,&
               mod_flow_unrouted mod_flow_routed obs_flow'
  do i = 1,sim_length
    if(unit_shape .gt. 0) then
      write(45,31) year(i),month(i),day(i),hour(i)+12,sneqv(i)*1000.,precip(i),precip(i)*scf,raim(i),&
                        eta(i),pet(i),tair(i),uztwc_dp(i),uzfwc_dp(i),lztwc_dp(i),lzfsc_dp(i),lzfpc_dp(i),&
                        adimc_dp(i),tci(i),route_tci(i),streamflow(i)
    else
      write(45,31) year(i),month(i),day(i),hour(i)+12,sneqv(i)*1000.,precip(i),precip(i)*scf,raim(i),&
                        eta(i),pet(i),tair(i),uztwc_dp(i),uzfwc_dp(i),lztwc_dp(i),lzfsc_dp(i),lzfpc_dp(i),&
                        adimc_dp(i),tci(i),route_tci(i),streamflow(i)
    endif
  enddo
  close(unit=45)

  if(write_restart .gt. 0) then

    !need to write out a snow17 state file for carryover variables and tprev
    call write_snow17_state(cs,tprev)
    !need to write a SAC-SMA state file
    call write_sac_state(uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp)
    !need to write Unit hydrograph state file
!    if(sim_length .lt. uh_length) then
!      call write_uh_state(tci(sim_length),unit_hydro,uh_state(sim_length+1:uh_length-1),uh_length)
    call write_uh_state(route_tci(sim_length:sim_length+uh_length),uh_state,sim_length,uh_length)
!    else
!      call write_uh_state(tci(sim_length),unit_hydro,uh_state(uh_length-2:uh_length-1),uh_length)
!      call write_uh_state(tci(sim_length),unit_hydro,uh_state(1:uh_length-1),uh_length)
!    endif
  endif


end program
