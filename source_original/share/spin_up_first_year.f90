
subroutine spin_up_first_year(a, spinup_crit, uztwc, uzfwc, lztwc, &
                              lzfsc, lzfpc, adimc)
  use nrtype
  use constants, only: sec_hour
  use gauge_calib, only: sfc_pressure, calc_pet_pt
  use snow17_sac

  implicit none

!input variables
  real(dp), intent(in)			:: spinup_crit
  real(sp), dimension(30), intent(inout) 	:: a

!output variables
  real(dp), intent(out)		:: uztwc,uzfwc,lztwc
  real(dp), intent(out)		:: lzfsc,lzfpc,adimc

!local variables

  integer(I4B) :: i,cnt

  logical :: spin_up_flag

!previous sac-sma state
  real(dp)              :: uztwc_prev,uzfwc_prev,lztwc_prev
  real(dp)              :: lzfsc_prev,lzfpc_prev,adimc_prev


!diff variables
  real(dp)              :: uztwc_diff,uzfwc_diff,lztwc_diff
  real(dp)              :: lzfsc_diff,lzfpc_diff,adimc_diff

!single precision state variables
  real(sp)		:: uztwc_sp,uzfwc_sp,lztwc_sp
  real(sp)		:: lzfsc_sp,lzfpc_sp,adimc_sp

!single precision forcing variables
  real(sp)		:: tair_sp
  real(sp)		:: precip_sp
  real(sp)		:: pet_sp

!snow-17 carry over variables
  real(sp) :: tprev				!carry over variable
  real(sp), dimension(19)  :: cs		!carry over variable

!snow-17 output variables
  real(sp)				:: snowh, sneqv, snow, raim_snow17

!sac-sma output variables
  real(sp)				:: qs,qg,eta,tci


  real(dp)		:: pa	!surface pressure for snow-17

!code


!set spin up flag
  spin_up_flag = .true.

  cnt = 0

!setup initial state sac-sma
  uztwc_sp = init_smois(1)
  uzfwc_sp = init_smois(2)
  lztwc_sp = init_smois(3)
  lzfsc_sp = init_smois(4)
  lzfpc_sp = init_smois(5)
  adimc_sp = init_smois(6)

!reset snow-17 carryover variables
  tprev = 0.0
  cs = 0.0
!set sac output variables to zero
  qs = 0.0
  qg = 0.0
  tci = 0.0
  eta = 0.0

!get snow17 surface pressure
  call sfc_pressure(elev,pa)

!need to set non-optimized parameters here again due to way sce passes things back out
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

!  print *,a

!first calc pet
  call calc_pet_pt(a)


!  print *,'starting while loop'
!print *,pet(200),a(20)
  do while (spin_up_flag)
    !run model combo over first year (use first full water year
    !first set previous state variables to initial state variables
    uztwc_prev = real(uztwc_sp,kind(dp))
    uzfwc_prev = real(uzfwc_sp,kind(dp))
    lztwc_prev = real(lztwc_sp,kind(dp))
    lzfsc_prev = real(lzfsc_sp,kind(dp))
    lzfpc_prev = real(lzfpc_sp,kind(dp))
    adimc_prev = real(adimc_sp,kind(dp))

!    print *,'running wy'
!  print *,'before',uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp
    !run first water year to spin things up
    do i = 1,365
    !set single precision inputs
      tair_sp   = real(tair(i),kind(sp))
      precip_sp = real(precip(i),kind(sp))
      pet_sp    = real(pet(i),kind(sp))
            !precip multiplier
!	precip_sp = precip_sp*1.9
!print *,i,pet(i),precip(i),tair(i)

      CALL EXSNOW19(int(dt),int(dt/sec_hour),day(i),month(i),year(i),&
	!SNOW17 INPUT AND OUTPUT VARIABLES
			  precip_sp,tair_sp,raim_snow17,sneqv,snow,snowh,&
	!SNOW17 PARAMETERS
!ALAT,SCF,MFMAX,MFMIN,UADJ,SI,NMF,TIPM,MBASE,PXTEMP,PLWHC,DAYGM,ELEV,PA,ADC
!			    alat,a(1),a(2),a(3),a(4),a(5),a(7),a(8),a(9),&
!			    a(6),a(10),a(11),elev,pa,adc(1),&
			  real(lat,kind(sp)),a(1),a(2),a(3),a(4),a(5),a(21),a(22),a(23),&
			  a(6),a(24),a(25),real(elev,kind(sp)),real(pa,kind(sp)),adc,&
	!SNOW17 CARRYOVER VARIABLES
			  cs,tprev) 

! print *,'here 4',raim_snow17,i,pet_sp
! print *,uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp
!print *,a(7),a(8),a(12),a(27),a(26),a(28),a(15), &
!		      a(16),a(9),a(11),a(10),a(14),a(13),a(17),&
!		      a(29),a(30)
!print *,real(dt)

      call exsac(1,real(dt),raim_snow17,tair_sp,pet_sp,&
	!SAC PARAMETERS
!UZTWM,UZFWM,UZK,PCTIM,ADIMP,RIVA,ZPERC, &
!REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE, &
!SIDE,RSERV, &
!a(12),a(13),a(18),a(23),a(17),a(25),a(21),a(22),a(14),a(16),a(15),a(20),a(19),a(24),a(26),a(27)
!			a(12),a(13),a(18),a(23),a(17),a(25),a(21),&
!			a(22),a(14),a(16),a(15),a(20),a(19),a(24),&
!			a(26),a(27),&
		      a(7),a(8),a(12),a(27),a(26),a(28),a(15), &
		      a(16),a(9),a(11),a(10),a(14),a(13),a(17),&
		      a(29),a(30), &
	!SAC State variables
			uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp,&
	!SAC OUTPUTS
			qs,qg,tci,eta)

    enddo !end model loop

    uztwc = real(uztwc_sp,kind(dp))
    uzfwc = real(uzfwc_sp,kind(dp))
    lztwc = real(lztwc_sp,kind(dp))
    lzfsc = real(lzfsc_sp,kind(dp))
    lzfpc = real(lzfpc_sp,kind(dp))
    adimc = real(adimc_sp,kind(dp))
!print *,'after',uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc
    !check for convergence
    !the units of all six state variables are mm

    !state variables are:
    !uztwc:	upper-zone tension water storage content
    !uzfwc:	upper-zone free water storage content
    !lztwc:	lower-zone tension water storage content
    !lzfpc:	lower-zone free primary water storage content
    !lzfsc:	lower-zone free secondary water storage content
    !adimc:	additional impervious area content

    uztwc_diff = abs(uztwc-uztwc_prev)
    uzfwc_diff = abs(uzfwc-uzfwc_prev)
    lztwc_diff = abs(lztwc-lztwc_prev)
    lzfsc_diff = abs(lzfsc-lzfsc_prev)
    lzfpc_diff = abs(lzfpc-lzfpc_prev)
    adimc_diff = abs(adimc-adimc_prev)

    cnt = cnt + 1
    !print *,cnt,uztwc_diff,uzfwc_diff,lztwc_diff,lzfsc_diff,lzfpc_diff,adimc_diff
    !print *,cnt,lztwc,lztwc_prev,lztwc_diff
    !print *,cnt,lzfpc,lzfpc_prev,lzfpc_diff
!      print *,cnt,uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc

    if(uztwc_diff .le. spinup_crit .and. uzfwc_diff .le. spinup_crit .and. lztwc_diff .le. spinup_crit .and. &
	lzfsc_diff .le. spinup_crit .and. lzfpc_diff .le. spinup_crit .and. adimc_diff .le. spinup_crit) then
      spin_up_flag = .false.

    endif

    if(cnt .gt. 50) then
      spin_up_flag = .false.
!	print *,'failed'
    endif
  enddo  !end while loop for spin up
!print *,cnt

  return
end subroutine spin_up_first_year
