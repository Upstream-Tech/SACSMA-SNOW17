subroutine write_snow17_state(cs,tprev)
  use nrtype
  use snow17_sac, only: snow17_state_out_file

  implicit none

  !input variables
  real(sp), intent(in) 			:: tprev				!carry over variable
  real(sp), dimension(:), intent(in)	:: cs					!carry over array

  !local variables
  integer(I4B)	:: i

  open(unit=95,FILE=trim(snow17_state_out_file),FORM='formatted')

  do i = 1,19
    write(95,*) cs(i)
  enddo
  write(95,*) tprev
  close(unit=95)

  return
end subroutine write_snow17_state

!ccccccccccccccccccccccccccccccc

subroutine write_sac_state(uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
  use nrtype
  use snow17_sac, only: sac_state_out_file
 
  implicit none

  !input variables
  real(sp), intent(in)	:: uztwc					!state variable
  real(sp), intent(in)	:: uzfwc					!state variable
  real(sp), intent(in)	:: lztwc					!state variable
  real(sp), intent(in)	:: lzfsc					!state variable
  real(sp), intent(in)	:: lzfpc					!state variable
  real(sp), intent(in)	:: adimc					!state variable

  open(unit=95,FILE=trim(sac_state_out_file),FORM='formatted')

  write(95,*) uztwc
  write(95,*) uzfwc
  write(95,*) lztwc
  write(95,*) lzfsc
  write(95,*) lzfpc
  write(95,*) adimc

  close(unit=95)

  return

end subroutine write_sac_state

!cccccccccccccccccccccccccccccccccccccccccc

subroutine read_uh_state(uh_flow,uh_length,sim_length)
  use nrtype
  use snow17_sac, only: uh_state_in_file

  implicit none

  !input variable
  integer(I4B), intent(in)				:: uh_length,sim_length

  !output variables
  real(sp), dimension(:), intent(out) 			:: uh_flow

  !local variables
  integer(I4B)		:: i,end_pt

!  if(sim_length .lt. uh_length) then
!    end_pt = sim_length
!  else
!    end_pt = uh_length
!  endif
  end_pt = uh_length-1
  open(unit=95,FILE=trim(uh_state_in_file),FORM='formatted',status='old')

  do i = 1,end_pt
    read(95,*) uh_flow(i)
  enddo

  close(unit=95)

  return
end subroutine read_uh_state

!cccccccccccccccccccccccccccccccccccccccccc

subroutine write_uh_state(route_tci,old_uh,sim_length,uh_length)
  use nrtype
  use snow17_sac, only: uh_state_out_file

  implicit none

  !input variables
  integer(I4B), intent(in)				:: uh_length
  integer(I4B), intent(in)				:: sim_length
  real(sp), dimension(:), intent(in) 			:: route_tci
  real(sp), dimension(:), intent(in) 			:: old_uh


  !local variables
  integer(I4B)				:: i,len_old,start_old
  real(sp),allocatable,dimension(:)	:: out_uh

  allocate(out_uh(uh_length))


  start_old = sim_length

!print *,'lere',start_old,uh_length,sim_length
  out_uh = 0.0

  if(uh_length-start_old .lt. 1) then
    do i = 1,uh_length
      out_uh(i) = route_tci(i)
    enddo
  else
    do i = 1,uh_length
      if(i .lt. uh_length-start_old) then
	out_uh(i) = route_tci(i) + old_uh(i+start_old-1)
!print *,old_uh(i+start_old),route_tci(i)
      else
	out_uh(i) = route_tci(i)
      endif
    enddo
  endif

  open(unit=95,FILE=trim(uh_state_out_file),FORM='formatted')

  !need to take the last time input channel flow from sac and route it
  !through uh and output
  !also need to consider old_uh and add any flow that occurs in the future from the end of the simulation
  do i = 2,uh_length
    write(95,*) out_uh(i)
  enddo

  close(unit=95)

  return
end subroutine write_uh_state

!ccccccccccccccccccccccccccccccc

subroutine read_snow17_state(cs,tprev)
  use nrtype
  use snow17_sac, only: snow17_state_in_file

  implicit none

  !input variables
  real(sp), intent(out) 			:: tprev				!carry over variable
  real(sp), dimension(:), intent(out)	:: cs					!carry over array

  !local variables
  integer(I4B)	:: i

  open(unit=95,FILE=trim(snow17_state_in_file),FORM='formatted',status='old')

  do i = 1,19
    read(95,*) cs(i)
  enddo
  read(95,*) tprev
  close(unit=95)

  return
end subroutine read_snow17_state

!ccccccccccccccccccccccccccccccc

subroutine read_sac_state(uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
  use nrtype
  use snow17_sac, only: sac_state_in_file
 
  implicit none

  !input variables
  real(sp), intent(out)	:: uztwc					!state variable
  real(sp), intent(out)	:: uzfwc					!state array
  real(sp), intent(out)	:: lztwc					!state array
  real(sp), intent(out)	:: lzfsc					!state array
  real(sp), intent(out)	:: lzfpc					!state array
  real(sp), intent(out)	:: adimc					!state array

  open(unit=95,FILE=trim(sac_state_in_file),FORM='formatted',status='old')

  read(95,*) uztwc
  read(95,*) uzfwc
  read(95,*) lztwc
  read(95,*) lzfsc
  read(95,*) lzfpc
  read(95,*) adimc

  close(unit=95)

  return
end subroutine read_sac_state

!ccccccccccccccccccccccccccccccc

subroutine get_sim_length(sim_length)
  use nrtype
  use snow17_sac, only: forcing_name, start_year,start_day,start_month, &
                        end_year,end_month,end_day

  !output variables
  integer(I4B),intent(out)	:: sim_length

  !local variables
  integer(I4B)			:: ios=0
  integer(I4B)			:: dum_int
  integer(I4B)			:: year,month,day
  integer(I4B)			:: read_flag

  character(len=64), parameter		:: read_format = "(I4.4, 3(1x,I2.2),1x,7(F10.2))"
  character(len = 1024)			:: dum_str
  real(dp)		:: dum_real

  

!code
  read_flag = 0
  sim_length = 0

  open (UNIT=50,file=trim(forcing_name),form='formatted',status='old')

  read (UNIT=50,FMT='(F7.2)') dum_real
  read (UNIT=50,FMT='(F7.2)') dum_real
  read (UNIT=50,FMT='(F11.0)') dum_real
  read (UNIT=50,FMT='(80A)') dum_str


  do while(ios .ge. 0)
    read (UNIT=50,FMT=read_format,IOSTAT=ios) year,month,day,dum_int,&
				dum_real,dum_real,dum_real,dum_real,dum_real,&
				dum_real,dum_real

    if(year .eq. start_year .and. month .eq. start_month .and. day .eq. start_day) then
      read_flag = 1
    end if

    if(read_flag .eq. 1) then
      sim_length = sim_length + 1
    end if

    if(year .eq. end_year .and. month .eq. end_month .and. day .eq. end_day) then
      read_flag = 0
    end if
  end do

  close(unit=50)

end subroutine get_sim_length

!ccccccccccccccccccccccccccccccc

subroutine read_areal_forcing(year,month,day,hour,tmin,tmax,vpd,dayl,swdown,precip)
  use nrtype
  use snow17_sac, only: forcing_name, start_year,start_day,start_month, &
                        end_year,end_month,end_day,lat,area_basin,elev

  implicit none

!output variables
  integer(I4B),dimension(:),intent(out)	:: year
  integer(I4B),dimension(:),intent(out)	:: month
  integer(I4B),dimension(:),intent(out)	:: day
  integer(I4B),dimension(:),intent(out)	:: hour
  real(dp),dimension(:),intent(out)	:: tmin
  real(dp),dimension(:),intent(out)	:: tmax
  real(dp),dimension(:),intent(out)	:: vpd
  real(dp),dimension(:),intent(out)	:: dayl
  real(dp),dimension(:),intent(out)	:: swdown
  real(dp),dimension(:),intent(out)	:: precip

!local variables
  integer(I4B)				:: i,ios=0
  integer(I4B)				:: yr,mnth,dy,hr
  integer(I4B)				:: read_flag

  character(len=64), parameter		:: read_format = "(I4.4, 3(1x,I2.2),1x,7(F10.2))"
  character(len = 1024)			:: dum_str
  real(DP)				:: swe
  real(DP)				:: dum_real

  real(DP)				:: dl,pcp,sw,tma,tmn,vp




!code
  i = 1
!read met file
  open (UNIT=50,file=trim(forcing_name),form='formatted',status='old')

  read (UNIT=50,FMT='(F7.2)') lat
  read (UNIT=50,FMT='(F7.2)') elev
  read (UNIT=50,FMT='(F11.0)') area_basin
  read (UNIT=50,FMT='(80A)') dum_str
  !read the rest of the input file
  !forcing_offset is the first day of the first complete water year in corresponding observed streamflow
  !this is the point at which we want to keep the forcing data
  !keep through the end of the sim_len or end of file (if validation period)
  !need to do this because forcing starts 01 Jan 1980
  !observed streamflow varies in its start date by gauge

  do while(ios .ge. 0)
    read (UNIT=50,FMT=*,IOSTAT=ios) yr,mnth,dy,hr,&
			    dl,pcp,sw,swe,tma,&
			    tmn,vp

    if(yr .eq. start_year .and. mnth .eq. start_month .and. dy .eq. start_day) then
      read_flag = 1
    end if

    if(read_flag .eq. 1) then
      year(i)	= yr
      month(i)	= mnth
      day(i)	= dy
      hour(i)	= hr
      dayl(i)	= dl
      precip(i)	= pcp
      swdown(i) = sw
      tmax(i)	= tma
      tmin(i)	= tmn
      vpd(i)	= vp
      i = i + 1

    end if

    if(yr .eq. end_year .and. mnth .eq. end_month .and. dy .eq. end_day) then
      read_flag = 0
    end if

  end do

  close(unit=50)

  return

end subroutine read_areal_forcing

!cccccccccccccccccccccccccccccccccccccccccc
!change & simplify
subroutine read_streamflow(streamflow)
  use nrtype
  use constants, only: sec_day,cfs_cms
  use snow17_sac, only: stream_name,area_basin,start_year,start_day,start_month, &
			end_year,end_month,end_day

  implicit none

!output variables
  real(dp),dimension(:),intent(out)	:: streamflow

!local variables
  integer(I4B) :: i
  integer(I4B) :: yr,mn,dy,gauge,ios=0,error
  real(dp)    :: flow

  logical	:: valid

!this subroutine assumes streamflow is daily data of the following format:
  character(len=64), parameter :: read_format = "(I8.8,1x,I4.4, 2(1x,I2.2),1x,1(F8.2))"


!code
  valid = .false.
  i = 1

! open streamflow file
  open (UNIT=50,file=stream_name,form='formatted',status='old')


  do while(ios .ge. 0)
    read (UNIT=50,FMT=read_format,IOSTAT=ios) gauge,yr,mn,dy,flow

    if(yr .eq. start_year .and. mn .eq. start_month .and. dy .eq. start_day) then
      valid = .true.
    end if

    if(valid) then
      streamflow(i)	= flow

      !convert flow to mm/day
      !convert streamflow (cfs) to cms
      streamflow(i) = streamflow(i)*cfs_cms  !now in cubic meters per second

      !need to convert to mm/day
      streamflow(i) = streamflow(i)*sec_day !now in cubic meters per day m^3 day^-1
                                           !m^3/day

      !1 cubic meter per day is 1000 mm per square m -> mm*m^2/day
      streamflow(i) = streamflow(i)*1000./area_basin  !now in mm/day
      
      i = i + 1
    end if

    if(yr .eq. end_year .and. mn .eq. end_month .and. dy .eq. end_day) then
      valid = .false.
    end if

  end do

  close(unit=50)

  return

end subroutine read_streamflow

!ccccccccccccccccccccccccccccccc
!
! subroutines for reading param files: sac, snow17, unit hydrograph & pet
!
!ccccccccccccccccccccccccccccccc

subroutine read_sac_params(param_name)
  use nrtype
  use snow17_sac, only: uztwm,uzfwm,uzk,pctim,adimp,zperc,rexp, &
			lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,riva,side,rserv

  implicit none
 
!input variables
  character(len=1024),intent(in)	:: param_name
 
!local variables
  character(len=50)		:: param
  
  real(sp)			:: value

  integer(I4B)			:: ios=0
 
  open(unit=50,file=trim(param_name))

  do while(ios .eq. 0)
    read(unit=50,FMT=*,IOSTAT=ios) param,value

    if(param == 'uztwm') then
      uztwm = value
    else if(param == 'uzfwm') then
      uzfwm = value
    else if(param == 'uzk') then
      uzk = value
    else if(param == 'pctim') then
      pctim = value
    else if(param == 'adimp') then
      adimp = value
    else if(param == 'zperc') then
      zperc = value
    else if(param == 'rexp') then
      rexp = value
    else if(param == 'lztwm') then
      lztwm = value
    else if(param == 'lzfsm') then
      lzfsm = value
    else if(param == 'lzfpm') then
      lzfpm = value
    else if(param == 'lzsk') then
      lzsk = value
    else if(param == 'lzpk') then
      lzpk = value
    else if(param == 'pfree') then
      pfree = value
    else if(param == 'riva') then
      riva = value
    else if(param == 'side') then
      side = value
    else if(param == 'rserv') then
      rserv = value
    end if
  end do
  close(unit=50)

  return
end subroutine read_sac_params

!ccccccccc
subroutine read_snow17_params(param_name)
  use nrtype
  use snow17_sac, only: scf,mfmax,mfmin,uadj,si,pxtemp,nmf,&
                        tipm,mbase,plwhc,daygm,adc

  implicit none
 
!input variables
  character(len=1024),intent(in)	:: param_name

!local variables
  character(len=50)		:: param
  
  real(sp)			:: value

  integer(I4B)			:: ios=0

  open(unit=50,file=trim(param_name))

  do while(ios .eq. 0)
    read(unit=50,FMT=*,IOSTAT=ios) param,value

    if(param == 'mfmax') then
      mfmax = value
    else if(param == 'mfmin') then
      mfmin = value
    else if(param == 'scf') then
      scf = value
    else if(param == 'uadj') then
      uadj = value
    else if(param == 'si') then
      si = value
    else if(param == 'pxtemp') then
      pxtemp = value
    else if(param == 'nmf') then
      nmf = value
    else if(param == 'tipm') then
      tipm = value
    else if(param == 'mbase') then
      mbase = value
    else if(param == 'plwhc') then
      plwhc = value
    else if(param == 'daygm') then
      daygm = value
    else if(param == 'adc1') then
      adc(1) = value
    else if(param == 'adc2') then
      adc(2) = value
    else if(param == 'adc3') then
      adc(3) = value
    else if(param == 'adc4') then
      adc(4) = value
    else if(param == 'adc5') then
      adc(5) = value
    else if(param == 'adc6') then
      adc(6) = value
    else if(param == 'adc7') then
      adc(7) = value
    else if(param == 'adc8') then
      adc(8) = value
    else if(param == 'adc9') then
      adc(9) = value
    else if(param == 'adc10') then
      adc(10) = value
    else if(param == 'adc11') then
      adc(11) = value
    end if
  end do

  close(unit=50)

  return
end subroutine read_snow17_params

!cccccccccccc
subroutine read_uhp_params(param_name)
  use nrtype
  use snow17_sac, only: unit_shape,unit_scale,pet_coef

  implicit none
 
!input variables
  character(len=1024),intent(in)	:: param_name

!local variables
  character(len=50)		:: param
  
  real(sp)			:: value

  integer(I4B)			:: ios=0

  open(unit=50,file=trim(param_name))
  
  do while(ios .eq. 0)
    read(unit=50,FMT=*,IOSTAT=ios) param,value

    if(param == 'unit_shape') then
      unit_shape = value
    else if(param == 'unit_scale') then
      unit_scale = value
    else if(param == 'pet_coef') then
      pet_coef = value
    end if
  end do
  close(unit=50)

  return
end subroutine read_uhp_params
