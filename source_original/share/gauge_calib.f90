module gauge_calib

  interface 

    subroutine julian_day(year,month,day,jday)
      use nrtype
      !input variables
      integer(I4B),dimension(:),intent(in) :: year
      integer(I4B),dimension(:),intent(in) :: month
      integer(I4B),dimension(:),intent(in) :: day

      !output variables
      integer(I4B),dimension(:),intent(out) :: jday

    end subroutine julian_day


    subroutine julianday_scalar(iyear,imonth,iday,jday_scalar)
      use nrtype

      integer(I4B), intent(in) :: iyear
      integer(I4B), intent(in) :: imonth
      integer(I4B), intent(in) :: iday

      integer(I4B),intent(out) :: jday_scalar
    end subroutine julianday_scalar


    subroutine read_namelist(namelist_name)
      use nrtype
      !input variable
      character(len=2000),intent(in)	:: namelist_name

    end subroutine read_namelist

    subroutine sfc_pressure(elev, pres)
      use nrtype
      use constants, only: sfc_pres_a,sfc_pres_b,sfc_pres_c,&
                       sfc_pres_d,sfc_pres_e

      real(dp), intent(in)		:: elev

      real(dp), intent(out)		:: pres
    end subroutine sfc_pressure


    subroutine calc_pet_pt(jday,tmax,tmin,tair,vpd,swdown,dayl,pet)
      use constants
      use nrtype
      use snow17_sac, only: elev,lat,pet_coef

      !input variable
      integer(I4B), dimension(:), intent(in) 	:: jday
      real(dp), dimension(:), intent(in) 	:: tair
      real(dp), dimension(:), intent(in) 	:: tmax
      real(dp), dimension(:), intent(in) 	:: tmin
      real(dp), dimension(:), intent(in) 	:: vpd
      real(dp), dimension(:), intent(in) 	:: swdown
      real(dp), dimension(:), intent(in) 	:: dayl

      !output variable
      real(dp),dimension(:),intent(out)	:: pet
    end subroutine calc_pet_pt

    subroutine get_sim_length(sim_length)
      use nrtype
      use snow17_sac, only: forcing_name, start_year,start_day,start_month, &
                        end_year,end_month,end_day

      integer(I4B), intent(out)	:: sim_length	
    end subroutine get_sim_length


    subroutine read_streamflow(streamflow)
      use nrtype
      use constants,  only: cfs_cms, sec_day
      use snow17_sac, only: stream_name,area_basin

      !output variables
      real(dp),dimension(:),intent(out)	:: streamflow
    end subroutine read_streamflow

    subroutine write_snow17_state(cs,tprev)
      use nrtype
      use snow17_sac, only: snow17_state_out_file

      !input variables
      real(sp), intent(in) 			:: tprev				!carry over variable
      real(sp), dimension(:), intent(in)	:: cs					!carry over array

    end subroutine write_snow17_state

    subroutine write_sac_state(uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
      use nrtype
      use snow17_sac, only: sac_state_out_file

      !input variables
      real(sp), intent(in)	:: uztwc					!state variable
      real(sp), intent(in)	:: uzfwc					!state array
      real(sp), intent(in)	:: lztwc					!state array
      real(sp), intent(in)	:: lzfsc					!state array
      real(sp), intent(in)	:: lzfpc					!state array
      real(sp), intent(in)	:: adimc					!state array
    end subroutine write_sac_state

    subroutine read_snow17_state(cs,tprev)
      use nrtype
      use snow17_sac, only: snow17_state_in_file

      !input variables
      real(sp), intent(out) 			:: tprev				!carry over variable
      real(sp), dimension(:), intent(out)	:: cs					!carry over array

    end subroutine read_snow17_state

    subroutine read_sac_state(uztwc,uzfwc,lztwc,lzfsc,lzfpc,adimc)
      use nrtype
      use snow17_sac, only: sac_state_in_file

      !input variables
      real(sp), intent(out)	:: uztwc					!state variable
      real(sp), intent(out)	:: uzfwc					!state array
      real(sp), intent(out)	:: lztwc					!state array
      real(sp), intent(out)	:: lzfsc					!state array
      real(sp), intent(out)	:: lzfpc					!state array
      real(sp), intent(out)	:: adimc					!state array
    end subroutine read_sac_state

    subroutine write_uh_state(route_tci,old_uh,sim_length,uh_length)
      use nrtype
      use snow17_sac, only: uh_state_out_file

      !input variables
      integer(I4B), intent(in)				:: uh_length
      integer(I4B), intent(in)				:: sim_length
      real(sp), dimension(:),intent(in) 		:: old_uh
      real(sp), dimension(:),intent(in)			:: route_tci
    end subroutine write_uh_state

    subroutine read_uh_state(uh_flow,uh_length,sim_length)
      use nrtype
      use snow17_sac, only: uh_state_in_file

      !input variable
      integer(I4B), intent(in)				:: uh_length,sim_length

      !output variables
      real(sp), dimension(:), intent(out) 			:: uh_flow
    end subroutine read_uh_state


    subroutine read_areal_forcing(year,month,day,hour,tmin,tmax,vpd,dayl,swdown,precip)
      use nrtype
      use snow17_sac, only: forcing_name, start_year,start_day,start_month, &
                        end_year,end_month,end_day,lat,area_basin,elev

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
    end subroutine read_areal_forcing


    subroutine read_sac_params(param_name)
      use nrtype
      use snow17_sac, only: uztwm,uzfwm,uzk,pctim,adimp,zperc,rexp, &
			    lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,riva,side,rserv

      !input variables
      character(len=1024),intent(in)	:: param_name

    end subroutine read_sac_params

    subroutine read_snow17_params(param_name)
      use nrtype
      use snow17_sac, only: scf,mfmax,mfmin,uadj,si,pxtemp,nmf,&
                        tipm,mbase,plwhc,daygm,adc

      !input variables
      character(len=1024),intent(in)	:: param_name

    end subroutine read_snow17_params

    subroutine read_uhp_params(param_name)
      use nrtype
      use snow17_sac, only: unit_shape,unit_scale,pet_coef

      !input variables
      character(len=1024),intent(in)	:: param_name

    end subroutine read_uhp_params


  end interface

end module gauge_calib