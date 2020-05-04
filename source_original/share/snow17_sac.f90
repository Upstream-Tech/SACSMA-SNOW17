module snow17_sac
  use nrtype

  implicit none

!initialization variables
!these are in the namelist bloc &INIT_CONTROL
  integer(I4B)			:: dt			!model time step (seconds)
  character(len = 1024) 	:: forcing_name		!name of forcing data file 
  character(len = 1024) 	:: stream_name		!name of observed streamflow data
  character(len = 1024) 	:: model_out		!base name for output files
  character(len = 1024) 	:: sac_param_file		!name for sac-sma parameters
  character(len = 1024) 	:: snow17_param_file		!name for snow17 parameters
  character(len = 1024) 	:: uhp_param_file		!name for uh and pet parameters
  character(len = 1024)		:: uh_state_out_file		!name for uh state output file
  character(len = 1024)		:: uh_state_in_file		!name for uh state input file
  character(len = 1024)		:: snow17_state_out_file	!name for snow17 state output file
  character(len = 1024)		:: sac_state_out_file	!name for sac state output file
  character(len = 1024)		:: snow17_state_in_file	!name for snow17 state input file
  character(len = 1024)		:: sac_state_in_file	!name for sac state input file
  integer(I4B) 			:: gage_id		!usgs gage id 
  integer(I4B)			:: start_month		!starting month 
  integer(I4B)			:: start_day		!starting day
  integer(I4B)			:: start_year		!starting year
  integer(I4B)			:: end_month		!ending month 
  integer(I4B)			:: end_day		!ending day
  integer(I4B)			:: end_year		!ending year
  integer(I4B)			:: restart_run		!restart run flag
  integer(I4B)			:: write_restart	!flag to tell current run to write restart files
  real(sp)			:: in_swe
  real(sp)			:: in_uztwc
  real(sp)			:: in_uzfwc
  real(sp)			:: in_lztwc
  real(sp)			:: in_lzfsc
  real(sp)			:: in_lzfpc
  real(sp)			:: in_adimc


!variables for snow17,pet,streamflow calculations
!read in from the forcing file
  real(dp) :: lat
  real(dp) :: elev
  real(dp) :: area_basin


!SAC_model params
!in the param file for SAC_SMA
  real(sp)	:: uztwm,uzfwm,uzk,pctim,adimp,zperc,rexp
  real(sp)	:: lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree
  real(sp)	:: riva,side,rserv

!uh and pet
  real(sp)	:: pet_coef		!coefficient for p-t pet calculation
  real(sp)	:: unit_shape,unit_scale  !unit hydrograph parameters

!Snow17_model params
!in the paramfile for SNOW_17
  real(sp)	:: scf,mfmax,mfmin,uadj,si,pxtemp
  real(sp)	:: nmf,tipm,mbase,plwhc,daygm
  real(sp), dimension(11)  :: adc

!namelists

!namelist for snow_17
!right now read in the following parameters with upper and lower bounds
!SCF, MFMAX,MFMIN,UADJ,SI, Areal depletion curve info, and PXTMP
!also read in: nmf,tipm,mbase,plwhc,daygm
!  namelist / SNOW_17 / scf,mfmax,mfmin,uadj,si,adc,nmf,tipm,&
!                       pxtemp,mbase,plwhc,daygm


!namelist for SAC-SMA
!read in the following parameters with upper and lower bounds
!UZTWM,UZFWM,UZK,PCTIM,ADIMP,ZPERC,REXP,LZTWM,LZFSM,LZFPM,LZSK,LZPK,PFREE
!  namelist / SAC_SMA / uztwm,uzfwm,uzk,pctim,adimp,zperc,rexp,&
!                       lztwm,lzfsm,lzfpm,lzsk,lzpk,pfree,riva,side,rserv,&
!                       unit_shape,unit_scale

  namelist / INIT_CONTROL / dt,forcing_name,stream_name,model_out,gage_id, &
                            start_day,start_month,start_year,end_year,end_month, &
                            end_day,in_swe,in_uztwc,in_uzfwc,in_lztwc,in_lzfsc, &
                            in_lzfpc,in_adimc,sac_param_file,snow17_param_file,uhp_param_file, &
			    uh_state_in_file,restart_run,write_restart, &
			    snow17_state_out_file,sac_state_out_file,snow17_state_in_file, &
			    sac_state_in_file,uh_state_out_file

  save
end module
