
!ccccccccccccccccccccccccccccccccc

subroutine calc_pet_pt(jday,tmax,tmin,tair,vpd,swdown,dayl,pet)
  use constants
  use nrtype
  use gauge_calib, only: sfc_pressure
  use snow17_sac, only: elev,lat,pet_coef

  implicit none

  !input variable
  integer(I4B), dimension(:), intent(in) 	:: jday    !day of year (1-365(6))
  real(dp), dimension(:), intent(in) 		:: tair    !daily average air temperature (deg C)
  real(dp), dimension(:), intent(in) 		:: tmax    !daily max air temperature (deg C)
  real(dp), dimension(:), intent(in) 		:: tmin    !daily min air temperature (deg C)
  real(dp), dimension(:), intent(in) 		:: vpd     !daily average vapor pressure (Pa)
  real(dp), dimension(:), intent(in) 		:: swdown  !daily average swdown (w m-2)
  real(dp), dimension(:), intent(in) 		:: dayl    !length of day (sun above horizon) seconds

  !output variable
  real(dp),dimension(:),intent(out)		:: pet     ! estimated PET (mm/day)

!local variables
  integer(I4B)          :: i
  integer(I4B)		:: end_pt
  
  real(DP)               :: albedo	!albedo for pet calculation
  real(DP)             	 :: apt		!p-t coefficient for aird regions...

  real(DP)               :: l		!latent heat of vaporization (MJ kg-1)
  real(DP)               :: g		!psychrometric constant
  real(DP)               :: s		!slope of the saturation vapour pressure-temperature relationship
  real(DP)               :: tavg	!daily average temperature ( deg C )
  real(DP)               :: r_net	!daily net radiation (MJ m-2 day-1)
  real(DP)               :: r_nl	!daily net longwave (MJ m-2 day-1)
  real(DP)               :: r_ns	!daily net shortwave (MJ m-2 day-1)
  real(DP)               :: pressure	!surface pressure (kPa)
  real(DP)               :: r_s		!daily estimated solar from forcing dataset (MJ m-2 day-1)
  real(DP)               :: r_a		!daily extraterrestrial radiation (MJ m-2 day-1)
  real(DP)               :: r_so	!daily clear sky solar shortwave (MJ m-2 day-1)
  real(DP)               :: d_r		!inverse relative earth-sun distance
  real(DP)               :: dec		!solar declination in radians
  real(DP)               :: lat_rad	!latitude in radians
  real(DP)               :: sha		!sunset hour angle
  real(DP)               :: e_a		!vapor pressure (kPa)
  real(DP)               :: e_s		!saturation vapor pressure (kPa)


!!set pet coefficient from namelist parameter
  apt = real(pet_coef,kind(dp))

!!set albedo to a constant
  albedo = 0.20_dp

!!calculate pressure from elevation using standard atmosphere (taken from Snow-17)
  call sfc_pressure(elev,pressure) !pressure in hPa
  pressure = pressure/10.0_dp !in kPa



!how many days to do?
  end_pt = size(jday)

!now lets get going on p-t pet
  do i = 1,end_pt
!    tavg = (tmax(i)+tmin(i))/2
    tavg = tair(i)

    l = l_v - tadj*tavg   !(MJ kg-1)
    s = slope_svpc_a*exp(slope_svpc_b*tavg)
    g = (c_p*pressure)/(e*l)

    e_s = e_sa*exp((e_sb*tavg)/(tavg+e_sc)) !in kPa
    e_a = vpd(i)/1000.0_dp

    !radiation terms
    d_r = 1 + sun_e_inv_d*cos(((2*pi_d)/365.0_dp) * jday(i))
    dec = sol_dec_a*sin( (((2*pi_d)/365.0_dp) * jday(i)) - sol_dec_b)
    lat_rad = (pi_d/180.0_dp)*lat
    sha = acos(-tan(lat_rad)*tan(dec))

    r_a = ((24.*60.)/pi_d)*gsc*d_r*((sha*sin(lat_rad)*sin(dec)) + (cos(lat_rad)*cos(dec)*sin(sha))) !in MJ m-2 day-1
    r_so = (clear_sky_a + clear_sky_b*elev)*r_a !in MJ m-2 day-1
    r_s = (swdown(i)*dayl(i)/86400.)*w_to_mj !in MJ m-2 day-1
    r_ns = (1.0_dp-albedo)*r_s  !assume a constant albedo of ~0.20    !in MJ m-2 day-1

    r_nl = sbc*((((tmax(i)+273.16)**4)+((tmin(i)+273.16)**4))/2)*(net_lw_a-net_lw_b*sqrt(e_a))*(net_lw_c*(r_s/r_so)-net_lw_d) !in MJ m-2 day-1

    r_net = r_ns - r_nl


    pet(i) = (1./l)*((s*r_net)/(s+g))*apt

  enddo

  return
end subroutine
