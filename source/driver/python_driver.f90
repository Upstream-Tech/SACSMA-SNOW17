 subroutine timestep_driver(dt, year, month, day,                               & ! timestep
    elev, lat                                                                   & ! characteristics                      
    PRECIP, TMAX, TMIN, VPD, SRAD, PA,                                          & ! forcings
    SFC, MFMAX, MFMIN, UADJ, SI, NMF, TIPM, MBASE, PXTEMP, PLWHC, DAYGM, ADC,   & ! snow17 parameters
    UZTWM, UZFWM, UZK, PCTIM, ADIMP, RIVA, ZPERC, REXP, LZTWM, LZFSM, LZFPM,    & ! sac parameters
    LZSK, LSPK, PFREE, SIDE, RSERV,                                             & ! sac parameters
    cs, tprev,                                                                  & ! snow17 states
    uztwc_sp,uzfwc_sp,lztwc_sp,lzfsc_sp,lzfpc_sp,adimc_sp,                      & ! sac states
    raim, sneqv, snow, snowh,                                                   & ! snow17 outputs
    qs, qg, tci, eta)                                                             ! sac outputs

    ! Day of Year
    call julian_day(year, month, day, jday)

    ! Compute PET
    tair = (tmax+tmin)/2.0_dp
    call calc_pet_pt(jday, tmax, tmin, tair, vpd, srad, dayl, pet)

    ! Compute surface pressure
    call sfc_pressure(elev, pa)

    ! Set single precision inputs
    tair_sp   = real(tair(i),kind(sp))
    precip_sp = real(precip(i),kind(sp))
    pet_sp    = real(pet(i),kind(sp))

    ! Call snow-17
    call exsnow19(int(dt), int(dt/sec_hour), day, month, year,                  &
               precip_sp, tair_sp,                                              & ! SNOW17 INPUTS
               raim, sneqv, snow, snowh,                                        & ! SNOW17 OUTPUTS
               real(lat,kind(sp)), scf, mfmax, mfmin, uadj, si, nmf,            & ! SNOW17 PARAMETERS
               tipm, mbase, pxtemp, plwhc, daygm, real(elev,kind(sp)),          & ! SNOW17 PARAMETERS
               real(pa,kind(sp)), adc,                                          & ! SNOW17 PARAMETERS
	            cs, tprev)                                                         ! SNOW17 CARRYOVER VARIABLES

    ! Call SAC-SMA           
    call exsac(1, real(dt), raim, tair_sp, pet_sp,                              &
	            uztwm, uzfwm, uzk, pctim, adimp, riva, zperc,                    & ! SAC PARAMETERS
	            rexp, lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree,                    & ! SAC PARAMETERS
	            side, rserv,                                                     & ! SAC PARAMETERS
               uztwc_sp, uzfwc_sp, lztwc_sp, lzfsc_sp, lzfpc_sp, adimc_sp,      & ! SAC State variables
	            qs, qg, tci, eta)                                                  ! SAC OUTPUTS
