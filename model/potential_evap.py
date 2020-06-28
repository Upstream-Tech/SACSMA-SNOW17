import numpy as np

LAMBDA = 2.45  # Kept constant, MJkg-1
ALPHA = 1.26  # Calibrated in CAMELS, here static
STEFAN_BOLTZMAN = 4.903e-09
PI = np.pi


# Routine plucked directly from NCAR Fortran code
# I don't know what any of these constants are
def calc_surface_pressure(elev):
  # constants
  sfc_pres_a = 33.86
  sfc_pres_b = 29.9
  sfc_pres_c = 0.335
  sfc_pres_d = 0.00022
  sfc_pres_e = 2.4

  # sfc pres in hPa
  sfc_pres = sfc_pres_a * (sfc_pres_b - (sfc_pres_c * (elev / 100)) + (sfc_pres_d * ((elev / 100) ** sfc_pres_e)))

  return sfc_pres


def priestley_taylor_pet(t_min: np.ndarray, t_max: np.ndarray, s_rad: np.ndarray, lat: float, elev: float,
                             doy: np.ndarray) -> np.ndarray:
  """Calculate PET as an approximation following the Priestley-Taylor equation.

    The ground head flux G is assumed to be 0 at daily time steps (see Newman et al. 2015). The
    equation follow FAO-56 (Allen et al. (1998))

    Parameters
    ----------
    t_min : np.ndarray
        Daily min temperature (degree C)
    t_max : np.ndarray
        Daily max temperature (degree C)
    s_rad : np.ndarray
        Solar radiation (Wm-2)
    lat : float
        Latitude in degree
    elev : float
        Elevation in m
    doy : np.ndarray
        Day of the year

    Returns
    -------
    np.ndarray
        Array containing PET estimates in mm/day
    """
  lat = lat * (PI / 180)  # degree to rad

  # Slope of saturation vapour pressure curve
  t_mean = 0.5 * (t_min + t_max)
  slope_svp = get_slope_svp_curve(t_mean)

  # incoming netto short-wave radiation
  s_rad = s_rad * 0.0864  # conversion Wm-2 -> MJm-2day-1
  in_sw_rad = get_net_sw_srad(s_rad)

  # outgoginng netto long-wave radiation
  sol_dec = get_sol_decl(doy)
  sha = get_sunset_hour_angle(lat, sol_dec)
  ird = get_ird_earth_sun(doy)
  et_rad = get_extraterra_rad(lat, sol_dec, sha, ird)
  cs_rad = get_clear_sky_rad(elev, et_rad)
  a_vp = get_avp_tmin(t_min)
  out_lw_rad = get_net_outgoing_lw_rad(t_min, t_max, s_rad, cs_rad, a_vp)

  # net radiation
  net_rad = get_net_rad(in_sw_rad, out_lw_rad)

  # gamma
  atm_pressure = get_atmos_pressure(elev)
  gamma = get_psy_const(atm_pressure)

  # PET MJm-2day-1
  pet = (ALPHA / LAMBDA) * (slope_svp * net_rad) / (slope_svp + gamma)

  # convert energy to evap
  pet = pet * 0.408

  return pet


def get_slope_svp_curve(t_mean: np.ndarray) -> np.ndarray:
  """Slope of saturation vapour pressure curve

    Equation 13 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    t_mean : np.ndarray
        Mean temperature (degree C)

    Returns
    -------
    np.ndarray
        Slope of the saturation vapor pressure curve in kPa/(degree C)
    """
  delta = 4098 * (0.6108 * np.exp((17.27 * t_mean) / (t_mean + 237.3))) / ((t_mean + 237.3)**2)
  return delta


def get_net_sw_srad(s_rad: np.ndarray, albedo: float = 0.23) -> np.ndarray:
  """Calculate net shortwave radiation

    Equation 38 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    s_rad : np.ndarray
        Incoming solar radiation (MJm-2day-1)
    albedo : float, optional
        Albedo, by default 0.23

    Returns
    -------
    np.ndarray
        Net shortwave radiation (MJm-2day-1)
    """
  net_srad = (1 - albedo) * s_rad
  return net_srad


def get_sol_decl(doy: np.ndarray) -> np.ndarray:
  """Get solar declination

    Equation 24 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    doy : np.ndarray
        Day of the year

    Returns
    -------
    np.ndarray
        Solar declination in rad
    """
  # equation 24 FAO Allen
  sol_dec = 0.409 * np.sin((2 * np.pi) / 365 * doy - 1.39)
  return sol_dec


def get_sunset_hour_angle(lat: float, sol_dec: np.ndarray) -> np.ndarray:
  """Sunset hour angle



    Parameters
    ----------
    lat : float
        Latitude in rad
    sol_dec : np.ndarray
        Solar declination in rad

    Returns
    -------
    np.ndarray
        Sunset hour angle in rad
    """
  term = -np.tan(lat) * np.tan(sol_dec)
  term[term < -1] = -1
  term[term > 1] = 1
  sha = np.arccos(term)
  return sha


def get_ird_earth_sun(doy: np.ndarray) -> np.ndarray:
  """Inverse relative distance between Earth and Sun

    Equation 23 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    doy : np.ndarray
        Day of the year

    Returns
    -------
    np.ndarray
        Inverse relative distance between Earth and Sun
    """
  ird = 1 + 0.033 * np.cos((2 * PI) / 365 * doy)
  return ird


def get_extraterra_rad(lat: float, sol_dec: np.ndarray, sha: np.ndarray, ird: np.ndarray) -> np.ndarray:
  """Extraterrestrial Radiation

    Equation 21 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    lat : float
        Lat in rad (pos for northern hemisphere)
    sol_dec : np.ndarray
        Solar declination in rad
    sha : np.ndarray
        Sunset hour angle in rad
    ird : np.ndarray
        Inverse relative distance of Earth and Sun

    Returns
    -------
    np.ndarray
        Extraterrestrial radiation MJm-2day-1
    """
  term1 = (24 * 60) / PI * 0.082 * ird
  term2 = sha * np.sin(lat) * np.sin(sol_dec) + np.cos(lat) * np.cos(sol_dec) * np.sin(sha)
  et_rad = term1 * term2
  return et_rad


def get_clear_sky_rad(elev: float, et_rad: np.ndarray) -> np.ndarray:
  """Clear sky radiation

    Equation 37 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    elev : float
        Elevation in m
    et_rad : np.ndarray
        Extraterrestrial radiation in MJm-2day-1

    Returns
    -------
    np.ndarray
        Clear sky radiation MJm-2day-1
    """
  cs_rad = (0.75 + 2 * 10e-5 * elev) * et_rad
  return cs_rad


def get_avp_tmin(t_min: np.ndarray) -> np.ndarray:
  """Actual vapor pressure estimated using min temperature

    Equation 48 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    t_min : np.ndarray
        Minimum temperature in degree C

    Returns
    -------
    np.ndarray
        Actual vapor pressure kPa
    """
  avp = 0.611 * np.exp((17.27 * t_min) / (t_min + 237.3))
  return avp


def get_net_outgoing_lw_rad(t_min: np.ndarray, t_max: np.ndarray, s_rad: np.ndarray, cs_rad: np.ndarray,
                            a_vp: np.ndarray) -> np.ndarray:
  """Net outgoing longwave radiation

    Expects temperatures in degree and does the conversion in kelvin in the function.

    Equation 49 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    t_min : np.ndarray
        Min temperature in degree C
    t_max : np.ndarray
        Max temperature in degree C
    s_rad : np.ndarray
        Measured or modeled solar radiation MJm-2day-1
    cs_rad : np.ndarray
        Clear sky radiation MJm-2day-1
    a_vp : np.ndarray
        Actuatal vapor pressure kPa

    Returns
    -------
    np.ndarray
        Net outgoing longwave radiation MJm-2day-1
    """
  term1 = ((t_max + 273.16)**4 + (t_min + 273.16)**4) / 2  # conversion in K in equation
  term2 = 0.34 - 0.14 * np.sqrt(a_vp)
  term3 = 1.35 * s_rad / cs_rad - 0.35
  net_lw = STEFAN_BOLTZMAN * term1 * term2 * term3
  return net_lw


def get_net_rad(sw_rad: np.ndarray, lw_rad: np.ndarray) -> np.ndarray:
  """Net radiation

    Equation 40 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    sw_rad : np.ndarray
        Net incoming shortwave radiation MJm-2day-1
    lw_rad : np.ndarray
        Net outgoing longwave radiation MJm-2day-1

    Returns
    -------
    np.ndarray
        [description]
    """
  return sw_rad - lw_rad


def get_atmos_pressure(elev: float) -> float:
  """Atmospheric pressure

    Equation 7 FAO-56 Allen et al. (1998)

    Parameters
    ----------
    elev : float
        Elevation in m

    Returns
    -------
    float
        Atmospheric pressure in kPa
    """
  temp = (293.0 - 0.0065 * elev) / 293.0
  return np.power(temp, 5.26) * 101.3


def get_psy_const(atm_pressure: float) -> float:
  """Psychometric constant

    Parameters
    ----------
    atm_pressure : float
        Atmospheric pressure in kPa

    Returns
    -------
    float
        Psychometric constant in kPa/(degree C)
    """
  return 0.000665 * atm_pressure


def srad_from_t(et_rad, cs_rad, t_min, t_max, coastal=False):
  # equation 50
  if coastal:
    adj = 0.19
  else:
    adj = 0.16

  sol_rad = adj * np.sqrt(t_max - t_min) * et_rad

  return np.minimum(sol_rad, cs_rad)
