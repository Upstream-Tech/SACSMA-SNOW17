import pandas as pd
import numpy as np
from model.potential_evap import priestley_taylor_pet
import sacsma_source.sac.ex_sac1 as sacsma
import sacsma_source.sac.duamel as unit_hydrograph
import sacsma_source.snow19.exsnow as snow17


def rmse_obj_fun(forcings: pd.DataFrame,
                 parameters: pd.Series,
                 observations: pd.Series,
                 latitude: float,
                 elevation: float):

    sac_fluxes, sac_states = run_sacsma(forcings, parameters, latitude, elevation)

    # concat and remove nans
    df = sac_fluxes['sacsma_uh_qq'].to_frame()
    df = df.join(observations)
    df = df.dropna(axis=0, how='any')

    return np.sqrt(np.mean((df['sacsma_uh_qq'] - df['QObs']) ** 2))


def run_sacsma(forcings: pd.Series,
               parameters: pd.Series,
               latitude: float,
               elevation: float):

    # Timestep in different units
    dt_seconds = int((forcings.index[1] - forcings.index[0]).total_seconds())
    dt_days = dt_seconds / 86400
    dt_hours = dt_seconds / (60 * 60)

    # Keys for parameters, states, fluxes
    sacsma_parameter_keys = [
        'uztwm', 'uzfwm', 'uzk', 'pctim', 'adimp', 'riva', 'zperc', 'rexp', 'lztwm', 'lzfsm', 'lzfpm', 'lzsk', 'lzpk',
        'pfree', 'side', 'rserv'
    ]
    snow17_parameter_keys = [
        'scf', 'mfmax', 'mfmin', 'uadj', 'si', 'nmf', 'tipm',
        'mbase', 'pxtemp', 'plwhc', 'daygm'
    ]
    hydrograph_parameter_keys = [
        'unit_shape', 'unit_scale'
    ]
    state_keys = [
        # 'sacsma_snow17_cs', 
        'sacsma_snow17_tprev',
        'sacsma_uztwc', 'sacsma_uzfwc', 'sacsma_lztwc', 'sacsma_lzfsc', 'sacsma_lzfpc', 'sacsma_adimc'
    ]
    flux_keys = [
        'sacsma_pet',
        'sacsma_snow17_raim', 'sacsma_snow17_sneqv', 'sacsma_snow17_snow', 'sacsma_snow17_snowh',
        'sacsma_surf', 'sacsma_grnd', 'sacsma_qq', 'sacsma_tet'
    ]

    # Set constant snow17 parameter
    adc = np.array([0.05, 0.15, 0.26, 0.45, 0.5, 0.56, 0.61, 0.65, 0.69, 0.82, 1.0]).astype('f4')

    # Initial states must be numpy scalars for the fortran intent(inout) to work.
    # Also must have trailing decimals, or fortran will treat as integers even if typed as real in the pyf interface.
    uztwc = np.array(200.)
    uzfwc = np.array(200.)
    lztwc = np.array(200.)
    lzfsc = np.array(200.)
    lzfpc = np.array(200.)
    adimc = np.array(200.)

    tprev = np.array(0.)
    cs = np.full(19, 0., dtype='f4')

    # Extract parameters as numpy array
    sacsma_parameters_np = parameters[sacsma_parameter_keys].values.copy()
    hydrograph_parameters_np = parameters[hydrograph_parameter_keys].values.copy()
    snow17_parameters_np = parameters[snow17_parameter_keys].values.copy()

    # Extract vectors as numpy arrays for speed
    precipitation = forcings['PRCP(mm/day)'].values
    temperature = 0.5 * (forcings['Tmax(C)'].values + forcings['Tmin(C)'].values)
    day = forcings['Day'].values
    month = forcings['Mnth'].values
    year = forcings['Year'].values
    # latitude = attributes['gauge_lat'].astype('f4')
    # elevation = attributes['elev_mean'].astype('f4')

    # Calculate potential evaporation as numpy array
    pet = priestley_taylor_pet(forcings['Tmin(C)'], forcings['Tmax(C)'], forcings['SRAD(W/m2)'],
                                                  latitude, elevation,
                                                  forcings.index.to_series().dt.dayofyear)

    # estimate surface pressure [hPa] as a function of elevation only
    surf_pres = calc_surface_pressure(elevation)

    # Init storage as numpy arrays
    states = np.full([forcings.shape[0], len(state_keys)], np.nan)
    fluxes = np.full([forcings.shape[0], len(flux_keys)], np.nan)

    # Time loop
    for t in range(forcings.shape[0]):

        # Run the model at one timestep
        raim, sneqv, snow, snowh = snow17.exsnow19(dt_seconds, dt_hours, day[t], month[t], year[t],
                                                   precipitation[t], temperature[t],
                                                   latitude,
                                                   *snow17_parameters_np,
                                                   elevation,
                                                   surf_pres,
                                                   adc,
                                                   cs, tprev)

        surf, grnd, qq, tet = sacsma.exsac(dt_seconds, raim, temperature[t], pet[t],
                                           *sacsma_parameters_np,  # 0.,
                                           uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)

        # Save states & outputs in time arrays
        states[t] = tprev, uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc
        fluxes[t] = pet[t], raim, sneqv, snow, snowh, surf, grnd, qq, tet

    # channel routing
    m_unit_hydro = 1000
    n_unit_hydro = forcings.shape[0] + m_unit_hydro
    hydrograph_qq = unit_hydrograph.duamel(fluxes[:, 7], *hydrograph_parameters_np, dt_days, n_unit_hydro, m_unit_hydro,
                                           1, 0)

    # Turn into dataframes
    states_df = pd.DataFrame(states, index=forcings.index, columns=state_keys)
    fluxes_df = pd.DataFrame(fluxes, index=forcings.index, columns=flux_keys)
    fluxes_df.insert(2, "sacsma_uh_qq", hydrograph_qq[:-m_unit_hydro], True)

    return fluxes_df, states_df # hydrograph_qq[:-m_unit_hydro] #


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
