import sacsma
import pandas as pd
import numpy as np
from tqdm import tqdm
import pet
from collections import OrderedDict


def run_sacsma(dates, forcings, parameters, attributes):

  # Number of time steps
  num_time_steps = forcings.shape[0]

  # Timestep (I believe this is in units 'days', but I am not 100% sure)
  assert type(dates) == pd.Series
  assert type(dates[0]) == pd.Timestamp
  assert all(dates.diff().values[1:] == dates.diff().values[1])
  dt = (dates[1]-dates[0]).total_seconds()/(60*60*24)
 
  # Initial states must be numpy scalars for the fortran intent(inout) to work. 
  # Also must have trailing decimals, or fortran will treat as integers even if typed as real in the interface.
  state_keys = ['uztwc', 'uzfwc', 'lztwc', 'lzfsc', 'lzfpc', 'adimc']
  uztwc = np.array(100.)
  uzfwc = np.array(100.)
  lztwc = np.array(100.)
  lzfsc = np.array(100.)
  lzfpc = np.array(100.)
  adimc = np.array(100.)

  # Calculate potential evaporation
  if not('PET(mm/day)' in forcings):
    forcings['PET(mm/day)'] = pet.get_priestley_taylor_pet(forcings['Tmin(C)'], forcings['Tmax(C)'], forcings['SRAD(W/m2)'], 
                                                           attributes['gauge_lat'], attributes['elev_mean'], dates.dt.dayofyear)    

  # # Init storage
  flux_keys = ['surf', 'grnd', 'tet']
  states = np.full([num_time_steps, len(state_keys)], np.nan)
  fluxes = np.full([num_time_steps, len(flux_keys)], np.nan)

  # Make parameters as numpy array
  parameters_np = parameters[['uztwm', 'uzfwm', 'uzk', 'pctim', 'adimp', 'riva', 'zperc', 'rexp', 'lztwm', 'lzfsm', 'lzfpm', 'lzsk', 'lzpk', 'pfree', 'side']].values.copy()

  # Time loop
  for t in tqdm(range(num_time_steps)):

    # Run the model at one timestep
    current_flux_np = sacsma.fland1(forcings['PRCP(mm/day)'][t], forcings['PET(mm/day)'][t], dt, uztwc, uzfwc, lztwc,
                                    lzfsc, lzfpc, adimc, *parameters_np, 0)

    # Save states & outputs
    states[t] = uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc
    fluxes[t] = current_flux_np

  # Turn into dataframes
  states_df = pd.DataFrame(states, index=dates)
  fluxes_df = pd.DataFrame(fluxes, index=dates)

  return fluxes_df, states_df



