import pandas as pd
import numpy as np
import potential_evap, sacsma, camels_utilities
from tqdm import tqdm

def run_sacsma(inputs: pd.DataFrame, gauge_id: str):

  # Load parameters
  parameters = camels_utilities.load_sacsma_parameters(gauge_id)

  # Load attributes
  attributes = camels_utilities.load_basin_attributes(gauge_id)

  # Timestep in [seconds]
  dt = int((inputs.index[1] - inputs.index[0]).total_seconds())

  # Keys for parameters, states, fluxes
  parameter_keys = [
      'uztwm', 'uzfwm', 'uzk', 'pctim', 'adimp', 'riva', 'zperc', 'rexp', 'lztwm', 'lzfsm', 'lzfpm', 'lzsk', 'lzpk',
      'pfree', 'side'
  ]
  state_keys = ['sacsma_uztwc', 'sacsma_uzfwc', 'sacsma_lztwc', 'sacsma_lzfsc', 'sacsma_lzfpc', 'sacsma_adimc']
  flux_keys = ['sacsma_pet', 'sacsma_surf', 'sacsma_grnd', 'sacsma_tet']

  # Initial states must be numpy scalars for the fortran intent(inout) to work.
  # Also must have trailing decimals, or fortran will treat as integers even if typed as real in the pyf interface.
  uztwc = np.array(10.)
  uzfwc = np.array(10.)
  lztwc = np.array(10.)
  lzfsc = np.array(10.)
  lzfpc = np.array(10.)
  adimc = np.array(10.)

  # Extract parameters as numpy array
  parameters_np = parameters[parameter_keys].values.copy()

  # Extract precipitation as numpy array
  precipitation = inputs['PRCP(mm/day)'].values

  # Calculate potential evaporation as numpy array
  pet = potential_evap.get_priestley_taylor_pet(inputs['Tmin(C)'], inputs['Tmax(C)'], inputs['SRAD(W/m2)'], attributes['gauge_lat'],
                                                attributes['elev_mean'],
                                                inputs.index.to_series().dt.dayofyear)

  # Init storage as numpy arrays
  states = np.full([inputs.shape[0], len(state_keys)], np.nan)
  fluxes = np.full([inputs.shape[0], len(flux_keys)], np.nan)

  # Time loop
  for t in tqdm(range(inputs.shape[0])):

    # Run the model at one timestep
    surf, grnd, tet = sacsma.fland1(precipitation[t], pet[t], dt, uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc,
                                    *parameters_np, 0)

    # Save states & outputs in time arrays
    states[t] = uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc
    fluxes[t] = pet[t], surf, grnd, tet

  # Turn into dataframes
  states_df = pd.DataFrame(states, index=inputs.index, columns=state_keys)
  fluxes_df = pd.DataFrame(fluxes, index=inputs.index, columns=flux_keys)

  return fluxes_df, states_df
