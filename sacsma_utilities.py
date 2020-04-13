import sacsma
import pandas as pd
import numpy as np
from tqdm import tqdm


def run_sacsma(dates, forcings, parameters):

  # Number of time steps
  num_time_steps = forcings.shape[0]

  # Timestep (I believe this is in units 'days', but I am not 100% sure)
  dt = (dates[1]-dates[0]).total_seconds()/(60*60*24)
 
  # Initial states
  # states = {'uztwc': 100., 'uzfwc': 100., 'lztwc': 100., 'lzfsc': 10., 'lzfpc': 10., 'adimc': 100.}
  uztwc = np.array(100.)
  uzfwc = np.array(100.)
  lztwc = np.array(100.)
  lzfsc = np.array(100.)
  lzfpc = np.array(100.)
  adimc = np.array(100.)

  # Init storage
  states = pd.DataFrame(index=dates, columns=['uztwc', 'uzfwc', 'lztwc', 'lzfsc', 'lzfpc', 'adimc'])
  outputs = pd.DataFrame(index=dates, columns=['surf', 'grnd', 'tet'])

  # Time loop
  for t in tqdm(range(num_time_steps)):

    # Run
    surf, grnd, tet = sacsma.fland1(forcings['PRCP(mm/day)'][t], forcings['PET(mm/day)'][t], dt, uztwc, uzfwc, lztwc, lzfsc, lzfpc,
                                    adimc, parameters['uztwm'], parameters['uzfwm'], parameters['uzk'],
                                    parameters['pctim'], parameters['adimp'], parameters['riva'], parameters['zperc'],
                                    parameters['rexp'], parameters['lztwm'], parameters['lzfsm'], parameters['lzfpm'],
                                    parameters['lzsk'], parameters['lzpk'], parameters['pfree'], parameters['side'], 0)

    # Save states
    states.loc[forcings['Date'][t], 'uztwc'] = uztwc
    states.loc[forcings['Date'][t], 'uzfwc'] = uzfwc
    states.loc[forcings['Date'][t], 'lztwc'] = lztwc
    states.loc[forcings['Date'][t], 'lzfsc'] = lzfsc
    states.loc[forcings['Date'][t], 'lzfpc'] = lzfpc
    states.loc[forcings['Date'][t], 'adimc'] = adimc

    # Save outputs
    outputs.loc[forcings['Date'][t], 'surf'] = surf
    outputs.loc[forcings['Date'][t], 'grnd'] = grnd
    outputs.loc[forcings['Date'][t], 'tet'] = tet

  return outputs, states



