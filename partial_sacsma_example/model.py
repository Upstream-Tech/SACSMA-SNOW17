import pandas as pd
import numpy as np
from potential_evap import get_priestley_taylor_pet
from sacsma_utilities import run_sacsma


class model(object):

    def __init__(self, 
                 forcings: pd.DataFrame, 
                 latitude: float, 
                 elevation: float, 
                 observations: pd.DataFrame, 
                 warmup_timesteps: int):

      # Calculate potential evaporation as numpy array
      forcings['PET(mm/day)'] = get_priestley_taylor_pet(forcings['Tmin(C)'], 
                                                        forcings['Tmax(C)'], 
                                                        forcings['SRAD(W/m2)'], 
                                                        latitude,
                                                        elevation,
                                                        forcings.index.to_series().dt.dayofyear)

      # Store data
      merged_df = pd.merge(forcings, observations, on='Date')
      self.precip = merged_df['PRCP(mm/day)'].values
      self.pet = merged_df['PET(mm/day)'].values
      self.observations = merged_df['OBS_RUN'].values

      # Dates for objective function
      self.dates = forcings.index
      self.eval_dates = merged_df.index
      #self.eval_dates = merged_df.iloc[warmup_timesteps:].index
      self.dt_seconds = int((forcings.index[1] - forcings.index[0]).total_seconds())

                       
    def run(self, args):
        return self._run(*args)


    def _run(self,
             uztwm=None,
             uzfwm=None,
             uzk=None,
             pctim=None,
             adimp=None,
             riva=None,
             zperc=None,
             rexp=None,
             lztwm=None,
             lzfsm=None,
             lzfpm=None,
             lzsk=None,
             lzpk=None,
             pfree=None,
             side=None,
             uztwc=None,
             uzfwc=None,
             lztwc=None,
             lzfsc=None,
             lzfpc=None,
             adimc=None):

      # Initial states
      uztwc = np.array(uztwc).astype(float)
      uzfwc = np.array(uzfwc).astype(float) 
      lztwc = np.array(lztwc).astype(float)
      lzfsc = np.array(lzfsc).astype(float)
      lzfpc = np.array(lzfpc).astype(float)
      adimc = np.array(adimc).astype(float)

      # Parameter data type
      uztwm = np.array(uztwm).astype(float)
      uzfwm = np.array(uzfwm).astype(float) 
      uzk   = np.array(uzk).astype(float)
      pctim = np.array(pctim).astype(float)
      adimp = np.array(adimp).astype(float)
      riva  = np.array(riva).astype(float)
      zperc = np.array(zperc).astype(float)
      rexp  = np.array(rexp).astype(float)
      lztwm = np.array(lztwm).astype(float)
      lzfsm = np.array(lzfsm).astype(float)
      lzfpm = np.array(lzfpm).astype(float)
      lzsk  = np.array(lzsk).astype(float)
      lzpk  = np.array(lzpk).astype(float)
      pfree = np.array(pfree).astype(float)
      side  = np.array(side).astype(float) 

      # Init storage 
      qq = np.full(self.precip.shape[0], np.nan)

      # Time loop
      for t in range(self.precip.shape[0]):
        surf, grnd, tet = sacsma.fland1(self.precip[t], self.pet[t], self.dt_seconds, 
                                        uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc,
                                        uztwm, uzfwm, uzk, pctim, adimp, riva, zperc, rexp, lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree, side, 0)
        qq[t] = surf + grnd
        #print(uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)
        assert not np.any(np.isnan([surf, grnd, tet, uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc]))

      return qq 




