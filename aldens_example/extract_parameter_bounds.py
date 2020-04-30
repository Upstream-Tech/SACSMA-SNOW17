import pickle as pkl
import pandas as pd
import numpy as np
from tqdm import tqdm

with open('results/all_camels_parameters_nldas.pkl', 'rb') as f:
  parameters = pkl.load(f)

basins = list(parameters.keys())
#parm_names = list(parameters[basins[0]].index)
parameter_names = ['uztwm', 'uzfwm', 'uzk', 'pctim', 'adimp', 'riva', 'zperc', 'rexp', 'lztwm', 'lzfsm', 'lzfpm', 'lzsk', 'lzpk', 'pfree', 'side']

parm_bounds = pd.DataFrame(np.nan, index=parm_names, columns=['min', 'max'])

for parm in parm_names:
  print(parm)
  for basin in tqdm(basins):
    basin_parm_min = np.min(parameters[basin].loc[parm].values.astype(float))
    basin_parm_max = np.max(parameters[basin].loc[parm].values.astype(float))
    if np.isnan(parm_bounds.loc[parm, 'min']):
      parm_bounds.loc[parm, 'min'] = basin_parm_min 
      parm_bounds.loc[parm, 'max'] = basin_parm_max 
    else:
      parm_bounds.loc[parm, 'min'] = np.min([parm_bounds.loc[parm, 'min'], basin_parm_min]) 
      parm_bounds.loc[parm, 'max'] = np.max([parm_bounds.loc[parm, 'max'], basin_parm_max]) 



print(parm_bounds.loc[parameter_keys])

