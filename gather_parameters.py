import pickle as pkl 
import numpy as np 
import camels_utilities
from tqdm import tqdm
import pandas as pd 

# Location of camels data
camels_dir = '/Users/grey/workspace/camels_data/'

# Which forcings
# forcing_type = 'nldas'
forcings = 'maurer'
#forcings = 'daymet'

# Load list of gauge IDs
with open('camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Grab parameter names
test_parameters = camels_utilities.load_sacsma_parameters(camels_dir, forcing_type, basins[0])
param_names = test_parameters.keys()

# Init storage
parameters = pd.DataFrame(index=basins, columns=param_names)

# Loop through all basins
for basin in tqdm(basins):
     basin_params = camels_utilities.load_sacsma_parameters(camels_dir, forcing_type, basin)
     parameters.loc[basin] = basin_params

# Pickle
fname = f'results/all_camels_parameters_{forcing_type}.pkl'
with open(fname, 'wb') as f:
    pkl.dump(f, parameters)


