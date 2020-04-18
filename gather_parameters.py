import pickle as pkl 
import numpy as np 
import camels_utilities
from tqdm import tqdm
import pandas as pd 

# Location of camels data
camels_dir = '/Users/grey/workspace/camels_data/'

# Which forcings
forcing_type = 'nldas'
# forcing_type = 'maurer'
# forcing_type = 'daymet'

# Load list of gauge IDs
with open('camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Grab parameter names
test_parameters = camels_utilities.load_sacsma_parameters(basins[0])
param_names = test_parameters.keys()

# Init storage
parameters = {}

# Loop through all basins
for basin in tqdm(basins):
    for param in range(10):

     # Load parameters
     basin_params = camels_utilities.load_all_sacsma_parameters(basin, forcing_type)
     
     # Store in dictionary
     parameters[basin] = basin_params

# Pickle
fname = f'results/all_camels_parameters_{forcing_type}.pkl'
with open(fname, 'wb') as f:
    pkl.dump(f, parameters)


