import pickle as pkl 
import numpy as np 
import camels_utilities
from tqdm import tqdm
import pandas as pd 

# Location of camels data
camels_dir = '/Users/grey/workspace/camels_data/'

# Load list of gauge IDs
with open('camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Grab parameter names
test_attributes = camels_utilities.load_basin_attributes(basins[0])
att_names = test_attributes.keys()

# Init storage
attributes = {}

# Loop through all basins
for basin in tqdm(basins):

    # Load parameters
    basin_atts = camels_utilities.load_basin_attributes(basin)
    
    # Store in dictionary
    attributes[basin] = basin_atts

# Pickle
fname = f'results/all_camels_attributes.pkl'
with open(fname, 'wb') as f:
    pkl.dump(attributes,f)


