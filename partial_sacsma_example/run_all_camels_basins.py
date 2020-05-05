import pandas as pd
import sacsma_utilities
import camels_utilities
import matplotlib.pyplot as plt
import pickle
from tqdm import tqdm

# Location of camels data
camels_dir = '/Users/grey/workspace/camels_data/'

# Which forcings
forcing_type = 'nldas'
#forcings = 'maurer'
#forcings = 'daymet'

# Load list of gauge IDs
with open('camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Loop through all basins
for gauge_id in tqdm(basins):

    # Load all necessary data
    parameters = camels_utilities.load_sacsma_parameters(camels_dir, forcing_type, gauge_id)
    forcings = camels_utilities.load_forcings(camels_dir, forcing_type, gauge_id)
    attributes = camels_utilities.load_basin_attributes(camels_dir, gauge_id)
    benchmark = camels_utilities.load_discharge(camels_dir, forcing_type, gauge_id)
    dates = forcings['Date']

    # Run SAC-SMA
    outputs, states = sacsma_utilities.run_sacsma(dates, forcings, parameters, attributes)

    # Save results
    fname = f'results/{gauge_id}.pkl'
    with open(fname, 'wb') as f:
        pickle.dump([outputs,states], f)

