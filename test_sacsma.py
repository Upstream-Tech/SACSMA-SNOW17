import pandas as pd
import sacsma_utilities
import camels_utilities
import matplotlib.pyplot as plt

# Location of camels data
camels_dir = '/Users/grey/workspace/camels_data/'

# Which forcings
forcing_type = 'nldas'
#forcings = 'maurer'
#forcings = 'daymet'

# Gauge ID to test
gauge_id = '09378170'

# Load all necessary data
parameters = camels_utilities.load_sacsma_parameters(camels_dir, forcing_type, gauge_id)
forcings = camels_utilities.load_forcings(camels_dir, forcing_type, gauge_id)
attributes = camels_utilities.load_basin_attributes(camels_dir, gauge_id)
benchmark = camels_utilities.load_discharge(camels_dir, forcing_type, gauge_id)
dates = forcings['Date']

# Run SAC-SMA
outputs, states = sacsma_utilities.run_sacsma(dates, forcings, parameters, attributes)

# Plot results
plt.plot(outputs['surf'] + outputs['grnd'], label='Our SAC-SMA')
plt.plot(benchmark['OBS_RUN'], label='Observations')
plt.plot(benchmark['MOD_RUN'], label='NCAR SAC-SMA')
plt.grid()
plt.show()

import pdb
pdb.set_trace()
