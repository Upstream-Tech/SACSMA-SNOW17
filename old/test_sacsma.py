import pandas as pd
import sacsma_utilities as sacsma
import camels.camels_utilities as camels
import matplotlib.pyplot as plt

# Load list of gauge IDs
with open('camels/camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Gauge ID to test
# basin = '06352000' 
basin = basins[300]
print('Basin', basin)

# Load all necessary data
forcings, area = camels.load_forcings(basin)
benchmark = camels.load_discharge(basin)

# Run SAC-SMA
fluxes, states = sacsma.run_sacsma(forcings, basin)

# test = pd.read_csv('/Users/grey/workspace/SACSMA-SNOW17/sacsma/SACSMA-Snow17/bin/test_model_output_1.txtb', sep='\s+')
# test['Date'] = pd.to_datetime(test['year'] * 10000 + test['month'] * 100 + test['day'], format='%Y%m%d')
# test = test.set_index('Date')

# Plot results
plt.plot(fluxes['sacsma_uh_qq'], lw=3, label='Our SAC-SMA')
# plt.plot(test['mod_flow_routed'], lw=1, label='Fortran SAC-SMA')
# plt.plot(test['obs_flow'], lw=3, label='Observation')
plt.plot(benchmark['OBS_RUN'], label='CAMELS Observations')
plt.plot(benchmark['MOD_RUN'], label='CAMELS SAC-SMA')
plt.plot(forcings['PRCP(mm/day)'], label='Precipitation', linewidth=0.1)
plt.legend()
plt.grid()
plt.show()
