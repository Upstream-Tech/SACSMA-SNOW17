import pandas as pd
import aldens_example.sacsma_utilities as sacsma
import camels.camels_utilities as camels
import matplotlib.pyplot as plt

# Load list of gauge IDs
with open('camels/camels_basin_id_list.txt', 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]

# Gauge ID to test
#basin = '02465493' 
basin = basins[200]

# Load all necessary data
forcings = camels.load_forcings(basin)
benchmark = camels.load_discharge(basin)

# Run SAC-SMA
fluxes, states = sacsma.run_sacsma(forcings, basin)

# Plot results
plt.plot(fluxes['sacsma_surf'] + fluxes['sacsma_grnd'], label='Our SAC-SMA')
plt.plot(benchmark['OBS_RUN'], label='Observations')
plt.plot(benchmark['MOD_RUN'], label='NCAR SAC-SMA')
plt.plot(forcings['PRCP(mm/day)'], label='Precipitation', linewidth=0.1)
plt.legend()
plt.grid()
plt.show()
