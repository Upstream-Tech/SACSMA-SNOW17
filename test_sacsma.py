import pandas as pd
import sacsma_utilities
import camels_utilities
import matplotlib.pyplot as plt

# Gauge ID to test
gauge_id = '09378170'

# Load all necessary data
forcings = camels_utilities.load_forcings(gauge_id)
benchmark = camels_utilities.load_discharge(gauge_id)

# Run SAC-SMA
fluxes, states = sacsma_utilities.run_sacsma(forcings, gauge_id)

# Plot results
plt.plot(fluxes['sacsma_surf'] + fluxes['sacsma_grnd'], label='Our SAC-SMA')
plt.plot(benchmark['OBS_RUN'], label='Observations')
plt.plot(benchmark['MOD_RUN'], label='NCAR SAC-SMA')
plt.legend()
plt.grid()
plt.show()
