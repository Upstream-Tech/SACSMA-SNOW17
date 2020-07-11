import spotpy
import camels_utilities as camels
from optimizer import spotpy_setup

basin = '01054200'

# load data
# parameters = camels.load_sacsma_parameters(basin)
attributes = camels.load_basin_attributes(basin)
forcings, area = camels.load_forcings(basin)
observations = camels.load_usgs(basin, area)
# benchmarks = camels.load_discharge(basin)

# initialize spotpy optimizer (and model)
optimizer = spotpy_setup(
                 forcings=forcings,
                 observations=observations['QObs'],
                 latitude=attributes['gauge_lat'],
                 elevation=attributes['elev_mean'],
                 warmup=365*2)

# initialize spotpy sampler
sampler=spotpy.algorithms.sceua(optimizer,
                                dbname='SCE',
                                dbformat='csv')
max_model_runs = 1e5
# sampler.sample(max_model_runs, ngs=len(optimizer.optimized_parameter_names))
sampler.sample(max_model_runs, ngs=20)

# run the optimization
results = sampler.getdata()





