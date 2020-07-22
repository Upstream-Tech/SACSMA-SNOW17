import pandas as pd
import numpy as np
from itertools import chain
import spotpy
import pickle as pkl
import model.camels_utilities as camels
from optimizer.optimizer import spotpy_setup

# function to run a single basin
def run_single_basin(basin, forcing_type, train_dates, algorithm, max_model_runs, dds_trials, out_dir_run):

    # training dates for this basin
    sd = train_dates['start_dates'][basin]
    ed = train_dates['end_dates'][basin]
    obj_fun_dates = pd.DataFrame(list(chain.from_iterable(pd.date_range(sdi, edi) for sdi, edi in zip(sd, ed))), columns = ('train_dates',))

    # load data
    mask_dates = obj_fun_dates['train_dates']
    attributes = camels.load_basin_attributes(basin)
    forcings, area = camels.load_forcings(basin, forcing_type)
    observations = camels.load_usgs(basin, area)

    # set up optimizer
    if algorithm == 'SCE':
      algorithm_minimize =  True
    elif algorithm == 'DDS':
      algorithm_minimize = False
    optimizer = spotpy_setup(forcings=forcings,
                             observations=observations['QObs'],
                             latitude=attributes['gauge_lat'],
                             elevation=attributes['elev_mean'],
                             algorithm_minimize=algorithm_minimize,
                             mask_dates=mask_dates)

    # SCE hyperparameters
    if algorithm == 'SCE':
      sampler=spotpy.algorithms.sceua(optimizer,
                                      dbname='SCE',
                                      dbformat='ram',
                                      parallel='seq',
                                      save_sim=False)
      sampler.sample(repetitions=int(max_model_runs), ngs=len(optimizer.optimized_parameter_names))
    # DDS hyperparameters
    elif algorithm == 'DDS':
      sampler=spotpy.algorithms.dds(optimizer, 
                                    dbname='DDS', 
                                    dbformat='ram',
                                    parallel='seq',
                                    save_sim=False)
      sampler.sample(repetitions=int(max_model_runs), trials=int(dds_trials))

    # get best parameters
    results = sampler.getdata()
    best_parameters = spotpy.analyser.get_best_parameterset(results,maximize=(not algorithm_minimize))
    best_parameters_df = pd.DataFrame(best_parameters)
    for key in best_parameters_df.keys():
        new_key = key.split('par')[-1]
        best_parameters_df = best_parameters_df.rename(columns={key: new_key})
    best_parameters_series = best_parameters_df.transpose()[0]

    # get simulation with best parameters
    parm_vector = best_parameters_series.loc[optimizer.optimized_parameter_names].values
    sim = optimizer.simulation(parm_vector)

    # grab best likelihood
    try:
        likelihoods=results['like']
    except ValueError:
        likelihoods=results['like1']
    if algorithm_minimize:
      best_likelihood=np.nanmin(likelihoods)
    else:
      best_likelihood=np.nanmax(likelihoods)
    index = np.where(likelihoods==best_likelihood)[0]

    # save output
    outfile = out_dir_run / f"{basin}.pkl"
    with open(outfile, 'wb') as f:
        pkl.dump([best_parameters_series, sim, likelihoods, best_likelihood], f)



