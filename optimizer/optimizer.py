import spotpy
from optimizer.model import model
import pandas as pd
import numpy as np
import model.camels_utilities as camels


class spotpy_setup(object):
    def __init__(self,
                 forcings: pd.DataFrame,
                 observations: pd.Series,
                 latitude: float,
                 elevation: float,
                 warmup: int = 365 * 2):

        self.observations = observations.values[warmup:]

        # Parameter bounds
        forcing_types = ['nldas', 'daymet', 'maurer']
        columns = pd.MultiIndex.from_product([['min', 'max'], forcing_types])
        parameter_df = pd.DataFrame(columns=columns)
        for forc in forcing_types:
            forc_df = camels.load_all_sacsma_parameters(forc)
            parameter_df[('min', forc)] = forc_df.min(axis=1)
            parameter_df[('max', forc)] = forc_df.max(axis=1)
        parameter_df['all_mean'] = parameter_df.mean(axis=1)
        parameter_df['all_min'] = parameter_df.min(axis=1)
        parameter_df['all_max'] = parameter_df.max(axis=1)
        parameter_df = parameter_df.drop('PT_COEF', axis=0)

        # Parameter sampling
        self.all_parameter_names = list(parameter_df.index)
        self.optimized_parameter_names = []
        self.params = []
        for parm in list(parameter_df.index):
            if np.abs(parameter_df.loc[parm, 'all_max'].values[0] - parameter_df.loc[parm, 'all_min'].values[0]) > 1e-5:
                self.params.append(spotpy.parameter.Uniform(parm, parameter_df.loc[parm, 'all_min'],
                                                            parameter_df.loc[parm, 'all_max']))
                self.optimized_parameter_names.append(parm)

        # Init model
        self.model = model(forcings=forcings,
                           latitude=latitude,
                           elevation=elevation,
                           default_parameters=parameter_df['all_mean'],
                           warmup=warmup)

    def parameters(self):
        parameters = spotpy.parameter.generate(self.params)
        return parameters

    def simulation(self, vector):
        simulations = self.model._run(vector, self.optimized_parameter_names)
        return simulations

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.model.eval_dates
        else:
            return self.observations

    def objectivefunction(self, simulation, evaluation):
        objectivefunction = spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction
