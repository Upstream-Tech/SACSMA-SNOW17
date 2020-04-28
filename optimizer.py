from model import model
import camels_utilities
import spotpy



class spotpy_setup(object):
    def __init__(self):

      # Gauge ID to test
      basin = '01013500'

      # Load all necessary data
      forcings = camels_utilities.load_forcings(basin)
      benchmark = camels_utilities.load_discharge(basin)
      attributes = camels_utilities.load_basin_attributes(basin)

      # Init model
      self.model = model(forcings=forcings,
                         latitude=attributes['gauge_lat'],
                         elevation=attributes['elev_mean'],
                         observations=benchmark,
                         warmup_timesteps=730)

      self.params = [spotpy.parameter.Uniform('uztwm', 1    , 1000 ),
                     spotpy.parameter.Uniform('uzfwm', 1    , 1000 ),
                     spotpy.parameter.Uniform('uzk'  , 0.1  , 0.7  ),
                     spotpy.parameter.Uniform('pctim', 0.005, 0.005),
                     spotpy.parameter.Uniform('adimp', 0    , 0    ),
                     spotpy.parameter.Uniform('riva' , 0    , 0    ),
                     spotpy.parameter.Uniform('zperc', 1    , 250  ),
                     spotpy.parameter.Uniform('rexp' , 1.1  , 6    ),
                     spotpy.parameter.Uniform('lztwm', 1    , 100  ),
                     spotpy.parameter.Uniform('lzfsm', 1    , 1000 ),
                     spotpy.parameter.Uniform('lzfpm', 1    , 1000 ),
                     spotpy.parameter.Uniform('lzsk' , 0.001, 0.25 ),
                     spotpy.parameter.Uniform('lzpk' , 0    , 0.03 ),
                     spotpy.parameter.Uniform('pfree', 0    , 1    ),
                     spotpy.parameter.Uniform('side' , 0    , 0    ),
                     spotpy.parameter.Uniform('uztwc', 0    , 10 ),
                     spotpy.parameter.Uniform('uzfwc', 0    , 10 ), 
                     spotpy.parameter.Uniform('lztwc', 0    , 10 ),
                     spotpy.parameter.Uniform('lzfsc', 0    , 10 ),
                     spotpy.parameter.Uniform('lzfpc', 0    , 10 ),
                     spotpy.parameter.Uniform('adimc', 0    , 10 )
                     ]

    def parameters(self):
        parameters = spotpy.parameter.generate(self.params)
       # print(parameters)
        return parameters

    def simulation(self, vector):
        simulations = self.model._run(uztwm = vector[0],
                                      uzfwm = vector[1],
                                      uzk   = vector[2],
                                      pctim = vector[3],
                                      adimp = vector[4],
                                      riva  = vector[5],
                                      zperc = vector[6],
                                      rexp  = vector[7],
                                      lztwm = vector[8],
                                      lzfsm = vector[9],
                                      lzfpm = vector[10],
                                      lzsk  = vector[11],
                                      lzpk  = vector[12],
                                      pfree = vector[13],
                                      side  = vector[14],
                                      uztwc = vector[15],
                                      uzfwc = vector[16],
                                      lztwc = vector[17],
                                      lzfsc = vector[18],
                                      lzfpc = vector[19],
                                      adimc = vector[20])

        return simulations

    def evaluation(self, evaldates=False):
        if evaldates:
            return self.model.eval_dates
        else:
            return self.model.observations

    def objectivefunction(self, simulation, evaluation):
        objectivefunction = spotpy.objectivefunctions.rmse(evaluation, simulation)
        return objectivefunction
