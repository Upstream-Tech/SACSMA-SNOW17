import pandas as pd
import numpy as np
from model.potential_evap import priestley_taylor_pet, calc_surface_pressure
import sacsma_source.sac.ex_sac1 as sacsma
import sacsma_source.sac.duamel as unit_hydrograph
import sacsma_source.snow19.exsnow as snow17


class model(object):

    def __init__(self,
                 forcings: pd.DataFrame,
                 latitude: float,
                 elevation: float,
                 default_parameters: pd.DataFrame,
                 warmup: int):

        self.default_parameters = default_parameters
        self.latitude = latitude.astype(float)
        self.elevation = elevation.astype(float)
        self.warmup = warmup

        # Timestep in different units
        self.dt_seconds = int((forcings.index[1] - forcings.index[0]).total_seconds())
        self.dt_days = self.dt_seconds / 86400
        self.dt_hours = self.dt_seconds / (60 * 60)

        # test period
        self.eval_dates = forcings.index.values[warmup:]

        # Keys for parameters, states, fluxes
        self.sacsma_parameter_keys = [
            'uztwm', 'uzfwm', 'uzk', 'pctim', 'adimp', 'riva', 'zperc', 'rexp', 'lztwm', 'lzfsm', 'lzfpm', 'lzsk',
            'lzpk',
            'pfree', 'side', 'rserv'
        ]
        self.snow17_parameter_keys = [
            'scf', 'mfmax', 'mfmin', 'uadj', 'si', 'nmf', 'tipm',
            'mbase', 'pxtemp', 'plwhc', 'daygm'
        ]
        self.snow17_adc_parameter_keys = [
            'adc1', 'adc2', 'adc3', 'adc4', 'adc5', 'adc6', 'adc7', 'adc8', 'adc9', 'adc10', 'adc11'
        ]
        self.hydrograph_parameter_keys = [
            'unit_shape', 'unit_scale'
        ]

        # Extract vectors as numpy arrays for speed
        self.precipitation = forcings['PRCP(mm/day)'].values
        self.temperature = 0.5 * (forcings['Tmax(C)'].values + forcings['Tmin(C)'].values)
        self.day = forcings['Day'].values
        self.month = forcings['Mnth'].values
        self.year = forcings['Year'].values

        # Calculate potential evaporation as numpy array
        self.pet = priestley_taylor_pet(forcings['Tmin(C)'], forcings['Tmax(C)'], forcings['SRAD(W/m2)'],
                                   self.latitude, self.elevation,
                                   forcings.index.to_series().dt.dayofyear)

        # estimate surface pressure [hPa] as a function of elevation only
        self.surf_pres = calc_surface_pressure(self.elevation)


    def run(self, args):
        return self._run(*args)

    def _run(self, test_parms, test_param_names):

        parameters = self.default_parameters.copy()
        for p, parm in enumerate(test_param_names):
            parameters[parm] = test_parms[p]

        # Extract parameters as numpy array
        sacsma_parameters_np = parameters[self.sacsma_parameter_keys].values.copy().astype(float)
        hydrograph_parameters_np = parameters[self.hydrograph_parameter_keys].values.copy().astype(float)
        snow17_parameters_np = parameters[self.snow17_parameter_keys].values.copy().astype(float)
        snow17_adc_parameters_np = parameters[self.snow17_adc_parameter_keys].values.copy().astype('f4')

        # Initial states must be numpy scalars for the fortran intent(inout) to work.
        # Also must have trailing decimals, or fortran will treat as integers even if typed as real in the pyf interface.
        uztwc = np.array(200.)
        uzfwc = np.array(200.)
        lztwc = np.array(200.)
        lzfsc = np.array(200.)
        lzfpc = np.array(200.)
        adimc = np.array(200.)

        tprev = np.array(0.)
        cs = np.full(19, 0., dtype='f4')

        # Init storage as numpy arrays
        qq = np.full(self.precipitation.shape[0], np.nan)

        # Time loop
        for t in range(self.precipitation.shape[0]):

            # Run the model at one timestep
            raim, sneqv, snow, snowh = snow17.exsnow19(self.dt_seconds, self.dt_hours, self.day[t], self.month[t], self.year[t],
                                                             self.precipitation[t], self.temperature[t],
                                                             self.latitude,
                                                             *snow17_parameters_np,
                                                             self.elevation,
                                                             self.surf_pres,
                                                             snow17_adc_parameters_np,
                                                             cs, tprev)

            surf, grnd, qq[t], tet = sacsma.exsac(self.dt_seconds, raim, self.temperature[t], self.pet[t],
                                                     *sacsma_parameters_np,  # 0.,
                                                     uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc)

        # channel routing
        m_unit_hydro = 1000
        n_unit_hydro = self.precipitation.shape[0] + m_unit_hydro
        hydrograph_qq = unit_hydrograph.duamel(qq, *hydrograph_parameters_np, self.dt_days, n_unit_hydro,
                                                 m_unit_hydro,
                                                 1, 0)

        return hydrograph_qq[self.warmup:-m_unit_hydro]






        # # Set constant snow17 parameter
        # adc = np.array([0.05, 0.15, 0.26, 0.45, 0.5, 0.56, 0.61, 0.65, 0.69, 0.82, 1.0]).astype('f4')

# # Parameter data type
    # mfmax = np.array(mfmax).astype(float)
    # mfmin = np.array(mfmin).astype(float)
    # uadj = np.array(uadj).astype(float)
    # si = np.array(si).astype(float)
    # scf = np.array(scf).astype(float)
    # pxtemp = np.array(pxtemp).astype(float)
      # uztwm = np.array(uztwm).astype(float)
      # uzfwm = np.array(uzfwm).astype(float)
      # lztwm = np.array(lztwm).astype(float)
      # lzfpm = np.array(lzfpm).astype(float)
      # lzfsm = np.array(lzfsm).astype(float)
      # uzk   = np.array(uzk).astype(float)
      # lzpk  = np.array(lzpk).astype(float)
      # lzsk  = np.array(lzsk).astype(float)
      # zperc = np.array(zperc).astype(float)
      # rexp  = np.array(rexp).astype(float)
      # pfree = np.array(pfree).astype(float)
      # ushape = np.array(ushape).astype(float)
      # uscale = np.array(uscale).astype(float)
