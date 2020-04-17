from datetime import datetime, timezone
import sys

import keras as k
from keras import Sequential
from keras.layers import LSTM, TimeDistributed, Dense
import matplotlib.pyplot as plt
import numpy as np
from scipy import optimize
from sklearn.metrics import mean_squared_error

from leaf_river import load_leaf_river, DAYS_PER_TIMESTEP
import sacsma

date_str = datetime.now(timezone.utc).isoformat()[:-10].replace(':', '-')
MODEL_SAVE_PATH_FORMAT = f'checkpoints/lstm_{date_str}_{{epoch:02d}}-{{loss:.2f}}.hdf5'
TENSORBOARD_LOG_DIR = f'tensorboard/lstm_{date_str}'

SAMPLE_LENGTH = 365 # days
BATCH_SIZE = 10
EPOCH_SIZE = 20
EPOCHS = 400

def get_model(input_shape):
  model = Sequential()
  model.add(LSTM(30, input_shape=input_shape, return_sequences=True))
  model.add(TimeDistributed(Dense(1)))
  model.compile(optimizer='adam', loss='mean_squared_error')
  model.summary()
  return model

def train_lstm():
  X_train, X_test, y_train, y_test = load_leaf_river()

  train_generator = generate_data(X_train, y_train, BATCH_SIZE, SAMPLE_LENGTH)
  test_generator = generate_data(X_test, y_test, BATCH_SIZE, SAMPLE_LENGTH)

  input_shape = (SAMPLE_LENGTH, X_train.shape[-1])

  model = get_model(input_shape)

  save_model_callback = k.callbacks.ModelCheckpoint(MODEL_SAVE_PATH_FORMAT,
      monitor='loss', verbose=0, save_best_only=False, save_weights_only=False, period=5)

  tensorboard_callback = k.callbacks.TensorBoard(log_dir=TENSORBOARD_LOG_DIR)

  model.fit_generator(train_generator, EPOCH_SIZE, EPOCHS,
      validation_data=test_generator, validation_steps=2,
      callbacks=[save_model_callback, tensorboard_callback])
  print(f'Finished training')

def predict(weights_path):
  X_train, X_test, y_train, y_test = load_leaf_river()

  # Calculate consistent value ranges for each plot
  viz_max_percentile = 99.9
  viz_max_precip = np.percentile(X_train[:, 0], viz_max_percentile)
  viz_max_pet = np.percentile(X_train[:, 1], viz_max_percentile)
  viz_max_flow = np.percentile(y_train, viz_max_percentile)
  plt.rcParams["figure.figsize"] = [10.0, 4.0]
  #plt.rcParams["figure.figsize"] = [10.0, 8.0]

  test_generator = generate_data(X_test, y_test, BATCH_SIZE, SAMPLE_LENGTH)

  # LSTM
  input_shape = (SAMPLE_LENGTH, X_train.shape[-1])
  model = get_model(input_shape)
  model.load_weights(weights_path)

  # SAC-SMA
  params = [12.521811802883539, 5.211564916610896, 0.7708752328297034, -0.1050904048188684, 0.06065184395299683, -0.0023806638772820793, 37.38798398791645, 37.223184062889004, 241.31540250716318, 926.1296101258185, 138.32163203343808, 0.24007759493580283, 0.015664888700619946, 0.3285593242050193, 0.014977676821340569, -0.0004886515232322797]
  state = [np.array(0.0) for _ in range(6)]

  for X_batch, y_true_batch in test_generator:
    y_pred_lstm_batch = model.predict(X_batch)

    for x, y_true, y_pred_lstm in zip(X_batch, y_true_batch, y_pred_lstm_batch):
      y_pred_sac = sacsma_predict(x, params, state) # state is reused across calls. that's good since we are iterating cronologically

      mse_sac = mean_squared_error(y_true, y_pred_sac)
      mse_lstm = mean_squared_error(y_true, y_pred_lstm)

      plt.subplot(311)
      plt.title('Precipitation')
      plt.plot(x[: , 0])
      plt.ylim([0, viz_max_precip])
      plt.xticks([])

      plt.subplot(312)
      plt.title('Potential ET (Energy)')
      plt.plot(x[: , 1])
      plt.ylim([0, viz_max_pet])
      plt.xticks([])

      plt.subplot(313)
      plt.title('Streamflow')
      plt.plot(y_true, label='Observed')
      plt.plot(y_pred_lstm, label=f'Predicted LSTM (mse: {mse_lstm:.2f})')
      plt.plot(y_pred_sac, label=f'Predicted SAC-SAM (mse: {mse_sac:.2f})')
      plt.ylim([0, viz_max_flow])
      plt.legend()

      plt.show()

def generate_data(X, Y, batch_size, sample_len):
  samples_in_batch = batch_size * sample_len
  i = 0
  while True:
    start = i % X.shape[0]
    end = (i+samples_in_batch) % X.shape[0]
    if start < end:
      x = X[start:end]
      y = Y[start:end]

      x_batch = x.reshape(batch_size, sample_len, -1)
      y_batch = y.reshape(batch_size, sample_len, -1)

      yield x_batch, y_batch

    i += samples_in_batch

def calibrate_sacsma():
  X_train, X_test, y_train, y_test = load_leaf_river()

  # SAC-SMA Parameters
  # Parameter bounds via Samuel, J. et al. "Assessing model state and forecasts variation in hydrologic data assimilation" (2014)
  # https://drive.google.com/open?id=1-Buzq6twSLYfW6ooxmGfgWMjU_Vq0S4p, Table 2
  # Note the paper gives no range for parameters RIVA, SIDE and SAVED, not sure why. TODO try giving some range
  PARAMETER_BOUNDS = [
    (1., 100.),     # uztwm  Upper-zone tension water maximum storage (mm)
    (1., 100.),     # uzfwm  Upper-zone free water maximum storage (mm)
    (0.01, 0.99),   # uzk    Upper-zone free water lateral depletion rate (per day)
    (0., 0.1),      # pctim  Impervious fraction of the watershed area (fraction)
    (0.01, 0.1),    # adimp  Additional impervious area (fraction)
    (0., 0.),       # riva   Riparian vegetation area (fraction I think)
    (1., 50.),      # zperc  Maximum percolation rate (no units?)
    (1., 50.),      # rexp   Exponent of the percolation equation
    (50., 500.),    # lztwm  Lower-zone tension water maximum storage (mm)
    (100., 1000.),  # lzfsm  Lower-zone free water supplemental maximum storage (mm)
    (100., 2000.),  # lzfpm  Lower-zone free water primary maximum storage (mm)
    (0.1, 0.5),     # lzsk   Lower-zone supplemental free water depletion rate (per day)
    (0.0001, 0.02), # lzpk   Lower-zone primary free water depletion rate (per day)
    (0., 0.5),      # pfree  Fraction percolating from upper- to lower-zone free water storage (fraction)
    (0., 0.),       # side   Ratio of deep recharge to channel base flow (fraction)
    (0., 0.)        # saved  Fraction of lower zone free water not transferable to tension water (fraction)
  ]
  # Start optimization with each parameter set to the middle of it's range
  initial_params = [np.mean(bounds) for bounds in PARAMETER_BOUNDS]

  # SAC-SMA State variables
  # uztwc Upper-zone tension water storage content (mm)
  # uzfwc Upper-zone free water storage content (mm)
  # lztwc Lower-zone tension water storage content (mm)
  # lzfpc Lower-zone free primary water storage content (mm)
  # lzfsc Lower-zone free secondary water storage content (mm)
  # adimc Additional impervious area content directly link to the stream network (mm)

  # Dumb initialization of all state variables to zero (aka everything's dry)
  state = [np.array(0.0) for _ in range(6)]

  # Calibrate the model parameters
  result = optimize.minimize(sacsma_mean_squared_error,
      initial_params,
      args=(state, X_train, y_train),
      method='Nelder-Mead',
      bounds=PARAMETER_BOUNDS)
      #callback=print)
  print('Calibration', 'success' if result.success else 'failure, optimization failed.')
  print('Parameters:')
  print(list(params))

def predict_sacsma():
  X_train, X_test, y_train, y_test = load_leaf_river()

  params = [12.521811802883539, 5.211564916610896, 0.7708752328297034, -0.1050904048188684, 0.06065184395299683, -0.0023806638772820793, 37.38798398791645, 37.223184062889004, 241.31540250716318, 926.1296101258185, 138.32163203343808, 0.24007759493580283, 0.015664888700619946, 0.3285593242050193, 0.014977676821340569, -0.0004886515232322797]

  state = [np.array(0.0) for _ in range(6)]

  # Calculate consistent value ranges for each plot
  viz_max_percentile = 99.99
  viz_max_precip = np.percentile(X_train[:, 0], viz_max_percentile)
  viz_max_pet = np.percentile(X_train[:, 1], viz_max_percentile)
  viz_max_flow = np.percentile(y_train, viz_max_percentile)
  plt.rcParams["figure.figsize"] = [8.0, 6.0]

  test_generator = generate_data(X_test, y_test, BATCH_SIZE, SAMPLE_LENGTH)

  for X_batch, y_true_batch in test_generator:
    for x, y_true in zip(X_batch, y_true_batch):
      y_pred = sacsma_predict(x, params, state) # state is reused across calls. that's good since we are iterating cronologically

      mse = mean_squared_error(y_true, y_pred)

      plt.subplot(311)
      plt.title('Precipitation')
      plt.plot(x[: , 0])
      plt.ylim([0, viz_max_precip])
      plt.xticks([])

      plt.subplot(312)
      plt.title('Potential ET (Energy)')
      plt.plot(x[: , 1])
      plt.ylim([0, viz_max_pet])
      plt.xticks([])

      plt.subplot(313)
      plt.title('Streamflow')
      plt.plot(y_true, label='Observed')
      plt.plot(y_pred, label=f'Predicted sac-sam (MSE: {mse:.2f})')
      plt.ylim([0, viz_max_flow])
      plt.legend()

      plt.show()

def sacsma_mean_squared_error(params, state, X, y_true):
  y_pred = sacsma_predict(X, params, state)
  assert y_true.shape == y_pred.shape
  return mean_squared_error(y_true, y_pred)

def sacsma_predict(X, params, state):
  """
  X.shape = (timesteps, 2) where [:, 0] is precipitation, [:, 1] is potential ET.
  returns y, 1D list of streamflow values of len timesteps
  """

  predicted_streamflow = np.full(X.shape[0], np.nan)

  # iterate over daily timesteps
  for i, x in enumerate(X):
    """
    sacsma.fland1 fortran function
    Performs calculations for the Sacramento rainfall-runoff model.
    Parameter list: All parameters are passed, here we brake them into how they are actually used.
    Inputs:
      pxv, rainfall per one time interval
      edmnd, potential evapotranspiration
      dt, simulation time interval, in days
    Parameters (calibrated)
      uztwm, uzfwm, uzk, pctim, adimp, riva, zperc, rexp, lztwm, lzfsm, lzfpm, lzsk, lzpk, pfree, side, saved, SAC-SMA parameters
    States
      uztwc, uzfwc, lztwc, lzfsc, lzfpc, adimc, SAC-SMA states
    Outputs
      surf, simulated fast runoff
      grnd, simulated slow runoff
      tet, simulated actual evapotranspiration
    More parameter explanations https://drive.google.com/open?id=1-Buzq6twSLYfW6ooxmGfgWMjU_Vq0S4p
    More conceptual SAC-SMA documentation http://www.nws.noaa.gov/oh/hrl/nwsrfs/users_manual/part2/_pdf/23sacsma.pdf
    """
    # Inputs
    precip = x[0]   # precipitation
    potental_et = x[1] # potential evapotranspiration

    ### Call Fortran Model. Advances one day. Updates state arguments by reference. ###
    outputs = sacsma.fland1(precip, potental_et, DAYS_PER_TIMESTEP, *state, *params)

    # Outputs (sacsma variable name)
    runoff_fast = outputs[0] # surface, simulated fast runoff (surf)
    runoff_slow = outputs[1] # ground, simulated slow runoff (grnd)
    et = outputs[2] # simulated actual evapotranspiration (tet)

    predicted_streamflow[i] = runoff_fast + runoff_slow

  return np.array(predicted_streamflow)

if __name__ == '__main__':
  assert len(sys.argv) > 1, 'Usage: $ {} train|predict [args]'.format(sys.argv[0])
  command = sys.argv[1]
  if command == 'train-lstm':
    train_lstm(*sys.argv[2:])
  elif command == 'calibrate-sacsma':
    calibrate_sacsma(*sys.argv[2:])
  elif command == 'predict':
    predict(*sys.argv[2:])
  elif command == 'predict-sacsma':
    predict_sacsma(*sys.argv[2:])
  else:
    print('Unrecognized command.')
    exit(1)
