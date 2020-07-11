import glob
from pathlib import Path
import pandas as pd

# Hard-code paths
DATA_DIR = '/home/gsnearing/projects/camels_data'
FORCING_TYPE = 'nldas'


def load_all_sacsma_parameters(forcing_type: str):

  # Construct file name from pieces
  filenames = glob.glob(f'{DATA_DIR}/model_output/{forcing_type}/**/*_model_parameters.txt')

  # loop through all files
  for i, filename in enumerate(filenames):

    # Load a dictionary of parameter values
    parameters = {}
    with open(filename) as f:
      for line in f:
        key, val = line.split()
        parameters[key] = float(val)

    # Convert to pandas series
    if not('parameters_df' in locals()):
      parameters_df = pd.DataFrame(index = parameters.keys())
    parameters_df[str(i)] = pd.Series(parameters)#.to_frame

  return parameters_df


def load_sacsma_parameters(gauge_id: str):

  # Construct file name from pieces
  filename = glob.glob(f'{DATA_DIR}/model_output/{FORCING_TYPE}/**/{gauge_id}_*_model_parameters.txt')[0]

  # Load a dictionary of parameter values
  parameters = {}
  with open(filename) as f:
    for line in f:
      key, val = line.split()
      parameters[key] = val

  # Convert to dataframe
  parameters_series = pd.Series(parameters)#.to_frame

  return parameters_series


def load_basin_attributes(gauge_id: str):

  # Add the attributes path extention
  attributes_dir = Path(DATA_DIR) / 'camels_attributes_v2.0'

  # Check for existence
  if not attributes_dir.exists():
    raise RuntimeError(f"Attribute folder not found at {attributes_dir}")

  # Grab all attributes files
  txt_files = attributes_dir.glob('camels_*.txt')

  # Read-in attributes into one big dataframe
  dfs = []
  for txt_file in txt_files:
    df_temp = pd.read_csv(txt_file, sep=';', header=0, dtype={'gauge_id': str})
    df_temp = df_temp.set_index('gauge_id')
    dfs.append(df_temp)

  df = pd.concat(dfs, axis=1)

  # convert huc column to double digit strings
  df['huc'] = df['huc_02'].apply(lambda x: str(x).zfill(2))
  df = df.drop('huc_02', axis=1)

  # Choose only the basin we want
  attributes = df.loc[gauge_id]

  return attributes


def load_forcings(gauge_id: str):

  # Grab the correct forcing file
  forcing_files = glob.glob(
      f'{DATA_DIR}/basin_dataset_public_v1p2/basin_mean_forcing/{FORCING_TYPE}/**/{gauge_id}_*_forcing_leap.txt')
  assert len(forcing_files) == 1
  forcing_file = forcing_files[0]

  # load area from header
  with open(forcing_file, 'r') as fp:
    content = fp.readlines()
    area = int(content[2])

  # Grab the data
  forcings = pd.read_csv(forcing_file, sep='\s+', header=3)

  # Datetime index
  forcings['Date'] = pd.to_datetime(
      forcings['Year'] * 10000 + forcings['Mnth'] * 100 + forcings['Day'], format='%Y%m%d')
  forcings = forcings.set_index('Date')

  return forcings, area


def load_discharge(gauge_id: str):

  # Grab the correct forcing file
  filename = glob.glob(f'{DATA_DIR}/model_output/{FORCING_TYPE}/**/{gauge_id}_*_model_output.txt')[0]

  # Grab the data
  output = pd.read_csv(filename, sep='\s+')

  # Datetime index
  output['Date'] = pd.to_datetime(output['YR'] * 10000 + output['MNTH'] * 100 + output['DY'], format='%Y%m%d')
  output = output.set_index('Date')

  return output


def load_usgs(gauge_id: str, area: int):

  # Grab the correct forcing file
  filename = glob.glob(f'{DATA_DIR}/basin_dataset_public_v1p2/usgs_streamflow/**/{gauge_id}_streamflow_qc.txt')[0]

  # Grab the data
  col_names = ['basin', 'Year', 'Mnth', 'Day', 'QObs', 'flag']
  obs = pd.read_csv(filename, sep='\s+', header=None, names=col_names)

  # unit conversion cfs --> mm/day
  obs.QObs = 28316846.592 * obs.QObs * 86400 / (area * 10 ** 6)

  # Datetime index
  obs['Date'] = pd.to_datetime(obs.Year.map(str) + "/" + obs.Mnth.map(str) + "/" + obs.Day.map(str))
  obs = obs.set_index('Date')

  return obs

