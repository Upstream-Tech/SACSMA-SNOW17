from pathlib import Path
#from batch_run_experiment import batch_run_experiment

# fraction of total cores to use
use_cores_frac = 1

# run directory
run_dir = Path('/home/gsnearing/projects/lstm_based_hydrology/extreme_year_runs/')
out_dir = Path('./results/')

# load config files
config_files = list(run_dir.glob('**/config.yml'))
print(f'There are {len(config_files)} experiments.')

# optimizer hypers
max_model_runs = 1e1 # 1e3 # 1e5
dds_trials = 1e2

# loop over experiments
for f, config_file in enumerate(config_files):
  print(config_file, max_model_runs, out_dir, use_cores_frac)
#  batch_run_experiment(config_file, max_model_runs, out_dir, use_cores_frac)    


