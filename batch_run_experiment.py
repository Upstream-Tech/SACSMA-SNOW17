import argparse
from pathlib import Path
import pandas as pd
from ruamel.yaml import YAML
import pickle as pkl
import os
import shutil
import multiprocessing
from joblib import Parallel, delayed
from optimize_single_basin import run_single_basin as run_basin

#def batch_run_experiment(config_file, max_model_runs, out_dir, use_cores_frac):

parser = argparse.ArgumentParser()
parser.add_argument('--config_file', type=str, required=True)
parser.add_argument('--max_model_runs', type=int, required=True)
parser.add_argument('--dds_trials', type=int, default=1)
parser.add_argument('--out_dir', type=str, default='results')
parser.add_argument('--use_cores_frac', type=float, default=1)
parser.add_argument('--algorithm', type=str, default='DDS')
#parser.add_argument('--forcing_type', type=str, required=True)
args = parser.parse_args()

# read config file
with Path(args.config_file).open('r') as fp:
    yaml = YAML(typ="safe")
    yaml.allow_duplicate_keys = True
    cfg = yaml.load(fp)  

# extract forcing type
assert(len(cfg['forcings']))==1
forcing_type = cfg['forcings'][0]
print('Using forcings of type: ', forcing_type)

# extract training dates
train_dates = pd.date_range(cfg['train_start_date'], cfg['train_end_date'])

# list all basins in this experiment    
with open(cfg['train_basin_file'], 'r') as f:
    basins = f.readlines()
basins = [basin.strip() for basin in basins]
assert len(basins) == 531

# create output directory
out_dir_run = Path(args.out_dir) / f"{str(args.config_file).split('/')[-1][:-4]}"
shutil.rmtree(out_dir_run, ignore_errors=True)
os.mkdir(out_dir_run)

# parallel loop over basins
num_cores = multiprocessing.cpu_count()
use_n_cores = int(num_cores*args.use_cores_frac)
print(f'Using {use_n_cores} cores of {num_cores} total.')
Parallel(n_jobs=use_n_cores)(delayed(run_basin)(basin, 
                                                forcing_type,
                                                train_dates,
                                                args.algorithm,
                                                args.max_model_runs,
                                                args.dds_trials,
                                                out_dir_run) 
                             for basin in basins)
#run_basin(basin=basins[0], 
#          forcing_type=forcing_type, 
#          train_dates=train_dates, 
#          algorithm=args.algorithm, 
#          max_model_runs=args.max_model_runs, 
#          dds_trials=args.dds_trials, 
#          out_dir_run=out_dir_run)



