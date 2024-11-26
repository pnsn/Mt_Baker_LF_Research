"""
:module: Mt_Baker_LF_Research/src/python/drivers/stack_mbs_traces.py
:auth: Benz Poobua
:email: spoobu@uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: This driver stacks traces using shifted clusters.

"""

import os, glob
from pathlib import Path
from obspy import Stream
from eqcorrscan.utils.stacking import linstack

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging
Logger = setup_terminal_logger(name=__name__)
Logger.addHandler(CriticalExitHandler(exit_code=1))
# Get absolute path for repo root directory
root = Path(__file__).parent.parent.parent.parent
# Safety check that it's running from root directory
if os.path.split(root)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {root}')
else:
    Logger.warning(f'running {__file__} from {root}')

# Get path for ClusteringTribe
processed_path = os.path.join(root,'processed_data','shifted_clusters','well_located_20km')
ctr_file_path = os.path.join(processed_path,'hyper_parameter_trials')
ctr_files = glob.glob(os.path.join(ctr_file_path, 'fv_mean_sl_1.0_ls_False.tgz'))

# Load the ClusteringTribe
filepath = ctr_files[0]
ctr = ClusteringTribe().read(filepath)
Logger.info('ClusteringTribe loaded')
ctr2 = ctr.copy()

# Construct input
st_list = [temp.st for temp in ctr2.templates]

# Stack
stack = linstack(st_list)