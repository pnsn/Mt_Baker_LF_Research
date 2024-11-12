"""
:module: Mt_Baker_LF_Research/src/python/drivers/create_mbs_tribe.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This driver clusters templates using waveform
    cross-correlation based grouping routines provided
    by EQcorrscan. These routines are re-packaged in the
    `eqcutil.core.clusteringtribe.ClusteringTribe` class.
"""

import os, glob
from pathlib import Path

from eqcorrscan import Template

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

# Get path for templates
tdir = os.path.join(root,'processed_data','templates','well_located_20km')
# Get list of template names
tnames = glob.glob(os.path.join(tdir,'*.tgz'))

# Major parameters
fill_value = 1
shift_len = 5.

ccckwargs = {'method': 'correlation_cluster',
             'replace_nan_distances_with': fill_value,
             'shift_len': shift_len,
             'corr_thresh': 0.7,
             'allow_individual_trace_shifts': False,
             'show': False,
             'cores': 'all',
             'save_corrmat': False}


# save directory
savename = os.path.join(root,'processed_data',
                        'clusters','well_located_20km',
                        f'fv_{fill_value}_sl_{shift_len}')

### PROCESSING SECTION ###
# Load Templates
Logger.info(f'Constructing ClusteringTribe from {len(tnames)} saved templates')
ctr = ClusteringTribe()
for tfile in tnames:
    Logger.debug(f'loading: {tfile}')
    itemp = Template().read(tfile)
    ctr.extend(itemp)
Logger.info(f'{len(ctr.templates)} of {len(tnames)} loaded')

# Run Clustering
Logger.info(f'Running waveform-based clustering analysis on ClusteringTribe')
ctr.cluster(**ccckwargs)

# Save to Disk
Logger.info(f'Writing results to: {savename}.tgz')
ctr.write(savename)

# Signal End
Logger.critical('SUCCESSFULLY COMPLETED PROCESS - EXITING')