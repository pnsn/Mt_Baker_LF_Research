"""
TODO: Nate fill in docstring header
"""


import os, sys, glob
from pathlib import Path

from eqcorrscan import Tribe, Template
from eqcorrscan.utils.clustering import cluster

# Local mappings
root = Path(__file__).parent.parent.parent.parent
pypath = os.path.join(root,'src','python')
# Add repo to path
sys.path.append(pypath)
# Load repo utilities
from eqc_utils.logger_utils import setup_standard_logger
from eqc_utils.template_utils import compose_template_list, save_template_clustering_output


# Processing Parameters
ckwargs = {'show': False,
           'corr_thresh': 0.4,
           'shift_len': 10.,
           'allow_individual_trace_shifts': False,
           'save_corrmat': True,
           'replace_nan_distances_with': 'mean',
           'cores': 'all',
           'method': 'single',
           'metric': 'euclidian',
           'optimal_ordering': False}

# Path to data
loaddir = os.path.join(root,'processed_data','templates','20km_catalog')
# Get list of files to load
flist = glob.glob(os.path.join(loaddir,'*.tgz'))
# Save Path
save_dir = Path().cwd() / 'results' / 'tables' / '20km_clustering' / 'nanmean_fill'
# save_dir = Path().cwd() / 'results' / 'tables' / 'wl_clustering' / 'nan_mean_fill'

### DRIVER PROCESSING PAST HERE ###

# Initialize Logging
Logger = setup_standard_logger(__name__)

# Create save directory if it doesn't exist
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
    Logger.warning(f'generated new directory: {save_dir}')
else:
    Logger.warning(f'will write to: {save_dir}')

# Create input parameter save file
with open(os.path.join(save_dir,'params.csv'), 'w') as p:
    p.write('name,value\n')
    for _k, _v in ckwargs.items():
        p.write(f'{_k},{_v}\n')

Logger.info('writing input parameters to disk')

# Load templates from disk
mytribe = Tribe()
for f in flist:
    Logger.info(f'loading template: {f}')
    mytemplate = Template().read(f)
    mytemplate.st.merge()
    mytribe += mytemplate

Logger.info('creating template list')
tlist = compose_template_list(mytribe)

Logger.info('running clustering')
groups = cluster(tlist, **ckwargs)
Logger.info('clustering finished')
save_template_clustering_output(save_dir, groups, savename='full_test')

