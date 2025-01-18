import os, logging, warnings, glob
from pathlib import Path

from eqcutil.util.logging import basic_logger_config, CriticalExitHandler
from eqcutil.core.clusteringtribe import ClusteringTribe

basic_logger_config(level=logging.WARNING)
Logger = logging.getLogger(os.path.split(__file__)[-1])
# ch = logging.StreamHandler()
Logger.addHandler(CriticalExitHandler())

# Ignore FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)


# MAP ABSOLUTE PATHS
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow"
# Template Save Directory
CTSD = PD_DIR / 'templates' / 'single_station' / 'xcc_test'

### CROSS CORRELATION CLUSERING KEY WORD ARGUMENTS ####

ccckwargs = {'method': 'xcc',
            'replace_nan_distances_with': 'mean',
            'shift_len': 10,
            'corr_thresh': 0.5,
            'allow_individual_trace_shifts': True,
            'cores': 'all'}

### PROCESSING SECTION ###

flist = glob.glob(str(CTSD / '*.tgz'))

for _f in flist:
    Logger.warning(f'reading {_f}')
    ctr = ClusteringTribe().read(_f)
    Logger.warning(f'Clustering')
    try:
        ctr.cluster(**ccckwargs)
    except:
        breakpoint()
    Logger.warning(f'Saving to same file')
    ctr.write(_f, compress=True)
