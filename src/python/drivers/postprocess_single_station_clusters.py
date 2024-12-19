import os, logging, glob
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from obspy import Catalog, read_inventory
from eqcorrscan import Template, Tribe
from obsplus import EventBank, WaveBank

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path for EventBank 
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Path to event metadata
AQMS = ROOT / 'data' / 'Events' / 'MtBaker_50km_radius_origins.csv'
# Path to templating directory
TMPD = ROOT / 'processed_data' / 'bank_templates' / 'd20km'
# Path to INV
INVF = ROOT / 'data' / 'XML' / 'INV' / 'station_inventory_50km_MBS_IRISWS_STATION.xml'

### CONNECT / LOAD ###
# Connect To EventBank
EBANK = EventBank(EBBP)
# Load Inventory
INV = read_inventory(str(INVF))

### GET ETYPES FROM AQMS METADATA ###
# Load Event Metadata
df_O = pd.read_csv(str(AQMS))
# Subset to evid and etype
df_O = df_O[['evid','etype']]
# Drop duplicates
df_O = df_O.drop_duplicates(keep='first')
# Generate matching index to CTR.clusters
df_O.index = [f'uw{_e}' for _e in df_O.evid]

# Get eventbank index
df_E = EBANK.read_index()
# Generate matching index to CTR.clusters
df_E.index = [os.path.splitext(os.path.split(row.path)[-1])[0] for _, row in df_E.iterrows()]

### LOAD SINGLE SITE(ish) CLUSTERS ###
flist = glob.glob(str(TMPD/'corr_cluster_*_FreeStep.tgz'))
tribes = {}
for _f in flist:
    if 'XSTA' in _f:
        continue
    elif 'Composite' in _f:
        continue
    path, fname = os.path.split(_f)
    parts = fname.split('_')
    _k = parts[2]
    Logger.info(f'loading {_f}')
    _ctr = ClusteringTribe().read(_f)
    _resave = False
    if 'etype' not in _ctr.clusters.columns:
        _resave = True
        Logger.info(f'appending etypes to clusters dataframe')
        _ctr.clusters = _ctr.clusters.join(df_O['etype'], how='left')
    if 'event_id' not in _ctr.clusters.columns:
        _resave = True
        Logger.info(f'appending EventBank summary to clusters dataframe')
        _ctr.clusters = _ctr.clusters.join(df_E, how='left')
    tribes.update({_k: _ctr})
    if _ctr.clusters.correlation_cluster.value_counts().index[0] != 1:
        _resave = True
        Logger.info('reindexing correlation_clusters')
        _ctr.reindex_columns()
    if _resave:
        Logger.info(f'resaving clustering tribe with appended event metadata')
        _ctr.write(_f)

# Concatenate groups with different thresholding values
df_C = pd.DataFrame()
# Iterate across thresholding values (10x corr_thresh)

# Iterate across clusteringtribes

for _t in [9, 8, 7, 6, 5, 4, 3, 2, 1]:
    _cct = _t*0.1
    Logger.info(f'rethresholding at NCC level: {_cct}')

    for _e, (_k, _ctr) in enumerate(tribes.items()):
        Logger.debug(f'running rethresholding on {_k}')
        # Run re-thresholding
        _series = _ctr.cct_regroup(corr_thresh=_cct, inplace=False)
        # Relabel
        _series.name = f'{_k}_{_t:02d}'
        # Concatenate
        df_C = pd.concat([df_C, _series], axis=1, ignore_index=False)

if 'etype' not in df_C.columns:
    df_C = df_C.join(df_O['etype'], how='left')
if 'event_id' not in df_C.columns:
    df_C = df_C.join(df_E, how='left')

df_C.to_csv(str(TMPD / 'multi_thresholded_clusters.csv'), header=True, index=True)
