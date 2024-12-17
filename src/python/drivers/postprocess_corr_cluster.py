import os, logging
from pathlib import Path

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
# Path to CTR Save File
CTRF = TMPD / 'corr_cluster.tgz'


### CONNECT / LOAD ###
# Connect To Banks
EBANK = EventBank(EBBP)
# Load Inventory
INV = read_inventory(str(INVF))
# Load ClusteringTribe
CTR = ClusteringTribe().read(str(CTRF))
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

# Join etype to clusters
CTR.clusters = CTR.clusters.join(df_O['etype'],how='left')
CTR.clusters = CTR.clusters.join(df_E, how='left')

df_C = CTR.clusters
