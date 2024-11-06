
import os, sys, logging
from pathlib import Path

import pandas as pd
from eqcorrscan import Template, Tribe

# Load EQcorrscan_utils library
from eqcorrscan_utils.logging import CriticalExitHandler, setup_terminal_logger
from eqcorrscan_utils.eventbank import EventBank

# Set up logging at debug level reporting
Logger = setup_terminal_logger(__name__,level=logging.DEBUG)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Define paths data sources
# Get root directory of the repository
root = Path(__file__).parent.parent.parent.parent
# Safety check that it's running from root directory
if os.path.split(root)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {root}')
else:
    Logger.warning(f'running {__file__} from {root}')

# Define connection to EventBank in repository
ebpath = os.path.join(root,'data','XML','QUAKE','BANK')
ebkwargs = {'base_path': ebpath,
            'path_structure': '{year}',
            'name_structure': 'uw{event_id_end}'}
# Define location of 20km radius catalog with phases
phase_csv = os.path.join(root,'data','Events','MtBaker_20km_radius_phases.csv')


##### PROCESSING SECTION #####

## CONNECT TO EVENT BANK
Logger.info(f'Connecting to EventBank with base_path: {ebpath}')
ebank = EventBank(**ebkwargs)
df_eb = ebank.read_index()
# Safety check that event bank connection is good
if len(df_eb) > 0:
    Logger.info(f'Successfully connected: {len(df_eb)} events in bank')
else:
    Logger.critical('Returned an empty event bank - exiting')

## GET 20km EVIDS
_df = pd.read_csv(phase_csv)
evid = [f'uw{eid}' for eid in _df.evid.value_counts().index.values]

## MATCH 20km EVIDs to event bank contents
df_eb = df_eb[df_eb.index.isin(evid)]
if len(df_eb)> 0:
    Logger.info(f'Subsampled EventBank contents to {len(df_eb)} events')
else:
    breakpoint()

