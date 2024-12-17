"""
:module: Mt_Baker_LF_Research/src/python/drivers/create_mbs_tribe.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This driver creates templates for a sub-catalog
    of PNSN events near Mt. Baker that satisfy pre-defined
    requirements for being "well-located" and augments the
    templates by expanding the channel coverage to 3-component
    data where available.
"""

import os, logging
from pathlib import Path

import pandas as pd
from obspy import Catalog
from eqcorrscan import Template, Tribe
from obsplus import EventBank, WaveBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.augment.template import rename_templates, augment_template
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.DEBUG)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path's for WaveBank and EventBank 
WBBP = ROOT / 'data' / 'WF' / 'BANK'
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Path to Event Subset
C20K = ROOT/'data'/'Events'/'MtBaker_20km_radius_phases.csv'
# Path to INV
INVF = ROOT / 'data' / 'XML' / 'INV' / 'station_inventory_50km_MBS_IRISWS_STATION.xml'



# Define Pick Filtering kwargs
pfkwargs = {'phase_hints': ['P'],
            'enforce_single_pick': 'earliest',
            'stations': ['MBW','MBW2','SHUK']}
            # ,
            #              'SAXON','RPW2','PASS',
            #              'MULN','RPW2']}

# Define template construction kwargs
ckwargs = {'method': 'from_client',
           'lowcut': 0.5,
           'highcut': 20.,
           'filt_order': 4,
           'samp_rate': 50.,
           'prepick': 10.,
           'length': 50.,
           'process_len': 300.,
           'min_snr': 3.,
           'parallel': True,
           'num_cores': 4, 
           'save_progress': False
           }


#### PROCESSING SECTION ###

# Define event_id set
Logger.info('Loading 20km event subset file')
df_20k = pd.read_csv(C20K)
unique_20km_evids = df_20k.evid.unique().astype(str)
Logger.debug(f'got {len(unique_20km_evids)} unique EVIDs')

# Connect to EventBank
ebank = EventBank(base_path=str(EBBP))
# Assert that the index exits
if len(ebank.read_index()) == 0:
    Logger.critical('EventBank shows zero entries')
else:
    Logger.info('Connected to EventBank')

# Connect to WaveBank
Logger.info(f'Initializing/connecting to WaveBank at {WBBP}')
wbank = WaveBank(base_path=str(WBBP))
# Assert that the index exits
if len(wbank.read_index()) == 0:
    Logger.critical('WaveBank shows zero entries')
else:
    Logger.info("Connected to WaveBank")


# Add WaveBank to template construction kwargs
ckwargs.update({'client_id': wbank})




# Define template augmentation kwargs
akwargs = {'client': wfclient,
           'padding': 120.}


# Define Save Directory
savedir = os.path.join(root, 'processed_data','templates','full_catalog')



##### PROCESSING SECTION #####
if not os.path.exists(savedir):
    Logger.info(f'making savedir: {savedir}')
    os.makedirs(savedir)
else:
    Logger.info(f'saving to existing savedir: {savedir}')


## CONNECT TO EVENT BANK
Logger.info(f'Connecting to EventBank with base_path: {ebpath}')
ebank = EventBank(**ebkwargs)
# Get eb index
df_eb = ebank.read_index()
# Update index as EVIDs
idx = [row.event_id.split('/') for _, row in df_eb.iterrows()]
idx = [f'{_e[-2].lower()}{_e[-1]}' for _e in idx]
df_eb.index = idx

# Safety check that event bank connection is good
if len(df_eb) > 0:
    Logger.info(f'Successfully connected: {len(df_eb)} events in bank')
else:
    Logger.critical('Returned an empty event bank - exiting')

## GET 20km EVIDS
_df = pd.read_csv(phase_csv)
evid20 = [f'uw{eid}' for eid in _df.evid.value_counts().index.values]

## GET WELL-LOCATED EVIDS
with open(well_located_evids, 'r') as file:
    lines = file.readlines()

evidwl = list(set([f'uw{int(eid)}' for eid in lines]))

## Event Bank Pre-Filtering
df_eb = df_eb[(df_eb.index.isin(evid20))]# & (df_eb.index.isin(evidwl))]
if len(df_eb)> 0:
    Logger.info(f'Subsampled EventBank contents to {len(df_eb)} events')
else:
    Logger.critical('No events selected from EventBank index subsampling')

## Get event catalog
cat = ebank.get_events(event_id=df_eb.event_id)
if len(cat) > 0:
    Logger.info(f'Fetched {len(cat)} event\'s metadata from the EventBank')
else:
    Logger.critical(f'Failed to fetch any event metadata from the EventBank')

## Apply phase hints
Logger.info('Applying phase hints')
cat = apply_phase_hints(cat)

## Filter Catalog
Logger.info(f'Filtering catalog picks for: {pfkwargs}')
pfkwargs.update({'catalog': cat})
cat = filter_picks(**pfkwargs)

sta_set = dict(zip(pfkwargs['stations'], [0]*len(pfkwargs['stations'])))
for _e in cat.events:
    for _p in _e.picks:
        sta_set[_p.waveform_id.station_code] += 1
Logger.info(f'Pick counts are: {sta_set}')
for event in cat.events:
    icat = Catalog(events=[event])
    ckwargs.update({'catalog': icat})
    Logger.info(f'constructing template for: {icat[0].resource_id}')
    try:
        tribe = Tribe().construct(**ckwargs)
    except:
        Logger.warning('failed to construct template - continuing to next')
        continue
    # Rename templates
    Logger.info(f'renaming templates with EVID')
    tribe = rename_templates(tribe)

    for template in tribe:
        # Augment Templates
        Logger.info(f'Augmenting waveform data on {template.name}')
        template = augment_template(template=template, **akwargs)    
        Logger.info(f'Saving template {tribe[0].name} to {savedir}')
        template.write(os.path.join(savedir, f'{template.name}.tgz'))


        


# ## Iterate across events & create templates
# for event in cat.events:
#     Logger.info('processing event')
# for evid, row in df_eb.iterrows():
#     event = ebank.get_events(event_id=row.event_id)

