"""
TODO: Nate fill in docstring header
"""

import os, sys
from pathlib import Path

import pandas as pd
from obspy import read_inventory
from obspy.clients.fdsn import Client
from eqcorrscan import Tribe
from pyrocko import obspy_compat

obspy_compat.plant()


# Define paths
root = Path(__file__).parent.parent.parent.parent
ebpath = os.path.join(root,'data','XML','QUAKE','BANK')
wlpath = os.path.join(root,'results','tables','well_located_mt_baker_evids_20km.txt')
evidpath = os.path.join(root,'data','Events','MtBaker_20km_radius_phases.csv')
invpath = os.path.join(root, 'data','XML','INV','station_inventory_50km_MBS_IRISWS_STATION.xml')

pypath = os.path.join(root,'src','python')
savedir = os.path.join(root,'processed_data','templates','20km_catalog')
sys.path.append(pypath)
# Import utility Scripts
from eqc_utils.catalog_utils import connect_to_eventbank, apply_phase_hints, filter_picks
from eqc_utils.template_utils import rename_templates, augment_template
from eqc_utils.logger_utils import setup_standard_logger

# Create savedir
if not os.path.exists(savedir):
    os.makedirs(savedir)

# Define Catalog filtering Parametes
pickfiltkw = {
    'stations': ['MBW','MBW2','SHUK','SAXON','MULN','PUBD'],
    'phase_hints': ['P'],
    'enforce_single_pick': 'earliest'
}

# Define Template Construction Parameters
constructkw = {
    'method': 'from_client',
    'lowcut': 0.5,
    'highcut': 20.,
    'samp_rate': 50.,
    'filt_order': 4,
    'length': 40.,
    'prepick': 5.,
    'process_len': 3600,
    'delayed': True,
    'all_horiz': False,
    'min_snr': 3.,
    'parallel': True,
    'num_cores': 6
}
# Define additional options for construct
# TODO: ensure that removing this section does not mess things up. - superceded by `augment_template`
construct_options = {
    'all_vertical': True,
    'all_horiz': False,
    'vertical_chans': ['Z','N','E','1','2','3'],
    'horizontal_chans': []
}

constructkw.update(construct_options)

### DRIVER SECTION PAST HERE ###
# Set Up Logging
Logger = setup_standard_logger(__name__)


## Connect to WaveBank
ebank= connect_to_eventbank(base_path=ebpath)
Logger.info(f'Connected to EventBank at {ebpath}')

## Connect to IRIS webservices client
client = Client('IRIS')
Logger.info(f'Initialized IRIS webservices client. Verifying connection...')
try:
    client.get_stations(network='UW', station='GPW', level='station')
except Exception as e:
    Logger.critical(f'{type(e).__name__}: {e}')
    sys.exit(1)
Logger.info('Success! Adding to Tribe().construct kwargs')

# Add client_id to construct keyword arguments 
constructkw.update({'client_id': client})
Logger.info('Added client to construct key-word arguments')

# Read EVIDs
## Load AQMS data
df = pd.read_csv(evidpath)

evids = df.evid.value_counts().index.values
evids = [f'uw{evid}' for evid in evids]
Logger.info(f'Selected {len(evids)} unique EVIDs for the 20km radius catalog')

# with open(wlpath, 'r') as _f:
#     evids = _f.read().splitlines()
# # Format to match ebank index
# evids = [f'uw{evid}' for evid in evids]
# Logger.info(f'Loaded {len(evids)} (e.g., {evids[0]}) from well-located subset')

# Subsample Events
df_evid = ebank.read_index()
df_subset = df_evid[df_evid.index.isin(evids)]
Logger.info(f'Filtered event bank contents from {len(df_evid)} to {len(df_subset)} well-located EVIDS ')

counter = 1
total_count = len(df_subset)
for evid, row in df_subset.iterrows():
    # Announce iteration
    Logger.warning(f'...running {evid} ({counter}/{total_count})')
    # Check if already exists
    savename = os.path.join(savedir, f'{evid}')
    if os.path.isfile(f'{savename}.tgz'):
        Logger.warning(f'{savename} exists - skipping to next')
        counter += 1
        continue
    else:
        pass
    # Fetch Catalog for single event
    cat = ebank.get_events(event_id = df_evid.loc[evid].event_id)
    if len(cat) == 1:
        Logger.info('...got event from EventBank')
    else:
        Logger.warning('...no event in EventBank, skipping to next')
        continue
    # Apply phase hints
    cat = apply_phase_hints(cat)
    Logger.info(f'...copying arrival phases to pick phase_hints')
    # Filter Catalog
    cat = filter_picks(cat, **pickfiltkw)
    if len(cat) == 0:
        Logger.warning('...filtering eliminated all picks. Skipping to next')
        counter += 1
        continue
    else:
        Logger.info(f'...filtered catalog: {len(cat.events[0].picks)} picks')
    
    # Update constructkw
    constructkw.update({'catalog': cat})
    # Create template

    # Try to generate Tribe
    try:
        itribe = Tribe().construct(**constructkw)
        rename_templates(itribe, prefix='uw')
    except Exception as e:
        Logger.warning(f'Could not constuct tribe: {type(e).__name__}: {e}')
        continue

    for template in itribe.templates:
        Logger.info(f'Attempting to augment {template.name}')
        template = augment_template(template, client)
        Logger.info(f'Saving Template {template.name} as {savename}')
        template.write(savename)
    counter += 1

Logger.critical('SUCCESSFULLY COMPLETED - ENDING')

    # if len(itribe) >= 1:
    #         Logger.info(f'Tribe of {len(itribe.templates)} generated!')
    #         rename_templates(itribe, prefix='uw')
    #         Logger.info(f'Saving to {savename}')
    #         itribe.write(filename=savename)
    #         breakpoint()
    
    # 


# # Get Events from Event Bank
# cat = ebank.get_events(event_id=df_wlevid.event_id.values)
# Logger.info(f'Loaded event catalog from EventBank: {len(cat)} events.')

# # Apply phase hints to cat
# cat = apply_phase_hints(cat)
# Logger.info(f'Applied phase hints to catalog')

# # Check that phase hints are applied (eventually shift to testing/method)
# matches = 0; mismatches = 0
# for _e in cat:
#     for _p in _e.picks:
#         if hasattr(_p, 'phase_hint'):
#             if isinstance(_p.phase_hint,str):
#                 if len(_p.phase_hint) == 1:
#                     matches += 1
#                 else:
#                     mismatches += 1
#             else:
#                 mismatches += 1
#         else:
#             mismatches += 1
# if mismatches > 0:
#     Logger.warning(f'Found {mismatches} picks without phase_hints of {mismatches + matches} total picks scanned')
# else:
#     Logger.info(f'All {matches} picks in catalog have single-character phase_hint values')


# # Filter for stations and select phases
# cat = filter_picks(cat, **pickfiltkw)

# Logger.info(f'Filtered catalog to {len(cat)} events with {sum([len(_e.picks) for _e in cat])} total picks - phase hints: {pickfiltkw["phase_hints"]}')

# # Try to create a tribe of 1 with an event object
# test_tribe = Tribe().construct(catalog=Catalog(events=[cat[-1]]),**constructkw)

# test_tribe.templates[0].st.snuffle(catalog=Catalog(events=[cat[-1]]), inventory=inv)