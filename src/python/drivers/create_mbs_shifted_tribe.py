"""
:module: Mt_Baker_LF_Research/src/python/drivers/create_mbs_shifted_tribe.py
:auth: Benz Poobua 
:email: spoobu@uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: This driver creates shifted templates for a sub-catalog
    of PNSN events near Mt. Baker that satisfy pre-defined
    requirements for being "well-located". The shift is refered to 
    `align_mbs_traces.py`
"""

import os, logging
from pathlib import Path

import pandas as pd
from obspy import Catalog, read_events
from obspy.clients.fdsn import Client
from eqcorrscan import Template, Tribe

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.augment.template import rename_templates, augment_template
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(__name__,level=logging.DEBUG)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### CONTROL PARAMETER DEFINITION SECTION ###

# Define paths data sources
# Get root directory of the repository
root = Path(__file__).parent.parent.parent.parent
# Safety check that it's running from root directory
if os.path.split(root)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {root}')
else:
    Logger.warning(f'running {__file__} from {root}')

# Path to shifted catalog
subcat_path = os.path.join(root,'processed_data','sub_catalog')

# Load shifted catalog
xml_files = [f for f in os.listdir(subcat_path) if f.endswith('.xml')]
subcat_filepath = os.path.join(subcat_path, xml_files[0])

# Define waveform client
wfclient = Client("IRIS")

# Define Pick Filtering kwargs
pfkwargs = {'phase_hints': ['P'],
            'enforce_single_pick': 'earliest',
            'stations': ['MBW','MBW2','SHUK',
                         'SAXON','RPW2','PASS',
                         'MULN','RPW2']}

# Define template construction kwargs
ckwargs = {'method': 'from_client',
           'client_id': wfclient,
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

# Define template augmentation kwargs
akwargs = {'client': wfclient,
           'padding': 120.}

# Define Save Directory
savedir = os.path.join(root, 'processed_data','shifted_templates','well_located_20km')

##### PROCESSING SECTION #####
if not os.path.exists(savedir):
    Logger.info(f'making savedir: {savedir}')
    os.makedirs(savedir)
else:
    Logger.info(f'saving to existing savedir: {savedir}')

## Get event catalog
cat = read_events(subcat_filepath)

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
        # Logger.info(f'Augmenting waveform data on {template.name}')
        # template = augment_template(template=template, **akwargs)    
        Logger.info(f'Saving template {tribe[0].name} to {savedir}')
        template.write(os.path.join(savedir, f'{template.name}.tgz'))


    

