import os, logging, warnings
from pathlib import Path

import numpy as np
import pandas as pd

from tqdm import tqdm
from obspy.core.event import Comment
from obspy.clients.fdsn import Client
from obsplus import EventBank, WaveBank
from eqcorrscan import Tribe

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.augment.catalog import filter_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Setup Logging
Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler())

# Ignore FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)


# MAP ABSOLUTE PATHS
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow"
# Augmented Event Bank Base Path
AEBBP = PD_DIR / "catalog" / "AUGMENTED_BANK"
# WaveBank Base Path
WBBP = ROOT / "data" / "WF" / "BANK"
# Preferred Event-Station File
PSEF = PD_DIR / 'catalog' / 'preferred_event_sta_picks.csv'
# Template Save Directory
CTSD = PD_DIR / 'templates' / 'single_station'

############ TEMPLATE CONSTRUCTION PARAMETERS ############

# Define location code aliases
ALIASES = {'01': ''}

# Define Template Construction Parameters
tckwargs = {
    'method':'from_client',
    'lowcut': 0.5,
    'highcut': 20.,
    'filt_order': 4,
    'samp_rate': 50.,
    'prepick': 5.,
    'length': 45.,
    'process_len': 300.,
    'min_snr': 3.,
    'parallel': True,
    'num_cores': 12,
    'save_progress': False
}

ccckwargs = {'method': 'correlation_cluster',
            'replace_nan_distances_with': 'mean',
            'shift_len': 5,
            'corr_thresh': 0.5,
            'allow_individual_trace_shifts': True,
            'show': False,
            'cores': 'all',
            'save_corrmat': False}
################################

### HELPER FUNCTIONS ###
def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))

def check_prepick_rms(trace, prepick):
    ntr = trace.copy().trim(endtime=trace.stats.starttime + prepick)
    amp = _rms(ntr.data)
    return amp


#### PROCESSING SECTION ####
if not os.path.exists(CTSD):
    os.makedirs(CTSD)

# Connect to Clients
# IRIS Webservices
IRIS = Client('IRIS')
# Local WaveBank
WBANK = WaveBank(WBBP)
# Pre-Filtered / Augmented EventBank
EBANK = EventBank(AEBBP)

# Read index
df_eb = EBANK.read_index()

# Get preferred station codes
STAS = set()
with open(PSEF, 'r') as _f:
    lines = _f.readlines()
for line in lines[1:]:
    parts = line.split(',')
    STAS.add(parts[4])
# Cleanup
del lines


# Iterate across event_ids
if Logger.level < 20:
    tqdm_disable = True
else:
    tqdm_disable = False

for sta in STAS:
    Logger.info(f'Processing events for station {sta}')
    npotential = 0
    ictr = ClusteringTribe()
    for event_id in tqdm(df_eb.event_id, disable=tqdm_disable):
        uwevid = ''.join(event_id.split('/')[-2:]).lower()
        Logger.debug(f'Constructing templates for event {event_id}')
        # Get catalog
        cat = EBANK.get_events(event_id=[event_id])
        # Subset for to station
        cat = filter_picks(cat, stations=[sta])
        # if there are no picks, continue to next
        if len(cat) == 0:
            continue
        else:
            npotential += 1
        # Try creating the templates from the local wavebank as default
        tckwargs.update({'client_id': WBANK, 'catalog': cat})
        try:
            itribe = Tribe().construct(**tckwargs)
            Logger.debug(f'Successful local build of {sta} {event_id} template')
        except:
            Logger.debug(f'Unsuccessful local build of {sta} {event_id} template')
            tckwargs.update({'client_id': IRIS})
            try:
                itribe = Tribe().construct(**tckwargs)
                Logger.debug(f'Successful webservice build of {sta} {event_id} template')
            except:
                Logger.debug(f'Unsuccessful webservice build of {sta} {event_id} template - Skipping')
                continue
        # If tribe constructed, but is empty somehow
        if len(itribe) == 0:
            Logger.debug(f'Construction of {event_id} unsuccessful')
            continue
        # If tribe construction made multiple templates
        elif len(itribe) > 1:
            breakpoint()
        else:
            tmp = itribe[0]
            # Rename template
            tmp.name = f'{sta}_{uwevid}'
            tmp.event.comments[-1].text = f'eqcorrscan_template_{tmp.name}'

            Logger.debug(f'Construction of template {tmp.name} successful')
            # Ensure single pick
            if len(tmp.event.picks) != 1:
                Logger.warning(f'template has {len(tmp.event.picks)} picks')
                breakpoint()
            # Ensure single trace
            tmp.st.merge(method=1, interpolation_samples=-1)
            if len(tmp.st) != 1:
                Logger.warning(f'template has {len(tmp.st)} traces despite merge')
                breakpoint()
            
            # Apply aliases (if applicable)
            loc = tmp.event.picks[0].waveform_id.location_code
            if loc in ALIASES.keys():
                Logger.warning(f'Applying location alias to {sta} {event_id} template')
                tmp.event.picks[0].waveform_id.location_code = ALIASES[loc]
                tmp.event.picks[0].comments.append(Comment(text=f'Location code aliased: "{loc}" to "{ALIASES[loc]}"'))
                tmp.st[0].stats.location = ALIASES[loc]
            
            # Append template to clustering tribe
            ictr += tmp
    # Add metadata to ictr.clusters
    holder = []
    for tmp in ictr:
        event = tmp.event
        prefor = event.preferred_origin()
        prefmag = event.preferred_magnitude()
        line = [tmp.name.split('_')[-1], event.comments[-2].text,prefor.time,
                tmp.st[0].stats.channel, event.picks[0].evaluation_mode,
                prefor.longitude, prefor.latitude, prefor.depth]
        if prefor.origin_uncertainty is None:
            line += [-9999, -9999]
        else:
            line += [prefor.origin_uncertainty.horizontal_uncertainty,
                    prefor.depth_errors.uncertainty]
        line += [prefmag.mag, prefmag.magnitude_type]
        holder.append(line)
    df = pd.DataFrame(holder, columns=['evid', 'etype',
                                       'channel','eval_mode','time'
                                       'longitude','latitude','depth',
                                       'horizontal_uncertainty','vertical_uncertainty',
                                       'magnitude','magtype'])
    df.index.name='id_no'
    # Join etype and pick evaluation mode
    ictr.clusters = ictr.clusters.join(df, on='id_no')
    Logger.info(f'Constructed {len(ictr)} templates for {sta} (of {npotential} possible)')

    # RUN CLUSTERING
    Logger.info(f'Running clustering for {sta}')
    ictr.cluster(**ccckwargs)
    # Save station-specific clustering tribes
    isavename = str(CTSD/sta)
    ictr.write(isavename, compress=True)

            # # Check if the pick location code has a known alias
            # if pick.waveform_id.location_code in ALIASES.keys():
            #     # Make doubly sure that the pick metadata matches
            #     if tmp.event.picks[0] != pick:
            #         breakpoint()
            #     # Update the location code in the template
            #     # itr

    

