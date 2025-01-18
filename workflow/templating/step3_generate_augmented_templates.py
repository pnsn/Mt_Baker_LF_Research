import os, logging, warnings
from pathlib import Path
import matplotlib.pyplot as plt

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
CTSD = PD_DIR / 'templates' / 'single_station' / 'xcc_test'

############ TEMPLATE CONSTRUCTION PARAMETERS ############

# Define location code aliases
LOC_ALIASES = {'01': ''}
MIN_PICK_SNR = 1.2
MIN_RMS = 3
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
    'min_snr': 2.,
    'parallel': True,
    'num_cores': 12,
    'save_progress': False
}

ccckwargs = {'method': 'xcc',
            'replace_nan_distances_with': 'mean',
            'shift_len': 2,
            'corr_thresh': 0.5,
            'allow_individual_trace_shifts': True,
            'cores': 'all'}
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

def get_template_pick_snr(tmp, scale=0.9):
    tr = tmp.st[0].copy()

    noise_RMS = _rms(tr.copy().trim(endtime=tr.stats.starttime + tmp.prepick*scale))
    onset_RMS = _rms(tr.copy().trim(starttime=tr.stats.starttime + tmp.prepick,
                                    endtime=tr.stats.starttime + (1. + scale)*tmp.prepick))
    # pick_snr = onset_RMS / noise_RMS
    return onset_RMS, noise_RMS

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
if Logger.level < 10:
    tqdm_disable = True
else:
    tqdm_disable = True

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
            tmp.name = uwevid
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
            # Double check SNR around pick
            onset_rms, noise_rms = get_template_pick_snr(tmp)
            pick_snr = onset_rms/noise_rms
            if pick_snr <= MIN_PICK_SNR:
                Logger.warning(f'Rejecting template {tmp.name} due to low SNR around pick ({pick_snr:.3f})')
                continue
            
            if noise_rms <= MIN_RMS and onset_rms <= MIN_RMS:
                Logger.warning(f'Rejecting template {tmp.name} due to low RMS on both sides of pick ({noise_rms:.3f} | {onset_rms:.3f})')
                continue

            # Append template to clustering tribe
            ictr += tmp

    # Add metadata to ictr.clusters
    holder = []
    for tmp in ictr:
        event = tmp.event
        prefor = event.preferred_origin()
        prefmag = event.preferred_magnitude()
        # Get essential origin information
        line = [tmp.name.split('_')[-1], event.comments[-2].text,
                tmp.st[0].stats.network, tmp.st[0].stats.station, tmp.st[0].stats.location,
                tmp.st[0].stats.channel, event.picks[0].evaluation_mode, prefor.time,
                prefor.longitude, prefor.latitude, prefor.depth]
        # Get origin uncertainties if provided
        if prefor.origin_uncertainty is None:
            line.append(-9999)
        else:
            line.append(prefor.origin_uncertainty.horizontal_uncertainty)
        try:
            sdep = prefor.depth_errors.uncertainty
            if not np.isfinite(sdep):
                sdep = -9999
            line.append(sdep)
        except:
            line.append(-9999)
        # Get magnitude and magnitude type
        line += [prefmag.mag, prefmag.magnitude_type]
        # Calculate pick SNR
        onset_rms, noise_rms = get_template_pick_snr(tmp)
        line += [onset_rms/noise_rms, onset_rms, noise_rms]
        holder.append(line)
    try:
        df = pd.DataFrame(holder, columns=['evid', 'etype','network','station','location',
                                        'channel','eval_mode','time',
                                        'longitude','latitude','depth',
                                        'horizontal_uncertainty','vertical_uncertainty',
                                        'magnitude','magtype','snr', 'noise_rms','onset_rms'])
    except:
        breakpoint()
    df.index.name='id_no'
    # Join etype and pick evaluation mode
    ictr.clusters = ictr.clusters.join(df, on='id_no')
    Logger.info(f'Constructed {len(ictr)} templates for {sta} (of {npotential} possible)')

    ## APPLY ALIASES TO TEMPLATES IN TRIBES WITH MORE THAN ONE CHANNEL CODE
    # Get unique set of trace IDs
    id_set = set([])
    for tmp in ictr:
        id_set.add(tmp.st[0].id)
    # If there is more than one NSLC
    if len(id_set) > 1:
        # Blind the band and instrument characters
        new_id_set = set([])
        for tmp in ictr:
            tmp.st[0].stats.channel='??Z'
            # Apply aliases (if applicable)
            loc = tmp.event.picks[0].waveform_id.location_code
            # If the location is in LOC aliases, also adjust that
            if loc in LOC_ALIASES.keys():
                Logger.warning(f'Applying location alias to {sta} {event_id} template')
                tmp.st[0].stats.location = LOC_ALIASES[loc]
            new_id_set.add(tmp.st[0].id)
    
    # # RUN CLUSTERING
    # Logger.info(f'Running clustering for {sta}')
    # ictr.cluster(**ccckwargs)
    # # Save station-specific clustering tribes
    isavename = str(CTSD/sta)
    ictr.write(isavename, compress=True)


    

