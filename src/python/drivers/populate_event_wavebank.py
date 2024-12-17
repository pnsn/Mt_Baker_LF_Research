import os, logging, glob
from pathlib import Path

import numpy as np
import pandas as pd

from obsplus import WaveBank, EventBank
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.client import FDSNNoDataException
from obspy import UTCDateTime, read_inventory, Stream

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
Logger.addHandler(CriticalExitHandler(1))

# Define absolute path to the repository root directory
ROOT = Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical('ROOT does not look like the root directory path of this repository')
C20K = ROOT/'data'/'Events'/'MtBaker_20km_radius_phases.csv'
EBBP = ROOT/'data'/'XML'/'QUAKE'/'BANK'
WBBP = ROOT/'data'/'WF'/'BANK'
WBPS = '{year}/{julday}'
WBNS = '{seedid}_{year}{julday}_{hour}{minute}{second}'
INVP = ROOT/'data'/'XML'/'RESP'
PLOG = WBBP.parent/'processing_log.csv'

# Define temporal padding
wlen = 300. # [sec] lead padding relative to event origin time
wpos = range(-2, 4)
channels = '[BHE]?[ZNE12]' # Channel codes to use for waveform queries

# Define event_id set
df_20k = pd.read_csv(C20K)
unique_20km_evids = df_20k.evid.unique().astype(str)

### PROCESSING SECTION ###

# Create output directory
if not os.path.exists(WBBP):
    os.makedirs(WBBP)

# Initialize waveform client
client = Client("IRIS")

# Load station inventories
inv_list = glob.glob(str(INVP/'*.xml'))
if inv_list == []:
    Logger.critical('Did not find any station inventories')
else:
    for _e, _f in enumerate(inv_list):
        Logger.debug(f'Reading inventory: {os.path.split(_f)[-1]}')
        if _e == 0:
            INV = read_inventory(_f, format='stationXML')
        else:
            INV += read_inventory(_f)
    INV = INV.select(channel=channels)
    _c = INV.get_contents()
    msg = f'Loaded {len(inv_list)} stationXML files with: '
    for _k, _v in _c.items():
        msg += f' {len(_v)} {_k}'
    Logger.info(msg)    

# Initialize EventBank Connection
ebank = EventBank(base_path=str(EBBP))
# Assert that the index exits
if len(ebank.read_index()) == 0:
    Logger.critical('EventBank shows zero entries')

# Initialize WaveBank Instance
Logger.info(f'Initializing/connecting to WaveBank at {WBBP}')
wbank = WaveBank(base_path=str(WBBP),
                 path_structure=str(WBPS),
                 name_structure=str(WBNS))


# Get summary data for all events
df_eb = ebank.read_index()
# Subset Events to 20km events
mask = [row.event_id.split('/')[-1] in unique_20km_evids for _, row in df_eb.iterrows()]
df_eb = df_eb[mask]
# Sort by time for sanity's sake
df_eb = df_eb.sort_values('time')

# Get a set of all net.sta strings
ns_set = set()
for net in INV.networks:
    for sta in net.stations:
        ns_set.add(f'{net.code}.{sta.code}')

# Create logging file if it doesn't exist
log_hdr = 'event_id,window_no,bulk_no,st_no,st_sta_no,build_mode,commit_status,t1,t2'
if not os.path.isfile(PLOG):
    with open(PLOG, 'w') as LOG:
        LOG.write(f'{log_hdr}\n')
    df_LOG = pd.DataFrame(data=[],columns=log_hdr.split(','))
# Otherwise, load log
else:
    df_LOG = pd.read_csv(PLOG)

_idx = 1
# Iterate across events
for _, row in df_eb.iterrows():
    # Re-open log to append new entries
    with open(PLOG, '+a') as LOG:
        Logger.info(f'Processing: {row.event_id} ({_idx}/{len(df_eb)})')
        _idx += 1
        # Safety catch for already-processed data
        if row.event_id in df_LOG.event_id:
            Logger.info(f'already processed - continuing to next')
            continue
        # Get origin time from index
        tO = UTCDateTime(row.time)
        # Get active stations
        current_inv = INV.select(time=tO)
        # QC if no current stations
        if len(current_inv) == 0:
            Logger.info('No active stations in inventory - continuing to next event')
            LOG.write(f'{row.event_id},,,,,EmptyInv,False,{tO},{tO}\n')
            continue
        # Get start & end time for bulk request
        st = Stream()
        # Iterate across position vector (makes chunks)
        for _x in wpos:
            t1 = tO + _x*wlen
            t2 = t1 + wlen
            # Compose bulk request
            bulk = [tuple(_e.split('.') + [t1, t2]) for _e in current_inv.get_contents()['channels']]
            if len(bulk) == 0:
                breakpoint()
            Logger.debug(f'requesting {len(bulk)} channel records')
            # First try bulk request
            try:
                ist = client.get_waveforms_bulk(bulk)
                st += ist
                _bm = 'BulkRequest'
            except FDSNNoDataException:
                Logger.warning('Bulk request raised FDSNNoDataException - attempting single request approach')
                for _e, _b in enumerate(bulk):
                    kwargs = dict(zip(['network','station','location','channel','starttime','endtime'], _b))
                    try:
                        st += client.get_waveforms(**kwargs)
                    except FDSNNoDataException:
                        continue
                _bm = 'IndividualRequest'
                Logger.info(f'Single request approach returned {len(st)} waveform segments for {len(bulk)} initial requests')

            st_sta_no = len(set([tr.id for tr in st]))
            if any(np.ma.is_masked(tr.data) for tr in st):
                st = st.split()

            if len(st) > 0:
                Logger.info('adding waveforms to wavebank')
                try:
                    wbank.put_waveforms(st)
                    commit_status = True
                except:
                    commit_status = False
            else:
                commit_status=False

            LOG.write(f'{row.event_id},{_x},{len(bulk)},{len(st)},{st_sta_no},{_bm},{commit_status},{t1},{t2}\n')

    
