import os, sys
from pathlib import Path
# Import from installed dependencies
from obspy import read_inventory
from obspy.clients.fdsn import Client
from eqcorrscan import Template, Tribe

# Import visualization tools for Pyrocko's Snuffler GUI
from pyrocko import obspy_compat
obspy_compat.plant()

# Map local paths
# Set repository root direcotry absolute path
root = Path(__file__).parent.parent.parent.parent
# Set ebank absolute path
ebpath = os.path.join(root,'data','XML','QUAKE','BANK')
# Set import path for python scripts
pypath = os.path.join(root, 'src','python')
sys.path.append(pypath)
# Import local methods into scope
from eqc_utils.catalog_utils import (connect_to_eventbank,
                                     apply_phase_hints,
                                     filter_picks)
from eqc_utils.wavebank_utils import connect_to_wavebank


# Connect to eventbank
ebank = connect_to_eventbank(base_path=ebpath)
print('successfully connected to local machine EventBank')

# Get Station Inventory
sta_XML = os.path.join(root, 'data','XML','INV','station_inventory_50km_MBS_IRISWS_STATION.xml')
inv = read_inventory(sta_XML)

# Connect to IRIS Webservices
client = Client('IRIS')
print('successfully connected to IRIS webservices client')
# client = connect_to_wavebank()
# print('successfully connected to local machine WaveBank')
# Select a single event from our EventBank
event_df = ebank.read_index()
eid = event_df.event_id.iloc[10]
cat = ebank.get_events(event_id=eid)
print(cat)

# Copy phase labels from associations/arrivals to picks
cat = apply_phase_hints(cat)

# Subset to MBW, MBW2, and SHUK picks with P phase hints
# cat_p = filter_picks(cat.copy(),
#                     phase_hints=['P','p'])

print(f'Original Event Has: {len(cat[0].picks)} picks')
# print(f'Filtered Event Has: {len(cat_p[0].picks)} picks')

# Pull 5 minutes of data for all stations & channels with picks
# Downsample to 50 Hz (if not already there)
# Reject traces that don't have enough data (skip_short_chans)
# Trim a 50 second window of waveform data for all components
# Window starts 5 sec before pick time & runs for 50 sec
# Bandpass filter all data at 0.15 - 10 Hz
mytribe = Tribe().construct(
    method='from_client',
    client_id=client,
    catalog=cat,
    samp_rate=50.,
    lowcut=0.15,
    highcut=10,
    filt_order=4.,
    length=50.,
    prepick=5.,
    process_len=300.,
    all_horiz=True,
    delayed=True,
    skip_short_chans=True
)

# Visualize template in Snuffler
mytribe.templates[0].st.snuffle(catalog=cat, inventory=inv)
