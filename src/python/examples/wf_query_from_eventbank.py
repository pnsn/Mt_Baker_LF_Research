import os, sys, logging
import pandas as pd
from pathlib import Path
from obspy.clients.fdsn import Client
from obspy.clients.fdsn.header import FDSNNoDataException

# Set repository root directory absolute path
root = Path(__file__).parent.parent.parent.parent

pypath = os.path.join(root, 'src', 'python')
sys.path.append(pypath)

import eqc_utils.catalog_utils as cutil 

# Import visualization tools for Pyrocko's Snuffler GUI
from pyrocko import obspy_compat
obspy_compat.plant()

# Add path to store waveforms (mseed)
bank_base_path = os.path.join(root, 'data', 'waveforms', 'BANK')
if not os.path.exists(bank_base_path):
    os.makedirs(bank_base_path)

# Establish connection to event bank
ebpath = os.path.join(root, 'data', 'XML', 'QUAKE', 'BANK')
evidpath = os.path.join(root, 'results', 'tables', 'well_located_mt_baker_evids_20km.csv')
ebank = cutil.connect_to_eventbank(base_path=ebpath)
df = ebank.read_index()
df_ = pd.read_csv(evidpath, header=None)
df_ = df_.drop(index=0) # the first row is not number 
df_ = df_[0].apply(lambda x: f'uw{int(x)}')
df = df[df.index.isin(df_.values)].event_id
cat = ebank.get_events(event_id = df.values)
print(f'connected to ebank: {type(ebank)}')

# Create an IRIS webservices client
client = Client('IRIS')

for event in cat:
    try:
        origin = event.preferred_origin()
        uw_evid = f"uw{str(event.resource_id).split('/')[-1]}"
        savename = os.path.join(bank_base_path, f'{uw_evid}.mseed')
        st = cutil.origin_bulk_waveform_request(origin, client, all_components=True)
        st.write(savename, format='MSEED')
        print(f'Successfully saved {savename}')
    except FDSNNoDataException:
        print(f'FDSN could not get data for {event.resource_id}')
        continue
    
# Below is how to use Snuffler with mseed files

# import subprocess
# from pathlib import Path
# example file = uw61280137.mseed
# mseed_file_path = Path('processed_data/uw61280137.mseed')
# subprocess.run(['snuffler', str(mseed_file_path)])
