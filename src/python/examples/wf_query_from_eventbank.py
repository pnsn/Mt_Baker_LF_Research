import os, sys, logging
from pathlib import Path

from obspy.clients.fdsn import Client

# Get file path on your machine for Mt_Baker_LF_Research/src/python
pyroot = Path(__file__).parent.parent
# Add that path to your current path
sys.path.append(pyroot)
# Import catalog
import src.python.eqc_utils.catalog_utils as cutil


# establish connection to event bank
ebank = cutil.connect_to_eventbank()
print(f'connected to ebank: {type(ebank)}')
# Create a wildcard based query for a set of events and include fixed status
df = ebank.read_index(event_id='*/UW/6204*', include_fixed_status=True)
print(f'subset query from wildcard search has {len(df)} entries')
# Add a temporary column to see if anything is fixed
df = df.assign(anyfix=[any([x.fepi,x.ftime,x.fdepth]) for _, x in df.iterrows()])

# Subset to just have non-fixed events
df_filt = df[~df.anyfix]
print(f'subset query for non-fixed events has {len(df_filt)} entries')
# Get all these events as a catalog from the EventBank
cat = ebank.get_events(event_id=df_filt.event_id)
print('events pulled from EventBank')

# Create a bulk data request from the first event return
event = cat[0]
origin = event.preferred_origin()

# Create an IRIS webservices client
client = Client('IRIS')
# Run query
st = cutil.origin_bulk_waveform_request(origin, client)
breakpoint()
# bulk = cutil.compose_origin_bulk_lines(origin, method='pick', lead_time=10, lag_time=60, all_components=True)



# print('fetching data')
# st = client.get_waveforms_bulk(bulk)

# st.plot(equal_scale=False, method='full')