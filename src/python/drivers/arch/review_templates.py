import os, glob, sys
from pathlib import Path

from eqcorrscan import Template, Tribe
from obspy import Stream, read_inventory
from pyrocko import obspy_compat

obspy_compat.plant()

root = Path(__file__).parent.parent.parent.parent
ebpath = os.path.join(root,'data','XML','QUAKE','BANK')
invpath = os.path.join(root, 'data','XML','INV','station_inventory_50km_MBS_IRISWS_STATION.xml')
tempath = os.path.join(root,'processed_data','templates','well_located')
pypath = os.path.join(root,'src','python')
savedir = os.path.join(root,'processed_data','templates','well_located')
sys.path.append(pypath)
# Import utility Scripts
from eqc_utils.catalog_utils import connect_to_eventbank

# Get list of template files
temp_files = glob.glob(os.path.join(tempath, '*.tgz'))

# Load templates
tribe = Tribe()
for tfile in temp_files:
    tribe += Template().read(tfile)

# Create a big stream
big_st = Stream()
for template in tribe:
    big_st += template.st

# Get all event_id's
evids = [template.name for template in tribe]

# Connect to EventBank
ebank = connect_to_eventbank(ebpath)

# Get EventBank index
df_evid = ebank.read_index()
# Filter down to events
df_wlevid = df_evid[df_evid.index.isin(evids)]

# Load Events
cat = ebank.get_events(event_id=df_wlevid.event_id.values)

# Load inventory
inv = read_inventory(invpath)

# Run snuffler with station and event metadata
big_st.snuffle(catalog=cat, inventory=inv)



