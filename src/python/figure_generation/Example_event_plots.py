import os
import matplotlib.pyplot as plt
import pandas as pd
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

# Map Data Sources
root = os.path.join('..','..','..')
data = os.path.join(root,'data')
events = os.path.join(data,'Events','MtBaker_50km_radius_origins.csv')
phases = os.path.join(data,'Events','MtBaker_20km_radius_phases.csv')
stations = os.path.join(data,'Sensors','Mt_Baker_1deg_MDA_Stations.csv')
# Mt. Baker Summit Position
mbs_lat = 48.7745
mbs_lon = -121.8172
# Search Parameters
sta_rad_m = 30e3        # [m] search radius for stations
eve_rad_m = 20e3        # [m]

# Initalize IRIS Client
client = Client('IRIS')

# Load event data
eq_df = pd.read_csv(events)
# Rename a few headers
eq_df = eq_df.rename(columns={'to_timestamp':'origin_time'})
# Convert all origin times to UTCDateTime
# eq_df.origin_time = eq_df.origin_time.apply(lambda x: UTCDateTime(x))

# # Load Phase Data
# ph_df = pd.read_csv(phases)
# ph_df = ph_df.rename(columns={'to_timestamp':'origin_time','to_timestamp.1':'arrival_time'})

select_evids = {'deep LF': [60493937, 10242483],
                'vt EQ': [6168203, 61915302],
                'shallow LF': [61958071]}
# Get origin times into UTCDateTime format
holder = {}
for key, value in select_evids.values():
    for evid in value:
        t0 = eq_df[eq_df.evid == evid].origin_time
        inv = client.get_stations(startbefore=t0 - 60, endafter=t0 + 300,latitude=mbs_lat, longitude=mbs_lon, maxradius=km2deg(sta_rad_m*1e-3))

