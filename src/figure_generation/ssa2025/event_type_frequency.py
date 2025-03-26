"""
Event Type Frequency and Network Evolution Time-Series for the Mount Baker area

Figure shows the monthly frequency of the four catalog event types (excluding one UnKnown event)
cataloged within 30 km of Mount Baker between 1980 and early 2025. 

TODO: Color coordinate event types across figures

"""

import os
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from obsplus import EventBank
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent.parent
# Absolute path to event bank base_path
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Absolute path to EVent IDentifier - Event TYPE association table
ETYPE = ROOT / 'data' / 'Events' / 'Mount_Baker_evids_etypes_10_JAN_2025.csv'
# Absolute path to catalog membership table output by Phase 1 Step 1 script
CATM = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'

# Reference Latitude / Longitude / Elevation [m] of Mount Baker's summit
LAT, LON, ELE = 48.7745,-121.8172, 3286

# Absolute Path and File name structure for saving this figure
SPATH = ROOT / 'processed_data' / 'SSA2025' / 'figures'
SNAME = 'event_frequency_figure_{dpi}dpi.{fmt}'
# Rendering parameters
DPI = 120       # Figure rendering Dots Per Inch
FMT = 'png'     # Save Format
issave = False  # Should a figure be saved to disk?
isshow = True   # Should the figure be displayed?


### LOADING SECTION ###
# Load etype styles shared formatting
_styles = pd.read_csv(Path(__file__).parent / 'etype_styles.csv', index_col=[0])
# Connect to event bank client
# EBANK = EventBank(EBBP)

# Get IRIS webservice client
IRIS = Client("IRIS")
# Get stations within 1 degree of Mount Baker
inv = IRIS.get_stations(longitude=LON,
                        latitude=LAT,
                        maxradius=75/111.2,
                        level='station',
                        network='UW,CN,TA')
inv30 = IRIS.get_stations(longitude=LON,
                        latitude=LAT,
                        maxradius=30/111.2,
                        level='station',
                        network='UW,CN,TA')

# # Read event bank index
# df = EBANK.read_index()
# df.index = df.event_id

### PLOTTING SECTION ###

# Append catalog membership
df = pd.read_csv(CATM, index_col='event_id', parse_dates=['prefor_time'])

# breakpoint()
# df = pd.concat([df, df_cat], axis=1, ignore_index=False)


df0 = df[df.CAT0]
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)

# Create frequency time-series for event types
# Flag to create inventory-based active station count time-series
_run_inv = True
for _e in ['eq','lf','su','px']:
    # Subset by event type
    _df = df0[df0.etype == _e]
    # Create holders for event frequency counts and reference times
    x = []
    y = []
    # For first run across times, 
    # create the active station count time series
    if _run_inv:
        z = []
        z30 = []
    # Iterate across years
    for _yr in range(1980, 2026, 1):
        # Iterate across months
        for _mo in range(1,12,1):
            # Generate time bounds
            t0 = pd.Timestamp(f'{_yr}-{_mo:02d}-01')
            t1 = pd.Timestamp(f'{_yr}-{_mo+1:02d}-01')
            # Get event count for a (t0, t1] range - don't double count event if it sits at t0
            ct = len(_df[(_df.prefor_time > t0) & (_df.prefor_time <= t1)])
            # Capture event count and timestamps
            y.append(ct)
            x.append(t1)
            # If creating an inventory count of active stations
            if _run_inv:
                # Get count for this month at larger radius
                z.append(len(inv.select(time=UTCDateTime(f'{_yr}-{_mo:02d}-01')).get_contents()['stations']))
                # Get count for this month at smaller radius
                z30.append(len(inv30.select(time=UTCDateTime(f'{_yr}-{_mo:02d}-01')).get_contents()['stations']))
    # Turn off inventory-based profiling
    _run_inv = False
    # Plot frequencies
    # ax.fill_between(x, [0]*len(y), y, label=f'"{_e.upper()}" Frequency (Total Ct.: {len(_df)})', alpha=0.75)
    ax.plot(x,y, color=_styles[_e].color, label=f'"{_e.upper()}" Frequency (Total Ct.: {len(_df)})')

# Plot station counts
ax.plot(x, z,'m',label='Stations Within 75 km')
ax.plot(x,z30, color='darkgoldenrod', label='Stations Within 30 km')
ax.plot(x, [4]*len(x), ':', color='grey', alpha=0.5)
# Annotations & Axis Labels
plt.legend(loc='best', ncols=2)
ax.set_xlabel('Year', fontsize=12)
ax.set_ylabel('Frequency (ct./month) | Count (ct.)', fontsize=12)
ax.set_xlim([min(x), max(x)])
ax.set_ylim([0, 45])
### SAVE / DISPLAY SECTION ###
if issave:
    try:
        os.makedirs(SPATH, existsok=False)
    except:
        pass
    fig.savefig(str(SPATH/SNAME).format(dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)
if isshow:
    plt.show()