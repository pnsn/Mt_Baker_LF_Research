"""
:script: DLP_Mt_Baker/src/Python_Scripts/figure_generation/Motivation_Fig_1.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPL v3
:purpose: Plot frequency, depth, and magnitude time-series of deep long-period events in the PNSN
    catalog within a specified radius of Mt. Baker volcano. Accompanies the Motivation section in
    this repository's READ_ME.md.
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth

# PATH MAPPING
# Relative Root Directory Path
rrpath = os.path.join('..','..','..')
# PNSN Catalog Events for Mt. Baker Area of Observation
pnsn_events = os.path.join(rrpath,'data','GIS','Points','PNSN_CAT_Mount_Baker_Region_Full_Catalog_Event_Summary.csv')
local_stations = os.path.join(rrpath,'data','GIS','Points','Mt_Baker_AO_MDA_Stations.csv')
# OUTPUT CONTROL
issave = True
isshow = False
DPI = 120
FMT = 'PNG'
# Save Path
sfpath = os.path.join(rrpath,'results','figures')
if not os.path.exists(sfpath):
    os.makedirs(sfpath)
savename = os.path.join(sfpath, 'motivations_fig_1_%ddpi.%s'%(int(DPI), FMT.lower()))

# SEARCH CONTROL
# Catalog Subset SearchParameters
eq_rad_m = 20e3         # [m] search radius to include events as "Mt. Baker" events
sta_rad_m = 30e3        # [m] inclusion radius for stations around Mt. Baker
mbs_lat = 48.7745       # [degrees N] latitude of Mt. Baker Summit
mbs_lon = -121.8172     # [degrees E] longitude of Mt. Baker Summit
dlp_min_depth = 10e3    # [m] minimum source depth to classify event as DLP (after Nichols et al., 2011)
include_etypes = ['lf'] # [str] AQMS event type strings to include in search
n11_cutoff = pd.Timestamp('2009-01-01') # Approximate cut-off date for Nichols et al. (2011) catalog

# Load CSV output by AQMS
eq_df = pd.read_csv(pnsn_events)
# Load CSV from IRIS MDA
sta_df = pd.read_csv(local_stations)
# Ensure timestamps are pandas Timestamp objects and convert from local time to UTC
eq_df.timestamp = eq_df.timestamp.apply(lambda x: pd.Timestamp(x).tz_convert(None))
# Rename lat lon columns
sta_df = sta_df.rename(columns={'Longitude': 'lon', 'Latitude': 'lat'})
# Format StartTime and EndTime
sta_df.StartTime = sta_df.StartTime.apply(lambda x: pd.Timestamp(x))
sta_df.EndTime = sta_df.EndTime.apply(lambda x: pd.Timestamp(x))

# Run geodetic calcs for eq's and stations
for ct, df in enumerate([eq_df, sta_df]):
    df_geod = pd.DataFrame(data = [gps2dist_azimuth(mbs_lat, mbs_lon, row.lat, row.lon) for _, row in df.iterrows()],
                        columns=['lateral_offset_m','mbs_eq_az','eq_mbs_az'])
    # # Append geodetic calcs to dataframe
    cols = df.columns.copy()
    df = pd.concat([df, df_geod], ignore_index=True, axis=1)
    # DEBUG: preserve column names internal to concat
    df = df.rename(columns=dict(zip(range(len(df.columns)), list(cols) + list(df_geod.columns))))
    if ct == 0:
        eq_df = df
    elif ct == 1:
        sta_df = df
# Apply Distance, Depth, and Event Type Filters
eq_dff = eq_df[(eq_df.lateral_offset_m <= eq_rad_m) &\
        (eq_df.etype.isin(include_etypes)) &\
        (eq_df.depth >= dlp_min_depth*1e-3)]

sta_dff = sta_df[(sta_df.lateral_offset_m <= sta_rad_m)]


time_index = pd.date_range(start=pd.Timestamp('1979-01-01'), end=pd.Timestamp('2026-01-01'), freq='1M')
sta_count = [len(sta_dff[(sta_dff.StartTime <= time_index[x]) &\
                         (sta_dff.EndTime >= time_index[x+1])])
             for x in range(len(time_index) - 1)]

# Create additional filter for if events should be in the catalog of Nichols et al. (2011)
n11_ind = (eq_dff.timestamp >= pd.Timestamp('1980-01-01')) & (eq_dff.timestamp < n11_cutoff)

# Generate Figure
# Initialize figure and axes
fig = plt.figure(figsize=(6,4))
ax = fig.add_subplot(111)
# Plot Nichols et al. (2011) and subsequent DLP's
eq_dff[n11_ind].plot.scatter(x='timestamp',y='depth',ax=ax,c='blue', zorder=10)
eq_dff[~n11_ind].plot.scatter(x='timestamp',y='depth',ax=ax,c='orange', zorder=10)
# Generate legend
plt.legend(['1980-2009: %d DLP'%(sum(n11_ind)),
            '2010-present: %d DLP'%(sum(~n11_ind))],
            loc='upper left')
# Format Axes
ax.set_ylim([40, 5])
ax.set_xlabel('Origin Time [year]')
ax.set_ylabel('Source Depth [km b.s.l.]')
# Twin x-axis to produce another y-axis
ax2 = ax.twinx()
# Plot monthly station counts 
ax2.plot(time_index[1:], sta_count,'r-', zorder=1, alpha=0.5)
# Add second y-axis label
ax2.set_ylabel(f'Active Stations within {int(sta_rad_m*1e-3)} km of Mt. Baker\nAnnual DLP Frequency',
               rotation=270, labelpad=20)
ax2.hist(eq_dff.timestamp,bins=pd.date_range(start=time_index[0], end=time_index[-1], freq='12M'),
         alpha=0.5, color='black')


if isshow:
    plt.show()

if issave:
    fig.savefig(savename, dpi=DPI, format=FMT)