"""
:script: DLP_Mt_Baker/src/Python_Scripts/figure_generation/Mt_Baker_Events_Stations_Timeseries.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPL v3

.. rubric:: Purpose
    Plot frequency, depth, and magnitude time-series of deep long-period events in the PNSN
    catalog within a specified radius of Mt. Baker volcano. Accompanies the Motivation section in
    this repository's READ_ME.md.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from obspy.geodetics.base import gps2dist_azimuth

# Define small geodetic function
def my_km_distance(lat0,lon0,lat1,lon1):
    offset_m, _, _ = gps2dist_azimuth(lat0, lon0, lat1, lon1)
    return offset_m*1e-3


# PATH MAPPING
# Relative Root Directory Path
rrpath = os.path.join('..','..','..')
# PNSN Catalog Events for Mt. Baker Area of Observation
events = os.path.join(rrpath,'data','Events','MtBaker_50km_radius_origins.csv')
stations = os.path.join(rrpath,'data','Sensors','Mt_Baker_1deg_MDA_Stations.csv')
mtbaker = os.path.join(rrpath,'data','POI','Mt_Baker_Summit_LLH.csv')
# OUTPUT CONTROL
issave = True
isshow = False
DPI = 120
FMT = 'PNG'
# Save Path
sfpath = os.path.join(rrpath,'docs','Figures')
if not os.path.exists(sfpath):
    os.makedirs(sfpath)
savename = os.path.join(sfpath, 'Mt_Baker_Catalog_and_Stations_Timeseries_%ddpi.%s'%(int(DPI), FMT.lower()))

# SEARCH CONTROL
# Catalog Subset SearchParameters
eq_rad_km = 10.         # [km] search radius to include events as "Mt. Baker" events
sta_rad_km = 30.        # [km] inclusion radius for stations around Mt. Baker
nowtime = pd.Timestamp.now() # Timestamp for RIGHT NOW! in UTC
n11_depth_min = 10e3    # [m] minimum depth requirement for LF == DLP per Nichols et al. (2011)
n11_cutoff = pd.Timestamp('2009-10-01') # Catalog cutoff for deep LF events published in Nichols et al. (2011)


# Create time indices for getting monthly counts of stations and yearly counts of event types
month_index = pd.date_range(start=pd.Timestamp('1979-01-01'), end=nowtime, freq='1M')
year_index = pd.date_range(start=pd.Timestamp('1979-01-01'), end=nowtime, freq='12M')
year_midpoints = [year_index[x:x+2].mean() for x in range(len(year_index) - 1)]
# Load location information for Mt. Baker summit
mbs_df = pd.read_csv(mtbaker)

## EVENT PROCESSING ##
# Load Events from AQMS query output
eq_df = pd.read_csv(events)
# Update column name
eq_df = eq_df.rename(columns={'to_timestamp':'origin_time'})
# Filter down to events within the radius specified
eq_df = eq_df[eq_df.mbs_distance_km <= eq_rad_km]
# Convert origin_times to UTC Timestamp objects
eq_df.origin_time = eq_df.origin_time.apply(lambda x: pd.Timestamp(x).tz_localize(None))

# Create monthly event counts
event_counts = {}
for _e in eq_df.etype.unique():
    # Subset by event type
    eq_dfe = eq_df[eq_df.etype==_e]
    e_count = [len(eq_dfe[(eq_dfe.origin_time >= year_index[x]) &\
                           (eq_dfe.origin_time < year_index[x+1])])
                for x in range(len(year_index) - 1)]
    event_counts.update({_e:e_count})
# Create dataframe with datetime index
event_counts_df = pd.DataFrame(event_counts, index=year_index[1:])

## STATION PROCESSING
# Load Stations from IRIS Metadata Aggregator Query
sta_df = pd.read_csv(stations)
# Make filter for stations by distance from Mt. Baker
sta_dist = [my_km_distance(mbs_df.lat[0], mbs_df.lon[0], row.Latitude, row.Longitude) <= sta_rad_km 
            for idx, row in sta_df.iterrows()]
             #sta_df.Network.isin(['UW']) & (sta_df.Station != 'MBW2')]
sta_df = sta_df[sta_dist]
# Convert StartTime and Endtime into Timestamp objects
sta_df.StartTime = sta_df.StartTime.apply(lambda x: pd.Timestamp(x))
sta_df.EndTime = sta_df.EndTime.apply(lambda x: pd.Timestamp(x))
# Convert future-ending stations into now-ending stations
sta_df.EndTime = sta_df.EndTime.apply(lambda x: nowtime if x > nowtime else x)
# Make filter for stations with runtime > 2 years or has a future-dated end-time
sta_runtime = [any(
                    [(row.EndTime - row.StartTime) > pd.Timedelta(2*365.24, unit='d'),
                     row.EndTime == nowtime]
                ) for _, row in sta_df.iterrows()]
# Apply runtime filter
# sta_df = sta_df[sta_runtime]
# Exclude AM network (Raspberry Shake) for now... 
sta_df = sta_df[sta_df.Network == 'UW']
# Sort times for station dataframe
sta_df = sta_df.sort_values('StartTime')
                
# Create monthly station counts
sta_count = [len(sta_df[(sta_df.StartTime <= month_index[x]) &\
                         (sta_df.EndTime >= month_index[x+1])])
             for x in range(len(month_index) - 1)]
sta_count_series = pd.Series(sta_count,index=month_index[1:])

### PLOTTING ###

fig = plt.figure(figsize=(6,4))
axes = [fig.add_subplot(111)]
axes.append(axes[0].twinx())
# gs = fig.add_gridspec(ncols=1, nrows=2)
# axes = [fig.add_subplot(gs[0])]
# axes.append(fig.add_subplot(gs[1],sharex=axes[0]))

# Plot Event Type Frequencies(stacked)
etype_list = [('orange','lf'),('red','eq'),('blue','su'),('black','px')]#,('purple','uk')]
for _c, _e in etype_list:
    idf = event_counts_df[_e]#.sum(axis=1)
    axes[0].fill_between(idf.index,idf.values,color=_c,label=_e.upper(), alpha=0.5)
# Plot station counts
axes[1].plot(sta_count_series.index, sta_count_series.values, 'm-', label='Station Count')

# Make legend
axes[0].legend(loc='upper left')
# Label Axes
axes[0].set_xlabel('Year')
axes[1].set_ylabel('Network UW stations within %d km of Mt. Baker\n Monthly Count'%(sta_rad_km),
                   rotation=270, labelpad=20, color='m')
axes[0].set_ylabel('PNSN catalog events within %d km of Mt. Baker\nAnnual Count'%(eq_rad_km))

# Manual sets to x and y limits
axes[0].set_xlim([pd.Timestamp('1980-01-01'), event_counts_df.index[-1]])
axes[0].set_ylim([0,80])
axes[1].set_ylim([0,10])

if issave:
    fig.savefig(savename, dpi=DPI, format=FMT)
if isshow:
    plt.show()
