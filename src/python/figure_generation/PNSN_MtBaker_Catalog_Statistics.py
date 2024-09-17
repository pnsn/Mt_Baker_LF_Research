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
from copy import deepcopy
import numpy as np

# PATH MAPPING
# Relative Root Directory Path
rrpath = os.path.join('..','..','..')
# PNSN Catalog Events for Mt. Baker Area of Observation
pnsn_events = os.path.join(rrpath,'data','Events','MtBaker_50km_radius_origin_query.csv')

# Catalog filters
search_radius_km = 20.

# OUTPUT CONTROL
issave = True
isshow = True
DPI = 120
FMT = 'PNG'
# Save Path
sfpath = os.path.join(rrpath,'results','figures')
if not os.path.exists(sfpath):
    os.makedirs(sfpath)

# Read CSV
df = pd.read_csv(pnsn_events)
# Format datetime as UTC
df.origin_time = df.origin_time.apply(lambda x: pd.Timestamp(x).tz_convert(None))
# Subset data by distance from MBS
df = df[df.calculated_distance_km <= search_radius_km]
# Set index to origin times
df.index = df.origin_time
# Sort by origin time
df = df.sort_index()

df_zetype = df.copy()
df_zetype = df_zetype.assign(zetype=['LF ge 10 km' if row.etype == 'lf' and row.depth >= 10
                                     else 'LF lt 10 km' if row.etype == 'lf' and row.depth < 10
                                     else row.etype.upper() for _, row in df_zetype.iterrows()])

time_index = pd.date_range(start=pd.Timestamp('1979-01-01'), end=pd.Timestamp('2026-01-01'), freq='12M')
# zet_count = [df_zetype[(df_zetype.origin_time >= time_index[x]) & (df_zetype.origin_time <= time_index[x+1])].zetype.value_counts().values
#              for x in range(len(time_index) - 1)]
# breakpoint()
# Figure 1 - time-series of etype frequency
fig = plt.figure(figsize=(12,8))
gs = fig.add_gridspec(nrows=3,ncols=5, hspace=0)
axes = [fig.add_subplot(gs[0,:3])]
axes += [fig.add_subplot(gs[_e+1, :3], sharex=axes[0]) for _e in range(2)]
axes.append(fig.add_subplot(gs[:,3:]))

colors = ['dodgerblue','red','grey','orange','magenta']
ordered_zet = ['LF ge 10 km','LF lt 10 km','EQ','SU']
for _z in df_zetype.zetype.unique():
    if _z not in ordered_zet:
        ordered_zet.append(_z)
for _e, zet in enumerate(ordered_zet):
    if _e == 0:
        pcounts = np.zeros(len(time_index) - 1)
    elif _e > 0:
        pcounts = deepcopy(tcounts)
    subset = df_zetype[df_zetype.zetype == zet]
    tcounts = [len(subset[(subset.index >= time_index[x]) &
                          (subset.index <= time_index[x + 1])])
            for x in range(len(time_index) - 1)]
    tcounts = np.array(tcounts)
    # Subplot A) Time-Frequency-Event Type Distribution
    axes[0].fill_between(time_index[1:], pcounts, tcounts + pcounts, label=zet, color=colors[_e])
    tcounts = tcounts +  pcounts
    axes[0].legend()
    axes[0].set_ylabel('Cumulative\nFrequency [$yr^{-1}$]')
    axes[0].set_title('PNSN Catalog Events Within %.1f km of Mt. Baker Summit'%(search_radius_km))
    # axes[0].set_xlabel('Origin Year')
    # Subplot B) 
    axes[1].plot(subset.index, subset.depth, 'o', label=zet, color=colors[_e], zorder=_e+3)
    axes[1].set_ylim([40,-2])
    axes[1].set_ylabel('Origin Depth [km b.s.l.]')
    # axes[1].set_xlabel('Origin Year')
    # axes[1].legend(ncols=2, loc='upper left')
    axes[1].fill_between([pd.Timestamp("1980-01-01"), pd.Timestamp('2009-01-01')], [40,40], [10,10], color='firebrick', alpha=0.05)
    axes[1].text(pd.Timestamp("1985-01-01"), 38, 'Nichols et al. (2011) DLP Catalog', color='black')

    # Subplot C) Time magnitude
    axes[2].plot(subset.index, subset.magnitude, 'o', label=zet, color=colors[_e], zorder=_e)
    # axes[2].set_ylim([40,-2])
    axes[2].set_ylabel('Magnitude')
    axes[2].set_xlabel('Origin Year')
    axes[2].legend(ncols=2, loc='lower left')

    # Subplot D) Depth-Magnitude
    axes[3].plot(subset.magnitude, subset.depth, 'o', label=zet, color=colors[_e], zorder=_e)
    axes[3].yaxis.set_ticks_position('right')
    axes[3].yaxis.set_label_position('right')
    axes[3].set_ylim([40, -2])
    axes[3].set_ylabel('Origin Depth [km b.s.l.]')
    axes[3].set_xlabel('Magnitude')    


    for _x in range(4):
        axes[_x].grid(visible=True, linestyle=':')

    axes[0].text(pd.Timestamp('1980-01-01'), 10, 'A', fontstyle='italic',fontweight='extra bold', fontsize=14)
    axes[1].text(pd.Timestamp("1980-01-01"), 5, 'B', fontstyle='italic',fontweight='extra bold', fontsize=14)
    axes[2].text(pd.Timestamp("1980-01-01"), 1, 'C', fontstyle='italic',fontweight='extra bold', fontsize=14)
    axes[3].text(-5, 38, 'D', fontstyle='italic',fontweight='extra bold', fontsize=14)



if isshow:
    plt.show()
if issave:
    savename = f'Mt_Baker_PNSN_Catalog_Statistics_Plots_{DPI}dpi.{FMT.lower()}'
    fig.savefig(os.path.join(sfpath, savename))