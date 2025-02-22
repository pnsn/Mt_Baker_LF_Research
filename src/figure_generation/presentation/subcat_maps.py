import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from obsplus import EventBank
import cartopy.crs as ccrs

from map_util import mount_baker_basemap, rad2llur

# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to catalog membership CSV
CATD = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'
# path to preferred event/channel pairs CSV
PESD = ROOT / 'processed_data' / 'catalog' / 'P1S2_Preferred_Sta_Event_Picks.csv'
# path to analyst review summary CSV
REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'

# SAVEPATH
SAVEPATH = ROOT / 'results' / 'figures' / 'seismolunch'

map_rad_km = 30.

bounds = rad2llur(rad = map_rad_km*1e3)

## PREPROCESSING SECTION ###
# Connect to Event Bank
EBANK = EventBank(EBBP)
try:
    os.makedirs(str(REVD), exist_ok=False)
except:
    pass

# Read EBANK Index
df_eb = EBANK.read_index()
# Merge with catalog membership
df_eb.index = df_eb.event_id
df_eb = df_eb.join(pd.read_csv(CATD, index_col='event_id'), how='left')

# Iterate across event types
# TODO: Lump together PX and EX
# TODO: Apply vmin-vmax for entire timeframe to all plots
for etype in df_eb.etype.unique():

    _df = df_eb[(df_eb.etype == etype) & (df_eb.CAT0)]

    fig = plt.figure(figsize=(7*1.1, 6.9*1.1)) 
    gs = fig.add_gridspec(ncols=4, nrows=4)
    axes = {}
    # Plot basemap as part of subaxis initialization
    axes['map'], map_attrib = mount_baker_basemap(
        fig=fig, sps=gs[:3,:3], radius_km=map_rad_km)
    gl = axes['map'].gridlines(
        draw_labels=True,
        zorder=1)
    gl.bottom_labels=False
    gl.right_labels=False
    axes['map'].scatter(_df.longitude, _df.latitude,
                        s=1 + 2**_df.magnitude, c=_df.time,
                        cmap='inferno_r',
                        transform=ccrs.PlateCarree(),zorder=3)
    axes['map'].scatter(_df.longitude, _df.latitude,
                        s=_df.horizontal_uncertainty*1e-3, c='red',
                        transform=ccrs.PlateCarree(),zorder=2)
    axes['map'].set_title(etype)

    # Initialize side plots
    axes['b'] = fig.add_subplot(gs[3,:3])
    axes['b'].scatter(_df.time, _df.depth*1e-3,
                       s=1 + 2**_df.magnitude, c=_df.time,
                       cmap='inferno_r')
    # axes['b'].set_xlim(bounds[:2])
    # TODO: Flip y axis
    # TODO: Insert frequency plot as subplot



    axes['c'] = fig.add_subplot(gs[:3, 3])
    ch = axes['c'].scatter(
        _df.vertical_uncertainty*1e-0,
        _df.latitude,
        c=_df.time,
        s=1 + 2**_df.magnitude,
        cmap='inferno_r', zorder=3)
    axes['c'].scatter(
        _df.vertical_uncertainty*1e-0,
        _df.latitude,
        s=_df.horizontal_uncertainty*1e-3,
        c='red',zorder=2)
    axes['c'].set_xscale('log')
    axes['c'].yaxis.set_ticks_position('right')
    axes['c'].yaxis.set_label_position('right')
    axes['c'].set_ylabel('Latitude [$^oN$]', rotation=270, labelpad=15)
    axes['c'].set_xlabel('Depth Uncertainty [m]')
    axes['c'].set_ylim(bounds[2:])
    # fig.colorbar(ch, use_gridspec=True, cax=fig.add_subplot(gs[-1,-1]))

# TODO - Summative figure with all etypes

plt.show()
    # breakpoint()