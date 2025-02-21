"""


:helpful links:
 - Open Source Terrain Imagery from AWS: https://registry.opendata.aws/terrain-tiles/
    - Specifically for Mapzen Terrain: https://docs.safe.com/fme/html/FME-Form-Documentation/FME-ReadersWriters/terraintilesaws/terraintilesaws.htm
 - PIL Image Conversion Reference for `desired_tile_format`: https://pillow.readthedocs.io/en/stable/handbook/concepts.html#concept-modes
    - look at the "Modes" section - here we use 8-bit pixels to match the native elevation valuse of Normal (PNG) Mapzen tiles
"""
import datetime
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize


import cartopy.crs as ccrs
from cartopy.io.img_tiles import OSM, GoogleTiles, StadiaMapsTiles

import pandas as pd
from obsplus import EventBank
from map_util import mount_baker_basemap, rad2llur

# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to catalog membership CSV
CATD = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'
# path to preferred event/channel pairs CSV
PESD = ROOT / 'processed_data' / 'catalog' / 'P1S2_Preferred_Sta_Event_Picks.csv'
# path to analyst review summary CSV
REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'

map_rad_km = 30.

bounds = rad2llur(rad = map_rad_km*1e3)


## PREPROCESSING SECTION ###
# Connect to Event Bank
EBANK = EventBank(EBBP)


# Read EBANK Index
df_eb = EBANK.read_index()
# Merge with catalog membership
df_eb.index = df_eb.event_id
df_eb = df_eb.join(pd.read_csv(CATD, index_col='event_id'), how='left')
# breakpoint()



## PLOTTING SECTION ##
data_crs = ccrs.PlateCarree()
# Initialize figure
fig = plt.figure(figsize=(8.5, 6.5))
# Initalize GridSpec
gs = fig.add_gridspec(ncols=4, nrows=3)
# Initialize Basemap & Axis
axmap, attribution = mount_baker_basemap(
    fig=fig, sps=gs[:2,:2], radius_km=map_rad_km, zoom=10,
    open_street_map=False)
axes = {'map': axmap}
# Create additional axes
# axes = [fig.add_subplot(gs[:2,:2], projection=imagery.crs)]
axes['b'] = fig.add_subplot(gs[2,:2])
axes['c'] = fig.add_subplot(gs[0,2:])
axes['d'] = fig.add_subplot(gs[1,2:])
axes['e'] = fig.add_subplot(gs[2,2:])




# Catalog Color Marker Edge Mapping
corder = ['maroon','goldenrod','navy','cyan']
# sorder = [2.75, 2.5, 2.25, 2]

# Event Type Shape and Face Color Mappings
shapemap = {'lf': 's', 'eq': '*', 'su': 'o', 'px': 'x','uk': '.'}
facemap = {'lf':'k', 'eq': 'r', 'su': 'c', 'px': 'g', 'uk': 'm'}
for _e, catname in enumerate(['CAT0','CAT1','CAT2','WC']):
    for _k, _v in shapemap.items():
        _dfe = df_eb[(df_eb[catname]) & (df_eb.etype == _k)]
        # Also plot excluded >30km radius, neighboring events
        if catname == 'CAT0':
            _dfx = df_eb[~df_eb[catname]]
            axes['map'].plot(_dfx.longitude.values, _dfx.latitude.values,
                        'k.', ms=.25, transform=ccrs.PlateCarree())
            axes['b'].plot(_dfx.longitude, _dfx.depth*1e-3, 'k.',
                        ms=0.25)

        _kw = {'color': facemap[_k], 'marker': _v, 'linestyle': '',
                'ms': 2, 'label': catname, 
                'markeredgecolor': corder[_e],
                'markeredgewidth': 0.3}

        axes['map'].plot(_dfe.longitude.values, _dfe.latitude.values,
                    transform=ccrs.PlateCarree(), **_kw)
        
        axes['b'].plot(_dfe.longitude, _dfe.depth*1e-3,**_kw)
        axes['c'].plot(_dfe.time, _dfe.depth*1e-3,**_kw)

        __dfe_r = _dfe.copy()
        __dfe_r.index = __dfe_r.time
        _r = __dfe_r.etype.resample(pd.Timedelta(1, unit='w')).count()
        if catname in ['CAT0','CAT1']:
            _x = 'd'
        else:
            _x = 'e'
        
        axes[_x].fill_between(_r.index, [0]*len(_r), _r.values,
                                color=facemap[_k], alpha=0.5)
        

# Flip depth axes
axes['b'].set_ylim([axes['b'].get_ylim()[1], -2])

axes['c'].set_ylim([axes['c'].get_ylim()[1], -2])

# Match scope of 'b' to 'map'
axes['b'].set_xlim([bounds[0], bounds[1]])

# Shift axes to right for c-e
for _c in 'cde':
    axes[_c].yaxis.set_ticks_position(position='right')
    axes[_c].yaxis.set_label_position(position='right')
# Label Axes
axes['b'].set_ylabel('Depth [km]')
axes['b'].set_xlabel('Longitude [$^oE$]')
axes['c'].set_ylabel('Depth [km]', rotation=270, labelpad=15)
axes['d'].set_ylabel('Events per week [no.]')
axes['e'].set_ylabel('Events per week [no.]')


# Add attribution to map
xlims = axes['map'].get_xlim()
ylims = axes['map'].get_ylim()
# axes['map'].text(xlims[0], ylims[1], '\n'.join(attribution), fontsize=6, ha='left',va='top')
# Show
plt.show()