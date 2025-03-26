import os
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
import numpy as np

from obsplus import EventBank
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER

from map_util import *

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
FMT = 'png'
DPI = 200

issave = True
isshow = False

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
# Create matplotlib timestamp-friendly times
df_eb = df_eb.assign(epoch=[x.timestamp() for x in df_eb.time])
# Get easting/northing
points = UTM10N.transform_points(
    x = df_eb.longitude.values,
    y = df_eb.latitude.values,
    src_crs = WGS84)
mE = points[:,0]
mN = points[:,1]
df_eb = df_eb.assign(mE = mE)
df_eb = df_eb.assign(mN = mN)

# Initialize Figure
fig = plt.figure(figsize=(10, 8))
gs = fig.add_gridspec(ncols=5, nrows=6)
# Initialize Map Subplot
axm, map_attrib = mount_baker_basemap(fig=fig, sps=gs[:-2, :3], radius_km=map_rad_km,
                                      open_street_map=False)
# Initialize depth / longitude plot
axe = fig.add_subplot(gs[:-2,3:])
# Initialize timeline
axt = fig.add_subplot(gs[-2:,:4])
# Initialize depth histogram
axz = fig.add_subplot(gs[-2:, -1])
gl = axm.gridlines(draw_labels=True, zorder=1)
gl.bottom_labels=False
gl.right_labels=False
gl.xlocator = mticker.FixedLocator(np.arange(-123, 121, 0.2))

# # Plot surrounding events w/o splits
# _dfx = df_eb[~df_eb.CAT0]
# axm.plot(_dfx.longitude, _dfx.latitude, 'k.', 
#          ms=2, alpha=0.5, transform=ccrs.PlateCarree())

# Iterate across event typesdf_
for _e, _c in [('eq','ko'), ('lf','rs'), ('su','bs'), ('px','co'), ('uk', 'mo')]:
    # Split marker designation for scatter calls
    sc = _c[0]
    sm = _c[1]
    # Subset event types
    _df = df_eb[(df_eb.etype == _e) & (df_eb.CAT0)]
    # Get number of events with incomplete/unavailable uncertainties
    _n = len(_df) - len(_df[(_df.vertical_uncertainty.notna()) & (_df.horizontal_uncertainty.notna())])
    # Number of events within 10 km
    _n = len(_df[_df.offset_km <= 10.])
    # Plot lat/lon/mag on map
    axm.scatter(_df.longitude, _df.latitude, c=sc, s=3**_df.magnitude,
                transform=ccrs.PlateCarree(), zorder=2, marker=sm, alpha=0.5)
    # Uncertainty Plot
    axe.plot(_df.horizontal_uncertainty*1e-3, _df.vertical_uncertainty*1e-3, _c,
             alpha=0.5, label=f'{_e.upper()}: {len(_df)} ({_n})')

    # Time\Depth\Mag\plot
    axt.scatter(_df.time, _df.depth*1e-3, c=sc,
                s=3**_df.magnitude, alpha=0.5, label=_e.upper())

    # Depth histogram
    axz.hist(_df.depth*1e-3, np.arange(-2.5,50,1.25), density=False,
             orientation='horizontal', color=sc,
             alpha=0.5)
    # _n = len(_df) - (_df[_df.horizontal_uncertainty.isna()])
    # _n
    # _n = sum(_df[(_df.vertical_uncertainty.isna()) |\
    #               (_df.horizontal_uncertainty.isna())])
    # breakpoint()


    # # axe.hist(_df.horizontal_uncertainty*1e-3, 50, log=True, orientation='vertical', density=True)
    # axe.hist(_df.vertical_uncertainty*1e-3, color=sc, alpha=0.5, orientation='horizontal')
    # axe.errorbar(_df.mE, _df.depth*1e-3, xerr=_df.horizontal_uncertainty,
    #              yerr=_df.vertical_uncertainty*1e-3, capsize=5, fmt='o', color=_c)

# Map Formatting
# Add attribution to map
xlims = axm.get_xlim()
ylims = axm.get_ylim()
axm.text(xlims[1] - 500, ylims[1] - 500, 
         '\n'.join(map_attrib), 
         fontsize=6, ha='right',va='top')


# Uncertainty Formatting
axe.set_xscale('log')
axe.set_yscale('log')
axe.grid(linestyle=':')
xlims = axe.get_xlim()
ylims = axe.get_ylim()
axe.plot(xlims, [31.2,31.2], 'm:', alpha=0.5)
axe.text(1.3e-2, 33, '31.2 km', color='m')
axe.set_xlim(xlims)
axe.set_ylim(ylims)
axe.yaxis.set_ticks_position('right')
axe.yaxis.set_label_position('right')
axe.set_ylabel('Vertical Uncertainty [km]', rotation=270, labelpad=15)
axe.xaxis.set_ticks_position('top')
axe.xaxis.set_label_position('top')
axe.set_xlabel('Horizontal Uncertainty [km]')
axe.legend()

# Timeline Formatting
axt.set_xlabel('Year')
axt.set_ylabel('Depth [km]')
ylims = list(axt.get_ylim())
ylims[1] = 45
axt.set_ylim(ylims[::-1])
axt.legend(ncols=5)

# Depth Histogram formatting
axz.set_ylim(ylims[::-1]) #Intentionally using axt's ylims
axz.set_xscale('log')
axz.set_xlabel('Event Counts [no]')
axz.yaxis.set_ticks_position('right')
axz.yaxis.set_label_position('right')
axz.set_ylabel('Depth [km]', rotation=270, labelpad=15)

if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'etype_map_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)

if isshow:
    plt.show()
