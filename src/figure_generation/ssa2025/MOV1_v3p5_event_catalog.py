"""

TODO: Change depth scale to <=0 to >= 30
TODO: Add in US-Canada border in main plot
TODO: Add in coastline with main plot
"""

import os
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt
import cartopy.crs as ccrs


import numpy as np
import pandas as pd

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client

import map_util as mutil

### PATH DEFS SECTION ###
ROOT = Path(__file__).parent.parent.parent.parent
CPCSV = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
# STYLE_FILE = Path(__file__).parent / 'etype_styles.csv'
SPATH = ROOT / 'results' / 'figures' / 'SSA2025' / 'catalog_frames'
SFNAME = 'frame_{fnumber:05d}_{dpi}dpi.{fmt}'

### PLOTTING RENDERING CONTROL SECTION ###
DPI = 120
FMT = 'png'
istestframe = False
issave = True
isshow = False
PP = mutil.pnsn_pallet()

# BASEMAP PARAMETERS
LAT = mutil.BAKER_LAT
LON = mutil.BAKER_LON
RADIUS_KM = 31.
BASEMAP_ZOOM = 10
BASEMAP_CMAP ='Greys_r'
BASEMAP_ALPHA=0.1

# TIME STEPPING PARAMETERS
T0 = pd.Timestamp('1980-01-01') # Start Datetime (first frame newest event)
T1 = pd.Timestamp('2025-03-01') # End Datetime (last frame newest event)
WLEN = 180       # Window length for events to show
WUNIT = 'day'     # Window length unit
SLEN = 10       # Slide length between frames
SUNIT = 'day'     # Slide length unit
HALFLIFE = pd.Timedelta(WLEN/5,unit=WUNIT)

# EVENT SYTLE FORMATTING
marker_map = {'su':'d', 'px':'*','lf':'s','eq':'o'}
label_map = {'eq':'EQ: Earthquake',
             'lf':'LF: Low Frequency Event',
             'su':'SU: Surface Event',
             'px':'PX: Probable Blast'}
mec_map = {'su':'b', 'px':'m', 'eq':'k','lf':'r'}
lw_map = {'su': 1, 'px': 1, 'eq': 1, 'lf': 1}
depth_cmap = 'Greens'
base_alpha=0.8
dbins = np.arange(0, 35, 5)

title_str = "Seismicity and Seismic Stations Near Mount Baker\n{datestr}"


def alpha_decay(timestamp, reference_time, halflife, hldelay=1):
    """Set time-dependent function for fading out marker opacity
    decay the "alpha" plotting parameter!
    """
    # Age of timestamp relative to reference_time
    dt = reference_time - timestamp
    # Number of halflives
    nhl = dt/halflife
    if nhl < hldelay:
        return 1
    else:
        return np.exp(-1*(nhl - hldelay))

def magscale(mag, base=4, offset=9):
    """Magnitude scaling function

    scale = offset + base**mag

    :param mag: magnitude measurement
    :type mag: scalar or array-like
    :param base: exponential base to use, defaults to 2.3
    :type base: float, optional
    :param offset: prefactor, defaults to 9
    :type offset: int, optional
    :return: scales
    :rtype: scalar or array-like (same as **mag**)
    """    
    return base**(mag) + offset

def depth_binner(depths, bins=dbins):
    return np.digitize(depths, bins=bins)


def plot_events(geoaxis, df):
    handles = []
    if 'alpha' not in df.columns:
        df = df.assign(alpha=[base_alpha for _ in range(len(df))])
    for _k,_v in marker_map.items():
        _df = df[df.etype==_k]
        if len(_df) == 0:
            continue

        _df = _df.sort_values(by='depth', ascending=True)
        _h = geoaxis.scatter(
            _df.lon, _df.lat, 
            c=depth_binner(_df.depth*1e-3),
            marker=_v, 
            edgecolors=mec_map[_k],
            linewidths=lw_map[_k],
            s=magscale(_df.mag), 
            alpha=_df.alpha,
            vmin=1, vmax=len(dbins),
            transform=mutil.ccrs.PlateCarree(),
            cmap=depth_cmap)
        handles.append(_h)
    return handles

def plot_stations(geoaxis, inv, reference_time=None, label_netsta=False, **options):
    if isinstance(reference_time, pd.Timestamp):
        rtime = UTCDateTime(reference_time.timestamp())
        _inv = inv.select(time=rtime)
    else:
        _inv = inv
    holder = []
    for net in _inv.networks:
        for sta in net.stations:
            line = [net.code, sta.code, sta.longitude, sta.latitude]
            holder.append(line)
    df_inv = pd.DataFrame(holder, columns=['net','sta','lon','lat'])
    hdl = geoaxis.scatter(df_inv.lon, df_inv.lat, s=6**2, c='k',
        marker='v', edgecolors='k', linewidths=0.5,
         transform=ccrs.PlateCarree(), **options)
    if label_netsta:
        for _, row in df_inv.iterrows():
            axm.text(row.lon+0.005, row.lat+0.005, f'{row.net}.{row.sta}', ha='left', va='bottom', transform=ccrs.PlateCarree())
    return hdl
    

# Load Event Metadata
df = pd.read_csv(CPCSV, index_col='prefor_time', parse_dates=['prefor_time'])
# Apply time filtering
df = df[(df.index >= T0) & (df.index <= T1)]

# Get station inventory
inv = Client('IRIS').get_stations(longitude=LON, latitude=LAT, maxradius=1,
                                  network='UW,CC,CN,TA,GS,NP', level='station')

#### PLOT INITIALIZATION SECTION ####

# Initialize Figure and GridSpec
fig = plt.figure(figsize=(10,7))
gs = fig.add_gridspec(nrows=14, ncols=21)

# Generate Basemap
axm, map_attr = mutil.mount_baker_basemap(
    fig=fig, sps=gs[:, :18], radius_km=RADIUS_KM, zoom=BASEMAP_ZOOM,
    latnudge=0, lonnudge=0,
    open_street_map=False, aws_add_image_kwargs={'cmap':BASEMAP_CMAP,
                                                 'alpha': BASEMAP_ALPHA})

# Add political boundaries
axm.add_feature(mutil.cfeature.BORDERS, linestyle='-')
# Label State/Province
axm.text(-122.21, 49.005, 'BC', fontsize=14, transform=ccrs.PlateCarree())
axm.text(-122.21, 48.97, 'WA', fontsize=14, transform=ccrs.PlateCarree())
# Add coastline
axm.add_feature(mutil.cfeature.COASTLINE)

# Add distance Rings
mutil.add_rings(axm,
                rads_km=[10,20,30],
                rads_colors=['k']*3,
                include_units=[True, True, True],
                label_pt=18,ha='left',va='top')

# Add Lat/Lon annotations
# axm.set_xticks([-122.2 + _e*0.2 for _e in range(5)])
gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
                   xlocs=[-122.2, -122, -121.8, -121.6, -121.4, -121.2],
                   ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
                   alpha=0)

## GENERATE LEGEND ELEMENTS
legend_handles = []
# Add Mount Baker
mbh = mutil.plot_baker(axm, zorder=100, color='orange')
legend_handles.append(mbh)

# Add Seismic Station
sth = axm.scatter(
    0,0,c='k', s=6**2, 
    marker='v', edgecolors='k', linewidths=0.5,
    zorder=90, label='Seismic Station')
legend_handles.append(sth)

# Add Event Types (and get colorbar handle as a kicker)
for _k, _v in label_map.items():
    _h = axm.scatter(
        0,0,
        c=2,
        s=36,
        marker=marker_map[_k],
        edgecolors=mec_map[_k],
        linewidths=lw_map[_k],
        alpha=base_alpha, vmin=1, vmax=len(dbins),
        cmap=depth_cmap,
        label=_v)
    legend_handles.append(_h)

# Generate Legend
plt.legend(handles=legend_handles, bbox_to_anchor=(1,0.98))

# Generate Depth Scale
cbh = plt.colorbar(_h,
                   cax=fig.add_subplot(gs[4:8, 18]),
                   orientation='vertical',
                   label='Event Depth [km]')
# cbh.set_ticks(np.arange(1, len(dbins) + 1)[::2],
#               labels=[f'{int(_b)}' for _b in dbins][::2])
cbh.set_ticks([1, 3, 5, 7], labels=['$\leq$0','10','20','$\geq$30'])
cbh.ax.invert_yaxis()

# Generate Magnitude Scale
magax = fig.add_subplot(gs[9:13, 17:])
xpos = 0
for _m in range(0, 6):
    for _e, _k in enumerate(label_map.keys()):
        if _m <= np.ceil(df[df.etype==_k].mag.max()):
            _s= magscale(_m)

            magax.scatter(_m, 5 - _e, s=_s, c=2, alpha=base_alpha,
                        vmin=1,vmax=len(dbins), cmap=depth_cmap,
                        linewidths=lw_map[_k],
                        marker=marker_map[_k],
                        edgecolors=mec_map[_k])
magax.set_xlim([-0.5, 6.25])
magax.set_ylim([1,6])
magax.yaxis.set_ticks_position('right')
magax.yaxis.set_label_position('right')
magax.set_xlabel('Magnitude')
magax.set_ylabel('Event Type', rotation=270, labelpad=15)
magax.yaxis.set_ticks([5, 4, 3, 2], labels=['EQ','LF','SU','PX'])
magax.xaxis.set_ticks([0,1,2,3,4,5],labels=['$\leq$0','1','2','3','4','5'])
# magax.text(4, 3, 'None\nObserved', ha='center', va='center')

# Generate Reference Map
axms = mutil.add_inset_map(
    fig=fig,pad_deg=7.5, extent=[0.648, 0.778, 0.1, 0.10])
mutil.plot_baker(axms, zorder=20, color='orange')

#### START OF FRAME GENERATION ####
try:
    os.makedirs(str(SPATH),exist_ok=False)
except:
    pass

_frame_number=0

#### FIRST FRAME WITH DISPOSABLE CONTENTS (OVERVIEW) ####
# Plot Events
tmp_handles = plot_events(axm, df)
# Plot Stations
tmp_sta_hdl = plot_stations(axm, inv, reference_time=None)

tstr = f"{T0.strftime('%b %Y')} to {T1.strftime('%b %Y')}"

axm.set_title(title_str.format(datestr=tstr),
              fontsize=14)


# Save Summary Frame
fig.savefig(str(SPATH/SFNAME).format(fnumber=_frame_number, dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)
# Depopulate Stations and Events
tmp_handles.append(tmp_sta_hdl)
for _h in tmp_handles:
    _h.remove()
_frame_number += 1

# Initialize moving window trackers
t0 = T0 - pd.Timedelta(WLEN, unit=WUNIT)
t1 = T0

total_frames = int((T1 - t1)/pd.Timedelta(SLEN, unit=SUNIT)) + 1


# Iterate across windows
while t1 < T1:
    # Update Title
    axm.set_title(title_str.format(datestr=f"{t1.strftime('%b %Y')}"))

    # Subset Events
    _df = df[(df.index > t0) & (df.index <= t1)].copy()

    if len(_df) > 0:
        # Attach alphas
        _df = _df.assign(alpha=[alpha_decay(idx, reference_time=t1, halflife=HALFLIFE) for idx in _df.index])
        # Plot events
        tmp_handles = plot_events(axm, _df)
    else:
        tmp_handles = []
    # Plot stations
    tmp_sta_hdl = plot_stations(axm, inv, reference_time=t1, label_netsta=False)
    # Gather handles
    tmp_handles.append(tmp_sta_hdl)


    # Save frame to file
    fig.savefig(str(SPATH/SFNAME).format(fnumber=_frame_number, dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)
    print(f'Frame {_frame_number:d}/{total_frames:d} ({t0} - {t1}) Complete. Advancing Indices')

    # Advance time indices
    t0 += pd.Timedelta(SLEN, unit=SUNIT)
    t1 += pd.Timedelta(SLEN, unit=SUNIT)
    _frame_number += 1
    # Depopulate Stations and Events
    for _h in tmp_handles:
        _h.remove()


