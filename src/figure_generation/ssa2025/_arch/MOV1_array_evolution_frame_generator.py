"""
This script generates a sequence of frame-number annotated frames
that can be assembled into a movie with QuickTime that shows the evolution
of permanent stations and seismicity around Mount Baker from an ObsPy Inventory
and the catalog profile output from Phase1 Step1

TODO: Terminate when t1 exceeds the timestamp of the last event

"""
import os
from pathlib import Path
from collections import defaultdict

import matplotlib.pyplot as plt

import pandas as pd

from obspy import UTCDateTime
from obspy.geodetics import gps2dist_azimuth
from obspy.clients.fdsn import Client

from map_util import *

### PATH DEFS SECTION ###
ROOT = Path(__file__).parent.parent.parent.parent
# CPCSV = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
# STYLE_FILE = Path(__file__).parent / 'etype_styles.csv'
SPATH = ROOT / 'results' / 'figures' / 'SSA2025' / 'network_evolution_frames'
SFNAME = 'frame_{fnumber:05d}_{dpi}dpi.{fmt}'


### PLOTTING RENDERING CONTROL SECTION ###
DPI = 120
FMT = 'png'
istestframe = False
issave = True
isshow = False
PP = pnsn_pallet()

### INVENTORY QUERY PARAMETER VALUES ###
# Reference Location - Mount Baker Summit
LAT, LON, ELE = 48.7745,-121.8172, 3286
# Network codes to include
NETS = 'UW,CC,CN,TA,GS'
# Radius limit for station query around Reference Location
INV_RADIUS = 100 # [km] Maximum station offset radius
netcolors = {'UW':PP['navy'],
             'CC':'royalblue',
             'CN':'firebrick',
             'GS':'olive',
             'TA':'darkgoldenrod'}
            #  'NP':'yellow'}
NETS = ','.join(list(netcolors.keys()))
stasize = 36 # Station marker size (used in plt.scatter)
smec = 'w'   # Station marker edge color (used in plt.scatter)
smlw = 0.5   # Station marker line width (used in plt.scatter)


### TIME WINDOWING PARAMETER SECTION ###
T0 = pd.Timestamp('1971-01-01')
T1 = pd.Timestamp('2025-03-01')
WLEN = 30       # Window length for events to show
WUNIT = 'day'     # Window length unit
SLEN = 30       # Slide length between frames
SUNIT = 'day'     # Slide length unit

slidercolor = 'cyan'
def title_fmt(t0, t1):
    return f"Active Stations\n{t1.strftime('%Y')}\n{t1.strftime('%B')}"

events = [[pd.Timestamp('1975-01-01'), pd.Timestamp('1980-01-01'), 'Unrest at Mount Baker Begins'],
          [pd.Timestamp('1980-01-01'), pd.Timestamp('1995-09-01'), 'Progressive UW and Strong Motion Expansion'],
          [pd.Timestamp('1995-09-01'), pd.Timestamp('1997-04-01'), 'NRCan Data Become Accessible'],
          [pd.Timestamp('2009-09-01'), pd.Timestamp('2010-09-01'), 'UW.SHUK Installed'],
          [pd.Timestamp('2013-01-01'), pd.Timestamp('2025-03-01'), 'ShakeAlert Build-Out']]
events = pd.DataFrame(events, columns=['start','end','text'])


# Get station inventory
client = Client("IRIS")
inv = client.get_stations(
    longitude=LON,
    latitude=LAT,
    maxradius=INV_RADIUS/111.2,
    level='channel',
    network=NETS)

# Form dataframe
holder = []
for _n in inv.networks:
    for _s in _n.stations:
        if _s.end_date is None:
            endtime = UTCDateTime(T1.timestamp())

        else:
            endtime = _s.end_date
        dt = (endtime - _s.start_date)/(3600*24*365.24)

        if _n.code in netcolors.keys():
            color = netcolors[_n.code]
        else:
            color = 'k'
        
        distm, seaz, esaz = gps2dist_azimuth(_s.latitude, _s.longitude, LAT, LON)
        line = [_n.code, _s.code, _s.longitude, _s.latitude, distm*1e-3, seaz, esaz,
                pd.Timestamp(_s.start_date.timestamp, unit='s'), 
                pd.Timestamp(endtime.timestamp, unit='s'),
                dt, color]
        holder.append(line)

df_inv = pd.DataFrame(holder, columns=['net','sta','lon','lat','distkm','seaz','esaz','start','end','dt','color'])



### PLOTTING SECTION ###
# Initialize Figure
fig = plt.figure(figsize=(7,7))
gs = fig.add_gridspec(nrows=1, ncols=1, hspace=0)

# Get time iterators
t0 = T0 - pd.Timedelta(WLEN, unit=WUNIT)
t1 = T0

### SUBPLOT A ###
# Render basemap once
# Initial basemap
axm, map_attrib = mount_baker_basemap(fig=fig, sps=gs[0], radius_km=INV_RADIUS + 5,
                                        open_street_map=False, zoom=10)

axm.add_feature(cfeature.COASTLINE, color='dodgerblue', linewidth=0.5, alpha=0.5)
axm.add_feature(cfeature.BORDERS)

# Add coordinate frames
gl = axm.gridlines(draw_labels=True, zorder=1)
gl.bottom_labels = False
gl.right_labels=False

# Add distance reticles
add_rings(axm, 
          rads_km=[25, 50, 75, 100],
          rads_colors=['k','k','k','k'],
          label_pt=-13,
          annotations=[None, None, None, None],
          va='top')

lhand = []
# Plot permanent Mount Baker marker
pmbh = plot_baker(axm, zorder=100)
lhand.append(pmbh)
# Spoof network markers

for net in NETS.split(','):
    icolor = netcolors[net]
    nmh = axm.scatter(0,0,c=icolor, s=stasize, zorder=3, marker='v',
                edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree(),
                label=net)
    lhand.append(nmh)


# Insert Legend
legend_labels = ['Mount Baker']
legend_labels += [f'{net}: {len(df_inv[df_inv.net==net]):2d}' for net in NETS.split(',')]
axm.legend(lhand, legend_labels, loc='upper right')

# Add attribution
xlims = axm.get_xlim()
ylims = axm.get_ylim()
axm.text(xlims[1] - 500, ylims[0] + 500, 
         '\n'.join(map_attrib), 
         fontsize=6, ha='right',va='bottom')


### THEN DO SUMMARY FIGURE 0th FRAME
handles = []
# Plot all station-to-mountain lines
for _, row in df_inv.iterrows():
    sshl = axm.plot([row.lon, LON],[row.lat, LAT],'-',
                    color=row.color, zorder=2, 
                    lw=1, alpha=(100-row.distkm)/100.,
                    transform=ccrs.PlateCarree())
    handles.append(sshl)
# Plot all stations
ssh = axm.scatter(df_inv.lon, df_inv.lat, c=df_inv.color, 
                  s=10 + 2*df_inv.dt, zorder=3, marker='v',
                  edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
handles.append(ssh)

# Add Title
axm.set_title(f"Mount Baker Area Seismic Stations\nNetworks: {', '.join(NETS.split(','))} \n{T0.strftime('%b %Y')} to {T1.strftime('%b %Y')}")


if issave:
    try:
        os.makedirs(str(SPATH))
        print(f'Created new save directory: {SPATH}')

    except:
        pass
    fig.savefig(str(SPATH/SFNAME).format(fnumber=0, dpi=DPI, fmt=FMT),
                dpi=DPI, format=FMT)

# Remove temporary items
for _v in handles:
    try:
        _v.remove()
    except:
        list(_v)[0].remove()


# Calcuate total number of frames
# Do not exceed the index of the last event (by too much)

total_frames = int((T1 - t1)/pd.Timedelta(SLEN, unit=SUNIT)) + 1
_fno = 1

ftimes = []
fvalues = defaultdict(list)

while t1 <= T1:
    # Capture leading edge time
    ftimes.append(t1)  
    
    ### DATA WINDOWING ###
    # Subset inventory df
    _df_inv = df_inv[(df_inv.start <= t1) & (df_inv.end >= t0)]

    handles = []
   
    ### WINDOW PLOTTING ###
    handles = []
    # Plot all station-to-mountain lines
    for _, row in _df_inv.iterrows():
        sshl = axm.plot([row.lon, LON],[row.lat, LAT],'-',
                        color=row.color, zorder=2, 
                        lw=1, alpha=(100-row.distkm)/100.,
                        transform=ccrs.PlateCarree())
        handles.append(sshl)
    # Plot all stations
    ssh = axm.scatter(_df_inv.lon, _df_inv.lat, c=_df_inv.color, 
                    s=[10 + 2*(t1 - start).days/365.24 for start in _df_inv.start], zorder=3, marker='v',
                    edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
    handles.append(ssh)
    
    # Set new title
    axm.set_title(title_fmt(t0, t1))
    
    # Update legend with counts
    legend_labels = ['Mount Baker']
    for net in NETS.split(','):
        legend_labels.append(f'{net}: {len(_df_inv[_df_inv.net==net]):2d}')
    
    axm.legend(lhand, legend_labels, loc='upper right')

    ### CLEANUP SECTION ###

    if istestframe:
        if _fno == 155:
            plt.show()  
            breakpoint()

    # Save figure
    if not istestframe:
        if issave:

            fig.savefig(str(SPATH/SFNAME).format(fnumber=_fno, dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)
    
    # Clear Handles
    for _v in handles:
        try:
            _v.remove()
        except:
            list(_v)[0].remove()

    print(f'Frame {_fno:d}/{total_frames:d} ({t0} - {t1}) Complete. Advancing Indices')
    
    # Increment up frame number
    _fno += 1
    # Increment up time window
    t0 += pd.Timedelta(SLEN, unit=SUNIT)
    t1 += pd.Timedelta(SLEN, unit=SUNIT)

    

### THEN DO SUMMARY FIGURE (AGAIN) AS LAST FRAME
handles = []
# Plot all station-to-mountain lines
for _, row in df_inv.iterrows():
    sshl = axm.plot([row.lon, LON],[row.lat, LAT],'-',
                    color=row.color, zorder=2, 
                    lw=1, alpha=(100-row.distkm)/100.,
                    transform=ccrs.PlateCarree())
    handles.append(sshl)
# Plot all stations
ssh = axm.scatter(df_inv.lon, df_inv.lat, c=df_inv.color, 
                  s=10 + 2*df_inv.dt, zorder=3, marker='v',
                  edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
handles.append(ssh)

# Add Title
axm.set_title(f"Mount Baker Area Seismic Stations\n{', '.join(NETS.split(','))} Networks\n{T0.strftime('%b %Y')} to {T1.strftime('%b %Y')}")

# Insert Legend
legend_labels = ['Mount Baker']
legend_labels += [f'{net}: {len(df_inv[df_inv.net==net]):2d}' for net in NETS.split(',')]
axm.legend(lhand, legend_labels, loc='upper right')

if issave:
    try:
        os.makedirs(str(SPATH))
        print(f'Created new save directory: {SPATH}')

    except:
        pass
    fig.savefig(str(SPATH/SFNAME).format(fnumber=_fno, dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)

# Remove temporary items
for _v in handles:
    try:
        _v.remove()
    except:
        list(_v)[0].remove()




