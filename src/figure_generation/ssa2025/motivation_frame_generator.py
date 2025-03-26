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
from obspy.clients.fdsn import Client

from map_util import *

### PATH DEFS SECTION ###
ROOT = Path(__file__).parent.parent.parent.parent
CPCSV = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
STYLE_FILE = Path(__file__).parent / 'etype_styles.csv'
SPATH = ROOT / 'results' / 'figures' / 'SSA2025' / 'pnsn_catalog_frames'
SFNAME = 'frame_{fnumber:05d}_{dpi}dpi.{fmt}'


### PLOTTING RENDERING CONTROL SECTION ###
DPI = 120
FMT = 'png'
istestframe = False
issave = True
isshow = False


### INVENTORY QUERY PARAMETER VALUES ###
# Reference Location - Mount Baker Summit
LAT, LON, ELE = 48.7745,-121.8172, 3286
# Network codes to include
NETS = 'UW,CN,TA'
# Radius limit for station query around Reference Location
INV_RADIUS = 100 # [km] Maximum station offset radius
netcolors = {'UW':'purple','CN':'firebrick','TA':'darkgoldenrod'}
stasize = 36 # Station marker size (used in plt.scatter)
smec = 'w'   # Station marker edge color (used in plt.scatter)
smlw = 0.5   # Station marker line width (used in plt.scatter)

### CATALOG SELECTION VALUES ###
CAT_RADIUS = 30. # [km] Maximum radius from 
ETYPES = ['eq','lf','su','px']

### TIME WINDOWING PARAMETER SECTION ###
T0 = pd.Timestamp('1979-01-01')
T1 = pd.Timestamp('2026-01-01')
WLEN = 365       # Window length for events to show
WUNIT = 'day'     # Window length unit
SLEN = 10       # Slide length between frames
SUNIT = 'day'     # Slide length unit
HALFLIFE = WLEN/5 # Marker opacity decay half-life
HLUNIT = 'day'    # Marker opacity decay half-life unit
slidercolor = 'cyan'

def title_fmt(t0, t1):
    # return f"{t0.strftime('%Y-%m-%d')} to {t1.strftime('%Y-%m-%d')}"
    # return f"{t0.strftime('%b %d %Y')} to {t1.strftime('%b %d %Y')}"
    return f"{t1.strftime('%b %d %Y')}"


### EVENT RENDERING HELPER FUNCTIONS ###
def alpha_decay(timestamp, t1, hl, hldelay=1):
    """Set time-dependent function for fading out marker opacity
    decay the "alpha" plotting parameter!
    """
    # Age since leading edge of window
    dt = t1 - timestamp
    # Number of halflives
    nhl = dt/hl
    if nhl < hldelay:
        return 1
    else:
        return np.exp(-1*(nhl - hldelay))

def magscale(mag, base=2.3, offset=9):
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


### PROCESSING SECTION ###

# Get maximum view radius
MAX_RADIUS = max([CAT_RADIUS, INV_RADIUS])

# Get station inventory
client = Client("IRIS")
inv = client.get_stations(
    longitude=LON,
    latitude=LAT,
    maxradius=INV_RADIUS/111.2,
    level='station',
    network=NETS)

# Load catalog profile
df = pd.read_csv(CPCSV, index_col='prefor_time', parse_dates=['prefor_time'])

# Subset catalog profile by min-max time
df = df[(df.index >= T0) & (df.index <= T1)]
# Subset to inside and outside CAT_RADIUS
df_in = df[df.offset_km <= CAT_RADIUS]
df_out = df[df.offset_km > CAT_RADIUS]

# Load event type style file
_styles = pd.read_csv(STYLE_FILE, index_col=[0]).T



### PLOTTING SECTION ###
# Initialize Figure
fig = plt.figure(figsize=(5.9,8.68))
gs = fig.add_gridspec(nrows=3, ncols=1, hspace=0)

# Get time iterators
t0 = T0 - pd.Timedelta(WLEN, unit=WUNIT)
t1 = T0

### SUBPLOT A ###
# Render basemap once
# Initial basemap
axm, map_attrib = mount_baker_basemap(fig=fig, sps=gs[:2], radius_km=MAX_RADIUS+5,
                                        open_street_map=False, zoom=10)
# Add coordinate frames
gl = axm.gridlines(draw_labels=True, zorder=1)
gl.bottom_labels = False
gl.right_labels=False

# Add distance reticles
add_rings(axm, 
          rads_km=[int(CAT_RADIUS), 50, int(INV_RADIUS)],
          rads_colors=['royalblue','dodgerblue', 'cornflowerblue'],
          label_pt=13,
          annotations=['Study Events','Nearby Events','Study\nStations'])
# Annotate distances


# Get axis limits
xlims = axm.get_xlim()
ylims = axm.get_ylim()

# Plot Mount Baker once on basemap
def plot_baker(geoaxis, zorder=2):
    handle = geoaxis.scatter(
        LON, LAT, marker='^',
        s=64, c='orange',
        edgecolors='k', zorder=zorder, 
        label='Mount Baker',
        transform=ccrs.PlateCarree())
    return handle

_ = plot_baker(axm)

# Spoof network markers
for net in NETS.split(','):
    icolor = netcolors[net]
    axm.scatter(0,0,c=icolor, s=stasize, zorder=3, marker='v',
                edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree(),
                label=net)


# Insert Legend
axm.legend(loc='upper left')

# Add attribution
axm.text(xlims[1] - 500, ylims[0] + 500, 
         '\n'.join(map_attrib), 
         fontsize=6, ha='right',va='bottom')



### SUBPLOT B ####
# Plot Time-Depth-Mag-Etype Figure Once
axf = fig.add_subplot(gs[-1])
axt = axf.twinx()

for _etype in ETYPES:
    _df = df_in[df_in.etype == _etype]
    axt.scatter(_df.index, _df.depth*1e-3, c=_styles[_etype].color, s=magscale(_df.mag),
                marker=_styles[_etype].marker, zorder=_styles[_etype].zorder, alpha=0.5,
                label=_etype.upper())
axt.legend(ncols=4, loc='upper center')
# Format axes
ylims = list(axt.get_ylim())
ylims[0] -= 6
axt.set_ylim([ylims[1], ylims[0]])
axt.set_xlim([T0, T1])
axt.set_ylabel('Depth [km]', rotation=270, labelpad=15)
axt.set_xlabel('Time')

axf.set_ylim([-1, 100])
if WLEN > 1:
    modifier = 's'
else:
    modifier = ''
axf.set_ylabel(f'Event Frequency\n[Ct./{WLEN} {WUNIT}{modifier}]', labelpad=-5)


### THEN DO SUMMARY FIGURE 0th FRAME LAST FRAME
handles = {}
ilon, ilat, inet, ista, icolor= [], [], [], [], []
for _n in inv.networks:
    for _s in _n.stations:
        ilon.append(_s.longitude)
        ilat.append(_s.latitude)
        ista.append(_s.code)
        inet.append(_n.code)
        icolor.append(netcolors[_n.code])
    
# Plot all nearby events
neh = axm.scatter(df_out.lon, df_out.lat, c='k', s=magscale(df_out.mag),
            zorder=2, marker='+', alpha=0.25, linewidths=0.5, transform=ccrs.PlateCarree())
handles.update({'nearby': neh})
# Plot all stations
ssh = axm.scatter(ilon, ilat, c=icolor, s=stasize, zorder=3, marker='v',
            edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
handles.update({'stations': ssh})
# Plot all study events by etype
for _etype in ETYPES:
    _df_in = df_in[df_in.etype == _etype]
    seh = axm.scatter(_df_in.lon, _df_in.lat, c=_styles[_etype].color, s=magscale(_df_in.mag),
                zorder=_styles[_etype].zorder, marker=_styles[_etype].marker,alpha=0.5,
                transform=ccrs.PlateCarree())
    handles.update({_etype: seh})

# Title
axm.set_title(f"Mount Baker cataloged seismicity\n{df_in.index.min().strftime('%b %Y')} to {df_in.index.max().strftime('%b %Y')}")

# Temporary Mount Baker
bh = plot_baker(axm, zorder=100)
handles.update({'tmpmb': bh})


# plt.show()
# breakpoint()

if issave:
    try:
        os.makedirs(str(SPATH))
        print(f'Created new save directory: {SPATH}')

    except:
        pass
    fig.savefig(str(SPATH/SFNAME).format(fnumber=0, dpi=DPI, fmt=FMT), dpi=DPI, format=FMT)

# Remove temporary items
for _v in handles.values():
    _v.remove()


# Calcuate total number of frames
# Do not exceed the index of the last event (by too much)
T1 = min([T1, df_in.index.max()])

total_frames = int((T1 - t1)/pd.Timedelta(SLEN, unit=SUNIT)) + 1
_fno = 1

ftimes = []
fvalues = defaultdict(list)

while t1 <= T1:
    # Capture leading edge time
    ftimes.append(t1)
    
    ### DATA WINDOWING ###
    # Subset events
    _df_in = df_in[(df_in.index > t0) & (df_in.index <= t1)]
    _df_out = df_out[(df_out.index > t0) & (df_out.index <= t1)]
    # Subset inventory
    _inv = inv.select(time=UTCDateTime(t1.timestamp()))

    # Get station locations & names
    ilon, ilat, inet, ista, icolor= [], [], [], [], []
    for _n in _inv.networks:
        for _s in _n.stations:
            ilon.append(_s.longitude)
            ilat.append(_s.latitude)
            ista.append(_s.code)
            inet.append(_n.code)
            icolor.append(netcolors[_n.code])
    
    handles = {}

    # # Get event delta times
    # _df = _df.assign(_dt=[(idx - t0)/pd.Timedelta(WLEN, unit=WUNIT) for idx in _df.index])
    # # Get event marker alphas
    # _df = _df.assign(_alpha=[np.exp(-1*(idx - t0)/pd.Timedelta(HALFLIFE, unit=HLUNIT)) for idx in _df.index])
    _df_in = _df_in.assign(_alpha=[alpha_decay(idx, t1, pd.Timedelta(HALFLIFE, unit=HLUNIT), hldelay=1) for idx in _df_in.index])
    _df_out = _df_out.assign(_alpha=[alpha_decay(idx, t1, pd.Timedelta(HALFLIFE, unit=HLUNIT), hldelay=1) for idx in _df_out.index])
   
   
    ### WINDOW PLOTTING ###
    # Plot Stations in basemap
    sh = axm.scatter(ilon, ilat, c=icolor, s=stasize, zorder=3, marker='v',
                     edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
    handles.update({'stations': sh})
    # Plot Events in basemap

    # Iterate across etypes
    for _etype in ETYPES:
        __df_in = _df_in[_df_in.etype == _etype]
        fvalues[_etype].append(len(__df_in))
        if len(__df_in) > 0:
            # Plot Events
            eh = axm.scatter(
                __df_in.lon, __df_in.lat, c=_styles[_etype].color, s=magscale(__df_in.mag),
                zorder=_styles[_etype].zorder, marker=_styles[_etype].marker, alpha=__df_in._alpha,
                transform=ccrs.PlateCarree())
            handles.update({_etype: eh})
    if len(_df_out) > 0:
        oh = axm.scatter(_df_out.lon, _df_out.lat, c='k', s=magscale(_df_out.mag),
                        zorder=2, marker='+', linewidths=0.5, alpha=_df_out._alpha,
                        transform=ccrs.PlateCarree())
        handles.update({'outside': oh})
    

    axm.set_title(title_fmt(t0, t1))

    # Plot sampling window in time-depth-mag-etype plot
    wh = axt.fill_between([t0, t1], [ylims[0]]*2, [ylims[1]]*2, color=slidercolor, alpha=0.25)
    # Plot rolling-window count of events
    handles.update({'window': wh})
    for _k, _v in fvalues.items():
        efh = axf.plot(ftimes, _v, color=_styles[_k].color)
        handles.update({f'{_k}_freq': list(efh)})
        # Plot rolling-window start points
        eflh = axf.scatter(ftimes[-1], _v[-1], color=_styles[_k].color,
                           marker=_styles[_k].marker, s=16)
        handles.update({f'{_k}_flead': eflh})
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
    for _v in handles.values():
        try:
            _v.remove()
        except:
            _v[0].remove()

    print(f'Frame {_fno:d}/{total_frames:d} ({t0} - {t1}) Complete. Advancing Indices')
    # Increment up frame number
    _fno += 1
    # Increment up time window
    t0 += pd.Timedelta(SLEN, unit=SUNIT)
    t1 += pd.Timedelta(SLEN, unit=SUNIT)

    

### THEN DO SUMMARY FIGURE (AGAIN) AS LAST FRAME
### THEN DO SUMMARY FIGURE 0th FRAME LAST FRAME
handles = {}
ilon, ilat, inet, ista, icolor= [], [], [], [], []
for _n in inv.networks:
    for _s in _n.stations:
        ilon.append(_s.longitude)
        ilat.append(_s.latitude)
        ista.append(_s.code)
        inet.append(_n.code)
        icolor.append(netcolors[_n.code])
    
# Plot all nearby events
neh = axm.scatter(df_out.lon, df_out.lat, c='k', s=magscale(df_out.mag),
            zorder=2, marker='+', alpha=0.25, linewidths=0.5, transform=ccrs.PlateCarree())
handles.update({'nearby': neh})
# Plot all stations
ssh = axm.scatter(ilon, ilat, c=icolor, s=stasize, zorder=3, marker='v',
            edgecolors=smec, linewidths=smlw, transform=ccrs.PlateCarree())
handles.update({'stations': ssh})
# Plot all study events by etype
for _etype in ETYPES:
    _df_in = df_in[df_in.etype == _etype]
    seh = axm.scatter(_df_in.lon, _df_in.lat, c=_styles[_etype].color, s=magscale(_df_in.mag),
                zorder=_styles[_etype].zorder, marker=_styles[_etype].marker,alpha=0.5,
                transform=ccrs.PlateCarree())
    handles.update({_etype: seh})

# Title
axm.set_title(f"Mount Baker cataloged seismicity\n{df_in.index.min().strftime('%b %Y')} to {df_in.index.max().strftime('%b %Y')}")

# Temporary Mount Baker
bh = plot_baker(axm, zorder=100)
handles.update({'tmpmb': bh})


# Plot rolling-window count of events
handles.update({'window': wh})
for _k, _v in fvalues.items():
    efh = axf.plot(ftimes, _v, color=_styles[_k].color)
    handles.update({f'{_k}_freq': list(efh)})

if issave:
    fig.savefig(str(SPATH/SFNAME).format(fnumber=_fno, dpi=DPI, fmt=FMT), format=FMT, dpi=DPI)

# Remove temporary items
for _v in handles.values():
    try:
        _v.remove()
    except:
        _v[0].remove()




