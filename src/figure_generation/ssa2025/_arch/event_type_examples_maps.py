from pathlib import Path

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from obspy.imaging.spectrogram import spectrogram
from obspy.clients.fdsn import Client

from obsplus import EventBank

import cartopy.crs as ccrs

import map_util as mutil

# Repository Root Absolute Path
ROOT = Path(__file__).parent.parent.parent.parent
# Catalog Event Bank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Catalog Profile CSV
CPCSV = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
# Save Path
SPATH = ROOT / 'results' / 'figures' / 'SSA2025'
# Rendering Settings
DPI = 120
FMT = 'png'
# Output control
issave = True
isshow = True


### SPECTROGRAM INPUTS ###
pref_netsta = ['UW','MBW']

event_id_fstr = 'quakeml:uw.anss.org/Event/UW/{evid}'

pref_evid = {'eq': '10243683',
             'lf': '61822721',
             'px': '10800098',
             'su': '61802931'}

### MAP PARAMETERS
# BASEMAP PARAMETERS
LAT = mutil.BAKER_LAT
LON = mutil.BAKER_LON
RADIUS_KM = 31.
BASEMAP_ZOOM = 10
BASEMAP_CMAP ='Greys_r'
BASEMAP_ALPHA=0.1

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

# Waveform Query Arguments
prepick=10
duration = 60

# Spectrogram kwargs
spect_kw = {'per_lap': 0.95,
            'wlen': 2.56,
            'log': True,
            'dbscale': False,
            'cmap': 'Spectral_r'}


# CONNECT TO CLIENTS
ebank = EventBank(EBBP)
client = Client('IRIS')

# Compose bulk waveform request
picks = {}
for _etype, _evid in pref_evid.items():
    cat = ebank.get_events(event_id=event_id_fstr.format(evid=_evid))
    for _p in cat[0].picks:
        if _p.waveform_id.id.split('.')[:2] == pref_netsta:
            if _etype not in picks.keys():
                picks[_etype] = _p
            elif _p.time < picks[_etype].time:
                picks[_etype] = _p

bulk = {_e: (_p.waveform_id.id.split('.') + [_p.time - prepick, _p.time + duration]) for _e, _p in picks.items()}

# Get waveforms
st = client.get_waveforms_bulk(list(bulk.values()))

# Associate traces to etypes and assess completeness
traces = {}
prepads = {}
durations = {}
for _etype, pick in picks.items():
    for tr in st:
        if tr.stats.starttime < pick.time < tr.stats.endtime:
            traces[_etype] = tr
            prepads[_etype] = picks[_etype].time - tr.stats.starttime
            durations[_etype] = tr.stats.endtime - picks[_etype].time


### PLOTTING SECTION ###
            
### SPECTROGRAMS ###
figs = {}
axes = {}
caxes = {}
for _etype, _tr in traces.items():
    __tr = _tr.copy().detrend('linear').normalize()
    evid = pref_evid[_etype]
    pick = picks[_etype]
    fig = plt.figure(figsize=(4,5.25))
    gs = fig.add_gridspec(ncols=1, nrows=2, hspace=0, wspace=0)
    axs = [fig.add_subplot(gs[0])]
    axs.append(fig.add_subplot(gs[1], sharex=axs[0]))

    figs[_etype] = fig
    axes[_etype] = axs

    axs[0].plot(__tr.times(reftime=__tr.stats.starttime),
                __tr.data,
                mec_map[_etype],
                lw=0.5)
    shdl = spectrogram(__tr.data, __tr.stats.sampling_rate, axes=axs[1], **spect_kw)
    caxes[_etype] = shdl

    # axs[0].xaxis.set_visible(False)
    axs[1].set_xticks([_e + prepads[_etype] for _e in range(-10, 70, 10)], labels=[f'{_e:d}' for _e in range(-10, 70, 10)])


    xmin = -5 + prepads[_etype] + 2.56
    xmax = 30 + prepads[_etype] + 2.56
    axs[0].set_xlim([xmin,xmax])
    axs[0].set_yticks([-1, 0, 1], labels=[-1, 0, 1])
    axs[0].set_ylim([-1.2, 1.2])
    axs[1].set_ylim([0.5,30])

    # Add Labels
    axs[0].set_title(f'uw{evid} | ETYPE: {_etype.upper()}\n{pick.time}')
    axs[0].set_ylabel('Normalized Amplitude [ - ]')
    axs[0].grid(alpha=0.5, which='major')
    axs[0].text(20, __tr.data.min(), pick.waveform_id.id, fontsize=12, ha='left', va='bottom')

    axs[1].set_yticks([0.5,1,2,3,5,10,20,30], [0.5,1,2,3,5,10,20,30])
    axs[1].set_ylabel('Frequency [Hz]',labelpad=0)
    axs[1].set_xlabel(f'Elapsed Time Since P-Wave Arrival [sec]')
    axs[1].grid(alpha=0.5)

            # req = _p.waveform_id.id.split('.') + [_p.time - prepad, _p.time + duration]
            # req = dict(zip(['network','station','location','channel','starttime','endtime'], req))
            # traces.update({_etype: client.get_waveforms(**req)})

if issave:
    for _etype, _fig in figs.items():
        _fig.savefig(str(SPATH/f'wave_spect_UW_MBW_{_etype}_uw{pref_evid[_etype]}_narrow_{DPI}dpi.{FMT}'), dpi=DPI, format=FMT)    



### TYPE MAPS
df_cat = pd.read_csv(CPCSV, parse_dates=['prefor_time'])
etype_dfs = {_etype: df_cat[df_cat.etype==_etype] for _etype in ['eq','lf','px','su']}
etype_dfs['all'] = df_cat




fig = plt.figure(figsize=(10, 7))
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
mbh = mutil.plot_baker(axm, color='orange')
legend_handles.append(mbh)

# Add Seismic Stations
inv = Client('IRIS').get_stations(longitude=LON, latitude=LAT, maxradius=1,level='station')
lllist = []
for net in inv.networks:
    for sta in net.stations:
        lllist.append([sta.longitude, sta.latitude, net.code, sta.code])
df_inv = pd.DataFrame(lllist, columns=['lon','lat','net','sta'])
sth = axm.scatter(
    df_inv.lon,df_inv.lat,c='k', s=6**2, 
    marker='v', edgecolors='k', linewidths=0.5,
    zorder=90, label='Seismic Station',
    transform=ccrs.PlateCarree())
legend_handles.append(sth)

# Highlight MBW
_df_mbw = df_inv[df_inv.sta=='MBW']
axm.scatter(_df_mbw.lon, _df_mbw.lat, c='cyan',
            s=6**2, marker='v', edgecolors='k',linewidths=0.5,
            zorder=91, transform=ccrs.PlateCarree())

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

cbh.set_ticks([1, 3, 5, 7], labels=['$\leq$0','10','20','$\geq$30'])
cbh.ax.invert_yaxis()

for _etype, _df in etype_dfs.items():
    # Plot events with scaling
    handles = mutil.plot_events(axm, _df, depth_cmap=depth_cmap)

    # Plot example event
    if _etype == 'all':
        pass
    else:
        _pref_event = _df[_df.event_id == event_id_fstr.format(evid=pref_evid[_etype])]
        # Place marker
        shdl = axm.scatter(_pref_event.lon, _pref_event.lat,
                           marker=marker_map[_etype],
                           c='cyan',
                           edgecolors=mec_map[_etype],
                           s=mutil.magscale(_pref_event.mag),
                           vmin=1, vmax=len(dbins),
                           transform=ccrs.PlateCarree(),
                           zorder=100)
        handles.append(shdl)
        # Place connector line to MBW
        rhdl = axm.plot([_pref_event.lon, _df_mbw.lon],
                        [_pref_event.lat, _df_mbw.lat],
                        'c-', transform=ccrs.PlateCarree(),
                        zorder=2)
        handles.append(rhdl)

    # Make relevant magnitude Legend
    if _etype == 'all':
        magax = fig.add_subplot(gs[9:13, 17:])
    else:
        magax = fig.add_subplot(gs[9:11, 17:])
    xpos = 0
    for _m in range(0, 6):
        if _etype == 'all':
            for _e, _k in enumerate(label_map.keys()):
                if _m <= np.ceil(_df[_df.etype==_k].mag.max()):
                    _s= mutil.magscale(_m)

                    magax.scatter(_m, 5 - _e, s=_s, c=2, alpha=base_alpha,
                                vmin=1,vmax=len(dbins), cmap=depth_cmap,
                                linewidths=lw_map[_k],
                                marker=marker_map[_k],
                                edgecolors=mec_map[_k])
        else:
            if _m <= np.ceil(_df.mag.max()):
                _s = mutil.magscale(_m)
                magax.scatter(_m, 1, s=_s, c=2, alpha=base_alpha,
                              vmin=1, vmax=len(dbins), cmap=depth_cmap,
                              linewidths=lw_map[_etype],
                              marker=marker_map[_etype],
                              edgecolors=mec_map[_etype])
    if _etype == 'all':
        magax.set_ylim([1,6])
        magax.yaxis.set_ticks_position('right')
        magax.yaxis.set_label_position('right')
        magax.set_ylabel('Event Type', rotation=270, labelpad=15)
        magax.yaxis.set_ticks([5, 4, 3, 2], labels=['EQ','LF','SU','PX'])
    else:
        magax.set_ylim([-1, 3])
        magax.yaxis.set_visible(False)
    

    magax.set_xlim([-0.5, 6.25])
    magax.set_xlabel('Magnitude')
    magax.xaxis.set_ticks([0,1,2,3,4,5],labels=['$\leq$0','1','2','3','4','5'])
    breakpoint()
    # IF saving
    if issave:
        fig.savefig(str(SPATH/f'overview_map_{_etype.upper()}_{DPI}dpi.{FMT}'), dpi=DPI, format=FMT)
    # Clearout magax
    magax.remove()
    # Remove previous events
    for _h in handles:
        try:
            _h.remove()
        except:
            list(_h)[0].remove()


if isshow:
    plt.show()