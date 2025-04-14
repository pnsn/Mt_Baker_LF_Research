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
# Shapefiles for mapping
SHPDIR = ROOT/'data'/'SHP'
NFSHP = SHPDIR/'USDA'/'Forest_Administrative_Boundaries_(Feature_Layer)'/'Forest_Administrative_Boundaries_(Feature_Layer).shp'
NWSHP = SHPDIR/'USDA'/'Mount_Baker_Wilderness'/'Mount_Baker.shp'

# Save Path
SPATH = ROOT / 'results' / 'figures' / 'SSA2025'


# Rendering Settings
DPI = 120
FMT = 'png'
plt.rcParams.update({'font.size':14})
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
depths = {}
lats = {}
lons = {}
mags = {}
for _etype, _evid in pref_evid.items():
    cat = ebank.get_events(event_id=event_id_fstr.format(evid=_evid))
    for _p in cat[0].picks:
        if _p.waveform_id.id.split('.')[:2] == pref_netsta:
            if _etype not in picks.keys():
                picks[_etype] = _p
            elif _p.time < picks[_etype].time:
                picks[_etype] = _p
    depths[_etype] = cat[0].preferred_origin().depth*1e-3
    lats[_etype] = cat[0].preferred_origin().latitude
    lons[_etype] = cat[0].preferred_origin().longitude
    mags[_etype] = cat[0].preferred_magnitude().mag
    # zerr[_etype]

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
# def plot_spect(fig, wf_sps, spect_sps, trace, spect, depth, prepad, duration, evid, color, spect_kw=spect_kw):
#     axwf = fig.add_subplot(wf_sps)
#     axsp = fig.add_subplot(spect_sps, sharex=axwf)

#     # Detrend and normalize trace
#     tr = trace.copy().detrend('linear').normalize()
#     # Plot trace
#     axwf.plot(tr.times(reftime=tr.stats.starttime),
#               tr.data,
#               color,
#               lw=0.5)
#     shdl = spectrogram(tr.data, tr.stats.sampling_rate, axes=axsp, **spect_kw)    
#     axsp.set_xticks[]     

### SPECTROGRAMS ###
figs = {}
axes = {}
caxes = {}
for _etype, _tr in traces.items():
    __tr = _tr.copy().detrend('linear').normalize()
    evid = pref_evid[_etype]
    pick = picks[_etype]
    fig = plt.figure(figsize=(4.8,5.25))
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
    axs[0].set_title(f'uw{evid} | Depth: {depths[_etype]:.1f} km\n{pick.time}')#{_etype.upper()}\n{pick.time}')
    axs[0].set_ylabel('Norm. Amp. [ - ]')
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
        _fig.savefig(str(SPATH/f'wave_spect_UW_MBW_{_etype}_uw{pref_evid[_etype]}_v6_{DPI}dpi.{FMT}'), dpi=DPI, format=FMT)    


# Make Maps of pairwise combinations with EQ

# Load National Forest Boundaries
gdf_mbsnf = mutil.gpd.read_file(NFSHP)
# Subset to Mount Baker-Snoqualmie National Forest ShapeFile contents
gdf_mbsnf = gdf_mbsnf[gdf_mbsnf.FORESTNAME=='Mt. Baker-Snoqualmie National Forest']
# Load Mount Baker Wilderness Shapefile
gdf_mbw = mutil.gpd.read_file(NWSHP)

# Get stations within 0.5 degrees
inv = Client('IRIS').get_stations(level='station',
                                  network='UW',
                                  station='MBW')
# Extract station location / net-sta codes
holder = []
for net in inv.networks:
    for sta in net.stations:
        line = [sta.longitude, sta.latitude, net.code, sta.code]
        holder.append(line)
df_inv = pd.DataFrame(holder, columns=['lon','lat','net','sta'])


# Load shapefiles into geopandas
# Load National Forest Boundaries
gdf_mbsnf = mutil.gpd.read_file(NFSHP)
# Subset to Mount Baker-Snoqualmie National Forest ShapeFile contents
gdf_mbsnf = gdf_mbsnf[gdf_mbsnf.FORESTNAME=='Mt. Baker-Snoqualmie National Forest']
# Load Mount Baker Wilderness Shapefile
gdf_mbw = mutil.gpd.read_file(NWSHP)

# Initialize Basemap Figure and Gridspec
fig = plt.figure(figsize=(5,5))
gs = fig.add_gridspec(ncols=1, nrows=1, hspace=0, wspace=0)
# Define SpecShape
map_sps = gs[0] # Left column, full

# Initialize geoaxes and get AWS hillshade basemap (very light)
axm, attr = mutil.mount_baker_basemap(
    fig=fig, sps=map_sps, radius_km=15,
    latnudge=df_inv.lat.iloc[0] - mutil.BAKER_LAT + 0.03, 
    lonnudge=df_inv.lon.iloc[0] - mutil.BAKER_LON,
    open_street_map=False,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})
# Add distance Rings
mutil.add_rings(axm,
                rads_km=[10,15,20],
                rads_colors=['k']*3,
                include_units=[True, True, True],
                label_pt=55,ha='left',va='bottom', fontsize=12)

gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
                   xlocs=[-122.0, -121.9, -121.8],
                   ylocs=[48.7, 48.8, 48.9, 49.],
                   xlabel_style={'fontsize':14},
                   ylabel_style={'fontsize':14,'rotation':270},
                   alpha=0)



# # Add Mount Baker-Snoqualmie NF boundary
hdl = mutil.plot_gdf_contents(axm, gdf_mbsnf, transform=None, alpha=0.1, edgecolor='forestgreen')



# Plot Mount Baker
axm.scatter([mutil.BAKER_LON, mutil.SHUKS_LON, mutil.STWIN_LON],
            [mutil.BAKER_LAT, mutil.SHUKS_LAT, mutil.STWIN_LAT],
            s=[144, 81, 81], 
            marker='^',
            facecolor='none',
            edgecolors='orange',
            linewidths=[2, 1, 1],
            zorder=30,
            transform=ccrs.PlateCarree())


# Plot station
axm.scatter(df_inv.lon, df_inv.lat, marker='v', c='cyan',edgecolors='k', s=7**2,
            transform=ccrs.PlateCarree(),zorder=5)
axm.text(df_inv.lon[0]+0.005, df_inv.lat[0]+0.005, 'UW.MBW', fontsize=14, ha='left', va='bottom',
         transform=ccrs.PlateCarree())

# Plot Events
for _et in ['eq','lf','su','px']:
    axm.scatter(lons[_et], lats[_et], s=3.2**mags[_et] + 36, marker=marker_map[_et],
                c=mec_map[_et], transform=ccrs.PlateCarree(), zorder=6)
    axm.plot([lons[_et], df_inv.lon[0]], [lats[_et], df_inv.lat[0]],
             color=mec_map[_et], alpha=0.5, zorder=1, transform=ccrs.PlateCarree())



if issave:
    fig.savefig(str(SPATH/f'Example_Event_Map_{DPI}dpi.{FMT}'),
                format=FMT, dpi=DPI, bbox_inches='tight')

if isshow:
    plt.show()
# Plot 