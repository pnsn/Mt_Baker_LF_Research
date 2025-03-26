import os
from pathlib import Path

import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees

import shapely.geometry as sgeom
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from map_util import *

# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to catalog membership CSV
CATD = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'
# path to preferred event/channel pairs CSV
PESD = ROOT / 'processed_data' / 'catalog' / 'P1S2_Preferred_Sta_Event_Picks.csv'

# SAVEPATH
SAVEPATH = ROOT / 'results' / 'figures' / 'seismolunch'
FMT = 'png'
DPI = 200

issave = True
isshow = True
map_rad_km = 100.
zoom = 9

RADII_KM = [10, 30, 50, 70, 90]
RADII_COL = ['black','firebrick','darkgreen','darkblue','purple']

# Network Code(s) to include
NETS = ['UW','CN','TA','PB','CC','UO']
# Channel Codes to include
CHANS = []
for _b in 'BHE':
    for _i in 'HN':
        for _c in 'ZNE12':
            CHANS.append(''.join([_b,_i,_c]))
INV_LEVEL = 'station'

PREFSTA = ['VDB','CMW','HDW','HTW','JCW','MBW','MULN','PASS','RPW','SAXON','SHUK']

## PREPROCESSING SECTION ###
# Connect to client
IRIS = Client('IRIS')
# Get station inventory
INV = IRIS.get_stations(latitude=BAKER_LAT,
                        longitude=BAKER_LON,
                        maxradius=map_rad_km/111.2,
                        network=','.join(NETS),
                        channel=','.join(CHANS),
                        level=INV_LEVEL)
# Create station table
holder = []
for net in INV.networks:
    for sta in net.stations:
        # Ignore backup telemetered station
        if sta.code == 'MB2':
            continue
        # If at least 2 years of data
        start = pd.Timestamp(sta.start_date.datetime)
        if sta.end_date is None:
            end = pd.Timestamp('2030-01-01')
            isactive = True
            delta = UTCDateTime() - sta.start_date
        else:
            isactive = sta.end_date > UTCDateTime() 
            delta = sta.end_date - sta.start_date
            end = pd.Timestamp(sta.end_date.datetime)
        dist = locations2degrees(BAKER_LAT, BAKER_LON, sta.latitude, sta.longitude)
        dist *= 111.2
        line = [net.code, sta.code, 
                start, end,
                sta.latitude, sta.longitude, sta.elevation,
                isactive, delta/(3600*24*365.24), sta.code in PREFSTA,
                dist]
        holder.append(line)
df_sta = pd.DataFrame(holder, columns=['net','sta','start','end','lat','lon','elev','isactive','runyears','ispref', 'offset_km'])            

df_sta = df_sta.sort_values(['runyears'], ascending=False)



# Initialize Figure
fig = plt.figure(figsize=(8, 8))
gs = fig.add_gridspec(ncols=6, nrows=6)
# Initialize Map Figure
axm, map_attrib = mount_baker_basemap(fig=fig, sps=gs[:,:], radius_km=map_rad_km,
                                      open_street_map=True, zoom=zoom,
                                      aws_add_image_kwargs={'cmap': 'Greys_r', 'alpha': 0.25})

# Plot stations by deploy date
sax = axm.scatter(
    df_sta.lon, df_sta.lat,
    c=[x.year for x in df_sta.start],
    s=1.7**round(4 + df_sta.runyears/10),
    vmin=1970, vmax=2030,
    marker='v',
    edgecolor = ['k' if x else 'r' for x in df_sta.isactive],
    linewidth=0.5,
    cmap='inferno_r', transform=ccrs.PlateCarree(),
    zorder=10)
# Station deployment year colorbar
cbinset = inset_axes(axm, width="25%", height='3%', loc='lower right', borderpad=1)
cbh = fig.colorbar(sax, cax=cbinset, orientation='horizontal', 
                   shrink=0.5, location='bottom', ticklocation='top')
cbh.ax.set_xlabel('Deployment Year')
# Overlay triangles for still-active stations
# _dfa = df_sta[~df_sta.isactive]
# axm.plot(_dfa.lon, _dfa.lat, 'v', ms=6, color=None, mec='c', linewidth=0.5, alpha=0.5, zorder=11, transform=ccrs.PlateCarree())
# Plot Mount Baker Location
axm.plot(BAKER_LON, BAKER_LAT, '^', ms=6, mec='k', mfc='w', transform=ccrs.PlateCarree())
axm.text(BAKER_LON + 0.05, BAKER_LAT, 'Mount\nBaker',va='top',transform=ccrs.PlateCarree())

# Label stations with > 10 years of records
for _, row in df_sta.iterrows():
    if row.ispref:
        axm.text(row.lon + 0.05, row.lat, f'{row.net}.{row.sta}', 
                 fontsize=6, va='bottom', transform=ccrs.PlateCarree())
    elif row.runyears > 20:
        axm.text(row.lon + 0.05, row.lat, f'{row.net}.{row.sta}',
                 fontsize=4, va='bottom', color='b', transform=ccrs.PlateCarree())

# Add attribution
xlim = axm.get_xlim()
ylim = axm.get_ylim()
axm.text(xlim[0], ylim[1]-1e3, '\n'.join(map_attrib), ha='left', va='top',
         fontsize=6)

# Add coordinate labels
gl = axm.gridlines(draw_labels=True, zorder=1)
gl.top_labels = False
gl.left_labels = False
gl.xlines = False
gl.ylines = False
# Plot radii
for _e, _r in enumerate(RADII_KM):
    # Plot circles
    mE, mN = radiusllsets(rad=_r*1e3)
    axm.plot(mE, mN, ':', color=RADII_COL[_e], alpha=0.667, transform=UTM10N)
    # Label Radii
    axm.text(mE[28], mN[28]+1e3, f'{_r:d} km', va='bottom', ha='center', transform=UTM10N)



# Initialize inset map
axs = fig.add_axes([0.725, 0.7, 0.15, 0.15], projection=ccrs.PlateCarree())
# axs = inset_axes(axm, width='30%', height="30%", loc='upper right', borderpad=1)
# axs = projection=ccrs.PlateCarree())
# axs.set_projection=ccrs.PlateCarree()
# axs = fig.add_subplot(gs[0,0], projection=ccrs.PlateCarree())
# extent = rad2llur(rad=5*111.2e3)
axs.set_extent([BAKER_LON-10, BAKER_LON + 10, BAKER_LAT - 10, BAKER_LAT + 10], crs=ccrs.PlateCarree())
axs.add_feature(cfeature.LAND)
axs.add_feature(cfeature.OCEAN)
axs.add_feature(cfeature.COASTLINE)
axs.add_feature(cfeature.BORDERS, linestyle='-')
axs.add_feature(cfeature.STATES, linestyle=':')
axs.add_feature(cfeature.LAKES, alpha=0.5)
ext = rad2llur(rad=map_rad_km*1e3)
extbox = sgeom.box(ext[0], ext[2], ext[1], ext[3])
axs.add_geometries([extbox], ccrs.PlateCarree(), facecolor='none', edgecolor='red', linewidth=2, alpha=0.5)
axs.plot(BAKER_LON, BAKER_LAT, '^', ms=6, mec='k', mfc='w', transform=ccrs.PlateCarree())

if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'subarray_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)



fig2 = plt.figure(figsize=(9,4))
ax = fig2.add_subplot(111)

for _e, _r in enumerate(RADII_KM):
    _df = df_sta[df_sta.offset_km <= _r]
    _df.index = _df.start
    _sr = _df.copy().resample(pd.Timedelta(1, unit='w')).count().net
    counts = []
    for _t in _sr.index:
        ct = len(_df[(_df.start <= _t) & (_df.end >= _t)])
        counts.append(ct)
    # counts = [len(_df[(_df.start <= _t) & (_df.end >= _t)]) for _t in _sr.index.values]
    ax.plot(_sr.index, counts, color=RADII_COL[_e], label=f'R $\leq${_r:d} km')

ax.legend(loc='upper left')
ax.set_xlabel('Year')
ax.set_ylabel('Active Stations Within Radius "R" [no.]')
ax.set_yscale('log')
ax.set_xlim([_sr.index[0], _sr.index[-1]])
ax.set_ylim([0.9, 101])
ax.grid(which='both', linestyle=':')


if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'station_counts_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)




if isshow:
    plt.show()

