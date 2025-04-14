from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import gps2dist_azimuth

import cartopy.crs as ccrs

import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
# Absolute paths to input data files
CPD = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
SHPDIR = ROOT/'data'/'SHP'
NFSHP = SHPDIR/'USDA'/'Forest_Administrative_Boundaries_(Feature_Layer)'/'Forest_Administrative_Boundaries_(Feature_Layer).shp'
NWSHP = SHPDIR/'USDA'/'Mount_Baker_Wilderness'/'Mount_Baker.shp'
# RGISHP = SHPDIR/'RGI7'/'RGI2000-v7'/'RGI2000-v7.0-G-02_western_canada_usa.shp'

SDIR = ROOT/'results'/'figures'/'SSA2025'

## OUTPUT CONTROL ##
# Figure output control
FMT = 'png'
issave = True
isshow = True
DPI = 120
# Create key word argument holder for all fig.savefig calls
sfkw = {'format': FMT, 'dpi': DPI, 'bbox_inches':'tight'}

## Font Size Control
plt.rcParams.update({'font.size': 14})

yr2s = 3600*24*365.24

radius = 100e3
alpha_radius = 120e3
# Load National Forest Boundaries
gdf_mbsnf = mutil.gpd.read_file(NFSHP)
# Subset to Mount Baker-Snoqualmie National Forest ShapeFile contents
gdf_mbsnf = gdf_mbsnf[gdf_mbsnf.FORESTNAME=='Mt. Baker-Snoqualmie National Forest']
# Load Mount Baker Wilderness Shapefile
gdf_mbw = mutil.gpd.read_file(NWSHP)


# Get stations within 100 km
inv = Client('IRIS').get_stations(latitude=mutil.BAKER_LAT, longitude=mutil.BAKER_LON,
                                  maxradius=radius/111.2e3, level='channel',
                                  network='UW,CN,CC,TA,GS',
                                  channel='HHZ,EHZ,ENZ,BHZ,HNZ')
# Extract station location / net-sta codes & on/off dates
holder = []

future_utc = UTCDateTime('2026-01-01')

for net in inv.networks:
    for sta in net.stations:
        if sta.end_date is None:
            _ed = future_utc
        else:
            _ed = sta.end_date
        line = [net.code, sta.code, sta.longitude, sta.latitude,
                pd.Timestamp(sta.start_date.timestamp, unit='s'), pd.Timestamp(_ed.timestamp, unit='s'),
                _ed - sta.start_date]
        holder.append(line)

df_ins = pd.DataFrame(holder,columns=['net','sta','lon','lat','firston','lastoff','duration'])


### PLOTTING SECTION ####

# Initialize Basemap Figure and Gridspec
fig = plt.figure(figsize=(7.56,7.4))
gs = fig.add_gridspec(ncols=1, nrows=31, wspace=0)
# Define SpecShape
map_sps = gs[:-1, 0] # Left column, full
# Z-error ShapeSpec
# zerr_sps = gs[:-1, 1]  # Right column, full

# Initialize Basemap
axm, attr = mutil.mount_baker_basemap(
    fig=fig, sps=map_sps, radius_km=103,
    latnudge=0, lonnudge=0,
    open_street_map=False,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})

# Canadia-USA border
# axm.plot([-124, -120], [49,49], 'k-', transform=ccrs.PlateCarree())
axm.text(-120.89, 49.05, 'BC', fontsize=14, transform=ccrs.PlateCarree())
axm.text(-120.89, 48.92, 'WA', fontsize=14, transform=ccrs.PlateCarree())

# Add NF boundaries
hdl = mutil.plot_gdf_contents(axm, gdf_mbsnf, transform=None, alpha=0.3, linewidth=1, edgecolor='forestgreen')

# Add coastlines
# axm.add_feature(mutil.cfeature.COASTLINE, color='navy')
axm.add_feature(mutil.cfeature.OCEAN, color='dodgerblue', alpha=0.3)

axm.add_feature(mutil.cfeature.BORDERS, color='k', linewidth=3, alpha=0.5)

# Plot Mount Baker, Mount Shuksan, South Twin
axm.scatter([mutil.BAKER_LON, mutil.SHUKS_LON, mutil.STWIN_LON],
            [mutil.BAKER_LAT, mutil.SHUKS_LAT, mutil.STWIN_LAT],
            s=[144, 81, 81], 
            marker='^',
            facecolor='none',
            edgecolors='orange',
            linewidths=[2, 1, 1],
            zorder=30,
            transform=ccrs.PlateCarree())

# Add distance Rings
mutil.add_rings(axm,
                rads_km=[30, 60, 90],
                rads_colors=['k']*3,
                include_units=[True, True, True],
                label_pt=-13,ha='left',va='top')

# Add Lat/Lon annotations
# axm.set_xticks([-122.2 + _e*0.2 for _e in range(5)])
gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
                   xlocs=[-123, -122, -121, -120],
                   ylocs=[48, 48.5, 49, 49.5],
                   alpha=0, xlabel_style={'fontsize':14},
                   ylabel_style={'fontsize':14})

# Plot area stations
_df_ins_a = df_ins[df_ins.lastoff > pd.Timestamp('2025-04-01')]
_df_ins_x = df_ins[df_ins.lastoff <= pd.Timestamp('2025-04-01')]
hld = axm.scatter(_df_ins_a.lon, _df_ins_a.lat,
                  marker='v',
                  c=[_e.year for _e in _df_ins_a.firston],
                  vmin=df_ins.firston.min().year,
                  vmax=df_ins.lastoff.max().year,
                  s=5+ 3*(_df_ins_a.duration/yr2s),
                  edgecolors='k',
                  transform=ccrs.PlateCarree(),
                  cmap='inferno_r', zorder=20,
                  label='Active')
hld = axm.scatter(_df_ins_x.lon, _df_ins_x.lat,
                  marker='v',
                  c=[_e.year for _e in _df_ins_x.firston],
                  vmin=df_ins.firston.min().year,
                  vmax=df_ins.lastoff.max().year,
                  s=5 + 3*(_df_ins_x.duration/yr2s),
                  edgecolors='r',
                  transform=ccrs.PlateCarree(),
                  cmap='inferno_r', zorder=15,
                  label='Decomissioned')
axm.legend(loc='upper right')
# Render colorbar
linyears = np.linspace(df_ins.firston.min().year,
                     df_ins.lastoff.max().year,
                     201)
cbax = fig.add_subplot(gs[-1,0])
cbax.scatter(linyears, np.ones(shape=linyears.shape), c=linyears, marker='s',s=8**2,cmap='inferno_r')
cbax.set_ylim([0.9999,1.0001])
cbax.set_xlim([1971, 2025])
cbax.set_xticks(np.arange(1975,2035,10), labels=np.arange(1975,2035,10, dtype=np.int32))
cbax.yaxis.set_visible(False)
cbax.set_xlabel('Deployment Year')
# Plot connecting lines
llllac = []
for _, row in df_ins.iterrows():
    line = [row.lon, mutil.BAKER_LON, row.lat, mutil.BAKER_LAT]
    distm, _, _ = gps2dist_azimuth(line[-1], line[-3], line[-2],line[-4])
    line.append(1 - distm/alpha_radius)
    if row.lastoff <= pd.Timestamp('2025-04-01'):
        line.append('r')
    else:
        line.append('k')
    llllac.append(line)

for line in llllac:
    axm.plot(line[:2], line[2:4], color=line[5], alpha=line[4], linewidth=0.5, transform=ccrs.PlateCarree())

# Annotate networks in upper left
axm.text(-123.18, 49.575, 'Included Networks:\nUW, CN, CC, TA, GS', ha='left', va='center', transform=ccrs.PlateCarree())

if issave:
    fig.savefig(str(SDIR/f'Network_Geometry_Overview_{DPI}dpi.{FMT}'), **sfkw)

if isshow:
    plt.show()