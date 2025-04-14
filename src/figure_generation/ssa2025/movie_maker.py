import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from obspy import UTCDateTime
from obspy.clients.fdsn import Client

import cartopy.crs as ccrs

import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
# Absolute paths to input data files
CPD = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
SHPDIR = ROOT/'data'/'SHP'
NFSHP = SHPDIR/'USDA'/'Forest_Administrative_Boundaries_(Feature_Layer)'/'Forest_Administrative_Boundaries_(Feature_Layer).shp'
NWSHP = SHPDIR/'USDA'/'Mount_Baker_Wilderness'/'Mount_Baker.shp'
RGISHP = SHPDIR/'RGI7'/'RGI2000-v7'/'RGI2000-v7.0-G-02_western_canada_usa.shp'
# NBSHP = SHPDIR/'CENSUS'/'cb_2018_us_nation_20m'/'cb_2018_us_nation_20m.shp'
# LMD = ROOT/'results'/'tables'/'SSA2025'/'linkmat_average_cct0.300.npy'
# TPK = ROOT/'results'/'tables'/'SSA2025'/'template_profile.csv'
# Save Directory for Figures
SDIR = ROOT/'results'/'figures'/'SSA2025'


## OUTPUT CONTROL ##
# Figure output control
FMT = 'png'
issave = True
DPI = 120
# Movie frame output control
makemovie = True
MDPI = 100
T0 = pd.Timestamp('1980-01-01')
T1 = pd.Timestamp('2025-03-01')
DTWIN = pd.Timedelta(180, unit='d')
DTSLIDE = pd.Timedelta(7, unit='d')


# Create key word argument holder for all fig.savefig calls
sfkw = {'format': FMT, 'dpi': DPI, 'bbox_inches':'tight'}

## Font Size Control
plt.rcParams.update({'font.size': 14})

# Marker Rendering by Event Type
eveoff = 6**2
evebase = 3.2
marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
color_map = {'px':'m','su':'b','lf':'r','eq':'k'}
ungrouped_color = 'xkcd:ugly brown'
unanalyzed_color = 'xkcd:gross green'

# Initialize Basemap Figure and Gridspec
fig = plt.figure(figsize=(6.55, 8.07))
gs = fig.add_gridspec(ncols=1, nrows=5, hspace=0, wspace=0)
# Define SpecShape
map_sps = gs[:4] # Left column, full
tz_sps = gs[4]  # Right column, upper


### SUPPORT FUNCTIONS ###
def plotmap(map_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, base=evebase, offset=eveoff):
    hdl = map_axis.scatter(df.lon, df.lat,
                           marker=marker,
                           c=colors,
                           s=mutil.magscale(df.mag, base=base, offset=offset),
                           alpha=alpha,
                           zorder=zorder,
                           label=label,
                           transform=ccrs.PlateCarree())
    handles = [hdl]
    return handles

def plottz(tz_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, **options):
    hdl = tz_axis.scatter(df.lon, df.lat,
                           marker=marker,
                           c=colors,
                           s=mutil.magscale(df.mag, **options),
                           alpha=alpha,
                           zorder=zorder,
                           label=label)
    return hdl


def plotfun(map_axis, tz_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, base=evebase, offset=eveoff):
    handles = []
    # Plot events on map
    hdl = map_axis.scatter(df.lon, df.lat,
                      marker=marker,
                      c=colors,
                      s=mutil.magscale(df.mag, base=base, offset=offset),
                      alpha=alpha,
                      zorder=zorder,
                      transform=ccrs.PlateCarree())
    handles.append(hdl)
    # Plot events on time-depth chart
    hdl = tz_axis.scatter(df.prefor_time, df.depth*1e-3,
                      c=colors,
                      s=mutil.magscale(df.mag, base=base, offset=offset),
                      marker=marker,
                      alpha=alpha,
                      zorder=zorder,
                      label=label)
    handles.append(hdl)
    return handles

def clear_handles(handles):
    for _h in handles:
        try:
            _h.remove()
        except:
            list(_h).remove()
    return []



### LOADING SECTION ###
# Load Catalog Profile
df_cat = pd.read_csv(CPD, parse_dates=['prefor_time'], index_col=[0])
# Subset to 30 km radius
df_cat = df_cat[df_cat.CAT0]
# Remove Unknown event
df_cat = df_cat[df_cat.etype != 'uk']

# Get stations within 0.5 degrees
inv = Client('IRIS').get_stations(latitude=mutil.BAKER_LAT, longitude=mutil.BAKER_LON,
                                  maxradius=0.5, level='station',
                                  network='UW,CN,CC,TA,GS')
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
# Load Census Bureau 2018 national boundaries
# gdf_cb = mutil.gpd.read_file(NBSHP)
# Load Randolph Glacier Inventory
gdf_rgi = mutil.gpd.read_file(RGISHP)
# Subset to glaciers within 0.2 degrees of Mount Baker & dont include WA grab-bag
gdf_rgi = gdf_rgi[(abs(gdf_rgi.cenlat - mutil.BAKER_LAT) < 0.2) & 
                  (abs(gdf_rgi.cenlon - mutil.BAKER_LON) < 0.2)]# &
                #   (gdf_rgi.glac_name != 'WA')]
# breakpoint()



########################
#### STATIC BASEMAP ####
########################

# Initialize geoaxes and get AWS hillshade basemap (very light)
axm, attr = mutil.mount_baker_basemap(
    fig=fig, sps=map_sps, radius_km=31,
    latnudge=0, lonnudge=0,
    open_street_map=False,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})

# Add US/Canada Border (49th parallel)
axm.plot([-123, -121], [49,49], 'k-', transform=ccrs.PlateCarree())
# Label State/Provinces
axm.text(-122.21, 49.005, 'BC', fontsize=14, transform=ccrs.PlateCarree())
axm.text(-122.21, 48.97, 'WA', fontsize=14, transform=ccrs.PlateCarree())

# Add Mount Baker-Snoqualmie NF boundary
hdl = mutil.plot_gdf_contents(axm, gdf_mbsnf, transform=None, alpha=0.1, edgecolor='forestgreen')
# Add Mount Baker Wilderness boundary
hdl = mutil.plot_gdf_contents(axm, gdf_mbw, transform=None, alpha=0.2, edgecolor='forestgreen')

# Add distance Rings
mutil.add_rings(axm,
                rads_km=[10,20,30],
                rads_colors=['k']*3,
                include_units=[True, True, True],
                label_pt=18,ha='left',va='bottom')

# Add Lat/Lon annotations
gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
                   xlocs=[-122.1, -121.8, -121.5],
                   ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
                   alpha=0)
                #    xlocs=[-122.2, -122, -121.8, -121.6, -121.4, -121.2],
                #    ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
                #    alpha=0)

# Plot area stations
sta_hdl = axm.scatter(
    df_inv.lon, df_inv.lat,
    marker='v', c='w',
    edgecolors='k',
    s=7**2, zorder=10,
    transform=ccrs.PlateCarree())

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

## SUBPLOT 2
axt = fig.add_subplot(tz_sps)
axt.yaxis.set_ticks_position('right')
axt.yaxis.set_label_position('right')
axt.set_ylim([45, -3])
axt.set_xlim([pd.Timestamp('1978-01-01'), pd.Timestamp('2026-01-01')])
axt.set_xlabel('Time [year]')
axt.set_ylabel('Depth [km]', rotation=270, labelpad=15)
## Create Movie
    

# breakpoint()
sfkw.update({'dpi': MDPI})

# Make initial window bounds
t0 = T0 - DTWIN
t1 = T0
# breakpoint()
# Get total frame count (for status reports to stdout)
total_frames = int(np.ceil((T1 - T0)/DTSLIDE))
# Set initial frame number
frame_no = 0
# March forward in time
hdls = []

# Remove semi-temporary stations
sta_hdl.remove()

# Plot event depths and times 
zo=1
for _pet in ['eq','su','lf','px']:
    _df = df_cat[df_cat.petype==_pet]
    # __df = _df[_df.tidy_group==-2]
    # hdl = axt.scatter(__df.prefor_time, __df.depth*1e-3,
    #                 marker=marker_map[_pet],
    #                 c=unanalyzed_color,
    #                 s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
    #                 alpha=0.6)
    # __df = _df[_df.tidy_group==-1]
    # hdl = axt.scatter(__df.prefor_time, __df.depth*1e-3,
    #                 marker=marker_map[_pet],
    #                 c=ungrouped_color,
    #                 s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
    #                 alpha=0.6)
    __df = _df[_df.tidy_group > -3]
    hdl = axt.scatter(__df.prefor_time, __df.depth*1e-3,
                    marker=marker_map[_pet],
                    c=color_map[_pet],
                    s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
                    alpha=0.6)
    # tz_hld = plottz(axt, _df[_df.tidy_group==-2],
    #                 marker_map[_pet],
    #                 unanalyzed_color,
    #                 zorder=1+zo)
    # tz_hdl = plottz(axt, _df[_df.tidy_group==-1],
    #                 marker_map[_pet],
    #                 ungrouped_color,
    #                 zorder=2+zo)
    # tz_hdl = plottz(axt, _df[_df.tidy_group > -1],
    #                 marker_map[_pet],
    #                 color_map[_pet],
    #                 zorder=3+zo)
    zo += 1

# breakpoint()
# Loop over time to create updating slider & event distributions
while t1 <= T1:
    print(f'frame {frame_no} of {total_frames:d}')

    # Plot slider on time figure
    shdl = axt.fill_between([t0, t1], [50,50], [-5,-5], color='cyan', alpha=0.25)
    hdls.append(shdl)
    axt.set_ylim([45, -3])

    # Subset all events for map
    _df = df_cat[(df_cat.prefor_time > t0) & (df_cat.prefor_time <= t1)]


    # Subset stations
    _inv = inv.select(time=UTCDateTime(t1.timestamp()))

    # Plot active stations
    ilats, ilons = [], []
    for net in _inv.networks:
        for sta in net.stations:
            ilats.append(sta.latitude)
            ilons.append(sta.longitude)
    ihdl = axm.scatter(ilons, ilats, marker='v', edgecolors='k', c='w', s=7**2,
                        transform=ccrs.PlateCarree(), zorder=4)
    hdls.append(ihdl)
    # breakpoint()
    # Plot Events
    for _pet in ['eq','su','lf','px']:
        # Plot unanalyzed
        # __df = _df[(_df.petype==_pet)&(_df.tidy_group==-2)]
        # if len(__df) > 0:
        #     __alphas = [1 - (t1 - row.prefor_time)/DTWIN for _, row in __df.iterrows()]
        #     hdl = axm.scatter(__df.lon, __df.lat,
        #                       marker=marker_map[_pet],
        #                       c=unanalyzed_color,
        #                       s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
        #                       alpha=__alphas,
        #                       transform=ccrs.PlateCarree())
        #     hdls.append(hdl)

        # # Plot ungrouped
        # __df = _df[(_df.petype==_pet)&(_df.tidy_group==-1)]
        # if len(__df) > 0:
        #     __alphas = [1 - (t1 - row.prefor_time)/DTWIN for _, row in __df.iterrows()]
        #     hdl = axm.scatter(__df.lon, __df.lat,
        #                       marker=marker_map[_pet],
        #                       c=ungrouped_color,
        #                       s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
        #                       alpha=__alphas,
        #                       transform=ccrs.PlateCarree())
        #     hdls.append(hdl)

        # Plot grouped
        __df = _df[(_df.petype==_pet)&(_df.tidy_group>-3)]
        if len(__df) > 0:
            __alphas = [1 - (t1 - row.prefor_time)/DTWIN for _, row in __df.iterrows()]
            hdl = axm.scatter(__df.lon, __df.lat,
                              marker=marker_map[_pet],
                              c=color_map[_pet],
                              s=mutil.magscale(__df.mag, base=evebase, offset=eveoff),
                              alpha=__alphas,
                              transform=ccrs.PlateCarree())
            hdls.append(hdl)


    # Save frame
    if issave:
        savename = SDIR/'reclassified_movie'/f'frame_{frame_no:05d}_{MDPI:d}dpi.{FMT}'
        fig.savefig(str(savename), **sfkw)
    # Advance timer
    t0 += DTSLIDE
    t1 += DTSLIDE

    # Advance frame number
    frame_no += 1
    # breakpoint()
    # Clear temporary items
    hdls = clear_handles(hdls)