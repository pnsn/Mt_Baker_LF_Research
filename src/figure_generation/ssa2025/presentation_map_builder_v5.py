"""
TODO: Add attribution
TODO: Fix marker rendering on "applied relabeling"
TODO: Move relabeling into preprocessing (should carry in with CPD file)
"""

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
makemovie = False
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
marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
color_map = {'px':'m','su':'b','lf':'r','eq':'k'}
ungrouped_color = 'xkcd:ugly brown'
unanalyzed_color = 'xkcd:gross green'

# Initialize Basemap Figure and Gridspec
fig = plt.figure(figsize=(15.12,7.4))
gs = fig.add_gridspec(ncols=2, nrows=2, hspace=0, wspace=0)
# Define SpecShape
map_sps = gs[:,0] # Left column, full
tz_sps = gs[0,1]  # Right column, upper





### SUPPORT FUNCTIONS ###
def plotmap(map_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, base=3.2, offset=25):
    hdl = map_axis.scatter(df.lon, df.lat,
                           marker=marker,
                           c=colors,
                           s=mutil.magscale(df.mag, base=base, offset=offset),
                           alpha=alpha,
                           zorder=zorder,
                           label=label,
                           transform=ccrs.PlateCarree())
    return hdl

def plottz(tz_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, **options):
    hdl = tz_axis.scatter(df.lon, df.lat,
                           marker=marker,
                           c=colors,
                           s=mutil.magscale(df.mag, **options),
                           alpha=alpha,
                           zorder=zorder,
                           label=label)
    return hdl


def plotfun(map_axis, tz_axis, df, marker, colors, alpha=0.6, zorder=1, label=None, base=3.2, offset=25):
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


# axm = fig.add_subplot(map_sps, projection=ccrs.PlateCarree())
# extent = mutil.rad2llur(rad=31e3)
# axm.set_extent(extent, ccrs.PlateCarree())
# Add US/Canada Border (49th parallel)
axm.plot([-123, -121], [49,49], 'k-', transform=ccrs.PlateCarree())
# axm.add_feature(mutil.cfeature.BORDERS, linestyle='-')
# axm.add_feature(mutil.cfeature.BOUNDARIES, linestyle='-')
# Label State/Provinces
axm.text(-122.21, 49.005, 'BC', fontsize=14, transform=ccrs.PlateCarree())
axm.text(-122.21, 48.97, 'WA', fontsize=14, transform=ccrs.PlateCarree())

# # Add Mount Baker-Snoqualmie NF boundary
hdl = mutil.plot_gdf_contents(axm, gdf_mbsnf, transform=None, alpha=0.1, edgecolor='forestgreen')
# hdl = mutil.plot_gdf_contents(axm, gdf_mbw, transform=None, alpha=0.2, edgecolor='forestgreen')
# hdl = mutil.plot_gdf_contents(axm, gdf_rgi, transform=gdf_mbsnf.crs, alpha=0.2, facecolor='dodgerblue', linewidth=0.5, edgecolor='white')
# Add Mount Baker Wilderness border
# hdl = mutil.add_shp(str(ROOT/'data'/'SHP'/'USDA'/'Mount_Baker_Wilderness'/'Mount_Baker.shp'),
#                     axm, facecolor='k')
# hdl = mutil.plot_gdf_contents(axm, gdf_cb, transform='geoaxis', alpha=1, linewidth=1)

# Add distance Rings
mutil.add_rings(axm,
                rads_km=[10,20,30],
                rads_colors=['k']*3,
                include_units=[True, True, True],
                label_pt=18,ha='left',va='bottom')

# Add Lat/Lon annotations
# axm.set_xticks([-122.2 + _e*0.2 for _e in range(5)])
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
    s=36, zorder=10,
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
# Insert basemap attribution

## SUBPLOT 2
axt = fig.add_subplot(tz_sps)
axt.yaxis.set_ticks_position('right')
axt.yaxis.set_label_position('right')
axt.set_ylim([45, -3])
axt.set_xlim([pd.Timestamp('1978-01-01'), pd.Timestamp('2026-01-01')])
axt.set_xlabel('Time [year]')
axt.set_ylabel('Depth [km]', rotation=270, labelpad=15)


#### TEMPORARY PLOTTED ELEMENTS ####

##############
#### DATA ####
##############
#### Overview Map
hdls = []
mountain_labels = zip(['Mount Baker','South Twin','Mount Shuksan'],
                      [(mutil.BAKER_LON, mutil.BAKER_LAT),
                       (mutil.STWIN_LON, mutil.STWIN_LAT),
                       (mutil.SHUKS_LON, mutil.SHUKS_LAT)])
for name, (lon, lat) in mountain_labels:
    hdl = axm.text(lon, lat, name, transform=ccrs.PlateCarree(), ha='left', va='bottom')
    hdls.append(hdl)

hdl = axm.text(-121.96, 48.61, 'Contributing Seismic Stations\n1980-2025 Composite',
               transform=ccrs.PlateCarree(), ha='center', va='center',
               color='k', fontweight='extra bold')
hdls.append(hdl)
hdl = axm.text(-121.8133, 48.905, 'Mount Baker-\nSnoqualmie\nNational\nForest',
               transform=ccrs.PlateCarree(), ha='center', va='center',
               color='forestgreen', fontweight='extra bold')
hdls.append(hdl)
hdl = axt.text(pd.Timestamp('2002-06-15'), 20, 'Event Depths and Origin Times Go Here\nMarkers Scaled by Magnitude (M0-M5)',
               ha='center',va='center', fontsize=16)
hdls.append(hdl)

# Save Figure
if issave:
    fig.savefig(str(SDIR/f'Map_and_TimeDepth_Overview_{DPI:d}dpi.{FMT}'),
                **sfkw)
else:
    plt.show()
    breakpoint()
# Clear handles
hdls = clear_handles(hdls)

#### PNSN Catalog Summary Figure ####
zo = 5 # Keep below 10 to underlie stations & baker
hdls = []
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    handles = plotfun(axm, axt, _df, marker_map[_et], color_map[_et],
                      zorder=zo, label=f'{_et.upper()}')#, base=4,offset=16)
                    #   zorder=zo,label=f'{_et.upper()}: {len(_df)}')
    hdls += handles
    zo -= 1
lhdl = axt.legend(handles=hdls[::-1], loc='lower center', ncols=4)
hdls.append(lhdl)

# Save Figure
if issave:
    fig.savefig(str(SDIR/f'PNSN_current_catalog_{DPI:d}dpi.{FMT}'),
                **sfkw)
else:
    plt.show()
    breakpoint()
# Clear temporary elements
clear_handles(hdls)

### PNSN Catalog Event Type Specific Catalogs
zo = 20 # Overlie stations, underlie baker
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    hdls = plotfun(axm, axt, _df, marker_map[_et], color_map[_et], 
                   zorder=zo,label=f'{_et.upper()}: {len(_df)}')
    lhdl = axt.legend(loc='lower center')
    hdls.append(lhdl)

    # Save Figure
    if issave:
        fig.savefig(str(SDIR/f'PNSN_current_catalog_{_et.upper()}_{DPI:d}dpi.{FMT}'),
                    **sfkw)
    # Clear temporary elements
    clear_handles(hdls)

#################
#### RESULTS ####
#################


### SHOW WHAT'S GROUPED ###
# Plot all groups
zo = df_cat.tidy_group.max() + 200 # Overlie stations & baker
hdls = []

# Iterate across event types for marker shape
for _et in ['eq','lf','su','px']:
    # Iterate across groups for marker color
    for _grp in df_cat.tidy_group.value_counts(ascending=True).index:
        _df = df_cat[(df_cat.tidy_group==_grp) & (df_cat.etype==_et)]
        if len(_df) > 0:
            if _grp == -2:
                color = unanalyzed_color
            elif _grp == -1:
                color = ungrouped_color
            else:
                color = _df.leaf_color.iloc[0]
        
            handles = plotfun(
                axm, axt, _df,
                marker_map[_et],
                color,
                zorder=zo,
            )
            hdls += handles
            zo -= 1
        else:
            continue

# Save
savename = SDIR/f'map_clustered_events_{DPI}dpi.{FMT}'
if issave:
    fig.savefig(str(savename), **sfkw)
#Cleanup
hdls = clear_handles(hdls)

### SHOW WHAT'S MISSED ###
# Plot Only unanalyzed Events with catalot etype markings
zo = 20
_df = df_cat[df_cat.tidy_group == -2]
for _et in ['px','lf','eq']:
    __df = _df[_df.etype==_et]
    handles = plotfun(
        axm, axt, __df,
        marker_map[_et],
        unanalyzed_color,
        zorder=zo, label=f'{_et.upper()}: {len(__df):d}')
    hdls += handles
    zo -= 1
lhdl = axt.legend(handles=hdls[::-1], loc='lower center', ncols=3)
hdls.append(lhdl)

if issave:
    fig.savefig(str(SDIR/f'map_unanalyzed_events_{DPI}dpi.{FMT}'), **sfkw)
hdls = clear_handles(hdls)
# Plot Only Ungrouped Events with catalog etype markings
zo = 20 # overlie stations, underlie baker
_df = df_cat[df_cat.tidy_group == -1]
# Iterate across etype
for _et in ['su','lf','px','eq']:
    __df = _df[_df.etype==_et] 
    handles = plotfun(
        axm, axt, __df, 
        marker_map[_et], 
        ungrouped_color, 
        zorder=zo, label=f'{_et.upper()}: {len(__df):d}')
    # Include labels with counts
    hdls += handles
    zo -= 1
lhdl = axt.legend(handles=hdls[::-1], loc='lower center', ncols=4)
hdls.append(lhdl)
# Save
savename = SDIR/f'map_missed_events_{DPI}dpi.{FMT}'
# breakpoint()
if issave:
    fig.savefig(str(savename), **sfkw)
# Cleanup
hdls = clear_handles(hdls)


### PROPOSE CHANGES ### Include missed events of _ret?
zo = 20 # Overlie stations, underlie baker
# Plot each etype (re)classification subset
for _ret in ['su','lf','px','eq']:
    # Get subset events based on relabel etype
    _df = df_cat[df_cat.petype == _ret]
    # Plot unanalyzed
    handles = plotfun(
        axm, axt, _df[_df.tidy_group == -2],
        marker_map[_ret],
        unanalyzed_color,
        zorder=zo,label=f'Unanalyzed {_ret.upper()}'
    )
    hdls += handles
    # Plot ungrouped
    handles = plotfun(
        axm, axt, _df[_df.tidy_group == -1],
        marker_map[_ret],
        ungrouped_color,
        zorder=zo,label=f'Ungrouped {_ret.upper()}'
    )
    hdls += handles
    # Plot grouped that are not changed
    handles = plotfun(
        axm, axt, _df[(_df.etype==_df.petype)&(_df.tidy_group > 0)],
        marker_map[_ret],
        color_map[_ret],
        zorder=zo,label=f'Grouped {_ret.upper()}'
    )
    hdls += handles

    # Reset z-ordering indexer for each map
    _zo = zo + 1
    # Iterate over each event type for shape and labeling
    for _et, _ret in _df[['etype','petype']].value_counts().index:
        if _et == _ret:
            continue
        else:
            handles = plotfun(
                axm, axt, _df[(_df.etype==_et) & (_df.petype==_ret)],
                marker_map[_ret],
                color_map[_et],
                zorder=_zo
            )
            hdls += handles
            _zo += 1
        # first = True

        # _df = df_cat[(df_cat.petype==_ret)&(df_cat.etype != df_cat.petype)&(df_cat.etype==_et)]
        # if len(_df) > 0:
        #     handles = plotfun(
        #         axm, axt, _df,
        #         marker_map[_ret],
        #         color_map[_et],
        #         zorder=_zo,
        #     )
        #     _zo += 1
    # Save
    savename = SDIR/f'map_{_ret}_proposed_relabeling_{DPI}dpi.{FMT}'
    if issave:
        fig.savefig(str(savename), **sfkw)
    else:
        breakpoint()
    # Cleanup
    hdls = clear_handles(hdls)


    # # # Iterate over each group for color
    # #     for _id in _ids:
    # #         __df = _df[_df.tidy_group==_id]
    # #         if len(__df) > 0:
    # #             if first:
    # #                 if _et == _ret:
    # #                     label = f'Grouped {_et.upper()}'
    # #                 else:
    # #                     label=f'{_et.upper()}➔{_ret.upper()}'
    # #                 first = False
    # #             else:
    # #                 label = None
    # #             handles = plotfun(
    # #                 axm, axt, __df,
    # #                 marker_map[_et],
    # #                 __df.leaf_color.iloc[0],
    # #                 zorder=_zo,
    # #                 label=label)
    # #             _zo += 1
    # #             hdls += handles
    # # lhdl = axt.legend(loc='lower center', ncols=5)
    # hdls.append(lhdl)



# ### ENACT CHANGES ###
# zo = 20
# # Plot each etype (re)classification subset
# for _ret, _ids in zip(['su','lf','px','eq'],[SU_set, LF_set, PX_set, EQ_set]):
#     # Plot ungrouped of target reassignment
#     _df = df_cat[(df_cat.tidy_group == -1) & (df_cat.etype==_ret)]
#     handles = plotfun(
#         axm, axt, _df,
#         marker_map[_ret],
#         ungrouped_color,
#         zorder=zo,label=f'Ungrouped {_ret.upper()}'
#     )
#     hdls += handles
#     _zo = zo + 1
#     # Iterate over each event type for shape and labeling
#     for _et in ['eq','lf','su','px']:
#         first = True
#         _df = df_cat[(df_cat.tidy_group.isin(_ids))&(df_cat.etype==_et)]
#     # Iterate over each group for color
#         for _id in _ids:
#             __df = _df[_df.tidy_group==_id]
#             if len(__df) > 0:
#                 if first:
#                     label=f'{_et.upper()}➔{_ret.upper()}'
#                     first = False
#                 else:
#                     label = None
#                 handles = plotfun(
#                     axm, axt, __df,
#                     marker_map[_et],
#                     __df.leaf_color[0],
#                     zorder=_zo,
#                     label=label)
#                 _zo += 1
#                 hdls += handles
#     lhdl = axt.legend(loc='lower center', ncols=5)
#     hdls.append(lhdl)
    
#     # Plot each etype (re)classification subset
# for _ret, _ids in zip(['su','lf','px','eq'],[SU_set, LF_set, PX_set, EQ_set]):
#     # Plot ungrouped of target reassignment
#     _df = df_cat[(df_cat.tidy_group == -1) & (df_cat.etype==_ret)]
#     handles = plotfun(
#         axm, axt, _df,
#         marker_map[_ret],
#         ungrouped_color,
#         zorder=2
#     )
#     hdls += handles
#     # Iterate over each group for color
#     _zo = 3
#     for _id in _ids:
#         _df = df_cat[df_cat.tidy_group == _id]
#         # Iterate over each catalog event type for shape
#         for _et in ['eq','lf','su','px']:
#             __df = _df[_df.etype==_et]
#             if len(__df) > 0:
#                 handles = plotfun(
#                     axm, axt, __df,
#                     marker_map[_et],
#                     color_map[_ret],
#                     zorder=_zo)
#                 _zo += 1
#                 hdls += handles
    # savename = SDIR/f'map_{_ret}_applied_relabeling_{DPI}dpi.{FMT}'
    # if issave:
    #     fig.savefig(str(savename), **sfkw)
    # hdls = clear_handles(hdls)


## Create Movie
    
if makemovie:
    sfkw.update({'dpi': MDPI})

    # Make initial window bounds
    t0 = T0 - DTWIN
    t1 = T0
    # Get total frame count (for status reports to stdout)
    total_frames = int(np.ceil((T1 - T0)/DTSLIDE))
    # Set initial frame number
    frame_no = 0
    # March forward in time
    ihdls = []
    ahdls = []
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

        # Subset all events for tz plot
        # Iterate across event types for shapes
        _zo = 2
        for _et, _set in zip(['eq','su','lf','px'],[EQ_set, SU_set, LF_set, PX_set]):
            # Get etype subset by assignment group
            # Get ungrouped events of target etype
            __df_x = _df[(_df.etype == _et)&(~_df.tidy_group.isin(_set))]
            # Get grouped events reassigned to etype
            __df_i = _df[_df.tidy_group.isin(_set)]
            if len(__df_x) != len(_df):
                breakpoint()
            _alphas_x = [1 - (t1 - row.prefor_time)/DTWIN for _, row in __df_x.iterrows()]
            _alphas_i = [1 - (t1 - row.prefor_time)/DTWIN for _, row in __df_i.iterrows()]
            if len(__df_x) > 0:
                # Plot ungrouped
                handles = plotfun(
                    axm, axt, __df_x,
                    marker_map[_et],
                    ungrouped_color,
                    alpha=_alphas_x,
                    zorder=1
                )
                hdls += handles
            if len(__df_i) > 0:
                # Plot grouped by intended class color
                handles = plotfun(
                    axm, axt, __df,
                    marker_map[_et],
                    color_map[_et],
                    alpha=_alphas_i,
                    zorder = _zo
                )
                _zo += 1
                hdls += handles
        
        # Save
        if issave:
            savename = SDIR/'reclassified_movie'/f'frame_{frame_no:05d}_{MDPI:d}dpi.{FMT}'
            fig.savefig(str(savename), **sfkw)
        # Advance timer
        t0 += DTSLIDE
        t1 += DTSLIDE

        # Advance frame number
        frame_no += 1
        # Clear temporary items
        hdls = clear_handles(hdls)

