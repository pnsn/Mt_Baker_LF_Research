from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from scipy.cluster.hierarchy import dendrogram

import cartopy.crs as ccrs

import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
# Absolute paths to input data files
CPD = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
LMD = ROOT/'results'/'tables'/'SSA2025'/'linkmat_average_cct0.300.npy'
TPK = ROOT/'results'/'tables'/'SSA2025'/'template_profile.csv'
# Save Directory for Figures
SDIR = ROOT/'results'/'figures'/'SSA2025'


## OUTPUT CONTROL ##
issave = True
DPI = 120
MDPI = 100
FMT = 'png'

# Marker Rendering by Event Type
marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
color_map = {'px':'m','su':'b','lf':'r','eq':'k'}

# Dendrogram hardcoded elements
CT = 0.3
at_color = 'k'


### SUPPORT FUNCTIONS ###
def plotfun(map_axis, tz_axis, df, marker, color, alpha=0.6, zorder=1, label=None):
    handles = []
    # Plot events on map
    hdl = map_axis.scatter(df.lon, df.lat,
                      marker=marker,
                      c=color,
                      s=mutil.magscale(df.mag),
                      alpha=alpha,
                      zorder=zorder,
                      transform=ccrs.PlateCarree())
    handles.append(hdl)
    # Plot events on time-depth chart
    hdl = tz_axis.scatter(df.prefor_time, df.depth*1e-3,
                      c=color,
                      s=mutil.magscale(df.mag),
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

# # Load template profile
df_tpk = pd.read_csv(TPK)

### CREATE SETS FOR QUICK INDEXING AGAINST tidy_group
# Hard-set group set interpretations (from visual analysis!)
superset = set(df_cat.tidy_group.unique())
# Ungrouped
UG_set = set([-1])
# Surface Event Superset (to be converted to SU)
to_SU_set = set([1])
# Groups judged to warrant conversion to all-EQ
to_EQ_set = (set([26]))
# Native all-(deep)-LF groups
is_LF_set = set([6, 8, 38, 51])
# Groups judged to warrant conversion to all-(deep)-LF
to_LF_set = set([20, 63])
# Native all-probable blast groups
is_PX_set = set([57])
# Groups judged to warrant conversion to all-PX
to_PX_set = set([52])

PX_set = is_PX_set.union(to_PX_set)
SU_set = to_SU_set
LF_set = is_LF_set.union(to_LF_set)
EQ_set = superset.difference(PX_set, SU_set, LF_set, UG_set)



### MAPS ###
###############
#### SETUP ####
###############
# Initialize Basemap
fig = plt.figure(figsize=(15.12,7.4))
gs = fig.add_gridspec(ncols=2, nrows=2, hspace=0, wspace=0)

axm, attr = mutil.mount_baker_basemap(
    fig=fig,sps=gs[:, 0], radius_km=31,
    latnudge=0, lonnudge=0,
    open_street_map=False,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.1})

# Add Borders
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

## SUBPLOT 2
axt = fig.add_subplot(gs[0, 1])
axt.yaxis.set_ticks_position('right')
axt.yaxis.set_label_position('right')
axt.set_ylim([45, -3])
axt.set_xlim([pd.Timestamp('1978-01-01'), pd.Timestamp('2026-01-01')])
axt.set_xlabel('Year')
axt.set_ylabel('Depth [km]')


#### TEMPORARY PLOTTED ELEMENTS ####

##############
#### DATA ####
##############
#### PNSN Catalog Summary Figure ####
zo = 5
hdls = []
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    handles = plotfun(axm, axt, _df, marker_map[_et], color_map[_et],
                      zorder=zo,label=f'{_et.upper()}: {len(_df)}')
    hdls += handles
    zo -= 1
lhdl = axt.legend(handles=hdls[::-1], loc='lower center', ncols=4)
hdls.append(lhdl)

# Save Figure
if issave:
    fig.savefig(str(SDIR/f'PNSN_current_catalog_{DPI:d}dpi.{FMT}'),
                dpi=DPI, format=FMT)
# Clear temporary elements
clear_handles(hdls)

### PNSN Catalog Event Type Specific Catalogs
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    hdls = plotfun(axm, axt, _df, marker_map[_et], color_map[_et], 
                   zorder=zo,label=f'{_et.upper()}: {len(_df)}')
    lhdl = axt.legend(loc='lower center')
    hdls.append(lhdl)

    # Save Figure
    if issave:
        fig.savefig(str(SDIR/f'PNSN_current_catalog_{_et.upper()}_{DPI:d}dpi.{FMT}'),
                    dpi=DPI, format=FMT)
    # Clear temporary elements
    clear_handles(hdls)

#################
#### RESULTS ####
#################


### SHOW WHAT'S GROUPED ###
# Plot all groups
zo = df_cat.tidy_group.max() + 10
hdls = []

# Iterate across event types for marker shape
for _et in ['eq','lf','su','px']:
    # Iterate across groups for marker color
    for _grp in df_cat.tidy_group.value_counts(ascending=True).index:
        _df = df_cat[(df_cat.tidy_group==_grp) & (df_cat.etype==_et)]
        if len(_df) > 0:
            handles = plotfun(
                axm, axt, _df,
                marker_map[_et],
                _df.leaf_color[0],
                zorder=zo,
            )
            hdls += handles
            zo -= 1
        else:
            continue

# Save
savename = SDIR/f'map_clustered_events_{DPI}dpi.{FMT}'
if issave:
    fig.savefig(str(savename), dpi=DPI, format=FMT)
#Cleanup
hdls = clear_handles(hdls)

### SHOW WHAT'S MISSED ###
# Plot Only Ungrouped Events with catalog etype markings/colors
zo = 5
_df = df_cat[df_cat.tidy_group.isin(UG_set)]
# Iterate across etype
for _et in ['su','lf','px','eq']:
    __df = _df[_df.etype==_et] 
    handles = plotfun(
        axm, axt, __df, 
        marker_map[_et], 
        'xkcd:ugly brown', 
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
    fig.savefig(str(savename), dpi=DPI, format=FMT)
# Cleanup
hdls = clear_handles(hdls)


### PROPOSE CHANGES ### Include missed events of _ret?
# Plot each etype (re)classification subset
for _ret, _ids in zip(['su','lf','px','eq'],[SU_set, LF_set, PX_set, EQ_set]):
    # Plot ungrouped of target reassignment
    _df = df_cat[(df_cat.tidy_group == -1) & (df_cat.etype==_ret)]
    handles = plotfun(
        axm, axt, _df,
        marker_map[_ret],
        'xkcd:ugly brown',
        zorder=2,label=f'Ungrouped {_ret.upper()}: {len(_df)}'
    )
    hdls += handles
    _zo = 3
    # Iterate over each event type for shape and labeling
    for _et in ['eq','lf','su','px']:
        first = True
        _df = df_cat[(df_cat.tidy_group.isin(_ids))&(df_cat.etype==_et)]
    # Iterate over each group for color
        for _id in _ids:
            __df = _df[_df.tidy_group==_id]
            if len(__df) > 0:
                if first:
                    label=f'{_et.upper()}: {len(_df)}'
                    first = False
                else:
                    label = None
                handles = plotfun(
                    axm, axt, __df,
                    marker_map[_et],
                    __df.leaf_color[0],
                    zorder=_zo,
                    label=label)
                _zo += 1
                hdls += handles
    lhdl = axt.legend(loc='lower center', ncols=5)
    hdls.append(lhdl)
    # Save
    savename = SDIR/f'map_{_ret}_proposed_relabeling_{DPI}dpi.{FMT}'
    if issave:
        fig.savefig(str(savename), dpi=DPI, format=FMT)
    # Cleanup
    hdls = clear_handles(hdls)


### ENACT CHANGES ###
# Plot each etype (re)classification subset
for _ret, _ids in zip(['su','lf','px','eq'],[SU_set, LF_set, PX_set, EQ_set]):
    # Plot ungrouped of target reassignment
    _df = df_cat[(df_cat.tidy_group == -1) & (df_cat.etype==_ret)]
    handles = plotfun(
        axm, axt, _df,
        marker_map[_ret],
        'xkcd:ugly brown',
        zorder=2,label=f'Ungrouped {_ret.upper()}: {len(_df)}'
    )
    hdls += handles
    _zo = 3
    # Iterate over each event type for shape and labeling
    for _et in ['eq','lf','su','px']:
        first = True
        _df = df_cat[(df_cat.tidy_group.isin(_ids))&(df_cat.etype==_et)]
    # Iterate over each group for color
        for _id in _ids:
            __df = _df[_df.tidy_group==_id]
            if len(__df) > 0:
                if first:
                    label=f'{_et.upper()}: {len(_df)}'
                    first = False
                else:
                    label = None
                handles = plotfun(
                    axm, axt, __df,
                    marker_map[_ret],
                    __df.leaf_color[0],
                    zorder=_zo,
                    label=label)
                _zo += 1
                hdls += handles
    lhdl = axt.legend(loc='lower center', ncols=5)
    hdls.append(lhdl)
    
#     # Plot each etype (re)classification subset
# for _ret, _ids in zip(['su','lf','px','eq'],[SU_set, LF_set, PX_set, EQ_set]):
#     # Plot ungrouped of target reassignment
#     _df = df_cat[(df_cat.tidy_group == -1) & (df_cat.etype==_ret)]
#     handles = plotfun(
#         axm, axt, _df,
#         marker_map[_ret],
#         'xkcd:ugly brown',
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
    savename = SDIR/f'map_{_ret}_applied_relabeling_{DPI}dpi.{FMT}'
    if issave:
        fig.savefig(str(savename), dpi=DPI, format=FMT)
    hdls = clear_handles(hdls)


## Create Movie
    
# Assign new EVIDs to events


T0 = pd.Timestamp('1980-01-01')
T1 = pd.Timestamp('2025-03-01')
DTWIN = pd.Timedelta(180, unit='d')
DTSLIDE = pd.Timedelta(7, unit='d')

t0 = T0 - DTWIN
t1 = T0
total_frames = int(np.ceil((T1 - T0)/DTSLIDE))
frame_no = 0
# March forward in time
while t1 <= T1:
    print(f'frame {frame_no} of {total_frames:d}')
    # Plot slider on time figure
    shdl = axt.fill_between([t0, t1], [50,50], [-5,-5], color='cyan', alpha=0.25)
    hdls.append(shdl)
    axt.set_ylim([45, -3])

    # Subset all events
    _df = df_cat[(df_cat.prefor_time > t0) & (df_cat.prefor_time <= t1)]
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
                'xkcd:ugly brown',
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
        fig.savefig(str(savename), dpi=MDPI, format=FMT)
    # Advance timer
    t0 += DTSLIDE
    t1 += DTSLIDE

    # Advance frame number
    frame_no += 1
    # Clear temporary items
    hdls = clear_handles(hdls)

