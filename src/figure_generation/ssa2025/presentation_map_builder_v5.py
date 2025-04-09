from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

import cartopy.crs as ccrs

import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
CPD = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
LMD = ROOT/'results'/'tables'/'SSA2025'/'linkmat_average_cct0.300.npy'
TPK = ROOT/'results'/'tables'/'SSA2025'/'template_profile.csv'
# Load Catalog Profile
df_cat = pd.read_csv(CPD, parse_dates=['prefor_time'], index_col=[0])
# Subset to 30 km radius
df_cat = df_cat[df_cat.CAT0]
# Load precalculated linkage matrix
linkmat = np.load(LMD)
# Load template profile
df_tpk = pd.read_csv(TPK)

marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
color_map = {'px':'m','su':'b','lf':'r','eq':'k'}

def plotfun(map_axis, tz_axis, df, marker, color, alpha=0.6, zorder=1, label=None):
    handles = []
    # Plot events on map
    hdl = axm.scatter(df.lon, df.lat,
                      marker=marker,
                      c=color,
                      s=mutil.magscale(df.mag),
                      alpha=alpha,
                      zorder=zorder,
                      transform=ccrs.PlateCarree())
    handles.append(hdl)
    # Plot events on time-depth chart
    hdl = axt.scatter(df.prefor_time, df.depth*1e-3,
                      c=color,
                      s=mutil.magscale(df.mag),
                      marker=marker,
                      alpha=alpha,
                      zorder=zorder,
                      label=label)
    handles.append(hdl)
    return handles


### PLOTTING SECTION ###

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

# Plot total catalog map
zo = 5
hdls = []
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    handles = plotfun(axm, axt, _df, marker_map[_et], color_map[_et], zorder=zo,label=_et.upper())
    hdls += handles
    zo -= 1
axt.legend(handles=hdls[::-1], loc='lower center', ncols=4)

# Save

# Clear Axes
for _h in hdls:
    try:
        _h.remove()
    except:
        list(_h)[0].remove()

# Plot individual event type catalog maps
for _et in ['px','su','lf','eq']:
    _df = df_cat[df_cat.etype==_et]
    hdls = plotfun(axm, axt, _df, marker_map[_et], color_map[_et], zorder=zo,label=_et.upper())
    axt.legend(loc='lower center')

    # Save

    # Clear Axes
    for _h in hdls:
        _h.remove()
    

# Hard-set group set interpretations (from visual analysis!)
superset = set(df_cat.tidy_group.unique())
SU_set = set([1])
UG_set = set([-1])
LF_set = set([6, 8, 20, 38, 51, 63])
PX_set = set([59])

EQ_set = superset.difference(PX_set, SU_set, UG_set)

# Plot ungrouped events
zo = 5
hdls = []
for _et in ['su','lf','px','eq']:
    _df = df_cat[(df_cat.tidy_group==-1) & (df_cat.etype==_et) & (df_cat.group.notna())]    
    handles = plotfun(axm, axt, _df, marker_map[_et], color_map[_et], zorder=zo, label=f'{_et.upper()}: {len(_df):d}')
    hdls += handles
    zo -= 1
axt.legend(handles=hdls[::-1], loc='lower center', ncols=4)
breakpoint()


for _grps in [[1], [6, 8, 15, 48], [59], ]
