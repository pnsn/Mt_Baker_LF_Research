"""
TODO: Use this to plot the main dendrogram
TODO: Incorporate sidecar showing ungrouped events of etype, in position - use ugly brown


"""
import os
from pathlib import Path
from glob import glob

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from obspy.clients.fdsn import Client

from eqcorrscan import Tribe
from eqcorrscan.utils.stacking import align_traces

import cartopy.crs as ccrs

import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
# Absolute paths to input data files
TDIR = ROOT/'results'/'tables'/'SSA2025'
CPD = TDIR/'catalog_profile.csv'
TPK = TDIR/'template_profile.csv'


# Save Directory for Figures
SDIR = ROOT/'results'/'figures'/'SSA2025'

issave = True
isshow = True
DPI = 120
FMT = 'png'
plt.rcParams.update({'font.size': 14})
sfkw = {'dpi': DPI, 'format': FMT, 'bbox_inches':'tight'}
plt.rcParams.update({'font.size': 18})


def load_traces_and_cross_corr(_df_tpk, su_include):
    _df_tpks = _df_tpk[su_include]


    trace_list = []
    shifts = []
    corrs = []
    for _grp in _df_tpks.tidy_group.unique():
        _df = _df_tpks[_df_tpks.tidy_group==_grp]
        _trace_list = []
        for _f in _df.file:
            tribe = Tribe().read(_f)
            tr = tribe[0].st[0]
            _trace_list.append(tr)
        # If correlation-based group
        if _grp > 0:
            # Run in-group correlationshifting
            _shifts, _corrs = align_traces(_trace_list, shift_len=int(5*50), master=False)
        # Otherwise provide NaN
        else:
            _shifts = [np.nan]*len(_df)
            _corrs = [np.nan]*len(_shifts)
        shifts += _shifts
        corrs += _corrs
        trace_list += _trace_list
    # Attach 
    _df_tpks = _df_tpks.assign(corr=corrs)
    _df_tpks = _df_tpks.assign(shifts=shifts)
    _df_tpks = _df_tpks.assign(trace=trace_list)
    return _df_tpks


### SUPPORTING FUNCTIONS ###

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

def plot_dendrogram(dend, df_cat, s=25):
    marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
    color_map = {'px':'m','su':'b','lf':'r','eq':'k'}
    fig = plt.figure(figsize=(15,7))
    ax = fig.add_subplot(111)
    # Plot dendrogram
    for ic, dc, cc in zip(dend['icoord'],dend['dcoord'],dend['color_list']):
        if cc == 'k':
            cc = 'xkcd:ugly brown'
        else:
            pass
        ax.plot(ic, dc, color=cc, linewidth=0.5, zorder=1)
    # terminate endpoints with event type markers
    for _et in ['eq','lf','su','px']:
        _df_cat = df_cat[df_cat.etype==_et]
        ax.scatter(_df_cat.leafpos*10 + 5, [0]*len(_df_cat),
                   c=color_map[_et],
                   marker=marker_map[_et],
                   s=s, alpha=0.6, zorder=2)
    
    ax.set_xlim([-15, df_cat.leafpos.max()*10 + 20])
    ax.set_ylim([-0.033, 1])
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[f'{_e/10:.1f}' for _e in range(10, -2, -2)])
    ax.set_ylabel('Average Coherence [ - ]')
    return (fig, ax)


def plot_group_waveforms(_df_cat, _df_tpk, dend):
    # TODO: pin group color marker to ends of each trace
    # TODO: Tidy up group splits in waveforms (see EQ)
    # TODO: Include really simple zone map + station location (more or less JCW/MBW)
    # Marker Rendering by Event Type
    marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
    color_map = {'px':'m','su':'b','lf':'r','eq':'k'}
    ungrouped_color = 'xkcd:ugly brown'
    unanalyzed_color = 'xkcd:gross green'

    # Initialize Figure
    fig = plt.figure(figsize=(15.12,7.4))
    gs = fig.add_gridspec(ncols=8,nrows=3, hspace=0)
    # Initialize Subplots
    axd = fig.add_subplot(gs[0,:6])
    axw = fig.add_subplot(gs[1:,:6])
    # axm, attr = mutil.mount_baker_basemap(
    #     fig=fig, sps=gs[:,6:],
    #     latnudge=0, lonnudge=0,
    #     open_street_map=False,
    #     aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})

    # extent = [mutil.BAKER_LON-0.5, mutil.BAKER_LON+0.5,
    #           mutil.BAKER_LAT-1, mutil.BAKER_LAT+0.4]
    # axm.set_extent(extent)
    # gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
    #                xlocs=[-122.1, -121.8, -121.5],
    #                ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
    #                alpha=0)
    # gl.left_labels=False
    # gl.right_labels=True
    # xlims = [float('inf'), 0]
    # _lm = 0
    # Plot base dendrogram for specific group
    for ic, dc, cc in zip(dend['icoord'],dend['dcoord'], dend['color_list']):

        if cc in _df_cat.leaf_color.values:
            axd.plot(ic, dc, color=cc, linewidth=0.5, zorder=2)
        else:
            axd.plot(ic, dc, color='xkcd:ugly brown', linewidth=0.5, zorder=1)
    # Plot etype markers
    for evid, row in _df_cat.iterrows():
        if evid in _df_tpk.index:
            _y = 0
        else:
            _y = 0.08
        axd.scatter(row.leafpos*10 + 5, _y,
                    s=25, 
                    c=color_map[row.etype], 
                    marker=marker_map[row.etype],
                    zorder=3,
                    alpha=0.6)
        # if _y == 0:
        #     _xrng = _df_cat.leafpos.max() - _df_cat.leafpos.min() - 6.5
        #     _dx = _xrng/len(_df_tpk)
        #     _x0 = _df_cat.leafpos.min() + 5

        #     axd.plot([row.leafpos*10 + 5, 
        #               (_x0 + _lm*_dx)*10 + 5],
        #               [_y, -0.03],
        #               color=color_map[row.etype],
        #               zorder=2, linewidth=0.33,
        #               alpha=0.6)
        #     _lm += 1
    xlims = list(axd.get_xlim())
    xlims[0] = 5 + 10*_df_cat.leafpos.min() - 15
    xlims[1] = 5 + 10*_df_cat.leafpos.max() + 15
    axd.set_xlim(xlims)
    axd.set_ylim([-0.033,1])
    axd.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=[f'{_e/10:.1f}' for _e in range(10, -2, -2)])
    axd.set_ylabel('Average\nCoherence [ - ]')
    axd.set_xticks([])

    # Plot waveforms with shifts descending
    _dx = 0
    for _e, (_, row) in enumerate(_df_tpk.iterrows()):
        _y = row.trace.times() + row.shifts - 5
        _x = row.trace.copy().normalize().data
        axw.plot(_x + _dx, _y, color=color_map[row.etype], linewidth=0.5, zorder=1)
        # Plot the color of this WF's group, but original shape
        axw.scatter(_x[0] + _dx, -0.5, s=25, c=row.leaf_color,
                    marker=marker_map[row.etype],
                    zorder=2,
                    alpha=0.6)
        # Increment to next position
        _dx += 1
        # If there are more traces to plot
        if _e < len(_df_tpk) - 1:
            # If the next trace to plot is in another group
            if row.tidy_group != _df_tpk.tidy_group.iloc[_e+1]:
                # Plot a dashed vertical line
                axw.plot([_dx, _dx],[-2,45],':', color='xkcd:ugly brown', alpha=0.5, linewidth=0.75)
                # Insert Padding
                _dx += 1
                # And
        # # If there are more samples
        # if _k < len(_df_tpk) - 1:
        #     # If the next sample has a different etype
        #     if row.tidy_group != _df_tpk.tidy_group.iloc[_k+1]:
        #         # Add padding between traces
        #         _e += 2
    # Format axes
    axw.set_ylim([45, -2])
    axw.set_xlim([-2, _dx + 2])
    axw.set_xticks([])
    axw.set_xlabel('Less Coherent        Relative Coherence        More Coherent')
    axw.set_ylabel('Elapsed Time Since\nP-wave Arrival [sec]')



    # # Add distance Rings
    # mutil.add_rings(axm,
    #                 rads_km=[10,30,50,90],
    #                 rads_colors=['k']*4,
    #                 include_units=[True, True, True, True],
    #                 label_pt=-50,ha='left',va='bottom')

    # # Add Lat/Lon annotations
    # # axm.set_xticks([-122.2 + _e*0.2 for _e in range(5)])
    # gl = axm.gridlines(draw_labels=['top','left'], zorder=1,
    #                 xlocs=[-122.1, -121.8, -121.5],
    #                 ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
    #                 alpha=0)
    #                 #    xlocs=[-122.2, -122, -121.8, -121.6, -121.4, -121.2],
    #                 #    ylocs=[48.5, 48.6, 48.7, 48.8, 48.9, 49],
    #                 #    alpha=0)


    # # Plot Mount Baker
    # axm.scatter([mutil.BAKER_LON, mutil.SHUKS_LON, mutil.STWIN_LON],
    #         [mutil.BAKER_LAT, mutil.SHUKS_LAT, mutil.STWIN_LAT],
    #         s=[144, 81, 81], 
    #         marker='^',
    #         facecolor='none',
    #         edgecolors='orange',
    #         linewidths=[2, 1, 1],
    #         zorder=30,
    #         transform=ccrs.PlateCarree())

    # # Populate Map
    # # Plot observing station
    # inv = Client('IRIS').get_stations(network=_df_tpk.net.iloc[0],
    #                                   station=_df_tpk.sta.iloc[0],
    #                                   level='station')
    # ista = inv.networks[0].stations[0]
    # axm.scatter(
    #     ista.longitude, ista.latitude,
    #     marker='v', c='w',
    #     edgecolors='k',
    #     s=36, zorder=10,
    #     transform=ccrs.PlateCarree()
    # )
    # axm.text(
    #     ista.longitude, ista.latitude,
    #     f'UW.{ista.code}',fontsize=14,
    #     ha='left',va='bottom',
    #     transform=ccrs.PlateCarree()
    # )
    # # Iterate across source ETYPE
    # for _et in _df_tpk.etype.unique():
    #     # Subset & plot using group color
    #     _df = _df_tpk[_df_tpk.etype==_et]
    #     _ = plotmap(axm, _df,
    #                 marker_map[_et],
    #                 _df.leaf_color)

    return (fig, axd, axw) #, axm)






#### LOAD TABLES ####
df_cat = pd.read_csv(CPD, parse_dates=['prefor_time'], index_col=[0])
df_tpk = pd.read_csv(TPK)

#### LOAD DENDROGAM STRUCTURES ####
dend = dict()
for _f in glob(str(TDIR/'*_.npy')):
    arr = np.load(_f)
    fname = os.path.splitext(os.path.split(_f)[-1])[0]
    key = '_'.join(fname.split('_')[1:])[:-1]
    dend.update({key:arr}) 



### PROCESS SU SUPER GROUP
# Subset to group 1 & sort by increasing leaf position
_df_cat_su = df_cat[df_cat.tidy_group == 1].sort_values('leafpos')

# Get template entries matching desired events and station
_df_tpk_su = df_tpk[(df_tpk.evid.isin(_df_cat_su.index))&
                 (df_tpk.pref_nslc)&
                 (df_tpk.sta=='SHUK')]

# Attach _df_tpk_su to _df_cat
_df_tpk_su.index = _df_tpk_su.evid
_df_tpk_su = _df_cat_su.join(_df_tpk_su, how='right')
# Sort by leafpos
_df_tpk_su = _df_tpk_su.sort_values('leafpos')

# Get subsample of the superset
su_include = []
for _e, (_, row) in enumerate(_df_tpk_su.iterrows()):
    # Inclued all PX and EQ to be reassigned
    if row.etype in ['px','eq','su']:
        su_include.append(True)
    # Include the "worst" 5 matches included
    elif _e <= 5:
        su_include.append(True)
    # Include the "best" 5 matches included
    elif _e >= len(_df_tpk_su) - 6:
        su_include.append(True)
    elif _e%7 == 0:
        su_include.append(True)
    else:
        su_include.append(False)



### PROCESS LF GROUPS ###
# Do things that will be mapped to LF
_df_cat_lf = df_cat[(df_cat.petype=='lf')&(df_cat.tidy_group > 0)]
# View from JCW - sees all
_df_tpk_lf = df_tpk[(df_tpk.evid.isin(_df_cat_lf.index))&
                    (df_tpk.pref_nslc)&
                    (df_tpk.sta=='JCW')]
                    # (df_tpk.sta.isin(['MBW','MBW2']))]
# Attach _df_tpk_lf to _df_cat
_df_tpk_lf.index = _df_tpk_lf.evid
_df_tpk_lf = _df_cat_lf.join(_df_tpk_lf, how='right')
# Sort by leafpos
_df_tpk_lf = _df_tpk_lf.sort_values('leafpos')

lf_include = [True]*len(_df_tpk_lf)




### PROCESS EQ GROUPS ###
_df_cat_eq = df_cat[(df_cat.petype=='eq')&(df_cat.tidy_group > 0)]
# Remove Group 18 due to bad signal on JCW
_df_cat_eq = _df_cat_eq[(df_cat.tidy_group != 18)]
# Limit to groups with 6+ events
eq_grp_cts = _df_cat_eq.tidy_group.value_counts()
eq_cids = eq_grp_cts[eq_grp_cts > 2].index.values
_df_cat_eq = _df_cat_eq[_df_cat_eq.tidy_group.isin(eq_cids)]
# View from JCW
_df_tpk_eq = df_tpk[(df_tpk.evid.isin(_df_cat_eq.index))&
                    (df_tpk.pref_nslc)&
                    (df_tpk.sta=='JCW')]
# Join CAT & TPK
_df_tpk_eq.index = _df_tpk_eq.evid
_df_tpk_eq = _df_cat_eq.join(_df_tpk_eq, how='right')
# Sort by leafpos
_df_tpk_eq = _df_tpk_eq.sort_values('leafpos')


eq_include = []
for _grp in _df_tpk_eq.tidy_group.unique():
    _df = _df_tpk_eq[_df_tpk_eq.tidy_group==_grp]
    for _e, (_, row) in enumerate(_df.iterrows()):
        # Include anything being relabeled
        if row.etype != 'eq':
            eq_include.append(True)
        # Include "worst" of each group
        elif _e == 0:
            eq_include.append(True)
        # Include "best" of each group
        elif _e == len(_df) - 1:
            eq_include.append(True)
        # Otherwise, include every 3th
        elif _e%3 == 0:
            eq_include.append(True)
        else:
            eq_include.append(False)
        
### PROCESS PX GROUPS ###
# Do things that will be mapped to LF
_df_cat_px = df_cat[(df_cat.petype=='px')&(df_cat.tidy_group > 0)]
# View from JCW - sees all
_df_tpk_px = df_tpk[(df_tpk.evid.isin(_df_cat_px.index))&
                    (df_tpk.pref_nslc)&
                    (df_tpk.sta=='MBW')&
                    (df_tpk.revstat=='manual')]
                    # (df_tpk.sta.isin(['MBW','MBW2']))]
# Attach _df_tpk_px to _df_cat
_df_tpk_px.index = _df_tpk_px.evid
_df_tpk_px = _df_cat_px.join(_df_tpk_px, how='right')
# Sort by leafpos
_df_tpk_px = _df_tpk_px.sort_values('leafpos')

px_include = [True]*len(_df_tpk_px)





### TRACE LOADING AND ALIGNMENT SECTION ###
_df_tpk_su = load_traces_and_cross_corr(_df_tpk_su, su_include)
_df_tpk_lf = load_traces_and_cross_corr(_df_tpk_lf, lf_include)
_df_tpk_eq = load_traces_and_cross_corr(_df_tpk_eq, eq_include)
_df_tpk_px = load_traces_and_cross_corr(_df_tpk_px, px_include)

### PLOTTING SECTION ###

### PLOT FULL DENDROGRAM WITH CATALOG EVENT LABELS
dend_handles = plot_dendrogram(dend, df_cat[df_cat.tidy_group > -2])
ax = dend_handles[1]
ax.set_xticks([])
xlims = ax.get_xlim()
ax.plot(ax.get_xlim(), [0.7, 0.7], ':', color='dodgerblue')
ax.set_xlim(xlims)

if issave:
    dend_handles[0].savefig(str(SDIR/f'full_dendrogram_{DPI}dpi.{FMT}'), **sfkw)

### PLOT RECLASSIFICATION FIGURE FOR SU
su_fig_axes = plot_group_waveforms(_df_cat_su, _df_tpk_su, dend)
if issave:
    su_fig_axes[0].savefig(str(SDIR/f'su_reclassified_cluster_waveforms_{DPI}dpi.{FMT}'), **sfkw)

### PLOT RECLASSIFICATION FIGURE FOR LF
lf_fig_axes = plot_group_waveforms(_df_cat_lf, _df_tpk_lf, dend)
if issave:
    lf_fig_axes[0].savefig(str(SDIR/f'lf_reclassified_cluster_waveforms_{DPI}dpi.{FMT}'), **sfkw)

### PLOT RECLASSIFICATION FIGURE FOR EQ
eq_fig_axes = plot_group_waveforms(_df_cat_eq, _df_tpk_eq, dend)
if issave:
    eq_fig_axes[0].savefig(str(SDIR/f'eq_reclassified_cluster_waveforms_{DPI}dpi.{FMT}'), **sfkw)

### PLOT RECLASSIFICATION FIGURE FOR PX
px_fig_axes = plot_group_waveforms(_df_cat_px, _df_tpk_px, dend)
if issave:
    px_fig_axes[0].savefig(str(SDIR/f'px_reclassified_cluster_waveforms_{DPI}dpi.{FMT}'), **sfkw)


if isshow:
    plt.show()
