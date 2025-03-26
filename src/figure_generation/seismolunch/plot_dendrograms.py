import os
from collections import defaultdict
from glob import glob
from pathlib import Path

import pandas as pd
from obspy.clients.fdsn import Client
from obsplus import EventBank
from eqcorrscan import Tribe
from eqcorrscan.utils.stacking import align_traces
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from sklearn.cluster import AgglomerativeClustering

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
# Coherence distance table
COHD = ROOT / 'results' / 'tables' / 'coherence_distance_table.csv'
# Template directory
TMPD = ROOT / 'processed_data' / 'template' / 'single_station'


def get_symmetric(df, i_field='event_i', j_field='event_j', k_field='coh', trace_value=1., aggfunc='mean'):
    """Get a symmetric matrix from a sparse representation of the upper triangle of 
    said matrix with positions designated by i_field, j_field and values designated
    by k_field. The trace is uniformly populated with trace_value and repeated (i,j,k)
    values are combined with the specified aggfunc (see pandas.pivot_table)

    :param df: _description_
    :type df: _type_
    :param i_field: _description_, defaults to 'event_i'
    :type i_field: str, optional
    :param j_field: _description_, defaults to 'event_j'
    :type j_field: str, optional
    :param k_field: _description_, defaults to 'coh'
    :type k_field: str, optional
    :param trace_value: _description_, defaults to 1.
    :type trace_value: _type_, optional
    :param aggfunc: _description_, defaults to 'mean'
    :type aggfunc: str, optional
    :return: _description_
    :rtype: _type_
    """    
    # Get all event names
    fullset = sorted(set(df[i_field]).union(set(df[j_field])))
    # Create pivot table of upper triangle
    cov = df.pivot_table(index=i_field, columns=j_field, values=k_field, aggfunc=aggfunc)
    # Pad out missing columns & rows
    cov = cov.reindex(index=fullset, columns=fullset, fill_value=np.nan)
    # Fill in lower triangle
    cov = cov.combine_first(cov.T)
    # Fill in trace
    np.fill_diagonal(cov.values, trace_value)
    return cov

def join_cov_df(df1, df2, aggfunc=np.nanmean, fill_value=np.nan):
    """Combine two symmetric, labeled arrays

    :param df1: _description_
    :type df1: _type_
    :param df2: _description_
    :type df2: _type_
    :param aggfunc: _description_, defaults to np.nanmean
    :type aggfunc: _type_, optional
    :param fill_value: _description_, defaults to np.nan
    :type fill_value: _type_, optional
    :return: _description_
    :rtype: _type_
    """    
    fullset = sorted(set(df1.index).union(df2.index))
    cov1_part = df1.reindex(index=fullset, columns=fullset, fill_value=fill_value)
    cov2_part = df2.reindex(index=fullset, columns=fullset, fill_value=fill_value)
    # Create 3-D array
    covstack = np.stack([cov1_part.values, cov2_part.values], axis=0)
    # Apply aggfunc across stack index
    covjoin = aggfunc(covstack, axis=0)
    # Convert back to dataframe
    cov_joined = pd.DataFrame(covjoin, index=fullset, columns=fullset)
    # Ensure symmetry
    cov_joined = cov_joined.combine_first(cov_joined.T)
    return cov_joined


CCT = 0.485

# SAVEPATH
SAVEPATH = ROOT / 'results' / 'figures' / 'seismolunch'
FMT = 'png'
DPI = 200

### PROCESSING SECTION ###
client = Client('IRIS')
# Iterate across templates to confirm if they had waveform data
tlist = glob(str(TMPD/'*'/'*'/'*'/'*.tgz'))
evid_sets = {'automatic': set(), 'manual': set()}
holder = []
rmap = {'manual': 1, 'automatic': -1}
for _f in tlist:
    parts = _f.split('/')
    evid = parts[-1].split('.')[0]
    year = int(parts[-2])
    revstat = parts[-3]
    stachan = parts[-4]
    evid_sets[revstat].add(evid)
    line = stachan.split('.') + [year, revstat, evid, rmap[revstat], _f]
    holder.append(line)

df_tpk = pd.DataFrame(holder, columns=['net','sta','loc','chan','year','revstat','evid','ohe','file'])

# Read Coherence Sparse Matrix
df_coh = pd.read_csv(COHD)

# Process individual coherence and shift matrices
coh_dict = {}
shift_dict = {}
for _k in df_coh.trace.unique():
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)

# Merge for all templates
coh_merge = pd.DataFrame()
for _v in coh_dict.values():
    coh_merge = join_cov_df(coh_merge, _v)

# Fill NaN entries with 0
coh_merge = coh_merge.replace(to_replace={np.nan: 0})

dmerge = 1. - coh_merge
# Run agglomerative clustering
model = AgglomerativeClustering(
    linkage='single',
    distance_threshold=1-CCT,
    n_clusters=None,
    metric='precomputed').fit(dmerge.values)

# Form Linkage matrix
counts = np.zeros(model.children_.shape[0])
n_samples = len(model.labels_)
for _e, merge in enumerate(model.children_):
    cc = 0
    for child_idx in merge:
        if child_idx < n_samples:
            cc += 1
        else:
            cc += counts[child_idx - n_samples]
    counts[_e] = cc
linkmat = np.column_stack([model.children_, model.distances_, counts]).astype(float)


EBANK = EventBank(EBBP)
df_eb = EBANK.read_index()
df_pp = pd.read_csv(CATD)

# Merge with catalog membership
df_eb.index = df_eb.event_id
df_eb = df_eb.join(pd.read_csv(CATD, index_col='event_id'), how='left')

# Update index for merging with uw####### evid format
df_eb.index = [f'{x.split("/")[-2].lower()}{x.split("/")[-1]}' for x in df_eb.index]

df_eb = pd.concat([df_eb, pd.Series(model.labels_, index=dmerge.index, name='aggc')], axis=1, ignore_index=False)

ngrps = sum(df_eb.aggc.value_counts() > 1)




## DENDROGRAM ALONE

# dcmap = plt.cm.get_cmap('nipy_spectral_r', ngrps)
# # Stagger colors
# colors = []
# _ee = 0
# for _e in range(ngrps):
#     _ee += 3
#     if _ee > ngrps - 1:
#         _ee -= ngrps
#     colors.append(mcolors.to_hex(dcmap(_ee)))

colors = [
    '#00a400',
    '#8eff00',
    '#00a2c1',
    '#00a8af',
    '#ffb200',
    '#00f600',
    '#f70000',
    '#830094',
    '#52005e',
    '#780089',
    '#cc3f3f',
    '#df0000',
    '#00de00',
    '#f6de00',
    '#00c700',
    '#7d008e',
    '#007ddd',
    '#00aaa2',
    '#ff4500',
    '#0062dd',
    '#0000af',
    '#0094dd',
    '#0089dd',
    '#cccccc',
    '#009a05',
    '#d2f700',
    '#0039dd',
    '#0000c1',
    '#00aa96',
    '#f0e900',
    '#0000d2',
    '#cd0000',
    '#00ea00',
    '#ff7900',
    '#00bb00',
    '#4dff00',
    '#29002f',
    '#fcd200',
    '#00aa8a',
    '#7f009a',
    '#0010dd',
    '#d80000',
    '#009dd2',
    '#009f34',
    '#d20000',
    '#cc8686',
    '#00d200',
    '#0dff00',
    '#00af00',
    '#5000a0',
    '#000000',
    '#00a562',
    '#ff1000',
    '#ffc300',
    '#e3f200',
    '#2100a6',
    '#eb0000',
    '#ffa000',
    '#c0fd00'
    ]

# # colors = [mcolors.to_hex(dcmap(_e)) for _e in ]
# colors = [mcolors.to_hex(dcmap(_e)) for _e in range(ngrps)]
# random.shuffle(colors)
fig = plt.figure(figsize=(12,5))

set_link_color_palette(colors)
gs = fig.add_gridspec(ncols=1, nrows=1)
ax = fig.add_subplot(gs[0])
labels = df_eb.loc[dmerge.index][['etype','depth']]
labels = [f'{x.depth*1e-3:.0f} | {x.etype}' for _, x in labels.iterrows()]
# labels = dmerge.index
dend_out = dendrogram(linkmat,ax=ax, color_threshold=1-CCT,
                      labels=labels , distance_sort=True, above_threshold_color='gray')#,
                    #   no_labels=True)
ax.set_ylabel('Correlation Distance (1 - mean correlation coefficient)')
ax.set_xlabel('Depths [km] | Event Types')




## GROUPS MAP
fig = plt.figure(figsize=(9.54,9.3))
gs = fig.add_gridspec(ncols=2, nrows=2, wspace=0, hspace=0)

## MAP PLOT
axm, mapattr = mount_baker_basemap(fig=fig, sps=gs[0,0],
                          open_street_map=False, radius_km=33.)
add_rings(axm, rads_km = [10,20,30], rads_colors=['k']*3, label_pt=10)
# Get df_eb clusters in leaf sort order
_df = df_eb.loc[dmerge.index].iloc[dend_out['leaves']]
# Attach color information
_df = _df.assign(ecolor=dend_out['leaves_color_list'])
# Scatter groups on map with colors
axm.scatter(_df.longitude, _df.latitude, c=_df.ecolor, s=9, alpha=0.667, transform=ccrs.PlateCarree())
# Add Gridlines
gl = axm.gridlines(draw_labels=True, zorder=1)
gl.bottom_labels = False
gl.right_labels = False
gl.xlines = False
gl.ylines = False

## Scatter lat/time in upper right
axe = fig.add_subplot(gs[0,1])
# axe.scatter(_df.depth*1e-3, _df.latitude, c=_df.ecolor, s=9, alpha=0.667)
axe.scatter(_df.time, _df.latitude, c=_df.ecolor, s=9, alpha=0.667)
# axe.set_xlim([-3, 40])
axe.yaxis.set_ticks_position('right')
axe.yaxis.set_label_position('right')
axe.xaxis.set_ticks_position('top')
axe.xaxis.set_label_position('top')
axe.set_xlabel('Year')
axe.set_ylabel('Latitude [$^o$E]', rotation=270, labelpad=15)
## Scatter  lon/depth in lower left
axd = fig.add_subplot(gs[1,0])
axd.scatter(_df.longitude, _df.depth*1e-3, c=_df.ecolor, s=9, alpha=0.667)
axd.set_ylim([40, -3])
axd.set_ylabel('Depth [km]')
axd.set_xlabel('Longitude [$^o$E]')

## Scatter time depth
axz = fig.add_subplot(gs[1,1])
axz.scatter(_df.time, _df.depth*1e-3, c=_df.ecolor, s=9, alpha=0.667)
axz.set_ylim([40, -3])
axz.yaxis.set_ticks_position('right')
axz.yaxis.set_label_position('right')
axz.set_ylabel('Depth [km]', rotation=270, labelpad=15)
axz.set_xlabel('Year')
# axz.axes.set(fontsize=8)
# axz.xaxis.set_ticks([pd.Timestamp(f'{x}-01-01') for x in [1980, 1990, 2000, 2010, 2020, 2030]])

# breakpoint()
# plt.show()

emark = {'eq': 'ro', 'su': 'kd', 'lf': 'bs', 'px': 'm*'}


## INDIVIDUAL GROUP COMPOUND PLOTS
for _e, (_gn, _ct) in enumerate(_df.aggc.value_counts().items()):
    if _ct < 5:
        continue
    fig = plt.figure(figsize=(12,9))
    gs = fig.add_gridspec(ncols=2, nrows=4)#, wspace=0, hspace=0)
    # In zone, in-group subset metadata
    _idf = _df[_df.aggc == _gn]
    # All other events not in this subset
    _xdf = _df[_df.aggc != _gn]

    # # Dendrogram subplot
    # axd = fig.add_subplot(gs[0,0])
    # _dend = dendrogram(linkmat, color_threshold=1 - CCT, ax=axd)


    # Map subplot
    axm, mapattr = mount_baker_basemap(fig=fig, sps=gs[:2,0],
                            open_street_map=False, radius_km=50.)
    # Longitude/depth subplot
    axz = fig.add_subplot(gs[0,1])
    # Timeline subplot
    axt = fig.add_subplot(gs[1,1])
    # Waveform subplot
    axwf = fig.add_subplot(gs[2:, :])

    axm.plot(_xdf.longitude, _xdf.latitude, 'xk', ms=2, alpha=0.15, transform=ccrs.PlateCarree())
    axz.plot(_xdf.longitude, _xdf.depth*1e-3, 'xk', ms=2, alpha=0.15)
    axt.plot(_xdf.time, _xdf.depth*1e-3, 'xk', ms=2, alpha=0.15)
    add_rings(axm, rads_km = [10,30,50], rads_colors=['k']*3)
    for _etype in _idf.etype.value_counts().index:
        __idf = _idf[_idf.etype == _etype]
        axm.plot(__idf.longitude, __idf.latitude, emark[_etype], ms=4, alpha=0.85,
                    transform=ccrs.PlateCarree(), label=f'{_etype.upper()} ({len(__idf)})')
        axz.plot(__idf.longitude, __idf.depth*1e-3, emark[_etype], ms=4, alpha=0.85)
        # axd.plot(__idf.latitude, __idf.depth*-1e-3, emark[_etype], ms=4, alpha=0.85)
        axt.plot(__idf.time, __idf.depth*1e-3, emark[_etype], ms=4, alpha=0.85)

    #TODO: swap out one or two of these plots for waveform plots
    evidset = set(_idf.index)
    trace_df = pd.DataFrame()
    index = []
    data = []
    _tpk = df_tpk[df_tpk.evid.isin(evidset)]
    rep_order = list(_tpk.sta.value_counts().index.values)
    __r = 0
    while evidset != set():
        # Get subset of template file metadata
        _tpk = df_tpk[df_tpk.evid.isin(evidset)]
        _to_load = _tpk[_tpk.sta == rep_order[__r]]
        # _topnetsta = '.'.join(_tpk[['net','sta']].value_counts().index[0])
        # rep_order.append(_topnetsta)
        # # Get file names for most present station
        # _to_load = _tpk[_tpk.sta == _topnetsta.split('.')[1]]
        # Load templates
        for _f in _to_load.file:
            tmp = Tribe().read(_f)[0]
            # Get trace
            tr = tmp.st[0]
            # attach trace to list
            # trace_lists[_topsta].append(tr)
            # Get EVID from template name
            _tevid = tmp.name
            index.append(_tevid)
            data.append([tr, tr.stats.network, tr.stats.station])

            # Remove EVID from evidset
            evidset.remove(_tevid)
        __r += 1
        if __r == len(rep_order):
            break
    # Append traces and net/sta metadata to temporary dataframe
    _idf = pd.concat([_idf, pd.DataFrame(data, index=index, columns=['trace', 'net','sta'])], axis=1, ignore_index=False)
    # breakpoint()
    # Iterate across topsta
    _offset = 2*len(_idf)
    _initial_offset = _offset
    ylims = [30, -5]
    for _sta in rep_order:
        # if len(_idf.sta.unique()) > 1:
        __idf = _idf[_idf.sta == _sta]
        tr_list = list(__idf.trace)
        if len(tr_list) < 2:
            continue
        shifts_corrs = align_traces(tr_list, 50)
        for _t, _tr in enumerate(tr_list):
            __tr = _tr.copy()
            __tr.detrend('linear')
            __tr.normalize()
            __tr.normalize(norm=__tr.std()*3)
            _ccsign = np.sign(shifts_corrs[1][_t])
            _yv = __tr.times() + shifts_corrs[0][_t] - tmp.prepick
            _xv = _ccsign*__tr.data/1. + _offset
            _offset -= 2
            # Plot waveform
            axwf.plot(_xv, _yv, color=emark[__idf.etype[_t]][0], lw=0.5)

        # axwf.fill_between(xlims, [_initial_offset]*2, [_offset]*2, alpha=0.25)
        _netsta = f'{_tr.stats.network}.{_tr.stats.station}'
        axwf.text(0.5*(_initial_offset + _offset),ylims[0] + 1, _netsta, ha='center', va='top')
        _initial_offset = _offset

        # Plot station location on map
        kw = dict(zip(['network','station','location','channel'], __tr.id.split('.')))
        kw.update({'level': 'station'})
        _inv = client.get_stations(**kw)
        for net in _inv.networks:
            for sta in net.stations:
                axm.plot(sta.longitude, sta.latitude, 'wv', mec='k', linewidth=0.5, ms=4, transform=ccrs.PlateCarree())
                axm.text(sta.longitude - 0.02, sta.latitude + 0.01, _netsta, transform=ccrs.PlateCarree(),
                         color='black', fontsize=8, ha='right', va='bottom')


    
    
    # Figure Formatting
    # Map
    axm.legend(loc='upper right', fontsize=10)
    gl = axm.gridlines(draw_labels=True, zorder=1)
    gl.bottom_labels = False
    gl.right_labels = False
    gl.xlines = False
    gl.ylines = False
    # Londepth
    axz.set_ylabel('Depth [km]')
    axz.set_ylim([40, -3])
    # Timedepth
    axt.set_ylim([40, -3])
    axt.set_ylabel('Depth [km]')
    axt.set_xlabel('Year')
    # axt.xaxis.set_ticks_position()
    # Waveforms
    axwf.set_ylim(ylims)
    axwf.xaxis.set_ticks([])
    axwf.set_xlim([0, 2*len(_idf) + 2])
    axwf.set_ylabel('Time Relative to Aligned Onset [sec]')
    # axwf.set_xlabel('')
    #     _df_tpk
    # _df_tpk = df_tpk[df_tpk.evid.isin(_idf.index)]
    # breakpoint()
    # _v1 = _df_tpk.sta.value_counts()
    # _to_load = _df_tpk[_df_tpk.sta==_v1.index[0]]
    # files_to_load = 
    # breakpoint()


plt.show()


