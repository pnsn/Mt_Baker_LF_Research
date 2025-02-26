import os, random
from pathlib import Path

import pandas as pd

from obsplus import EventBank

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from sklearn.cluster import AgglomerativeClustering

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
# Coherence distance table
COHD = ROOT / 'results' / 'tables' / 'coherence_distance_table.csv'



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
fig = plt.figure(figsize=(10,10))
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
set_link_color_palette(colors)
gs = fig.add_gridspec(ncols=1, nrows=1)
axes = {'dend': fig.add_subplot(gs[0,:])}
labels = df_eb.loc[dmerge.index][['etype','depth']]
labels = [f'{x.depth*1e-3:.0f} | {x.etype}' for _, x in labels.iterrows()]
# labels = dmerge.index
dend_out = dendrogram(linkmat,ax=axes['dend'], color_threshold=1-CCT,
                      labels=labels , distance_sort=True, above_threshold_color='gray')#,
                    #   no_labels=True)

## GROUPS MAP
fig = plt.figure(figsize=(8,8))
gs = fig.add_gridspec(ncols=2, nrows=2, wspace=0, hspace=0)
axm, mapattr = mount_baker_basemap(fig=fig, sps=gs[0,0],
                          open_street_map=False, radius_km=33.)
add_rings(axm, rads_km = [10,20,30], rads_colors=['k']*3)
# Get df_eb clusters in leaf sort order
_df = df_eb.loc[dmerge.index].iloc[dend_out['leaves']]
# Scatter with colors
axm.scatter(_df.longitude, _df.latitude, c=dend_out['leaves_color_list'], s=9, alpha=0.667, transform=ccrs.PlateCarree())

axe = fig.add_subplot(gs[0,1])
axe.scatter(_df.latitude, _df.depth*-1e-3, c=dend_out['leaves_color_list'], s=9, alpha=0.667)


axd = fig.add_subplot(gs[1,0])
axd.scatter(_df.longitude, _df.depth*-1e-3, c=dend_out['leaves_color_list'], s=9, alpha=0.667)

axz = fig.add_subplot(gs[1,1])
axz.scatter(_df.time, _df.depth*-1e-3, c=dend_out['leaves_color_list'], s=9, alpha=0.667)

emark = {'eq': 'ro', 'su': 'c^', 'lf': 'bs', 'px': 'm*'}

for _gn, _ct in _df.aggc.value_counts().items():
    if _ct > 5:
        fig = plt.figure(figsize=(8,8))
        gs = fig.add_gridspec(ncols=2, nrows=4, wspace=0, hspace=0)
        _idf = _df[_df.aggc == _gn]
        _xdf = _df[_df.aggc != _gn]
        axm, mapattr = mount_baker_basemap(fig=fig, sps=gs[:2,0],
                                open_street_map=False, radius_km=33.)
        axz = fig.add_subplot(gs[2,0])
        # axd = fig.add_subplot(gs[0,1])
        axt = fig.add_subplot(gs[3,0])

        axm.plot(_xdf.longitude, _xdf.latitude, 'xk', ms=2, alpha=0.15, transform=ccrs.PlateCarree())
        axz.plot(_xdf.longitude, _xdf.depth*-1e-3, 'xk', ms=2, alpha=0.15)
        axt.plot(_xdf.time, _xdf.depth*-1e-3, 'xk', ms=2, alpha=0.15)
        add_rings(axm, rads_km = [10,20,30], rads_colors=['k']*3)
        for _etype in _idf.etype.unique():
            __idf = _idf[_idf.etype == _etype]
            axm.plot(__idf.longitude, __idf.latitude, emark[_etype], ms=4, alpha=0.85, transform=ccrs.PlateCarree())
            axz.plot(__idf.longitude, __idf.depth*-1e-3, emark[_etype], ms=4, alpha=0.85)
            # axd.plot(__idf.latitude, __idf.depth*-1e-3, emark[_etype], ms=4, alpha=0.85)
            axt.plot(__idf.time, __idf.depth*-1e-3, emark[_etype], ms=4, alpha=0.85)

    #TODO: swap out one or two of these plots for waveform plots
            




plt.show()


