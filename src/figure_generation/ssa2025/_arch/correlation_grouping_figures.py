from glob import glob
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

import pandas as pd
from scipy.cluster.hierarchy import dendrogram, set_link_color_palette
from sklearn.cluster import AgglomerativeClustering

from obspy.clients.fdsn import Client
from obsplus import EventBank
from eqcorrscan import Tribe
from eqcorrscan.utils.stacking import align_traces

import cartopy.crs as ccrs
from map_util import *



# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to catalog membership CSV
CATP = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
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

####### PROCESSING SECTION ########


## COHERENCE MATRIX PROSTPROCESSING ##
# Scrape metadata from templates (without loading them) to inform waveform loading for group plots
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


## EVENT CATALOG METADATA LOAD ##
# Read catalog profile
df_cat = pd.read_csv(CATP, parse_dates=['prefor_time'])

# Update index for merging group labels with uw####### evid format
df_cat.index = [f'{x.split("/")[-2].lower()}{x.split("/")[-1]}' for x in df_cat.event_id]
# Conduct concatenation
df_cat = df_cat.join(pd.Series(model.labels_, index=dmerge.index, name='group'), how='right')
# Create "tidy group index" for display purposes
tidy_group = np.ones(len(df_cat))*-1
for _e, (_gn, _gc) in enumerate(df_cat.group.value_counts().items()):
    if _gc > 1:
        tidy_group[df_cat.group==_gn] = _e + 1

df_cat = df_cat.assign(tidy_group=tidy_group)

## STATION LOCATION FETCH ##
# Connect to client
client = Client('IRIS')
# Get station locations
inv = client.get_stations(network='UW,CN',station='MBW,MBW2,RPW,RPW2,JCW,SHUK,MULN,VDB,DEB,CMW', level='station')


### GENERATE DENDROGAM ELEMENTS WITH COLOR FORMATTING
ngrps = sum(df_cat.group.value_counts() > 1)
dcmap = plt.cm.get_cmap('jet', ngrps+4)

colors = [mcolors.to_hex(dcmap(_e)) for _e in range(ngrps)]
# Create dendrogram outputs & structure
dend = dendrogram(
    linkmat,
    labels=df_cat.index,
    distance_sort=False, count_sort='ascending',
    color_threshold=1-CCT, above_threshold_color='k',
    leaf_rotation=0, no_plot=True,
    get_leaves=True)

# Get leaf coordinates associated back to df_cat
df_leaf = pd.DataFrame(data={'leafpos': range(len(dend['leaves'])),'leaf': dend['leaves'],'leaf_color':dend['leaves_color_list']}, 
                    index=dend['ivl'])
df_cat = df_cat.join(df_leaf, how='left')

## PLOTTING SECTION ##


# random.shuffle(colors)
# 3 step shuffle
# colors = colors[::2] + colors[1::2][::-1]

fig = plt.figure(figsize=(12,5))

set_link_color_palette(colors)
gs = fig.add_gridspec(ncols=1, nrows=1)
ax = fig.add_subplot(gs[0])


                    #   no_labels=True)
# Plot thinned dendrogram elements

for ic, dc, cc in zip(dend['icoord'], dend['dcoord'], dend['color_list']):
    # Adjust first position to 0
    rescaled_ic = np.array(ic) - np.min(np.min(dend['icoord']))
    # scale S.T. maximum value becomes 1
    rescaled_ic /= (np.max(np.max(dend['icoord'])) - np.min(np.min(dend['icoord'])))
    # scale S.T. maximum value becomes number of events - 1
    rescaled_ic *= len(dend['leaves']) - 1
    # PLOT PRETTY & SCALED LINKS
    if cc == 'k':
        ax.plot(rescaled_ic,dc, lw=0.5, color=cc)
    else:
        ax.plot(rescaled_ic, dc, lw=0.75, color=cc)

# Format Y-axis for Coherence labeling
ax.set_ylabel('Coherence Linkage Level')
ax.set_yticks(np.arange(0,0.9, 0.1), [f'{_e:.1f}' for _e in np.arange(1.0,0.1, -0.1)])
ax.set_ylim([0,ax.get_ylim()[1]])

# Format X-axis
ticks = []
labels = []
for _gn, _gc in df_cat.group.value_counts().items():
    if _gc > 1:
        _df = df_cat[df_cat.group==_gn]
        _lx = _df.leafpos.mean()
        _ll = int(_df.tidy_group.mean())
        ticks.append(_lx)
        labels.append(f'{_ll:d}')
ax.set_xticks(ticks, labels=labels, fontsize=4)

ax.set_xlim([-1, len(df_cat)])



# Label groups with counts greater than X
# for _gn , _gc in df_cat.group.value_counts().items():
#     if _gc > 10:
#         _df = df_cat[df_cat.group==_gn]
#         _lx = _df.leafpos.mean()
#         _ll = f'Grp {int(_df.tidy_group.values[0]):d}'
#         _lc = _df.leaf_color.values[0]

#         ax.text(_lx, 0.8, _ll, color=_lc)



