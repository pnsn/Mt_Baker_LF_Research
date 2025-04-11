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
COHD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'coh_shift_table.csv'
# Template directory
TMPD = ROOT / 'processed_data' / 'template' / 'single_station'

# Coherence Clustering Params

CT = 0.3
LINKAGE = 'average'
usepref = True

# Define preferred net-sta-loc-chans
PREF_NSLC = ['UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
             'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
             'UW.JCW..EHZ',
             'UW.CMW..EHZ',
             'UW.SHUK..BHZ','UW.SHUK..HHZ',
             'CN.VDB..EHZ','CN.VDEB..HHZ',
             'UW.MULN..HHZ'] 

## SOLUTION/OUTPUT CONTROL
issave = True
isshow = True


### SUPPORTING METHODS ###

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



####### PROCESSING SECTION ########

# Scrape metadata from templates (without loading them) to inform waveform loading
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
    line = stachan.split('.') + [year, revstat, evid, rmap[revstat], _f, stachan in PREF_NSLC]
    holder.append(line)

# Form template profile dataframe
df_tpk = pd.DataFrame(holder, columns=['net','sta','loc','chan','year','revstat','evid','rstat_ohe','file','pref_nslc'])

## COHERENCE MATRIX PROSTPROCESSING ##

# Read Coherence Sparse Matrix
df_coh = pd.read_csv(COHD)

# Reconstitute individual coherence and shift matrices for preferred templates
coh_dict = {}
shift_dict = {}
for _k in df_coh.trace.unique():
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)

# Merge for preferred templates
coh_merge = pd.DataFrame()
for _k, _v in coh_dict.items():
    print(f'merging {_k}')
    if _k in PREF_NSLC and usepref:
        coh_merge = join_cov_df(coh_merge, _v)
    elif usepref:
        continue
    else:
        coh_merge = join_cov_df(coh_merge, _v)


# Fill NaN entries with 0
coh_merge = coh_merge.replace(to_replace={np.nan: 0})

# Convert to distance matrix
dmerge = 1. - coh_merge

# Run agglomerative clustering
model = AgglomerativeClustering(
    linkage=LINKAGE,
    distance_threshold=1-CT,
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
dend_labels = coh_merge.index.values


## EVENT CATALOG METADATA LOAD ##
# Read catalog profile
df_cat = pd.read_csv(CATP, parse_dates=['prefor_time'])

# Update index for merging group labels with uw####### evid format
df_cat.index = [f'{x.split("/")[-2].lower()}{x.split("/")[-1]}' for x in df_cat.event_id]

# Strip high granularity QC columns
df_cat = df_cat[['etype','lat','lon','depth','prefor_time','zerr','herr','mag','magtype','offset_km','CAT0','CAT1','CAT2','WC']]

# Add CAT3 - events with viable templates based on the linkage work done above
df_cat = df_cat.assign(CAT3=df_cat.index.isin(coh_merge.index))
# Conduct concatenation with clustering groups
df_cat = df_cat.join(
    pd.Series(model.labels_, index=dmerge.index, name='group'),
    how='left')

### GENERATE DENDROGAM ELEMENTS WITH COLOR FORMATTING
ngrps = sum(df_cat.group.value_counts() > 1)
dcmap = plt.cm.get_cmap('jet_r', ngrps+4)
colors = [mcolors.to_hex(dcmap(_e)) for _e in range(ngrps)]
# Create dendrogram outputs & structure
set_link_color_palette(colors)
dend = dendrogram(
    linkmat,
    labels=dend_labels,
    distance_sort=False, count_sort='ascending',
    color_threshold=1-CT, above_threshold_color='k',
    leaf_rotation=0, no_plot=True,
    get_leaves=True)


# Create tidy_group field
tidy_group = df_cat.group.copy().values
_e = 1
for _grp, _count in df_cat.group.value_counts(ascending=False).items():
    if _count > 1:
        tidy_group[df_cat.group == _grp] = _e
        _e += 1
    elif _count == 1:
        tidy_group[df_cat.group == _grp] = -1

# Assign non-analyzed events as -2?
tidy_group[df_cat.group.isna()] = -2

df_cat = df_cat.assign(tidy_group=tidy_group)

link_level = []
leaf_coord = []
link_end = []
for ic, dc in zip(dend['icoord'], dend['dcoord']):
    # If start of U presents as a leaf
    if ic[0]%5 == 0 and dc[0] == 0:
        link_level.append(dc[1])
        leaf_coord.append(ic[0])
        link_end.append(ic[2])
    # If end of U presents as a leaf 
    if ic[-1]%5 == 0 and dc[-1] == 0:
        link_level.append(dc[2])
        leaf_coord.append(ic[-1])
        link_end.append(ic[-3])

# Get leaf coordinates and colors associated back to df_cat
df_leaf = pd.DataFrame(
    data={'leafpos': range(len(dend['leaves'])),'leaf_color':dend['leaves_color_list'],
          'link_level': link_level, 'leaf_coord': leaf_coord, 'link_end': link_end}, 
                    index=dend['ivl'])
df_cat = df_cat.join(df_leaf, how='left')



### A PRIORI RECLASSIFICATION STRUCTURE
# This must be consistent with input "CPD" file and is based
# on prior analysis of grouping results. No logical checks are done
# in this script - strictly for plotting!
# All group numbers correspond to the "tidy_group" field.
# "Ungrouped" Group
UG_set = set([-1])
# Native all-LF groups
is_LF_set = set([6, 8, 38, 51])
# Native all-probable blast group(s)
is_PX_set = set([57])
# Group(s) judged to warrant conversion to all-SU
to_SU_set = set([1, 80])
# Group(s) judged to warrant conversion to all-EQ
to_EQ_set = (set([26]))
# Groups judged to warrant conversion to all-LF
to_LF_set = set([20, 63])
# Groups judged to warrant conversion to all-PX
to_PX_set = set([52])

# Using interpreted group reassignments, create proposed etype colum (petype)
# Hard-set group set interpretations (from visual analysis!)
superset = set(df_cat.tidy_group.unique())
# Form sets by unions
PX_set = is_PX_set.union(to_PX_set)
SU_set = to_SU_set
LF_set = is_LF_set.union(to_LF_set)
# Everything else is EQ
EQ_set = superset.difference(PX_set, SU_set, LF_set, UG_set)
petype = df_cat.etype.copy()
petype.name='petype'
for _ret, _set in zip(['su','eq','lf','px'], [to_SU_set, to_EQ_set, to_LF_set,to_PX_set]):
    petype[df_cat.tidy_group.isin(_set)] = _ret

df_cat = pd.concat([df_cat, petype], axis=1, ignore_index=False)



if isshow:
    fig = plt.figure()
    ax = fig.add_subplot()
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
    plt.title(f'CT: {CT} | METHOD: {LINKAGE} | PREFERRED?: {str(usepref)}')

    ax.set_ylabel('Coherence Linkage Level')
    ax.set_yticks(np.arange(0,0.9, 0.1), [f'{_e:.1f}' for _e in np.arange(1.0,0.1, -0.1)])
    ax.set_ylim([0,ax.get_ylim()[1]])

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

    ax.set_xlim([-1, len(df_cat[df_cat.CAT3])])

if issave:
    df_tpk.to_csv(ROOT/'results'/'tables'/'SSA2025'/'template_profile.csv', header=True,index=False)
    np.save(ROOT/'results'/'tables'/'SSA2025'/f'linkmat_{LINKAGE}_cct{CT:.3f}.npy', linkmat)
    df_cat.to_csv(ROOT/'results'/'tables'/'SSA2025'/f'catalog_profile.csv')
    for _k, _v in dend.items():
        np.save(ROOT/'results'/'tables'/'SSA2025'/f'dend_{_k}_.npy', _v)


if isshow:
    plt.show()
