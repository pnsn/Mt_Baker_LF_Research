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
coh_merge_pref = pd.DataFrame()
for _k, _v in coh_dict.items():
    print(f'merging {_k}')
    if _k in PREF_NSLC:
        coh_merge_pref = join_cov_df(coh_merge_pref, _v)
    coh_merge = join_cov_df(coh_merge, _v)


# Fill NaN entries with 0
coh_merge = coh_merge.replace(to_replace={np.nan: 0})
coh_merge_pref = coh_merge_pref.replace(to_replace={np.nan: 0})
# Convert to distance matrix
dmerge = 1. - coh_merge
dmerge_pref = 1. - coh_merge_pref

# Run agglomerative clustering
model = AgglomerativeClustering(
    linkage=LINKAGE,
    distance_threshold=1-CT,
    n_clusters=None,
    metric='precomputed')


model_all = model.fit(dmerge.values)
model_pref = model.fit(dmerge_pref.values)

if usepref:
    model = model_pref
else:
    model = model_all

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



if issave:
    df_tpk.to_csv(ROOT/'results'/'tables'/'SSA2025'/'template_profile.csv', header=True,index=False)
    np.save(ROOT/'results'/'tables'/'SSA2025'/f'linkmat_{LINKAGE}_cct{CT:.3f}.npy', linkmat)
    df_cat.to_csv(ROOT/'results'/'tables'/'SSA2025'/f'catalog_profile.csv')



## ALL OF THIS TO JIT
# Conduct concatenation with clustering groups
df_cat = df_cat.join(pd.Series(model.labels_, index=dmerge.index, name='group'), how='left')



# Create "tidy group index" for display purposes
tidy_group = np.ones(len(df_cat))*-1
for _e, (_gn, _gc) in enumerate(df_cat.group.value_counts().items()):
    # # Remove not included
    # if not np.isfinite(_gn):
    #     tidy_group[df_cat.group==_gn] = np.nan
    if _gc > 1:
        tidy_group[df_cat.group==_gn] = _e + 1

df_cat = df_cat.assign(tidy_group=tidy_group)

if isshow:
plt.figure()
_df = df_cat[df_cat.group.notna()]
dend = dendrogram(linkmat, color_threshold = 1 - CT, orientation='right', 
                    labels=[f'{row.mag:.1f}|{row.etype}|{row.depth*1e-3:.1f}|{row.tidy_group:.0f}|{idx}' for idx, row in _df.iterrows()],
                    count_sort='descending',
                    above_threshold_color='k')
plt.title(f'CT: {CT} | METHOD: {LINKAGE} | PREFERRED?: {str(usepref)}')

if isshow:
    plt.show()
