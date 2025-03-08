"""
:module: src/phase3_grouping_analysis/step2_assess_ensemble_grouping.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose:
    This script seeks to find optimal correlation coefficient thresholds
    for station-ensembled waveform cross-correlation based grouping of
    events around Mount Baker based on the following metrics:
    manual - user-specified cutoff for "what looks good"
    nmi - normalized mutual info score using correlation-based groups
        as "true" labels and event "etype" labels as the "predicted"
        labels
        see :meth:`~sklearn.metrics.normalized_mutual_info_score`
    ami - adjusted mutual info score using same assignments for inputs as `nmi`
            see :meth:`~sklearn.metrics.adjusted_mutual_info_score`
    ari - adjusted rand index using same assignments for inputs as `nmi`
            see :meth:`~sklearn.metrics.adjusted_rand_score`
    igh - `intra-group-homogeneity`. An ad-hoc metric.
        Assesses how many "predicted" labels differ from the dominant 
        "predicted" label within each non-singleton (1 member) group
        in "true" labels. Similar to  the methods above, `igh` approaches 1 
        with perfect homogeneity within each non-singleton group and approaches
        0 with highly heterogeneous "predicted" labelings within each "true"
        labeling group. In its current form this metric never reaches 0.
        Unlike the metrics above, `igh` is non-symmetric.
"""

import logging
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

from eqcutil.util.logging import setup_terminal_logger

from hypothesis_utils import *

# Setup logger
Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# # Reviewer / AQMS classes
# REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Event Bank Base path
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Event Space/Time cluster table
XTCD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'dist_table.csv'
# Coherence distance table
COHD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'coh_shift_table.csv'
# EVID ETYPE Reference table
EVETD = ROOT / 'data' / 'Events' / 'MtBaker_EVID_ETYPE.csv'
# Save Directory
SAVEDIR = ROOT / 'results' / 'tables'


# Read precomputed AQMS hypocenter distance table
df_dxt = pd.read_csv(XTCD)
# Read precomputed correlation coherence event-event-station table
df_coh = pd.read_csv(COHD)

# Read evid to etype mapping table
df_ee = pd.read_csv(EVETD, index_col=['uw_evid'])

# Iterate across NSLC and reconstitute an array for each
coh_dict = {}
shift_dict = {}

# Reconstitute in descending event count order
for _k in df_coh.trace.value_counts().index:
    Logger.info(f'reconstituting {_k}')
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)

# Compose total merge of coh matrices
coh_merge = pd.DataFrame()
for _v in coh_dict.values():
    coh_merge = join_cov_df(coh_merge, _v)
# Fill unpopulated values
coh_merge.fillna(0, inplace=True)

# Convert to distances
X = 1. - coh_merge.values


adkw = {'n_clusters': None,
                'linkage': 'single',
                'metric': 'precomputed'}

## GROUPING 0: "By Eye" Manual Correlation Threshold Grouping
# Run manaual preferred clustering threshold grouping
cct_pref_manual = 0.45
model_manual = AgglomerativeClustering(
    **adkw,
    distance_threshold=1. - cct_pref_manual).fit(X)

# Initialize group labeling dataframe starting with manually assigned 
df_pref = pd.DataFrame(data=model_manual.labels_, index=coh_merge.index, columns=['manual'])

## GROUPING 1: Multi-Metric Optimized Threshold Grouping
SEARCH_RANGE = 0.1
SEARCH_STEP = 0.001
# Run "group label homogeneity" optimizing sweep"
best_scores = dict(zip(['nmi','ari','ami','igh','mean'],[0]*5))
best_scores_cct = best_scores.copy()
etypes = df_ee.loc[coh_merge.index].etype
cct_search = np.arange(
    cct_pref_manual - SEARCH_RANGE,
    cct_pref_manual + SEARCH_RANGE + SEARCH_STEP,
    SEARCH_STEP
)

for _cct in cct_search:
    # Fit to new threshold
    _fit = AgglomerativeClustering(**adkw, distance_threshold=1. - _cct).fit(X)
    # Check by a selection of metrics
    nmi = normalized_mutual_info_score(_fit.labels_, etypes)
    ari = adjusted_rand_score(_fit.labels_, etypes)
    ami = adjusted_mutual_info_score(_fit.labels_, etypes)
    igh = intragroup_homogeneity(_fit.labels_, etypes)
    mean = np.mean([nmi, ari, ami, igh])
    # Update best-fit scores individually by metric
    for _k, _v in zip(best_scores.keys(),[nmi,ari,ami,igh,mean]):
        if _v > best_scores[_k]:
            best_scores[_k] = _v
            best_scores_cct[_k] = _cct



# Take mean of best score correlation thresholds
cct_post_mean = 0
for _k in ['nmi','ari','ami','igh']:
    cct_post_mean += best_scores_cct[_k]
cct_post_mean /= 4.
mod_pref = AgglomerativeClustering(
    **adkw, distance_threshold=1. - cct_post_mean).fit(X)

_ser = pd.Series(mod_pref.labels_, index=coh_merge.index, name='post_mean')
df_pref = pd.concat([df_pref, _ser], axis=1, ignore_index=False)

# Iterate across individual metric's "optimal" grouping and save grouping
for _k, _v in best_scores_cct.items():
    _fit = AgglomerativeClustering(
        **adkw, distance_threshold=1. - _v).fit(X)
    _ser = pd.Series(index=coh_merge.index, data=_fit.labels_, name=_k)
    df_pref = pd.concat([df_pref, _ser], axis=1, ignore_index=False)
    



## OUTPUT SECTION
#  attach mean and manual values to best_scores_cct
best_scores_cct.update({'post_mean': cct_post_mean, 'manual': cct_pref_manual})
# Attach cct values to labelings
df_out = pd.concat([df_pref.T, pd.Series(best_scores_cct, name='cct')], axis=1).T.sort_index()
df_out.to_csv(str(SAVEDIR/'ensemble_cct_groupings.csv'), header=True, index=True)




    
    
# # Make a series of ensemble labels
# ser_e = pd.Series(model_manual.labels_, index=coh_merge.index, name='ensemble')
# cct_ind = np.arange(0.3, 0.8, 0.025)

# # Iterate across individual channels
# holders = defaultdict(list)
# for _k, _v in coh_dict.items():
#     _X = 1. - _v.values
#     _ids = _v.index
#     # Iterate across a band of correlation thresholds
#     scores = []
#     for _cct in cct_ind:
#         _mod = AgglomerativeClustering(
#             n_clusters=None,
#             linkage='single',
#             metric='precomputed',
#             distance_threshold=1. - _cct
#         ).fit(_X)
#         _ser = pd.Series(_mod.labels_, index=_ids, name=f'cct{_cct:.2f}')
#         nmi = normalized_mutual_info_score(_ser.values, ser_e.loc[_ser.index])
#         scores.append(nmi)
#     holder.append(scores)

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ch = ax.pcolor(holder, vmin=0.5, vmax=1.)
# ax.set_yticks(np.array(range(len(coh_dict))) + 0.5, [f'{_k}({len(_v)})' for _k, _v in coh_dict.items()])
# ax.set_xticks(np.array(range(len(cct_ind)))[::2] + 0.5, [f'{_c:.3f}' for _c in cct_ind[::2]])
# fig.colorbar(ch)
    


# # # Iterate across thresholds and test against consistency of etype labeling
# # for cthr in np.arange(.3, .8, 0.01):
# #     cmod = AgglomerativeClustering(
# #         distance_threshold=cthr,
# #         n_clusters=None,
# #         metric='precomputed',
# #         linkage='single').fit(X)
# #     _ser_ct = pd.Series(cmod.labels_, index=coh_merge.index, name=f'cc{cthr:.2f}')
# Z = get_linkage_matrix(mod_pref)
# fig = plt.figure(figsize=(6,6))
# ax = fig.add_subplot(111)
# outs = dendrogram(Z, color_threshold=1. - cct_pref, labels=df_ee.loc[coh_merge.index.values].etype.values)


# plt.show()