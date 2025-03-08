"""
:module: src/phase3_grouping_analysis/step3_assess_single_station_grouping.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose:
    This script seeks to find optimal correlation coefficient threshold
    for each observing station that best matches labelings from the 


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
# Preferred Ensemble Labelings & Correlation Coefficient Thresholds (cct)
PREF_E = SAVEDIR / 'ensemble_cct_groupings.csv'


# Cross Correlation Threshold search parameters
CCT_SEARCH_RANGE = 0.1
CCT_SEARCH_INCR = 0.001
REF_CCT_TYPE = 'mean'
# AgglomerativeClustering shared parameters
adkw = {'n_clusters': None,
        'linkage': 'single',
        'metric': 'precomputed'}


## PROCESSING SECTION ##

# Read precomputed correlation coherence event-event-station sparse table
df_coh = pd.read_csv(COHD)

# Read preferred labeling csv
df_pref = pd.read_csv(PREF_E, index_col=[0])
df_pref_cct = df_pref.iloc[0,:]
df_pref_lbl = df_pref.iloc[1:,:]

# Get starting CCT
CCT_CENTER = df_pref_cct[REF_CCT_TYPE]
# Populate range of CCT to search
CCT_RANGE = np.arange(CCT_CENTER - CCT_SEARCH_RANGE,
                      CCT_CENTER + CCT_SEARCH_RANGE + CCT_SEARCH_INCR,
                      CCT_SEARCH_INCR)


# Reconstitute by NSLC code in order of descending event count
df_bests = pd.DataFrame()
for _k in df_coh.trace.value_counts().index:
    Logger.info(f'reconstituting {_k} coherence matrix from sparse')
    # Get subset
    _df = df_coh[df_coh.trace==_k]
    # Rebuild coherence matrix from sparse representation
    coh = get_symmetric(_df, k_field='coh', trace_value=1.)
    # Get subset labeling from ensemble
    _df_ref = df_pref_lbl.loc[coh.index]
    best_scores = dict(zip(['nmi','ari','ami','igh','mean'],[0]*5))
    best_cct_lower = best_scores.copy()
    best_cct_upper = best_scores.copy()
    # Iterate across CCT
    Logger.info('scanning across correlation coefficient thresholds')
    X = 1. - coh.values
    for _cct in CCT_RANGE:
        _fit = AgglomerativeClustering(
            **adkw, distance_threshold=1. - _cct
        ).fit(X)
        # Iterate across each ensemble metric
        for _c in _df_ref.columns:
            _ser_ref = _df_ref[_c]
            nmi = normalized_mutual_info_score(_ser_ref.values, _fit.labels_)
            ami = adjusted_mutual_info_score(_ser_ref.values, _fit.labels_)
            ari = adjusted_rand_score(_ser_ref.values, _fit.labels_)
            igh = intragroup_homogeneity(_ser_ref.values, _fit.labels_)
            mean = np.mean([nmi,ami,ari,igh])
            for _l, _v in zip(best_scores.keys(), [nmi, ari, ami, igh, mean]):
                if _v > best_scores[_l]:
                    best_scores[_l] = _v
                    best_cct_lower[_l] = _cct
                    best_cct_upper[_l] = _cct
                elif _v == best_scores[_l]:
                    best_cct_upper[_l] = _cct
    _df_bests = pd.DataFrame(
        {'trace': _k,
         'count': len(_fit.labels_),
         'score': best_scores,
         'cct_l': best_cct_lower,
         'cct_u': best_cct_upper})
    _df_bests = _df_bests.assign(metric=_df_bests.index.values)
    df_bests = pd.concat([df_bests, _df_bests], axis=0, ignore_index=True)
    # df_bests = pd.concat([df_bests, _df_bests], axis=0, ignore_index=)

# SAVE RAW RESULT
df_bests.to_csv(str(SAVEDIR/'single_station_cct_opt_values.csv'), index=False, header=True)



