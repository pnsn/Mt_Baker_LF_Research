"""
SAME FEATURE COMPARISON TESTS
-----------------------------

EXACT MATCH - ad hoc method assessing # exact matches of labels i & j for same feature k
   divided by the number of features k
   EM = sum(ik == jk) / k

FLEXIBLE MATCH - ad hoc method assessing # matches of white-space delimited, compound labels
    e.g., 'lf su' -> ['lf','su'] vs 'lf' -> ['lf'] with
    FM = sum(intersection(parts1_k, parts2_k)/union(parts1_k, parts2_k))/ count(k)

MUTUAL INFORMATION, NORMALIZED MI, ADJUSTED MI
 - MI measures the agreement between the labels for the same feature and places
  difference in precise label in the context of their label grouping with other features

PAIRWISE FEATURE COMPARISON TESTS - Have I grouped my data in a similar manner?
---------------------------------
RAND_INDEX and ADJUSTED_RAND_INDEX (note Rand is the statistician's name)
 - symmetric -> can be thought of as concesus measurements
 - proportional to the number of samples pairs that are both the same
    or both different for sample pairs i & j with i != j
 - require knowledge of the ground truth classes, but can still be useful for
   testing concensus for "unsupervised" groupings

ADJUSTED_RAND_INDEX is adapted to approach 0.0 for any number of clusters and samples
    that appear to be drawn from a random labeling. As such, the ARI is in [-0.5, 1]

My Best Pass at plain-language questions we ask with each of these tests

RAND_INDEX - Did this labeling group my features in a similar structure as my baseline?
    
ADJUSTED_RAND_INDEX - Did this labeling group my features in a similar structure 
    to my baseline and in a way that is significantly different from a random labeling?

EXACT_MATCH - How identically did this labeling group my features compared to my baseline?

FLEXIBLE_MATCH - Given multiple 

"""

import logging
from pathlib import Path

import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from eqcorrscan.utils.clustering import handle_distmat_nans

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_mutual_info_score, rand_score, adjusted_rand_score
from sklearn.cluster import AgglomerativeClustering
from eqcutil.util.logging import setup_terminal_logger
from hypothesis_utils import *

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# Reviewer / AQMS classes
REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Event Space/Time cluster table
XTCD = ROOT / 'results' / 'tables' / 'event_distance_table.csv'
# Coherence distance table
COHD = ROOT / 'results' / 'tables' / 'coherence_distance_table.csv'

Logger = setup_terminal_logger(name=__name__, level=logging.INFO)


Logger.info('reading analyst reviewed event table')
# Read reviewed event results
df_rev = pd.read_csv(REVD, index_col=[0])
# Encode Etypes
mapping = {}
_e = int(0)
for _, row in df_rev.iterrows():
    for _, etype in row.items():
        if etype not in mapping.keys():
            # Ignore nan for up-indexing
            if isinstance(etype, str):
                mapping[etype] = _e
                _e += int(1)
            # Include nan to nan mapping
            else:
                mapping[etype] = etype
df_rev_encoded = df_rev.copy().replace(mapping)


Logger.info('reading hypocenter distance table')
# Read precomputed AQMS hypocenter distance table
df_dxt = pd.read_csv(XTCD)
Logger.info('reading coherence matrix table')
# Read precomputed correlation coherence event-event-station table
df_coh = pd.read_csv(COHD)


# Subset df_coh and df_dxt by accepted values
wcset = set(df_rev.index)
df_coh = df_coh[(df_coh.event_i.isin(wcset)) & (df_coh.event_j.isin(wcset))]
df_dxt = df_dxt[(df_dxt.event_i.isin(wcset)) & (df_dxt.event_j.isin(wcset))]

# Process individual coherence and shift matrices
coh_dict = {}
shift_dict = {}
for _k in df_coh.trace.unique():
    Logger.info(f'reconstituting {_k}')
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)

# Process distance tables into distance matrices for each metric
dist_dict = {}
for _k in df_dxt.columns:
    if _k not in ['event_i','event_j']:
        dist_dict[_k] = get_symmetric(df_dxt, k_field=_k, trace_value=0.)


# HYPOTHESIS 0: Analysts agree with AQMS
df_H0 = pd.DataFrame()
for _r in range(4):
    revname = f'R{_r}'
    testscores = assess_labeling(df_rev, 'AQMS', revname)
    df_H0 = pd.concat([df_H0, pd.DataFrame(testscores, index=[_r])],axis=0,ignore_index=False)

print(df_H0[['exact','flexible','RI','ARI','count','unique_pairs','left','right']])

# HYPOTHESIS 1: Analysts concensus
df_H1 = pd.DataFrame()
for _r in range(3):
    for _s in range(_r+1, 4):
        rname = f'R{_r}'
        sname = f'R{_s}'
        testscores = assess_labeling(df_rev, rname, sname)
        df_H1 = pd.concat([df_H1, pd.DataFrame(testscores, index=[0])], ignore_index=True, axis=0)

# print(df_H1[['exact','flexible','']])


# CLUSTERING TESTING
"""
I want to ascertain if there is a distinct distance threshold and correlation threshold pair
that optimizes location based groupings S.T. we can isolate
"""
# Iterate across agglomeration method
for lnkg in ['single','average','complete']:
    # Iterate across coarse range of clustering distances
    for hthr in np.arange(500, 10000, 500):
        # Iterate across coarse range of correlation thresholds
        for cthr in np.arange(0.1,0.8, 0.05):
            # Iterate across stations
            for _k, _v in coh_dict.items():
                # Subset event_distances
                _dfh = dist_dict['delh_ij_m']
                _dfh = _dfh.loc[_v.index, _v.index]
                _dfz = dist_dict['delz_ij_m']
                _dfz = _dfz.loc[_v.index, _v.index]
                _dfd = (_dfh**2 + _dfz**2)**0.5
                dmodel = AgglomerativeClustering(distance_threshold=hthr, n_clusters=None, metric='precomputed', linkage=lnkg)
                dmodel = dmodel.fit(_dfd.values)
                _dgrp = pd.Series(dmodel.labels_, index=_dfd.index)
                cmodel = AgglomerativeClustering(distance_threshold=cthr, n_clusters=None, metric='precomputed', linkage=lnkg)
                cmodel = cmodel.fit(_v.values)
                _cgrp = pd.Series(cmodel.labels_, index=_v.index)
                _df_grp = pd.concat([_dgrp, _cgrp], axis=1, ignore_index=False)
                breakpoint()
                # TODO: Run testsuite and update best score iteratively