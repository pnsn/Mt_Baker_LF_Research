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
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from pyproj import Proj

from eqcorrscan.utils.clustering import handle_distmat_nans
from obsplus import EventBank

from sklearn import metrics
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import AgglomerativeClustering, DBSCAN, OPTICS

from eqcutil.util.logging import setup_terminal_logger

from hypothesis_utils import *

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# Reviewer / AQMS classes
REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Event Bank Base path
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Event Space/Time cluster table
XTCD = ROOT / 'results' / 'tables' / 'event_distance_table.csv'
# Coherence distance table
COHD = ROOT / 'results' / 'tables' / 'coherence_distance_table.csv'

Logger = setup_terminal_logger(name=__name__, level=logging.INFO)
# Connect to eventbank
EBANK = EventBank(EBBP)

# Define Lat Lon to UTM converter
myproj = Proj(proj='utm', zone=10, ellips='WGS84', preserve_units=False)

Logger.info('reading analyst reviewed event table')
# Read reviewed event results
df_rev = pd.read_csv(REVD, index_col=[0])
# # Encode Etypes
# mapping = {}
# _e = int(0)
# for _, row in df_rev.iterrows():
#     for _, etype in row.items():
#         if etype not in mapping.keys():
#             # Ignore nan for up-indexing
#             if isinstance(etype, str):
#                 mapping[etype] = _e
#                 _e += int(1)
#             # Include nan to nan mapping
#             else:
#                 mapping[etype] = etype
# df_rev_encoded = df_rev.copy().replace(mapping)


Logger.info('reading hypocenter distance table')
# Read precomputed AQMS hypocenter distance table
df_dxt = pd.read_csv(XTCD)
Logger.info('reading coherence matrix table')
# Read precomputed correlation coherence event-event-station table
df_coh = pd.read_csv(COHD)

Logger.info('reading eventbank index')
df_eb = EBANK.read_index()

# Convert into UTM Zone 10N
x, y = myproj(df_eb.longitude.values, df_eb.latitude.values)
df_eb = df_eb.assign(mE=x)
df_eb = df_eb.assign(mN=y)
# Convert into XYZ array with evid annotations
df_XYZ = pd.DataFrame(
    df_eb[['mE','mN','depth','horizontal_uncertainty','vertical_uncertainty']].values,
    columns=['mE','mN','mZ','mSH','mSZ'],
    index=[f'uw{_id.split("/")[-1]}' for _id in df_eb.event_id]
    )

# Subset df_coh and df_dxt by accepted values
wcset = set(df_rev.index)
# df_coh = df_coh[(df_coh.event_i.isin(wcset)) & (df_coh.event_j.isin(wcset))]
# df_dxt = df_dxt[(df_dxt.event_i.isin(wcset)) & (df_dxt.event_j.isin(wcset))]

# Process individual coherence and shift matrices
coh_dict = {}
shift_dict = {}
for _k in df_coh.trace.unique():
    Logger.info(f'reconstituting {_k}')
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)

coh_merge = pd.DataFrame()
for _v in coh_dict.values():
    coh_merge = join_cov_df(coh_merge, _v)

# Fill remaining NaN values
coh_merge = pd.DataFrame(
    data = handle_distmat_nans(coh_merge.values, replace_nan_distances_with='mean'),
    index=coh_merge.index,
    columns=coh_merge.columns)                     


# breakpoint()

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


# CLUSTERING EXPERIMENTS
"""
I want to ascertain if there is a distinct distance threshold and correlation threshold pair
that optimizes location based groupings and correlation based groupings. In this first experiment
I'm testing out different agglomerative clustering algorithms and recording the correlation and
distance threshold values that maximize individual concensus metrics.

From initial examination, it appears that this approach favors reducing down to an unreasonably
small number of clusters that maximize the MI and RI metrics, as well as the ad hoc methods
I've included for one to one labeling comparisons in the analyst concensus analysis above
"""
run_clustering_experiment = 1
if run_clustering_experiment == 0:
    bestscores = {}
    # Iterate across agglomeration method
    for lnkg in ['single','average','complete']:
        Logger.info(lnkg)
        # Iterate across coarse range of clustering distances [m]
        for hthr in np.arange(500, 10000, 500):
            # Logger.info(hthr)
            # Iterate across coarse range of correlation thresholds
            for cthr in np.arange(0.2, 0.9, 0.025):            
                _dfh = dist_dict['delh_ij_m']
                _dfh = _dfh.loc[coh_merge.index, coh_merge.index]
                _dfz = dist_dict['delz_ij_m']
                _dfz = _dfz.loc[coh_merge.index, coh_merge.index]
                _dfd = (_dfh**2 + _dfz**2)**0.5
                # Cluster on calculated 3D distances
                dmodel = AgglomerativeClustering(distance_threshold=hthr, n_clusters=None, metric='precomputed', linkage=lnkg)
                dmodel = dmodel.fit(_dfd.values)
                _dgrp = pd.Series(dmodel.labels_, index=_dfd.index, name='distance')
                cmodel = AgglomerativeClustering(distance_threshold=cthr, n_clusters=None, metric='precomputed', linkage=lnkg)
                # Cluster on coh_dist = 1 - coh to have distance matrix
                cmodel = cmodel.fit(1. - coh_merge.values)
                _cgrp = pd.Series(cmodel.labels_, index=coh_merge.index, name='correlation')
                _df_grp = pd.concat([_dgrp, _cgrp], axis=1, ignore_index=False)
                testscores = assess_labeling(_df_grp, 'distance','correlation')
                if lnkg not in bestscores.keys():
                    bestscores[lnkg] = {}
                idx = bestscores[lnkg]
                for _l, _w in testscores.items():
                    if _l in ['exact','flexible','VM','RI','ARI','MI','AMI','NMI']:
                        pass
                    else:
                        continue
                    if _l not in idx.keys():
                        idx[_l] = [_w, lnkg, hthr, cthr]
                    elif idx[_l][0] < _w:
                        Logger.info(f'New improvement: {_k}, {_l}, {_w}, {lnkg}, {hthr}, {cthr}, {len(_dgrp.unique())}, {len(_cgrp.unique())}')
                        idx[_l] = [_w, lnkg, hthr, cthr]
                    else:
                        continue

if run_clustering_experiment == 1:
    """
    In this experiment I am approaching comparative grouping from another angle. We have
    event location estimates (and uncertainties for most locations) and there are sets of
    density based clustering methods that can reject outliers. This is a similar structure
    as observed in single linkage clustering of cross correlation coherences, so here
    we test a set of denstity based clusterings on location estimates and compare them
    to agglomeratively clustered correlation coherence "distances" over a range of threshold values
    """
    uxyz = df_XYZ[df_XYZ.index.isin(coh_merge.index)][['mE','mN','mZ']].values
    sxyz = df_XYZ[df_XYZ.index.isin(coh_merge.index)][['mSH','mSH','mSZ']].values
    # Create standardized geometry representation
    uxyz_ss = StandardScaler().fit_transform(uxyz)

df_dbl = pd.DataFrame()
df_dbs = pd.DataFrame()
Logger.info('Ru')
for eps in np.arange(0.05, 1., 0.05):
    for ms in range(2,6):
        # Calculate clustering
        db = DBSCAN(eps=eps, min_samples=ms).fit(uxyz_ss)
        # Annotate group labels with EVIDs
        _ser_db = pd.Series(db.labels_, index=coh_merge.index, name=f'db{eps:.2f}_ms{ms:d}')
        # Capture label set in summary DF
        df_dbl = pd.concat([df_dbl, _ser_db], axis=1, ignore_index=False)
        # Calculate Metrics
        n_clust = len(set(db.labels_)) - (1 if -1 in db.labels_ else 0)
        n_outly = list(db.labels_).count(-1)
        sil_coef = metrics.silhouette_score(uxyz_ss, db.labels_)
        _ser_dbs = pd.Series({'no_clust': n_clust, 'no_outliers': n_outly, 'sil_coef': sil_coef, 'eps': eps, 'ms': ms},
                                name=f'db{eps:.2f}_ms{ms:d}')
        df_dbs = pd.concat([df_dbs, _ser_dbs], axis=1, ignore_index=False)

df_acl = pd.DataFrame()
df_acs = pd.DataFrame()
for cthr in np.arange(0.4,0.8, 0.01):
    cmodel = AgglomerativeClustering(
        distance_threshold=cthr,
        n_clusters=None,
        metric='precomputed',
        linkage='single').fit(1. - coh_merge.values)
    # Form series
    _ser_ct = pd.Series(cmodel.labels_, index=coh_merge.index, name=f'cc{cthr:.2f}')
    # Calculate Metrics
    # Get indices of singletons (outliers)
    _outlier_mask = _ser_ct.value_counts() == 1
    # Calculate pre-masking metrics
    n_outly = sum(_outlier_mask)
    n_clust = sum(_ser_ct.value_counts() > 1)

    df_acl = pd.concat([df_acl, _ser_ct], axis=1, ignore_index=False)
    
    try:
        sil_coef_masked = metrics.silhouette_score(uxyz_ss, _ser_ct.values)
    except ValueError:
        sil_coef_masked = np.nan
    _ser_acs = pd.Series({'no_clust': n_clust, 'no_outliers': n_outly, 'sil_coef': sil_coef_masked, 'cthr': cthr},
                            name=f'cc{cthr:.2f}')
    df_acs = pd.concat([df_acs, _ser_acs], axis=1, ignore_index=False)

# # Run intercomparison
holder = []
spare_index_start = 10000
for dbpar, drow in df_dbl.T.iterrows():
    _drow = drow.copy()
    parts = dbpar.split('_')
    # Extract DBSCAN hyper parameters
    ms = int(parts[1][2:])
    eps = float(parts[0][2:])
    for _e, _v in _drow.items():
        if _v == -1:
            _drow[_e] = spare_index_start
            spare_index_start += 1
    
    for acpar, crow in df_acl.T.iterrows():
        # Excract corr threshold
        cct = float(acpar[2:])

        _crow = crow.copy()
        # join on EVIDS
        _dfc = pd.concat([_crow, _drow], axis=1, ignore_index=False)
        _dfc = _dfc[(_dfc[dbpar].notna())&(_dfc[acpar].notna())]

        # Calculate metrics
        nclust = len(_crow.value_counts()[_crow.value_counts() >= ms])
        exact = sum(_dfc[dbpar] == _dfc[acpar])/len(_dfc)
        ami = metrics.adjusted_mutual_info_score(_dfc[dbpar], _dfc[acpar])
        nmi = metrics.normalized_mutual_info_score(_dfc[dbpar], _dfc[acpar])
        ari = metrics.adjusted_rand_score(_dfc[dbpar], _dfc[acpar])    
        line = [ms, eps, cct, nclust, exact, ami, nmi, ari]
        holder.append(line)

        # # Mask outliers as defined by min membership
        # msk = crow.value_counts() < ms
        # mapper = {_k: -1 for _k in crow.value_counts()[msk].index}
        # nclust = 0
        # if len(mapper) > 0:
        #     _crow = crow.replace(mapper)
        #     nclust -= 1
        # else:

# FIXME: ami/nmi/ari zero out when there is a negative label. Need to convey non-grouped labels as unique integers

    # for eps, drow in df_dbl.T.iterrows():
    #     _dfc = pd.concat([crow, drow], axis=1, ignore_index=False)
    #     _dfc = _dfc[(_dfc[cct].notna()) & (_dfc[eps].notna())]
    #     AMI = metrics.adjusted_mutual_info_score(_dfc[cct], _dfc[eps])
    #     NMI = metrics.normalized_mutual_info_score(_dfc[cct], _dfc[eps])
    #     ARI = metrics.adjusted_rand_score(_dfc[cct], _dfc[eps])
    #     line = [cct, eps, AMI, NMI, ARI]
    #     holder.append(line)



# # Iterate across DBSCAN parameter rows
# for dbpar, dblbl in df_dbs.T.iterrows():
#     parts = dbpar.split('_')
#     # Get minimum membership for cluster
#     ms = int(parts[1][2:])
#     # Iterate across AGGCLUST parameter rows
#     for cct, agclbl in df_sla.T.iterrows():


