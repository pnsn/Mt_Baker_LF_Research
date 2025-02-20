"""
:purpose:
The goal of this analysis is to find an optimum inter-event clustering
spatial distance and correlation distance paired with sets of analyst / catalog event type labels to 


"""
import logging
from pathlib import Path

import numpy as np

from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from eqcorrscan.utils.clustering import handle_distmat_nans

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_mutual_info_score, rand_score, adjusted_rand_score
from eqcutil.util.logging import setup_terminal_logger
from hypothesis_utils import *

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# Reviewer / AQMS classes
REVD = ROOT / 'processed_data' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Event Space/Time cluster table
XTCD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'event_distance_table.csv'
# Coherence distance table
COHD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'coherence_distance_table.csv'

Logger = setup_terminal_logger(name=__name__, level=logging.INFO)


Logger.info('reading analyst reviewed event table')
# Read reviewed event results
df_rev = pd.read_csv(REVD, index_col=[0])
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

breakpoint()
# def to_symmetric(df, i_field, j_field, k_field, trace_value=1.):
#     # Create pivot table to sample from
#     pt = df.pivot_table(columns=i_field, index=j_field, values=k_field, aggfunc='mean')

#     iunique = df[i_field].unique()
#     iunique.sort()
#     junique = df[j_field].unique()
#     if set(iunique) != set(junique):
#         raise AttributeError
        
#     values = np.full(shape=(len(iunique), len(iunique)), fill_value=np.nan)
#     for ii, iv in enumerate(iunique):
#         for jj, jv in enumerate(iunique):
#             if ii == jj:
#                 values[ii, jj] = trace_value
#             else:
#                 values[ii, jj] = pt.loc[iv, jv]
#                 values[]

#     return pd.DataFrame(values, columns=iunique, index=iunique)
            
    



#     for ii, iname in enumerate(iunique):
#         for jj, jname in enumerate(junique):
#             val = df[()]





# # Reconstitute lower triangle
# _df = df_dxt.copy()
# _df = _df.rename(columns={'event_i': 'event_j',
#                           'event_j': 'event_i',
#                           'etype_i': 'etype_j',
#                           'etype_j': 'etype_i'})
# # Reconstitute diagnoal
# iparts = df_dxt[['event_i','etype_i']].value_counts().index.values
# holder = []
# for evid, etype in iparts:
#     holder.append([evid, evid, 0, 0, 0, etype,etype])
# # Reassemble into tabular format
# df_dxt = pd.concat([df_dxt, _df, pd.DataFrame(holder, columns=df_dxt.columns)],
#                    axis=0, ignore_index=True)

# # Read precomputed correlation coherence table
# df_coh = pd.read_csv(COHD)
# # Reconstitude lower triangle
# _df = df_coh.copy()
# _df = _df.rename(columns={'event_i': 'event_j',
#                           'event_j': 'event_i'})
# # Reconstitude diagonal
# iparts = df_coh[['trace','event_i']]
# holder = []
# for trid, evid in iparts.values:
#     holder.append([trid, evid, evid, 1, 0])
# df_coh = pd.concat([df_coh, _df, pd.DataFrame(holder, columns=df_coh.columns)],
#                     axis=0, ignore_index=True)




# EXPERIMENT SET 1: Well Constrained Events'
# Grouping is entirely explained by difference in source

# Get well constrained event IDs
wce = df_rev.index
# Subset distance/time table
df_dxt_wc = df_dxt[(df_dxt.event_i.isin(wce)) & (df_dxt.event_j.isin(wce))]
# Subset coherence table
df_coh_wc = df_coh[(df_coh.event_i.isin(wce)) & (df_coh.event_j.isin(wce))]

# EXPERIMENT 1A - All Stations

# Get the mean coherence matrix for all stations observing well constrained events
coh_table = df_coh_wc.pivot_table(
    index='event_i',columns='event_j',
    values='coh', aggfunc='mean')
# Fill nans
cdist_arr = handle_distmat_nans(
    1 - coh_table.values,
    replace_nan_distances_with='mean')
# Vectorize
cdist_vect = squareform(cdist_arr)
# Get linkage matrix
COH_Z = linkage(cdist_vect)


# Get inter-event distance matrix
dh_table = df_dxt_wc.pivot_table(
    index='event_i', columns='event_j',
    values='delh_ij_km', aggfunc='sum')
# dh_arr = handle_distmat_nans(
#     dh_table.values,
#     replace_nan_distances_with='mean')

# Vectorize
dh_vect = squareform(dh_table.values)
# Get linkage matrix
DHD_Z = linkage(dh_vect)

results = []

# Iterate across analyst labels
for src, labels in df_rev.T.iterrows():
    # Iterate across correlation thresholds
    for _cct in np.arange(0.1, 1, 0.05):
        comp = f'cct_{_cct:.2f}'
        # Run linkage clustering at set threshold
        coh_grps = fcluster(COH_Z, t = 1 - _cct, criterion='distance')
        # Label data
        coh_grp_series = pd.Series(data=coh_grps, index=coh_table.index, name=comp)
        # Get overlapping
        _pairs = pd.concat([labels, coh_grp_series], axis=1)
        _pairs = _pairs[(_pairs[src].notna()) & (_pairs[comp].notna())]
        nmi = normalized_mutual_info_score(_pairs[src].values, _pairs[comp].values)
        ami = adjusted_mutual_info_score(_pairs[src].values, _pairs[comp].values)
        result = [src, comp, nmi, ami, len(_pairs)]
        results.append(result)
df_result = pd.DataFrame(results)

# for rad in np.arange(1, 20.5, 0.5): 
#     for _cct in np.arange


# EXPERIMENT 1B - Individual Stations


# EXPERIMENT 1C - 

# H0: Event Labels Alone Are Sufficient To Explain All Grouping

# H1: Event Labels + Locations are Sufficient To Explain All Grouping



# EXPERIMENT SET 2: Whole Catalog