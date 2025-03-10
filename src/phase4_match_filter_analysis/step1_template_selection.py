import logging, glob
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.cluster import AgglomerativeClustering

from eqcutil import ClusteringTribe

ROOT = Path(__file__).parent.parent.parent
CTRDIR = ROOT / 'processed_data' / 'cluster' / 'single_station'

PREFNSLC = ['UW.JCW..EHZ',
            'UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.CMW..EHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'UW.PASS..BHZ','UW.PASS..HHZ',
            'CN.VDB..EHZ','CN.VDEB..HHZ',
            ]

MIN_SNR = 3.

### SUPPORTING METHODS ###
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


def get_linkage_matrix(model):
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
    Z = np.column_stack([model.children_, model.distances_, counts]).astype(float)
    return Z

### PROCESSING SECTION ###
ctr_dict = {}
coh_merged = pd.DataFrame()
for nslc in PREFNSLC:
    fname = CTRDIR / f'{nslc}_clustered.tgz'
    print(f'reading {fname}')
    # Load clustering tribe
    _ctr = ClusteringTribe().read(str(fname))
    ctr_dict[nslc] = _ctr
    # Subset by SNR
    _ictr = _ctr.get_subset(
        _ctr._c[_ctr._c.mean_snr_dB >= MIN_SNR].index
    )
    # Create labeled dist_mat
    df_dist = pd.DataFrame(data=_ictr.dist_mat, index=_ictr._c.index, columns=_ictr._c.index)
    # Convert to coherence score and merge
    coh_merged = join_cov_df(coh_merged, 1 - df_dist)

model = AgglomerativeClustering(
    n_clusters=None,
    metric='precomputed',
    linkage='single',
    distance_threshold=1 - 0.445
)
X = coh_merged.values
