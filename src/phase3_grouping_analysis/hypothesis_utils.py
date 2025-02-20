import pandas as pd
import numpy as np
from sklearn.metrics import (
    rand_score, adjusted_rand_score,
    mutual_info_score, normalized_mutual_info_score, adjusted_mutual_info_score,
    silhouette_samples, v_measure_score)

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
    cov_joined = cov_joined.combine_first(cov_joined.t)
    return cov_joined

def exact_match_score(labeling1, labeling2):
    if len(labeling1) != len(labeling2):
        raise AttributeError
    arr = []
    for _ii in range(len(labeling1)):
        arr.append(labeling1[_ii] == labeling2[_ii])
    val = sum(arr)/len(arr)
    return val

def flex_label_score(labeling1, labeling2, delimiter=' '):
    if len(labeling1) != len(labeling2):
        raise AttributeError
    matches = 0
    for _ii in range(len(labeling1)):
        parts1 = str(labeling1[_ii]).split(delimiter)
        parts2 = str(labeling2[_ii]).split(delimiter)
        ninter = len(set(parts1).intersection(set(parts2)))
        ntot = len(set(parts1).union(set(parts2)))
        matches += ninter/ntot
    val = matches/len(labeling1)
    return val


def assess_labeling(df, labeling1, labeling2):
    """Apply a series of labeling assessments on two labelings of the same data
    Methods

    exact - sum(exact_matches)/sum(total_features)
        exact matches are those that satisfy featureX[labeling1] == featureX[labeling2]
        see :meth:`~exact_match_score`

    flexible - count(partial_matches)/sum(total_features)
        partial_matches are defined as any match between elements of space delimited label strings
        e.g., 'lf su' could partial match with 'lf' or 'su'
        see :meth:`~.flex_label_score`

    rand - Rand Index :meth:`~sklear.metrics.rand_score`
        RI = # agreeing pairs / # pairs
        where # agreeing pairs = # True Positives + # True Negatives

    arand - Adjusted Rand Index :meht:`~sklearn.metrics.adjusted_rand_score`   

    

    :param df: _description_
    :type df: _type_
    :param labeling1: _description_
    :type labeling1: _type_
    :param labeling2: _description_
    :type labeling2: _type_
    :return: _description_
    :rtype: _type_
    """    
    _df = df[(df[labeling1].notna()) & (df[labeling2].notna())]
    _ii = _df[labeling1].values
    _jj = _df[labeling2].values
    out = {
        'exact': exact_match_score(_ii, _jj),
        'flexible': flex_label_score(_ii, _jj),
        'VM': v_measure_score(_ii, _jj),
        'RI': rand_score(_ii, _jj),
        'ARI': adjusted_rand_score(_ii, _jj),
        'MI': mutual_info_score(_ii, _jj),
        'AMI': adjusted_mutual_info_score(_ii, _jj),
        'NMI': normalized_mutual_info_score(_ii, _jj),
        'count': len(_df),
        'unique_pairs': len(_df[[labeling1,labeling2]].value_counts()),
        'left': labeling1,
        'right': labeling2
    }
    return out