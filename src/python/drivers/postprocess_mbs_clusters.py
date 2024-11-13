"""
:module: Mt_Baker_LF_Research/src/python/drivers/postprocess_mbs_clusters.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This driver conducts the following post-processing steps
    on outputs from `cluster_mbs_tribe.py`

     - Load saved ClusteringTribe data from hyperparameter trials
     - Associate select AQMS origin metadata to clusteringtribe.clusters tables
     - Save associated `clusters` table to disk
     - Determine the best-fit correlation coefficient to explain selected thresholds
"""

import os, glob, logging
from pathlib import Path

import numpy as np
import pandas as pd
from obsplus import EventBank
from sklearn.metrics import normalized_mutual_info_score

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger


### PARAMETER SETTING SECTION ###
# Seed random number generator
np.random.seed(62323)

kwargs = {'n_depth_bins': 4,
          'average_method': 'arithmetic',
          'etype_wgt': 10}

cct_range = np.arange(0.005, 1., 0.005)

### MONTE CARLO SIMULATION SAMPLER
def generate_perturbed_value_bins(locs, scales, nbins=4, distribution='normal'):
    """Given the expected values (locs) and characteristic scales of a vector
    of data, generate perturbed data samples for a prescribed distribution type
    and bin them into equal sample-count bins

    :param locs: expected parameter values
    :type locs: vector-like
    :param scales: characteristic scales of parameter uncertainties (see **distribution**)
    :type scales: vector-like
    :param nbins: number of equal-sample scaled bins to produce, defaults to 4
        see :meth:`~pandas.qcut` - passed to parameter `q`
    :type nbins: int, optional
    :param distribution: name of the distribution to use, defaults to 'normal'
        Supported values:
            'normal' - locs are the means, scales are the standard deviations
            'uniform' - min is locs - scales, max is locs + scales
    :type distribution: str, optional
    :return:
     - **groups** (*pandas.Categorical*) - unlabeled bin designations
        of perturbed MC sample
    """    
    if distribution.lower() == 'normal':
        pert_vals = np.random.normal(loc=locs, scale=scales)
    elif distribution.lower() == 'uniform':
        pert_vals = np.random.uniform(low=locs - scales, high=locs + scales)
    else:
        raise ValueError(f'distribution {distribution} not supported')
    groups = pd.qcut(pert_vals, q=nbins, labels=False)
    return groups
        
def evaluate_score_mc(ref_grouping, locs, scales, nbins=4, ndraw=30, 
                      distribution='normal', average_method='arithmetic'):
    """Conduct a Monte Carlo simulation to assess the normalized mutual information
    score between a reference grouping and a vector of random variables defined by
    data values, uncertainties, and assumed underlying distribution. 

    :param ref_grouping: reference grouping
    :type ref_grouping: vector-like
    :param locs: expected values of random variable vector members
    :type locs: vector-like
    :param scales: uncertainty/range of random variable vector members
    :type scales: vector-like
    :param nbins: number of even-sample-count bins to group the random 
        variable by, defaults to 4
    :type nbins: int, optional
    :param ndraw: number of sample sets to draw for this simulation, 
        defaults to 30
    :type ndraw: int, optional
    :param distribution: name of the distribution to use, defaults to 'normal'
        Supported: 'normal','uniform'
    :type distribution: str, optional
    :param averaging_method: averaging method for the NMI score, defaults to 'arithmetic'
    :type averaging_method: str, optional
    :return: 
     - **scores** (*numpy.ndarray*) -- vector of scores from each MC simulation draw
    """    
    scores = []
    for _n in ndraw:
        sample = generate_perturbed_value_bins(locs, scales, nbins=nbins, distribution=distribution)
        score = normalized_mutual_info_score(ref_grouping, sample, average_method=average_method)
        scores.append(score)
    return np.array(scores)

def objective(ctr, cct, n_depth_bins=4, n_depth_draw=30, average_method='arithmetic', etype_wgt=2.):
    acdf = ctr.clusters
    locs = acdf.depth.values
    scales = acdf.vertical_uncertainty.values
    dscores = evaluate_



### EVALUATION FUNCTION ###
def evaluate_corr_thresh(clustertribe, threshold, n_depth_bins=8, average_method='arithmetic', etype_wgt=2.):
    # Get clusters as a dataframe
    clusters = clustertribe._cct_regroup(threshold)
    # df_c = pd.DataFrame(clustertribe._cct_regroup(threshold))

    # Left-join etype and depth 
    df_c = pd.concat([clusters, clustertribe.clusters[['etype','depth']]], axis=1, ignore_index=False)
    # Get cluster values
    cvals = df_c['correlation_cluster'].values
    escore = normalized_mutual_info_score(cvals, df_c.etype.values, average_method=average_method)
    dscore = normalized_mutual_info_score(cvals, pd.qcut(df_c.depth.values, n_depth_bins, labels=False), average_method=average_method)

    return (escore*etype_wgt + dscore)/(etype_wgt + 1.)

def perturb_depths(clustertribe):


def evaluate_nmi()


# Setup Logging
Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Define Root Directory
root = Path(__file__).parent.parent.parent.parent
# Safety check that it's running from root directory
if os.path.split(root)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {root}')
else:
    Logger.info(f'running {__file__} from {root}')

# Define path to origins metadata that contains event type & fixed statuses
etype_file = os.path.join(root, 'data','Events','MtBaker_50km_radius_origins.csv')
fields = ['evid','etype','fdepth','fepi','ftime']

# Define connection to EventBank in repository
ebpath = os.path.join(root,'data','XML','QUAKE','BANK')
ebkwargs = {'base_path': ebpath,
            'path_structure': '{year}',
            'name_structure': 'uw{event_id_end}'}

# Get get path to ClusteringTribe file names
processed_path = os.path.join(root,'processed_data','clusters','well_located_20km')
ctr_file_path = os.path.join(processed_path,'hyper_parameter_trials')
ctr_files = glob.glob(os.path.join(ctr_file_path, '*.tgz'))

# Save directory
save_path = os.path.join(processed_path,'assoc_groupings')

#### PROCESSING SECTION ####

# Create save directory
if not os.path.exists(save_path):
    Logger.warning(f'creating new output directory: {save_path}')
    os.makedirs(save_path)
else:
    Logger.info(f'writing to {save_path}')

# Load etype_file
df_etype = pd.read_csv(etype_file)
idx = [f'uw{row.evid}' for _, row in df_etype.iterrows()]
# Format index for joining to ctr.clusters
df_etype.index = idx
# Deduplicate based on EVID
df_etype = df_etype[fields].drop_duplicates(keep='first', subset=['evid'])


# Connect to ebank
ebank = EventBank(**ebkwargs)
# Get index
df_eb = ebank.read_index()
# Update index with COMCAT evids
idx = [eid.split('/') for eid in df_eb.event_id]
idx = [f'{eid[-2].lower()}{eid[-1]}' for eid in idx]
df_eb.index = idx

bests = {}
results = {}
# Iterate across clustering trials
for filepath in ctr_files:
    Logger.info(f'processing clustering trial: {filepath}')
    path, file = os.path.split(filepath)
    filename, ext = os.path.splitext(file)
    # LOAD
    ctr = ClusteringTribe().read(filepath)
    Logger.info('ClusteringTribe loaded')
    # Join etype
    ctr.clusters = ctr.clusters.join(df_etype, how='left')
    # Join eventbank metadata
    ctr.clusters = ctr.clusters.join(df_eb, how='left')
    Logger.info(f'Running scoring for {len(cct_range)} cross_correlation threshold values')
    scores = []
    for _e, thresh in enumerate(cct_range):
        Logger.debug(f'{_e+1} of {len(cct_range)}')
        scores.append(evaluate_corr_thresh(ctr, thresh, **kwargs))
    results.update({filename: scores})
    best_thresh = cct_range[np.argmax(scores)]
    bests.update({filename: (best_thresh, np.max(scores))})

    Logger.info(f'{filename}: {best_thresh} ({np.max(scores)})')

results = pd.DataFrame(results, index=cct_range)
breakpoint()
    

    