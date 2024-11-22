"""
:module: src/python/drivers/presentation_cct_sweep.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: Trimmed-down version of the `postprocess_mbs_clusters.py` for
    a presentation on 13. NOV 2024 to the PNSN SP Team meeting.

    This driver sweeps across a set of cross-correlation threshold (CCT)
    values and generates figure files for the purely correlation based
    clustering of 
"""
import os, glob, logging
from pathlib import Path

import numpy as np
import pandas as pd
from obsplus import EventBank
import matplotlib.pyplot as plt

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger

# Define objective function
# def objective(cct, etypes, )

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
ctr_files = glob.glob(os.path.join(ctr_file_path, 'fv_mean_sl_10.0_ls_False.tgz'))

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

cct_range = np.linspace(0.025, 1, 21)

filepath = ctr_files[0]
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

ctr2 = ctr.copy()
ctr2.select_template_traces(station='MBW*')
meth = 'correlation_cluster'
ctr2.cluster(meth, **ctr2.cluster_kwargs[meth])
breakpoint()
    # Set clustering to specified threshold

# Create Dendrogram For everything
fig = plt.figure(figsize=(9,4))
ax = fig.add_subplot(111)
ctr.dendrogram(corr_thresh=0.6, ax=ax, title='All Events\n')


# Create Dendrogram for LF
fig = plt.figure(figsize=(9,4))
ax = fig.add_subplot(111)
lf_ctr = ctr.get_subset(ctr.clusters[ctr.clusters.etype == 'lf'].index)
lf_ctr.dendrogram(xlabelfield='etype', corr_thresh=0.6, ax=ax, title='All LF Events\n')

# Create Dendrogram for EQ
fig = plt.figure(figsize=(9,4))
ax = fig.add_subplot(111)
eq_ctr = ctr.get_subset(ctr.clusters[ctr.clusters.etype == 'eq'].index)
eq_ctr.dendrogram(xlabelfield='index', corr_thresh=0.6, ax=ax, title='All EQ Events\n')

# Create Dendrogram for LF with depth - 2*sdepth > 0
fig = plt.figure(figsize=(9,4))
ax = fig.add_subplot(111)
dlf_filt = (ctr.clusters.etype=='lf') & (ctr.clusters.depth - ctr.clusters.vertical_uncertainty*2 > 0)
dlf_ctr = ctr.get_subset(ctr.clusters[dlf_filt].index)
dlf_ctr.dendrogram(xlabelfield='depth', corr_thresh=0.6, ax=ax, title='Subsurface (depth - 2*sdep > 0) LF Events\n')


fig = plt.figure(figsize=(9,4))
ax = fig.add_subplot(111)
dlf_filt = (ctr.clusters.etype=='lf') & (ctr.clusters.depth - ctr.clusters.vertical_uncertainty*2 < 0)
slf_ctr = ctr.get_subset(ctr.clusters[dlf_filt].index)
slf_ctr.dendrogram(xlabelfield='depth', corr_thresh=0.6, ax=ax, title='Shallow (depth - 2*sdep < 0) LF Events\n')

plt.show()

    




    # Logger.info(f'Running scoring for {len(cct_range)} cross_correlation threshold values')
    # scores = []
    # for _e, thresh in enumerate(cct_range):
    #     Logger.debug(f'{_e+1} of {len(cct_range)}')
    #     scores.append(evaluate_corr_thresh(ctr, thresh, **kwargs))
    # results.update({filename: scores})
    # best_thresh = cct_range[np.argmax(scores)]
    # bests.update({filename: (best_thresh, np.max(scores))})
