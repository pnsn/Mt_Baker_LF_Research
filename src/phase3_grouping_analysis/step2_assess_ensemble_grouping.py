import logging
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



# Read precomputed AQMS hypocenter distance table
df_dxt = pd.read_csv(XTCD)
# Read precomputed correlation coherence event-event-station table
df_coh = pd.read_csv(COHD)

# Read evid to etype mapping table
df_ee = pd.read_csv(EVETD, index_col=['uw_evid'])

# Iterate across NSLC and reconstitute an array for each
coh_dict = {}
shift_dict = {}

for _k in df_coh.trace.unique():
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

# Run "by eye" preferred clustering threshold
model_eyeball = AgglomerativeClustering(
    n_clusters=None,
    linkage='single',
    metric='precomputed',
    distance_threshold=0.55).fit(X)
Z = get_linkage_matrix(model_eyeball)
fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(111)
outs = dendrogram(Z, color_threshold=0.55, labels=df_ee.loc[coh_merge.index.values].etype.values)

# Iterate across thresholds and test against consistency of etype labeling
for cthr in np.arange(.3, .8, 0.01):
    cmod = AgglomerativeClustering(
        distance_threshold=cthr,
        n_clusters=None,
        metric='precomputed',
        linkage='single').fit(X)
    _ser_ct = pd.Series(cmod.labels_, index=coh_merge.index, name=f'cc{cthr:.2f}')


plt.show()