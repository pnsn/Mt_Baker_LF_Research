"""
:purpose:
The goal of this analysis is to quantify in the range [0, 1] how similar pairwise
combinations of event classifications are for the current PNSN catalog (AQMS),

These data have a mix of categorical (event type interpretation and cross correlation
grouping) and qualitative (event location solution, correlation threshold) elements.

As such, we can use a grid-search approach to optimize the consistency of groupings
wherein we iterate across the radius with which we define event clusters and
the correlation coefficient with which we threshold correlation base groupings.


"""

import os, glob
from pathlib import Path

from obspy import Catalog
from eqcorrscan.utils.clustering import catalog_cluster

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score, adjusted_mutual_info_score

from eqcutil import ClusteringTribe

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# Reviewer / AQMS classes
REVD = ROOT / 'processed_data' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Preclustered template root
TDIR = ROOT / 'results' / 'cluster' / 'single_station'

# Read reviewed event results
df_rev = pd.read_csv(REVD, index_col=[0])
r_evids = set(df_rev.index.values)
# Get file names of pre-clustered ClusteringTribes
files = glob.glob(str(TDIR/'*'))
for _f in files:
    # Load pre-clustered tribe
    ctr = ClusteringTribe().read(_f)
    # Get intersection of EVIDs for CTR (what data passed) and reviewed events
    c_evids = set(ctr._c.index.values)
    q_set = c_evids.intersection(r_evids)
    # Subset clustering tribe by events included in the review
    ctr.get_subset(list(q_set))
    # Create catalog from included templates
    cat = Catalog([t.event for t in ctr])
    # Run spatial clustering on catalog
    groups = catalog_cluster(cat, metric='distance', thresh=6, show=False)
    # reconstitute evid-groupno pairs
    pairs = []
    for _e, _sc in enumerate(groups):
        for _event in _sc:
            eparts = _event.resource_id.id.split('/')[-2:]
            evid = eparts[0].lower() + eparts[1]
            grp = _e
            pairs.append((evid, grp))
    # ctr.cluster(method='space_cluster', threshold=3000)
    breakpoint()
