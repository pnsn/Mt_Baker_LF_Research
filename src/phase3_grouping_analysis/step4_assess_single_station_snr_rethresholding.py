import logging, glob, os
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

from eqcutil import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger

from hypothesis_utils import *

# Setup logger
Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# ClusteringTribe Directory
CTRD = ROOT / 'processed_data' / 'cluster' / 'single_station'


# # Event Space/Time cluster table
# XTCD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'dist_table.csv'
# # Coherence distance table
# COHD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'coh_shift_table.csv'
# EVID ETYPE Reference table
EVETD = ROOT / 'data' / 'Events' / 'MtBaker_EVID_ETYPE.csv'
# Save Directory
SAVEDIR = ROOT / 'results' / 'tables'


flist = glob.glob(str(CTRD / '*_clustered.tgz'))
ctr_dict = {}
# Iterate across 
for _f in flist:
    _, namex = os.path.split(_f)
    name, ext = os.path.splitext(namex)
    nslc = name.split('_')[0]
    Logger.info(f'reading {nslc}')
    ctr_dict.update({nslc: ClusteringTribe().read(_f)})

