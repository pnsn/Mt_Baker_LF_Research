
import logging
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from obsplus import EventBank
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

from eqcutil.util.logging import setup_terminal_logger

from hypothesis_utils import *

# Setup logger
Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Absolute path to repository root directory
ROOT = Path(__file__).parent.parent.parent
# Event Bank Base path
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Coherence distance table
COHD = ROOT / 'processed_data' / 'cluster' / 'tables' / 'coh_shift_table.csv'
# EVID ETYPE Reference table
EVETD = ROOT / 'data' / 'Events' / 'MtBaker_EVID_ETYPE.csv'
# Save Directory
SAVEDIR = ROOT / 'results' / 'tables'

# Define preferred / well-performing NSLC's
PREFNSLC = ['UW.JCW..EHZ',
            'UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.CMW..EHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'UW.PASS..BHZ','UW.PASS..HHZ',
            'CN.VDB..EHZ','CN.VDEB..HHZ',
            ]

CCT_PREF = 0.425

MIN_SRN = 2.

### PROCESSING SECTION ###
# Connect to event bank and read index
EBANK = EventBank(EBBP)
df_eb = EBANK.read_index()
# Assign UW EVIDs to index
df_eb.index = [_e.split('/')[-1].split('.')[0] for _e in df_eb.path]

# Read precomputed correlation coherence event-event-station table
df_coh = pd.read_csv(COHD)

# Read evid to etype mapping table
df_ee = pd.read_csv(EVETD, index_col=['uw_evid'])
# Attach etype to event bank index
df_eb = df_eb.join(df_ee, how='left')
# Create depth annotated etypes
aetypes = []
for _, row in df_eb.iterrows():
    if np.isfinite(row.vertical_uncertainty):
        if row.depth + row.vertical_uncertainty < 1e4:
            aetypes.append(f's{row.etype}')
        elif row.depth - row.vertical_uncertainty > 1e4:
            aetypes.append(f'd{row.etype}')
        else:
            aetypes.append(f'i{row.etype}')
    else:
        aetypes.append(f'q{row.etype}')

df_eb = df_eb.assign(aetype=aetypes)


breakpoint()
# Iterate across NSLC and reconstitute an array for each
coh_dict = {}
shift_dict = {}

# Reconstitute individual matrices in descending event-count order
for _k in df_coh.trace.value_counts().index:
    Logger.info(f'reconstituting {_k}')
    _df = df_coh[df_coh.trace==_k]
    coh_dict[_k] = get_symmetric(_df, k_field='coh', trace_value=1.)
    shift_dict[_k] = get_symmetric(_df, k_field='shift', trace_value=0.)



# DIMENSION 0: NSLC Catalog Labeling Subsetting
