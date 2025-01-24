# Imports
import logging, glob, os

from pathlib import Path

import pandas as pd
import numpy as np

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import basic_logger_config

basic_logger_config(level=logging.INFO)
# Map absolute paths
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow" / "templates" / "single_station"
# Template Save Directory
CTSD = PD_DIR / 'xcc_test'
# Grouping Summary Save File
OUTFILE = PD_DIR / 'xcc_threshold_sweep.csv'

Logger = logging.getLogger(__file__)


cct_vector = np.arange(0.1,1.,0.025)


# Get file names for templates
flist = glob.glob(str(CTSD/'*.tgz'))

holder = []
# Iterate across stations
for file in flist:
    Logger.info(f'loading: {file}')
    ctr = ClusteringTribe().read(file)
    sta = ctr[0].st[0].stats.station
    # Iterate across cross correlation threshold values
    for cct in cct_vector:
        Logger.info(f'Threshold Value: {cct}')
        # Apply new cross correlation threshold
        ctr.cct_regroup(corr_thresh=cct, inplace=True)
        # Reindex
        ctr.reindex_columns()
        # Accumulate data
        for event_id, row in ctr._c.iterrows():
            line = [sta, event_id, row.xcc, cct, row.eval_mode]
            holder.append(line)

# Create dataframe
df = pd.DataFrame(data=holder, columns=['station','event_id','xcc_group','cc_threshold','pick_evaluation_mode'])
# Write to disk
df.to_csv(str(OUTFILE), header=True, index=False)
breakpoint()

    # "Peel the onion" - apply high CC thresholds to get at sub-groups
