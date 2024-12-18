import logging, os
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

from sklearn.metrics import normalized_mutual_info_score

from eqcutil.util.pandas import reindex_columns
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
DATA = ROOT / 'processed_data' / 'bank_templates' / 'd20km' / 'multi_thresholded_clusters.csv'

df = pd.read_csv(DATA, index_col=[0], parse_dates=['time','creation_time'])

cols = []
for _cct in range(1,10):
    cols += list(df.filter(like=f'{_cct:02d}').columns)

reindex_columns(df, cols, inplace=True)

# Get a count of the number of independent clusters observing a given evid
mutual_obs_index = df.filter(like='_07').notna().sum(axis=1).sort_values(ascending=False)
# Get a count of the number of linkages for events (number of independently observing clusters)
mutual_obs_count = mutual_obs_index.value_counts().sort_index()

## ASSESS PAIRWISE NMIS
# Iterate across threshold level

pairwise_nmis = []
unique_pairs = []
pair_counts = []
row_idx = []
col_idx = []
for _thr in range(1, 10):
    _df = df.filter(like=f'_{_thr:02d}')
    row_idx.append(_thr)
    # Iterate across "true" labels
    np_line = []   
    up_line = []
    pc_line = []
    for _p, pcol in enumerate(_df.columns):
        for _s, scol in enumerate(_df.columns):
            if _p > _s:
                # Subset to mutually represented entries
                __df = _df[(_df[pcol].notna()) & (_df[scol].notna())][[pcol, scol]]
                _nmis = normalized_mutual_info_score(__df[pcol].values, __df[scol].values)
                np_line.append(_nmis)
                up_line.append(len(__df.value_counts()))
                pc_line.append(len(__df))
                if (pcol[:-3], scol[:-3]) not in col_idx:
                    col_idx.append((pcol[:-3], scol[:-3]))
    pairwise_nmis.append(np_line)
    unique_pairs.append(up_line)
    pair_counts.append(pc_line)

df_nmis = pd.DataFrame(data=pairwise_nmis, index=row_idx, columns=col_idx)
df_up = pd.DataFrame(data=unique_pairs, index=row_idx, columns=col_idx)
df_pc = pd.DataFrame(data=pair_counts, index=row_idx, columns=col_idx)

