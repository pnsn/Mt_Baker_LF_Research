import logging, os, glob
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import normalized_mutual_info_score

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'bank_templates' / 'd20km'
DATA = TMPD / 'multi_thresholded_clusters.csv'
# Load Summary
df = pd.read_csv(DATA, index_col=[0], parse_dates=['time','creation_time'])

# Load tribes
flist = glob.glob(str(TMPD/'corr_cluster_*_FreeStep.tgz'))
tribes = {}
for _f in flist:
    if 'XSTA' in _f:
        continue
    elif 'Composite' in _f:
        continue
    path, fname = os.path.split(_f)
    parts = fname.split('_')
    _k = parts[2]
    Logger.info(f'loading {_f}')
    _ctr = ClusteringTribe().read(_f)
    tribes.update({_k: _ctr})
# cols = []
# for _cct in range(1,10):
#     cols += list(df.filter(like=f'{_cct:02d}').columns)

# reindex_columns(df, cols, inplace=True)

# Get a count of the number of independent clusters observing a given evid
mutual_obs_index = df.filter(like='_07').notna().sum(axis=1).sort_values(ascending=False)
# Get a count of the number of linkages for events (number of independently observing clusters)
mutual_obs_count = mutual_obs_index.value_counts().sort_index()

## ASSESS PAIRWISE NMIS
# Iterate across threshold level
# "How similar is the groupings between stations for a given correlation threshold?"

pairwise_nmis = []
unique_pairs = []
pair_counts = []
cct_vals = []
sta_pairs = []
for _thr in range(1, 10):
    _df = df.filter(like=f'_{_thr:02d}')
    cct_vals.append(_thr)
    # Collector for Normalized Mutual Info Scores
    np_line = []   
    # Collector for Unique Group Pairings
    up_line = []
    # Collector for number of shared events
    pc_line = []

    # Iterate across "true" labels
    for _p, pcol in enumerate(_df.columns):
        # Iterate across "modeled" labels
        for _s, scol in enumerate(_df.columns):
            # Only do upper triangle
            if _p > _s:
                # Subset to mutually represented entries
                __df = _df[(_df[pcol].notna()) & (_df[scol].notna())][[pcol, scol]]
                _nmis = normalized_mutual_info_score(__df[pcol].values, __df[scol].values)
                np_line.append(_nmis)
                up_line.append(len(__df.value_counts()))
                pc_line.append(len(__df))
                if (pcol[:-3], scol[:-3]) not in sta_pairs:
                    sta_pairs.append((pcol[:-3], scol[:-3]))
    pairwise_nmis.append(np_line)
    unique_pairs.append(up_line)
    pair_counts.append(pc_line)

df_nmis = pd.DataFrame(data=pairwise_nmis, index=cct_vals, columns=sta_pairs).T
df_up = pd.DataFrame(data=unique_pairs, index=cct_vals, columns=sta_pairs).T
df_pc = pd.DataFrame(data=pair_counts, index=cct_vals, columns=sta_pairs).T

### RUN NMIS ANALYSIS FOR EACH TEMPLATE
# "How sensitive is grouping for each station to the correlation threshold level?"
tribe_nmis = []
tribe_up = []
tribe_pc = []
cct_pairs = []
sta_names = []
for _name in tribes.keys():
    _df = df.filter(like=_name)
    _df = _df[_df[f'{_name}_07'].notna()]
    np_line = []
    up_line = []
    pc_line = []
    sta_names.append(_name)
    for _ii in range(1,10):
        for _jj in range(1,10):
            if _jj > _ii:
                idf = _df[f'{_name}_{_ii:02d}']
                jdf = _df[f'{_name}_{_jj:02d}']
                _nmis = normalized_mutual_info_score(idf.values, jdf.values)
                np_line.append(_nmis)
                up_line.append(len(idf.value_counts()))
                pc_line.append(len(idf))
                if (_ii, _jj) not in cct_pairs:
                    cct_pairs.append((_ii,_jj))
    tribe_nmis.append(np_line)
    tribe_up.append(up_line)
    tribe_pc.append(pc_line)

df_tribe_nmis = pd.DataFrame(tribe_nmis, index=sta_names, columns=cct_pairs).T
df_tribe_up = pd.DataFrame(tribe_up, index=sta_names, columns=cct_pairs).T
df_tribe_pc = pd.DataFrame(tribe_pc, index=sta_names, columns=cct_pairs).T

# Iterate over event ID's and look at station pairs' classification consistency
# "How stable is the grouping of an event across station pairs for identical correlation thresholds?"

event_nmis = []
event_up = []
event_pc = []
evids = []
sta_pairs = []
for evid, row in df.iterrows():
    np_line = []
    up_line = []
    pc_line = []
    evids.append(evid)
    for _ii, _iname in enumerate(tribes.keys()):
        for _jj, _jname in enumerate(tribes.keys()):
            # Upper triangle only
            if _jj > _ii:
                idf = row.filter(like=f'{_iname}_')
                jdf = row.filter(like=f'{_jname}_')
                if (_iname, _jname) not in sta_pairs:
                    sta_pairs.append((_iname,_jname))
                # Only if pairs match for this event
                if all(idf.notna()) and all(jdf.notna()):
                    try:
                        _nmis = normalized_mutual_info_score(idf.values, jdf.values)
                    except:
                        breakpoint()
                    np_line.append(_nmis)
                    up_line.append(len(idf.value_counts()))
                    pc_line.append(len(idf))
                else:
                    np_line.append(np.nan)
                    up_line.append(np.nan)
                    pc_line.append(np.nan)
    event_nmis.append(np_line)
    event_up.append(up_line)
    event_pc.append(pc_line)
# The lower the value of NMIS in this array, the less consistent grouping is
# for a given 
df_event_nmis = pd.DataFrame(event_nmis, index=evids, columns=sta_pairs).T.astype(pd.SparseDtype("float", np.nan))
df_event_up = pd.DataFrame(event_up, index=evids, columns=sta_pairs).T.astype(pd.SparseDtype("int", np.nan))
df_event_pc = pd.DataFrame(event_pc, index=evids, columns=sta_pairs).T.astype(pd.SparseDtype("int", np.nan))