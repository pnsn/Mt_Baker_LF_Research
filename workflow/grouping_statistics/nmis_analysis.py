import logging
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import normalized_mutual_info_score

from eqcutil.util.logging import basic_logger_config

plt.close('all')

# Set up Logging
basic_logger_config(level=logging.INFO)
Logger = logging.getLogger()

# Map absolute paths
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow" / "templates" / "single_station"
# Template Save Directory
CTSD = PD_DIR / 'xcc_test'
# Grouping Summary Save File
INFILE = PD_DIR / 'xcc_threshold_sweep.csv'

df = pd.read_csv(str(INFILE))

holder = []
for cct in df.cc_threshold.unique():
    # Subset by correlation threshold
    idf = df[df.cc_threshold==cct]
    # Do pairwise comparisons of stations' groupings
    ustas = idf.station.unique()
    for ii, ista in enumerate(ustas):
        for jj, jsta in enumerate(ustas):
            if ii < jj:
                # Subset by correlation threshold and station code
                idf = df[(df.cc_threshold == cct)&(df.station==ista)]
                jdf = df[(df.cc_threshold == cct)&(df.station==jsta)]
                # Subset by mutually shared EVIDS
                ij_evids = set(idf.event_id.values).intersection(set(jdf.event_id.values))
                if len(ij_evids) == 0:
                    continue
                # Get groups
                igrp = idf[idf.event_id.isin(ij_evids)].xcc_group
                jgrp = jdf[jdf.event_id.isin(ij_evids)].xcc_group
                # Get counts of unique elements
                iunique = len(igrp.unique())
                junique = len(jgrp.unique())
                # Get event #unique/#events
                i_uer = len(ij_evids)/iunique
                j_uer = len(ij_evids)/junique
                # Set NMIS for completely ungrouped to 0
                if len(ij_evids) == iunique:
                    nmis = 0
                elif len(ij_evids) == junique:
                    nmis = 0
                # Set NMIS to 0 for completely grouped
                elif iunique == junique and iunique == 1:
                    nmis = 0
                else:
                    nmis = normalized_mutual_info_score(igrp, jgrp)

                if cct == 0.1 and nmis == 1:
                    breakpoint()
                result = [ista, jsta, len(ij_evids), iunique, junique, i_uer, j_uer, cct, nmis]
                holder.append(result)

df_out = pd.DataFrame(holder, columns=['istation','jstation','no_evids','i_grp_ct','j_grp_ct','i_uer','j_uer','cc_threshold','nmis'])




fig = plt.figure()
gs = fig.add_gridspec(nrows=3, ncols=1, wspace=0, hspace=0)
axes = [fig.add_subplot(gs[_e]) for _e in range(3)]
for _e, _fld in enumerate(['nmis','i_uer','j_uer']):
    pt = pd.pivot_table(df_out, values=_fld, index=['istation','jstation'], columns=['cc_threshold'])
    # Plot individual station pairs
    for staij, row in pt.iterrows():
        if _fld == 'nmis':
            axes[_e].plot(row.index, row.values, 'k-', alpha=0.2)
        else:
            axes[_e].semilogy(row.index, row.values, 'k-', alpha=0.2)

    # Plot IQR bounds and median
    quants = pt.quantile(q=[0.25, 0.5, 0.75], axis=0)
    for _q, row in quants.iterrows():
        if _q == 0.5:
            fmt = 'r--'
        else:
            fmt = 'r:'
        if _fld == 'nmis':
            axes[_e].plot(row.index, row.values, fmt, alpha=0.5)
        else:
            axes[_e].semilogy(row.index, row.values, fmt, alpha=0.5)


    # axes[_e].set_ylabel(_fld.upper())
axes[0].set_ylabel('Station Pair NMIS\n-- Increasing Grouping Similarity -->')
axes[1].set_ylabel('Station i Event Count / Group Count')
axes[2].set_ylabel('Station j Event Count / Group Count')
axes[-1].set_xlabel('Cross Correlation Threshold')

