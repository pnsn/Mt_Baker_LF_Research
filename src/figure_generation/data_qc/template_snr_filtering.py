import logging, os, glob
from pathlib import Path

import matplotlib.pyplot as plt

import numpy as np
import pandas as pd

from eqcutil import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger

Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Absolute path to repo root directory
ROOT = Path(__file__).parent.parent.parent.parent
# Preclustered single-station clusters
CTRD = ROOT / 'processed_data' / 'cluster' / 'single_station'

clist = glob.glob(str(CTRD/'*.tgz'))
ctrd = {}
df_snr = pd.DataFrame()
df_eval = pd.DataFrame()
# Load all single-station, pre-clustered templates
for _c in clist:
    name = os.path.splitext(_c.split('/')[-1])[0]
    nslc = name.split('_')[0]
    Logger.info(f'Loading {_c} as key {nslc}')
    _ctr = ClusteringTribe().read(_c)
    ctrd.update({nslc: _ctr})

    _ser = _ctr._c.mean_snr_dB
    _ser.name = nslc
    df_snr = pd.concat([df_snr, _ser], axis=1, ignore_index=False)
    _ser = _ctr._c.pick_status
    _ser.name = nslc
    df_eval = pd.concat([df_eval, _ser], axis=1, ignore_index=False)

# Sort columns and index
df_snr = df_snr.sort_index().T.sort_index().T

fig = plt.figure(figsize=(6,6))
ax1 = fig.add_subplot(211)
ch = ax1.pcolor(df_snr.T, cmap='tab20c', vmin=0, vmax=10)
cbh = fig.colorbar(ch, orientation='horizontal', fraction=0.05)
cbh.set_label('Average Pick-Centered SNR [dB]')
ax1.set_yticks(np.array(range(len(df_snr.columns))) + 0.5, df_snr.columns)
ax1.set_xticks(np.array(range(len(df_snr.index))) + 0.5, df_snr.index)
ax1.set_ylabel('Channel Code')
ax1.set_xlabel('Event ID [ - ]')
ax2 = fig.add_subplot(212, sharex=ax1)
for _e in [1,2,5,8]:
    ax2.fill_between(df_snr.index, [0]*len(df_snr), df_snr[(df_snr >= _e)].count(axis=1), label=f'SNR $\geq$ {_e}')
ax2.legend()
ax2.set_ylabel('Included Tempates [Count]')
ax2.set_xlabel('Event ID [ - ]')

plt.show()