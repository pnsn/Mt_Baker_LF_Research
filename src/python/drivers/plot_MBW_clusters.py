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
flist = glob.glob(str(TMPD/'corr_cluster_MBW_FreeStep.tgz'))
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
    _ctr.cct_regroup(0.45, inplace=True)
    _ctr.reindex_columns();
    tribes.update({_k: _ctr})


fig = plt.figure()
gs = fig.add_gridspec(nrows=2, ncols=3)
dax = fig.add_subplot(gs[0,:])
axes = [fig.add_subplot(gs[1,_e]) for _e in range(3)]

_ctr = tribes['MBW']
grps = _ctr.clusters.correlation_cluster.unique()
grps.sort()
_ctr.dendrogram(xlabels=['correlation_cluster','etype'],
                scalar=[False, False], ax=dax)
pairs = [('longitude','latitude'),('longitude','depth'),('latitude','depth')]
for _f, (x,y) in enumerate(pairs):
    xs = 1
    if y=='depth':
        ys = -1e-3
    else:
        ys = 1
    for _e in grps:
        _subset = _ctr.clusters[_ctr.clusters.correlation_cluster==_e]
        _zmed = _subset.depth.median()
        _vc= _subset.etype.value_counts()
        _top_etype = _vc.index[0]
        if len(_vc) == 1:
            _top_etype = f'{_top_etype}!'
        elif len(_vc) > 1:
            _top_etype = f'~{_top_etype}'
        if len(_subset) > 5:
            axes[_f].plot(_subset[x]*xs, _subset[y]*ys, '*',
                     label=f'grp{_e} {_top_etype} z :{_zmed*1e-3:.2f} km',
                     zorder=5, ms=6)
        elif len(_subset) > 1:
            axes[_f].scatter(_subset[x]*xs, _subset[y]*ys, s=len(_subset),
                     label=f'grp{_e} {_top_etype} z:{_zmed*1e-3:.2f} km',
                     zorder=10)
        else:
            axes[_f].plot(_subset[x]*xs, _subset[y]*ys, 'sk', alpha=0.25)
        axes[_f].set_xlabel(x)
        axes[_f].set_ylabel(y)
    axes[0].legend(bbox_to_anchor=[-0.05, 1], fontsize=8)

plt.show()