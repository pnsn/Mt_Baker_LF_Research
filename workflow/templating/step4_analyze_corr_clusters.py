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
ROOT= Path(__file__).parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'workflow' / 'templates'

print('loading')
_f = TMPD/'step3_clustered_5sec_FreeStep.tgz'
_ctr = ClusteringTribe().read(str(_f))
_ctr.cct_regroup(0.45, inplace=True)
_ctr.reindex_columns();
_ctr._c.time = [pd.Timestamp(row.time) for _, row in _ctr._c.iterrows()]


def makeplot(_ctr, minmemb=4):

    fig = plt.figure()
    gs = fig.add_gridspec(nrows=5, ncols=4)
    dax = fig.add_subplot(gs[:2,:])
    tax = fig.add_subplot(gs[2,:-1])
    axes = [tax] + [fig.add_subplot(gs[3:,_e]) for _e in range(3)]

    grps = _ctr.clusters.correlation_cluster.unique()
    stas = set()
    for tmp in _ctr.templates:
        for tr in tmp.st:
            stas.add(tr.stats.station)
    tstas = ', '.join(list(stas))
    grps.sort()
    _ctr.dendrogram(xlabels=['correlation_cluster','etype'],
                    scalar=[False, False], ax=dax, title=f'{tstas}\n')
    daxlims = dax.get_xlim()
    dax.plot(daxlims,[1 - _ctr.cluster_kwargs['correlation_cluster']['corr_thresh']]*2, 'r:')
    # # Splotch to sync colorbars with 
    # dax.plot(-1,-1,'k.')
    dax.set_xlim(daxlims)
    pairs = [('time', 'depth'), ('longitude','latitude'),('longitude','depth'),('latitude','depth')]
    mcycle = ['s','^','.']

    for _f, (x,y) in enumerate(pairs):
        _c = 0
        _t = 0

        xs = 1
        if y=='depth':
            ys = -1e-3
        else:
            ys = 1
        for _e in grps:
            _subset = _ctr._c[_ctr._c.correlation_cluster==_e]
            if x != 'time':
                _x = _subset[x]*xs
            else:
                _x = _subset[x]
            _y = _subset[y]*ys
            _zmed = _subset.depth.median()
            _vc= _subset.etype.value_counts()
            _top_etype = _vc.index[0]
            if _top_etype == 'lf':
                _mkr = '*'
            elif _top_etype == 'eq':
                _mkr = '.'
            elif _top_etype == 'px':
                _mkr = 'x'
            elif _top_etype == 'su':
                _mkr = 's'
            else:
                _mkr = 'o'
        
            if len(_vc) == 1:
                _top_etype = f'{_top_etype}!'
                _mec = 'k'
            elif len(_vc) > 1:
                _top_etype = f'~{_top_etype}'
                _mec = 'c'
                
            # Plotting for large enough groups
            if len(_subset) >= minmemb:
                axes[_f].plot(_x, _y, _mkr,
                        label=f'grp{_e} {_top_etype} z :{_zmed*1e-3:.2f} km',
                        zorder=_e, ms=6, mec=_mec, mew=0.5)
            else:
                _et = _subset.etype.values[0]
                ms = 2
                if _et == 'lf':
                    ms = 6

                _mkr = _mkr + 'k'

                axes[_f].plot(_x, _y, _mkr,
                              alpha=0.25, zorder=_e, ms=ms)
            axes[_f].set_xlabel(x)
            axes[_f].set_ylabel(y)
        axes[-1].legend(bbox_to_anchor=[1.05, 1], ncols=2, fontsize=8)
    return fig, [dax] + axes



if __name__ == 'main':

    fig, axes = makeplot(_ctr)

    plt.show()