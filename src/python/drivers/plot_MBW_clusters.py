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

print('loading')
_f = TMPD/'corr_cluster_XSTA_FreeStep.tgz'
_ctr = ClusteringTribe().read(str(_f))
_ctr.cct_regroup(0.45, inplace=True)
_ctr.reindex_columns();
_ctr._c = [pd.Timestamp(row.time) for _, row in _ctr._c.iterrows()]


def makeplot(_ctr):

    fig = plt.figure()
    gs = fig.add_gridspec(nrows=3, ncols=4)
    dax = fig.add_subplot(gs[0,:])
    tax = fig.add_subplot(gs[1,:])
    axes = [tax] + [fig.add_subplot(gs[2,_e]) for _e in range(3)]

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

            # if x == 'latitude' and y == 'depth':
            #     axes[_f].errorbar(_subset[x]*xs, _subset[y]*ys,
            #                       xerr=_subset['horizontal_uncertainty']/111.2,
            #                       yerr=_subset['vertical_uncertainty']*ys,
            #                       capsize=2, color='k', alpha=0.25)
            if len(_vc) == 1:
                _top_etype = f'{_top_etype}!'
                _mec = 'k'
            elif len(_vc) > 1:
                _top_etype = f'~{_top_etype}'
                _mec = 'c'
            if len(_subset) > 5:
                axes[_f].plot(_x, _y, '*',
                        label=f'grp{_e} {_top_etype} z :{_zmed*1e-3:.2f} km',
                        zorder=_e, ms=6, mec=_mec, mew=0.5)
            elif len(_subset) > 1:
                # axes[_f].scatter(_x, _y, s=len(_subset)*2,
                #          label=f'grp{_e} {_top_etype} z:{_zmed*1e-3:.2f} km',
                #          zorder=10)
                axes[_f].plot(_x, _y, mcycle[_c],
                            label=f'grp{_e} {_top_etype} z:{_zmed*1e-3:.2f} km',
                            zorder=_e, ms=4, mec=_mec, mew=0.5)
                _t += 1
                if _t%10 == 0:
                    _c +=1
            else:
                _et = _subset.etype.values[0]
                ms = 2
                if _et == 'lf':
                    mkr = '*k'
                    ms = 6
                elif _et == 'eq':
                    mkr = 'ok'
                elif _et == 'px':
                    mkr = '.k'
                elif _et == 'su':
                    mkr = 'sc'
                else:
                    mkr = 'm.'
                axes[_f].plot(_x, _y, mkr,
                              alpha=0.25, zorder=_e, ms=ms)
            axes[_f].set_xlabel(x)
            axes[_f].set_ylabel(y)
        axes[-1].legend(bbox_to_anchor=[1.05, 1], ncols=2, fontsize=8)
    return fig, [dax] + axes



if __name__ == 'main':

    fig, axes = makeplot(_ctr)

    plt.show()