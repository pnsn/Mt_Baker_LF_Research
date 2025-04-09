from pathlib import Path
import numpy as np

import matplotlib.pyplot as plt

from map_util import pnsn_pallet

ROOT = Path(__file__).parent.parent.parent


cij1 = np.array([[1, 0.9, 0.5],
                 [0.9, 1, 0.2],
                 [0.5, 0.2, 1]])

cij2 = np.array([[1, 0.7, 0.3],
                 [0.7, 1, 0.1],
                 [0.3, 0.1, 1]])

cij3 = np.array([[1, 0.6, np.nan],
                 [0.6, 1, np.nan],
                 [np.nan, np.nan, np.nan]])


def plot_cijk(ax, cijk, events, cmap='Greens'):
    ax.pcolor(cijk, cmap=cmap, vmin=0, vmax=1)
    for _i in range(len(cijk)):
        for _j in range(len(cijk)):
            if not np.isfinite(cijk[_i,_j]):
                color = pnsn_pallet()['lime']
            elif _i == _j:
                color = 'w'
            else:
                color = 'k'
            ax.text(_i+0.5, _j+0.5, f'{cijk[_i, _j]:.2f}',
                    ha='center',va='center', fontsize=12,
                    fontweight='extra bold', color=color)
            

    ax.set_xticks(np.arange(len(cijk)) + 0.5, labels=events)
    ax.set_yticks(np.arange(len(cijk)) + 0.5, labels=events)

fig = plt.figure(figsize=(9,9))
gs = fig.add_gridspec(ncols=3, nrows=3)
axes = [fig.add_subplot(gs[_e]) for _e in range(9)]


# Simple averaging
plot_cijk(axes[0],cij1, ['E0','E1','E2'])
plot_cijk(axes[1],cij2, ['E0','E1','E2'])
plot_cijk(axes[2],np.mean(np.stack([cij1, cij2], axis=-1), axis=-1),
          ['E0','E1','E2'])

plot_cijk(axes[3],cij1, ['E0','E1','E2'])
plot_cijk(axes[4],cij3, ['E0','E1','E2'])
plot_cijk(axes[5],np.nanmean(np.stack([cij1, cij3], axis=-1), axis=-1),
          ['E0','E1','E2'])


plot_cijk(axes[6],cij1, ['E0','E1','E2'])
plot_cijk(axes[7],cij2, ['E0','E1','E3'])

cij4 = np.full(shape=(4,4), fill_value=np.nan)
cij4[:2, :2] = np.mean(np.stack([cij1[:2,:2], cij2[:2,:2]], axis=-1), axis=-1)
cij4[2,:] = cij1[2,:]
cij4[:,2] = cij1[:,2]
cij4[]
for _i in range(4):
    for _j in range(4):
        

plot_cijk(axes[8],np.nanmean(np.stack([cij1, cij3], axis=-1), axis=-1),
          ['E0','E1','E2'])

plt.show()


