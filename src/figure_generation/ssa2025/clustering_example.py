from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram

ROOT = Path(__file__).parent.parent.parent.parent
SDIR = ROOT/'results'/'figures'/'SSA2025'

isshow = True
issave = True
DPI = 120
FMT = 'png'
plt.rcParams.update({'font.size': 14})
sfkw = {'dpi':DPI, 'format':FMT, 'bbox_inches':'tight'}

#124
coh0 = np.array(
    [[1, 0.7, 0.6],
     [0.6, 1, 0.3],
     [0.7, 0.3, 1]])
#234
coh1 = np.array(
    [[1, 0.4, 0.3],
     [0.4, 1, 0.9],
     [0.3, 0.9, 1]])

coh3 = np.array(
    [[1., .7, .6, 0., 0.],
     [.7, 1., .3, .4, .3],
     [.6, .3, 1., 0., 0.],
     [0., .4, 0., 1., .9],
     [0., .3, 0., .9, 1.]])
#     [[1.0,0.9,0.0,0.4],
#      [0.9,1.0,0.8,0.45],
#      [0.0,0.8,1.0,0.7],
#      [0.4,0.45,0.7,1.0]]
# )

model = AgglomerativeClustering(
    linkage='average',
    distance_threshold=0.5,
    n_clusters=None,
    metric='precomputed'
)

def model2linkmat(model):
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for _e, merge in enumerate(model.children_):
        cc = 0
        for child_idx in merge:
            if child_idx < n_samples:
                cc += 1
            else:
                cc += counts[child_idx - n_samples]
        counts[_e] = cc
    linkmat = np.column_stack(
        [
            model.children_,
            model.distances_,
            counts
        ]
    ).astype(float)
    return linkmat


def label_axes(ax, x, xloc='top',yloc='left',labels=['E0','E1','E2']):
    ax.xaxis.set_ticks_position(xloc)
    ax.xaxis.set_label_position(xloc)
    ax.yaxis.set_ticks_position(yloc)
    ax.yaxis.set_label_position(yloc)
    ax.set_xticks(np.arange(0.5, len(labels) + 0.5),
                labels=labels)
    ax.set_yticks(np.arange(0.5, len(labels) + 0.5),
                labels=labels)

    for _i in range(x.shape[0]):
        for _j in range(x.shape[1]):
            ax.text(_i+0.5, _j+0.5,
                    f'{x[_i,_j]:.1f}',
                    ha='center',
                    va='center')
    return ax

def plot_dendrogram(model, ax, coh, labels, threshold=0.5, above_threshold_color='xkcd:ugly brown'):
    linkmat = model2linkmat(model.fit(1 - coh))
    dend = dendrogram(linkmat, ax=ax, color_threshold=threshold,
                      above_threshold_color=above_threshold_color,
                      labels=labels)
    xlims = ax.get_xlim()
    ax.plot(xlims, [threshold]*2, ':', color='dodgerblue')
    ax.set_xlim(xlims)


    return (linkmat, dend)

### SINGLE COHERENCE MATRIX EXAMPLE
fig = plt.figure(figsize=(9, 7))
gs = fig.add_gridspec(ncols=2, nrows=4, wspace=0.5)
axc = fig.add_subplot(gs[1:3,0])
axc.pcolor(coh0, cmap='Greens', vmin=0., vmax=1.5)
axc = label_axes(axc, coh0)

axd = fig.add_subplot(gs[:,1])
linkmat, dend = plot_dendrogram(
    model, axd, coh0,
    labels=['E0','E1','E2'])

ylims = axd.get_ylim()
axd.set_yticks(np.linspace(0,1,11), np.linspace(1,0,11))
axd.set_ylabel('Average Coherence [ - ]')
axd.set_ylim(ylims)
axd.text(np.mean(axd.get_xlim()) - 3, .51, 'Threshold', ha='center', va='bottom')
axc.set_xlabel('Observed Events [E#]')
axc.set_ylabel('Observed Events [E#]')
axc.set_title('Station-Channel 0')

if issave:
    fig.savefig(str(SDIR/f'clustering_example_part_1_{DPI}dpi.{FMT}'), **sfkw)


## HETEROGENEOUS DATA EXAMPLE
fig = plt.figure(figsize=(15,6))
gs = fig.add_gridspec(ncols=3, nrows=4, wspace=0.5)
ax0 = fig.add_subplot(gs[:2,0])
ax1 = fig.add_subplot(gs[2:,0])
axe = fig.add_subplot(gs[1:3,1])
axd = fig.add_subplot(gs[:,-1])

ax0.pcolor(coh0, cmap='Greens', vmin=0, vmax=1.5)
ax0 = label_axes(ax0, coh0, labels=['E0','E1','E2'])
ax1.pcolor(coh1, cmap='Greens', vmin=0, vmax=1.5)
ax1 = label_axes(ax1, coh1, labels=['E1','E3','E4'], xloc='bottom')
axe.pcolor(coh3, cmap='Greens', vmin=0, vmax=1.5)
axe = label_axes(axe, coh3, labels=['E0','E1','E2','E3','E4'])

(dend, linkmat) = plot_dendrogram(
    model, axd, coh3, labels=['E0','E1','E2','E3','E4'])

linkmat = model2linkmat(model.fit(1-coh3))
dend = dendrogram(linkmat, ax=axd, color_threshold=0.5,
                  above_threshold_color='xkcd:ugly brown',
                  labels=[f'E{_d}' for _d in range(5)])

# Window dresssing
ax0.set_xlabel('Station-Channel 0')
# ax0.set_xlabel('Observed Events [E#]')
ax1.set_xlabel('Station-Channel 1')
# ax1.set_xlabel('Observed Events [E#]')

ylims = axd.get_ylim()
axd.set_yticks([0,0.2,0.4, 0.6, 0.8,1], [1.0, 0.8, 0.6, 0.4, 0.2, 0.0])
axd.set_ylabel('Average Coherence [ - ]')
axd.set_ylim(ylims)
axd.text(np.mean(axd.get_xlim()), .51, 'Threshold', ha='center', va='bottom')

axe.set_xlabel('Combined\nCoherence Matrix')
axd.set_xticks(axd.xaxis.get_ticklocs(), axd.xaxis.get_ticklabels())

if issave:
    fig.savefig(str(SDIR/f'clustering_example_part_2_{DPI}dpi.{FMT}'), **sfkw)

if isshow:
    plt.show()