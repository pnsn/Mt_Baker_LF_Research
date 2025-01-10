import logging, os, glob
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import normalized_mutual_info_score

from eqcorrscan.utils.stacking import align_traces


from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler



#### SETUP ####

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'workflow' / 'templates'

# USER INPUTS
min_members = 2
cct = 0.7
_f = TMPD/'step_6_clustered_UW.MBW..EHZ_5sec_cct0.45.tgz'




#### PROCESSING ###

print('loading')
ctr = ClusteringTribe().read(str(_f))
ctr.cct_regroup(0.45, inplace=True)
ctr.reindex_columns();
ctr._c.time = [pd.Timestamp(row.time) for _, row in ctr._c.iterrows()]
ctr.cct_regroup(corr_thresh=cct,inplace=True)
ctr.reindex_columns()


### PLOTTING ROUTINE ###

def plot_summary(
        ctr,
        axes=None,
        title='',
        render_legend=True,
        plotting_included=True,
        calculate_shifts=True,
        positive_shift=False):
    if axes is None:
        fig = plt.figure()
        gs = fig.add_gridspec(ncols=3,nrows=3)
        axes = [fig.add_subplot(gs[0,:2]),
                fig.add_subplot(gs[1,:2]),
                fig.add_subplot(gs[2,0]),
                fig.add_subplot(gs[2,1]),
                fig.add_subplot(gs[:2,2]),
                fig.add_subplot(gs[2,2])]
    
    if len(axes) != 6:
        raise AttributeError
    
    emap = {'eq': 'k*', 'lf':'r*', 'px':'mx', 'su':'bo'}
    xfmt = {'fmt':'b.','alpha': 0.1}
    # Plot Dendrogram
    if plotting_included:
        ctr.dendrogram(ax=axes[0], xlabels=['depth','etype'],scalar=[1e-3,None], title=f'{title}\n')
        axes[0].set_xlabel('Catalog Depth [km]')

    # Plot Scatters
    unit = {'time': '', 'depth': '[km]', 'longitude': '[$^o$E]', 'latitude': '[$^o$N]'}
    for _ax,_xf,_yf,_xs,_ys in [(1,'time','depth',False,1e-3),(2,'longitude','latitude',False,False),(3,'longitude','depth',False, 1e-3)]:
        # Overlay with everything in i_ctr
        if plotting_included:
            for _etype in ctr._c.etype.unique():
                j_df_c = ctr.clusters[ctr._c.etype == _etype]
                xv = j_df_c[_xf].copy().values
                if _xs:
                    xv *= _xs
                yv = j_df_c[_yf].copy().values
                if _ys:
                    yv *= _ys
                # PLOT
                axes[_ax].plot(xv, yv, emap[_etype], label=f'{_etype} (ct: {len(j_df_c)})', zorder=2)
        else:
            xv = ctr.clusters[_xf].copy().values
            if _xs:
                xv *= _xs
            yv = ctr.clusters[_yf].copy().values
            if _ys:
                yv *= _ys
            axes[_ax].plot(xv, yv, xfmt['fmt'], alpha=xfmt['alpha'], label='Unaffiliated', zorder=1)
    
        if _yf == 'depth':
            ylims = axes[_ax].get_ylim()
            if ylims[0] < ylims[1]:
                axes[_ax].set_ylim([ylims[1],ylims[0]])
        
        # Label Axes
        axes[_ax].set_xlabel(_xf[0].upper() + _xf[1:] + unit[_xf])
        axes[_ax].set_ylabel(_yf[0].upper() + _yf[1:] + unit[_yf])

    if plotting_included:
        # Get Shifts & Plot shifted traces
        tr_list = [tmp.st[0].copy() for tmp in ctr]
        if calculate_shifts:
            shifts_corrs = align_traces(tr_list, shift_len=50, positive=positive_shift)
            display_cc = True
        else:
            shifts_corrs = np.full(shape=(2,len(tr_list)), fill_value=1)
            display_cc = False


        for _e, tr in enumerate(tr_list):
            _tr = tr.copy()
            _tr.detrend('linear')
            _tr.normalize()
            _tr.normalize(norm=_tr.std()*3)
            _ccsign = np.sign(shifts_corrs[1][_e])
            axes[4].plot(_tr.times() + shifts_corrs[0][_e] - ctr[_e].prepick,
                            _ccsign*_tr.data + _e*2,
                            emap[ctr._c.etype[_e]][0],
                            lw=0.5)
        if display_cc:
            if positive_shift:
                vkwargs = {'vmin': 0, 'vmax': 1}
            else:
                vkwargs = {'vmin': -1, 'vmax': 1}
            ch = axes[4].scatter(-1*np.ones(len(tr_list)), np.arange(len(tr_list))*2, c=shifts_corrs[1], cmap='seismic', **vkwargs)
            plt.colorbar(ch, location='bottom', orientation='horizontal', fraction=0.05)

        axes[4].yaxis.tick_right()
        axes[4].set_ylabel('-- Towards Present Day -->')
        axes[4].set_yticks(np.arange(len(ctr))*2.,ctr._c.index)
        axes[4].tick_params('y',labelsize=6)
        if calculate_shifts:
            axes[4].set_xlabel('Time Relative to\nXCorr Adjusted Arrival [sec]')
        else:
            axes[4].set_xlabel('Time Relative to\nCatalog Arrivals [sec]')
        
        if render_legend:
            axes[1].legend(ncol=2)

    return axes


#### ITERATIVE PLOT GENERATION SECTION ###

grouped = []
orphans = []

        
# Get net/station set for this clustering tribe
sta_set = set()
for tmp in ctr:
    for tr in tmp.st:
        sta_set.add(f'{tr.stats.network}.{tr.stats.station}')

stastr = ''
for _sta in sta_set:
    stastr += f'{_sta}, '
stastr = stastr[:-2]


# Start by plotting everything
axes = plot_summary(ctr, title=f'Stations: {stastr}',
                    axes=None, plotting_included=True, calculate_shifts=False)


# Get alias to clusters dataframe
df_c = ctr._c.copy()
# Get group membership counts
g_counts = df_c.correlation_cluster.value_counts()

# Iterate across groups
for _gn, _ct in g_counts.items():
    # If insufficient membership, log EVID(s) to orphans
    if _ct < min_members:
        orphans += list(df_c[df_c.correlation_cluster==_gn].index.values)
    # If sufficient membership, log EVIDs to grouped
    else:
        # Get this group subset and add to grouped EVIDs list
        grouped_subset = df_c[df_c.correlation_cluster==_gn].index.values
        grouped += list(grouped_subset)
        # Get the clusering tribe of included events
        i_ctr = ctr.get_subset(grouped_subset)
        # Get the clustering tribe of excluded events
        x_ctr = ctr.get_subset(df_c[df_c.correlation_cluster!=_gn].index.values)

        # Initialize Figure and Plot excluded events
        axes = plot_summary(x_ctr, axes=None, plotting_included=False, render_legend=False)
        # Plot included events
        axes = plot_summary(i_ctr, axes=axes, title=f'Stations: {stastr} Group: {_gn:.0f}')


o_ctr = ctr.get_subset(orphans)
for _etype in o_ctr._c.etype.unique():
    e_orphans = o_ctr._c[o_ctr._c.etype == _etype].index.values
    oe_ctr = o_ctr.get_subset(e_orphans)
    xe_ctr = ctr.get_subset(list(e_orphans) + grouped)
    axes = plot_summary(xe_ctr, axes=None, plotting_included=False, render_legend=False)
    axes = plot_summary(oe_ctr, axes=axes, title=f'Stations: {stastr} Orphaned: {_etype}')


#         # Plot Dendrogram
#         i_ctr.dendrogram(ax=axes[0], xlabels=['depth','etype'],scalar=[1e-3,None], title=f'Stations: {stastr} Group {_gn:.0f}\n')
#         axes[0].set_xlabel('Catalog Depth [km]')

#         # Plot Scatters
#         for _ax,_xf,_yf,_xs,_ys in [(1,'time','depth',False,1e-3),(2,'longitude','latitude',False,False),(3,'longitude','depth',False, 1e-3)]:
            
#             # Plot everything not in i_ctr
#             xv = x_ctr._c[_xf].copy().values
#             if _xs:
#                 xv *= _xs
#             yv = x_ctr._c[_yf].copy().values
#             if _ys:
#                 yv *= _ys
#             # PLOT
#             axes[_ax].plot(xv, yv,'c.',alpha=0.1,zorder=1)

#             # Overlay with everything in i_ctr
#             for _etype in i_ctr._c.etype.unique():
#                 j_df_c = i_ctr.clusters[i_ctr._c.etype == _etype]
#                 xv = j_df_c[_xf].copy().values
#                 if _xs:
#                     xv *= _xs
#                 yv = j_df_c[_yf].copy().values
#                 if _ys:
#                     yv *= _ys
#                 # PLOT
#                 axes[_ax].plot(xv, yv, emap[_etype], label=f'{_etype} (ct: {len(j_df_c)})')
        
#             if _yf == 'depth':
#                 ylims = axes[_ax].get_ylim()
#                 axes[_ax].set_ylim([ylims[1],ylims[0]])
            
#             axes[_ax].set_xlabel(_xf)
#             axes[_ax].set_ylabel(_yf)
#             if _ax == 1:
#                 axes[_ax].legend(ncol=2,loc='lower left')
#         # axes[1].set_xlabel('UTC Date/Time')
#         # axes[1].set_ylabel('Catalog Depth [km]')
#         # ylims = axes[1].get_ylim()
#         # axes[1].set_ylim([ylims[1], ylims[0]])
#         # axes[1].legend(loc='lower left')

#         # Get Shifts & Plot shifted traces
#         tr_list = [tmp.st[0].copy() for tmp in i_ctr]
#         shifts_corrs = align_traces(tr_list, shift_len=50, positive=True)
#         # times = []
#         # positions = []
#         # amps = []

#         for _e, tr in enumerate(tr_list):
#             _tr = tr.copy()
#             _tr.detrend('linear')
#             _tr.normalize()
#             _tr.normalize(norm=_tr.std()*3)
#             # times.append(_tr.times() + shifts_corrs[0][_e] - i_ctr[_e].prepick)
#             # amps.append(_tr.data)
#             # positions.append(np.ones(_tr.count())*2*_e)
        
#         # times = np.array(times)
#         # amps = np.array(amps)
#         # positions = np.array(positions)
#         # axes[4].pcolor(times, positions, amps, cmap='seismic',vmin=-2, vmax=2)
#         #     # axes[4].scatter(_tr.times() + shifts_corrs[0][_e] - i_ctr[_e].prepick,
#         #     #                 _e*2*np.ones(_tr.count()),
#         #     #                 c=_tr.data,cmap=ecmap[i_ctr._c.etype[_e]],s=4)
#             axes[4].plot(_tr.times() + shifts_corrs[0][_e] - i_ctr[_e].prepick,
#                          _tr.data + _e*2,
#                          emap[i_ctr._c.etype[_e]][0],
#                          lw=0.5)
#         axes[4].set_title('Shifted Templates')
#         axes[4].yaxis.tick_right()
#         axes[4].set_ylabel('-- Towards Present Day -->')
#         axes[4].set_yticks(np.arange(len(i_ctr))*2.,i_ctr._c.index)
#         axes[4].tick_params('y',labelsize=6)
#         axes[4].set_xlabel('Time Relative to XCorr Adjusted Arrival [sec]')


# # Create Orphaned ClusteringTribe
# o_ctr = _ctr.get_subset(orphans)

# for _etype in o_ctr._c.etype.unique():
#     # Get list of Orphaned Templates of Given Event Type
#     oe_list = o_ctr._c[o_ctr._c.etype == _etype].index.values
#     i_ctr = _ctr.get_subset(oe_list)
#     # Get everything else that isnt in the orphaned events of a given event type
#     x_list = set(grouped).union(set(orphans).difference(set(oe_list)))
#     x_ctr = _ctr.get_subset(df_c[df_c.correlation_cluster!=_gn].index.values)
#     fig = plt.figure()
#     gs = fig.add_gridspec(ncols=3,nrows=3)
#     axes = [fig.add_subplot(gs[0,:2]), fig.add_subplot(gs[1,:2]), fig.add_subplot(gs[2,0]), fig.add_subplot(gs[2,1]),fig.add_subplot(gs[:,2]), ]
#     # Plot Dendrogram
#     i_ctr.dendrogram(ax=axes[0], xlabels=['depth','etype'],scalar=[1e-3,None], title=f'Stations: {stastr} Orphaned {_etype}\n')
#     axes[0].set_xlabel('Catalog Depth [km]')

#     # Plot Scatters
#     for _ax,_xf,_yf,_xs,_ys in [(1,'time','depth',False,1e-3),(2,'longitude','latitude',False,False),(3,'longitude','depth',False, 1e-3)]:
        
#         # Plot everything not in i_ctr
#         xv = x_ctr._c[_xf].copy().values
#         if _xs:
#             xv *= _xs
#         yv = x_ctr._c[_yf].copy().values
#         if _ys:
#             yv *= _ys
#         # PLOT
#         axes[_ax].plot(xv, yv,'c.',alpha=0.1,zorder=1)

#         # Overlay with everything in i_ctr
#         for _etype in i_ctr._c.etype.unique():
#             j_df_c = i_ctr.clusters[i_ctr._c.etype == _etype]
#             xv = j_df_c[_xf].copy().values
#             if _xs:
#                 xv *= _xs
#             yv = j_df_c[_yf].copy().values
#             if _ys:
#                 yv *= _ys
#             # PLOT
#             axes[_ax].plot(xv, yv, emap[_etype], label=f'{_etype} (ct: {len(j_df_c)})')
    
#         if _yf == 'depth':
#             ylims = axes[_ax].get_ylim()
#             axes[_ax].set_ylim([ylims[1],ylims[0]])
        
#         axes[_ax].set_xlabel(_xf)
#         axes[_ax].set_ylabel(_yf)

#         # Subplot Specific Formatting
#         if _ax == 1:
#             axes[_ax].legend(ncol=2,loc='lower left')

#     # Get Shifts & Plot shifted traces
#     tr_list = [tmp.st[0].copy() for tmp in i_ctr]
#     # tr_list = tr_list[1:]
#     # shifts_corrs = align_traces(tr_list, shift_len=250, positive=True)
#     for _e, tr in enumerate(tr_list):
#         _tr = tr.copy()
#         _tr.detrend('linear')
#         _tr.normalize()
#         _tr.normalize(norm=_tr.std()*3)
#         axes[4].plot(_tr.times() - i_ctr[_e].prepick,
#                         _tr.data + _e*2,
#                         emap[i_ctr._c.etype[_e]][0],
#                         lw=0.5)
#     axes[4].set_title('Unadjusted Templates')
#     axes[4].yaxis.tick_right()
#     axes[4].set_ylabel('-- Towards Present Day -->')
#     axes[4].set_yticks(np.arange(len(i_ctr))*2.,i_ctr._c.index)
#     axes[4].tick_params('y',labelsize=6)
#     axes[4].set_xlabel('Time Relative to Catalog Pick Time [sec]')



