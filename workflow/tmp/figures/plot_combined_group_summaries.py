import logging, os, glob, warnings
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from obspy.clients.fdsn import Client

from eqcorrscan.utils.stacking import align_traces

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Ignore FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)

#### SETUP ####

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'workflow' / 'templates' / 'single_station' / 'xcc_test'

##### USER INPUTS #####
min_members = 4
file = TMPD/'JCW.tgz'
cct = 0.45
##### END OF USER INPUT SECTION #####

#---------------------------------------------------------------------------#

#### METHODS ####

### PLOTTING ROUTINE ###
def plot_summary(
        ctr,
        axes=None,
        title='',
        render_legend=True,
        plotting_grouped=True,
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
                # fig.add_subplot(gs[2,2])]
    
    if len(axes) != 6:
        raise AttributeError
    
    emap = {'eq': 'kv', 'lf':'ro', 'px':'gx', 'su':'b.','master': '*'}
    xfmt = {'fmt':'c','alpha': 0.2}
    # Plot Dendrogram
    if plotting_grouped:
        # Create/update tmp_id_no
        ctr.clusters = ctr.clusters.assign(tmp_id_no = list(range(len(ctr))))
        R = ctr.dendrogram(ax=axes[0], xlabels=['tmp_id_no'], scalar=[None], title=f'{title}\n')
        # Plot threshold level
        axes[0].plot(axes[0].get_xlim(), [1 - cct]*2, 'r:')
        # Plot etype for each leaf
        xcoords = axes[0].get_xticks()
        for _e, _l in enumerate(R['leaves']):
            row = ctr._c.iloc[_l,:]
            axes[0].plot(xcoords[_e], 0.01, emap[row.etype], zorder=10)
        # Label Axis
        axes[0].set_xlabel('Event Index Number')
    else:
        R = None

    # Plot Scatters
    unit = {'time': '', 'depth': '[km]', 'longitude': '[$^o$E]', 'latitude': '[$^o$N]'}
    for _ax,_xf,_yf,_xs,_ys in [(1,'time','depth',False,1e-3),(2,'longitude','latitude',False,False),(3,'longitude','depth',False, 1e-3)]:
        # Overlay with everything in i_ctr
        for _etype in ctr._c.etype.unique():
            j_df_c = ctr.clusters[ctr._c.etype == _etype]
            xv = j_df_c[_xf].copy().values
            if _xs:
                xv *= _xs
            yv = j_df_c[_yf].copy().values
            if _ys:
                yv *= _ys
            # PLOT
            label = f'{_etype} (ct: {len(j_df_c)})'
 
            if plotting_grouped:
                fmt = emap[_etype]
                label = f'Grouped {label}'
                zorder=2
                alpha=1
            else:
                fmt = xfmt['fmt']+emap[_etype][1]
                label = f'Unaffiliated {label}'
                zorder=1
                alpha=0.2

            axes[_ax].plot(xv, yv, fmt, label=label, zorder=zorder, alpha=alpha)

        # if _yf == 'depth':
        #     ylims = axes[_ax].get_ylim()
        #     if ylims[0] < ylims[1]:
        #         axes[_ax].set_ylim([ylims[1],ylims[0]])
        
        # Label Axes
        axes[_ax].set_xlabel(_xf[0].upper() + _xf[1:] + unit[_xf])
        axes[_ax].set_ylabel(_yf[0].upper() + _yf[1:] + unit[_yf])

    # If intending to plot grouped events
    if plotting_grouped:

        ## TRACE SUBPLOT ##

        # Get Shifts & Plot shifted traces
        tr_list = [tmp.st[0].copy() for tmp in ctr]
        # NEW: Use pick SNR to pick master, rather than overall trace RMS
        snr_list = ctr._c.snr
        _emaster = np.argmax(snr_list)
        master_tr = tr_list[_emaster]
        if calculate_shifts:
            shifts_corrs = align_traces(tr_list, shift_len=50, master=master_tr, positive=positive_shift)
            display_cc = True
        else:
            shifts_corrs = np.full(shape=(2,len(tr_list)), fill_value=1)
            display_cc = False
        _act= 0
        _mct = 0
        axes[4].plot([1,1], [-1, len(tr_list)+1], 'k--', alpha=0.5)
        for _e, tr in enumerate(tr_list):
            _row = ctr._c.iloc[_e,:]
            _tr = tr.copy()
            _tr.detrend('linear')
            _tr.normalize()
            _tr.normalize(norm=_tr.std()*3)
            _ccsign = np.sign(shifts_corrs[1][_e])
            _xv = _tr.times() + shifts_corrs[0][_e] - ctr[_e].prepick
            _yv = _ccsign*_tr.data/3. + _e
            # Plot trace
            axes[4].plot(_xv, _yv, emap[_row.etype][0], lw=0.5, zorder=2)
            # Mark master trace
            if _e == _emaster:
                axes[4].plot(_xv[0], _yv[0], emap['master'], color='goldenrod', zorder=4)
                axes[4].plot(_xv, _yv, color='goldenrod', lw=2, zorder=1)
            # Mark catalog picks
            if _row.eval_mode == 'manual':
                if _mct == 0:
                    label='Catalog Picks'
                else:
                    label=None
                axes[4].plot(_xv[0], _yv[0], 'md', zorder=3, label=label)
                _mct += 1
            # Mark modeled picks
            else:
                if _act == 0:
                    label='Modeled Picks'
                else:
                    label=None
                axes[4].plot(_xv[0], _yv[0], 'cs', zorder=3, label=label)
                _act += 1

        # If displaying correlation coefficients on traces
        if display_cc:
            if positive_shift:
                vkwargs = {'vmin': 0, 'vmax': 1}
            else:
                vkwargs = {'vmin': 0, 'vmax': 1}
            ch = axes[4].scatter(-1*np.ones(len(tr_list)), np.arange(len(tr_list)), c=np.abs(shifts_corrs[1]),
                                 s=16, cmap='inferno', zorder=3, **vkwargs)
            # plt.colorbar(ch, location='bottom', orientation='horizontal', fraction=0.05)

        # Annotate SNR
        ch = axes[4].scatter(-0.5*np.ones(len(ctr)),
                        ctr._c.tmp_id_no, s=16,
                        c=ctr._c.snr, cmap='cool')
        plt.colorbar(ch, ax=axes[5], location='top', orientation='horizontal')
        axes[4].yaxis.tick_right()
        axes[4].set_ylabel('-- Towards Present Day -->')
        # axes[4].set_yticks(np.arange(len(ctr))*2.,ctr._c.index)
        axes[4].tick_params('y',labelsize=6)
        if calculate_shifts:
            axes[4].set_xlabel('Time Relative to\nXCorr Adjusted Arrival [sec]')
        else:
            axes[4].set_xlabel('Time Relative to\nCatalog Arrivals [sec]')
        


        ## CORRELATION COEFFICIENT SUBPLOT ###
        XX, YY = np.meshgrid(np.arange(0,len(ctr)), np.arange(0,len(ctr)))
        # XX, YY = np.meshgrid(ctr._c.id_no, ctr._c.id_no)
        ch = axes[5].pcolor(XX, YY, 1 - ctr.dist_mat, cmap='inferno', vmin=0, vmax=1)
        plt.colorbar(ch, ax=axes[5], orientation='horizontal', location='top')
        # Plot etypes
        for etype in ctr._c.etype.unique():
            rows = ctr._c[ctr._c.etype==etype]
            axes[5].plot(rows.tmp_id_no, rows.tmp_id_no, emap[etype], ms=4, zorder=3)
        axes[5].set_xlabel('Event Index Number')
        
        ## OVERLAY MASTER EVENT ON ALL OTHER PLOTS
        # Get the dendrogram x-coordinate for the master event
        mdxc = axes[0].get_xticks()[np.argwhere(R['leaves']==_emaster)]
        # Get the row entry corresponding to the master event
        mrow = ctr._c.iloc[_emaster, :]
        # Dendrogram
        axes[0].plot(mdxc, 0.01, emap['master'], color='goldenrod', zorder=20, ms=4)
        # Depth/time scatter
        axes[1].plot(mrow.time, mrow.depth*1e-3, emap['master'], color='goldenrod', zorder=20, ms=4, label='Reference Event')
        # Lat Lon scatter
        axes[2].plot(mrow.longitude, mrow.latitude, emap['master'], color='goldenrod', zorder=20, ms=4)
        # Lon Depth scatter
        axes[3].plot(mrow.longitude, mrow.depth*1e-3, emap['master'], color='goldenrod', zorder=20, ms=4)
        # Correlation Matrix Scatter
        axes[5].plot(_emaster, _emaster, emap['master'], color='goldenrod', ms=2, zorder=4)
        
    # Finally, if rendering legend, do so on figure 1
    if render_legend:
        axes[1].legend(ncol=2, loc='lower left', fontsize=8)
        axes[4].legend(ncol=2, loc='lower left', fontsize=8)
    return axes, R




#### LOADING ###

def load_data(file):
    print('loading')
    ctr = ClusteringTribe().read(str(file))
    return ctr

#### RE-THRESHOLDING ####
def rethreshold(ctr, cct=cct):
    ctr.cct_regroup(cct, inplace=True)
    ctr.reindex_columns('xcc');
    ctr._c.time = [pd.Timestamp(row.time) for _, row in ctr._c.iterrows()]
    ctr.cct_regroup(corr_thresh=cct,inplace=True)
    ctr.reindex_columns()
    return ctr


#### ITERATIVE PLOT GENERATION ###
def synopsis(ctr):
    # Get net/station set for this clustering tribe
    nslc_set = set()
    for tmp in ctr:
        for tr in tmp.st:
            nslc_set.add(tr.id)

    nslcstr = ''
    for _nslc in nslc_set:
        nslcstr += f'{_nslc}, '
    nslcstr = nslcstr[:-2]

    # Start by plotting everything
    axes, R = plot_summary(ctr, title=f'Stations: {nslcstr}',
                        axes=None, plotting_grouped=True, calculate_shifts=False)
    # Flip depth axes to have positive down
    for _x in [1, 3]:
        ylim1 = axes[_x].get_ylim()
        axes[_x].set_ylim(ylim1[1], ylim1[0])

def grouped_plots(ctr, min_members=min_members):
    grouped = []
    orphans = []

    # Get net/station set for this clustering tribe
    nslc_set = set()
    for tmp in ctr:
        for tr in tmp.st:
            nslc_set.add(tr.id)

    nslcstr = ''
    for _nslc in nslc_set:
        nslcstr += f'{_nslc}, '
    nslcstr = nslcstr[:-2]     

    # # Apply group number labels to the dendrogram
    # xcoords = axes[0].get_xticks()
    # id_nos = ctr._c.tmp_id_no
    # for _gn, _ct in ctr._c.correlation_cluster.value_counts().items():
    #     if _ct >= min_members:
    #         _xt = [xcoords[_x] for _x in ctr._c[ctr._c.correlation_cluster == _gn]]

    # Get alias to clusters dataframe
    df_c = ctr._c.copy()
    # Get group membership counts
    g_counts = df_c.xcc.value_counts()

    # Iterate across groups
    for _gn, _ct in g_counts.items():
        # If insufficient membership, log EVID(s) to orphans
        if _ct < min_members:
            orphans += list(df_c[df_c.xcc==_gn].index.values)
        # If sufficient membership, log EVIDs to grouped
        else:
            # Get this group subset and add to grouped EVIDs list
            grouped_subset = df_c[df_c.xcc==_gn].index.values
            grouped += list(grouped_subset)
            # Get the clusering tribe of included events
            i_ctr = ctr.get_subset(grouped_subset)
            # Get the clustering tribe of excluded events
            x_ctr = ctr.get_subset(df_c[df_c.xcc !=_gn].index.values)
            # Plot grouped events
            axes, Ri = plot_summary(i_ctr, axes=None, title=f'Stations: {nslcstr} Group: {_gn:.0f}', render_legend=False)
            # Plot unaffiliated events
            axes, Rx = plot_summary(x_ctr, axes=axes, plotting_grouped=False, render_legend=True)
            for _x in [1, 3]:
                ylim1 = axes[_x].get_ylim()
                axes[_x].set_ylim(ylim1[1], ylim1[0])

    # Iterate across etypes for orphans
    o_ctr = ctr.get_subset(orphans)
    for _etype in o_ctr._c.etype.unique():
        e_orphans = o_ctr._c[o_ctr._c.etype == _etype].index.values
        if len(e_orphans) < 2:
            continue
        oe_ctr = ctr.get_subset(e_orphans)
        xe_ctr = ctr.get_subset(list(e_orphans) + grouped)
        axes, Roe = plot_summary(oe_ctr, axes=None, title=f'Stations: {nslcstr} Orphaned: {_etype}', render_legend=False)
        axes, Rxe = plot_summary(xe_ctr, axes=axes, plotting_grouped=False, render_legend=True)
        for _x in [1, 3]:
            ylim1 = axes[_x].get_ylim()
            axes[_x].set_ylim(ylim1[1], ylim1[0])



### RUN AS MAIN ###
if __name__ == '__main__':
    plt.close('all')
    ctr = load_data(file)
    ctr = rethreshold(ctr, cct=cct)
    synopsis(ctr)
    grouped_plots(ctr, min_members=min_members)

    plt.show()

#