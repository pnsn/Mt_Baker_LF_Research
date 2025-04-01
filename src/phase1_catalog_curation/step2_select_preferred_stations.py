
import os, logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from obspy.clients.fdsn.client import Client
from obsplus import EventBank

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Absolute path to response inventory directory
INVD = ROOT / "data" / "XML" / "RESP"
# Aboslute path to processed data / catalog directory
PDDIR = ROOT / "processed_data" 
# Catalog Membership CSV File (From phase1 step1)
CAT_MEMBERSHIP = PDDIR / "catalog" / "P1S1_Catalog_Profile.csv"
# 
SAVEPATH = PDDIR / "catalog" 
# Save file for station-event-pick lines
SAVESEP = SAVEPATH / "P1S2_Preferred_Sta_Event_Picks.csv"
# Save file for preferred station codes
SAVESTAPREF = SAVEPATH / 'P1S2_Preferred_Stations.csv'
# STATION SELECTION PARAMETERS
# Network Code(s) to include
NETS = ['UW','CN','TA','PB','CC','UO']
# Channel Codes to include
CHANS = []
for _b in 'BHE':
    for _i in 'HN':
        for _c in 'ZNE12':
            CHANS.append(''.join([_b,_i,_c]))
# Geographic limits
LAT_REF = 48.7745   # [deg N] Reference point latitude
LON_REF= -121.8172  # [deg E] Reference point longitude
RAD_LIM_KM = 100.   # [km] Radius limit for station query
INV_LEVEL = 'channel' # Inventory query depth
MINSTAPICK = 50
MINEVEPICK = 1
def main():
    # Load catalog membership CSV
    df_cat = pd.read_csv(CAT_MEMBERSHIP, index_col=[0], parse_dates=['prefor_time'])
    # Connect to IRIS webservice client
    IRIS = Client('IRIS')
    # Connect to Event Bank
    EBANK = EventBank(EBBP)
    # Load inventory from webservices
    INV = IRIS.get_stations(latitude=LAT_REF,
                            longitude=LON_REF,
                            maxradius=RAD_LIM_KM/111.2,
                            network=','.join(NETS),
                            channel=','.join(CHANS),
                            level=INV_LEVEL)

    # Subset by CAT1 - Within 30 km of MB & after 1980-01-01T00:00:00
    df_cat1 = df_cat[df_cat.CAT1]
    # Iterate across CAT1 events
    pick_status = defaultdict(list)
    for _e, (evid, row) in enumerate(df_cat1.iterrows()):
        # Log progress
        Logger.info(f'Processing {evid} ({_e+1}/{len(df_cat1)})')
        # Get prefor from EBANK
        cat = EBANK.get_events(event_id=evid)
        prefor = cat[0].preferred_origin()
        # Subset inventory by origin time (active stations during event)
        inv = INV.select(time=prefor.time)
        # Warn if no picks are available
        if len(prefor.arrivals) == 0:
            Logger.warning(f'No associated picks for {evid}')
            continue
        # Iterate across arrivals
        for arr in prefor.arrivals:
            # Get associated pick
            pick = arr.pick_id.get_referred_object()
            # Get NSLC channel code
            nslc = pick.waveform_id.id
            # Split into parts
            parts = nslc.split('.')
            # Make sure network code is acknowledged
            if parts[-1] in CHANS:
                # Append to defaultdict accumulator
                pick_status['event_id'].append(evid)
                pick_status['phase'].append(arr.phase)
                pick_status['nslc'].append(nslc)
                pick_status['time'].append(pick.time)
                pick_status['resource_id'].append(pick.resource_id.id)
                pick_status['etype'].append(row.etype)
                
                for _x, _y in zip(['network','station','location','channel'],parts):
                    pick_status[_x].append(_y)
            else:
                continue
    # Form selected picks dataframe
    df_picks = pd.DataFrame(pick_status)
    
    # Conduct Pick Continuity/Availability Analysis
    # Get top P-picked stations
    sta_P_pick_counts = df_picks[df_picks.phase=='P'].station.value_counts()
    sta_S_pick_counts = df_picks[df_picks.phase=='S'].station.value_counts()
    sta_E_pick_counts = df_picks.station.value_counts()
    # Get top P-picked stations
    top_P_picked_stas = sta_P_pick_counts[sta_P_pick_counts >= MINSTAPICK].index.values
    top_S_picked_stas = sta_S_pick_counts[sta_S_pick_counts >= MINSTAPICK].index.values
    top_picked_stas = sta_E_pick_counts[sta_E_pick_counts >= MINSTAPICK].index.values
    
    # Pivot table for just P picks
    ptable = df_picks[(df_picks.phase=='P') & (df_picks.station.isin(top_P_picked_stas))][['event_id','station']]
    ptable = ptable.pivot_table(index=['event_id'], columns=['station'], aggfunc=lambda x: 1, fill_value=0)
    # Translate and sort channels by pick count
    ptable = ptable.T.loc[top_P_picked_stas]

    # Pivot table for just S picks
    stable = df_picks[(df_picks.phase=='S') & (df_picks.station.isin(top_S_picked_stas))][['event_id','station']]
    stable = stable.pivot_table(index=['event_id'], columns=['station'], aggfunc=lambda x: 1, fill_value=0)
    # Translate and sort channels by pick count
    stable = stable.T.loc[top_S_picked_stas]
     
    # Pivot table for all picks
    etable = df_picks[df_picks.station.isin(top_picked_stas)][['event_id','station']]
    etable = etable.pivot_table(index=['event_id'], columns=['station'], aggfunc=lambda x: 1, fill_value=0)
    # Translate and sort channels by pick count
    etable = etable.T.loc[top_picked_stas]


    # Analyze for how many stations to include to cover all events with N picks
    include = set()
    event_pick_sums = ptable.cumsum(axis=0)
    for sta, row in event_pick_sums.iterrows():
        if any(_e < 1 for _e in row):
            include.add(sta)
            continue
        else:
            include.add(sta)
            break
    
    # Save preferred station list to disk
    df_out = df_picks[
        (df_picks.phase=='P')&\
        (df_picks.station.isin(include))&\
        (df_picks.event_id.isin(df_cat1.index))]
    df_out.to_csv(SAVESEP, header=True, index=False)
    return df_out


if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
    Logger.addHandler(CriticalExitHandler(exit_code=1))

    # RUN MAIN 
    df_out = main()
    ser_cts = df_out[['network','station']].value_counts()
    ser_cts.to_csv(SAVESTAPREF)



    # if SHOW:
    #     plt.show()

            # for event_id in tqdm(df_eb_sub.event_id, disable = _dtqdmf):
    #     # Increment up indexer
    #     _e += 1
    #     Logger.debug(f"{event_id} ({_e+1} of {newcount})")
    #     # Fetch event metadata
    #     cat = EBANK.get_events(event_id=event_id)
    #     # Get preferred origin metadata
    #     prefor = cat[0].preferred_origin()

    #     # Assess if this is a "well-constrained event"
    #     wcstatus = True


    # # Create Picked Dataframe    
    # df_picked = pd.DataFrame(pick_info, columns=['event_id','phase','nslc','net','sta','loc','chan'])
    # # Subset to only channels in inventory
    # df_picked = df_picked[df_picked.nslc.isin(INV.get_contents()['channels'])]
    # # Subset to only P picked stations
    # ser_top_p = df_picked[df_picked.phase=='P'].sta.value_counts()
    # # Subset to stations with GE MINPICK P picks
    # top_p = ser_top_p[ser_top_p >= MINPICK].index.values

    # # Create Availability array
    # df_p_cont = df_picked[(df_picked.phase=='P') & (df_picked.sta.isin(top_p))][['event_id','sta']]
    # df_p_cont = df_p_cont.pivot_table(index=['event_id'], columns=['sta'], aggfunc=lambda x: 1, fill_value=0)
    # # Translate and sort channels by number of picks
    # df_p_cont = df_p_cont.T.loc[top_p]


    # ### FIGURE RENDERING SECTION ###
    # well_located_evids = pd.DataFrame(well_located_evids, columns=['event_id','depth_type'])
    # # From TASK1 Figure - Effects of time, distance, and origin quality filtering
    # fig= plt.figure(figsize=(6,6))
    # gs = fig.add_gridspec(ncols=2, nrows=2)
    # axes = [fig.add_subplot(gs[_e]) for _e in range(4)]

    # for _lbl in ['Catalog','Radius','Radius+Time','R+T+Well-Constrained']:
    #     if _lbl == 'Catalog':
    #         df_e = df_eb
    #         fmt = 's'
    #         ms = 1
    #         color='black'
    #     elif _lbl == 'Radius':
    #         df_e = df_eb[df_eb.radius_offset_km <= RAD_LIM_KM]
    #         fmt = 'v'
    #         ms = 2
    #         color='firebrick'
    #     elif _lbl == 'Radius+Time':
    #         df_e = df_eb[(df_eb.radius_offset_km <= RAD_LIM_KM) & 
    #                     (df_eb.time >= STARTDATE)]
    #         fmt = 'x'
    #         ms = 3
    #         color='dodgerblue'
    #     elif _lbl == 'R+T+Well-Constrained':
    #         df_e = df_eb[df_eb.event_id.isin(well_located_evids.event_id)]
    #         fmt = '*'
    #         ms = 4
    #         color='goldenrod'

    #     plotfields = [('time','depth'),     ('longitude','latitude'),
    #                 ('time','magnitude'), ('time','vertical_uncertainty')]
                    
    #     for _b, (xfld, yfld) in enumerate(plotfields):
    #         if yfld in ['depth','vertical_uncertainty']:
    #             _yv = df_e[yfld].values*1e-3
    #             if yfld == 'depth':
    #                 _yv *= -1
    #         else:
    #             _yv = df_e[yfld].values
    #         _xv = df_e[xfld].values
    #         axes[_b].plot(_xv, _yv, fmt, ms=ms, color=color, label=f'{_lbl} ({len(_xv)})')
    #         axes[_b].set_xlabel(xfld)
    #         axes[_b].set_ylabel(yfld)
    #     axes[0].legend()

    # plt.savefig(str(OUTD / f'step1_event_selection_{DPI:d}dpi.{FMT.lower()}'),
    #             dpi=DPI,format=FMT, **plot_save_kwargs)

    # # Render Figure for Pick Availability & Data Informed Station Selection
    # fig = plt.figure(figsize=(8,6))
    # gs = fig.add_gridspec(ncols=1, nrows=2, wspace=0, hspace=0)
    # axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
    # axes[0].pcolor(df_p_cont)
    # axes[0].set_yticks(np.arange(len(df_p_cont.index))+0.5, df_p_cont.index.values)
    # axes[0].set_title(f'Pick Availability for Stations with $\geq$ {MINPICK} Picks')
    # for _e in range(3, len(top_p) + 1):
    #     nobs = df_p_cont.iloc[:_e,:].sum(axis=0)
    #     nobs_bins = nobs.value_counts()
    #     if 0 in nobs_bins.index:
    #         n0 = nobs.value_counts()[0]
    #     else:
    #         n0 = 0
    #     axes[1].fill_between(list(range(len(nobs))), [0]*len(nobs), nobs.values, label=f'Top {_e} ({n0} no-obs)',zorder=20-_e)
    #     if n0 > 0:
    #         for _j, _n in enumerate(nobs.values):
    #             if _n == 0:
    #                 axes[0].plot([_j, _j],[0,_e], 'r-', alpha=0.2)
    #     if n0 == 0:
    #         axes[0].fill_between([0, len(nobs)],[_e]*2, [len(df_p_cont)]*2,color='k',alpha=0.5)
    #         axes[0].text(len(nobs)//2, np.mean([_e, len(df_p_cont)]), 'Redundant', color='w')
    #         # If this is the first station that completes event converage, break iteration
    #         break
    # axes[1].set_xlabel('Event Index (--Towards Present-->)')
    # axes[0].set_ylabel('Station Code (<--Increasing Total Picks--)')
    # axes[1].set_ylabel('Number of Analyst Picks [ct.]')
    # axes[1].legend(ncols=1)
    # axes[1].set_xlim(axes[0].get_xlim())

    # plt.savefig(str(OUTD / f'step1_station_selection_{DPI:d}dpi.{FMT.lower()}'),
    #             dpi=DPI,format=FMT, **plot_save_kwargs)

    # # Save df_picked for non-redundant stations (or full set if there are still gaps)
    # df_out = df_picked[df_picked.sta.isin(top_p[:_e])]
    # # Annotate if this entry corresponds to a well-constrained event
    # df_out = df_out.assign(WCE=[x in well_located_evids.event_id.values for x in df_out.event_id])
    # # Save to disk
    # df_out.to_csv(str(OUTD/'preferred_event_sta_picks.csv'), header=True, index=False)


    # if display_after:
    #     plt.show()
