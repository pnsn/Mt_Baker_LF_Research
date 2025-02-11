
import os, logging, glob
from pathlib import Path

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
CAT_MEMBERSHIP = PDDIR / "catalog" / "P1S1_Event_ID_Catalog_Membership.csv"


# STATION SELECTION PARAMETERS
# Network Code(s) to include
NETS = ['UW']      
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

def main():
    # Load catalog membership CSV
    df_cat = pd.read_csv(CAT_MEMBERSHIP, index_col=[0], parse_dates=['time'])
    # Connect to IRIS webservice client
    IRIS = Client('IRIS')
    # Connect to Event Bank
    EBANK = EventBank(EBBP)

    INV = IRIS.get_stations(latitude=LAT_REF,
                            longitude=LON_REF,
                            maxradius=RAD_LIM_KM/111.2,
                            network=','.join(NETS),
                            channel=','.join(CHANS))
    breakpoint()

    # Read Index
    df_eb = EBANK.read_index()
    for evid in df_eb.event_id:
        # If not included in CAT0, continue to next
        if not df_cat.loc[evid,'CAT0']:
            continue

    

    # Load inventory, subsetting for desired channels
    INV = Inventory()
    for _e, _f in enumerate(glob.glob(str(INVD/'*.xml'))):
        inv = read_inventory(_f)
        for net in NETS:
            INV += inv.select(network=net, channel=CHANS)




    #     # Get active channels for this event
    #     inv = INV.select(time=prefor.time)
    #     # Check if it has any associated picks
    #     if len(prefor.arrivals) == 0:
    #         Logger.warning(f'No associated picks for {event_id} - skipping to next')
    #         continue
    #     for arr in prefor.arrivals:
    #         pick = arr.pick_id.get_referred_object()
    #         nslc = pick.waveform_id.id
    #         line = [event_id, arr.phase, nslc] + nslc.split('.')
    #         pick_info.append(line)

    
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


if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
    Logger.addHandler(CriticalExitHandler(exit_code=1))

    # RUN MAIN 
    main()

    if SHOW:
        plt.show()