import os, logging, glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
from obsplus import EventBank
from obspy.geodetics import locations2degrees
from obspy import read_inventory, Inventory, UTCDateTime

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to station response files
INVD = ROOT / 'data' / 'XML' / 'RESP'
# Output Path
OUTD = ROOT / 'processed_data' / 'workflow' / 'catalog'
# Preferred Station Picked Events File (save file)
PSPEF = OUTD / 'preferred_station_picked_event_ids.csv'


EBANK = EventBank(EBBP)

LATMBS, LONMBS, HBMS = 48.7745,-121.8172, 3286.
RAD_LIM_KM = 30.
STARTDATE = pd.Timestamp('2001-01-01T00:00:00')
MINPICK=50

NETS = ['UW']
# Manually identified stations to use
STAS = ['MBW','MBW2','SHUK','RPW','RPW2','JCW','CMW','SAXON','MULN','PASS']
CHANS = '[BHE][HN][ZNE123]'
# Load inventory, subsetting for desired channels
INV = Inventory()
for _e, _f in enumerate(glob.glob(str(INVD/'*.xml'))):
    inv = read_inventory(_f)
    for net in NETS:
        INV += inv.select(network=net, channel=CHANS)


# Load event summary
df_eb = EBANK.read_index()
df_eb = df_eb.sort_values(by='time')
Logger.info(f'EventBank has a total of {len(df_eb)} events')
# Calculate distances from Mt. Baker Summit
df_eb = df_eb.assign(mbs_dist_km=[
    111.2* locations2degrees(LATMBS, LONMBS, row.latitude, row.longitude)
    for _, row in df_eb.iterrows()])

### ANNOUNCE SUBSETTING IMPACTS ON EVENT COUNTS ###
# Test Subset Events by Date
newcount = len(df_eb[df_eb.time >= STARTDATE])
Logger.info(f'Subsetting to {STARTDATE} to PRESENT retains {newcount} events')

# Test Subset Events by Distance
newcount = len(df_eb[df_eb.mbs_dist_km <= RAD_LIM_KM])
Logger.info(f'Subsetting to {RAD_LIM_KM} km distance retains {newcount} events')

# Apply Subsetting
df_eb_sub = df_eb[(df_eb.time >= STARTDATE) & (df_eb.mbs_dist_km <= RAD_LIM_KM)]
newcount = len(df_eb_sub)
Logger.info(f'Subsetting by both retains {newcount} events')

### CONDUCT CONTINUITY/COVERAGE CHECK
if Logger.level < 20:
    _dtqdmf = True
else:
    _dtqdmf = False

picked_channels = set()
pick_info = []
_e = -1
for event_id in tqdm(df_eb_sub.event_id, disable = _dtqdmf):
    # Increment up indexer
    _e += 1
    Logger.debug(f"{event_id} ({_e+1} of {newcount})")
    # Fetch event metadata
    cat = EBANK.get_events(event_id=event_id)
    # Get preferred origin metadata
    prefor = cat[0].preferred_origin()
    # Get active channels for this event
    inv = INV.select(time=prefor.time)
    # Check if it has any associated picks
    if len(prefor.arrivals) == 0:
        Logger.warning(f'No associated picks for {event_id} - skipping to next')
        continue
    for arr in prefor.arrivals:
        pick = arr.pick_id.get_referred_object()
        nslc = pick.waveform_id.id
        line = [event_id, arr.phase, nslc] + nslc.split('.')
        pick_info.append(line)

# Create Picked Dataframe    
df_picked = pd.DataFrame(pick_info, columns=['event_id','phase','nslc','net','sta','loc','chan'])
# Subset to only channels in inventory
df_picked = df_picked[df_picked.nslc.isin(INV.get_contents()['channels'])]
# Subset to only P picked stations
ser_top_p = df_picked[df_picked.phase=='P'].sta.value_counts()
# Subset to stations with GE MINPICK P picks
top_p = ser_top_p[ser_top_p >= MINPICK].index.values

# Create Availability array
df_p_cont = df_picked[(df_picked.phase=='P') & (df_picked.sta.isin(top_p))][['event_id','sta']]
df_p_cont = df_p_cont.pivot_table(index=['event_id'], columns=['sta'], aggfunc=lambda x: 1, fill_value=0)
# Translate and sort channels by number of picks
df_p_cont = df_p_cont.T.loc[top_p]

# Render Figure for  Pick Availability & Data Informed Station Selection
fig = plt.figure()
gs = fig.add_gridspec(ncols=1, nrows=2, wspace=0, hspace=0)
axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
axes[0].pcolor(df_p_cont)
axes[0].set_yticks(np.arange(len(df_p_cont.index))+0.5, df_p_cont.index.values)
axes[0].set_title(f'Pick Availability for Stations with $\geq$ {MINPICK} Picks')
for _e in range(3, len(top_p) + 1):
    nobs = df_p_cont.iloc[:_e,:].sum(axis=0)
    nobs_bins = nobs.value_counts()
    if 0 in nobs_bins.index:
        n0 = nobs.value_counts()[0]
    else:
        n0 = 0
    axes[1].fill_between(list(range(len(nobs))), [0]*len(nobs), nobs.values, label=f'Top {_e} ({n0} no-obs)',zorder=20-_e)
    if n0 > 0:
        for _j, _n in enumerate(nobs.values):
            if _n == 0:
                axes[0].plot([_j, _j],[0,_e], 'r-', alpha=0.2)
    if n0 == 0:
        axes[0].fill_between([0, len(nobs)],[_e]*2, [len(df_p_cont)]*2,color='k',alpha=0.5)
        axes[0].text(len(nobs)//2, np.mean([_e, len(df_p_cont)]), 'Redundant', color='w')
        break
axes[1].set_xlabel('Event Index (--Towards Present-->)')
axes[0].set_ylabel('Station Code (<--Increasing Total Picks--)')
axes[1].set_ylabel('Number of Analyst Picks [ct.]')
axes[1].legend(ncols=1)
axes[1].set_xlim(axes[0].get_xlim())

# Save df_picked for non-redundant stations
df_picked[df_picked.sta.isin(top_p[:_e])].to_csv(str(OUTD/'preferred_event_sta_picks.csv'), header=True, index=False)

plt.show()

# ax.set_xticks(np.arange(len(df_p_cont.index))+0.5, df_p_cont.index.values)
breakpoint()
# # Convert to boolean array of present/not present
# df['multi_index']


# # Construct P-pick mapping for preferred_origins to each station
# picked_stations = {}
# picking_continuity = {}

# _e = -1

# aa

# # Get pick counts for all events' preferred origins
# picked_times = []
# picked_evids = []
# # Iterate across events in chronologic order
# for event_id in tqdm(df_eb_sub.event_id, disable=_dtqdmf):
#     # Increment up indexer
#     _e += 1
#     Logger.debug(f"{event_id} ({_e+1} of {newcount})")
#     # Fetch event record
#     cat = EBANK.get_events(event_id=event_id)
#     # Get preferred origin object
#     prefor = cat[0].preferred_origin()
#     # Skip prefor with no associated picks
#     if len(prefor.arrivals) == 0:
#         _e -= 1
#         continue
#     # Iterate across prefor arrivals
#     for arr in prefor.arrivals:
#         # If the arrival is a P-wave
#         if arr.phase == 'P':
#             # Get the associated pick
#             pick = arr.pick_id.get_referred_object()
#             # Get teh channel code for that pick
#             chan = pick.waveform_id.id

#             ## CONTINUITY LOGGING
#             # If channel ID is new
#             if chan not in picked_stations.keys():
#                 # Populate the picked station & include the event_id
#                 picked_stations.update({chan: [event_id]})
#                 # pad with 0's on initialization plus 1 for shoing a pick was present
#                 picking_continuity.update({chan: [0]*(_e - 1) + [1]})
#             # If the picked channel is active at this point
#             else:
#                 picked_stations[chan].append(event_id)
#                 picking_continuity[chan].append(1)

#     # Iterate across all channels present in picking_continuity
#     for _k, _v in picking_continuity.items():
#         # if the channel was active during this event
#         if chan in INV.select(time=prefor.time).get_contents()['channels']:    
#             # If the continuity vector is missing an element
#             if len(_v) < _e + 1:
#                 # append a missed pick value
#                 picking_continuity[_k].append(-1)
#             # DEBUG: Safety catch of too-long
#             elif len(_v) > _e + 1:
#                 breakpoint()
#         # If the channel is not active
#         else:
#             # And it shows a missing last pick
#             if len(_v) < _e + 1:
#                 # Append a non-active pick value
#                 picking_continuity[_k].append(0)

#         #     if len(_v) != _e + 1:
#         #         picking_continuity[_k].append(-1)
#         # # If the channel is inactive
#         # else:
#         #     if len(_v) != _e + 1:
#         #         picking_continuity[_k].append(0)
#     picked_times.append(prefor.time)
#     picked_evids.append(event_id)


# # Construct pick continuity dataframe
# df_cont = pd.DataFrame(picking_continuity).T
# df_cont = df_cont.astype(pd.SparseDtype('int', 0))
# # Get number of picks minus number of missed picks
# npossible = (df_cont**2).sum(axis=1)
# # Get number of picks
# npicks = df_cont[df_cont > 0].sum(axis=1).values
# # Number of picks minus number of missed picks
# sums = df_cont.sum(axis=1).values

# df_cont = df_cont.assign(summed=sums)
# df_cont = df_cont.assign(npicks=npicks)
# df_cont = df_cont.sort_values(['npicks','summed'], ascending=False)
# df_cont.pop('summed');
# df_cont.pop('npicks')
# df_cont = df_cont.T
# df_cont.index = picked_evids

# # Extract unique station codes
# auto_pref_sta = set()
# for col in df_cont.iloc[:,:20].columns:
#     n, s, l, c = col.split('.')
#     auto_pref_sta.add(s)

# # Construct pick count series
# ser_picks = pd.Series({_k: len(_v) for _k, _v in picked_stations.items()}).sort_values(ascending=False)
# event_id_set = set(df_eb.event_id.values)
# for _cha, _val in ser_picks.items():
#     Logger.info(f'{_cha} - {_val}')
#     if _val < 50:
#         break


# # Get the set of events observed by preferred stations
# pref_sta_evid_set = set()
# pref_sta_set = set()
# for _k, _v in picked_stations.items():
#     for _psn in STAS:
#         if _psn in _k:
#             pref_sta_evid_set = pref_sta_evid_set.union(set(_v))
#             pref_sta_set.add(_k)

# # REPORT TO LOG
# Logger.info(f'preferred station channels and pick counts:')
# pref_sta_list = list(pref_sta_set)
# pref_sta_list.sort()
# for _pst in pref_sta_list:
#     Logger.info(f'\t{_pst} - {ser_picks[_pst]}')
# # Get the 
# missed_event_ids = event_id_set.difference(pref_sta_evid_set)
# fmissed = len(missed_event_ids)/newcount
# if fmissed > 0.1:
#     Logger.warning(f'preferred stations list results in more than 10% event loss')
# else:
#     Logger.info(f'preferred stations list results in a {fmissed*100:.3f}% event loss ({len(missed_event_ids)})')

# # Write preferred station evid set to disk
# with open(str(PSPEF), 'w') as LOG:
#     LOG.write('event_id\n')
#     for evid in pref_sta_evid_set:
#         LOG.write(f'{evid}\n')

# # Write continuity dataframe to disk
# # Subset to preferred stations
# df_cont = df_cont.filter(regex='|'.join(STAS))
# # Sort by station name
# df_cont = df_cont.T.sort_index().T
# df_cont.to_csv(str(OUTD/'preferred_stachan_pick_continuity.csv'),
#                header=True, index=True)

# # Plot station and event distribution, showing other suggested stations
# plt.figure()
# df_ebf = df_eb[df_eb.mbs_dist_km <= RAD_LIM_KM]
# plt.plot(df_eb.longitude, df_eb.latitude, 'k.', ms=2, zorder=1)

# plt.scatter(df_ebf.longitude, df_ebf.latitude, c=df_ebf.mbs_dist_km, s=4,zorder=2)
# plt.colorbar()
# for sta in auto_pref_sta:
#     inv = INV.select(station=sta)
#     for net in inv.networks:
#         for sta in net.stations:
#             if sta.code in STAS:
#                 marker = 'rv'
#             else:
#                 marker = 'kv'
#             plt.plot(sta.longitude, sta.latitude, marker)
#             if '2' in sta.code:
#                 adj = -0.01
#                 ha = 'right'
#             else:
#                 adj = 0.01
#                 ha ='left'
            
#             plt.text(sta.longitude+adj, sta.latitude+ adj, f'{net.code}.{sta.code}',
#                      color='r', ha=ha, fontweight='extra bold')
#             break
# plt.axis('square')   

# # Plot Pick Continuity for preferred stations

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.pcolor(df_cont.values.T)
# ax.set_yticks(np.arange(len(df_cont.columns))+0.5, df_cont.columns.values)
# ax.set_xlabel('Event Index [no.] Increasing in Time')
# ylims = ax.get_ylim()
# for _e, ind in enumerate(df_cont.index):
#     _s = df_cont.loc[ind]
#     if all(_s <= 0):
#         plt.plot([_e]*2, [0, len(df_cont.columns)], 'r-')
# ax.set_ylim(ylims)

# # for _e, (_k, _v) in enumerate(picking_continuity.items()):
# #     plt.plot(_v, [_e]*len(_v), 'k-')
    
# # Plot Pick continuity

# plt.show()