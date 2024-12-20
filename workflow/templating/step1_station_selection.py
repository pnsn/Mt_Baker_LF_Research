import os, logging, glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
from obsplus import EventBank
from obspy.geodetics import locations2degrees
from obspy import read_inventory

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

PSN = ['MBW','MBW2','SHUK','RPW','RPW2','JCW','CMW','SAXON','MULN']

# Load inventory
for _e, _f in enumerate(glob.glob(str(INVD/'*.xml'))):
    if _e == 0:
        INV = read_inventory(_f)
    else:
        INV += read_inventory(_f)

# Load event summary
df_eb = EBANK.read_index()
df_eb = df_eb.sort_values(by='time')
Logger.info(f'EventBank has a total of {len(df_eb)} events')
# Calculate distances from Mt. Baker Summit
df_eb = df_eb.assign(mbs_dist_km=[
    111.2* locations2degrees(LATMBS, LONMBS, row.latitude, row.longitude)
    for _, row in df_eb.iterrows()])

# Construct P-pick mapping for preferred_origins to each station
picked_stations = {}
picking_continuity = {}

_e = -1
if Logger.level < 20:
    _dtqdmf = True
else:
    _dtqdmf = False

# Subset Events by Distance
newcount = len(df_eb[df_eb.mbs_dist_km <= RAD_LIM_KM])
Logger.info(f'Subsetting to {RAD_LIM_KM} km distance retains {newcount} events')

# Get pick counts for all events' preferred origins
picked_times = []
picked_evids = []
for event_id in tqdm(df_eb[df_eb.mbs_dist_km <= RAD_LIM_KM].event_id, disable=_dtqdmf):
    _e += 1
    Logger.debug(f"{event_id} ({_e+1} of {newcount})")
    cat = EBANK.get_events(event_id=event_id)
    prefor = cat[0].preferred_origin()
    # Skip no-pick prefor
    if len(prefor.arrivals) == 0:
        _e -= 1
        continue
    for arr in prefor.arrivals:
        if arr.phase == 'P':
            pick = arr.pick_id.get_referred_object()
            chan = pick.waveform_id.id
            if chan not in picked_stations.keys():
                picked_stations.update({chan: [event_id]})
                # pad with 0's on initialization
                picking_continuity.update({chan: [0]*(_e - 1) + [1]})
            else:
                picked_stations[chan].append(event_id)
                picking_continuity[chan].append(1)
    # If channel wasn't picked, but has been in the past, assign -1 to penalize
    for _k, _v in picking_continuity.items():
        if len(_v) != _e + 1:
            picking_continuity[_k].append(-1)

    picked_times.append(prefor.time)
    picked_evids.append(event_id)


# Construct pick continuity dataframe
df_cont = pd.DataFrame(picking_continuity).T
df_cont = df_cont.astype(pd.SparseDtype('int', 0))
# Get number of picks minus number of missed picks
npossible = (df_cont**2).sum(axis=1)
# Get number of picks
npicks = df_cont[df_cont > 0].sum(axis=1).values
# Number of picks minus number of missed picks
sums = df_cont.sum(axis=1).values

df_cont = df_cont.assign(summed=sums)
df_cont = df_cont.assign(npicks=npicks)
df_cont = df_cont.sort_values(['npicks','summed'], ascending=False)
df_cont.pop('summed');
df_cont.pop('npicks')
df_cont = df_cont.T
df_cont.index = picked_evids

# Extract unique station codes
auto_pref_sta = set()
for col in df_cont.iloc[:,:20].columns:
    n, s, l, c = col.split('.')
    auto_pref_sta.add(s)

# Construct pick count series
ser_picks = pd.Series({_k: len(_v) for _k, _v in picked_stations.items()}).sort_values(ascending=False)
event_id_set = set(df_eb.event_id.values)
for _cha, _val in ser_picks.items():
    Logger.info(f'{_cha} - {_val}')
    if _val < 50:
        break


# Get the set of events observed by preferred stations
pref_sta_evid_set = set()
pref_sta_set = set()
for _k, _v in picked_stations.items():
    for _psn in PSN:
        if _psn in _k:
            pref_sta_evid_set = pref_sta_evid_set.union(set(_v))
            pref_sta_set.add(_k)
Logger.info(f'preferred station channels and pick counts:')
pref_sta_list = list(pref_sta_set)
pref_sta_list.sort()
for _pst in pref_sta_list:
    Logger.info(f'\t{_pst} - {ser_picks[_pst]}')
# Get the 
missed_event_ids = event_id_set.difference(pref_sta_evid_set)
fmissed = len(missed_event_ids)/newcount
if fmissed > 0.1:
    Logger.warning(f'preferred stations list results in more than 10% event loss')
else:
    Logger.info(f'preferred stations list results in a {fmissed*100:.3f}% event loss ({len(missed_event_ids)})')

# Write preferred station evid set to disk
with open(str(PSPEF), 'w') as LOG:
    LOG.write('event_id\n')
    for evid in pref_sta_evid_set:
        LOG.write(f'{evid}\n')

# Write continuity dataframe to disk
# Subset to preferred stations
df_cont = df_cont.filter(regex='|'.join(PSN))
# Sort by station name
df_cont = df_cont.T.sort_index().T
df_cont.to_csv(str(OUTD/'preferred_stachan_pick_continuity.csv'),
               header=True, index=True)

# Plot station and event distribution, showing other suggested stations
plt.figure()
df_ebf = df_eb[df_eb.mbs_dist_km <= RAD_LIM_KM]
plt.plot(df_eb.longitude, df_eb.latitude, 'k.', ms=2, zorder=1)

plt.scatter(df_ebf.longitude, df_ebf.latitude, c=df_ebf.mbs_dist_km, s=4,zorder=2)
plt.colorbar()
for sta in auto_pref_sta:
    inv = INV.select(station=sta)
    for net in inv.networks:
        for sta in net.stations:
            if sta.code in PSN:
                marker = 'rv'
            else:
                marker = 'kv'
            plt.plot(sta.longitude, sta.latitude, marker)
            if '2' in sta.code:
                adj = -0.01
                ha = 'right'
            else:
                adj = 0.01
                ha ='left'
            
            plt.text(sta.longitude+adj, sta.latitude+ adj, f'{net.code}.{sta.code}',
                     color='r', ha=ha, fontweight='extra bold')
            break
plt.axis('square')   

# Plot Pick Continuity for preferred stations

fig = plt.figure()
ax = fig.add_subplot(111)
ax.pcolor(df_cont.values.T)
ax.set_yticks(np.arange(len(df_cont.columns))+0.5, df_cont.columns.values)
ax.set_xlabel('Event Index [no.] Increasing in Time')
ylims = ax.get_ylim()
for _e, ind in enumerate(df_cont.index):
    _s = df_cont.loc[ind]
    if all(_s <= 0):
        plt.plot([_e]*2, [0, len(df_cont.columns)], 'r-')
ax.set_ylim(ylims)

# for _e, (_k, _v) in enumerate(picking_continuity.items()):
#     plt.plot(_v, [_e]*len(_v), 'k-')
    
# Plot Pick continuity