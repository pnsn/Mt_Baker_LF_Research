import os, logging
from pathlib import Path

import pandas as pd
from obsplus import EventBank

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"

EBANK = EventBank(EBBP)


PSN = ['MBW','RPW','JCW','SHUK','SAXON']

df_eb = EBANK.read_index()
# Construct P-pick mapping for preferred_origins to each station
picked_stations = {}
for _e, event_id in df_eb.event_id.items():
    Logger.info(f"{event_id} ({_e+1} of {len(df_eb)})")
    cat = EBANK.get_events(event_id=event_id)
    prefor = cat[0].preferred_origin()
    for arr in prefor.arrivals:
        if arr.phase == 'P':
            pick = arr.pick_id.get_referred_object()
            chan = pick.waveform_id.id
            if chan not in picked_stations.keys():
                picked_stations.update({chan: [event_id]})
            else:
                picked_stations[chan].append(event_id)

# Construct pick count series
pick_counts = pd.Series({_k: len(_v) for _k, _v in picked_stations.items()}).sort_values(ascending=False)
event_id_set = set(df_eb.event_id.values)

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
    Logger.info(f'\t{_pst} - {pick_counts[_pst]}')
# Get the 
missed_event_ids = event_id_set.difference(pref_sta_evid_set)
fmissed = len(missed_event_ids)/len(df_eb)
if fmissed > 0.1:
    Logger.warning(f'preferred stations list results in more than 10% event loss')
else:
    Logger.info(f'preferred stations list results in a {fmissed*100:.3f}% event loss ({len(missed_event_ids)})')