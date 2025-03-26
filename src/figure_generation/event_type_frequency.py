from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from obsplus import EventBank
from obspy import UTCDateTime
from obspy.clients.fdsn import Client

ROOT = Path(__file__).parent.parent.parent
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
ETYPE = ROOT / 'data' / 'Events' / 'Mount_Baker_evids_etypes_10_JAN_2025.csv'
CATM = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'
EBANK = EventBank(EBBP)

LAT, LON, ELE = 48.7745,-121.8172, 3286

_styles = pd.read_csv(Path(__file__).parent / 'etype_styles.csv')
breakpoint()
# Get IRIS webservice client
IRIS = Client("IRIS")
# Get stations within 1 degree of Mount Baker
inv = IRIS.get_stations(longitude=LON,
                        latitude=LAT,
                        maxradius=75/111.2,
                        level='station',
                        network='UW,CN,TA')
inv30 = IRIS.get_stations(longitude=LON,
                        latitude=LAT,
                        maxradius=30/111.2,
                        level='station',
                        network='UW,CN,TA')

# Read event bank index
df = EBANK.read_index()
df.index = df.event_id
# Append catalog membership
df_cat = pd.read_csv(CATM, index_col='event_id')
df = pd.concat([df, df_cat], axis=1, ignore_index=False)

# # Assign time indices
# df.index = df.time
# # Sort by time
# df = df.sort_index()

df0 = df[df.CAT0]
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)

# Create frequency time-series for event types
_run_inv = True
for _e in ['eq','lf','su','px','uk']:
    _df = df0[df0.etype == _e]
    x = []
    y = []
    if _run_inv:
        z = []
        z30 = []
    for _yr in range(1980, 2026, 1):
        for _mo in range(1,12,1):
            t0 = pd.Timestamp(f'{_yr}-{_mo:02d}-01')
            t1 = pd.Timestamp(f'{_yr}-{_mo+1:02d}-01')
            ct = len(_df[(_df.time >= t0) & (_df.time <= t1)])
            y.append(ct)
            x.append(t1)
            if _run_inv:
                z.append(len(inv.select(time=UTCDateTime(f'{_yr}-{_mo:02d}-01')).get_contents()['stations']))
                z30.append(len(inv30.select(time=UTCDateTime(f'{_yr}-{_mo:02d}-01')).get_contents()['stations']))
    _run_inv = False
    ax.fill_between(x, [0]*len(y), y, label=f'"{_e.upper()}" Frequency', alpha=0.75)
    # _df = df0[df0.etype == _e].etype.rolling(pd.Timedelta(365, unit='days')).count()
    # ax.plot(_df.index, _df, label=_e)
# ax2 = ax.twinx()
ax.plot(x, z,'m',label='Station Count Within 75 km')
ax.plot(x,z30, 'k', label='Station Count Within 30 km')

plt.legend(loc='upper left')
