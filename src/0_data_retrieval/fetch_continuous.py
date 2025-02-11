import os, logging, glob
from pathlib import Path

from obspy import read_inventory, UTCDateTime, Stream
import pandas as pd
from obspy.clients.fdsn import Client
from obsplus import WaveBank


from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

Logger = setup_terminal_logger(os.path.split(__file__)[-1])
Logger.addHandler(CriticalExitHandler())
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to wavebank
WBBP = ROOT/'data'/'WF'/'BANK'
WBPS = '{year}/{julday}'
WBNS = '{seedid}_{year}{julday}_{hour}{minute}{second}'


DT = 300
DAY = 24*3600
t0 = UTCDateTime('2007-05-01T00:00:00')
tf = UTCDateTime('2007-09-01T00:00:00')
STAS = ['SHUK','MBW','MBW2','RPW','RPW2','JCW','SAXON']

client = Client("IRIS")
wbank = WaveBank(base_path=WBBP, path_structure=WBPS, name_structure=WBNS)

# breakpoint()
for _e, sta in enumerate(STAS):
    if _e == 0:
        INV = client.get_stations(network='UW', station=sta, location='*', channel='??Z', level='channel')
    else:
        INV += client.get_stations(network='UW', station=sta, location='*', channel='??Z', level='channel')
breakpoint()
while t0 + DT < tf:
    Logger.info(f'{t0}')
    # Iterate over stations

    inv = INV.select(time=t0)
    channels = inv.get_contents()['channels']
    for seed in channels:
        Logger.info(f'{seed}')
        bulk = []
        # Iterate over one day to build bulk
        for _e in range(int(DAY//DT)):
            t1 = t0 + _e*DT
            t2 = t1 + DT
            line = tuple(seed.split('.') + [t1, t2])
            bulk.append(line)
        
        st = Stream()
        try:
            st = client.get_waveforms_bulk(bulk)
            Logger.info('bulk success')
        except:
            Logger.debug('trying single request')
            for _b in bulk:
                try:
                    st += client.get_waveforms(*_b)
                    Logger.info('success')
                except:
                    Logger.debug("fail")
                    continue
        if len(st) > 0:
            st.split()
            try:
                wbank.put_waveforms(st)
                Logger.info('successful put')
            except:
                Logger.debug('trying single put')
                for tr in st:
                    try:
                        wbank.put_waveforms(tr)
                        Logger.info('successful single put')
                    except:
                        Logger.debug('single put fail')
                        continue
        else:
            Logger.debug('empty stream')
            pass
    t0 += DAY


    