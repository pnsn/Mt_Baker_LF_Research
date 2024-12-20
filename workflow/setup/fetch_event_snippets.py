import os, logging
from pathlib import Path

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obsplus import WaveBank, EventBank

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to wavebank
WBBP = ROOT/'data'/'WF'/'BANK'
WBPS = '{year}/{julday}'
WBNS = '{seedid}_{year}{julday}_{hour}{minute}{second}'

# State which station and channel(s) will be queried
NET = 'UW'
STA = 'JCW'
LOC = '*'
CHA = '??Z'
DT = 300.

# Connect to EventBank
EBANK = EventBank(base_path=EBBP)
# Connect to WaveBank
WBANK = WaveBank(base_path=WBBP,
                 path_structure=WBPS,
                 name_structure=WBNS)
# Connect to IRIS DMC Client
client = Client("IRIS")

# Get Event Metadata
df_eb = EBANK.read_index()

for _e, row in df_eb.iterrows():
    # Get origin time
    t0 = UTCDateTime(row.time.isoformat())
    Logger.info(f'{t0} - ({_e+1} of {len(df_eb)})')
    # Compose bulk request
    bulk = []
    for _e in range(-2,4):
        t1 = t0 + _e*DT
        t2 = t1 + DT
        bulk.append((NET, STA, LOC, CHA, t1, t2))
    # Create empty holder stream
    st = Stream()
    # Try to run bulk request
    try:
        st += client.get_waveforms_bulk(bulk)
        Logger.info('successful bulk query')
    # Failing that, try to run single entry request
    except:
        Logger.warning('attempting individual query')
        for _b in bulk:
            try:
                st += client.get_waveforms(*_b)
            except:
                continue
    # If data were queried
    if len(st) > 0:
        # Split to remove gaps
        st.split()
        # Submit to wavebank
        try:
            WBANK.put_waveforms(st)
            Logger.info('successful bulk submit to WaveBank')
        # Failing that, try to submit individual traces
        except:
            Logger.warning('attempting individual submit to WaveBank')
            for tr in st:
                try:
                    WBANK.put_waveforms(tr)
                except:
                    continue
    # If no data were queried, continue
    else:
        Logger.warning('empty stream')
        continue




    

