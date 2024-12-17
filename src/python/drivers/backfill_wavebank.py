import os, logging, argparse
from pathlib import Path

from obsplus import WaveBank, EventBank
from obspy.clients.fdsn.client import Client, FDSNNoDataException
from obspy import UTCDateTime

from eqcutil.util.logging import *

Logger = setup_terminal_logger(os.path.split(Path(__file__))[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(1))

### SET PATHS
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path's for WaveBank and EventBank 
WBBP = ROOT / 'data' / 'WF' / 'BANK'
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'




# Connect To Clients
client = Client("IRIS")
EBANK = EventBank(EBBP)
WBANK = WaveBank(WBBP)

