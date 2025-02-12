import os, logging, warnings
from pathlib import Path
from collections import defaultdict

import pandas as pd
from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.core.event import Comment
from obsplus import EventBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.catalog.model_phases import model_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# FIXME: Ignore FutureWarning s
warnings.simplefilter(action='ignore', category=FutureWarning)

# Absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Absolute path to processed data directory
PDDIR = ROOT / "processed_data"
# Absolute path to preferred event-station-phase CSV
PSEDAT = PDDIR / "catalog" / "P1S2_Preferred_Sta_Event_Picks.csv"
# Absolute path to establish augmented catalog
AEBBP = PDDIR / "catalog" / "AUGMENTED_BANK"


LAT_REF = 48.7745   # [deg N] Reference point latitude
LON_REF= -121.8172  # [deg E] Reference point longitude
RAD_LIM_KM = 100.   # [km] Radius limit for station query
INV_LEVEL = 'channel' # Inventory query depth

# Pick Transfer Mapping
SHIFTS = {'N': 'Z',
          'E': 'Z'}

# Pick Modeling Parameters
MODPHZ = ['P','p']
VMOD = 'P5'
MODCHAN = '??Z'
ITYPE_PRIORITY = ['HH','BH','EH','HN','EN']


def main():
    try:
        os.makedirs(AEBBP)
    except FileExistsError:
        pass

    # Connect to source Event Bank
    EBANK = EventBank(EBBP)
    # Connect to IRIS webservice client
    IRIS = Client('IRIS')
    # Initialize Augmented Event Bank
    ABANK = EventBank(base_path=AEBBP,
                      path_structure='{year}',
                      name_structure='uw{event_id_end}')

    # Load preferred event-station-pick file
    df_pref = pd.read_csv(PSEDAT)

    # Get unique NSLC components
    unique_nslc = defaultdict(set)
    for _, row in df_pref.iterrows():
        for _f in ['network','station','channel']:
            unique_nslc[_f].add(row[_f])

    qkwargs = {_k: ','.join(list(_v)) for _k, _v in unique_nslc.items()}
    qkwargs.update({'longitude': LON_REF,
                    'latitude': LAT_REF,
                    'maxradius': RAD_LIM_KM/111.2,
                    'level':INV_LEVEL})
    # Query webservice for metadata
    INV = IRIS.get_stations(**qkwargs)

    # Get unique event_ids
    pref_evid = df_pref.event_id.unique()
    for evid in pref_evid:
        # Load event
        cat = EBANK.get_events(event_id=evid)
        # Ensure event is non-empty
        if len(cat) != 1:
            Logger.critical(f'fetched catalog with {len(cat)} events')
        else:
            event = cat[0]

        # Apply phase hints to Pick objects
        cat = apply_phase_hints(cat)
        # Filter picks 
        cat = filter_picks(
            cat,
            phase_hints=MODPHZ,
            enforce_single_pick='preferred',
            stations=list(unique_nslc['station']))
        breakpoint()

if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
    Logger.addHandler(CriticalExitHandler(exit_code=1))

    # RUN MAIN 
    df_out = main()