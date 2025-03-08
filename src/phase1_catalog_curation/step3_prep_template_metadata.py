"""

TODO:
Update this method to populate a pick for all valid channels per INV metadata
"""


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
VMOD = 'P4'
MODCHAN = '??Z'
ITYPE_PRIORITY = ['HH','BH','HN','EH','EN']


def main():
    # Make sure augmented bank base path directory exists
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

    STAS = unique_nslc['station']
    # Compile key-word arguments for Inventory query
    qkwargs = {_k: ','.join(list(_v)) for _k, _v in unique_nslc.items()}
    qkwargs.update({'longitude': LON_REF,
                    'latitude': LAT_REF,
                    'maxradius': RAD_LIM_KM/111.2,
                    'level':INV_LEVEL})
    # Query webservice for metadata
    INV = IRIS.get_stations(**qkwargs)

    for event_id, etype in df_pref[['event_id','etype']].value_counts().index:
        # Load event
        cat = EBANK.get_events(event_id=event_id)
        # Ensure event is non-empty
        if len(cat) != 1:
            Logger.critical(f'fetched catalog with {len(cat)} events')

        # Apply phase hints to Pick objects
        cat = apply_phase_hints(cat)
        # Filter picks 
        cat = filter_picks(
            cat,
            phase_hints=MODPHZ,
            enforce_single_pick='preferred',
            stations=list(unique_nslc['station']))
        
        event = cat[0]
        # Get picked channels (and apply component shift if applicable)
        # picksta = set()
        for _p in event.picks:
            _sta = _p.waveform_id.station_code
            _cha = _p.waveform_id.channel_code
            if _cha[-1] in SHIFTS.keys():
                newcode = _cha[:-1] + SHIFTS[_cha[-1]]
                _p.waveform_id.channel_code = newcode
                Logger.warning(f'Shifted {_p.waveform_id.id} channel to {newcode}')
            # picksta.add(_sta)

        # Insert ETYPE as comment on event
        event.comments.append(Comment(text=etype))
        # Get metadata from retained picks & apply channel shifts
        picksta = set([])
        picknslc = set([])
        kprid = []
        Logger.debug(f'Filtering preserved {len(event.picks)} picks')
        for _p in event.picks:
            kprid.append(_p.resource_id)
            picksta.add(_p.waveform_id.station_code)
            picknslc.add(_p.waveform_id.id)
            ccode = _p.waveform_id.channel_code
            # Apply component shifts as needed
            if ccode[-1] in SHIFTS.keys():
                scode = ccode[:-1] + SHIFTS[ccode[-1]]
                comment = Comment(text=f'Originally on channel {ccode}')
                Logger.warning(f'Shifting pick on {_p.waveform_id.id} to component {scode}')
                _p.waveform_id.channel_code = scode
                _p.comments.append(comment)
                
        # Station magnitudes cleanup
        smk = []
        for smag in event.station_magnitudes:
            if smag.waveform_id.station_code in STAS:
                smk.append(smag)
        Logger.debug(f'Removing {len(event.station_magnitudes) - len(smk)} unreferenced station_magnitudes')

        event.station_magnitudes = smk

        # Amplitudes cleanup
        amk = []
        for amp in event.amplitudes:
            if amp.waveform_id.station_code in STAS:
                amk.append(amp)
        Logger.debug(f'Removing {len(event.amplitudes) - len(amk)} unreferenced amplitudes')
        event.amplitudes = amk
        prefor = event.preferred_origin()

        # Arrivals cleanup
        ark = []
        for arr in prefor.arrivals:
            if arr.pick_id in kprid:
                ark.append(arr)
        Logger.debug(f'Removing {len(prefor.arrivals) - len(ark)} unreferenced arrivals')
        prefor.arrivals = ark

        # MODEL ARRIVAL TIMES
        # Get active station list
        inv = INV.select(time=prefor.time, channel=MODCHAN)
        # Get list of active NSLC codes at the time of this event
        active = inv.get_contents()['channels']
        # # Create a set of station codes to model using set logic
        # astas = set([_c.split('.')[1] for _c in active])
        # # Get stations to model for
        # modsta = list(set(astas).difference(picksta))
        # # Announce stations being modeled
        # msg = 'Modelling arrival times for station(s):'
        # for ms in modsta:
        #     msg += f' {ms},'
        # msg = msg[:-1]
        # Logger.debug(msg)

        modnslc = list(set(active).difference(picknslc))
        msg = 'Modelling arrival times for channel(s):'
        for mi in modnslc:
            msg += f' {mi},'
        msg = msg[:-1]
        Logger.debug(msg)
        

        modeled_first_arrivals = []
        # Iterate across each active station that is not picked
        # for sta in modsta:
        for nslc in modnslc:
            n,s,l,c = nslc.split('.')
            # Get the subset inventory for that station
            # sinv = inv.select(station=sta)
            # # Iterate across instrument/band priorities
            # for itype in ITYPE_PRIORITY:
            #     # Check sub-sub-inventory
            #     sbinv = sinv.select(channel=f'{itype}?')
            #     # If we find a component that has one represented channel
            #     # Break and continue to modeling
            #     if len(sbinv.get_contents()['channels']) == 1:
            #         break
            # # If we've gone across all instrument types and get nothing
            # if len(sbinv.get_contents()['channels']) != 1:
            #     # Continue to next station
            #     continue

            # Subset inventory to this active channel
            sbinv = inv.select(network=n, station=s, location=l, channel=c)

            # Model picks 
            # TODO: this can be updated to run for all NSLC to pick with recent update to eqcutil
            picks_hat = model_picks(
                prefor,
                sbinv,
                model_name=VMOD,
                phases=MODPHZ
            )
            # Take earliest pick
            t0 = UTCDateTime()
            for _ph in picks_hat:
                if _ph.time < t0:
                    t0 = _ph.time
                    earliest = _ph

            if len(picks_hat) > 0:
                modeled_first_arrivals.append(earliest)
            else:
                continue
    
        Logger.info(f'Adding {len(modeled_first_arrivals)} modeled first arrivals to {event_id}')
        # Append modeled picks to event (should go back to cat)
        pre_append_phz = len(event.picks)
        event.picks += modeled_first_arrivals
        if len(cat[0].picks) != pre_append_phz + len(modeled_first_arrivals):
            breakpoint()
        # try:
        ABANK.put_events(cat)
        Logger.info('Successfully added to Augmented Event Bank')
        # except:
            # breakpoint()   


    # for _p in event.picks:
    #     nslc = _p.waveform_id.id
    #     line = [event_id, _p.phase_hint, nslc, _p.evaluation_mode] + nslc.split('.')
    #     pick_info.append(line)
    # Try committing updated event to AUGMENTED EVENT BANK
    
if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
    Logger.addHandler(CriticalExitHandler(exit_code=1))

    # RUN MAIN 
    df_out = main()