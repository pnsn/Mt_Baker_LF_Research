import os, logging, glob
from pathlib import Path

import numpy as np
import pandas as pd

from obspy import read_inventory, Inventory
from obspy.clients.fdsn import Client
from obspy.core.event import ResourceIdentifier, WaveformStreamID, Comment
from eqcorrscan import Tribe
from eqcorrscan.utils.catalog_utils import filter_picks
from obsplus import EventBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.catalog.model_phases import model_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
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
# Event Type Metadata
ETDATA = ROOT / 'data' / 'Events' / 'Mount_Baker_evids_etypes_10_JAN_2025.csv'
### USER CONTROLS ###

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(exit_code=1))

STAS = ['MBW','MBW2','SHUK','RPW','RPW2','JCW','CMW','SAXON','MULN','PASS']
min_snr = 1.1
min_noise_RMS = 3
phases=['P','p']
velocity_model = 'P4'

SHIFTS = {'UW.SHUK..HHN':'UW.SHUK..HHZ',
          'UW.MBW2..HHE':'UW.MBW2..HHZ',
          'UW.MBW2..ENZ':'UW.MBW2..HHZ',
          'UW.MULN..HHN':'UW.MULN..HHZ',
          'CN.VDB..SHZ':'CN.VDB..EHZ'}
# Location codes to remove in order to permit records at the same
# site from different instruments to be cross-correlated
ALIASES = {'UW.MBW.01.EHZ':'UW.MBW..EHZ',
           'UW.RPW.01.EHZ':'UW.RPW..EHZ'}

pick_filt_kwargs = {'phase_hints': phases,
                    'enforce_single_pick': 'preferred',
                    'stations': STAS}


def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))

def main():
    # INITIALIZE CLIENTS
    EBANK = EventBank(EBBP)
    # Ensure that augmented event bank has a home
    if not os.path.exists(str(OUTD/'AUGMENTED_BANK')):
        os.makedirs(str(OUTD/'AUGMENTED_BANK'))
    # Initialize Augmented event bank
    ABANK = EventBank(base_path=OUTD/'AUGMENTED_BANK',
                      path_structure='{year}',
                      name_structure='uw{event_id_end}')

    # Load inventory
    INV = Inventory()
    for _f in glob.glob(str(INVD/'*.xml')):
        # Read the given inventory
        inv = read_inventory(_f)
        # Include only vertical channels for preferred stations
        for _sta in STAS:
            INV += inv.select(station=_sta,channel='[BHE][NH][Z]')

    # Read preferred event IDs
    df_evids = pd.read_csv(str(PSPEF))
    df_evids = df_evids.assign(evid=[int(_e.split('/')[-1]) for _e in df_evids.event_id])
    # Load etype metadata
    df_meta = pd.read_csv(ETDATA, index_col=[0])
    # Join dataframes
    df_evids = df_evids.join(df_meta, on='evid', how='left')
    # Assign index
    df_evids.index = df_evids.event_id

    # Iterate across events
    for evid, meta in df_evids.iterrows():
        Logger.info(f'Processing {evid} ({meta.etype})')
        # Fetch event from event bank
        cat = EBANK.get_events(event_id=evid)
        # Apply phase hints
        cat = apply_phase_hints(cat)
        # Subset to selected stations and phases only
        cat = filter_picks(cat, **pick_filt_kwargs)
        # Create a comment for etype
        comment = Comment(text=meta.etype)
        cat[0].comments.append(comment)
        # # Apply channel shifts
        # for event in cat:
        #     for pick in event.picks:
        #         if pick.waveform_id.id in SHIFTS.keys():
        #             pick.waveform_id = WaveformStreamID(
        #                 seed_string=SHIFTS[pick.waveform_id.id]
        #             )
        # Get active channels
        active = INV.select(time=cat[0].preferred_origin().time)
        # Get picked channels
        picked = []
        kept_arrivals = []
        for _arr in cat[0].preferred_origin().arrivals:
            pick = _arr.pick_id.get_referred_object()
            if not pick is None:
                picked.append(pick.waveform_id.id)
                kept_arrivals.append(_arr)
        if kept_arrivals == []:
            breakpoint()
        cat[0].preferred_origin().arrivals = kept_arrivals

        
        # Get unpicked channels
        unpicked = set(active.get_contents()['channels']).difference(set(picked))
        Logger.info(f'{len(picked)} picked channels')
        Logger.info(f'{len(unpicked)} unpicked channels')
        Logger.info(f'Modeling for phases {phases}')
        # Model arrival times for unpicked channels
        for nslc in unpicked:
            n, s, l, c = nslc.split('.')
            iactive = active.select(network=n,station=s, location=l, channel=c)
            picks_hat = model_picks(cat[0].preferred_origin(),
                                    iactive,
                                    model_name=velocity_model,
                                    phases=phases)
            for ph in picks_hat:
                cat[0].picks.append(ph)
                Logger.debug(f'Modeled pick "{ph.phase_hint}" for stream "{ph.waveform_id.id}"')
        # Save subset, augmented catalog to EVENT BANK
        Logger.info('Submitting processed/augmented event to new eventbank')
        ABANK.put_events(cat)
        









if __name__ == '__main__':
    main()