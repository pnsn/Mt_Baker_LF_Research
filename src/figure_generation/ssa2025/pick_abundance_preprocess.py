import os
from pathlib import Path

import pandas as pd
import numpy as np

from tqdm import tqdm
from obsplus import EventBank

from map_util import *

ROOT = Path(__file__).parent.parent.parent.parent
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
AEBBP = ROOT / 'processed_data' / 'catalog' / 'AUGMENTED_BANK'
TDIR = ROOT / 'processed_data' / 'template' / 'single_station'

# Connect to EventBanks
OBANK = EventBank(EBBP)
ABANK = EventBank(AEBBP)

event_ids = ABANK.read_index().event_id
# Iterate across all events
holder = []
for event_id in tqdm(event_ids):

    evid = f'uw{os.path.split(event_id)[1]}'

    # Get augmented event
    augmented_event = ABANK.get_events(event_id=event_id)[0]    
    # Get original event
    original_event = OBANK.get_events(event_id=event_id)[0]
    # Iterate across picks
    for pick in augmented_event.picks:
        emode = pick.evaluation_mode
        nslc = pick.waveform_id.id
        net, sta, loc, chan = nslc.split('.')
        oinst = any(_op.waveform_id.id[:-1] == nslc[:-1] for _op in original_event.picks)
        ochan = any(_op.waveform_id.id == nslc for _op in original_event.picks)

        year = pick.time.year
        time = pd.Timestamp(pick.time.timestamp, unit='s')
        viable = os.path.exists(TDIR/nslc/emode/str(year)/f'{evid}.tgz')
        
        # Pick exists
        code = 0
        # Pick made a template
        if viable:
            code += 1
            # # Template from a manual pick
            # if emode =='manual':
            #     code += 1
            # A manual pick on the original instrument
            if oinst:
                code += 1
                # # A manual pick on the original channel
                # if ochan:
                #     code += 1
                        
        # # Pick was a manual pick
        # if emode == 'manual':
        #     code += 1
        # # Pick was a manual pick on the original instrument
        # if oinst:
        #     code += 1
        # # Pick was the analyst picks
        # if ochan:
        #     code += 1

        line = [evid, time, nslc, net, sta, loc, chan, oinst, ochan, emode, viable, code]
        holder.append(line)

# Construct dataframe from accumulated metadata
df = pd.DataFrame(holder, columns=['evid','time','nslc','net','sta','loc','chan','original_inst','original_chan','eval_mode','made_template','emode'])
df.to_csv(ROOT/'results'/'tables'/'SSA2025'/'padded_template_emode.csv')

