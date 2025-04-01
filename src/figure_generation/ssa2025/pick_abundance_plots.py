import os
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm
from obsplus import EventBank

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

        oinst = any(_op.waveform_id.id[:-1] == nslc[:-1] for _op in original_event.picks)
        ochan = any(_op.waveform_id.id == nslc for _op in original_event.picks)

        year = pick.time.year
        datetime = pd.Timestamp(pick.time.timestamp, unit='s')
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
                # A manual pick on the original channel
                if ochan:
                    code += 1
                        
        # # Pick was a manual pick
        # if emode == 'manual':
        #     code += 1
        # # Pick was a manual pick on the original instrument
        # if oinst:
        #     code += 1
        # # Pick was the analyst picks
        # if ochan:
        #     code += 1

        line = [evid, datetime, nslc, oinst, ochan, emode, viable, code]
        holder.append(line)

#
df = pd.DataFrame(holder, columns=['evid','datetime','nslc','original_inst','original_chan','eval_mode','made_template','reality_code'])

pt = df[df.reality_code > 0].pivot_table(index=['nslc'], columns=['evid'], values=['reality_code'])

fig = plt.figure()
gs = fig.add_gridspec(ncols=1, nrows=3, wspace=0, hspace=0)
ax1 = fig.add_subplot(gs[:-1])
ax2 = fig.add_subplot(gs[-1], sharex=ax1)

ax1.pcolor(pt, cmap='cividis_r')
labels = []
for _e, _id in enumerate(pt.index):
    if _e == 0:
        labels.append(_id)
        _s = _id.split('.')[1]
    else:
        n,s,l,c = _id.split('.')
        if s == _s:
            labels.append(c)
        else:
            labels.append(_id)
        _s = s
ax1.set_yticks(np.arange(len(pt)) + 0.5, labels=labels)
ax1.xaxis.set_visible(False)

ax2.plot((pt > 0).sum(axis=0).values, 'k', linewidth=1)
ax2.fill_between(np.arange(len(pt.columns)), (pt > 0).sum(axis=0).values, color='#fde725')
ax2.fill_between(np.arange(len(pt.columns)), (pt > 1).sum(axis=0).values, color='#3b528b')
ax2.legend(['Total','Modeled','Catalog'])
ax2.yaxis.set_label_position('right')
ax2.yaxis.set_ticks_position('right')
ax2.set_ylabel('Template Channels [Ct.]', rotation=270, labelpad=15)
ax2.set_xlabel('Event Index (In Chronologic Order) [No.]')


# pt_exists = df.pivot_table(index=['nslc'], columns=['evid'], aggfunc=lambda x: 1, fill_value=0)
# pt_qcfail = df.pivot_table(index=['nslc'], columns=['evid'], values=['made_template'], aggfunc=lambda x: )