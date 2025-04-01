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
        net, sta, loc, chan = nslc.split('.')
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

        line = [evid, datetime, nslc, net, sta, loc, chan, oinst, ochan, emode, viable, code]
        holder.append(line)

#
df = pd.DataFrame(holder, columns=['evid','datetime','nslc','net','sta','loc','chan','original_inst','original_chan','eval_mode','made_template','reality_code'])
df = df.sort_values(by='datetime')
breakpoint()
evid_chrono = df.evid.unique()
evid_datetimes = [df[df.evid == evid].datetime.mean() for evid in evid_chrono]
evid_dates = [_t.strftime('%Y-%b') for _t in evid_datetimes]
# Create initial pivot table
pt = df.pivot_table(index=['nslc'], columns=['evid'], values=['reality_code'])
# Strip extraneous "reality code" from columns
pt = pt.T.loc['reality_code']
# Sort evids by chronology
pt = pt.loc[evid_chrono]
# Translate back
pt = pt.T
# Convert all 0's (did not make it) to NaN
pt[pt==0] = np.nan

prefnslc = ['UW.JCW..EHZ','UW.MBW..EHZ','UW.RPW..EHZ','UW.CMW..EHZ',
            'UW.MBW.01.EHZ','UW.RPW.01.EHZ','UW.SHUK..BHZ','CN.VDB..EHZ',
            'UW.MBW2..HHZ','UW.MBW2..ENZ','UW.RPW2..HHZ','UW.SHUK..HHZ','CN.VDEB..HHZ','UW.MULN..HHZ']

prefnlsc = ['UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.JCW..EHZ',
            'UW.CMW..EHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'CN.VDB..EHZ','CN.VDEB..HHZ',
            'UW.MULN..HHZ']

# net_order = ['UW','CN']
# sta_order = ['JCW','MBW','MBW2','RPW','RPW2','CMW','SHUK','SAXON','MCW','PASS','MULN',
#              'HTW','HDW','MBKE','VDEB','VDB','HOPB','PUBD','PLBD','SMSH']
# loc_order = ['','01']
# cha_order = ['HHZ','EHZ','BHZ','ENZ','HNZ']
# order = []
# labels = []
# for sta in sta_order:
#     for net in net_order:
#         for loc in loc_order:
#             for cha in cha_order:
#                 nslc = f'{net}.{sta}.{loc}.{cha}'

#                 if nslc in pt.index.values:
#                     if len(order) > 0:
#                         if order[-1].split('.')[1] == sta:
#                             labels.append(f'.{loc}.{cha}')
#                         else:
#                             labels.append(nslc)
#                     else:
#                         labels.append(nslc)
                
#                     order.append(nslc)

# year_index = []
# evids = df.evid.unique()
# evids.sort()
# for evid in evids:
#     dt = pd.Timestamp(df[df.evid==evid].datetime.values[0])
#     year_index.append(dt.strftime('%Y-%b'))

# pt = pt.loc[order]

# pt2 = pt.loc[prefnslc]

def plot_coverage(pivoted, year_index, labels=None):
    pt = pivoted

    fig = plt.figure(figsize=(10,8))
    gs = fig.add_gridspec(ncols=1, nrows=3, wspace=0, hspace=0)
    ax1 = fig.add_subplot(gs[:-1])
    ax2 = fig.add_subplot(gs[-1])

    ax1.pcolor(pt.values[::-1,:], cmap='seismic_r',vmin=0, vmax=4)
    if labels is None:
        ax1.set_yticks(np.arange(len(pt)) + 0.5, labels=pt.index.values[::-1])
    else:
        ax1.set_yticks(np.arange(len(pt)) + 0.5, labels=labels[::-1])
    ax1.xaxis.set_ticks_position('top')
    ax1.xaxis.set_label_position('top')
    tl = np.arange(0,100*(len(pt.T)//100 + 1),100, dtype=np.int32)
    ax1.set_xticks(tl, labels=[year_index[_t] for _t in tl])
    ax1.set_xlabel('Reference Date (Nonlinear Scaling)')
    ax1.grid(linestyle='-', alpha=0.2)

    ax2.plot((pt > 0).sum(axis=0).values, 'k', linewidth=1,label='Total')
    ax2.fill_between(np.arange(len(pt.columns)), (pt > 0).sum(axis=0).values, color='r',label='Modeled/Cloned')
    ax2.fill_between(np.arange(len(pt.columns)), (pt > 1).sum(axis=0).values, color='b',label='Catalog')
    ax2.legend()
    ax2.set_ylabel('Template Channels [Ct.]', rotation=90, labelpad=5)
    ax2.set_xlabel('Event Order [No.]\n Towards Past <----> Towards Present')
    ax2.set_ylim([0, 24])
    ax2.set_xlim([0, len(pt.columns)])
    ax2.grid(linestyle='-', alpha=0.2)

    return (fig, gs, ax1, ax2)

# pt_exists = df.pivot_table(index=['nslc'], columns=['evid'], aggfunc=lambda x: 1, fill_value=0)
# pt_qcfail = df.pivot_table(index=['nslc'], columns=['evid'], values=['made_template'], aggfunc=lambda x: )