from pathlib import Path

import matplotlib.pyplot as plt

import pandas as pd
import numpy as np

from obspy.imaging.spectrogram import spectrogram
from obspy.clients.fdsn import Client

from obsplus import EventBank
# Repository Root Absolute Path
ROOT = Path(__file__).parent.parent.parent.parent
# Catalog Event Bank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Catalog Profile CSV
CPCSV = ROOT / 'processed_data' / 'catalog' / 'P1S1_Catalog_Profile.csv'
# Save Path
SPATH = ROOT / 'results' / 'figures' / 'SSA2025'
# Rendering Settings
DPI = 120
FMT = 'png'
# Output control
issave = True
isshow = True


### SPECTROGRAM INPUTS ###
pref_netsta = ['UW','MBW']

event_id_fstr = 'quakeml:uw.anss.org/Event/UW/{evid}'

pref_evid = {'eq': '10243683',
             'lf': '61822721',
             'px': '10800098',
             'su': '61802931'}

# EVENT SYTLE FORMATTING
marker_map = {'su':'d', 'px':'*','lf':'s','eq':'o'}
label_map = {'eq':'EQ: Earthquake',
             'lf':'LF: Low Frequency Event',
             'su':'SU: Surface Event',
             'px':'PX: Probable Blast'}
mec_map = {'su':'b', 'px':'m', 'eq':'k','lf':'r'}
lw_map = {'su': 1, 'px': 1, 'eq': 1, 'lf': 1}
depth_cmap = 'Greens'
base_alpha=0.8
dbins = np.arange(0, 35, 5)

# Waveform Query Arguments
prepick=10
duration = 60

# Spectrogram kwargs
spect_kw = {'per_lap': 0.95,
            'wlen': 2.56,
            'log': True,
            'dbscale': False,
            'cmap': 'Spectral_r'}


# CONNECT TO CLIENTS
ebank = EventBank(EBBP)
client = Client('IRIS')


# Compose bulk waveform request
picks = {}
for _etype, _evid in pref_evid.items():
    cat = ebank.get_events(event_id=event_id_fstr.format(evid=_evid))
    for _p in cat[0].picks:
        if _p.waveform_id.id.split('.')[:2] == pref_netsta:
            if _etype not in picks.keys():
                picks[_etype] = _p
            elif _p.time < picks[_etype].time:
                picks[_etype] = _p

bulk = {_e: (_p.waveform_id.id.split('.') + [_p.time - prepick, _p.time + duration]) for _e, _p in picks.items()}

# Get waveforms
st = client.get_waveforms_bulk(list(bulk.values()))

# Associate traces to etypes and assess completeness
traces = {}
prepads = {}
durations = {}
for _etype, pick in picks.items():
    for tr in st:
        if tr.stats.starttime < pick.time < tr.stats.endtime:
            traces[_etype] = tr
            prepads[_etype] = picks[_etype].time - tr.stats.starttime
            durations[_etype] = tr.stats.endtime - picks[_etype].time


### PLOTTING SECTION ###
            
### SPECTROGRAMS ###
figs = {}
axes = {}
caxes = {}
for _etype, _tr in traces.items():
    __tr = _tr.copy().detrend('linear').normalize()
    evid = pref_evid[_etype]
    pick = picks[_etype]
    fig = plt.figure(figsize=(5.25,5.25))
    gs = fig.add_gridspec(ncols=1, nrows=2, hspace=0, wspace=0)
    axs = [fig.add_subplot(gs[0])]
    axs.append(fig.add_subplot(gs[1], sharex=axs[0]))

    figs[_etype] = fig
    axes[_etype] = axs

    axs[0].plot(__tr.times(reftime=__tr.stats.starttime),
                __tr.data,
                mec_map[_etype],
                lw=0.5)
    shdl = spectrogram(__tr.data, __tr.stats.sampling_rate, axes=axs[1], **spect_kw)
    caxes[_etype] = shdl

    # axs[0].xaxis.set_visible(False)
    axs[1].set_xticks([_e + prepads[_etype] for _e in range(-10, 70, 10)], labels=[f'{_e:d}' for _e in range(-10, 70, 10)])


    xmin = -5 + prepads[_etype] + 2.56
    xmax = 30 + prepads[_etype] + 2.56
    axs[0].set_xlim([xmin,xmax])
    axs[0].set_yticks([-1, 0, 1], labels=[-1, 0, 1])
    axs[0].set_ylim([-1.2, 1.2])
    axs[1].set_ylim([0.5,30])

    # Add Labels
    axs[0].set_title(f'uw{evid} | ETYPE: {_etype.upper()}\n{pick.time}')
    axs[0].set_ylabel('Normalized Amplitude [ - ]')
    axs[0].grid(alpha=0.5, which='major')
    axs[0].text(20, __tr.data.min(), pick.waveform_id.id, fontsize=12, ha='left', va='bottom')

    axs[1].set_yticks([0.5,1,2,3,5,10,20,30], [0.5,1,2,3,5,10,20,30])
    axs[1].set_ylabel('Frequency [Hz]',labelpad=0)
    axs[1].set_xlabel(f'Elapsed Time Since P-Wave Arrival [sec]')
    axs[1].grid(alpha=0.5)

            # req = _p.waveform_id.id.split('.') + [_p.time - prepad, _p.time + duration]
            # req = dict(zip(['network','station','location','channel','starttime','endtime'], req))
            # traces.update({_etype: client.get_waveforms(**req)})

if issave:
    for _etype, _fig in figs.items():
        _fig.savefig(str(SPATH/f'wave_spect_UW_MBW_{_etype}_uw{pref_evid[_etype]}_v5_{DPI}dpi.{FMT}'), dpi=DPI, format=FMT)    
