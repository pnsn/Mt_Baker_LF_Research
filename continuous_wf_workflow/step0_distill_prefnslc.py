import os
from pathlib import Path
from collections import defaultdict
from tqdm import tqdm
from glob import glob
from eqcutil import ClusteringTribe
from obspy.core.event import Comment

from eqcorrscan import Template
from eqcorrscan.utils.stacking import align_traces, linstack

import pandas as pd


ROOT = Path(__file__).parent.parent
PCSV = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
SSCT = ROOT/'processed_data'/'cluster'/'single_station'
# Set perferred NSLC's
PREFNSLC = ['UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'UW.JCW..EHZ']
MINSNR = 1.2
# etype_mapper = {'lf': 'low frequency'}

xcckwargs = {
    'shift_len': 5.,
    'allow_individual_trace_shifts': False,
    'corr_thresh': 0.2}

# Load catalog profile
df_cp = pd.read_csv(PCSV, parse_dates=['prefor_time'], index_col=[0])

df_cp_lfsu = df_cp[(df_cp.petype.isin(['lf','su']))&(df_cp.prefor_time >= pd.Timestamp('2002-01-01'))]
# Load single-channel clustering tribes
tribes = {}
for _k in PREFNSLC:
    print(f'loading {_k}')
    tribes[_k] = ClusteringTribe().read(str(SSCT/f'{_k}_clustered.tgz'))

# Get all deep LF events'
templates_LF = defaultdict(list)
templates_SU = defaultdict(list)

# Assemble by evid & apply SNR filtering
for _k, _ct in tribes.items():
    olen = len(_ct)
    clusters = _ct.clusters
    # Get intermediate length filtering for SU/LF & post 2001
    ilen = len(clusters.index.isin(df_cp_lfsu.index))
    csub = clusters[(clusters.mean_snr_dB >= MINSNR) & (clusters.index.isin(df_cp_lfsu.index))]
    _ctf = _ct.get_subset(csub.index)
    print(f'filtered {_k} for post 2001 SU/LF and SNR g.e. {MINSNR}: {olen} --> {ilen} -> {len(_ctf)}')
    for _t in _ctf:
        petype = df_cp_lfsu.loc[_t.name].petype
        if petype == 'lf':
            templates_LF[_t.name].append(_t)
        elif petype == 'su':
            templates_SU[_t.name].append(_t)

# Compose multi-channel templates
newtribes = defaultdict(ClusteringTribe)
for etype, tmpset  in zip(['lf','su'],[templates_LF, templates_SU]):
    for evid, templates in tmpset.items():
        for _e, _tmp in enumerate(templates):
            if _e == 0:
                new_tmp = _tmp.copy()
                if new_tmp.event.comments[0].text != etype:
                    new_tmp.event.comments.append(Comment(text=etype))
                    # new_tmp.event.event_type = etype_mapper[etype]
            else:
                new_tmp.event.picks += _tmp.event.picks.copy()
                new_tmp.event.origins[0].arrivals += _tmp.event.origins[0].arrivals.copy()
                new_tmp.st += _tmp.st.copy()
        if len(new_tmp.st) > 0:
            newtribes[etype] += new_tmp
    # Run Clustering
    newtribes[etype].cluster('xcc', **xcckwargs)


        


