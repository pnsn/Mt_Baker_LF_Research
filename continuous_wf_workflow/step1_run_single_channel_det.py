import os
from pathlib import Path
from glob import glob
from eqcutil import ClusteringTribe

from obspy import UTCDateTime

from msiclient.client import WaveformClient

import numpy as np

import pandas as pd


ROOT = Path(__file__).parent.parent
PCSV = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
SSCT = ROOT/'processed_data'/'cluster'/'single_station'

client = WaveformClient(str(ROOT/'data'/'WF'/'baker_mseed.sqlite'))

# Set perferred NSLC's
PREFNSLC = ['UW.MBW..EHZ']#,'UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            # 'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            # 'UW.SHUK..BHZ','UW.SHUK..HHZ',
            # 'UW.JCW..EHZ']
MINSNR = 5
NSU = 4
# etype_mapper = {'lf': 'low frequency'}


# Load catalog profile
df_cp = pd.read_csv(PCSV, parse_dates=['prefor_time'], index_col=[0])


# Load single-channel clustering tribes
tribes = {}
for _k in PREFNSLC:
    print(f'loading {_k}')
    tribes[_k] = ClusteringTribe().read(str(SSCT/f'{_k}_clustered.tgz'))


df_cp_lf = df_cp[(df_cp.petype=='lf') &
                 (df_cp.prefor_time >= pd.Timestamp('2002-01-01'))]

df_cp_su = df_cp[(df_cp.petype=='su') & 
                 (df_cp.prefor_time >= pd.Timestamp('2002-01-01'))]

ctr = tribes['UW.MBW..EHZ']
cdf = ctr._c
# Filter & Sort
cdf_lf = cdf[(cdf.index.isin(df_cp_lf.index)) & (cdf.mean_snr_dB >= MINSNR)].sort_values('mean_snr_dB')
cdf_su = cdf[(cdf.index.isin(df_cp_su.index)) & (cdf.mean_snr_dB >= MINSNR)].sort_values('mean_snr_dB')

intrange = np.round(np.linspace(0, len(cdf_su)-1, NSU)).astype(int)

evids = list(cdf_lf.index.values)
for _e in intrange:
    evids.append(cdf_su.index.values[_e])

# Subset to all LFs and a few SUs with a range of SNRs
ctr_filt = ctr.get_subset(evids)

party = ctr_filt.client_detect(
    client, threshold=10,
    threshold_type='MAD',
    trig_int=10.,
    starttime=UTCDateTime('2007-01-08T00:00:00'),
    endtime=UTCDateTime('2007-01-12T12:00:00'), 
    save_progress=True, overlap='calculate')
breakpoint()

# # Filter for date, SNR, and etype
# df_cp_filt = df_cp[(df_cp.petype.isin(['lf','su']))&(df_cp.prefor_time >= pd.Timestamp('2002-01-01'))]    
# filt_tribes = {}
# for _k, _ct in tribes.items():
#     _ctf = _ct.get_subset(_ct._c[(_ct._c.index.isin(df_cp_filt.index)) & (_ct._c.mean_snr_dB >= MINSNR)].index)
#     filt_tribes[_k] = _ctf


