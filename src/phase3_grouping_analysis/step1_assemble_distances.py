import logging, os, glob
from pathlib import Path

import numpy as np
import pandas as pd

from obsplus import EventBank
from obspy.geodetics import locations2degrees, degrees2kilometers

from eqcutil import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger

Logger = setup_terminal_logger(name=__name__, level=logging.INFO)

# Absolute path to repo root directory
ROOT = Path(__file__).parent.parent.parent
# Preclustered single-station clusters
CTRD = ROOT / 'processed_data' / 'cluster' / 'single_station'
# Base Path for Event Bank
EBBP = ROOT / 'processed_data' / 'catalog' / 'AUGMENTED_BANK'
# Save Path
SAVEPATH = ROOT / 'processed_data' / 'cluster' / 'tables'

EBANK = EventBank(EBBP)

clist = glob.glob(str(CTRD/'*.tgz'))
ctrd = {}
# Load all single-station, pre-clustered templates
for _c in clist:
    name = os.path.splitext(_c.split('/')[-1])[0]
    Logger.info(f'Loading {name}')
    ctrd.update({name: ClusteringTribe().read(_c)})

# Compose distance / corrmax / shift / timedelta dataframe
df_coh = pd.DataFrame()
evid_etype = {}
for _t, _ctr in ctrd.items():
    Logger.info(_t)
    # Get coherence and shift upper triangles
    coh = 1. - _ctr.dist_mat[np.triu_indices(len(_ctr), k=1)]
    shi = _ctr.shift_mat[np.triu_indices(len(_ctr), k=1)]
    # coh = 1. - _ctr.dist_mat.ravel()#[np.triu_indices(len(_ctr), k=1)]
    # shi = _ctr.shift_mat.ravel()#[np.triu_indices(len(_ctr), k=1)]
    # # Get a vector of trace labels
    labels = [_t.split('_')[0]]*len(shi)
    event_i, event_j = [], []
    # Create the combination of upper triangle
    for ii, evid_i in enumerate(_ctr._c.index):
        # evid_etype.update({evid_i: _ctr._c.loc[evid_i, 'etype']})
        for jj, evid_j in enumerate(_ctr._c.index):
            if ii < jj:
                event_i.append(evid_i)
                event_j.append(evid_j)
            
    df_hold = pd.DataFrame({'trace': labels, 'event_i':event_i, 'event_j': event_j, 'coh': coh, 'shift': shi})
    df_coh = pd.concat([df_coh, df_hold], axis=0, ignore_index=True)


# df_evid = pd.DataFrame(evid_etype, index=['etype']).T
# Get Time/distance table from events
df_eb = EBANK.read_index()
# Assign shorthand event IDs
df_eb = df_eb.assign(evid=[f"uw{x.split('/')[-1]}" for x in df_eb.event_id])

# Subset by present event IDs
included = df_coh.event_i.unique()
df_eb = df_eb[df_eb.evid.isin(included)]
df_eb.index = df_eb.evid

# df_eb = df_eb.join(df_evid, how='left')

df_del = pd.DataFrame()

# THIS CREATES THE UPPER TRIANGLE OF DISTANCE MATRICES
for ii, (evid_i, row_i) in enumerate(df_eb.iterrows()):
    Logger.info(f'{ii + 1} / {len(df_eb) - 1}')
    if ii + 1 < len(df_eb):
        lats = df_eb.latitude[ii+1:]
        herr = df_eb.horizontal_uncertainty[ii+1:]
        zerr = df_eb.vertical_uncertainty[ii+1:]
        lons = df_eb.longitude[ii+1:]
        times = df_eb.time[ii+1:]
        # etype_j = df_eb.etype[ii+1:]
        dx = degrees2kilometers(locations2degrees(row_i.latitude, row_i.longitude, lats, lons))*1e3
        dz = (row_i.depth - df_eb.depth[ii+1:])
        if np.isfinite(row_i.horizontal_uncertainty):
            dxs = herr**2 + row_i.horizontal_uncertainty**2
        else:
            dxs = herr**2
        if np.isfinite(row_i.vertical_uncertainty):
            dzs = zerr**2 + row_i.vertical_uncertainty**2
        else:
            dzs = zerr**2
        deltime = row_i.time - times
        event_j = deltime.index.values
        event_i = [evid_i]*len(event_j)
        dt = [x.total_seconds() for _, x in deltime.items()]
        df_hold = pd.DataFrame({'event_i': event_i,
                                'event_j': event_j,
                                'delt_ij_sec': dt,
                                'delh_ij_m': dx,
                                'sigh_ij_m2': dxs,
                                'delz_ij_m': dz,
                                'sigz_ij_m2': dzs})
        df_del = pd.concat([df_del, df_hold], axis=0, ignore_index=True)


# Write tables to disk
df_coh.to_csv(str(SAVEPATH/'coh_shift_table.csv'), header=True, index=False)
df_del.to_csv(str(SAVEPATH/'dist_table.csv'), header=True, index=False)
    # for jj, (evid_j, row_j) in enumerate(df_eb.iterrows()):
    #     if ii < jj:
    #         dist_x = degrees2kilometers(
    #             locations2degrees(
    #                 row_i.latitude,
    #                 row_i.longitude,
    #                 row_j.latitude,
    #                 row_j.longitude
    #                 )
    #             )
    #         dist_t = (row_j.time - row_i.time).total_seconds()
    #         line = [evid_i, evid_j, dist_x, dist_t]
    #         holder.append(line)

# Write coherence table to disk
# df_coh.to_csv(COH_SAVE, header=True, index=False)


#     for evid_i, row_i  in _ctr._c.iterrows():
#         Logger.info(evid_i)
#         for evid_j, row_j in _ctr._c.iterrows():
#             # Do upper triangle
#             if row_i.id_no < row_j.id_no:
#                 # Get max correlaton coefficient amplitude
#                 cohmax = 1. - _ctr.dist_mat[row_i.id_no, row_j.id_no]
#                 # Get max correlation coefficient amplitude shift
#                 shiftmax = _ctr.shift_mat[row_i.id_no, row_j.id_no]
#                 # Get event separation distance
#                 distkm = degrees2kilometers(locations2degrees(row_i.latitude, row_i.longitude, row_j.latitude, row_j.longitude))
#                 # Get event origint time difference
#                 distsec = abs(UTCDateTime(row_i.time) - UTCDateTime(row_j.time))
#                 # Trace ID, EVID_i, EVID_j, coh, dt_shift, distkm, dt_orig
#                 line = [_t, evid_i, evid_j, cohmax, shiftmax, distkm, distsec]
#             elif row_i.id_no == row_j.id_no:
#                 line = [_t, evid_i, evid_j, 1., 0., 0., 0.]
#             else:
#                 continue
#             holder.append(line)

# df_coh = pd.DataFrame(holder, columns=['trace','evid_i','evid_j','cohmax','shiftsec','distkm','distsec'])