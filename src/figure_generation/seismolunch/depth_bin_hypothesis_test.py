import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

from obspy.geodetics import locations2degrees, degrees2kilometers

from obsplus import EventBank
from sklearn.cluster import AgglomerativeClustering


from map_util import UTM10N, WGS84

# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to catalog membership CSV
CATD = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'
# path to preferred event/channel pairs CSV
PESD = ROOT / 'processed_data' / 'catalog' / 'P1S2_Preferred_Sta_Event_Picks.csv'
# path to analyst review summary CSV
REVD = ROOT / 'results' / 'survey' / 'S1_extracted_reviewer_classes.csv'

# SAVEPATH
SAVEPATH = ROOT / 'results' / 'figures' / 'seismolunch'
FMT = 'png'
DPI = 200

issave = True
isshow = False
## PREPROCESSING SECTION ###
# Connect to Event Bank
EBANK = EventBank(EBBP)
try:
    os.makedirs(str(REVD), exist_ok=False)
except:
    pass

# Read EBANK Index
df_eb = EBANK.read_index()
# Merge with catalog membership
df_eb.index = df_eb.event_id
df_eb = df_eb.join(pd.read_csv(CATD, index_col='event_id'), how='left')
# Create matplotlib timestamp-friendly times
df_eb = df_eb.assign(epoch=[x.timestamp() for x in df_eb.time])
# Get easting/northing
points = UTM10N.transform_points(
    x = df_eb.longitude.values,
    y = df_eb.latitude.values,
    src_crs = WGS84)
mE = points[:,0]
mN = points[:,1]
df_eb = df_eb.assign(mE = mE)
df_eb = df_eb.assign(mN = mN)

# Subset to 30 km radius
df_eb = df_eb[df_eb.CAT0]

df_del = pd.DataFrame()
# THIS CREATES THE UPPER TRIANGLE OF DISTANCE MATRICES
for ii, (evid_i, row_i) in enumerate(df_eb.iterrows()):
    # Logger.info(f'{ii + 1} / {len(df_eb) - 1}')
    if ii + 1 < len(df_eb):
        lats = df_eb.latitude[ii+1:]
        herr = df_eb.horizontal_uncertainty[ii+1:]
        zerr = df_eb.vertical_uncertainty[ii+1:]
        lons = df_eb.longitude[ii+1:]
        times = df_eb.time[ii+1:]
        etype_j = df_eb.etype[ii+1:]
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

dh_km = np.abs(df_del.delh_ij_m)*1e-3
sh_km = df_del.sigh_ij_m2**0.5*1e-3
ssmd_h = dh_km/sh_km
ssmd_h_finite = ssmd_h[np.isfinite(ssmd_h)]

dz_km = np.abs(df_del.delz_ij_m)*1e-3
sz_km = df_del.sigz_ij_m2**0.5*1e-3
ssmd_z = dz_km/sz_km
ssmd_z_finite = ssmd_z[np.isfinite(ssmd_z)]

dx_km = np.sqrt(df_del.delh_ij_m**2 + df_del.delz_ij_m**2)*1e-3
sx_km = np.sqrt(df_del.sigh_ij_m2 + df_del.sigz_ij_m2)*1e-3
ssmd_x = dx_km/sx_km
ssmd_x_finite = ssmd_x[np.isfinite(ssmd_x)]

fig = plt.figure(figsize=(1.25*5.5, 1.25*4.5))
ax = fig.add_subplot(111)
_n = sum(ssmd_z_finite>=2)/(len(ssmd_z_finite))
_m = 1 - len(ssmd_z_finite)/len(ssmd_z)
outs = ax.hist(ssmd_z, 3000, density=False, alpha=0.33, color='blue',
               label=f'Depth SSMD ({_n*100:.2f}% $\geq$ 2)\n{_m*100:.2f}% pairs omitted (no $\sigma_z$)')
_n = sum(ssmd_h_finite>=2)/(len(ssmd_h_finite))
_m = 1 - len(ssmd_h_finite)/len(ssmd_h)
ax.hist(ssmd_h, bins=outs[1], density=False, alpha=0.33, color='red',
        label=f'Epicentral SSMD ({_n*100:.2f}% $\geq$ 2)\n{_m*100:.2f}% pairs omitted (no $\sigma_h$)')
_n = sum(ssmd_x_finite>=2)/(len(ssmd_x_finite))
_m = 1 - len(ssmd_x_finite)/len(ssmd_x) 
ax.hist(ssmd_x, bins=outs[1], density=False, alpha=0.33, color='black',
        label=f'Hypocetral SSMD ({_n*100:.2f}% $\geq$ 2)\n{_m*100:.2f}% pairs omitted (no $\sigma_z$ | $\sigma_h$)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend(loc='upper right')
ax.set_ylabel('Frequency [count]')
ax.set_xlabel('Strictly Scaled Mean Difference (SSMD) [km/km]')

if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'catalog_ssmd_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)


# fig = plt.figure(figsize=(12, 8))
# ax1 = fig.add_subplot(311)
# ax2 = fig.add_subplot(312)
# ax3 = fig.add_subplot(313)

# ax1.hist2d(sh_km[np.isfinite(sh_km)], dh_km[np.isfinite(sh_km)], bins=(50,50), norm='log')
# ax2.hist2d(sz_km[np.isfinite(sz_km)], dz_km[np.isfinite(sz_km)], bins=(50,50), norm='log')
# ax3.hist2d(sx_km[np.isfinite(sx_km)], dx_km[np.isfinite(sx_km)], bins=(50,50), norm='log')

# for _e in bins:
#     for ax in [ax1, ax2, ax3]:
#         ax.plot([0,ax.get_xlim()[1]], [0, _e*ax.get_xlim()[1]], 'r-', linewidth=0.5)

# fig = plt.figure()
#         # ax.text(, f'{_e}-$\sigma$')

# ax1.scatter(np.abs(df_del.delh_ij_m)*1e-3, )
# ax2.scatter(np.abs(df_del.delz_ij_m)*1e-3,)
# ax3.scatter(, )
# ax1.set_xlabel('Epicentral Distance [km]')
# ax2.set_xlabel('Depth Difference [km]')
# ax3.set_xlabel('Hypocentral Distance [km]')

# ax1.set_ylabel('Epicentral Distance Uncertainty [km]')
# ax2.set_ylabel('Depth Difference Uncertainty [km]')
# ax3.set_ylabel('Hypocentral Distance [km]')
# ax2.set_ylabel('Epicentral Distance Uncertainty [km]')
plt.show()



# df_acsl = pd.DataFrame()
# for _d in range(1,11):
#     acsl = AgglomerativeClustering(linkage='complete', distance_threshold=_d*1e3, n_clusters=None)
#     acsl = acsl.fit(XYZ.values)
#     _ser = pd.Series(data=acsl.labels_, index=XYZ.index, name=f'd{_d:d}_km')
#     df_acsl = pd.concat([df_acsl, _ser], axis=1, ignore_index=False)

# df_acsl = pd.concat([df_acsl, XYZ], axis=1, ignore_index=False)