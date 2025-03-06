import os
from pathlib import Path
from glob import glob

import matplotlib.pyplot as plt

import pandas as pd
from obsplus import EventBank

import cartopy.crs as ccrs

from map_util import *


# Absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.parent
# path to eventbank
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# path to template dir
TMPD = ROOT / 'processed_data' / 'template' / 'single_station'
# path to catalog membership CSV
CATD = ROOT / 'processed_data' / 'catalog' / 'P1S1_Event_ID_Catalog_Membership.csv'

RAD_KM = 45.
ZOOM = 9

isshow = True
issave = False
SAVEPATH = ROOT / 'results' / 'figures' / 'seismolunch'
FMT = 'png'
DPI = 200



### PROCESSING SECTION ###
EBANK = EventBank(EBBP)
df_eb = EBANK.read_index()
df_pp = pd.read_csv(CATD)

# Merge with catalog membership
df_eb.index = df_eb.event_id
df_eb = df_eb.join(pd.read_csv(CATD, index_col='event_id'), how='left')

# Update index for merging with uw####### evid format
df_eb.index = [f'{x.split("/")[-2].lower()}{x.split("/")[-1]}' for x in df_eb.index]

# Iterate across templates to confirm if they had waveform data
tlist = glob(str(TMPD/'*'/'*'/'*'/'*.tgz'))
evid_sets = {'automatic': set(), 'manual': set()}
holder = []
rmap = {'manual': 1, 'automatic': -1}
for _f in tlist:
    parts = _f.split('/')
    evid = parts[-1].split('.')[0]
    year = int(parts[-2])
    revstat = parts[-3]
    stachan = parts[-4]
    evid_sets[revstat].add(evid)
    line = stachan.split('.') + [year, revstat, evid, rmap[revstat]]
    holder.append(line)

df_tpk = pd.DataFrame(holder, columns=['net','sta','loc','chan','year','revstat','evid','ohe'])

strictly_manual = evid_sets['manual'].difference(evid_sets['automatic'])
strictly_augmented = evid_sets['automatic'].difference(evid_sets['manual'])
mixed_representation = evid_sets['manual'].intersection(evid_sets['automatic'])
all_templated = evid_sets['manual'].union(evid_sets['automatic'])

df_eb = df_eb.assign(only_manual=df_eb.index.isin(strictly_manual))
df_eb = df_eb.assign(only_modeled=df_eb.index.isin(strictly_augmented))
df_eb = df_eb.assign(mixed=df_eb.index.isin(mixed_representation))
df_eb = df_eb.assign(hastemplate=df_eb.index.isin(all_templated))

### PLOTTING SECTION ###

## PLOT DIFFERENT CATALOGS
fig = plt.figure(figsize=(8,8))
ax, attr = mount_baker_basemap(fig=fig, radius_km = RAD_KM, zoom=ZOOM, latnudge=0, lonnudge=0)
add_rings(ax, rads_km=[10, 20, 30, 40], rads_colors=['k','k','k','k'])
mark_mount_baker(ax, labeled=False)

# CATn1 - 50 km
ax.plot(df_eb.longitude, df_eb.latitude, 'm.', ms=7, alpha=0.5, transform=ccrs.PlateCarree(), label='PNSN Catalog Events')

# CAT0 - 30km
df = df_eb[df_eb.CAT0]
ax.plot(df.longitude, df.latitude, 'k.', ms=6, transform=ccrs.PlateCarree(), label=f'& R $\leq$ 30 km ({len(df)})')
# CAT1 - 30km + 1980-onward
# df = df[(df.time >= pd.Timestamp('1980-01-01'))]
# ax.plot(df.longitude, df.latitude, '.', ms=5, color='violet', transform=ccrs.PlateCarree())
# CAT2 - CAT1 + generates template
df = df[df.hastemplate]
ax.plot(df.longitude, df.latitude, '.', ms=3, color='firebrick', transform=ccrs.PlateCarree(), label=f'& has waveforms ({len(df)})')
# CAT3 - CAT2 + has uncertainties
df = df[(df.vertical_uncertainty.notna()) & (df.horizontal_uncertainty.notna())]
ax.plot(df.longitude, df.latitude, '.', ms=4, color='goldenrod', transform=ccrs.PlateCarree(), label=f'& has uncertainties ({len(df)})')

# CAT4 - CAT3 + well constrained
df = df[df.WC]
ax.plot(df.longitude, df.latitude, '.', ms=2, color='cyan', transform=ccrs.PlateCarree(), label=f'& well constrained ({len(df)})')

ax.legend(loc='upper right')

gl = ax.gridlines(draw_labels=True, zorder=1)
gl.top_labels = False
gl.left_labels = False
gl.xlines = False
gl.ylines = False

if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'subarray_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)

## PLOT PICK DENSITY
fig = plt.figure(figsize=(8,4))
ax = fig.add_subplot(111)
pt = df_tpk.pivot_table(columns='evid', index=['net','sta'], values='ohe', fill_value=0)
ch = ax.pcolor(pt.values, cmap='seismic_r', vmin=-2, vmax=2)
ax.set_yticks(np.array(range(len(pt))) + 0.5, ['.'.join(x) for x in pt.index])
ax.set_xlabel('Event Index Value [No.]')

## 
if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'subarray_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)


## PLOT DIFFERENT CATALOGS
fig = plt.figure(figsize=(8,8))
ax, attr = mount_baker_basemap(fig=fig, radius_km = 100, zoom=8, latnudge=0, lonnudge=0)
add_rings(ax, rads_km=[10, 30, 50, 70, 90], rads_colors=['k','k','k','k','k'])
mark_mount_baker(ax, labeled=False)

# CAT0 - 30km
df = df_eb[df_eb.CAT0]
# CAT2 - CAT1 + generates template
df = df[df.hastemplate]
ax.plot(df.longitude, df.latitude, '.', ms=3, color='firebrick', transform=ccrs.PlateCarree(), label=f'& has waveforms ({len(df)})')

# ax.legend(loc='upper right')

from obspy.clients.fdsn import Client
client = Client('IRIS')

pts = np.abs(pt)
pts= pt.sum(axis=1)
for ind, val in pts.items():
    inv= client.get_stations(network=ind[0], station=ind[1], location='*', channel='*', level='station')
    for net in inv.networks:
        for sta in net.stations:
            ax.scatter(sta.longitude, sta.latitude, c='c', s=val**0.8, marker='v', edgecolor='k', cmap='inferno_r', vmin=0, vmax=550,
                       transform=ccrs.PlateCarree(), zorder=20)
            ax.text(sta.longitude - 0.02, sta.latitude + 0.02, '.'.join(ind), ha='right', va='bottom', transform=ccrs.PlateCarree())


gl = ax.gridlines(draw_labels=True, zorder=1)
gl.top_labels = False
gl.left_labels = False
gl.xlines = False
gl.ylines = False

if issave:
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    SAVENAME = SAVEPATH / f'prefstaeve_{DPI}dpi.{FMT}'
    plt.savefig(str(SAVENAME), format=FMT, dpi=DPI)


if isshow:
    plt.show()