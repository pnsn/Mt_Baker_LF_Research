from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client


from map_util import *

ROOT = Path(__file__).parent.parent.parent.parent
INPUT = ROOT/'results'/'figures'/'SSA2025'/'padded_template_emode.csv'



prefnslc = ['UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.JCW..EHZ',
            'UW.CMW..EHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'CN.VDB..EHZ','CN.VDEB..HHZ',
            'UW.MULN..HHZ']

df = pd.read_csv(INPUT, parse_dates=['time'])
# Sort chronologically
df = df.sort_values(by='time')
# Get list of EVIDs in chronologic order
evid_chrono = df.evid.unique()
# Get average times of EVID data (we only use year and month here, so seconds don't matter)
evid_times = [df[df.evid == evid].time.mean() for evid in evid_chrono]
# Convert into Year-Mo. strings for plotting
evid_dates = [_t.strftime('%Y-%b') for _t in evid_times]
# Create initial pivot table for channel/evid/
pt = df.pivot_table(index=['nslc'], columns=['evid'], values=['emode'])
# Strip extraneous "emode" multi-index layer from columns
pt = pt.T.loc['emode']
# Sort evids by time using ordered list from earlier
pt = pt.loc[evid_chrono]
# Translate back to columns being EVIDs
pt = pt.T
# Convert all 0's (did not make it) to NaN
pt[pt==0] = np.nan
# Subset to the preferred channels
pt2 = pt.loc[prefnslc]

### PLOT PICK CONTINUITY
fig = plt.figure(figsize=(10,8))
gs = fig.add_gridspec(ncols=1, nrows=3, wspace=0, hspace=0)
ax1 = fig.add_subplot(gs[:-1])
ax2 = fig.add_subplot(gs[-1])

mycmap, norm = make_pnsn_cmap(pallet_names=['lime','navy'], discretization=[1,1.5,2])
pallet = pnsn_pallet()
ax1.pcolor(pt2.values[::-1,:], cmap=mycmap,norm=norm)
ax1.set_yticks(np.arange(len(pt2)) + 0.5, labels=prefnslc[::-1])
ax1.xaxis.set_ticks_position('top')
ax1.xaxis.set_label_position('top')
tl = np.arange(0,100*(len(pt2.T)//100 + 1),100, dtype=np.int32)
ax1.set_xticks(tl, labels=[evid_dates[_t] for _t in tl])
ax1.set_xlabel('Reference Date (Nonlinear Scaling)')
ax1.grid(linestyle='-', alpha=0.2)

ax2.plot((pt2 > 0).sum(axis=0).values, 'k', linewidth=1,label='Total',zorder=4)
ax2.fill_between(np.arange(len(pt2.columns)),
                 (pt2 > 1).sum(axis=0).values,
                 color=pallet['navy'],
                 label='Catalog',
                 zorder=3)
ax2.fill_between(np.arange(len(pt2.columns)),
                 (pt2 > 0).sum(axis=0).values,
                 color=pallet['lime'],
                 label='Modeled/Cloned',
                 zorder=2)

ax2.legend()
ax2.set_ylabel('Template Channels [Ct.]', rotation=90, labelpad=5)
ax2.set_xlabel('Event Order [No.]\n Towards Past <----> Towards Present')
ax2.set_ylim([0, 9])
ax2.set_xlim([0, len(pt2.columns)])
ax2.grid(linestyle='-', alpha=0.2)


### PLOT STATIONS WITH PICK ABUNDANCE UNDERLAYS
inv = Client('IRIS').get_stations(latitude=BAKER_LAT, longitude=BAKER_LON,
                                  maxradius=1,level='channel',
                                  channel='BHZ,HHZ,EHZ,HNZ,ENZ')

# Spool up reference df

holder = []
for net in inv.networks:
    for sta in net.stations:
        for cha in sta.channels:
            nslc_list = [net.code, sta.code, cha.location_code, cha.code]
            nslc = '.'.join(nslc_list)
            line = [nslc] + nslc_list + [sta.latitude, sta.longitude]
            holder.append(line)
df_inv = pd.DataFrame(holder, columns=['nslc','net','sta','loc','chan','lat','lon'])

df_inv.drop_duplicates(keep='first', inplace=True)

PP = pnsn_pallet()
netcolors = {'UW':PP['navy'],
             'CC':'royalblue',
             'CN':'firebrick',
             'GS':'olive',
             'TA':'darkgoldenrod'}
            #  'NP':'yellow'}
NETS = ','.join(list(netcolors.keys()))

fig = plt.figure(figsize=(7,7))
gs = fig.add_gridspec(ncols=1,nrows=1)

axm, map_attrib = mount_baker_basemap(
    fig=fig, sps=gs[0], radius_km=75, open_street_map=False,
    zoom=8)

# Plot Baker
plot_baker(axm, zorder=100)

# Plot stations
df_out = df_inv[~df_inv.net.isin(netcolors.keys())]
axm.scatter(df_out.lon, df_out.lat, marker='v', c=pnsn_pallet()['lime'], s=4,
            linewidths=0.5, transform=ccrs.PlateCarree())

df_np = df_inv[(df_inv.net.isin(netcolors.keys())) &
               (~df_inv.nslc.isin(prefnslc))]
axm.scatter(df_np.lon, df_np.lat, marker='v', c=pnsn_pallet()['evergreen'], 
            s=9,
            linewidths=0.5, transform=ccrs.PlateCarree())

df_in = df_inv[(df_inv.nslc.isin(prefnslc))]
axm.scatter(df_in.lon, df_in.lat, marker='v', c=pnsn_pallet()['forest green'],
            s= 25, linewidths=0.5, transform=ccrs.PlateCarree())

# Underlay pie charts
# For each channel included in this study
for _, row in df_in.iterrows():
    # Subset pivot table
    nobs = pt.loc[row.nslc].count()
    

# pt_exists = df.pivot_table(index=['nslc'], columns=['evid'], aggfunc=lambda x: 1, fill_value=0)
# pt_qcfail = df.pivot_table(index=['nslc'], columns=['evid'], values=['made_template'], aggfunc=lambda x: )



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
