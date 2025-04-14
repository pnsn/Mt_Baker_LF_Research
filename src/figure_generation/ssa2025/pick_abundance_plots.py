from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from obspy.clients.fdsn import Client

import cartopy.crs as ccrs
import map_util as mutil

ROOT = Path(__file__).parent.parent.parent.parent
INPUT = ROOT/'results'/'tables'/'SSA2025'/'padded_template_emode.csv'
CPD = ROOT/'results'/'tables'/'SSA2025'/'catalog_profile.csv'
SHPDIR = ROOT/'data'/'SHP'
NFSHP = SHPDIR/'USDA'/'Forest_Administrative_Boundaries_(Feature_Layer)'/'Forest_Administrative_Boundaries_(Feature_Layer).shp'


SDIR = ROOT/'results'/'figures'/'SSA2025'
issave = True
isshow = True
DPI = 120
FMT = 'png'
plt.rcParams.update({'font.size':14})
# Create key word argument holder for all fig.savefig calls
sfkw = {'format': FMT, 'dpi': DPI, 'bbox_inches':'tight'}

# Marker Rendering by Event Type
marker_map = {'px':'*','su':'d','lf':'s','eq':'o'}
color_map = {'px':'m','su':'b','lf':'r','eq':'k'}
ungrouped_color = 'xkcd:ugly brown'
unanalyzed_color = 'xkcd:gross green'

prefnslc = ['UW.MBW..EHZ','UW.MBW.01.EHZ','UW.MBW2..HHZ','UW.MBW2..ENZ',
            'UW.RPW..EHZ','UW.RPW.01.EHZ','UW.RPW2..HHZ',
            'UW.JCW..EHZ',
            'UW.CMW..EHZ',
            'UW.SHUK..BHZ','UW.SHUK..HHZ',
            'CN.VDB..EHZ','CN.VDEB..HHZ',
            'UW.MULN..HHZ']



# Load National Forest Boundaries
gdf_mbsnf = mutil.gpd.read_file(NFSHP)
# Subset to Mount Baker-Snoqualmie National Forest ShapeFile contents
gdf_mbsnf = gdf_mbsnf[gdf_mbsnf.FORESTNAME=='Mt. Baker-Snoqualmie National Forest']
# Load preprocessed catalog CSV
df_cat = pd.read_csv(CPD, parse_dates=['prefor_time'], index_col=[0])

# Load precompiled pick inventory
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
fig = plt.figure(figsize=(8,6))
axc = fig.add_subplot(111)



# Plot continuity
# axc = fig.add_subplot(sps_coverage)
axc.pcolor(pt2.values[::-1,:], cmap='cool_r', vmin=0.5, vmax=2.5)
axc.set_yticks([])
# axc.set_yticks(np.arange(len(pt2)) + 0.5, labels=prefnslc[::-1], fontsize=12, direction='in')
lbls = prefnslc[::-1]
for _e, _y in enumerate(np.arange(len(pt2)) + 0.5):
    axc.text(20, _y, lbls[_e], fontsize=14, ha='left', va='center', fontweight='extra bold')
axc.tick_params(axis='y',direction='in')
axc.set_xlabel('Subcatalog Event Number in Sequence')

if issave:
    fig.savefig(str(SDIR/f'Selected_Stachan_Survival_Plot_{DPI}dpi.{FMT}'),
                **sfkw)



fig = plt.figure(figsize=(15,8))
gs = fig.add_gridspec(ncols=2, nrows=1)
sps_event_map = gs[0]
sps_station_map = gs[1]

# Plot selected events
# Initialize geoaxes and get AWS hillshade basemap (very light)
axem, attre = mutil.mount_baker_basemap(
    fig=fig, sps=sps_event_map, radius_km=31,
    latnudge=0, lonnudge=0,
    open_street_map=False,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})

# Add gridliner
gl = axem.gridlines(draw_labels=['top','left'], zorder=1,
                   alpha=0, ylabel_style={'fontsize':12,'rotation':90},
                   xlabel_style={'fontsize':12})

# Add distance Rings
mutil.add_rings(axem,
                rads_km=[10,20,30],
                rads_colors=['k','k','xkcd:hot purple']*3,
                include_units=[True, True, True],
                label_pt=18,ha='left',va='top', fontsize=12)

axem.plot([-123, -121], [49,49], 'k-', transform=ccrs.PlateCarree())

# # Add Mount Baker-Snoqualmie NF boundary
hdl = mutil.plot_gdf_contents(axem, gdf_mbsnf, transform=None, alpha=0.2, linewidth=2, edgecolor='forestgreen')

# Plot Mount Baker & friends
axem.scatter([mutil.BAKER_LON, mutil.SHUKS_LON, mutil.STWIN_LON],
            [mutil.BAKER_LAT, mutil.SHUKS_LAT, mutil.STWIN_LAT],
            s=[144, 81, 81], 
            marker='^',
            facecolor='none',
            edgecolors='orange',
            linewidths=[2, 1, 1],
            zorder=30,
            transform=ccrs.PlateCarree())

# Plot events over stations for this map
_zo = 11
for _et in ['eq','lf','su','px']:
    _dfc = df_cat[(df_cat.tidy_group > -2)&(df_cat.etype==_et)]
    axem.scatter(_dfc.lon, _dfc.lat,
                 marker=marker_map[_et],
                 c=color_map[_et],
                 s=(3**_dfc.mag) + 6,
                 alpha=0.6,
                 zorder=_zo,
                 transform=ccrs.PlateCarree())
    _zo += 1    





# Initialize geoaxes and get AWS hillshade basemap (very light)
axsm, attrs = mutil.mount_baker_basemap(
    fig=fig, sps=sps_station_map, radius_km=61,
    open_street_map=False,
    latnudge=-0.15,
    aws_add_image_kwargs={'cmap':'Greys_r', 'alpha':0.05})

# Add gridliner
gl = axsm.gridlines(draw_labels=['top','left'], zorder=1,
                   alpha=0, ylabel_style={'fontsize':12,'rotation':90},
                   xlabel_style={'fontsize':12})

# axsm.add_feature(mutil.cfeature.BORDERS, linewidth=2, alpha=0.5)
axsm.add_feature(mutil.cfeature.OCEAN, alpha=0.5, color='dodgerblue')

# Add distance Rings
mutil.add_rings(axsm,
                rads_km=[10,30,50,70],
                rads_colors=['k','xkcd:hot purple','k','k'],
                include_units=[True, True, True, True],
                label_pt=60,ha='left',va='top', fontsize=12)
axsm.plot([-123, -121], [49,49], 'k-', transform=ccrs.PlateCarree())

# # Add Mount Baker-Snoqualmie NF boundary
hdl = mutil.plot_gdf_contents(axsm, gdf_mbsnf, transform=None, 
                              alpha=0.2, linewidth=2, edgecolor='forestgreen')

# Plot Mount Baker & friends
axsm.scatter([mutil.BAKER_LON, mutil.SHUKS_LON, mutil.STWIN_LON],
            [mutil.BAKER_LAT, mutil.SHUKS_LAT, mutil.STWIN_LAT],
            s=[144, 81, 81], 
            marker='^',
            facecolor='none',
            edgecolors='orange',
            linewidths=[2, 1, 1],
            zorder=30,
            transform=ccrs.PlateCarree())

inv = Client("IRIS").get_stations(network='UW,CN',station='MBW,MBW2,RPW,RPW2,JCW,CMW,SHUK,VDB,VDEB,MULN',level='station')
_holder = []
for net in inv.networks:
    for sta in net.stations:
        _line = [net.code, sta.code, sta.latitude, sta.longitude]
        _holder.append(_line)
df_inv = pd.DataFrame(_holder, columns=['net','sta','lat','lon'])


axsm.scatter(df_inv.lon, df_inv.lat, marker='v', edgecolors='k', s=7**2,
             c='cyan', zorder=50, transform=ccrs.PlateCarree())
for _, row in df_inv.iterrows():
    if row.sta == 'MBW':
        axsm.text(row.lon-0.01, row.lat-0.01, f'{row.net}.{row.sta}',
            fontsize=12, transform=ccrs.PlateCarree(),
            ha='right', va='top')
    elif row.sta == 'MBW2':
        axsm.text(row.lon-0.01, row.lat+0.01, f'{row.net}.{row.sta}',
            fontsize=12, transform=ccrs.PlateCarree(),
            ha='right', va='bottom')
    elif row.sta == 'VDB':
        axsm.text(row.lon-0.01, row.lat+0.01, f'{row.net}.{row.sta}',
            fontsize=12, transform=ccrs.PlateCarree(),
            ha='right')
    elif row.sta == 'RPW2':
        axsm.text(row.lon+0.01, row.lat-0.01, f'{row.net}.{row.sta}',
                fontsize=12, transform=ccrs.PlateCarree(),
                va='top')
    else:    
        axsm.text(row.lon+0.01, row.lat+0.01, f'{row.net}.{row.sta}',
                fontsize=12, transform=ccrs.PlateCarree())
if issave:
    fig.savefig(str(SDIR/f'Selected_Stations_And_Events_{DPI}dpi.{FMT}'), **sfkw)

if isshow:
    plt.show()
# tl = np.arange(0,150*(len(pt2.T)//200 + 1),200, dtype=np.int32)
# axc.set_xticks(tl, labels=[evid_dates[_t] for _t in tl])
# axc.set_xlabel('Reference Date (Nonlinear Scaling)')
# axc.grid(linestyle='-', alpha=0.2)
# xlims = axc.get_xlim()
# ylims = axc.get_ylim()
# axc.fill_between([-5,-4],[-4,-4],[-3,-3], color='b',label='Catalog Picks')
# axc.fill_between([-5,-4],[-4,-4],[-3,-3], color='r',label='Modeled Picks')
# axc.set_xlim(xlims)
# axc.set_ylim(ylims)
# axc.set_legend(loc='lower left')

# ax1 = fig.add_subplot(gs[:-1])
# ax2 = fig.add_subplot(gs[-1])






# mycmap, norm = make_pnsn_cmap(pallet_names=['xkcd:highligher green','navy'], discretization=[1,1.5,2])
# pallet = pnsn_pallet()
# # ax1.pcolor(pt2.values[::-1,:], cmap=mycmap,norm=norm)
# ax1.pcolor(pt2.values[::-1,:], cmap='seismic_r', vmin=0.5, vmax=2.5)
# ax1.set_yticks(np.arange(len(pt2)) + 0.5, labels=prefnslc[::-1])
# ax1.xaxis.set_ticks_position('top')
# ax1.xaxis.set_label_position('top')
# tl = np.arange(0,100*(len(pt2.T)//100 + 1),100, dtype=np.int32)
# ax1.set_xticks(tl, labels=[evid_dates[_t] for _t in tl])
# ax1.set_xlabel('Reference Date (Nonlinear Scaling)')
# ax1.grid(linestyle='-', alpha=0.2)

# ax2.plot((pt2 > 0).sum(axis=0).values, 'k', linewidth=1,label='Total',zorder=4)
# ax2.fill_between(np.arange(len(pt2.columns)),
#                  (pt2 > 1).sum(axis=0).values,
#                  color=pallet['navy'],
#                  label='Catalog',
#                  zorder=3)
# ax2.fill_between(np.arange(len(pt2.columns)),
#                  (pt2 > 0).sum(axis=0).values,
#                  color=pallet['xkcd:highligher green'],
#                  label='Modeled/Cloned',
#                  zorder=2)

# ax2.legend()
# ax2.set_ylabel('Template Channels [Ct.]', rotation=90, labelpad=5)
# ax2.set_xlabel('Event Order [No.]\n Towards Past <----> Towards Present')
# ax2.set_ylim([0, 9])
# ax2.set_xlim([0, len(pt2.columns)])
# # ax2.grid(linestyle='-', alpha=0.2)


# if issave:
#     breakpoint()


# if isshow:
#     plt.show()



# ### PLOT STATIONS WITH PICK ABUNDANCE UNDERLAYS
# inv = Client('IRIS').get_stations(latitude=BAKER_LAT, longitude=BAKER_LON,
#                                   maxradius=1,level='channel',
#                                   channel='BHZ,HHZ,EHZ,HNZ,ENZ')

# # Spool up reference df

# holder = []
# for net in inv.networks:
#     for sta in net.stations:
#         for cha in sta.channels:
#             nslc_list = [net.code, sta.code, cha.location_code, cha.code]
#             nslc = '.'.join(nslc_list)
#             line = [nslc] + nslc_list + [sta.latitude, sta.longitude]
#             holder.append(line)
# df_inv = pd.DataFrame(holder, columns=['nslc','net','sta','loc','chan','lat','lon'])

# df_inv.drop_duplicates(keep='first', inplace=True)

# PP = pnsn_pallet()
# netcolors = {'UW':PP['navy'],
#              'CC':'royalblue',
#              'CN':'firebrick',
#              'GS':'olive',
#              'TA':'darkgoldenrod'}
#             #  'NP':'yellow'}
# NETS = ','.join(list(netcolors.keys()))

# fig = plt.figure(figsize=(7,7))
# gs = fig.add_gridspec(ncols=1,nrows=1)

# axm, map_attrib = mount_baker_basemap(
#     fig=fig, sps=gs[0], radius_km=75, open_street_map=False,
#     zoom=8)

# # Plot Baker
# plot_baker(axm, zorder=100)

# # Plot stations
# df_out = df_inv[~df_inv.net.isin(netcolors.keys())]
# axm.scatter(df_out.lon, df_out.lat, marker='v', c=pnsn_pallet()['xkcd:highligher green'], s=4,
#             linewidths=0.5, transform=ccrs.PlateCarree())

# df_np = df_inv[(df_inv.net.isin(netcolors.keys())) &
#                (~df_inv.nslc.isin(prefnslc))]
# axm.scatter(df_np.lon, df_np.lat, marker='v', c=pnsn_pallet()['evergreen'], 
#             s=9,
#             linewidths=0.5, transform=ccrs.PlateCarree())

# df_in = df_inv[(df_inv.nslc.isin(prefnslc))]
# axm.scatter(df_in.lon, df_in.lat, marker='v', c=pnsn_pallet()['forest green'],
#             s= 25, linewidths=0.5, transform=ccrs.PlateCarree())

# # Underlay pie charts
# # For each channel included in this study
# for _, row in df_in.iterrows():
#     # Subset pivot table
#     nobs = pt.loc[row.nslc].count()
    

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
