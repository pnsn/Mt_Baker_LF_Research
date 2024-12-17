import os, logging
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd
from obspy import Catalog, read_inventory
from eqcorrscan import Template, Tribe
from obsplus import EventBank, WaveBank

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path for EventBank 
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Path to event metadata
AQMS = ROOT / 'data' / 'Events' / 'MtBaker_50km_radius_origins.csv'
# Path to templating directory
TMPD = ROOT / 'processed_data' / 'bank_templates' / 'd20km'
# Path to INV
INVF = ROOT / 'data' / 'XML' / 'INV' / 'station_inventory_50km_MBS_IRISWS_STATION.xml'
# Path to CTR Save File
CTRF = TMPD / 'corr_cluster_XSTA_LockStep.tgz'


### CONNECT / LOAD ###
# Connect To Banks
EBANK = EventBank(EBBP)
# Load Inventory
INV = read_inventory(str(INVF))
# Load ClusteringTribe
CTR = ClusteringTribe().read(str(CTRF))
# Load Event Metadata
df_O = pd.read_csv(str(AQMS))
# Subset to evid and etype
df_O = df_O[['evid','etype']]
# Drop duplicates
df_O = df_O.drop_duplicates(keep='first')
# Generate matching index to CTR.clusters
df_O.index = [f'uw{_e}' for _e in df_O.evid]

# Get eventbank index
df_E = EBANK.read_index()
# Generate matching index to CTR.clusters
df_E.index = [os.path.splitext(os.path.split(row.path)[-1])[0] for _, row in df_E.iterrows()]

# Join etype to clusters
CTR.clusters = CTR.clusters.join(df_O['etype'],how='left')
CTR.clusters = CTR.clusters.join(df_E, how='left')

# Recluster at higher threshold
CTR.cct_regroup(0.6, inplace=True)

# Alias clusters
df_C = CTR.clusters

gcounts = df_C.correlation_cluster.value_counts()
min_members = 1
min_grp = 5
fig1 = plt.figure()
gs1 = fig1.add_gridspec(ncols=2, nrows=2)
axes = [fig1.add_subplot(gs1[_e]) for _e in range(4)]

for grp, count in gcounts.items():
    if count >= min_members:
        idf_C = df_C[df_C.correlation_cluster == grp]

        if count >= min_grp:
            marker = '.'
            zorder=count
            label=f'G{grp} ({count})'
        else:
            marker = 'xk'
            zorder=1
            label = None
        # axes[0].errorbar(idf_C.latitude, idf_C.longitude,
        #                  xerr=idf_C.horizontal_uncertainty.values/111.2e3,
        #                  yerr=idf_C.horizontal_uncertainty.values/111.2e3,
        #                  marker=marker, zorder=zorder, label=label)
        axes[0].plot(idf_C.longitude, idf_C.latitude, marker, label=label)
        axes[1].plot(idf_C.depth*1e-3, idf_C.latitude, marker)
        axes[2].plot(idf_C.longitude, idf_C.depth*1e-3, marker)
        axes[3].plot(idf_C.time, idf_C.depth*1e-3, marker)
        # axes[1].plot(idf_C.depth*1e-3, idf_C.vertical_uncertainty, marker)
        # axes[2].plot(idf_C.latitude, idf_C.longitude, marker)
        # axes[3].plot(idf_C.longitude,idf_C.depth,marker)

axes[0].legend()
    
fig2 = plt.figure()
ax = fig2.add_subplot(111)
CTR.dendrogram(xlabelfield='etype')

# fig3 = plt.figure()
# gs3 = fig3.add_gridspecc

