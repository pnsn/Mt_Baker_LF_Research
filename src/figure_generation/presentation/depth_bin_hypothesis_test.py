import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

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
offset_max_km = 30.
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

# Subset events to 30 km 
df = df_eb[df_eb.CAT0]
# Further subset to have vertical uncertainties
df2 = df[df.vertical_uncertainty.notna()]
# Further subset to have both uncertainties
df3 = df2[df2.horizontal_uncertainty.notna()]
# Get meter-scaled positions
XYZ = df[['mE','mN','depth']]

df_acsl = pd.DataFrame()
for _d in range(1,11):
    acsl = AgglomerativeClustering(linkage='single', distance_threshold=_d, n_clusters=None)
    acsl = acsl.fit(XYZ.values())
    _ser = pd.Series(data=acsl.labels_, index=XYZ.index, name=f'{_d:d}_km')
    df_acsl = pd.concat([df_acsl, _ser], axis=1, ignore_index=False)
