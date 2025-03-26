import os
from pathlib import Path

import pandas as pd

import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

from obspy import UTCDateTime
from obspy.clients.fdsn import Client
from obspy.geodetics import locations2degrees

import shapely.geometry as sgeom
import cartopy.feature as cfeature
import cartopy.crs as ccrs
from map_util import *

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