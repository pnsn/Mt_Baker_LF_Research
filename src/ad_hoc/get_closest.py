from glob import glob
from pathlib import Path

from obspy.clients.fdsn import Client
from obsplus import EventBank

from eqcutil import ClusteringTribe

import pandas as pd
# Get absolute path to repo root
ROOT = Path(__file__).parent.parent.parent.absolute()
# Directory path to unclustered templates
TEMD = ROOT / 'processed_data' / 'template' / 'single_station'
# Base Path for Augmented Event Bank
AEBBP = ROOT / 'processed_data' / 'catalog' / 'AUGMENTED_BANK'
# ETYPE CSV from AQMS
ETYPF = ROOT / 'data' / 'Events' / 'MtBaker_EVID_ETYPE.csv'


# Load ETYPE-EVID mapping CSV
df_ee = pd.read_csv(str(ETYPF))
# Connect to event bank
EBANK = EventBank(AEBBP)
# load event bank index
df_eb = EBANK.read_index()
# attach etype

# filter down to 30 km radius (This should be in CAT0)

# get xcorr filtering from preferred stations
 
"""
Analysis strategy
-----------------
We seek to have the (approximate) arrival time of waves from an
event at the closest observing station if that station is operative
as determined by successful template construction via QC measures from
PHASE1 and PHASE2/step1.

For Mount Baker this will be UW.MBW2 and UW.SHUK in the current record, but may be 
a splice of UW.MWB.[--,01].EHZ, UW.CMW..EHZ, UW.JCW..EHZ for deep cuts

"""
# Iterate across events
for event_id in df_eb.event_id:
    # Load QuakeML
    cat = EBANK.read_index( event_id=event_id)
    