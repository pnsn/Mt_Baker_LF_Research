"""
:program: Mt_Baker_LF_Research/workflow/templating/step1_station_selection.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:

    This program conducts the following tasks
    
    TASK 1
    Subset earthquake metadata hosted in an ObsPlus ventBank based on their epicentral
    distance from a reference location (Mount Baker) and origin time
    
    TASK 2
    Iterate across events and log which stations have picks for each event. 

    TASK 3
    Filter stations by the number of P-picks they have on any channel and save this
    information to a CSV

    TASK 4
    Assess which stations necessary to have at least one catalog pick per event by
    successively including stations in decending pick-count order until all events have
    at least one observation or the list from TASK 3 is exhausted. Generate a figure
    displaying these results.

    Outputs are:
     - preferred_event_sta_picks.csv: CSV of  event/NSLC/phase-type for
        station/event combinations that meet the 
     - step1_station_selection_###dpi.png: Figure from Task 4

"""

import os, logging, glob
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
import pandas as pd
from obsplus import EventBank
from obspy.geodetics import locations2degrees
from obspy import read_inventory, Inventory

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

### START OF USER INPUT SECTION ###
# Catalog filtering parameters
LAT_REF = 48.7745   # [deg N] Reference point latitude
LON_REF= -121.8172  # [deg E] Reference point longitude
ELEV_REF = 3286.    # [m ASL] Reference point elevation (not used)
RAD_LIM_KM = 30.    # [km] Selection radius
MINPICK=50          # [ct] Minmum number of picks to consider for continuity analysis
# Start date time-based filtering
STARTDATE = pd.Timestamp('2001-01-01T00:00:00')

# Well Located Criteria
no_fixed = True
max_loc_err_m = 10e3
max_rms_s = 1.
min_obs = 6



# Inventory subsetting parameters
NETS = ['UW']
CHANS = '[BHE][HN][ZNE123]'

# Figure rendering parameters
DPI = 120             # [dots per inch] saved figure resolution
FMT = 'PNG'           # File format
plot_save_kwargs = {} # Additional key-word arguments to add to matplotlib.pyplot.savefig
display_after = True
### END OF USER INPUT SECTION ###

### START OF PROCESSING SECTION ###
Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to station response files
INVD = ROOT / 'data' / 'XML' / 'RESP'
# Output Path
OUTD = ROOT / 'processed_data' / 'workflow' / 'catalog'
# Preferred Station Picked Events File (save file)
PSPEF = OUTD / 'preferred_station_picked_event_ids.csv'


# Load inventory, subsetting for desired channels
INV = Inventory()
for _e, _f in enumerate(glob.glob(str(INVD/'*.xml'))):
    inv = read_inventory(_f)
    for net in NETS:
        INV += inv.select(network=net, channel=CHANS)

## INITIALIZE EVENT BANK
EBANK = EventBank(EBBP)


# Load event summary
df_eb = EBANK.read_index()
df_eb = df_eb.sort_values(by='time')
Logger.info(f'EventBank has a total of {len(df_eb)} events')
# Calculate distances from Mt. Baker Summit
df_eb = df_eb.assign(radius_offset_km=[
    111.2* locations2degrees(LAT_REF, LON_REF, row.latitude, row.longitude)
    for _, row in df_eb.iterrows()])

### ANNOUNCE SUBSETTING IMPACTS ON EVENT COUNTS ###
# Test Subset Events by Date
newcount = len(df_eb[df_eb.time >= STARTDATE])
Logger.info(f'Subsetting to {STARTDATE} to PRESENT retains {newcount} events')

# Test Subset Events by Distance
newcount = len(df_eb[df_eb.radius_offset_km <= RAD_LIM_KM])
Logger.info(f'Subsetting to {RAD_LIM_KM} km distance retains {newcount} events')

# Apply Subsetting
df_eb_sub = df_eb[(df_eb.time >= STARTDATE) & (df_eb.radius_offset_km <= RAD_LIM_KM)]
newcount = len(df_eb_sub)
Logger.info(f'Subsetting by both retains {newcount} events')

### CONDUCT CONTINUITY/COVERAGE CHECK
if Logger.level < 20:
    _dtqdmf = True
else:
    _dtqdmf = False

picked_channels = set()
pick_info = []
well_located_evids = []

_e = -1
for event_id in tqdm(df_eb_sub.event_id, disable = _dtqdmf):
    # Increment up indexer
    _e += 1
    Logger.debug(f"{event_id} ({_e+1} of {newcount})")
    # Fetch event metadata
    cat = EBANK.get_events(event_id=event_id)
    # Get preferred origin metadata
    prefor = cat[0].preferred_origin()

    # Assess if this is a "well-constrained event"
    wcstatus = True
    row = df_eb[df_eb.event_id==event_id]
    # Origin solution fixed value checks
    if no_fixed:
        if prefor.time_fixed:
            if wcstatus:
                Logger.debug(f'{event_id} - fixed time')
            wcstatus &= False
            
        if prefor.epicenter_fixed:
            if wcstatus:
                Logger.debug(f'{event_id} - fixed epi')
            wcstatus &= False

        if prefor.depth_type == 'from location':
            pass
        else:
            if wcstatus:
                Logger.debug(f'{event_id} - not "from location" depth')
            wcstatus &= False
    # Origin horizontal uncertainty checks
    if row.horizontal_uncertainty.values[0] > max_loc_err_m:
        if wcstatus:
            Logger.debug(f'{event_id} - horizontal uncertainty')
        wcstatus &= False
    elif not np.isfinite(row.horizontal_uncertainty.values[0]):
        if wcstatus:
            Logger.debug(f'{event_id} - NaN horizontal uncertainty')
        wcstatus &= False
    # Origin vertical uncertainty checks
    if row.vertical_uncertainty.values[0] > max_loc_err_m:
        if wcstatus:
            Logger.debug(f'{event_id} - vertical uncertainty')
        wcstatus &= False
    elif not np.isfinite(row.vertical_uncertainty.values[0]):
        if wcstatus:
            Logger.debug(f'{event_id} - NaN vertical uncertainty')
        wcstatus &= False
    # Time uncertainty checks
    if prefor.time_errors.uncertainty is None:
        if wcstatus:
            Logger.debug(f'{event_id} - no time errors - keeping')
        # wcstatus &= False
    elif prefor.time_errors.uncertainty > max_rms_s:
        if wcstatus:
            Logger.debug(f'{event_id} - min_rms')
        wcstatus &= False
    # Observation count check
    if row.used_phase_count.values[0] < min_obs:
        if wcstatus:
            Logger.debug(f'{event_id} - min_obs')
        wcstatus &= False


    if wcstatus:
        well_located_evids.append([event_id, prefor.depth_type])
        

    # Get active channels for this event
    inv = INV.select(time=prefor.time)
    # Check if it has any associated picks
    if len(prefor.arrivals) == 0:
        Logger.warning(f'No associated picks for {event_id} - skipping to next')
        continue
    for arr in prefor.arrivals:
        pick = arr.pick_id.get_referred_object()
        nslc = pick.waveform_id.id
        line = [event_id, arr.phase, nslc] + nslc.split('.')
        pick_info.append(line)


# Create Picked Dataframe    
df_picked = pd.DataFrame(pick_info, columns=['event_id','phase','nslc','net','sta','loc','chan'])
# Subset to only channels in inventory
df_picked = df_picked[df_picked.nslc.isin(INV.get_contents()['channels'])]
# Subset to only P picked stations
ser_top_p = df_picked[df_picked.phase=='P'].sta.value_counts()
# Subset to stations with GE MINPICK P picks
top_p = ser_top_p[ser_top_p >= MINPICK].index.values

# Create Availability array
df_p_cont = df_picked[(df_picked.phase=='P') & (df_picked.sta.isin(top_p))][['event_id','sta']]
df_p_cont = df_p_cont.pivot_table(index=['event_id'], columns=['sta'], aggfunc=lambda x: 1, fill_value=0)
# Translate and sort channels by number of picks
df_p_cont = df_p_cont.T.loc[top_p]


### FIGURE RENDERING SECTION ###
well_located_evids = pd.DataFrame(well_located_evids, columns=['event_id','depth_type'])
# From TASK1 Figure - Effects of time, distance, and origin quality filtering
fig= plt.figure(figsize=(6,6))
gs = fig.add_gridspec(ncols=2, nrows=2)
axes = [fig.add_subplot(gs[_e]) for _e in range(4)]

for _lbl in ['Catalog','Radius','Radius+Time','R+T+Well-Constrained']:
    if _lbl == 'Catalog':
        df_e = df_eb
        fmt = 's'
        ms = 1
        color='black'
    elif _lbl == 'Radius':
        df_e = df_eb[df_eb.radius_offset_km <= RAD_LIM_KM]
        fmt = 'v'
        ms = 2
        color='firebrick'
    elif _lbl == 'Radius+Time':
        df_e = df_eb[(df_eb.radius_offset_km <= RAD_LIM_KM) & 
                     (df_eb.time >= STARTDATE)]
        fmt = 'x'
        ms = 3
        color='dodgerblue'
    elif _lbl == 'R+T+Well-Constrained':
        df_e = df_eb[df_eb.event_id.isin(well_located_evids.event_id)]
        fmt = '*'
        ms = 4
        color='goldenrod'

    plotfields = [('time','depth'),     ('longitude','latitude'),
                  ('time','magnitude'), ('time','vertical_uncertainty')]
                   
    for _b, (xfld, yfld) in enumerate(plotfields):
        if yfld in ['depth','vertical_uncertainty']:
            _yv = df_e[yfld].values*1e-3
            if yfld == 'depth':
                _yv *= -1
        else:
            _yv = df_e[yfld].values
        _xv = df_e[xfld].values
        axes[_b].plot(_xv, _yv, fmt, ms=ms, color=color, label=f'{_lbl} ({len(_xv)})')
        axes[_b].set_xlabel(xfld)
        axes[_b].set_ylabel(yfld)
    axes[0].legend()

plt.savefig(str(OUTD / f'step1_event_selection_{DPI:d}dpi.{FMT.lower()}'),
            dpi=DPI,format=FMT, **plot_save_kwargs)

# Render Figure for Pick Availability & Data Informed Station Selection
fig = plt.figure(figsize=(8,6))
gs = fig.add_gridspec(ncols=1, nrows=2, wspace=0, hspace=0)
axes = [fig.add_subplot(gs[0]), fig.add_subplot(gs[1])]
axes[0].pcolor(df_p_cont)
axes[0].set_yticks(np.arange(len(df_p_cont.index))+0.5, df_p_cont.index.values)
axes[0].set_title(f'Pick Availability for Stations with $\geq$ {MINPICK} Picks')
for _e in range(3, len(top_p) + 1):
    nobs = df_p_cont.iloc[:_e,:].sum(axis=0)
    nobs_bins = nobs.value_counts()
    if 0 in nobs_bins.index:
        n0 = nobs.value_counts()[0]
    else:
        n0 = 0
    axes[1].fill_between(list(range(len(nobs))), [0]*len(nobs), nobs.values, label=f'Top {_e} ({n0} no-obs)',zorder=20-_e)
    if n0 > 0:
        for _j, _n in enumerate(nobs.values):
            if _n == 0:
                axes[0].plot([_j, _j],[0,_e], 'r-', alpha=0.2)
    if n0 == 0:
        axes[0].fill_between([0, len(nobs)],[_e]*2, [len(df_p_cont)]*2,color='k',alpha=0.5)
        axes[0].text(len(nobs)//2, np.mean([_e, len(df_p_cont)]), 'Redundant', color='w')
        # If this is the first station that completes event converage, break iteration
        break
axes[1].set_xlabel('Event Index (--Towards Present-->)')
axes[0].set_ylabel('Station Code (<--Increasing Total Picks--)')
axes[1].set_ylabel('Number of Analyst Picks [ct.]')
axes[1].legend(ncols=1)
axes[1].set_xlim(axes[0].get_xlim())

plt.savefig(str(OUTD / f'step1_station_selection_{DPI:d}dpi.{FMT.lower()}'),
            dpi=DPI,format=FMT, **plot_save_kwargs)

# Save df_picked for non-redundant stations (or full set if there are still gaps)
df_out = df_picked[df_picked.sta.isin(top_p[:_e])]
# Annotate if this entry corresponds to a well-constrained event
df_out = df_out.assign(WCE=[x in well_located_evids.event_id.values for x in df_out.event_id])
# Save to disk
df_out.to_csv(str(OUTD/'preferred_event_sta_picks.csv'), header=True, index=False)


if display_after:
    plt.show()