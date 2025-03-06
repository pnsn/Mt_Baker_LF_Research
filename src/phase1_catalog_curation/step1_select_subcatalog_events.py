"""
:program: Mt_Baker_LF_Research/src/phase1_catalog_curation/step1_select_stations.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:

    This program conducts the following tasks
    
    TASK 1
    Subset earthquake metadata hosted in an ObsPlus EventBank based on their epicentral
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

import os, logging
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
from obsplus import EventBank
from obspy.geodetics import locations2degrees

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

### START OF USER INPUT SECTION ###

 # Absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get event type metadata
ETDATA = ROOT / "data" / "Events" / "MtBaker_EVID_ETYPE.csv"
# Aboslute path to processed data / catalog directory
PDDIR = ROOT / "processed_data" 
# Catalog Membership Output
SAVENAME = PDDIR / "catalog" / "P1S1_Event_ID_Catalog_Membership.csv"

# Catalog filtering parameters
LAT_REF = 48.7745   # [deg N] Reference point latitude
LON_REF= -121.8172  # [deg E] Reference point longitude
ELEV_REF = 3286.    # [m ASL] Reference point elevation (not used)
RAD_LIM_KM = 30.    # [km] Selection radius

# Start-date-based filtering
DIGITAL_STARTDATE = pd.Timestamp('1980-01-01')
CONTINUOUS_STARTDATE = pd.Timestamp('2001-01-01')

# Well Constrained Event Criteria
no_fixed = True         # No fixed status flags on origin components
max_loc_err_m = 10e3    # [m] maximum location error (horizontal or vertical)
max_rms_s = 1.          # [seconds] maximum RMS misfit
min_obs = 6             # [no.] minimum number of phase observations for the preferred origin
min_sta = 4.            # [no.] minimum number of unique stations observing the preferred origin
max_closest_m = 10e3    # [m] largest closest station distance to accept as well constrained

def main():
    """MAIN PROCESS - DEFINE PATH PARAMETERS AND LOGGER AS GLOBALS"""

    ## INITIALIZE EVENT BANK CONNECTION
    EBANK = EventBank(EBBP)

    ## Event Metadata Processing
    # Load event summary
    df_eb = EBANK.read_index()
    # Sort values by origin time
    df_eb = df_eb.sort_values(by='time')
    # Assign event_id as index
    df_eb.index = df_eb.event_id
    # Calculate distances from Mt. Baker Summit
    df_eb = df_eb.assign(radius_offset_km=[
        111.2* locations2degrees(LAT_REF, LON_REF, row.latitude, row.longitude)
        for _, row in df_eb.iterrows()])
    
    # Load event metadata
    df_meta = pd.read_csv(ETDATA, index_col=[0])
    # Attach etype to df_eb
    df_eb = df_eb.assign(evid=[int(row.event_id.split('/')[-1]) for _, row in df_eb.iterrows()])
    df_eb = df_eb.join(df_meta, on='evid', how='left')
    # Iterate across events
    _e = -1
    catalog_status = defaultdict(list)
    for event_id, row in df_eb.iterrows():
        _e += 1
        Logger.debug(f'{event_id} ({_e+1} of {len(df_eb)})')
        # Get event_id
        catalog_status['event_id'].append(event_id)
        # Get event type
        catalog_status['etype'].append(row.etype)
        # Get preferred origin time (used for filtering)
        catalog_status['prefor_time'].append(row.time)
        # Get radius (derived value used for filtering)
        catalog_status['offset_km'].append(row.radius_offset_km)
        # CAT0 Membership (within 30 km of Mount Baker)
        CAT0 = row.radius_offset_km <= RAD_LIM_KM
        catalog_status['CAT0'].append(CAT0)
        # CAT1 Membership (CAT0 + after 1980-01-01)
        CAT1 = CAT0 and row.time >= DIGITAL_STARTDATE
        catalog_status['CAT1'].append(CAT1)
        # CAT2 Membership (CAT1 + post 2001-01-01)
        CAT2 = CAT1 and row.time >= CONTINUOUS_STARTDATE
        catalog_status['CAT2'].append(CAT2)
        # Get preferred origin object
        cat = EBANK.get_events(event_id=event_id)
        prefor = cat[0].preferred_origin()
        # Evaluate well constrained status
        # Fixed hypocentral parameter flags
        catalog_status['time_fixed'].append(prefor.time_fixed)
        catalog_status['epi_fixed'].append(prefor.epicenter_fixed)
        catalog_status['depth_fixed'].append(prefor.depth_type != 'from location')
        # Horizontal uncertainties too large?
        if row.horizontal_uncertainty > max_loc_err_m:
            herr_large = True
        elif not np.isfinite(row.horizontal_uncertainty):
            herr_large = True
        else:
            herr_large = False
        catalog_status['herr_large'].append(herr_large)
        # Vertical uncertainties too large?
        if row.vertical_uncertainty > max_loc_err_m:
            zerr_large = True
        elif not np.isfinite(row.vertical_uncertainty):
            zerr_large = True
        else:
            zerr_large = False
        catalog_status['zerr_large'].append(zerr_large)
        # Time uncertainties too large?
        if prefor.time_errors.uncertainty is None:
            terr_large = False
        elif prefor.time_errors.uncertainty > max_rms_s:
            terr_large = True
        else:
            terr_large = False
        catalog_status['terr_large'].append(terr_large)
        # Observation count
        if row.used_phase_count < min_obs:
            low_phase_count = True
        else:
            low_phase_count = False
        catalog_status['low_phase_count'].append(low_phase_count)
        # max distance to closest observing station
        if prefor.quality is None:
            closest_sta_too_far = True
        elif prefor.quality.minimum_distance is None:
            closest_sta_too_far = True
        elif not np.isfinite(prefor.quality.minimum_distance):
            breakpoint()
            closest_sta_too_far = True
        else:
            closest_sta_too_far = prefor.quality.minimum_distance > max_closest_m
        catalog_status['closest_sta_too_far'].append(closest_sta_too_far)
        # minimum number of unique observing stations
        if len(prefor.arrivals) < min_sta:
            too_few_sta = True
        else:
            unique = set()
            for arr in prefor.arrivals:
                pick = arr.pick_id.get_referred_object()
                nslc = pick.waveform_id.id
                n,s,l,c = nslc.split('.')
                unique.add('.'.join([n,s]))
            too_few_sta = len(unique) < min_sta
        catalog_status['too_few_sta'].append(too_few_sta)

        # Finally evaluate Well Constrained Status
        # Initially set to true, then check if anything fails (a True result for individual tests)
        wc_status = True
        # Iterate across keys
        for _f in catalog_status.keys():
            # Ignore non evaluation keys
            if _f in ['CAT0','CAT1','CAT2','event_id','prefor_time','offset_km','etype']:
                continue
                # ['time_fixed','epi_fixed','depth_fixed','herr_large','zerr_large','terr_large','low_phase_count','closest_sta_too_far','too_few_sta']:
            # If the last entry (this iteration) returns a True (True means fails test)
            elif catalog_status[_f][-1]:
                # set Well Constrained Status as False
                wc_status = False
                # Break iteration, one failure means total failure
                break
        # Append result
        catalog_status['WC'].append(wc_status)

    # Format for output
    evid = catalog_status.pop('event_id')
    df_out = pd.DataFrame(catalog_status, index=evid)
    df_out.index.name = 'event_id'
    df_out.to_csv(SAVENAME, header=True, index=True)
    return df_out


if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.DEBUG)
    Logger.addHandler(CriticalExitHandler(exit_code=1))
    # Make sure savepath exists
    SAVEPATH = os.path.split(SAVENAME)[0]
    try:
        os.makedirs(SAVEPATH)
    except FileExistsError:
        pass
    # RUN MAIN, RETURN df_out in case this is run inside (i)python
    df_out = main()


 