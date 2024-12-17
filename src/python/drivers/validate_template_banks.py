"""
:module: Mt_Baker_LF_Research/src/python/drivers/validate_template_banks.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This driver script uses data and metadata from EventBank and
    WaveBank construction to determine if templates can be generated
    from on-disk data.

    The script scans over the wavebank commit log to find which events
    have any waveform data associated with them and then  attempts to
    generate :class:`~eqcorrscan.Template` objects for each event using
    preferred origin picks. The parameters for generating the templates
    are written to a PARAM file and metadata from the attempted template 
    generation are written to a STATUS file.
"""

import os, logging
from pathlib import Path

import pandas as pd
from obspy import Catalog, read_inventory
from eqcorrscan import Template, Tribe
from obsplus import EventBank, WaveBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path's for WaveBank and EventBank 
WBBP = ROOT / 'data' / 'WF' / 'BANK'
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Path to Event Subset
C20K = ROOT/'data'/'Events'/'MtBaker_20km_radius_phases.csv'
# Path to INV
INVF = ROOT / 'data' / 'XML' / 'INV' / 'station_inventory_50km_MBS_IRISWS_STATION.xml'
# Path to WaveBank processing log file
WBLF = WBBP.parent / 'processing_log.csv'
# Save Path
SAVE = ROOT / 'processed_data' / 'bank_templates' / 'd20km'

####################################
### PROCESSING PARAMETER SECTION ###

# Define Catalog/Pick Subsetting Kwargs
PHZ = ['P','S']
STA = ['MBW','MBW2','SHUK']

# vel_model = 'P4'


# Define template construction kwargs
PARFILENAME='shortened_params.csv'
STATUSFILE=SAVE/'templating_status.csv'
overwrite_protect=True
ckwargs = {'method': 'from_client',
           'lowcut': 0.5,
           'highcut': 20.,
           'filt_order': 4,
           'samp_rate': 50.,
           'prepick': 5.,
           'length': 50.,
           'process_len': 300.,
           'min_snr': 3.,
           'parallel': True,
           'num_cores': 6, 
           'save_progress': False
           }

##########################
### PROCESSING SECTION ###

### Directory Creation
if not os.path.exists(SAVE):
    Logger.info(f'Creating save directory: {str(SAVE)}')
    os.makedirs(SAVE)

if not os.path.exists(SAVE/'param'):
    os.makedirs(SAVE/'param')

# Save template parameters to parameter file (sans WaveBank object)
PFNAME = SAVE/'param'/PARFILENAME
with open(str(PFNAME), 'w') as PAR:
    PAR.write('param,value\n')
    for _k, _v in ckwargs.items():
        PAR.write(f'{_k},{_v}\n')

# Initialize log file header if the file doesn't exist

if not os.path.isfile(STATUSFILE):
    overwrite_protect = False
else:
    if not overwrite_protect:
        overwrite_protect = False

if not overwrite_protect:
    with open(str(STATUSFILE), 'w') as LOG:
        LOG.write('event_id,build_summary,template_generates,all_picks_present,unpicked_stations,param_file_name,error_type\n')

### CONNECTION / LOADING ###
# Connect to EventBank
EBANK = EventBank(base_path=str(EBBP))
# Assert that the index exits
if len(EBANK.read_index()) == 0:
    Logger.critical('EventBank shows zero entries')
else:
    Logger.info('Connected to EventBank')

# Connect to WaveBank
Logger.info(f'Initializing/connecting to WaveBank at {WBBP}')
WBANK = WaveBank(base_path=str(WBBP))
# Assert that the index exits
if len(WBANK.read_index()) == 0:
    Logger.critical('WaveBank shows zero entries')
else:
    Logger.info("Connected to WaveBank")
    # Update Template Construction Kwargs with WaveBank client
    ckwargs.update({'client_id': WBANK})

# Get event_id set that loaded into WaveBank
Logger.info('Loading WaveBank submission log file')
df_WBLF = pd.read_csv(WBLF)
# Subset by successful commits to WaveBank
df_WBLF = df_WBLF[df_WBLF.commit_status]
# Get event_ids
REF_EVIDS = df_WBLF.event_id.unique()
Logger.info(f'WaveBank has {len(REF_EVIDS)} unique EVIDs from successful loads')

# # Load Inventory
INV = read_inventory(INVF)

# Get event subset
cat = EBANK.get_events(event_id = REF_EVIDS)
# Apply phase hints
cat = apply_phase_hints(cat)
# Filter to exclued multiple picks
cat = filter_picks(cat, evaluation_mode='all', enforce_single_pick='preferred')
# Attempt to make templates one event at a time
pfnameshort = os.path.split(PFNAME)[-1]
for _e, event in enumerate(cat):
    # Re-open Status Log file with append
    with open(str(STATUSFILE), '+a') as LOG:
        # Get EVID
        evid = event.resource_id.id
        Logger.info(f'Processing {evid} ({_e+1}/{len(cat)})')
        # Get active stations
        inv = INV.select(time=event.preferred_origin().time)
        Logger.info(f'Active Station Count: {len(inv.get_contents()["stations"])}')
        Logger.info(f'Pick Count: {len(event.picks)} | Assoc Count: {len(event.preferred_origin().arrivals)}')
        # Update event in construction kwargs
        ckwargs.update({'catalog': Catalog(events=[event])})
        # Try to run template creation
        try:
            tribe = Tribe().construct(**ckwargs)
            bld_status = True
            err_status = False
        # If this raises an error, capture the error in the output
        except Exception as e:
            Logger.warning(rich_error_message(e))
            sum_status = 'Error'
            err_status = type(e).__name__
            bld_status = False
            app_status = False
            ups_status = True

        # If tribe is built
        if bld_status:
            # Check if templates is empty
            if len(tribe.templates) == 0:
                Logger.warning(f'{evid} produced an empty Tribe')
                sum_status = 'Empty'
                app_status = False
                ups_status = True
            else:
                # Set assumed values
                app_status = True
                sum_status = 'Full'
                ups_status = False
                # Verify if all picks are present
                for pick in event.picks:
                    if len(tribe.templates[0].st.select(id=pick.waveform_id.id)) == 0:
                        app_status = False
                        sum_status='Partial'
                        break
                # Verify if all active stations are present
                for sta in inv.get_contents()['stations']:
                    if len(tribe.templates[0].st.select(id=f'{sta}*')) == 0:
                        ups_status = True
                        break
        # Write results to log
        LOG.write(f'{evid},{sum_status},{bld_status},{app_status},{ups_status},{pfnameshort},{err_status}\n')
        # Make sure tribe doesn't bleed over into other iterations
        del tribe


# Signal End Of Process
Logger.critical('SUCCESSFULLY COMPLETED PROCESS - EXITING')