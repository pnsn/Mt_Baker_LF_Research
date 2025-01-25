import os, logging, glob
from pathlib import Path

from obspy import read_inventory
from obspy.core.event import WaveformStreamID
from obspy.clients.fdsn import Client
import pandas as pd
from obsplus import EventBank, WaveBank
from eqcorrscan import Tribe

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.augment.template import rename_templates, augment_template
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

Logger = setup_terminal_logger(os.path.split(__file__)[-1])
Logger.addHandler(CriticalExitHandler())
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to wavebank
WBBP = ROOT / "data" / "WF" / "BANK"
# Get RESP file directory
INVD = ROOT / 'data' / 'XML' / 'RESP'
# Path to event metadata from AQMS
AQMS = ROOT / 'data' / 'Events' / 'MtBaker_50km_radius_origins.csv'
# Get result from step1
S1INFO = ROOT / 'processed_data' / 'workflow' / 'catalog' / 'preferred_stachan_pick_continuity.csv'

# Output directory absolute path
DOUT = ROOT / "processed_data" / "workflow" / "templates"
# Status file save
STATUSFILE = DOUT / "template_construction_status.csv"
# Parameter file for filter_pick (Pick Filtering) parameter documentation
PFNAME = DOUT / "catalog_filtering_kwargs.csv"
# Parameter file for construction (Template Construction) parameter documentation
TCNAME = DOUT / "template_construction_kwargs.csv"
# Name to save the constructed ClusteringTribe to
CTRIBE = DOUT / "clustering_tribe_v3.tgz"

### LOAD SUBSET EVENTS and CHANNELS
df_SE = pd.read_csv(S1INFO, index_col=[0])
# Get unique station names
STAS = {col.split('.')[1] for col in df_SE.columns}

# Channel Shifts - Move picks from horizontals so that
# templates will only use vertical data.
# Vertical traces still have to pass the SNR minimum!
SHIFTS = {'UW.SHUK..HHN':'UW.SHUK..HHZ',
          'UW.MBW2..HHE':'UW.MBW2..HHZ',
          'UW.MBW2..ENZ':'UW.MBW2..HHZ',
          'UW.MULN..HHN':'UW.MULN..HHZ',
          'CN.VDB..SHZ':'CN.VDB..EHZ'}
# Location codes to remove in order to permit records at the same
# site from different instruments to be cross-correlated
ALIASES = {'UW.MBW.01.EHZ':'UW.MBW..EHZ',
           'UW.RPW.01.EHZ':'UW.RPW..EHZ'}


## filter_picks key-word arguments
pick_filt_kwargs = {'phase_hints': ['P'],
                    'enforce_single_pick': 'preferred',
                    'stations': STAS}

## Construct Tribe key-word arguments
ckwargs = {'method': 'from_client',
           'lowcut': 0.5,
           'highcut': 20.,
           'filt_order': 4,
           'samp_rate': 50.,
           'prepick': 5.,
           'length': 45.,
           'process_len': 300.,
           'min_snr': 3.,
           'parallel': True,
           'num_cores': 12, 
           'save_progress': False
           }

### PROCESSING ###
Logger.info(f'Will save to {str(DOUT)}')
if not os.path.exists(str(DOUT)):
    os.makedirs(str(DOUT))

### GET INVENTORY ###
flist = glob.glob(str(INVD/'*.xml'))
for _e, _f in enumerate(flist):
    if _e == 0:
        INV = read_inventory(_f)
    else:
        INV += read_inventory(_f)

### SAVE PICK FILTERING PARAMETERS
with open(str(PFNAME), 'w') as PAR:
    PAR.write('param,value\n')
    for _k, _v in ckwargs.items():
        PAR.write(f'{_k},"{_v}"')

### SAVE TEMPLATE CONSTRUCTION PARAMETERS
with open(str(TCNAME), 'w') as PAR:
    PAR.write('param,value\n')
    for _k, _v in ckwargs.items():
        PAR.write(f'{_k},{_v}\n')

# Connect to eventbank
Logger.info('Connecting to eventbank')
EBANK = EventBank(EBBP)
df_eb = EBANK.read_index()
if len(df_eb) == 0:
    Logger.critical(f'Empty EventBank - check path: {str(EBBP)}')
# Create an EVID column 
df_eb = df_eb.assign(evid=[int(os.path.split(row.event_id)[-1]) for _, row in df_eb.iterrows()])
df_eb = df_eb.sort_values(by='time')

# Connect to wavebank
Logger.info('Connecting to wavebank')
WBANK = WaveBank(WBBP)

# Update wavebank as client for template construction
ckwargs.update({'client_id': WBANK})

# Initialize ClusteringTribe object
CTR = ClusteringTribe()

# Write STATUSFILE HEADER
with open(str(STATUSFILE), 'w') as LOG:
    LOG.write('event_id,build_summary,template_generates,all_picks_present,unpicked_stations,error_type,pick_filter_par_name,construct_par_name\n')

# Iterate across events
for _e, event_id in enumerate(df_SE.index):
    # Re-open statusfile to append the status line for this event
    with open(str(STATUSFILE), '+a') as LOG:
        Logger.info(f'Processing: {event_id} ({_e+1} of {len(df_SE)})')
        # Get event
        cat = EBANK.get_events(event_id=event_id)
        event = cat[0].copy()
        for event in cat.events:
            for pick in event.picks:
                if pick.waveform_id.id in SHIFTS.keys():
                    pick.waveform_id = WaveformStreamID(seed_string=SHIFTS[pick.waveform_id.id])
        # Apply phase hints to picks
        cat = apply_phase_hints(cat)
        # Filter picks
        cat = filter_picks(cat, **pick_filt_kwargs)
        # update ckwargs with catalog
        ckwargs.update({'catalog':cat})
        # Get preferred origin time
        prefor = event.preferred_origin()
        tO = prefor.time
        # Get current inventory
        inv = INV.select(time=tO)
        # Try to construct template
        try:
            tribe = Tribe().construct(**ckwargs)
            # Flag as successful build
            bld_status = True
            # Flag as no error status
            err_status = False
            tribe = rename_templates(tribe)
        # If construction fails
        except Exception as e:
            # Send error to the command line log
            Logger.warning(rich_error_message(e))
            # Update status fields for this line
            sum_status = 'Error'
            # Catch error type for statusfile
            err_status = type(e).__name__
            # Flag Build as a failure
            bld_status = False
            # Flag All Picks Present false
            app_status = False
            # Flag unpicked stations as true (a little counter-intuitive...?)
            ups_status = True

        # If tribe is built
        if bld_status:
            # Check if templates is empty
            if len(tribe.templates) == 0:
                Logger.warning(f'{event_id} produced an empty Tribe')
                sum_status = 'Empty'
                app_status = False
                ups_status = True
            # If not empty
            else:
                Logger.info(f'ClusteringTribe has {len(CTR)} members ({_e+1} attempted)')
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
                    else:
                        continue
                # Apply any aliases to trace ID's
                for tmp in tribe:
                    for tr in tmp.st:
                        if tr.id in ALIASES.keys():
                            tr.stats.location = ALIASES[tr.id].split('.')[-2]
                    # Merge traces
                    tmp.st.merge(method=1, interpolation_samples=-1)
                # Append tribe to ClusteringTribe  
                CTR += tribe.copy()
        else:
            pass
        # Write results to log
        LOG.write(f'{event_id},{sum_status},{bld_status},{app_status},{ups_status},{err_status},{str(PFNAME)},{str(TCNAME)}\n')
        # Make sure tribe doesn't bleed over into other iterations
        try:
            del tribe
        except NameError:
            continue

# Append EBANK summary and ETYPE to CTR.clusters
df_eb.index = df_eb.evid.apply(lambda x: f'uw{x}')
# Load Event Metadata
df_O = pd.read_csv(str(AQMS))
# Subset to evid and etype
df_O = df_O[['evid','etype']]
# Drop duplicates
df_O = df_O.drop_duplicates(keep='first')
# Generate matching index to CTR.clusters
df_O.index = [f'uw{_e}' for _e in df_O.evid]

CTR.clusters = CTR._c.join(df_eb, how='left')
CTR.clusters = CTR._c.join(df_O['etype'], how='left')
CTR.write(str(DOUT/'clustering_tribe.tgz'))