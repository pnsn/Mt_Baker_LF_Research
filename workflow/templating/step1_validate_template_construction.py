import os, logging, glob
from pathlib import Path

from obspy import read_inventory
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
INVD = ROOT / 'XML' / 'RESP'
# Get absolute path to phase file
PHZF = ROOT / "data" / "Events" / "MtBaker_20km_radius_phases.csv"
# Output directory absolute path
DOUT = ROOT / "processed_data" / "workflow" / "templates"
# Status file save
STATUSFILE = DOUT / "template_construction_status.csv"
PFNAME = DOUT / "catalog_filtering_kwargs.csv"
# Parameter file to save
TCNAME = DOUT / "template_construction_kwargs.csv"

# Preferred Stations
STAS = ['MBW','MBW2','SHUK','RPW','RPW2','JCW','SAXON']

pick_filt_kwargs = {'phase_hints': ['P'],
                    'enforce_single_pick': 'preferred',
                    'stations': STAS}

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

### PROCESSING ###
Logger.info(f'Will save to {str(DOUT)}')
if not os.path.exists(str(DOUT)):
    os.makedirs(str(DOUT))


### GET INVENTORY ###
flist = glob.glob(str(INVD/'*.xml'))

### SAVE PARAMETERS
with open(str(TCNAME), 'w') as PAR:
    PAR.write('param,value\n')
    for _k, _v in ckwargs.items():
        PAR.write(f'{_k},{_v}\n')

with open(str(PFNAME), 'w') as PAR:
    PAR.write('param,value\n')
    for _k, _v in ckwargs.items():
        PAR.write(f'{_k},"{_v}"')

# Connect to eventbank
Logger.info('Connecting to eventbank')
EBANK = EventBank(EBBP)
df_eb = EBANK.read_index()
if len(df_eb) == 0:
    Logger.critical(f'Empty EventBank - check path: {str(EBBP)}')
df_eb = df_eb.assign(evid=[int(os.path.split(row.event_id)[-1]) for _, row in df_eb.iterrows()])

# Connect to wavebank
Logger.info('Connecting to wavebank')
WBANK = WaveBank(WBBP)

ckwargs.update({'client_id': WBANK})

# Get phases
Logger.info("Loading phase file")
df_phz = pd.read_csv(str(PHZF))
df_phz_pf = df_phz[df_phz.sta.isin(STAS)]

# Get evids
if len(df_phz.evid.unique()) != len(df_phz_pf.evid.unique()):
    Logger.info(f'Preferred station list returns {len(df_phz_pf.evid.unique())} events')

# Subset preferred events
df_eb_pref = df_eb[df_eb.evid.isin(df_phz_pf.evid.unique())]

# Get event catalog
CTR = ClusteringTribe()

with open(str(STATUSFILE), 'w') as LOG:
    LOG.write('event_id,build_summary,template_generates,all_picks_present,unpicked_stations,param_file_name,error_type\n')

for _e, event_id in enumerate(df_eb_pref.event_id):
    with open(str(STATUSFILE), '+a') as LOG:
        Logger.info(f'Processing: {event_id} ({_e+1} of {len(df_eb_pref)})')
        cat = EBANK.get_events(event_id=event_id)
        cat = rename_templates(cat)
        cat = apply_phase_hints(cat)
        cat = filter_picks(cat, **pick_filt_kwargs)
        event = cat[0]
        try:
            tribe = Tribe().construct(**ckwargs)
            bld_status = True
            err_status = False
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
                    Logger.warning(f'{event_id} produced an empty Tribe')
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
            LOG.write(f'{event_id},{sum_status},{bld_status},{app_status},{ups_status},{pfnameshort},{err_status}\n')
            # Make sure tribe doesn't bleed over into other iterations
            del tribe  