import os, logging
from pathlib import Path

import pandas as pd
from obspy import read_inventory
from obsplus import EventBank, WaveBank

from eqcutil.core.load_from_bank import generate_clustering_tribe_from_banks
from eqcutil.augment.template import rename_templates
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler, rich_error_message

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
# base_path's for WaveBank and EventBank 
WBBP = ROOT / 'data' / 'WF' / 'BANK'
EBBP = ROOT / 'data' / 'XML' / 'QUAKE' / 'BANK'
# Path to templating directory
TMPD = ROOT / 'processed_data' / 'bank_templates' / 'd20km'
# templating STATUS file
TMPF = TMPD / 'templating_status.csv'
# Path to INV
INVF = ROOT / 'data' / 'XML' / 'INV' / 'station_inventory_50km_MBS_IRISWS_STATION.xml'
# Save Path
SAVE = TMPD

### FILTERING PARAMETER DEFINITIONS ###
# Channel Mapping Kwargs
TMAP = {'UW.MBW.01.EHZ':'UW.MBW..EHZ',
        'UW.MBW2..HHE':'UW.MBW2..HHZ',
        'UW.MBW2..ENZ':'UW.MBW2..HHZ',
        'UW.MULN..HHN':'UW.MULN..HHZ',
        'UW.RPW.01.EHZ':'UW.RPW..EHZ',
        'UW.SHUK..HHN':'UW.SHUK..HHZ',
        'UW.PASS..BHN':'UW.PASS..BHZ'}


STA_LISTS = [['MBW','MBW2'],['RPW','RPW2'],['SHUK'],['SAXON'],['PASS'],['MULN']]
for slist in STA_LISTS:
    # Filter Picks Kwargs
    FPKW = {'stations': slist,
            'phase_hints': ['P'],
            'enforce_single_pick': 'preferred'}
    # Template Build Status Filters
    TBS = ['Full','Partial']

    # Clustering Kwargs
    ccckwargs = {'method': 'correlation_cluster',
                'replace_nan_distances_with': 'mean',
                'shift_len': 10,
                'corr_thresh': 0.6,
                'allow_individual_trace_shifts': True,
                'show': False,
                'cores': 'all',
                'save_corrmat': False}
    # Save Clustering parameters
    with open(str(SAVE / 'clustering_params.csv'), 'w') as PAR:
        PAR.write('param,value\n')
        for _k, _v in ccckwargs.items():
            PAR.write(f'{_k},{_v}\n')    

    ### CONNECT / LOAD ###
    # Load template STATUS file
    df_status = pd.read_csv(TMPF)
    # Filter tempalate status file by template build status
    df_status = df_status[df_status.build_summary.isin(TBS)]
    # FIXME: Shorten list for debugging
    # df_status = df_status.iloc[20:100,:]

    # Get name(s) of unique param files
    param_files = df_status.param_file_name.unique()

    if len(param_files) != 1:
        Logger.critical('Found more than one parameter file in filtered STATUS')
    else:
        df_pf = pd.read_csv(str(TMPD/'param'/param_files[0]))
        # Load raw
        params = {_r.param: _r.value for _, _r in df_pf.iterrows()}
        # Parse types
        for _k, _v in params.items():
            if _v in ['True','False']:
                params[_k] = bool(_v)
            elif '.' in _v:
                params[_k] = float(_v)
            elif _v.isnumeric():
                params[_k] = float(_v)

    # Connect To Banks
    EBANK = EventBank(EBBP)
    WBANK = WaveBank(WBBP)
    # Load Inventory
    INV = read_inventory(INVF)

    ### PROCESSING ###
    # Generate Templates from Banked Data/Metadata
    ctr = generate_clustering_tribe_from_banks(
        WBANK, EBANK,
        df_status.event_id,
        transfer_mapping=TMAP,
        pick_filt_kwargs=FPKW,
        creation_kwargs=params)
    Logger.info('ClusteringTribe generated from Wave-/Event-Bank')
    # Rename Templates
    ctr = rename_templates(ctr)
    Logger.info('Templates renamed')

    # Run Cross-Correlation Clustering
    Logger.info(f'Running clustering for {slist[0]}')
    ctr.cluster(**ccckwargs)
    Logger.info(f'saving to disk')
    ctr.write(str(SAVE/f'corr_cluster_{slist[0]}_FreeStep.tgz'))
