"""
:program: Mt_Baker_LF_Research/workflow/templating/step2_augment_catalog.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose:
    This program iterates across events and stations selected in "step1_station_selection.py" 
    and does the following filtering/augmentations:
        1) Sub-sample event metadata to include only selected stations and phase types
        2) Models the arrival times of desired phase type(s) at stations reportedly active
            during a catalog event (i.e., station metadata says it was running)
        3) Appends the PNSN AQMS event-type classification as a 2-character comment on
            each obspy Event object's **comments** attribute
    
    filtered/augmented event metadata are then saved to a new "Augmented Event Bank" that
    is used as the input for template generation in "step3_generate_templates.py"
"""

import os, logging, glob, warnings
from pathlib import Path

import pandas as pd
from tqdm import tqdm

from obspy import read_inventory, Inventory, UTCDateTime
from obspy.core.event import Comment

# from obspy.clients.fdsn import Client
from obsplus import EventBank #, WaveBank

from eqcutil.augment.catalog import apply_phase_hints, filter_picks
from eqcutil.catalog.model_phases import model_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler
from eqcutil.viz.eqc_compat import plant

plant()

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))
# Ignore FutureWarning
warnings.simplefilter(action='ignore', category=FutureWarning)

### USER CONTROLS ###
# Modeled Pick Processing Parameters
phases=['P','p']            # classic phase labels for Pyrocko/Cake to model
velocity_model = 'P4'       # PNSN velocity model name to use for ray-tracing
MODCHAN = '??Z'             # Channel code to use for modeled picks

# Prioritization of instrument / band codes for assigning model picks
# I.e., if there are multiple active sensors at a station being modeled, which should be
# prioritized for retrieving trace data and auto-picking
ITYPE_PRIORITY_ORDER = ['HH','BH','EH','HN','EN'] 

# Developmental - auto_picker
# auto_pick_pad = [300, 300]    # [sec] front and back padding around modeled pick for autopicking
# RMS_win = 3.                # [sec] window on either side of a pick for calculating RMS amplitudes
# min_snr = 1.1               # [dB] Minimum acceptable SNR for final pick
# min_noise_RMS = 3           # [dB] Minimum acceptable SNR for noise segment in front of final modeled pick


# Manual Pick Processing Parameters
pick_filt_kwargs = {'phase_hints': phases,
                    'enforce_single_pick': 'preferred'}
# Mappings for transferring picks to other channels on the same instrument
# Syntax: {FROM: TO}
SHIFTS = {'N':'Z',
          'E':'Z'}

# Figure rendering parameters
FMT = 'PNG'
DPI = 120
plot_save_kwargs = {}

### END OF USER CONTROLS ###

### PATH DEFINITIONS ###

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to station response files
INVD = ROOT / 'data' / 'XML' / 'RESP'
# Output Path
OUTD = ROOT / 'processed_data' / 'workflow' / 'catalog'
# Preferred Station Picked Events File (save file)
PSPEF = OUTD / 'preferred_event_sta_picks.csv'
# Event Type Metadata
ETDATA = ROOT / 'data' / 'Events' / 'Mount_Baker_evids_etypes_10_JAN_2025.csv'
# WaveBank Base Path
WBBP = ROOT / "data" / "WF" / "BANK"


### PROCESSING SECTION ###
# Load preferred_event_sta_picks
df_evids = pd.read_csv(str(PSPEF), index_col=[0])
# Filter by phases
df_evids = df_evids[df_evids.phase.isin(phases)]
df_evids = df_evids.assign(evid=[int(_e.split('/')[-1]) for _e in df_evids.index])
# Load etype metadata
df_meta = pd.read_csv(ETDATA, index_col=[0])
# Join dataframes
df_evids = df_evids.join(df_meta, on='evid', how='left')
# Get preferred stations
STAS = list(df_evids.sta.unique())
# Get preferred event_id's
EVIDS = list(df_evids.index.unique())
EVIDS.sort()
pick_filt_kwargs.update({'stations': STAS})

# INITIALIZE CLIENTS
EBANK = EventBank(EBBP)
# Ensure that augmented event bank has a home
if not os.path.exists(str(OUTD/'AUGMENTED_BANK')):
    os.makedirs(str(OUTD/'AUGMENTED_BANK'))

# Initialize Augmented event bank (we're writing to this bank)
ABANK = EventBank(base_path=OUTD/'AUGMENTED_BANK',
                    path_structure='{year}',
                    name_structure='uw{event_id_end}')

# WBANK = WaveBank(base_path=WBBP)

# IRIS = Client("IRIS")

# Load inventory for preferred stations
INV = Inventory()
for _f in glob.glob(str(INVD/'*.xml')):
    # Read the provided stationXML
    inv = read_inventory(_f)
    # Filter for desired stations
    for _sta in STAS:
        INV += inv.select(station=_sta,channel='[BHE][NH][ZNE]')


df_seb = EBANK.read_index()
df_aeb = ABANK.read_index()
# Iterate across event_ids
if Logger.level < 30:
    tqdm_disable = True
else:
    tqdm_disable = False

pick_info = []
for event_id in tqdm(EVIDS, disable=tqdm_disable):
    # Confirm that event_id is not already in AUGMENTED
    if event_id in df_aeb.event_id.values:
        Logger.info(f'Event {event_id} already in Augmented Event Bank - Skipping')
        continue

    Logger.info(f'Processing {event_id}')
    # Get event
    cat = EBANK.get_events(event_id=[event_id])
    # Ensure event is non-empty
    if len(cat) != 1:
        Logger.critical('fetched non-singleton event')
    else:
        event=cat[0]

    ### FILTER EXISTING ANALYST DATA ###
    # Ensure that phase hints are applied
    cat = apply_phase_hints(cat)
    # Filter event
    cat = filter_picks(cat, **pick_filt_kwargs)
    # Ensure event is non-empty
    if len(cat) == 0:
        Logger.critical('filtering resulted in an empty catalog')
    
    
    event=cat[0]
    ### ADD EVENT TYPE IN EVENT COMMENTS ###
    etype = df_evids.loc[event_id,'etype']
    if not isinstance(etype, str):
        etype = etype[0]
    if len(etype) != 2:
        breakpoint()
    event.comments.append(Comment(text=etype))
    # Get metadata from retained picks & apply channel shifts
    picksta = set([])
    kprid = []
    Logger.debug(f'Filtering preserved {len(event.picks)} picks')
    for _p in event.picks:
        kprid.append(_p.resource_id)
        picksta.add(_p.waveform_id.station_code)
        ccode = _p.waveform_id.channel_code
        if ccode[-1] in SHIFTS.keys():
            scode = ccode[:-1] + SHIFTS[ccode[-1]]
            comment = Comment(text=f'Originally on channel {ccode}')
            Logger.warning(f'Shifting pick on {_p.waveform_id.id} to component {scode}')
            _p.waveform_id.channel_code = scode
            _p.comments.append(comment)
            
    # Station magnitudes cleanup
    smk = []
    for smag in event.station_magnitudes:
        if smag.waveform_id.station_code in STAS:
            smk.append(smag)
    Logger.debug(f'Removing {len(event.station_magnitudes) - len(smk)} unreferenced station_magnitudes')

    event.station_magnitudes = smk

    # Amplitudes cleanup
    amk = []
    for amp in event.amplitudes:
        if amp.waveform_id.station_code in STAS:
            amk.append(amp)
    Logger.debug(f'Removing {len(event.amplitudes) - len(amk)} unreferenced amplitudes')
    event.amplitudes = amk
    prefor = event.preferred_origin()

    # Arrivals cleanup
    ark = []
    for arr in prefor.arrivals:
        if arr.pick_id in kprid:
            ark.append(arr)
    Logger.debug(f'Removing {len(prefor.arrivals) - len(ark)} unreferenced arrivals')
    prefor.arrivals = ark

    ### MODEL ARRIVAL TIMES ###

    # Get active station list
    inv = INV.select(time=prefor.time, channel=MODCHAN)
    active = inv.get_contents()['channels']
    ASTAS = set([_c.split('.')[1] for _c in active])
    modsta = list(set(ASTAS).difference(picksta))
    # Model arrivals for stations not picked
    msg = 'Modelling arrival times for station(s):'
    for ms in modsta:
        msg += f' {ms},'
    msg = msg[:-1]
    Logger.debug(msg)

    # Model by individual station using P4 velocity model
    modeled_first_arrivals = []
    for sta in modsta:
        sinv = inv.select(station=sta)
        # Apply instrument type priority order
        for itype in ITYPE_PRIORITY_ORDER:
            sbinv = sinv.select(channel=f'{itype}?')
            if len(sbinv.get_contents()['channels']) == 1:
                break
        # If there is more than one channel in this inventory - go to debug
        if len(sbinv.get_contents()['channels']) != 1:
            breakpoint()

        # Model pick times using the P4 velocity model
        picks_hat = model_picks(
            prefor, sbinv, 
            model_name=velocity_model,
            phases=phases)
        
        # Identify earliest modeled arrival for each channel
        t0 = UTCDateTime()
        # Iterate across modeled picks
        for _ph in picks_hat:
            if _ph.time < t0:
                t0 = _ph.time
                earliest = _ph

        # # Run auto-picker approach for "earliest"
        # request = earliest.waveform_id.id.split('.') + [earliest.time - auto_pick_pad[0], earliest.time + auto_pick_pad[1]]
        # # Try to get waveforms from local client first
        # try:
        #     st = WBANK.get_waveforms(*request)
        #     Logger.debug('successful WaveBank request')
        # except:
        #     Logger.warning(f'WaveBank request failed for {request} - trying IRIS')
        #     try:
        #         st = IRIS.get_waveforms(*request)
        #     except:
        #         Logger.warning('IRIS request failed - abandoning modeled pick')
        #         continue

        # # Run picker
        # breakpoint()

         # If there were any modeled arrivals, append the earliest
        if len(picks_hat) > 0:
            modeled_first_arrivals.append(earliest)
        # If there were no modeled arrivals, proceed to next station
        # for modeling picks
        else:
            continue

    Logger.info(f'Adding {len(modeled_first_arrivals)} modeled first arrivals to {event_id}')
    event.picks += modeled_first_arrivals
    for _p in event.picks:
        nslc = _p.waveform_id.id
        line = [event_id, _p.phase_hint, nslc, _p.evaluation_mode] + nslc.split('.')
        pick_info.append(line)
    # Try committing updated event to AUGMENTED EVENT BANK
    try:
        ABANK.put_events(cat)
        Logger.info('Successfully added to Augmented Event Bank')
    except:
        breakpoint()



df_picked = pd.DataFrame(
    pick_info,
    columns=['event_id','phase','nslc','eval_mode','net','sta','loc','chan']
)
# # Get count of station-event picks
# ser_top_picked = df_picked.sta.value_counts()
# df_p_cont = df_picked[['event_id','sta']]
# breakpoint()
# df_p_cont = df_p_cont.pivot_table(
#     index=['event_id'],
#     columns=['sta'],
#     aggfunc=lambda x: 1, fill_value=0)
# for evid, row in df_p_cont.iterrows():
#     if 
                         

# Show new pick availability 

    # # Cleanup non-referred items
    # cat[0].amplitudes=[]
    # cat[0].station_magnitudes=[]
    # breakpoint()
    # # for event in cat:
    # #     for 
    # # Get origin time
    # prefor = cat[0].preferred_origin()
    # otime = prefor.time
    # # Get subset inventory
    # inv = INV.select(time=otime)

#     # Iterate across event_picks
#     for event_id, row in df_evids.iterrows():
#         Logger.info(f'Processing {row.uw}')
#         # Iterate across stations
#         for sta in 
    

#     # Iterate across events
#     for evid, meta in df_evids.iterrows():
#         Logger.info(f'Processing {evid} ({meta.etype})')
#         # Fetch event from event bank
#         cat = EBANK.get_events(event_id=evid)
#         # Apply phase hints
#         cat = apply_phase_hints(cat)
#         # Subset to selected stations and phases only
#         cat = filter_picks(cat, **pick_filt_kwargs)
#         # Create a comment for etype
#         comment = Comment(text=meta.etype)
#         cat[0].comments.append(comment)
#         # # Apply channel shifts
#         # for event in cat:
#         #     for pick in event.picks:
#         #         if pick.waveform_id.id in SHIFTS.keys():
#         #             pick.waveform_id = WaveformStreamID(
#         #                 seed_string=SHIFTS[pick.waveform_id.id]
#         #             )
#         # Get active channels
#         active = INV.select(time=cat[0].preferred_origin().time)
#         # Get picked channels
#         picked = []
#         kept_arrivals = []
#         for _arr in cat[0].preferred_origin().arrivals:
#             pick = _arr.pick_id.get_referred_object()
#             if not pick is None:
#                 picked.append(pick.waveform_id.id)
#                 kept_arrivals.append(_arr)
#         if kept_arrivals == []:
#             breakpoint()
#         cat[0].preferred_origin().arrivals = kept_arrivals

        
#         # Get unpicked channels
#         unpicked = set(active.get_contents()['channels']).difference(set(picked))
#         Logger.info(f'{len(picked)} picked channels')
#         Logger.info(f'{len(unpicked)} unpicked channels')
#         Logger.info(f'Modeling for phases {phases}')
#         # Model arrival times for unpicked channels
#         for nslc in unpicked:
#             n, s, l, c = nslc.split('.')
#             iactive = active.select(network=n,station=s, location=l, channel=c)
#             picks_hat = model_picks(cat[0].preferred_origin(),
#                                     iactive,
#                                     model_name=velocity_model,
#                                     phases=phases)
#             for ph in picks_hat:
#                 cat[0].picks.append(ph)
#                 Logger.debug(f'Modeled pick "{ph.phase_hint}" for stream "{ph.waveform_id.id}"')
#         # Save subset, augmented catalog to EVENT BANK
#         Logger.info('Submitting processed/augmented event to new eventbank')
#         ABANK.put_events(cat)
        









# if __name__ == '__main__':
#     main()