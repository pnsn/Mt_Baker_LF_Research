
import os, logging, glob
from pathlib import Path

import numpy as np

from obspy import read_inventory
from obspy.core.event import ResourceIdentifier
from eqcorrscan import Tribe
from obsplus import WaveBank

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.catalog.model_phases import model_picks
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))


def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'workflow' / 'templates'
SAVE_F = TMPD / 'step_6_clustered_{id}_5sec_cct0.45.tgz'
INVD = ROOT / 'data' / 'XML' / 'RESP' / '*.xml'
# Get absolute path to wavebank
WBBP = ROOT / "data" / "WF" / "BANK"

min_snr = 1.1
min_noise_RMS = 3

### SET CORRELATION PARAMETERS ###
ccckwargs = {'method': 'correlation_cluster',
            'replace_nan_distances_with': 'mean',
            'shift_len': 5,
            'corr_thresh': 0.45,
            'allow_individual_trace_shifts': True,
            'show': False,
            'cores': 'all',
            'save_corrmat': False}


### CONNECT TO WAVEBANK
WBANK = WaveBank(WBBP)

### LOAD INVENTORY
for _e, _f in enumerate(glob.glob(str(INVD))):
    if _e == 0:
        INV = read_inventory(_f)
    else:
        INV += read_inventory(_f)

### LOAD PRECLUSTERED TRIBE ###
_f = TMPD/'step3_clustered_5sec_FreeStep.tgz'
print(f'Loading {_f}')
_ctr = ClusteringTribe().read(str(_f))

### Generate New, Single-Station Tribes
ctribes = {}
xtribes = {}
# Create EVID lists keyed to stream names
for evid, row in _ctr._c.iterrows():
    for tr in _ctr.select(evid).st:
        if tr.id not in ctribes.keys():
            ctribes.update({tr.id: [evid]})
        else:
            ctribes[tr.id].append(evid)
xtribes = {}

evid_set = set(_ctr._c.copy().index.values)
for _k, _v in ctribes.items():
    unpicked_set = evid_set.difference(set(_v))
    # Confirm that the intersection is the null set
    if set(_v).intersection(unpicked_set) == set([]):
        # Convert to list
        unpicked_list = list(unpicked_set)
        # Sort
        unpicked_list.sort()
        xtribes.update({_k: unpicked_list})
    else:
        Logger.critical()
    

msg = 'Template Counts Per Channel\nSEED\t\tpicked_count\tto_model'
for _k, _v in ctribes.items():
    msg += f'\n{_k}\t{len(_v)}\t\t{len(xtribes[_k])}'
Logger.info(msg)

# iterate across stream evid lists to form new template objects
for _k, _v in ctribes.items():
    nslc = _k.split('.')
    Logger.info(f'subsetting {_k}')
    _ictr = ClusteringTribe()
    # Create single-channel templates from existing templates
    for _evid in _v:
        # Copy the original trace
        _itmp = _ctr.select(_evid).copy()
        # Get only the specified channel
        _itmp.st = _itmp.st.select(id=_k)
        # Append new single-channel template to new tribe
        _ictr += _itmp
    # Safety check on trace counts
    if not all(len(_tmp.st) == 1 for _tmp in _ictr):
        Logger.critical(f'Stream {_k} produced non-single-trace template(s)')
    # Try to create single-channel templates using ray-tracing to get arrival times
    _x = xtribes[_k]
    for _evid in _x:
        # Get template copy
        _itmp = _ctr.select(_evid).copy()
        # Get origin
        _iorig = _itmp.event.preferred_origin()
        # Check if station exists
        _iinv = INV.select(network=nslc[0],
                           station=nslc[1],
                           channel=nslc[3],
                           time=_iorig.time)
        # If channel does not exist during this event, skip to the next
        if len(_iinv) == 0:
            continue
        # FIXME: Development safety check
        elif len(_iinv) > 1:
            breakpoint()
        # If there is a single entry, model arrival time    
        else:
            pass
        # Model picks
        picks_hat = model_picks(
            _itmp.event.preferred_origin(),
            _iinv,
            phases=['P']
        )
        # Get earliest pick
        for _e, _p in enumerate(picks_hat):
            if _e == 0:
                _pt = _p.time
                _pick_hat = _p
            else:
                if _p.time < _pt:
                    _pt = _p.time
                    _pick_hat = _p
        # # TODO: Use Template construction API to run all of the next steps
        # _itribe = Tribe().construct(
        #     method='from_client',
        #     client_id=WBANK,
        #     lowcut=_itmp.lowcut,
        #     highcut=_itmp.highcut,
        #     filt_order=_itmp.filt_order,
        #     length=45.,
        #     prepick=5.,
        #     samp_rate=_itmp.samp_rate,
        #     parallel=True,
        #     min_snr=1.1,
        #     save_progress=False,
        #     num_cores=12,
        #     process_length=300.

        # )
        # Get data windowed around pick time
        st = WBANK.get_waveforms(
            network=nslc[0],
            station=nslc[1],
            location='*',
            channel=nslc[-1],
            starttime = _pt - _itmp.process_length/2,
            endtime = _pt + _itmp.process_length/2
        )
        # Check that data are present
        if len(st) == 0:
            continue
        nsmp = 0
        for _tr in st:
            nsmp += _tr.count()
        if nsmp < 3000:
            Logger.warning('Loaded data have less than 3000 samples - skipping')
            breakpoint()
            continue
        # Ensure data are merged
        else:
            st.split()
            st.merge(method=1, fill_value='interpolate',interpolation_samples=10)
        # FIXME: If more than one trace is present - safety check
        if len(st) > 1:
            breakpoint()
        if np.ma.is_masked(st[0].data):
            breakpoint()

        # Apply pre-processing
        # Filter
        st.filter('bandpass', freqmin=_itmp.lowcut, freqmax=_itmp.highcut)
        # Resample
        st.resample(_itmp.samp_rate)
        # Trim
        st.trim(starttime=_pt - 5.,
                endtime=_pt + 40. - 1./_itmp.samp_rate)
        # QC that enough data are present
        if st[0].count()/2250 < 0.8:
            Logger.warning(f'{_k} augmentation for {_evid} has insufficient data - skipping')
            continue

        # QC that pick has sufficient SNR
        # FIXME: Rigid implementation
        # Noise in leading 4.75 seconds of 
        noise_amp = _rms(st[0].data[:240])
        onset_amp = _rms(st[0].data[240:400])
        if noise_amp < min_noise_RMS:
            Logger.warning(f'{_k} augmentation for {_evid} has Noise RMS Amp < {min_noise_RMS} - likely dead trace - skipping')
            continue
        if onset_amp/noise_amp < min_snr:
            Logger.warning(f'{_k} augmentation for {_evid} has SNR < {min_snr} ({onset_amp/noise_amp}) - skipping')
            continue
        # Add pick to event
        _itmp.event.picks.append(_pick_hat)
        # Overwrite template stream
        _itmp.st = st
        # # Update template name with appended 'x'
        # _itmp.event.resource_id = ResourceIdentifier(id=_itmp.event.resource_id.id + 'x')
        # _itmp.name += 'x'
        # Append to _ictr
        _ictr += _itmp

    # Map pick quality to clusters attributes
    # _ictr.clusters = _ictr.clusters.assign(manual_pick=[not 'x' in idx for idx in _ictr._c.index])
    _ictr.clusters = _ictr.clusters.assign(manual_pick=['T' if idx in _v else 'F' for idx in _ictr._c.index])
    # Append metadata from _ctr
    _ictr.clusters = _ictr.clusters.join(_ctr.clusters[['etype','time','latitude','longitude','depth','horizontal_uncertainty','vertical_uncertainty']], how='left')
    # Overwrite EVID list with new clusteringtribe
    ctribes[_k] = _ictr
    # Run XCORR
    Logger.info(f'Running correlations on {_k} ({len(_v)} templates)')
    ctribes[_k].cluster(**ccckwargs)
    # Save to disk
    ctribes[_k].write(str(SAVE_F).format(id=_k))


breakpoint()