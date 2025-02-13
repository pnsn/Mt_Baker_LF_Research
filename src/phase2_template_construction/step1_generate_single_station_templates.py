import logging, os
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.signal import iirfilter

from obspy import Stream
from obspy.clients.fdsn.client import Client, FDSNNoDataException
from obsplus import EventBank, WaveBank
from eqcorrscan import Template
from eqcorrscan.utils import pre_processing

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler


# Repository root absolute path
ROOT = Path(__file__).parent.parent.parent
# Path to Augmented Event Bank
AEBBP = ROOT / 'processed_data' / 'catalog' / 'AUGMENTED_BANK'
# Wave Bank Base Path
WBBP = ROOT / 'data' / 'WF' / 'BANK'
# Wave Bank Path Structure
WBPS = '{year}/{julday}'
# Wave Bank Name Structure
WBNS = '{network}.{station}.{location}.{channel}_{year}{julday}_{hour}{minute}{second}.mseed'
# Save Path for Templates
SAVEPATH = ROOT / 'processed_data' / 'template' / 'single_station'
# Template Construction Logging File
CONST_LOG = ROOT / 'processed_data' / 'template' / 'single_station' / 'construct.log'

update_wavebank = False
write_protect = True
last_event_id = None #'quakeml:uw.anss.org/Event/UW/10679373'

# Trace Retrieval Parameters
PREPICK = 5.
DATA_PAD = 15.
PROCESS_LEN = 90.
HIGHCUT = 20.
LOWCUT = 0.5
SAMP_RATE = 50.
FILT_ORDER = 4
LENGTH = 50.
DELAYED = True

# Post Generation Pick SNR QC - Used for noise trace detection/QC
PERCENT_COVERAGE = 75  # Fraction of window that needs viable data
MIN_SNR = 1.2            # Minimum SNR using EQcorrscan's QC method (max amplitude / pre-pick noise RMS)
MIN_DATA_RMS = 3        # Minimum RMS amplitude for the whole trace
MIN_PICK_SNR = 1.2      # Minimum SNR calculated using pick-bounding windows
PICK_SNR_WINDOW = 5.    # [sec] length of pick-bounding windows for SNR calculation

# Post Generation Location Code Aliasing
LOC_ALIASES = {'01':''}

def _rms(array):
    """
    Calculate RMS of array.

    :type array: numpy.ndarray
    :param array: Array to calculate the RMS for.

    :returns: RMS of array
    :rtype: float
    """
    return np.sqrt(np.mean(np.square(array)))



def main():
    if PROCESS_LEN < DATA_PAD*2 + LENGTH:
        Logger.critical('Specified PROCESS_LEN < DATA_PAD + LENGTH + DATA_PAD')

    # Ensure output directory exists
    try:
        os.makedirs(SAVEPATH)
    except FileExistsError:
        pass
    # Connect to clients
    # Augmented Event Bank
    AEBANK = EventBank(base_path=AEBBP)
    # Wave Bank
    WBANK = WaveBank(base_path=WBBP,
                     path_structure=WBPS,
                     name_structure=WBNS)
    # IRIS Webservices Client
    IRIS = Client('IRIS')
    # Read event bank index
    index = AEBANK.read_index()
    # Update index index to event_id
    index.index = index.event_id
    index.index.name = 'event_id'
    # Sort by time
    index = index.sort_values(by='time')
    # Make holder for clustering tribes
    grouper = defaultdict(ClusteringTribe)
    _e = 0
    hit_last = False
    if last_event_id is None:
        hit_last = True
    for event_id, row in index.iterrows():
        _e += 1
        Logger.info(f'processing {event_id} ({_e}/{len(index)})')
        if event_id == last_event_id:
            hit_last = True
        if not hit_last:
            Logger.info(f'yet to hit {last_event_id} - skipping ({_e}/{len(index)})')
            continue
        # Read event
        cat = AEBANK.get_events(event_id=event_id)

        # Create common name for event-named templates
        name = f"{row.agency_id.lower()}{row.event_id.split('/')[-1]}" 
        # Iterate across picks
        _f = 0
        for _p in cat[0].picks:
            _f += 1
            # Get station code
            chanid = _p.waveform_id.id
            nslc = chanid.split('.')
            _year = _p.time.year
            # Get review status
            _rstatus = _p.evaluation_mode

            # Compose target save location
            full_savepath = SAVEPATH / chanid / _rstatus / str(_year)

            if os.path.isfile(str(full_savepath/f'{name}.tgz')):
                if write_protect:
                    Logger.warning(f'{full_savepath}/{name}.tgz already exits - skipping to next')
                    continue
                else:
                    Logger.warning(f'{full_savepath}/{name}.tgz already exits - may overwrite')

            # Create single-pick event
            event = cat[0].copy()
            event.picks = [_p]
            prefor = event.preferred_origin()
            # Cleanout unassociated arrivals
            arr_keep = []
            for arr in prefor.arrivals:
                if arr.pick_id == _p.resource_id:
                    arr_keep.append(arr)
            prefor.arrivals = arr_keep
            # Cleanout unassociated station magnitudes
            stm_keep = []
            for sm in event.station_magnitudes:
                if sm.waveform_id.id == chanid:
                    stm_keep.append(sm)
            event.station_magnitudes = stm_keep

            # Cleanout unassociated amplitudes
            amp_keep = []
            for amp in event.amplitudes:
                if amp.waveform_id.id == chanid:
                    amp_keep.append(amp)
            event.amplitudes = amp_keep
            
            # Compose data request
            if DELAYED:
                t1 = _p.time - PREPICK - DATA_PAD
                t2 = _p.time - PREPICK + LENGTH + DATA_PAD
            else:
                t1 = cat[0].preferred_origin().time - PREPICK - DATA_PAD
                t2 = cat[0].preferred_origin().time - PREPICK + LENGTH + DATA_PAD
            req = nslc + [t1, t2]
            # Try request with wavebank first
            try:
                st = WBANK.get_waveforms(*req)
            except Exception as e:
                Logger.debug(f'Wavebank retrieval failed - trying IRIS')
                st = Stream()
                
            # If empty, try IRIS
            if len(st) == 0: 
                try:
                    st = IRIS.get_waveforms(*req)
                # If no data, continue & log result to file
                except FDSNNoDataException:
                    # CLOG.info(f'{event_id}, FDSNNoDataError')
                    Logger.warning(f'{event_id} {chanid} - Wavebank & IRIS retrievals failed - skipping to next pick')
                    continue
                # Provide option to save IRIS-fetched waveforms locally
                # if update_wavebank and len(st) > 0:
                #     WBANK.put_waveforms(st)

            # EARLY STAGE RESAMPLING FOR ANALOG / OBS STATIONS
            # TO ENSURE ALL COMPONENTS AND TRACE CHUNKS HAVE
            # A SINGLE, INTEGER SAMPLING RATE
            srates = set()
            for tr in st:
                srates.add(tr.stats.sampling_rate)
            # If there are multiple sampling rates
            if len(srates) > 1:
                # Get the rounded mean sampling rate
                mean_srate = np.round(np.mean(list(srates)))
                # FIXME If any of these diverge from mean by more than 1 Hz, break
                if any(np.abs(_e - mean_srate) > 1 for _e in srates):
                    breakpoint()
                # Iterate across traces and resample to the rounded mean sampling rate
                for tr in st:
                    if tr.stats.sampling_rate != mean_srate:
                        Logger.warning(f'Resampling {tr.id} ({tr.stats.sampling_rate} -> {mean_srate})')
                        tr.resample(mean_srate)
            # If there is a single, non-integer sampling rate, resample
            elif list(srates)[0] != np.round(list(srates)[0]):
                for tr in st:
                    tr.resample(np.round(list(srates)[0]))
            # If there is a single, integer sampling rate, do nothing
            else:
                pass

            # THEN Merge traces - use method 1 to take interpolated values for overlaps
            try:
                st.merge(method=1)
            except:
                breakpoint()

            # Breakpoint if more than one trace made it through the merge
            if len(st) != 1:
                breakpoint()
            else:
                tr = st[0]
            
            # Pre-pre-processing data quality control checks
                
            # Run overall data coverage check (EQcorrscan method)
            frac_cover = ((tr.count() - 1)*tr.stats.delta)/(t2 - t1)
            pct_cover = np.round(frac_cover*100, decimals=1)
            if pct_cover < PERCENT_COVERAGE:
                Logger.info(f'Insufficient data for {event_id} {tr.id}: {pct_cover:.1f}% < {PERCENT_COVERAGE:.1f}%')
                continue
            
            # Run mostly zeros check (EQcorrscan method)
            if not pre_processing._check_daylong(tr.data):
                Logger.info(f'{event_id} {tr.id} is mostly zeros - skipping')
                continue
            
            # Test RMS/SNR values
            rms_trace = _rms(tr.copy().data)
            rms_noise = _rms(tr.copy().trim(
                starttime = _p.time - PICK_SNR_WINDOW,
                endtime = _p.time).data)
            rms_onset = _rms(tr.copy().trim(
                starttime = _p.time,
                endtime = _p.time + PICK_SNR_WINDOW).data)
            eqc_snr = max(tr.data) / rms_noise
            pick_snr = rms_onset/rms_noise
            # Run EQcorrscan's SNR QC check for all candidate traces
            if eqc_snr < MIN_SNR:
                Logger.info(f'{event_id} {tr.id} EQcorrscan SNR {eqc_snr} < {MIN_SNR}')
                continue
            else:
                pass
            # Run additional checks on auto-picks
            if _rstatus == 'automatic':
                # Then new Trace RMS check on whole trace, mostly for byte-noise detection
                if rms_trace < MIN_DATA_RMS:
                    Logger.info(f'{event_id} {tr.id} whole trace RMS amplitude too low ({rms_trace} < {MIN_DATA_RMS})')
                    continue
                # Then new Pick SNR on pick-bounded SNR calculations
                elif pick_snr < MIN_PICK_SNR:
                    Logger.info(f'{event_id} {tr.id} pick-centered SNR {pick_snr} < {MIN_PICK_SNR}')
                    continue
                # Passing everything, propagate trace
                else:
                    pass
            # Apply EQCorrscan preprocessing pipeline - TODO Check if starttime and endtime are needed..
            tr = pre_processing.multi_process(
                tr, LOWCUT, HIGHCUT, FILT_ORDER, SAMP_RATE,
                starttime=t1 + DATA_PAD, endtime=t2 - DATA_PAD)
            
            # Alias locations if applicable
            if tr.stats.location in LOC_ALIASES.keys():
                tr.stats.location = LOC_ALIASES[tr.stats.location]
                Logger.warning(f'{event_id} Aliasing location code for {tr.id}')

            # Construct template

            temp = Template(
                name=name,
                st=Stream([tr]),
                lowcut=LOWCUT,
                highcut=HIGHCUT,
                samp_rate=SAMP_RATE,
                filt_order=FILT_ORDER,
                process_length=PROCESS_LEN,
                prepick=PREPICK,
                event=event)
            Logger.info(f"TEMPLATE BUILD SUCCESS - {name} {chanid} {_rstatus} ({_f} / {len(cat[0].picks)})")
            Logger.debug(f'Writing to disk')
            # full_savepath = SAVEPATH / chanid / _rstatus / _year
            # Make full savepath if it doesnt already exist
            try:
                os.makedirs(str(full_savepath))
            except FileExistsError:
                pass
            temp.write(str(full_savepath/name), format='tar')
            # grouper[nslc[1]] += temp
            # breakpoint()

    # for _k, _ctr in grouper.keys():
    #     try:
    #         os.makedirs(SAVEPATH/_k)
    #     except FileExistsError:
    #         pass
    #     _ctr.write(str(SAVEPATH/_k/f'{_k}_unclustered'), compress=True)
    
    # return grouper


if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
    Logger.addHandler(CriticalExitHandler(exit_code=1))
    # Set up to-file logger
    # CLOG =

    main()           
    #         # Replicate eqcorrscan.util.pre_processing.multi_process in serial
    #         # 0. Enforce double-precision dtype
    #         if not tr.data.dtype == np.float64:
    #             tr.data = tr.data.astype(np.float64)
    #         # 1. Fill gaps with 0's
    #         if isinstance(tr.data, np.ma.MaskedArray):
    #             mask = tr.data.mask
    #             tr = tr.split()
    #             tr = tr.detrend().merge(fill_value=0)[0]
    #         else:
    #             mask = None
    #         # 2. Check for zeros (should have been done already?)
    #         # 3. detrend
    #         tr.detrend('simple')
    #         # 4. zero-pad to full processing length
    #         tr.trim(starttime=t1, endtime=t2, pad=True, fill_value=0)
    #         # 5. Resample
    #         tr.resample(SAMP_RATE)
    #         # 6. Filter
    #         if LOWCUT and HIGHCUT:
    #             tr.filter('bandpass', freqmin=LOWCUT, freqmax=HIGHCUT, corners=FILT_ORDER)
    #         elif LOWCUT:
    #             tr.filter('highpass', freq=LOWCUT, corners=FILT_ORDER)
    #         elif HIGHCUT:
    #             tr.filter('lowpass', freq=HIGHCUT, corners=FILT_ORDER)
    #         else:
    #             pass
    #         # 7. Reapply zeros to padding samples - these get trimmed in 8...?
    #         # 9. Zero Pad Gaps
    #         if mask is not None:
    #             tr.data = np.ma.MaskedArray(
    #                 data = tr.data,
    #                 mask= mask,
    #                 fill_value=0
    #             )
    #             tr = tr.split()
    #             tr = tr.detrend()
    #             tr = tr.merge(fill_value=0)[0]

    #         # 8. Trim to template length ("Recheck length") & pad with zeros
    #         tr.trim(starttime=t1 + DATA_PAD,
    #                 endtime=t2 - DATA_PAD,
    #                 pad=True,
    #                 fill_value=0.)
    #         if len()


    #         # Add to event stream for multiprocessing
    #         st_event += st
        
    #     # RUN MULTI THREADED EQCORRSCAN PREPROCESSING PIPELINE
    #     st_pp = pre_processing.multi_process(
    #         st_event, LOWCUT, HIGHCUT, FILT_ORDER, SAMP_RATE
    #     )
            


           
            
            
    #         # Enforce double-precision


    #         # Preprocess trace(s)
    #         if LOWCUT and HIGHCUT:
    #             st1.filter('bandpass', freqmin=LOWCUT, freqmax=HIGHCUT, corners=FILT_ORDER)
    #         elif LOWCUT:
    #             st1.filter('highpass', freq=LOWCUT, corners=FILT_ORDER)
    #         elif HIGHCUT:
    #             st1.filter('lowpass', freq=HIGHCUT, corners=FILT_ORDER)
    #         else:
    #             pass

    #         if any(tr.stats.sampling_rate != SAMP_RATE for tr in st1):
    #             for tr in st1:
    #                 if tr.stats.sampling_rate != SAMP_RATE:
    #                     tr.resample(SAMP_RATE)

    #             # Check for mostly-zero traces
            
    #         breakpoint()
    #         # Fetch Data from client
            
    #         # try:
    #         #     st = WBANK(*nslc.split('.'),
    #         #                _p.time - prepick - pad,)


    #         # # Filter copy of event for just those picks
    #         # _cat = filter_picks(cat.copy(), stations=[sta], enforce_single_pick='earliest')
    #         # # _kcat += _cat
    #         # gen_kwargs.update({'client_id': WBANK,
    #         #         'catalog': _cat})
    #         # # Try template construction with wavebank
    #         # # try:
    #         # template = template_gen(**gen_kwargs)
    #         # grouper[sta] += template
    #         # # except:
    #             # breakpoint()
    #         # breakpoint()    
    # breakpoint()


