import os, logging
from pathlib import Path

from obspy import UTCDateTime, Stream
from obspy.clients.fdsn import Client
from obsplus import WaveBank, EventBank

from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
Logger.addHandler(CriticalExitHandler(exit_code=1))

# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to wavebank
WBBP = ROOT/'data'/'WF'/'BANK'
WBPS = '{year}/{julday}'
WBNS = '{seedid}_{year}{julday}_{hour}{minute}{second}'

# State which station and channel(s) will be queried
NET = 'UW'
STA = 'MBW2'
LOC = '*'
CHA = '??Z'
DT = 300.

sampling_rate_precision = 10

# Connect to EventBank
EBANK = EventBank(base_path=EBBP)
# Connect to WaveBank
WBANK = WaveBank(base_path=WBBP,
                 path_structure=WBPS,
                 name_structure=WBNS)
# Connect to IRIS DMC Client
client = Client("IRIS")

# Get Event Metadata
df_eb = EBANK.read_index()

for _e, row in df_eb.iterrows():
    # Get origin time
    t0 = UTCDateTime(row.time.isoformat())
    Logger.info(f'{t0} - ({_e+1} of {len(df_eb)})')
    # Compose bulk request
    bulk = []
    for _e in range(-2,4):
        t1 = t0 + _e*DT
        t2 = t1 + DT
        # Try to get data from WBANK first
        st = WBANK.get_waveforms(NET, STA, LOC, CHA, t1, t2)
        # If the wavebank has some data for this request
        if len(st) > 0:
            Logger.info('Found waveforms for a bulk request in WaveBank')
            # De-float sampling rates with sync interpolation
            for _tr in st:
                _sr0 = _tr.stats.sampling_rate
                if _sr0 != int(_sr0):
                    _sr1 = round(_sr0/sampling_rate_precision)*sampling_rate_precision
                    Logger.info(f'Sampling rate for {tr.id} is {tr.stats.sampling_rate:0.4f} - interpolating to set to {_sr1}')
                    _tr.interpolate(sampling_rate=_sr1)
            # Make sure traces are merged
            st.merge()
            # Iterate across traces
            for tr in st:
                # Map contiguous intervals
                intervals = []
                # Split each trace to identify gaps in a given stream
                for itr in tr.split():
                    intervals.append([itr.stats.starttime, itr.stats.endtime])
                # If one interval (potentially incomplete data)
                if len(intervals) == 1:
                    # If full coverage, skip
                    if intervals[0][0] == t1 and intervals[0][1] == t2:
                        Logger.info(f'WAVEBANK already covers {t1} - {t2} for {tr.id} - skipping')
                        continue
                    # If missing leading samples
                    if intervals[0][0] > t1:
                        bulk.append(tr.id.split('.') + [t1, intervals[0][0]])
                    # Missing tailing samples
                    if intervals[0][1] < t2:
                        bulk.append(tr.id.split('.') + [intervals[0][1], t2])
                # If multiple intervals (gappy data)
                elif len(intervals) > 1:
                    for _e, interval in enumerate(intervals):
                        if _e == 0:
                            if interval[0] > t1:
                                bulk.append(tr.id.split('.') + [t1, interval[0]])
                        if _e == len(intervals) - 1:
                            if interval[1] < t2:
                                bulk.append(tr.id.split('.') + [interval[1], t2])
                        if _e < len(intervals) - 1:
                            bulk.append(tr.id.split('.') + [interval[1], intervals[_e+1][0]])
        # If no data are present, request full stretch
        else:
            Logger.info('No data in WAVEBANK - composing reqest for webservices')
            bulk.append((NET, STA, LOC, CHA, t1, t2))
    # Create empty holder stream
    st = Stream()
    # Try to run bulk request
    if len(bulk) == 0:
        Logger.info('empty bulk request - skipping to next')
        continue
    try:
        st += client.get_waveforms_bulk(bulk)
        Logger.info('successful bulk query')
    # Failing that, try to run single entry request
    except:
        Logger.warning('attempting individual query')
        for _b in bulk:
            try:
                st += client.get_waveforms(*_b)
            except:
                continue
    # If data were queried
    if len(st) > 0:
        # Split to remove gaps
        st.split()
        # Submit to wavebank
        try:
            WBANK.put_waveforms(st)
            Logger.info('successful bulk submit to WaveBank')
        # Failing that, try to submit individual traces
        except:
            Logger.warning('Stream submit to wavebank unsuccessful - attempting individual submit to WaveBank')
            for tr in st:
                try:
                    WBANK.put_waveforms(tr)
                except:
                    Logger.warning(f'unsuccessful single submit to wavebank for {tr.id} - skipping to next')
                    continue
    # If no data were queried, continue
    else:
        Logger.warning('empty stream')
        continue




    

