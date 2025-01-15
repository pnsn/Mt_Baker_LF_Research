import os, logging
from pathlib import Path

from eqcorrscan import Template

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging at debug level reporting
Logger = setup_terminal_logger(os.path.split(__file__)[-1],level=logging.INFO)
# Add a critical -> sysexit handler custom class
Logger.addHandler(CriticalExitHandler(exit_code=1))

### SET PATHS ###
ROOT= Path(__file__).parent.parent.parent
if os.path.split(ROOT)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {ROOT}')
TMPD = ROOT / 'processed_data' / 'workflow' / 'templates'
SAVE_F = TMPD / 'step_5_clustered_{id}_5sec_cct0.45.tgz'
### SET CORRELATION PARAMETERS ###
ccckwargs = {'method': 'correlation_cluster',
            'replace_nan_distances_with': 'mean',
            'shift_len': 5,
            'corr_thresh': 0.45,
            'allow_individual_trace_shifts': True,
            'show': False,
            'cores': 'all',
            'save_corrmat': False}

### LOAD PRECLUSTERED TRIBE ###
_f = TMPD/'step3_clustered_5sec_FreeStep.tgz'
print(f'Loading {_f}')
_ctr = ClusteringTribe().read(str(_f))

### Generate New, Single-Station Tribes
ctribes = {}
# Create EVID lists keyed to stream names
for evid, row in _ctr._c.iterrows():
    for tr in _ctr.select(evid).st:
        if tr.id not in ctribes.keys():
            ctribes.update({tr.id: [evid]})
        else:
            ctribes[tr.id].append(evid)
msg = 'Template Counts Per Channel\nSEED\t\tcount'
for _k, _v in ctribes.items():
    msg += f'\n{_k}\t{len(_v)}'
Logger.info(msg)

# iterate across stream evid lists to form new template objects
for _k, _v in ctribes.items():
    Logger.info(f'subsetting {_k}')
    _ictr = ClusteringTribe()
    for _evid in _v:
        # Copy the original trace
        _itmp = _ctr.select(_evid).copy()
        # Get only the specified channel
        _itmp.st = _itmp.st.select(id=_k)
        # Keep only the informing pick
        _pk = False
        for _p in _itmp.event.picks:
            if _p.waveform_id.id == _k:
                _pk = _p
                break
            else:
                continue
        # if _pk:       
        #     _itmp.event.picks = [_pk]
        # # Safety catch
        # else:
        #     breakpoint()
        # Append new single-channel template to new tribe
        _ictr += _itmp
    # Safety check on trace counts
    if not all(len(_tmp.st) == 1 for _tmp in _ictr):
        Logger.critical(f'Stream {_k} produced non-single-trace template(s)')
    # Append metadata from _ctr
    _ictr.clusters = _ictr.clusters.join(_ctr.clusters[['etype','time','latitude','longitude','depth','horizontal_uncertainty','vertical_uncertainty']], how='left')
    # Overwrite EVID list with new clusteringtribe
    ctribes[_k] = _ictr
    # Run XCORR
    Logger.info(f'Running correlations on {_k} ({len(_v)} templates)')
    ctribes[_k].cluster(**ccckwargs)
    # Save to disk
    ctribes[_k].write(str(SAVE_F).format(id=_k))



