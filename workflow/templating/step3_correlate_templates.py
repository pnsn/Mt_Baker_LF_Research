import os, logging, glob
from pathlib import Path

import pandas as pd

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1])
Logger.addHandler(CriticalExitHandler())
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent

# Working directory absolute path
WDIR = ROOT / "processed_data" / "workflow" / "templates"
# ClusteringTribe name
CTRIBE = WDIR / 'step2_unclustered_tribe.tgz'
# Output Name
ONAME = WDIR / 'step3_clustered_5sec_FreeStep.tgz'

ccckwargs = {'method': 'correlation_cluster',
            'replace_nan_distances_with': 'mean',
            'shift_len': 5,
            'corr_thresh': 0.5,
            'allow_individual_trace_shifts': True,
            'show': False,
            'cores': 'all',
            'save_corrmat': False}


Logger.info(f'loading {CTRIBE}')
CTR = ClusteringTribe().read(str(CTRIBE))
# Ensure that all template traces are merged!
for tmp in CTR.templates:
    lenst = len(tmp.st)
    tmp.st.merge(method=1, interpolation_samples=1)
    lenstm = len(tmp.st)
    if lenstm < lenst:
        Logger.info(f'merged {lenst - lenstm} traces in {tmp.name}')
Logger.info('starting clustering')
CTR.cluster(**ccckwargs)
Logger.info(f'saving to disk as {ONAME}')
if os.path.exists(str(ONAME)):
    Logger.warning(f'{ONAME} already exists! Breakpoint!')
    breakpoint()
    
CTR.write(str(ONAME))


