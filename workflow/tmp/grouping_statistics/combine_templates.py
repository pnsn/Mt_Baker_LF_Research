# Imports
import logging, glob, os

from pathlib import Path

import pandas as pd
import numpy as np

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import basic_logger_config

basic_logger_config(level=logging.INFO)
Logger = logging.getLogger(os.path.split(__file__)[-1])

# Map absolute paths
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow" / "templates" 
# Template Save Directory
CTSD = PD_DIR / "single_station" / 'xcc_test'
# Output File (also temporary save file)
OUTFILE = PD_DIR / 'MERGED.tgz'
FINAL_OUTFILE = PD_DIR / 'MERGED_CLUSTERED.tgz'

ccckwargs = {'method': 'xcc','replace_nan_distances_with': 'mean','shift_len': 2.,'corr_thresh': 0.5,'allow_individual_trace_shifts': True,'cores': 'all'}


Logger.info(f'Will read from {str(CTSD)}')
file_list = glob.glob(str(CTSD/'*.tgz'))
nfiles = len(file_list)
CTR = ClusteringTribe()
# Create placeholder for distmat
distmat = np.full(shape=(700,700, nfiles), fill_value=np.nan)
names = set()
for _e, file in enumerate(file_list):
    Logger.info(f'loading {os.path.split(file)[-1]}')
    # Load saved single-station clustering tribe
    ctr = ClusteringTribe().read(file)
    Logger.info(f'single station CTR has {len(ctr)} templates')
    # Iterate across templates
    for other in ctr:
        # If this is a new template name, use iadd
        if other.name not in names:
            # Append to CTR
            CTR += other
            # Add names to name set
            names.add(other.name)
            
        # Otherwise, splice traces and picks
        else:
            initial = CTR[other.name]
            initial.st += other.st
            for pick in other.event.picks:
                initial.event.picks.append(pick)

    # Join id_no's to get mapping for distance matrix elements        
    ijmap = CTR._c.join(ctr._c.id_no, how='right', rsuffix='_j')
    # Iterate across all indices
    for ii, irow in ijmap.iterrows():
        ii_from = irow.id_no_j
        ii_to = irow.id_no
        for jj, jrow in ijmap.iterrows():
            jj_from = jrow.id_no_j
            jj_to = jrow.id_no
            # put value in new distmat
            try:
                distmat[ii_to, jj_to, _e] = ctr.dist_mat[ii_from, jj_from]
            except:
                breakpoint()
    


    Logger.info(f'amalgamated CTR now has {len(CTR)} templates')
    
    Logger.info(f'writing to disk: {str(OUTFILE)}')
    CTR.write(str(OUTFILE), compress=True)


# Run Clustering
CTR.cluster(**ccckwargs)
CTR.write(str(FINAL_OUTFILE), compress=True)

breakpoint()
distmat = distmat[:len(CTR), :len(CTR), :]
# Calculate the nanmean distance 
distmat_mean = np.nanmean(distmat, axis=2)
# Calculate the nanmedian distances
distmat_median = np.nanmedian(distmat, axis=2)
