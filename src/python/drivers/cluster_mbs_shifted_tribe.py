"""
:module: Mt_Baker_LF_Research/src/python/drivers/cluster_mbs_tribe.py
:auth: Benz Poobua
:email: spoobu@uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: This driver clusters templates using shifted waveform
    cross-correlation based grouping routines provided
    by EQcorrscan. These routines are re-packaged in the
    `eqcutil.core.clusteringtribe.ClusteringTribe` class.

    It iterates across a set of hyper-parameters that influence
    the calculation of pairwise similarity "distance" values for
    template pairs and saves the results of each clustering analysis
    as a *.tgz file named with the hyper-parameter values

:Hyper Parameters:
    Fill Value (fv) - value used to pad un-paired channels between
        templates. Short-hand name for `replace_nan_distances_with`

        Currently assesses: fv = 1., 'mean'

        Value Explainer
        1 - Apply the maximum penalty to representative correlation
            strength for each missing channel pair.

            "Maximum penality - correlation strength of a missing pair is 0".

        0 - Apply the minimum penalty to representative correlation
            strength for each missing channel pair. WARNING: may tend
            to over-emphasize the importance of correlations between
            templates with very few matching channels.

        'mean' - Fill missing values with the average pair-wise
            maximum correlation amplitude of existing channel pairs.

            "Do not penalize overall correlation strength due to missing pairs"
            
        'min' - Fill missing values with the existing pairs' maximum correlation
            strength (minimum distance)
            "Accentuate the importance of existing channel pairs'
            correlation strength"


        DEBUG: 'min' appears to have some bugs at the eqcorrscan distribution
            level. 
    
    Shift Length (sl) - maximum correlation shift length permitted
        for finding a correlation amplitude maximum for each channel
        pair between templates. `shift_len`. Units of seconds

        Currently assesses: sl = 1., 5., 10.
    
    Lockstep Shifts (ls) - alias for the opposite of `allow_individual_trace_shifts`
        Currently assesses: ls = True, False
        
        Value Explainer
            ls = True - the maximum correlation amplitude (minimum distance)
                is calculated for a single, uniform shift of all traces in one
                template
            ls = False - the maximum correlation amplitude (minimum distance)
                between two templates is calculated from the channel-wise corelation
                amplitude maxima, which may have different shift values
    

"""

import os, glob
from pathlib import Path

from eqcorrscan import Template

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Set up logging
Logger = setup_terminal_logger(name=__name__)
Logger.addHandler(CriticalExitHandler(exit_code=1))
# Get absolute path for repo root directory
root = Path(__file__).parent.parent.parent.parent
# Safety check that it's running from root directory
if os.path.split(root)[-1] != 'Mt_Baker_LF_Research':
    Logger.critical(f'root path does not appear to look like the repo root: {root}')
else:
    Logger.warning(f'running {__file__} from {root}')

# Get path for templates
tdir = os.path.join(root,'processed_data','shifted_templates','well_located_20km')
# Get list of template names
tnames = glob.glob(os.path.join(tdir,'*.tgz'))

# Major parameters
fill_values = [1, 'mean']
shift_lens = [1., 10.]
locksteps = [True, False]
overwrite_protect = False

save_path = os.path.join(root,'processed_data',
                         'shifted_clusters','well_located_20km',
                         'hyper_parameter_trials')

### PROCESSING SECTION ###

# Load Templates
Logger.info(f'Constructing ClusteringTribe from {len(tnames)} saved templates')
ctr = ClusteringTribe()
for tfile in tnames:
    Logger.debug(f'loading: {tfile}')
    itemp = Template().read(tfile)
    ctr.extend(itemp)
Logger.info(f'{len(ctr.templates)} of {len(tnames)} loaded')

for fill_value in fill_values:
    for shift_len in shift_lens:
        for lockstep in locksteps:

            ccckwargs = {'method': 'correlation_cluster',
                        'replace_nan_distances_with': fill_value,
                        'shift_len': shift_len,
                        'corr_thresh': 0.7,
                        'allow_individual_trace_shifts': not lockstep,
                        'show': False,
                        'cores': 'all',
                        'save_corrmat': False}

            Logger.info(f'running fill_value {fill_value} | shift_len {shift_len} sec | lockstep shifts {lockstep}')

            # dynamically generate save directory/file name
            savename = os.path.join(save_path,
                                    f'fv_{fill_value}_sl_{shift_len}_ls_{lockstep}')

            # Overwrite Protection Clause
            if os.path.isfile(f'{savename}.tgz'):
                if overwrite_protect:
                    Logger.info(f'{savename}.tgz already exists - skipping to next')
                    continue
                else:
                    Logger.warning(f'overwriting existing file {savename}.tgz')

            # Run Clustering
            Logger.info(f'Running waveform-based clustering analysis on ClusteringTribe')
            ctr.cluster(**ccckwargs)

            # Save to Disk
            Logger.info(f'Writing results to: {savename}.tgz')
            ctr.write(savename)

# Signal End
Logger.critical('SUCCESSFULLY COMPLETED PROCESS - EXITING')