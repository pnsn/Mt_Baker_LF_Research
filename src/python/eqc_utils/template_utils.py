
"""
:module: BANJO-DETECT/src/python/utils/template_utils.py
:auth: Nathan T. Stevens & Barrett Johnson
:email: ntsteven@uw.edu; bnjo@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose: This module contains helper scripts for EQcorrscan Template
    management in addition to those native to EQcorrscan

:attribution: Based on scripts from EQcorrscan documentation and 
    example scripts by B. Johnson.
"""
import copy, logging, os
from pathlib import Path

import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform
from eqcorrscan.utils.clustering import handle_distmat_nans

tu_logger = logging.getLogger('eqc_utils.template_utils')
tu_logger.setLevel(logging.INFO)

def rename_templates(tribe, prefix=None, inplace=True):
    """Rename templates in a tribe with their PNSN EVID's
    with the option to enact the changes on a deepcopy of th e
    :param tribe: tribe to relabel
    :type tribe: :class:`~eqcorrscan.core.match_filer.Tribe
    :return:
     - **tribe** (*eqcorrscan.core.match_filter.Tribe*) - tribe of renamed templates
    """
    # If not inplace, make a deepcopy of the tribe within the
    # scope of this method
    if not inplace:
        tribe = copy.deepcopy(tribe)
    for template in tribe:
        event = template.event
        oldname = template.name
        name = str(event.resource_id).split(sep='/')[-1]
        if isinstance(prefix, str):
            name = f'{prefix}{name}'
        tu_logger.info(f'renamed template {oldname} -> {name}')
        template.name = name
    return tribe

def augment_template(template, client, padding= 120., min_ncomponents=3):
    """Retrieve additional waveform data for missing channels  for each 
    instrument in a EQcorrscan :class:`~eqcorrscan.core.match_filter.template.Template`
    object's **st** attribute from the provided **client**. If there are
    new trace ID's (NSLC codes) retrieved, they are pre-processed to match
    the sampling rate, bandpass filtering, and timing of the first matching
    trace for the relevant instrument code (NSLC code minus the component character)

    .. rubric:: Explainer
        Say a template has one entry for 'UW.SHUK..BHZ' and the following
        pre-processing steps:
        - 50 Hz sampling rate
        - 0.5 - 20 Hz bandpass filtering (4th order)

        this method will fetch data from the client for 'UW.SHUK..BH?' for 
        the start and end time of the 'UW.SHUK..BHZ' trace--padded by **padding**
        seconds--and apply pre-processing in the following order
        - resample 
            - downsampling uses :meth:`~obspy.core.trace.Trace.resample` with no_filter=False,
            - upsampling uses :meth:`~obspy.core.trace.Trace.interpolate` with method='lanczos')
        - filter
        - trim

    Parameters
    ----------
    :param template: template to augment, warning: modifications are made in-place
    :type template: eqcorrscan.core.match_filter.template.Template
    :param client: client object with a :meth:`~obspy.clients.fdsn.Client.get_waveforms`-type method
    :type client: obspy.clients.fdsn.Client or similar (e.g., obsplus.bank.wavebank.WaveBank)
    :param padding: amount of padding in seconds to add to each retrieved trace
        to capture and exclude filtering edge-effects, defaults to 120.
    :type padding: float-like, optional
    :param min_ncomponents: minimum number of expected components for 
        instruments, defaults to 3
    :type min_ncomponents: int, optional
    :return: augmented template
    :rtype: eqcorrscan.core.match_filter.template.Template
    """    
    # Get unique instrument codes
    insts = {tr.id[:-1] for tr in template.st}
    # Iterate across codes and see if there are missing components
    for inst in insts:
        sub_st = template.st.select(id=inst + '?')
        t0 = sub_st[0].stats.starttime
        t1 = sub_st[0].stats.endtime
        if len(sub_st) < min_ncomponents:
            n, s, l, c = inst.split('.')
            # Fetch waveforms
            ist = client.get_waveforms(network=n, station=s, location=l, channel=c + '?',
                                       starttime=t0 - padding, endtime=t1 + padding)
            # If there are new waveforms
            if len(ist) > len(sub_st):
                # Iterate across all fetched waveforms
                for tr in ist:
                    # If the trace ID is not in sub_st IDs
                    if tr.id not in [tr.id for tr in sub_st]:
                        # Confirm that S/R matches
                        if tr.stats.sampling_rate > template.samp_rate:

                            tu_logger.info(f'downsampling new trace to match template: {tr.id} {tr.stats.sampling_rate} -> {template.samp_rate}')
                            tr.resample(template.samp_rate, no_filter=False)
                        elif tr.stats.sampling_rate < template.samp_rate:
                            tu_logger.info(f'upsampling new trace to match template: {tr.id} {tr.stats.sampling_rate} -> {template.samp_rate}')
                            tr.interpolate(method='laczos')
                        else:
                            pass
                        tu_logger.info(f'Bandpass filtering new trace to match template preprocessing')
                        tr.filter('bandpass',
                                  freqmin = template.lowcut,
                                  freqmax=template.highcut,
                                  corners=template.filt_order)
                        tr.trim(starttime=t0, endtime=t1)
                        template.st += tr
    return template

def compose_template_list(tribe):
    """Compose a template list formatted as an input
    for :meth:`~eqcorrscan.util.clustering.cluster`

    template_list has the format:
     [(template.st, template.name), ...]
     where stream is template.st, and index is a unique
     integer key.

     i.e., a list of 2-tuples each containing: 
            [0] an ObsPy stream
            [1] an unique identifier key
    
    :param tribe: tribe containing templates
    :type tribe: eqcorrscan.core.match_filter.tribe.Tribe
    :return: tlist
    :rtype: list
    """    
    tlist = [(template.st, template.name) for template in tribe]
    return tlist

def save_template_clustering_output(save_dir, groups, savename='clusters'):

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    with open(os.path.join(save_dir, f'{savename}.csv'), 'w') as ofile:
        ofile.write('template_name,group_id\n')
        for _e, group in enumerate(groups):
            for entry in group:
                ofile.write(f'{entry[1]},{_e}\n')

    src = os.path.join(Path().cwd(), 'dist_mat.npy')
    dest = os.path.join(save_dir, f'{savename}_dist_mat.npy')
    if os.path.isfile(src):
        os.rename(src,dest)
    
def reconstitute_dendrogram(distance_matrix_file, corr_thresh=0.4, distance=False, fill_value=1, **kwargs):
    dist_mat = np.load(distance_matrix_file)
    dist_mat = handle_distmat_nans(dist_mat, replace_nan_distances_with=fill_value)
    dist_vect = squareform(dist_mat)
    Z = linkage(dist_vect, **kwargs)
    if distance:
        color_threshold = 1 - corr_thresh
    else:
        color_threshold = corr_thresh

    dendrogram(Z, color_threshold=color_threshold,
               distance_sort='ascending')
    

