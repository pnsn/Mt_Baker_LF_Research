"""
:module: Mt_Baker_LF_Research/src/python/drivers/align_mbs_traces.py
:auth: Benz Poobua
:email: spoobu@uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: This driver aligns traces relative to one another based on their cross-correlation value.

"""

import os, glob
from pathlib import Path
from obspy import Catalog
from eqcorrscan.utils.stacking import align_traces

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

# Get path for ClusteringTribe
processed_path = os.path.join(root,'processed_data','clusters','well_located_20km')
ctr_file_path = os.path.join(processed_path,'hyper_parameter_trials')
ctr_files = glob.glob(os.path.join(ctr_file_path, 'fv_mean_sl_1.0_ls_False.tgz'))

# Load the ClusteringTribe
filepath = ctr_files[0]
ctr = ClusteringTribe().read(filepath)
Logger.info('ClusteringTribe loaded')
ctr2 = ctr.copy() # 155 templates

# Set the ccv threshold of the alignment result
ccv_threshold = 0.7

# Initialize an empty catalog (shifted)
all_events = []

# Save path for subcatalog 
subcat_path = os.path.join(root,'processed_data','sub_catalog')
os.makedirs(subcat_path, exist_ok=True)
filename = f'sub_cat_ccv_threshold_{ccv_threshold}.xml'
subcat_filepath = os.path.join(subcat_path, filename)

for group_no in ctr2.clusters.correlation_cluster.unique():
    sub_ctr = ctr2.select_cluster(method='correlation_cluster', index=group_no)
    if len(sub_ctr) >= 3: # the number of templates
        trace_list = {}
        trace_list_metadata = {}
        for template in sub_ctr.templates:
            for tr in template.st:
                # if tr.id not in trace_list.keys():
                #     trace_list.update({tr.id: [tr]})
                #     trace_list_metadata.update({tr.id: [template.name]})
                if tr.id not in trace_list.keys():
                    trace_list[tr.id] = [tr]
                    trace_list_metadata[tr.id] = [template.name]
                else:
                    trace_list[tr.id].append(tr)
                    trace_list_metadata[tr.id].append(template.name)
        keys_to_remove = [] # Prevent RuntimeError (dictionary changed size)
        for _k, _v in trace_list.items(): # Iterate across trace lists and remove lists with fewer than 3 traces
            if len(_v) < 3:
                keys_to_remove.append(_k)
            else:
                realignment_results = align_traces(_v, shift_len=1, 
                                                   master=False, positive=True)
                shifts, ccvs = realignment_results  # Unpack realignment_results into two lists
                for name, shift, ccv in zip(trace_list_metadata[_k], shifts, ccvs):
                    if ccv >= ccv_threshold:
                        # Get a relevant template
                        template = sub_ctr.select(name) 
                        # Get the template's event
                        event = template.event 
                        # Iterate across picks to find matching
                        # No need to augment the waveforms in `create_mbs_tribe.py`
                        # Access the picks attribute of the Event
                        for pick in event.picks:
                            if pick.waveform_id: # Ensure waveform_id exists 
                                pick.time += shift  # Apply the shift to the pick's time

    # Remove the keys (less than 3 traces)
    for key in keys_to_remove:
        trace_list.pop(key, None)
        trace_list_metadata.pop(key, None)
        
    # Add the events from sub_ctr to all_events
    all_events.extend([template.event for template in sub_ctr.templates])
    # sub_cat = Catalog(events=[temp.event for temp in sub_ctr.templates])

# Construct catalog after looping
sub_cat = Catalog(events=all_events)

# Save the sub catalog
sub_cat.write(subcat_filepath, format="QUAKEML")





# # Example to check if the shift is applied successfully
# # Extract the event ID from cat (e.g., '61138827')
# cat_event_id = cat.events[10].resource_id.id.split('/')[-1]  # Extracts '61138827'

# # Find the corresponding event in sub_cat
# found_event_in_sub_cat = None
# for sub_cat_event in sub_cat:
#     sub_cat_event_id = sub_cat_event.resource_id.id.split('/')[-1]  # Extract event ID from sub_cat
    
#     if sub_cat_event_id == cat_event_id:
#         found_event_in_sub_cat = sub_cat_event
#         break

# if found_event_in_sub_cat:
#     print(f"Comparing pick times for Event ID: {cat_event_id}")
    
#     # Print pick times from the original catalog (cat)
#     print("\nPick times from cat (original):")
#     for pick in cat.events[10].picks:
#         print(f"  {pick.phase_hint} pick time: {pick.time}")
    
#     # Print pick times from the sub catalog (sub_cat)
#     print("\nPick times from sub_cat (shifted):")
#     for pick in found_event_in_sub_cat.picks:
#         print(f"  {pick.phase_hint} pick time: {pick.time}")
# else:
#     print(f"No matching event found for Event ID {cat_event_id} in sub_cat.")

# ctr2.select_template_traces(station='MBW*')
# ctr[1].st.select(station='MBW')