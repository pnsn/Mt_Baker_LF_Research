"""
:module: Mt_Baker_LF_Research/src/python/classes/eventbank2.py
:auth: Nathan T. Stevens
:email: ntsteven@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose:
    This module extends functionalites of the ObsPlus :class:`~obsplus.bank.eventbank.EventBank` class
    in a manner consistent with their GPLv3 licensing.
    see: https://github.com/niosh-mining/obsplus

    
"""
import fnmatch
import obspy
import pandas as pd
from obsplus import EventBank
from obsplus.bank.eventbank import compose_docstring
from obsplus.constants import get_events_parameters

def _parse_event_fix_status(event: obspy.core.event.Event) -> dict:
    """Parse the fixed solution status of an :class:`~obspy.core.event.Event` object
    and return a dictionary formatted with standard AQMS fields plus 'event_id'
    
    :param event: Event object to parse
    :type event: obspy.core.event.Event
    :return: 
     - **out** (*dict*) - dictionary with 
    """    
    event_id = str(event.resource_id)
    try:
        origin = event.preferred_origin()
    except AttributeError:
        origin = event.origins[0]
    if origin.depth_type == 'from location':
        fdepth = False
    else:
        fdepth = True
    ftime = origin.time_fixed
    fepi = origin.epicenter_fixed
    out = {'event_id': event_id,'fepi': fepi,'fdepth':fdepth,'ftime': ftime}
    return out

def _parse_event_min_sta_dist(event: obspy.core.event.Event, unit='m') -> dict:
    if unit.lower() in ['m','meter','km','kilometer','d','deg','degree']:
        unit = unit.lower()
    else:
        raise ValueError(f'unit value {unit} not supported.')
    
    event_id = str(event.resource_id)
    try: 
        origin = event.preferred_origin()
    except AttributeError:
        origin = event.origins[0]
    # Start with unrealistic degree distance
    mindist = 999
    for _arr in origin.arrivals:
        if hasattr(_arr, 'distance'):
            if isinstance(_arr.distance, float):
                if _arr.distance < mindist:
                    mindist = _arr.distance
                    pick = _arr.pick_id.get_referenced_object()
                    phz = _arr.phase
                    seed = pick.waveform_id.get_seed_string()
    if unit in ['km','kilometer']:
        mindist *= 111.2
        unit = 'km'
    elif unit in ['m','meter']:
        mindist *= 111.2e3
        unit = 'm'
    else:
        unit ='deg'
    out = {'event_id': event_id, f'mindist_{unit}': mindist,
           'mindist_phz': phz, 'mindist_chan': seed}
    return out

    
    


class EventBank2(EventBank):
    """An extension of the :class:`~obsplus.bank.eventbank.EventBank` class
    that adds some specific functionalities relevant to event location quality
    and other helper functions

    Parameters
    ----------
    :param EventBank: _description_
    :type EventBank: _type_
    """    
    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    @compose_docstring(get_events_params=get_events_parameters)
    def read_index(self, **options):
        """An augmentation of the :meth:`~obsplus.bank.eventbank.EventBank.read_index`
        method (and by extension the :meth:`~.get_event_summary` method).
        
        Augmentations
        -------------
        - Index is assigned to numeric EVID's that are taken as the last
          `/`-delimited element of **event_id** for each event in the EventBank
          and the agency_id.lower()
          e.g., 'quakeml:uw.anss.org/Event/UW/60593167'
           becomes: 'uw60593167'
          This conforms to the nomenclature in the US Geological Survey
          Comprehensive Catalog (COMCAT).

        Parameters
        ----------
        {get_events_params}

        # TODO - make a nearest station distance option extension (meter units)
        """
        # Restrict access to COMCAT_ID as a search parameter & provide hint
        if 'COMCAT_ID' in options:
            raise ValueError('use `event_id` column entries to fetch subset events by EVID')
        # If using new kwarg, pop off this kwarg from options to allow for expected inputs
        # to original read_index() method. If not using this, set internal var to False
        if 'include_fixed_status' in options.keys():
            ifs = options.pop('include_fixed_status')
        else:
            ifs = False
        
        if 'include_mindist' in options.keys():
            mds = options.pop('include_mindist')
        else:
            mds = False
        if 'unit' in options.keys():
            unit = options.pop('unit')
        else:
            unit = 'm'

        # Enable wildcard search for event_id
        # If event_id is in **options
        if 'event_id' in options:
            # If event_id's value is a string
            if isinstance(options['event_id'], str):
                # if that string has wildcard character(s) present
                if '?' in options['event_id'] or '*' in options['event_id']:
                    # pop the event_id entry off options
                    wild_entry = options.pop('event_id')
                    # run everyting else as normal
                    df = super().read_index(**options)
                    event_id=fnmatch.filter(df.event_id.values, wild_entry)
                    options.update({'event_id': event_id})

        # Run as normal
        df = super().read_index(**options)

        # Handle fetching fixed status entries
        if len(df) > 0:
            entries = []
            for eid in df.event_id:
                breakpoint()
                event = self.get_events(event_id=eid)[0]
                e_row = {}
                # Aggregate
                if ifs:
                    e_row.update(_parse_event_fix_status(event))
                if mds:
                    e_row.update(_parse_event_fix_status(event, unit=unit))
                entries.append(e_row)
            # Create dataframe
            entries = pd.DataFrame(entries)
            # Merge dataframe
            df = pd.merge(df, entries, on='event_id', how='left')
        
       
        # Update index with COMCAT_ID
        df.index = [row.agency_id.lower() + row.event_id.split('/')[-1] for _, row in df.iterrows()]
        # Update index name
        df.index.name = 'COMCAT_ID'
                    
        return df
    

    