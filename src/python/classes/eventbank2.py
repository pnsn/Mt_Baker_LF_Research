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
from obsplus import EventBank
from obsplus.bank.eventbank import compose_docstring
from obsplus.constants import get_events_parameters


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
        """
        # Run as normal
        df = super().read_index(**options)
        # Update index
        df.index = [row.agency_id.lower() + row.event_id.split('/')[-1] for _, row in df.iterrows()]
        # Update index name
        df.index.name = 'COMCAT_ID'
        return df
    
    
    