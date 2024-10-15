"""
:module: BANJO-DETECT/src/python/utils/catalog_utils.py
:auth: Nathan T. Stevens & Barrett Johnson
:email: ntsteven@uw.edu; bnjo@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose: This module contains helper scripts for earthquake catalog
    curation in preparation for EQcorrscan template generation

:attribution: Based on scripts from EQcorrscan documentation and 
    example scripts by B. Johnson.
"""

import logging, os, glob
import obspy
from obsplus import EventBank
import eqcorrscan.utils.catalog_utils as xutil
from eqcorrscan.utils.clustering import catalog_cluster
import pandas as pd

Logger = logging.getLogger(__name__)

# Catalog Construction Methods
def build_ebank(qml_dir, **options):
    """Very simple wrapper for ObsPlus EventBank.__init__

    :param qml_dir: directory containing desired QuakeML files
    :type qml_dir: str
    :return: 
     - **ebank** (*obsplus.EventBank*) - composed event bank
    """    
    ebank = EventBank(base_path=qml_dir, **options)
    return ebank

def assemble_catalog(qml_dir):
    """Assemble all QuakeML files in a specified directory
    into an ObsPy :class:`~obspy.core.event.Catalog` object

    :param qml_dir: directory holding QuakeML files
    :type qml_dir: str
    :return:
     - **cat** (*obspy.core.event.Catalog*) -- composed catalog
    """    
    qml_files = glob.glob(os.path.join(qml_dir, '*.xml'))
    qml_files.sort()
    for _e, _f in enumerate(qml_files):
        if _e == 0:
            cat = obspy.read_events(_f)
        else:
            cat += obspy.read_events(_f)
    return cat

# Catalog Curation Methods
def apply_phase_hints(catalog):
    """Apply phase hints from arrivals to picks in this catalog
    if they are used in the preferred origin of each event in
    this catalog.

    Note: Changes are made in-place, altering the input. If you
    want to preserve your data, use :meth:`~obspy.core.event.Catalog.copy()`
    on your input.

    :param catalog: input catalog
    :type catalog: :class:`~obspy.core.event.Catalog
    :return:
     - **catalog** (*obspy.core.event.Catalog*) -- updated catalog
    """    
    if not isinstance(catalog, obspy.core.event.Catalog):
        raise TypeError('catalog must be type obspy.core.event.Catalog')
    for event in catalog.events:
        # try: 
        origin = event.preferred_origin()
        # except ????: Make sure this isnt hanging
        #     event = event.origins[0]
        arrivals = origin.arrivals
        # Iterate across arrivals and ensure
        for _arr in arrivals:
            # Get 
            phz = _arr.phase
            # Get pick
            pick = _arr.pick_id.get_referred_object()
            # If the pick is present
            if pick != None:
                pick.phase_hint = str(phz)
            else:
                msg = f'Missing pick for ORID: {origin.resource_id} & ARID: {_arr.resource_id}'
                Logger.warning(msg)
    return catalog

def filter_picks(catalog, stations=None, channels=None, networks=None,
                 locations=None, top_n_picks=None, evaluation_mode='all',
                 phase_hints=None, enforce_single_pick=False,
                 min_delta=0, max_delta=180):
    """Augments the EQcorrscan :meth:`~eqcorrscan.utils.catalog_utils.filter_picks` method.
    Filter events in the catalog based on a number of parameters.

    :param catalog: Catalog to filter.
    :type catalog: obspy.core.event.Catalog
    :param stations: List for stations to keep picks from.
    :type stations: list
    :param channels: List of channels to keep picks from.
    :type channels: list
    :param networks: List of networks to keep picks from.
    :type networks: list
    :param locations: List of location codes to use
    :type locations: list
    :param top_n_picks: Filter only the top N most used station-channel pairs.
    :type top_n_picks: int
    :param evaluation_mode:
        To select only manual or automatic picks, or use all (default).
    :type evaluation_mode: str
    :param phase_hints: List of retained phase hints, or None to use all
    :type phase_hints: list
    :param enforce_single_pick:
        Method to enforce using only one pick of each phase-hint per
        station or False to leave all. Can be {False, "earliest"}
    :type enforce_single_pick: str
    :param min_delta: minmum great-circle distance for source-receiver separation
        in degrees, defaults to 0.
    :type min_delta: float, optional.
    :param max_delta: maximum great-circle distances for source-receiver separation
        in degrees, defaults to 180.
    :return:
     - **out** (*obspy.core.event.Catalog*) -- filtered catalog
    """
    # Run distance filtering first
    if min_delta > 0 or max_delta < 180:
        for event in catalog:
            picks = []
            origin = event.preferred_origin()
            if len(event.picks) == 0:
                continue
            for arrival in origin.arrivals:
                if min_delta <= arrival.distance <= max_delta:
                    pick = arrival.waveform_id.get_referred_object()
                    picks.append(pick)
            event.picks = picks

    # Then run the rest of the filtering
    catalog = xutil.filter_picks(
        catalog,stations=stations, channels=channels, networks=networks,
        locations=locations, top_n_picks=top_n_picks,
        evaluation_mode=evaluation_mode, phase_hints=phase_hints,
        enforce_single_pick=enforce_single_pick)

    return catalog


def catalog_cluster(catalog, thresh, metric='distance', show=False):
    """Alias to EQcorrscan :meth:`~eqcorrscan.utils.clustering.catalog_cluster`
    that does spatial or clustering on events in an ObsPy :class:`~obspy.core.event.Catalog`
    object. 

    https://eqcorrscan.readthedocs.io/en/latest/submodules/autogen/eqcorrscan.utils.clustering.catalog_cluster.html#eqcorrscan.utils.clustering.catalog_cluster

    :param catalog: catalog with events to cluster
    :type catalog: :class:`~obspy.core.event.Catalog`
    :param thresh: threshold value (km for *metric="distance"*, sec for *metric="time"*)
    :type thresh: float
    :param metric: threshold metric, either 'distance' or 'time', defaults to 'distance'
    :type metric: str, optional
    :param show: should the clusters be plotted? Defaults to False
    :type show: bool, optional
    :return:
     - **clusters** (*list* of *obspy.core.event.Catalog*) - sub-catalogs containing
        clustered events.
    """    
    clusters = catalog_cluster(catalog, thresh, metric=metric, show=show)
    return clusters

# Bulk Waveform Request Methods
def compose_pick_bulk_line(pick, lead_time=5., lag_time=30., all_components=False):
    """Compose a bulk data request line from an ObsPy :class:`~obspy.core.event.Pick`
    object and a specified amount of lead and lag time around the pick time.

    Provde an option to turn the component code to '?' to get all channels for a
    given station.

    :param pick: pick object
    :type pick: :class:`~obspy.core.event.Pick`
    :param lead_time: seconds of waveform before the pick to include, defaults to 5
    :type lead_time: float, optional
    :param lag_time: seconds of waveform after the pick to include, defaults to 30
    :type lag_time: float, optional
    :param all_components: Should the component code of the pick be replaced with
        a ``"?"`` wildcard? Defaults to False
    :type all_components: bool, optional
    :return:
     - **line** (*tuple*) -- tuple with bulk data reqest formatting
        (Network, Station, Location, Channel, Starttime, Endtime)
    """    
    t0 = pick.time - lead_time
    t1 = pick.time + lag_time
    n,s,l,c = pick.waveform_id.get_seed_string().split('.')
    if all_components:
        c = c[:2]+'?'
    line = (n, s, l, c, t0, t1)
    return line

def compose_origin_bulk_lines(origin, method='origin',lead_time=5,
                              lag_time=30, all_components=False):
    """Compose a bulk waveform request input from an ObsPy
    :class:`~obspy.core.event.Origin` object that has associated 
    :class:`~obspy.core.event.Arrival` and :class:`~obspy.core.event.Pick`
    objects contained within it. 
    
    method = 'origin':
    Time padding referenced to the origin time, using the picks to
    determine the stations to query data from

    method = 'pick':
    Time padding referenced to the timing of each pick. If using
    **all_components**=True, the first occurrence of a channel code
    is used 

    TODO: See if get_waveforms_bulk handles this already....?

    :param origin: origin object
    :type origin: :class:`~obspy.core.event.Origin`
    :param method: method for time referencing, defaults to 'origin'
    :type method: str, optional
    :param lead_time: seconds of waveform before the pick to include, defaults to 5
    :type lead_time: float, optional
    :param lag_time: seconds of waveform after the pick to include, defaults to 30
    :type lag_time: float, optional
    :param all_components: Should the component code of the pick be replaced with
        a ``"?"`` wildcard? Defaults to False
    :type all_components: bool, optional
    :return:
     - **bulk** (*list* of *tuple*) - input for a client bulk waveform request
    """    
    bulk = []
    stachans = []
    if method == 'origin':
        t0 = origin.time - lead_time
        t1 = origin.time + lag_time
        for _arr in origin.arrivals:
            pick = _arr.pick_id.get_referred_object()
            seed = pick.waveform_id.get_seed_string()
            if all_components:
                seed[-1] = '?'
            if seed not in stachans:
                stachans.append(seed)
                n,s,l,c = seed.split()
                bulk.append((n,s,l,c,t0,t1))
    elif method == 'pick':
        for _arr in origin.arrivals:
            pick = _arr.pick_id.get_referred_object()
            line = compose_pick_bulk_line(
                pick,
                lead_time=lead_time,
                lag_time=lag_time,
                all_components=all_components)
            seed = '.'.join(line[:4])
            if seed not in stachans:
                stachans.append(seed)
                bulk.append(line)
    return bulk


