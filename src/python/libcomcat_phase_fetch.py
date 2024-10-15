"""
:module: Mt_Baker_LF_Research/src/python/libcomcat_phase_fetch.py
:auth: Nathan T. Stevens
:email: ntsteven (at) uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU v3
:attribution: This script is based on the libcomcat DetailEvent.ipynb
 that can be found here: https://code.usgs.gov/ghsc/esi/libcomcat-python/-/blob/main/notebooks/DetailEvent.ipynb

:purpose: This script iterates across provided PNSN/UW EVID's and attempts
to retrieve
"""

import os, sys, logging, argparse
import pandas as pd
from obspy import Catalog
from obspy.io.quakeml.core import Unpickler
from libcomcat.search import get_event_by_id

Logger = logging.getLogger('libcomcat_phase_fetch.py')

def get_libcomcat_from_evid(evid, version='preferred'):
    """Core process - given an Event ID (EVID) from the PNSN catalog,
    fetch origin information for the specified version (see supported 
    values) and attempt to fetch phase data using the U.S. Geological
    Survey Earthquake COMCAT (Comprehensive Catalog) database. Return
    the resultant event metadata as a :class:`~obspy.core.event.Catalog`
    object

    :param evid: the event ID for a PNSN catalog event
    :type evid: int or str
    :param version: _description_, defaults to 'preferred'
    :type version: str, optional
    :return:
     - **cat** (*obspy.core.event.Catalog*) - Catalog object containing
        event, origin, and pick (meta)data
    """
    if isinstance(evid, int):
        evid = f'uw{evid}'
    elif isinstance(evid, str):
        if evid[:2] == 'uw':
            pass
        else:
            evid = f'uw{evid}'
    else:
        raise TypeError(f'evid type {type(evid)} not supported. Must be type int or str')
    
    if version not in ['first','last','latest','preferred']:
        raise ValueError(f'version "{version}" not supported.')

    # Initialize catalog & unpickler
    cat = Catalog()
    unpickler = Unpickler()
    try:
        detail = get_event_by_id(evid)
    except:
        Logger.error(f'could not get event detail for {evid}')
        origins = []
    try:
        origins = detail.getProducts('phase-data', source='uw', version=version)
    except:
        if args.exclude_phaseless:
            Logger.error(f'could not get phase-data for event {evid} - skipping')
            origins = []
        else:
            Logger.error(f'could not get phase-data for event {evid} - just getting origin information')
            origins = detail.getProducts('origin', source='uw', version=version)
    for origin in origins:
        qbytes, url = origin.getContentBytes('xml')
        catalog = unpickler.loads(qbytes)
        cat += catalog
    return cat

def main(args):
    """Main process - loads specified CSV file,
    filters by minimum and maximum offsets from Mt. Baker's
    summit, and interatively fetches phase dta from 

    :param args: _description_
    :type args: _type_
    """    
    # Load CSV
    df = pd.read_csv(args.input_csv, index_col='evid')
    # Filter by distance from Mt. Baker
    df = df[(df.mbs_distance_km <= args.max_radius_km) &\
            (df.mbs_distance_km >= args.min_radius_km)]
    # Make sure all rows have unique EVIDs (in case the phase file is read)
    evids = df.index.unique()
    if len(evids) != len(df):
        Logger.critical('not all rows in input CSV have unique EVIDs - exiting on 1')
        sys.exit(1)
    # Send starting report
    msg = f'AQMS dataframe filtered for offsets {args.min_radius_km} - '
    msg += f' {args.max_radius_km} km offset: {len(evids)} EVIDS'
    Logger.warning(msg)

    # Initialize catalog
    cat = Catalog()
    # Get output directory
    write_dir, _ = os.path.split(args.output_xml)
    # Create directory if it doesn't exist
    if not os.path.exists(write_dir):
        Logger.warning(f'creating output directory: {write_dir}')
        os.makedirs(write_dir)
    # Save input args to output directory
    with open(os.path.join(write_dir,'input_args.txt'), 'w') as _l:
        _l.write(str(args))
    
    for _e, evid, row in enumerate(df.iterrows()):
        icat = get_libcomcat_from_evid(evid)
        if len(icat) > 0:
            cat += icat
            Logger.debug(f'EVID: {evid} done - {_e+1} of {len(evids)}')
        # Enact incremental save-point
        if (_e+1) % args.nsavepoint == 0:
            cat.write(os.path.join(write_dir,'cat_savepoint.xml'), format='QUAKEML')
            Logger.info(f'Savepoint at {_e+1} of {len(evids)} complete.')
    cat.write(args.output_xml, format='QUAKEML')
    os.remove(os.path.join(write_dir,'cat_savepoint.xml'))


if __name__ == '__main__':
    ch = logging.StreamHandler()
    fmt = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    ch.setFormatter(fmt)
    Logger.addHandler(ch)

    parser = argparse.ArgumentParser(
        prog='libcomcat_search.py',
        description='Fetch event phase pick data for PNSN/UW events using libcomcat and a CSV of EVIDs')
    
    parser.add_argument(
        '-i',
        '--input',
        action = 'store',
        dest='input_csv',
        default= os.path.join('data','Events','MtBaker_50km_radius_origins.csv'),
        help='input CSV file with desired EVIDs'
    )

    parser.add_argument(
        '-o',
        '--output',
        action = 'store',
        dest='output_xml',
        default=os.path.join('data','XML','QUAKE','libcomcat_events.xml'),
        help='output QuakeML file'
    )

    parser.add_argument(
        '-R',
        '--max_radius_km',
        action = 'store',
        dest='max_radius_km',
        default=10,
        type=float,
        help='maximum origin from Mt. Baker summit in km'
    )

    parser.add_argument(
        '-r',
        '--min_radius_km',
        action = 'store',
        dest='min_radius_km',
        default=-1,
        type=float,
        help='minimum origin distance from Mt. Baker summit in km'
    )

    parser.add_argument(
        '-v',
        '--verbose',
        action = 'store_true',
        help='This turns logging level to INFO'
    )

    parser.add_argument(
        '-vv',
        '--extra_verbose',
        action='store_true',
        help='This turns logging level to DEBUG'
    )

    parser.add_argument(
        '-n',
        '-nsavepoint',
        action='store',
        dest='nsavepoint',
        type=int,
        default=10,
        help='maximum number of EVIDs to process'
    )

    parser.add_argument(
        '-x',
        '--exclude_phaseless',
        action='store_true',
        help='Use this flag to exclude saving events that do not return phase data'
    )

    args = parser.parse_args()

    if args.extra_verbose:
        Logger.setLevel(logging.DEBUG)
    elif args.verbose:
        Logger.setLevel(logging.INFO)
    else:
        Logger.setLevel(logging.ERROR)

    main(args)