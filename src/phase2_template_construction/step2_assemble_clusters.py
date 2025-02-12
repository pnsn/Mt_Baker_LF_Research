import logging, os, glob
from pathlib import Path

from eqcorrscan import Tribe
from eqcutil import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

# Absolute path to repository root
ROOT = Path(__file__).parent.parent.parent
# Absolute path to templates
TDIR = ROOT / 'processed_data' / 'template' / 'single_station'
# Absolute path to results subdirectory
RDIR = ROOT / 'results' / 'cluster' / 'single_station'

ALIAS_INSTRUMENTS = True
ccckwargs = {'method': 'xcc',
            'replace_nan_distances_with': 'mean',
            'shift_len': 10,
            'corr_thresh': 0.5,
            'allow_individual_trace_shifts': True,
            'cores': 'all'}


def main():
    # Ensure main save directory exists
    if not RDIR.is_dir():
        os.makedirs(str(RDIR))
    # Get channel names from 
    chan_dirs = glob.glob(str(TDIR/'*'))
    chanids = [_c.split('/')[-1] for _c in chan_dirs]
    chanids.sort()

    processed_chanids = set()

    for chanid in chanids:
        # Skip location 01's that will 
        if chanid.split('.')[2] == '01':
            Logger.warning(f'Skipping {chanid} - has 01 location code that will be merged with -- location code')
            continue
        Logger.info(f'Processing {chanid}')
        n,s,l,c = chanid.split('.')
        clist = glob.glob(str(TDIR/f'{n}.{s}.*.*'))
        # Check if there are Band / Instrument codes to alias
        if len(clist) > 1 and ALIAS_INSTRUMENTS:
            bi = 'XX'
            Logger.warning('Mixed instruments found, will alias band & instrument codes on traces to XX')
            chanid = f'{n}.{s}.{l}.{bi}{c[-1]}'
        else:
            bi = False
        if chanid in processed_chanids:
            Logger.warning(f'CHANID {chanid} already processed - skipping')
            continue
        tlist = glob.glob(str(TDIR/f'{n}.{s}.*.*'/'*'/'*'/'*.tgz'))
        ctr = ClusteringTribe()
        for _e, _t in enumerate(tlist):
            Logger.debug(f'Loading {_t}')
            # Update every 10
            if _e %10 == 0 and _e > 0:
                Logger.info(f'({_e+1}/{len(tlist)})')
            # Load Template
        #     template = Tribe().read(_t)[0]
        #     # If band/instrument are aliased
        #     if bi:
        #         # "Iterate" across trace
        #         for tr in template.st:
        #             Logger.warning(f'Alias {tr.stats.channel} --> {bi}{tr.stats.component}')
        #             tr.stats.channel = f'{bi}{tr.stats.component}'
        #     # add to clustering tribe
        #     ctr += template
        # # Apply alias to clustering tribe save name

        # # Run clustering
        # Logger.warning(f'Running clustering for {chanid}')
        # ctr.cluster(**ccckwargs)
        # # Populate metadata
        # ctr.populate_event_metadata()
        # # Write to disk
        # Logger.info(f'Writing to disk')
        # SAVENAME = RDIR / chanid
        # ctr.write(str(SAVENAME))
        # # Mark chanid as processed
        processed_chanids.add(chanid)

    breakpoint()


if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
    Logger.addHandler(CriticalExitHandler(exit_code=1))
    # Set up to-file logger
    # CLOG =

    main()   