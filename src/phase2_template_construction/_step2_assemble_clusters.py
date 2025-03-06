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

WRITE_PROTECT = False
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
    # Get channel names from directory names
    chan_dirs = glob.glob(str(TDIR/'*'))
    chanids = [_c.split('/')[-1] for _c in chan_dirs]
    chanids.sort()
    # Iterate across channel IDs
    processed = set()
    for chanid in chanids:
        mask_loc = False
        mask_bi = False
        n,s,l,c = chanid.split('.')
        # Group multiple location codes
        if len(glob.glob(str(TDIR/f'{n}.{s}.*.{c}'))) > 1:
            chanid = f'{n}.{s}.*.{c}'
            name = f'{n}.{s}.XX.{c}'
            mask_loc = True
        # Group multiple sensors/sampling bands
        elif len(glob.glob(str(TDIR/f'{n}.{s}.{l}.*'))) > 1:
            chanid = f'{n}.{s}.{l}.??{c[-1]}'
            name = f'{n}.{s}.{l}.XX{c[-1]}'
            mask_bi = True
        else:
            name = chanid
        if name in processed:
            continue
        # ACTUALLY get file list
        flist = glob.glob(str(TDIR/chanid/'*'/'*'/'*.tgz'))
        # Make savename
        SAVEFP = RDIR / f'{name}.tgz'
        if SAVEFP.is_file() and WRITE_PROTECT:
            continue
        
        ctr = ClusteringTribe()
        for _e, _f in enumerate(flist):
            # Update every 10
            if _e %10 == 0 and _e > 0:
                Logger.info(f'{name} - ({_e+1}/{len(flist)})')
            if mask_loc:
                Logger.debug(f'{name} - ALIASED LOCATIONS')
            if mask_bi:
                Logger.debug(f'{name} - ALIASED BAND AND INSTRUMENT CODES')
            Logger.debug(f'Loading {_f}')
            template = Tribe().read(_f)[0]
            if mask_loc or mask_bi:
                for tr in template.st:
                    if mask_loc:
                        Logger.debug(f'aliasing {tr.id} location')
                        tr.stats.location = 'XX'
                    if mask_bi:
                        Logger.debug(f'aliasing {tr.id} band and instrument codes')
                        tr.stats.channel = f'XX{tr.stats.component}'
            ctr += template
        Logger.info(f'RUNNING CLUSTERING FOR {name}')
        ctr.cluster(**ccckwargs)
        Logger.info(f'POPULATING EVENT METADATA')
        ctr.populate_event_metadata()
        Logger.info(f'WRITING TO DISK')
        ctr.write(str(RDIR/name))
        processed.add(name)



if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
    Logger.addHandler(CriticalExitHandler(exit_code=1))
    # RUN MAIN
    main()   