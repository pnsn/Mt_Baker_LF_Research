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
RDIR = ROOT / 'processed_data' / 'cluster' / 'single_station'

WRITE_PROTECT = True

ccckwargs = {'method': 'xcc',
            'replace_nan_distances_with': 'mean',
            'shift_len': 10,
            'corr_thresh': 0.45,
            'allow_individual_trace_shifts': True,
            'cores': 'all'}


def main():
    ctr_dict = {}
    try:
        os.makedirs(str(RDIR), exist_ok=False)
    except:
        pass
    dir_list = glob.glob(str(TDIR/'*'))
    for _dir in dir_list:
        # Ignore CSV last complete evid files
        if _dir[-3:] == 'csv':
            continue
        Logger.info(f'Reading templates from {_dir}')
        # NSLC_DIR/eval_mode/year/uwEVID.tgz
        ipath = Path(_dir) / '*' / '*' / '*.tgz'
        f_list = glob.glob(str(ipath))

        _, nslc = os.path.split(_dir)
        SAVEFILE = RDIR / f'{nslc}_clustered.tgz'
        
        if WRITE_PROTECT:
            if os.path.isfile(str(SAVEFILE)):
                Logger.warning(f'{SAVEFILE} already exists - skipping')
                continue


        ctr = ClusteringTribe()
        if len(f_list) > 1:
            pass
        else:
            Logger.warning(f'Only 1 tempate for {_dir}')
            continue
        Logger.info(f'...assembling {len(f_list)} templates...')
        for _f in f_list:
            tmp = Tribe().read(_f)[0]
            if tmp.process_length is None:
                tmp.process_length = 90.
            if tmp.prepick is None:
                tmp.prepick = 5.
            ctr += tmp
        Logger.info('...clustering...')
        ctr.cluster(**ccckwargs)
        Logger.info('...complete...')

        Logger.info('...assembling metadata...')
        ctr.populate_event_metadata()
        
        ctr_dict.update({nslc: ctr})
        Logger.info("...writing to disk.")
        ctr.write(str(SAVEFILE))
        
    return ctr_dict
        # fn = os.path.split(_tmp)
        # sn, ext = os.path.splitext(fn)
        # Logger.info('...overwriting unclustered file...')
        # ctr.write(_tmp)






if __name__ == '__main__':
    # Setup Logging
    Logger = setup_terminal_logger(os.path.split(__file__)[-1], level=logging.INFO)
    Logger.addHandler(CriticalExitHandler(exit_code=1))
    # RUN MAIN
    main()   