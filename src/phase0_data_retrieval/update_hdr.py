import os
from glob import glob
from pathlib import Path

import pandas as pd
from obspy import read


def main():
    ROOT = Path(__file__).parent.parent.parent
    WWU = ROOT/'data'/'WWU'
    NET = 'WW'
    LOC = ''
    CBI = 'HH'

    for YEAR in ['2007','2008','2009']:
        RAW_STA_DIRS = glob(str(WWU/'RAW'/YEAR/'*'))
        breakpoint()
        for STA in [os.path.split(_e)[-1] for _e in RAW_STA_DIRS]:
            RAWZ = WWU/'RAW'/YEAR/STA/'*z*'/f'{YEAR}*.sac'
            RAWN = WWU/'RAW'/YEAR/STA/'*n*'/f'{YEAR}*.sac'
            RAWE = WWU/'RAW'/YEAR/STA/'*e*'/f'{YEAR}*.sac'
            DGPS = WWU/'DATA'/YEAR/STA/'GPS'/'parsed_clock_data.csv'
            DATA = WWU/'DATA'/YEAR/STA/'SAC'

            df = pd.read_csv(DGPS, index_col=[0], parse_dates=['utcdatetime'])
            stla = df.lat.median()
            stlo = df.lon.median()
            stel = df.elev.median()

            try:
                os.makedirs(str(DATA), exist_ok=False)
            except:
                pass

            for _R, _C in [(RAWZ,'Z'),(RAWN,'N'),(RAWE,'E')]:
                # Get list of channel-specific files
                flist = glob(str(_R))
                flist.sort()
                # Load files
                for _f in flist:
                    try:
                        st = read(_f)
                    except:
                        continue
                    if len(st) != 1:
                        breakpoint()
                    else:
                        tr = st[0]
                    # Update NSLC
                    tr.stats.network = NET
                    tr.stats.station = STA
                    tr.stats.location = LOC
                    tr.stats.channel = CBI+_C
                    tr.stats.sac.update({'stla': stla, 'stlo': stlo, 'stle': stel})
                    savename = f"{NET}.{STA}.{LOC}.{CBI}{_C}_{tr.stats.starttime.strftime('%Y%m%d_%H%M%S')}.sac"
                    print(f'Writing: {DATA/savename}')
                    tr.write(str(DATA/savename), format='SAC')            


if __name__ == '__main__':
    main()