import os, warnings
import pandas as pd
from pathlib import Path

ROOT = Path(__file__).parent.parent
SAVE = ROOT/'data'/'MUSTANG'

try:
    os.makedirs(str(SAVE), exist_ok=False)
except:
    pass

base_url = 'https://service.iris.edu/'
service = 'mustang/measurements/1/'
interface = 'query?'
metrics = ['percent_availability','sample_min','max_range','sample_unique','num_spikes','num_gaps']
network='UW'
stations=['MBW','MBW2','JCW','SHUK','RPW','RPW2']
channel='??Z'
format='xml'

startdate = f'2002-01-01'
enddate = f'2026-01-01'

general_url = f'{base_url}{service}{interface}network={network}&channel={channel}&format={format}'
for station in stations:
    for metric in metrics:

        print(f'Fetching {metric} for {network}.{station} | {startdate} - {enddate}')
        url = general_url + f'&metric={metric}&station={station}&timewindow={startdate},{enddate}&nodata=404'
        try:
            df = pd.read_xml(url)
            try:
                os.makedirs(str(SAVE/f'{network}.{station}'), exist_ok=False)
            except:
                pass
            df = df.rename(columns={'value':metric})
            df.to_csv(SAVE/f'{network}.{station}'/f'{metric}_{startdate}_{enddate}.csv', header=True, index=False)
        except Exception as e:
            warnings.warn(f'{network}.{station} | {metric} | {startdate} - {enddate} | {e}')
            continue
