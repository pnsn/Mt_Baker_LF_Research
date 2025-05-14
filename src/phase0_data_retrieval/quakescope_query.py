import datetime, time, logging, os
from pathlib import Path
import pandas as pd

logger = logging.getLogger('quakescope_query')

def main():


    rotator = pd.DataFrame()

    current_datetime = start_date
    next_datetime = start_date + datedelta
    while next_datetime <= end_date:
        tstr = current_datetime.strftime('%Y-%m-%d')
        estr = next_datetime.strftime('%Y-%m-%d')
        url = f'{base_url}?tid={sta_id}&start_time={tstr}&end_time={estr}&limit={qlimit:d}'
        
        try:
            if service == 'picks':
                payload = pd.read_csv(url, delimiter='|', parse_dates=['start_time','peak_time','end_time'])
            elif service == 'classifies':
                payload = pd.read_csv(url, delimiter='|', parse_dates=['start_time'])
            else:
                raise NotImplementedError(f'service "{service}" not supported')
        except Exception as e:
            print(f'Error fetching {service} data for {sta_id} | {tstr} - {estr} | {e}')
            current_datetime = next_datetime
            next_datetime += datedelta
            continue

    
        print(f'{sta_id} | {tstr} - {estr} | {len(payload)}')
        if payload.columns[0] == 'No data found!':
            logger.info(f'{sta_id} | {tstr} - {estr} | no data')
            current_datetime = next_datetime
            next_datetime += datedelta
            continue
        # If picks dont exceed query limit
        elif len(payload) < qlimit:
            current_datetime = next_datetime
            next_datetime += datedelta
        else:
            if service == 'picks':
                current_datetime = payload.end_time.max()
                next_datetime = payload.end_time.max() + datedelta
            elif service == 'classifies':
                current_datetime = payload.start_time.max()
                next_datetime = payload.start_time.max() + datedelta
        rotator = pd.concat([rotator, payload], ignore_index=True)
        rotator.drop_duplicates(keep='first', ignore_index=True, inplace=True)
        # # Update time window
        # current_datetime = rotator.end_time.max()
        # next_datetime = rotator.end_time.max() + datedelta
        
        if len(rotator) >= cut_limit:
            if service == 'picks':
                savename = SAVEPATH/f'{sta_id}_{rotator.start_time.min().strftime("%Y-%m-%d")}_{rotator.end_time.max().strftime("%Y-%m-%d")}.csv'
            elif service == 'classifies':
                savename = SAVEPATH/f'{sta_id}_{rotator.start_time.min().strftime("%Y-%m-%d")}_{rotator.start_time.max().strftime("%Y-%m-%d")}.csv'
            print(f'saving to {savename} ({len(rotator)} {service})')
            rotator.to_csv(savename, header=True, index=False)
            rotator = pd.DataFrame()
        
        time.sleep(0.1)

    # Do final save
    if len(rotator) > 0:
        if service == 'picks':
            savename = SAVEPATH/f'{sta_id}_{rotator.start_time.min().strftime("%Y-%m-%d")}_{rotator.end_time.max().strftime("%Y-%m-%d")}.csv'
        elif service == 'classifies':
            savename = SAVEPATH/f'{sta_id}_{rotator.start_time.min().strftime("%Y-%m-%d")}_{rotator.start_time.max().strftime("%Y-%m-%d")}.csv'
        print(f'saving to {savename} ({len(rotator)} picks)')
        rotator.to_csv(savename, header=True, index=False)



if __name__ == '__main__':
    logger = logging.getLogger()
    logger.setLevel(logging.DEBUG)
    # service = 'picks' # 'classifies'
    services = {'etype': 'classifies','pick':'picks'}
    service = services['etype']

    # Base query URL
    base_url = f"https://dasway.ess.washington.edu/quakescope/service/{service}/query"
    sta_id = 'UW.MBW2.'

    # Time Domain
    deep_start = pd.Timestamp('1980-01-01')
    cont_start = pd.Timestamp('2002-01-01')
    # start_date = pd.Timestamp('2014-10-26')
    start_date = cont_start
    end_date = pd.Timestamp('2025-04-01')


    # end_date = start_date
    # start_date = deep_start
    # Time step
    datedelta = pd.Timedelta(60, unit='days')

    # Batching limits
    qlimit = 10000
    cut_limit = 25000


    ROOT = Path(__file__).parent.parent.parent
    SAVEPATH = ROOT/'data'/'QuakeScope'/service.upper()/sta_id
    try:
        os.makedirs(str(SAVEPATH), exist_ok=False)
    except:
        pass

    main()

