import os
from pathlib import Path
from glob import glob
import pandas as pd

def parse_gpgga_line(line):
    parts = line.split(',')
    # Message ID Header
    mid = parts[0]
    if mid != '$GPGGA':
        raise KeyError
    tstr = parts[1]
    # # UTC Observation Timestamp
    # hh = parts[1][:2]
    # mm = parts[1][2:4]
    # ss = parts[1][4:]
    # utctime=pd.Timestamp(f'{datestr}T{hh}:{mm}:{ss}')
    # LATITUDE
    lat = float(parts[2][:2]) + float(parts[2][2:])/60.

    if parts[3] == 'N':
        lat = lat
    elif parts[3] == 'S':
        lat = -1.*lat
    else:
        raise ValueError
    # LONGITUDE
    lon = float(parts[4][:3]) + float(parts[4][3:])/60.
    if parts[5] == 'W':
        lon = -1.*lon
    elif parts[5] == 'E':
        lon = lon
    else:
        raise ValueError
    # Position Fix Indicator
    pfind = int(parts[6])
    # Satellites Used
    nsat = int(parts[7])
    # HDOP Horizontal Dilution of Precision
    hdop = float(parts[8])
    # MLS ALT - Mean Sea Level Altitude
    elev = float(parts[9])
    # ELEV units 
    eunit = parts[10]
    if eunit == 'M':
        pass
    else:
        raise NotImplementedError(f'non-(M)eters united ({eunit}) elevations not supported')
    # Geoid separation
    if len(parts[11]) > 0:
        gsep = int(parts[11])
    else:
        gsep = None
    # Geoid separation unit
    gsu = parts[12]
    # age of differential
    if len(parts) == 16:
        age = int(parts[13])
        dcorr = int(parts[14])
        cksum = parts[15]
    elif len(parts) == 15:
        age = None
        dcorr = None
        cksum = parts[14]

    oline = [tstr, lon, lat, elev, pfind, hdop, nsat, gsep, gsu, age, dcorr]
    return oline

def parse_gpzda_line(line):
    """
    Parse a GPZDA line that transfers UTC Time and Date information
    """
    parts = line.split(',')
    mid = parts[0]
    hh = parts[1][:2]
    mm = parts[1][2:4]
    ssff = parts[1][4:]
    tstr = parts[1]
    DD = parts[2]
    MM = parts[3]
    YYYY = parts[4]
    tz = parts[5]
    if len(parts) == 8:
        tzm = parts[6]
        cksum = parts[7]
    else:
        tzm = None
        cksum = parts[6]
    return [YYYY, MM, DD, tstr, hh, mm, ssff, tz, tzm, cksum]
    # return [pd.Timestamp(f'{YYYY}-{MM}-{DD}T{hh}:{mm}:{ssff}'), tz, tzm, cksum]

if __name__ == '__main__':

    ROOT = Path(__file__).parent.parent
    for STA in ['BLDR','COUG','CRIM','PAUL','WELL']:
        # STA = 'RAIL'
        YEAR = '2009'
        RAW = ROOT/'RAW'/YEAR/STA/'*00'
        OUT = ROOT/'DATA'/YEAR/STA/'GPS'

        flist = glob(str(RAW/f'{YEAR}*.txt'))
        flist.sort()
        holder = []
        index = []
        for _f in flist:
            print(f'reading: {_f}')
            try:
                with open(_f, 'r') as f:
                    for line in f:
                        if '$GPZDA' in line:
                            _gpzdaline = parse_gpzda_line(line)
                        elif '$GPGGA' in line:
                            try:
                                _gpggaline = parse_gpgga_line(line)
                            except:
                                continue
                            if _gpggaline[0] == _gpzdaline[3]:
                                utctime = pd.Timestamp(f'{_gpzdaline[0]}-{_gpzdaline[1]}-{_gpzdaline[2]}T{_gpzdaline[4]}:{_gpzdaline[5]}:{_gpzdaline[6]}')
                                line = _gpggaline[1:]
                                holder.append(line)
                                index.append(utctime)
            except:
                continue
        
        df = pd.DataFrame(
            holder, 
            index=index,
            columns=['lon','lat','elev','pfind','hdop','nsat','gsep','gsu','age','dcorr'])
        df.index.name = 'utcdatetime'
        df = df.sort_index()
        try:
            os.makedirs(OUT, exist_ok=False)
        except:
            pass
        df.to_csv(OUT/'parsed_clock_data.csv', header=True, index=True)
