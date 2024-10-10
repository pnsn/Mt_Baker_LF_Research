"""
:module: src/python/pull_wf_by_evid.py
:auth: Benz Poobua & Nate Stevens
:email: spoobu (at) uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: This module is used to fetch events from AQMS and fetch subset waveform pickings from IRIS web services.

"""
import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from obspy import read
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

sys.path.append(os.path.join('..', 'python'))

from util import unix_to_UTCDateTime

def format_aqms_csv(filename, colmapping={'datetime': 'origin_datetime', 'datetime.1': 'arrival_datetime'}):
    """This method reads a raw output CSV file from an AQMS file, renames columns 
    as specified by **colmapping**, and formats POSIX datetimes to UTCDateTimes.

    :param filename: name of the CSV file to load
    :type filename: str
    :param colmapping: column name relabeling mapping, defaults to {'datetime': 'origin_datetime', 'datetime.1': 'arrival_datetime'}.
    :type colmapping: dict, optional

    :return:
     - **df** (*pandas.dataframe.DataFrame*) -- loaded, formatted
    """    
    # relative path
    rrpath = os.path.join('..', '..')
    phases = os.path.join(rrpath, 'data', 'Events', filename)
    phases_df = pd.read_csv(phases)

    # rename the columns
    phases_df.rename(columns=colmapping, inplace=True)

    # Delete rows that contain NaN in these two columns
    phases_df = phases_df.dropna(subset=['origin_datetime','arrival_datetime'])

    # convert float to int
    phases_df.origin_datetime = phases_df.origin_datetime.astype(int)
    phases_df.arrival_datetime = phases_df.arrival_datetime.astype(int)

    # employ unix_to_UTCDateTime function to convert epoch times to UTC datetime
    phases_df.origin_datetime = phases_df.origin_datetime.apply(lambda x: unix_to_UTCDateTime(x))
    phases_df.arrival_datetime = phases_df.arrival_datetime.apply(lambda x: unix_to_UTCDateTime(x))

    return phases_df

def get_phase_entries(dataframe, evid):
    """Get a subset of phase pick rows from input **dataframe** that match the specified **evid**

    :param dataframe: The DataFrame containing phase pick data.
    :type dataframe: pandas.DataFrame
    :param evid: The event ID to filter by.
    :type evid: int
    :return: A DataFrame containing rows that match the specified evid.
    :rtype: pandas.DataFrame
    """    
    if 'evid' not in dataframe.columns:
        raise ValueError("The DataFrame does not contain an 'evid' column.")

    filtered_entries = dataframe[dataframe.evid == evid]
    return filtered_entries

def get_waveforms_from_phases(client, pick, front_pad_sec=30, back_pad_sec=60):
    """Fetch waveform data from IRIS web services for the provided pick
    and pad pick times by specified seconds before (front) & after (back).

    The function will not terminate if it encounters a No-Data-Available exception.

    :param client: The IRIS web service client.
    :type client: str
    :param pick: The specific pick entry to fetch waveforms for.
    :type pick: pandas.Series
    :param front_pad_sec: Seconds to pad before the pick, defaults to 30.
    :type front_pad_sec: int, optional
    :param back_pad_sec: Seconds to pad after the pick, defaults to 60.
    :type back_pad_sec: int, optional
    """

    # request waveform data from IRIS web service 
    if pick is not None:
        client_ = Client(client)
        time_ = UTCDateTime(pick['arrival_datetime'])
        starttime = time_ - front_pad_sec
        endtime = time_ + back_pad_sec
        net = pick.net
        sta = pick.sta
        loc = pick.location
        chan = pick.seedchan

        try:
            Z = client_.get_waveforms(net,sta,loc,chan, starttime, endtime, attach_response=True)
            Z.remove_response(output='VEL')
            Z.merge()
            Z.detrend(type='linear')
            Z.taper(max_percentage=0.05)
            #Z.filter('bandpass',freqmin=0.5, freqmax=25)
            Z.plot();
            return Z
        
        except Exception as e:
            print(f"Error fetching waveforms: {e}")
            return None

    else: 
        print("No valid pick found for the specified evid.")
        return None

def get_waveforms(filename, evid, client, pick, savepath='../../results/waveforms'):
    """Composite function to get waveforms based on a CSV file and event ID.

    :param filename: Name of the CSV file to load phase picks.
    :type filename: str
    :param evid: Event ID to filter phase picks.
    :type evid: int or str
    :param client: The IRIS web service client.
    :type client: str
    :param savepath: Directory to save the waveform file, defaults to '../../results/waveforms'.
    :type savepath: str, optional
    :return: The Stream object containing waveforms.
    :rtype: obspy.Stream
    """
    
    # ensure the output directory exists
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    # load and filter phases
    phases_df = format_aqms_csv(filename)
    filtered_events = get_phase_entries(phases_df, evid)

    # get waveforms
    st = get_waveforms_from_phases(client, pick)

    # construct the output filename
    if st is not None and len(st) > 0:
        out_file_name = f'{st[0].id}.{evid}.mseed'
        save_name = os.path.join(savepath, out_file_name)

        # save the waveforms
        st.write(save_name, format='MSEED')
        print(f'Waveforms saved to {save_name}')
    else:
        print('No waveforms to save.')

    return st

