"""
:module: src/python/pull_wf_by_evid.py
:auth: Benz Poobua & Nate Stevens
:email: spoobu (at) uw.edu; ntsteven (at) uw.edu
:org: University of Washington; Pacific Northwest Seismic Network
:license: GNU GPLv3
:purpose: This module contains scripts that read and format CSVs output from AQMS
    and uses phase pick metadata to retrieve waveform picks from a web client.

"""
import os,sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np 
from obspy import read
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client
from obspy import Stream

sys.path.append(os.path.join('..', 'python'))

from util import unix_to_UTCDateTime
from well_located_events import *

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

def get_phase_entries_list(dataframe, evids):
    """Get a subset of phase pick rows from input **dataframe** that match the specified **evids**

    :param dataframe: The DataFrame containing phase pick data.
    :type dataframe: pandas.DataFrame
    :param evids: A list of event IDs to filter by.
    :type evids: list of int
    :return: A DataFrame containing rows that match the specified evids.
    :rtype: pandas.DataFrame
    """    
    if 'evid' not in dataframe.columns:
        raise ValueError("The DataFrame does not contain an 'evid' column.")
    
    if not isinstance(evids, list):
        raise ValueError("evids must be a list of integers.")

    filtered_entries = dataframe[dataframe.evid.isin(evids)]
    return filtered_entries

# define pick parameter used to select the interested station (e.g., station 'MBW)
# picks = phases_df_well_located[phases_df_well_located.sta == 'MBW'].iloc[:,:] if not phases_df_well_located.empty else None

def get_waveforms_from_phases(client, picks, filtered_entries, evids, front_pad_sec=10, back_pad_sec=30):
    """
    Fetch waveform data from IRIS (or other) web services for the provided picks
    and pad pick times by specified seconds before (front) and after (back).

    This function retrieves waveforms for multiple events based on their evidence IDs,
    applying necessary preprocessing steps such as detrending, filtering, and resampling.
    It collects all retrieved streams in a list and returns them.

    The function will not terminate if it encounters a No-Data-Available exception,
    and it handles errors gracefully by printing error messages.

    :param client: The IRIS web service client or equivalent service client.
    :type client: str
    :param picks: The specific pick entries to fetch waveforms for, usually in a DataFrame.
    :type picks: pandas.DataFrame
    :param filtered_entries: DataFrame containing filtered event entries with arrival times.
    :type filtered_entries: pandas.DataFrame
    :param evids: List of event IDs for which to retrieve waveforms.
    :type evids: list of int or str
    :param front_pad_sec: Seconds to pad before the pick time, defaults to 10 seconds.
    :type front_pad_sec: int, optional
    :param back_pad_sec: Seconds to pad after the pick time, defaults to 30 seconds.
    :type back_pad_sec: int, optional
    :return: A list of Stream objects containing the retrieved waveforms.
    :rtype: list of obspy.core.stream.Stream
    """

    if picks is None or picks.empty:
        print("No valid picks found for the specified station.")
        return

    # group by 'evid' and aggregate to find the earliest and latest arrival times
    earliest_latest_df = filtered_entries.groupby('evid').agg(
    earliest_arrival=('arrival_datetime', 'min'),
    latest_arrival=('arrival_datetime', 'max')
    ).reset_index()

    client = Client(client)
    streams = [] # collect all Stream objects
    for evid_ in evids:

        row = earliest_latest_df[earliest_latest_df['evid']==evid_]

        if not row.empty:
            earliest_arrival = row['earliest_arrival'].values[0]
            latest_arrival = row['latest_arrival'].values[0]

            starttime = earliest_arrival - pd.Timedelta(seconds=front_pad_sec)
            endtime = latest_arrival + pd.Timedelta(seconds=back_pad_sec)

            net = picks['net'].iloc[0]
            sta = picks['sta'].iloc[0]
            #chan = picks.seedchan
            bulk_ = [(str(net), str(sta), '*', '*', starttime, endtime)]

            try: 
                st = client.get_waveforms_bulk(bulk_)
                #st.remove_response(output='VEL')
                st.merge()
                st.detrend(type='linear')
                st.filter('bandpass',freqmin=0.1, freqmax=5)
                st.resample(10)
                st.taper(max_percentage=0.05)
                print(f"Evid: {evid_}, Retrieved Traces: {len(st)}")

                streams.append(st)

            except Exception as e:
                print(f'Error retrieving waveforms for evid {evid_}: {e}')
        
        else:
            print(f'No valid pick found for the specified evids')
    
    return streams
            
def get_waveforms(filename_events, filename_phases, evids_list, client, picks, savepath='../../results/waveforms'):
    """
    Composite function to retrieve waveforms based on provided event and phase pick data.

    This function loads event data from a CSV file, filters phase picks, and retrieves waveforms 
    from a web service (e.g., IRIS). It then saves the retrieved waveforms to specified files.

    :param filename_events: Name of the CSV file to load event data.
    :type filename_events: str
    :param filename_phases: Name of the CSV file to load phase picks.
    :type filename_phases: str
    :param evids_list: List of event IDs to filter phase picks.
    :type evids_list: list of int or str
    :param client: The IRIS (or other) web service client.
    :type client: str
    :param savepath: Directory to save the waveform files, defaults to '../../results/waveforms'.
    :type savepath: str, optional
    :return: None
    """
 
    # define the relative root directory path
    rrpath = os.path.join('..', '..')
    processed_path = os.path.join(rrpath, 'data', 'Processed Datasets')
    os.makedirs(processed_path, exist_ok=True)

    # ensure the output directory exists
    if not os.path.exists(savepath):
        os.makedirs(savepath)

    # load and curate the events
    input_file1 = filename_events
    cleaned_df = prep_data(input_file1)

    input_file2 = cleaned_df
    curated_df = curate_events(input_file2)

    curated_df_evid_list = curated_df['evid'].tolist()
    
    # load and filter phases
    phases_df = format_aqms_csv(filename_phases)
    phases_df_well_located = get_phase_entries_list(phases_df, curated_df_evid_list)

    # get waveforms
    for evid in evids_list:
        st_list = get_waveforms_from_phases(client, picks, phases_df_well_located, evids_list) 

        if st_list:  # if the list is not empty
            for st in st_list:  # iterate through each stream returned
                print(st)
                out_file_name = f'uw.{evid}.mseed'
                save_name = os.path.join(savepath, out_file_name)
                print(f"Attempting to save to: {save_name}") 

                try:
                    st.write(save_name, format='MSEED')
                    print(f'Waveforms saved to {save_name}')
                except Exception as e:
                    print(f"Error saving waveforms: {e}")
                
                st.plot();  
        else:
            print(f'No waveforms to save for evid {evid}.')




