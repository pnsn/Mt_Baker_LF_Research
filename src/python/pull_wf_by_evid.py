"""
:module: src/python/pull_wf_by_evid.py
:auth: Benz Poobua & Nate Stevens
:email: spoobu (at) uw.edu
:org: University of Washington
:license: GNU GPLv3
:purpose: 

"""

def format_aqms_csv(filename, colmapping={'datetime': 'origin_datetime', 'datetime.1': 'arrival_datetime'}):
    """This method reads a raw output CSV file from an AQMS file, renames columns 
    as specified by **colmapping**, and formats POSIX datetimes to UTCDateTimes and necessary
    corrections to make boolean entries.

    :param filename: name of the CSV file to load
    :type filename: str
    :param colmapping: column name relabeling mapping, defaults to {'datetime': 'origin_datetime', 'datetime.1': 'arrival_datetime'}.
    :type colmapping: dict, optional

    :return:
     - **df** (*pandas.dataframe.DataFrame*) -- loaded, formatted
    """    

def get_phase_entries(dataframe, evid):
    """Get a subset of phase pick rows from input **dataframe** that match the specified **evid**

    :param dataframe: _description_
    :type dataframe: _type_
    :param evid: _description_
    :type evid: _type_
    """    

def get_waveforms_from_phases(dataframe, client, front_pad_sec=30, back_pad_sec=60):
    """Fetch waveform data from IRIS web services for picks provided in **dataframe**
    and pad pick times by specified seconds before (front) & after (back).

    Make sure the function dosen't terminate if it hits a No-Data-Available type exception.

    :param dataframe: _description_
    :type dataframe: _type_
    :param client: _description_
    :type client: _type_
    :param front_pad_sec: _description_, defaults to 30
    :type front_pad_sec: int, optional
    :param back_pad_sec: _description_, defaults to 60
    :type back_pad_sec: int, optional
    """


def get_waveforms(filename, evid, savepath='.'):
    """Composite function of the above.

    :param filename: _description_
    :type filename: _type_
    :param evid: _description_
    :type evid: _type_
    :param savepath: _description_, defaults to '.'
    :type savepath: str, optional
    :return: _description_
    :rtype: _type_
    """



    out_file_name = f'uw{evid}.mseed'
    save_name = os.path.join(savepath, out_file_name)
    st.write(out_file_name)
    return st