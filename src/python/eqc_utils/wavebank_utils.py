import os, glob
from pathlib import Path
root = Path(__file__).parent.parent.parent.parent
from obspy import read
from obsplus import WaveBank

bank_base_path = os.path.join(root,'data','waveforms','BANK')

def initialize_wavebank(mseed_files=[],
                        base_path=bank_base_path,
                        path_structure='{year}',
                        name_structure='{seedid}.{time}',
                        **kwargs):
    """TODO: Fill out docstring

    :param mseed_files: _description_, defaults to []
    :type mseed_files: list, optional
    :param base_path: _description_, defaults to bank_base_path
    :type base_path: _type_, optional
    :param path_structure: _description_, defaults to '{year}'
    :type path_structure: str, optional
    :param name_structure: _description_, defaults to '{seedid}.{time}'
    :type name_structure: str, optional
    :return: _description_
    :rtype: _type_
    """
    if not os.path.exists(bank_base_path):
        os.makedirs(bank_base_path)

    wbank = WaveBank(base_path=base_path,
                     path_structure=path_structure,
                     name_structure=name_structure,
                     **kwargs)
    
    for msfile in mseed_files:
        st = read(msfile)
        wbank.put_waveforms(st)
    return wbank


def connect_to_wavebank(base_path=bank_base_path,
                        path_structure='{year}',
                        name_structure='{seedid}.{time}',
                        **kwargs):
    """TODO: Fill out docstring

    :param base_path: _description_, defaults to bank_base_path
    :type base_path: _type_, optional
    :param path_structure: _description_, defaults to '{year}'
    :type path_structure: str, optional
    :param name_structure: _description_, defaults to '{seedid}.{time}'
    :type name_structure: str, optional
    :return: _description_
    :rtype: _type_
    """    
    wbank = WaveBank(base_path=base_path,
                     path_structure=path_structure,
                     name_structure=name_structure,
                     **kwargs)
    return wbank
