
"""
:module: BANJO-DETECT/src/python/utils/template_utils.py
:auth: Nathan T. Stevens & Barrett Johnson
:email: ntsteven@uw.edu; bnjo@uw.edu
:org: Pacific Northwest Seismic Network
:license: GPLv3
:purpose: This module contains helper scripts for EQcorrscan Template
    management in addition to those native to EQcorrscan

:attribution: Based on scripts from EQcorrscan documentation and 
    example scripts by B. Johnson.
"""
import copy

def rename_templates(tribe, inplace=True):
    """Rename templates in a tribe with their PNSN EVID's
    with the option to enact the changes on a deepcopy of th e
    :param tribe: tribe to relabel
    :type tribe: :class:`~eqcorrscan.core.match_filer.Tribe
    :return:
     - **tribe** (*eqcorrscan.core.match_filter.Tribe*) - tribe of renamed templates
    """
    # If not inplace, make a deepcopy of the tribe within the
    # scope of this method
    if not inplace:
        tribe = copy.deepcopy(tribe)
    for template in tribe:
        event = template.event
        name = str(event.resource_id).split(sep='/')[-1]
        template.name = name
    return tribe