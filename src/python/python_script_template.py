"""
:module: python_script_template.py
:author: Nathan T. Stevens
:email: ntsteven (at) uw.edu
:org: Pacific Northwest Seismic Network
:license: GNU GPLv3

.. rubric:: Purpose
    This module provides a coding style example that contributions to this repository should follow.
    It closely matches/adheres to the ObsPy Coding Style Guide (https://docs.obspy.org/coding_style.html).
    All document strings (docstrings) are written in reStructuredText format to allow automatic generation
    of code documentation using Sphinx (https://sphinx-rtd-tutorial.readthedocs.io/en/latest/).

    This docstring is the module header and it must contain the module, author, email, org, and license
    lines, as well as this Purpose section that provides a brief overview of what the module does.
"""


# In-line comment - I don't show up in auto-generated documentation
def im_a_method(arg, kwarg=None):
    """im_a_method docstring

    This block shows up in auto-generated documentation and also comes up
    as the method help documentation.

    :param arg: a positional argument and it is required. 
    :type arg: obj
    :param kwarg: a key-word argument, defaults to None.
    :type kwarg: obj, optional.
    :return: this is the single output of this method. It returns a tuple with the values (arg, kwarg)
    :rtype: tuple
    """    
    print('Python calls functions "methods", it took me some time to get used to...')
    return (arg, kwarg)

