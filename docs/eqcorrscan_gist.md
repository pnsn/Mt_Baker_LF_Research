# EQcorrscan Gist
EQcorrscan (Chamberlain et al., 2017) is a popular seismic waveform template matching toolbox that provides
code implementation of methods from key template matching (also called matched filter) methods published
in the scienfic literature.

Template matching fundamentally uses the cross-correlation of template and continuous waveform time-series.
The method produces a detection function where increasing commonality between the two time-series' features result in
larger positive values in the detection function.

The general work-flow of a template matching analysis is:  
1) Candidate Template Generation  
    a) Determine a subset of already-characterized seismic events that have a reasonable likelihood of coming from periodically excitable sources.  
    b) Acquire an accurately picked catalog of phase arrival data for candidate template events and unprocessed waveform data that adequately encompass wavetrains of interest for each event.
2) Template Selection
    a) Which stations should be included?
    b) What is the minimum amount of signal processing needed to accentuate desired waveform features?
