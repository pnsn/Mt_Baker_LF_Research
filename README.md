# Mt_Baker_LF_Research: (Re)analysis of Low Frequency Seismic Sources at Mt. Baker Volcano
This repository hosts PNSN analysis codes and metadata for analysis of low-frequency seismic events at Mt. Baker volcano in Washington state, USA.

## Background (in brief)
Although Mt. Baker produces relatively few volcano-tectonic earthquakes, it produces a significant amount of low-frequency seismicity (LFs) from magmatic and glaciologic processes.
Glacier and volcano sourced LFs from Mt. Baker have similar waveform characteristics, typically requiring accurate event location solutions to differentiate LF source type.
From 1972 to 2009 Mt. Baker was continuously monitored by only one seismic station (UW.MBW) within 30 km of the summit, in the past 15 years 6 additional permanent seismic stations have been added
and UW.MBW was upgraded (site UW.MBW2). The current sub-network around Mt. Baker now allows detection of smaller seismic events and improved location accuracy, providing an opportunity to better constrain
differences between magmatic and glaciologic LF sources and use these insights to review and enhance the PNSN LF catalog at Mt. Baker, and potentially other glacier-clad volcanoes.

Also see the USGS synopsis on [Mt. Baker seismic monitoring](https://www.usgs.gov/volcanoes/mount-baker/science/earthquake-monitoring-mount-baker).

## Research Questions
1) What are the source processes associated with LFs at Mt. Baker?
2) Can we identify waveform data features that reliably differentiate LF source types?  
3) Can we enhance the PNSN catalog of LFs at Mt. Baker using new, single-station analytic methods?  
4) Can we gain new insights on Mt. Baker's magmatic (and glaciologic) processes?

## License  
![image](./docs/Figures/gplv3-with-text-136x68.png)  
Original works contained in this repository are distributed under the attached GNU General Public License v3.

## Repository Structure  
 - **data** - unprocessed or minimally processed data files that **should be treated as readonly**
    - Events - Data files for seismic catalog metadata  
    - Sensors - Data files for geophysical instrument metadata  
 - **docs** - supporting files for repository documentation  
 - **GIS** - Geographic Information Systems files  
    - Mt_Baker_LF_Maps.qgz - QGIS project file  
 - **processed_data** - intermediate analysis data files **all of which will be .gitignore'd for the online repository**
 - **results** - final analysis output files
    - figures - Rendered figure files  
    - tables - Final tabular data files
 - **src** - source code  
    - notebooks - Jupyter Notebook files  
    - PostgreSQL - SQL database query scripts  
    - python - python scripts
 - **todo** - holds todo lists for collaborators  
 - environment.yml - conda environment definition for scripts in `src/python`
 - LICENSE - distribution terms  
 - README.md - you are here!  

## Environment Installation
Source code living in the `src/python` and `src/notebooks` will likely make up the bulk of this repository. You can create a `conda` environment for required (and anticipated) dependencies
using the `environment.yml` file included in the root directory of this repository. Once you have `conda` installed, you can create the `baker_lf` environment using:  
```conda env create -f environment.yml```

## Collaboration
### How to Collaborate on this Repository
Collaborative work on this repository should roughly follow the [ObsPy contribution guidelines](https://docs.obspy.org/contributing.html). Detailed
guidance on general GitHub collaboration can be found [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests).
In a nutshell:  
1) Fork this repository
2) Make a new branch
3) Work on something
4) Push to your fork
5) [Submit a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request)
6) Wait for repo manager review of pull request

### People
Nathan T. Stevens (PNSN Seismologist/Developer) - Repository Manager / Research Mentor (ntsteven@uw.edu)  
Renate Hartog (PNSN Network Manager) - PNSN GitHub Organization Owner / Research Supervisor  
Benz Poobua (ESS Undergraduate Researcher) - Research Mentee  
Steve Malone (PNSN Director Emeritus)  
Alex Hutko (PNSN Research Seismologist)  
Amy Wright (PNSN Lead Seismic Analyst)  
Barrett Johnson (PNSN Seismic Analyst)  
Wes Thelen (CVO Research Geophysicist)  
...likely more to come   

### Organizations  
 - [PNSN Scientific Products Team](https://pnsn.org)  
 - [UW ESS Department Researchers](https://ess.uw.edu)    
 - [USGS Cascade Volcano Observatory Researchers](https://www.usgs.gov/observatories/cvo)   

### Collaborator Documents Repository  
A repository of written documents and references for permissioned collaborators is available on the PNSN GoogleDrive.  
