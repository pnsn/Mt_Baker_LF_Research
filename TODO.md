# General To-Do/De-Bug/Fix-Me Items for this repository (Triggers TodoTree)

## Tasks for Nate
### Initial Data
FIXME: Radius search is not returning the same population of events the search box grabbed. This needs to be reconciled.
FIXME: FOR README - Remove the use of DLP in the motivation to the PNSN standard classification of LF, and clarify the 10 km threshold comes from Nichols et al. (2011) and we use it here for comparison purposes.

### Collaboration Tools / Setup  
TODO: add note on command-line compile needed for `conda env create -f environment.yml` to README
### Code Examples  
TODO: Create template generation example from Stevens SGGS SSA repository  
TODO: Create template matching exaple from Stevens SGGS SSA repository  

## Tasks for Benz  
### Collaboration Tools Setup
TODO: Add "ToDo Tree" plugin to your VS Code installation  
TODO: Install the current version of `QGIS` on your local machine (if you haven't already)  
TODO: Install a conda environment using the `environment.yml`  
TODO: Pull a copy of this repository to your local machine, remove this To-Do item, and push the change. Then let Nate know and standby.  

### Collect Contextual Data
TODO: Find the land-surface elevations for station locations in Table 5 from Frank et al. (1977) (Table_5_Frank_etal_1977_Mt_Baker...csv). This may be through publications, digital elevation models, checking in with Steve Malone, etc.  
TODO: Get a subset of the Randolph Glacier Inventory 6.0 for Washington State, save as a Shape File (actually generates several files) in `data/GIS/Polygons`. You will need to get a (free) Earthdata Login for the NASA hosted website/data. Let Nate know if you have any issues with this.  
TODO: Get a CSV that has station information for every seismometer within a 30 km radius of Mt. Baker. Suggest using the IRIS Metadata Aggregator or an ObsPy client (method `obspy.core.client.Client.get_stations`).  

### Make Figures for the Project Overview Document  
TODO: Using data files in this repository and descriptions in the [Project Overview Document](https://docs.google.com/document/d/1d9cfkK4y1T6hd_p6u7PEp3NhYLGydLFzR-DxrLNWt7c/edit?usp=share_link), create first drafts of Figures 1A, 1B, 2A, and 2B using QGIS and `matplotlib`