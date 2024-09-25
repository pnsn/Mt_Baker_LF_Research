# To-Do/De-Bug/Fix-Me Items for this repository (Triggers TodoTree) for Benz Poobua

## Assigned by Nate Stevens
### Collaboration Tools Setup
TODO: Add "ToDo Tree" plugin to your VS Code installation  
TODO: Install the current version of `QGIS` on your local machine (if you haven't already)  
TODO: Install a conda environment using the `environment.yml`  

### Metadata Field Descriptions
https://ncedc.org/db/Documents/NewSchemas/PI/v1.6.4/PI.1.6.4/index.htm
Notes:
 - to_timestamp is actually origin.datetime in MtBaker_50km_radius_origins.csv
 - the first to_timestamp in MtBaker_20km_radius_phases.csv is the origin.datetime, the second to_timestamp is arrival.datetime


### Repository collaboration setup
TODO: Create a `benz-dev` branch on this repository
TODO: Clone the `benz-dev` branch to your local machine, remove this To-Do item, and push the change.

### Analysis Tasks: Catalog Exploration & Template Creation
 - Determine the frequency of different event classes at Mt. Baker.
 - Assess trends in Mt. Baker event source properties and solution qualities
 - Cross-reference the PNSN catalog with the REDPy catalog (see `README.md` for hyperlinks)
 - For each event class present at Mt. Baker, identify the event for each of these classes that has the most observations and make a script that isolates the phase picks for those events from `data/Events/MtBaker_20km_radius_phases.csv`.
 - Using the EQcorrscan API and example for generating templates from a `client`, create a template for each of the events identified in the previous task. Use P and S-wave pick times and fetch 3-C component data where available.

### Collect Contextual Data
TODO: Find the land-surface elevations for station locations in Table 5 from Frank et al. (1977) (Table_5_Frank_etal_1977_Mt_Baker...csv). This may be through publications, digital elevation models, checking in with Steve Malone, etc.  
TODO: Get a subset of the Randolph Glacier Inventory 6.0 for Washington State, save as a Shape File (actually generates several files) in `data/GIS/Polygons`. You will need to get a (free) Earthdata Login for the NASA hosted website/data. Let Nate know if you have any issues with this.  

