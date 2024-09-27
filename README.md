# Mt_Baker_LF_Research
## PNSN (Re)analysis of Low Frequency Seismic Sources at Mt. Baker Volcano, Washington, USA
**This repository hosts analysis codes and metadata for research on low-frequency seismic events at Mt. Baker volcano, Washington, USA.**

## Motivation

Mt. Baker does not produce many volcano-tectonic earthquakes, but it does produce a significant number of low-frequency seismic events (LFs) associated with either magmatic or glaciologic processes. Glacier and volcano sourced LFs from Mt. Baker have similar waveform characteristics, typically requiring accurate event location solutions to differentiate LF source type (Malone, 1977; Weaver & Malone, 1979; Caplan-Auerbach et al., 2009; Crider et al., 2011; Nichols et al., 2011; Thelen et al., 2013).


![image](./docs/Figures/Mt_Baker_Catalog_and_Stations_Timeseries_100dpi.png)  
***Figure 1** (left axis) Time-series of annual seismic event frequency at Mt. Baker (within 10 km of the summit) in the PNSN catalog. Events are split by event type: LF = low frequency, EQ = earthquake, SU = surface event, PX = probable blast (key). (right axis) Time-series of PNSN-operated seismic station counts within 30 km of Mt. Baker, sampled monthly. Station metadata from the IRIS MetaData Aggregator.*


From 1972 to 2001 Mt. Baker was continuously monitored by only one seismic station (UW.MBW) - a vertical component, analog station - within 30 km of the summit, making it difficult or impossible to locate lower magnitude events. Since 2001 seven new permanent stations were added to this region, and in 2023 UW.MBW2 was installed 2 km from UW.MBW while UW.MBW continued to operate until December 2023 (Fig. 2). Seismic monitoring of Mt. Baker has been periodically augmented by temporary seismic deployments on the volcano itself, including one from 1975 to 1977 (Malone, 1977; Weaver & Malone, 1979) and three seasonal deployments during 2007, 2008, and 2009 (Caplan-Auerbach et al., 2009). These temporary arrays enabled improved LF detection, location, and characteization, revealing that a substantial portion LF seismicity at Mt. Baker arises from glacier motion (Weaver & Malone, 1979; Caplan-Auerbach et al., 2009). The current permanent sub-network around Mt. Baker also provides similar benefits, resulting in a substantial increase in the occurrence of LFs in the post-2009 PNSN catalog (Fig. 1). Due to these observational biases, it is difficult to extract meaningful trends from the current catalog of LFs at Mt. Baker.


![image](./docs/Figures/Overview_Maps_120dpi.png)  
***Figure 2** Regional (left) and local (right) maps of PNSN stations within 40 km of Mt.Baker and PNSN catalog seismic events within 10 km of Mt.Baker. (key) Seismic station marker colors indicate the station's installation year, whereas seismic event epicenter colors indicate their classification. Station label colors indicate if the station is currently active (white) or decomissioned (orange).The domain of the local map is shown as a black outline in the regional map.Figure generated using QGIS 3.32, basemap imagery from Mapzen Global Terrain (both) and OpenStreetMap (local map).*


If we can identify waveform-based features that distinguish LF source types using these well-characterized events, we can use computationally fast, single-station analytic methods (e.g., Chamberlain et al., 2017; Münchmeyer et al., 2024) to reanalyze the entire continuous waveform archive for stations around Mt. Baker to enhance the PNSN LF catalog, address observational biases in the catalog, and potentially gain new insights on magmatic and glaciologic processes at Mt. Baker.

Also see the USGS summary on [Mt. Baker seismic monitoring](https://www.usgs.gov/volcanoes/mount-baker/science/earthquake-monitoring-mount-baker).

## Research Questions
1) What are the source processes associated with LFs at Mt. Baker?
2) Can we identify waveform data features that reliably differentiate LF source types?  
3) Can we enhance the PNSN catalog of LFs at Mt. Baker using new, single-station analytic methods?  
4) Can we gain new insights on Mt. Baker's magmatic (and glaciologic) processes?

## License  
![image](./docs/Figures/gplv3-with-text-136x68.png)  
Original works contained in this repository are distributed under the attached GNU General Public License v3.

## Repository Structure  
 - **data** - unprocessed or minimally processed data files **that should be treated as readonly files**
    - Events - PNSN seismic catalog metadata  
    - Sensors - geophysical instrument metadata  
 - **docs** - supporting files for repository documentation  
 - **GIS** - Geographic Information Systems files  
    - Mt_Baker_LF_Maps.qgz - QGIS project file  
 - **processed_data** - intermediate analysis data files **all of which is .gitignore'd in this GitHub repository**
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

## Environment & Relevant Software Versions
### Conda Environment for Python Scripts
Source code in the `src/python` and `src/notebooks` directories will likely make up the bulk of this repository. 
To create a `conda` environment `baker_lf`, use the included `environment.yml` file as such:  

```conda env create -f environment.yml```

### QGIS
The QGIS project in this repository was created using QGIS version 3.32 "Lima". QGIS can be downloaded [here](https://www.qgis.org).

### PostgreSQL + AQMS
SQL queries in the `src/PostgreSQL` directory were composed for PostgreSQL version 14.X databases used by AQMS (see Hartog et al., 2020). The schema for AQMS databases follows the CISN/ANSS Parametric Information Schema documented [here](https://ncedc.org/db/Documents/NewSchemas/PI/v1.6.4/PI.1.6.4/index.htm).

## Collaboration
### People
Nathan T. Stevens (PNSN Seismologist/Developer) - Repository Manager / Research Mentor (ntsteven@uw.edu)  
Renate Hartog (PNSN Network Manager) - PNSN GitHub Organization Owner / Research Supervisor  
Benz Poobua (ESS Undergraduate Researcher) - Research Mentee   

### Organizations  
 - [PNSN Scientific Products Team](https://pnsn.org)  
 - [UW Earth and Space Sciences Department](https://ess.uw.edu)    
 - [USGS Cascade Volcano Observatory](https://www.usgs.gov/observatories/cvo)   

### Document Repository  
A repository titled `Mt Baker LF Research` for written documents, presentations, and reference materials is available on the PNSN GoogleDrive account for permissioned collaborators.  

### Coding Style Guide
The ObsPy development team provides a great [Coding Style Guide](https://docs.obspy.org/coding_style.html) in their documentation. Please adhere to their guidelines.

### How to Collaborate on this Repository
Collaborative work on this repository follows the [Fork and Pull](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/about-collaborative-development-models#fork-and-pull-model) development model. This model is widely used for open-source software development: e.g., ObsPy, see their contribution guide [here](https://docs.obspy.org/contributing.html)).

Detailed guidance on collaborative development using GitHub is provided [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests).

Collaboration Steps In a Nutshell:  
1) Fork this repository
2) Make a new branch
3) Work on something locally
4) Push changes to your forked repository
5) [Submit a pull request to the main repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) to the primary repository.
6) Wait for review of pull request by an admin.

Best practices for pull requests are detailed [here](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/getting-started/best-practices-for-pull-requests).


## References

- [Caplan-Auerbach, J., Thelen, W.A., and Moran, S.C. (2009) An Unusual Cluster of Low-Frequency Earthquakes at Mount Baker, Washington, as Detected by a Local Broadband Network. EOS Trans. AGU 89, Fall Meeting Suppl., Paper number V23D-2111](https://ui.adsabs.harvard.edu/abs/2009AGUFM.V23D2111C/abstract)

- [Chamberlain, C.J., Hopp, C.J., Boese, C.M., Warren-Smith, E., Chambers, D., Chu, S.X., Michailos, K., and Townsend, J. (2017) EQcorrscan: Repeating and Near-Repeating Earthquake Detection and Analysis in Python. SRL 89(1): 173-181](https://doi.org/10.1785/0220170151)

- [Crider, J.G., Frank, D., Malone, S.D., Poland, M.P., werner, C., and Caplan-Auerbach, J. (2011). Magma at depth: a retrospective analysis of the 1975 unrest at Mt. Baker, Washington, USA. Bull. Volcanol., 73, 175-189.](https://doi.org/10.1007/s00445-010-0441-0)

- [Hartog, J.R., Friberg, P.A., Kress, V.C., Bodin, P., and Bhadha, R. (2020) Open-Source ANSS Quake Monitoring System Software. Seismol. Res. Lett. 91(2A). 677-686.](https://doi.org/10.1785/0220190219)

- [Malone, S.D. (1977) Summary of Seismicity and Gravity, pg. 19-25 in Frank, D., Meier, M.F., and Swanson, D., Assessmentof Increased THermal Activity at Mount Baker, Washington March 1975-March 1976. U.S. Geological Survey Professional Paper 1022-A](https://pubs.usgs.gov/pp/1022a/report.pdf)

- [Münchmeyer, J., Giffard-Roisin, S., Malfante, M., Frank, W.B., Poli, P., Marsan, D., and Socquet, A. (2024) Deep learning detects uncataloged low-frequency earthquakes across regions. Seismica , 3(1)](https://doi.org/10.26443/seismica.v3i1.1185)

- [Nichols, M.L., Malone, S.D., Moran, S.C., Thelen, W.A., and Vidale, J.E. (2011) Deep long-period earthquakes beneath Washington and Oregon volcanoes. J. Volc. and Geotherm. Res. 200(3–4), 116–128.](https://doi.org/10.1016/j.jvolgeores.2010.12.005)

- [Thelen, W.A., Allstadt, K., De Angelis, S., Malone, S.D., Moran, S.C., and Vidale, J. (2013) Shallow repeating seismic events under an alpine glacier at Mount Rainier, Washington, USA. J. Glac., 59(214), 345-356.](https://doi.org/10.3189/2013Jog12J111)

- [Weaver, C.S., and Malone, S.D. (1979) Seismic evidence for discrete glacier motion at the rock-ice interface. J. Glaciol 23:171–184](https://doi.org/10.3189/S0022143000029816)