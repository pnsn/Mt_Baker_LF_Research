# DLP_Mt_Baker: Distinguishing Long-Period Seismic Sources at Mt. Baker Volcano
This repository hosts code and metadata for analyzing long-period seismic events at Mt. Baker volcano in Washington state, USA.

## Motivation  
Mt. Baker has been seismically monitored since 1972, but has produced relatively few observable volcano-tectonic earthquakes, even during its period of unrest in 1975 (Crider et al., 2011). However, Mt. Baker produces Deep Long Period seismicity (DLP)--which are defined as earthquakes with source solutions deeper than 10 km b.s.l. and dominant periods in the 0.2 - 1 second band--with a relative abundance, accounting for 31 of 60 DLPs catalogged between 1980 and 2009 (Nichols et al., 2011). These events are associated with magma and/or volatile migration at depth and have been characterized at many of the Cascade volcanoes and other volcanoes globally (Nichols et al., 2011; and references therein). Subsequent monitoring by the PNSN adds 21 further DLPs to the Baker catalog between 2010 and 2024 (Figure 1). This percieved rise in DLP occurrence may be an observational artifact, likely arising from a combination of improvements to instrumentation, increased station density around Mt. Baker (3 stations added in 2009, 2019, and 2023), the nature of PNSN automated detection algorithms, and revised analyst practice motivated by findings from temporary deployments on Mt. Baker in 1975-76 (Malone, 1977; Weaver & Malone, 1979) and 2009 (Caplan-Auerbach et al., 2009) and detailed observation at other Cascades volcanoes (e.g., Nichols et al., 2011).

FIGURE 1 NEAR HERE (time-series of DLPs at Mt. Baker from PNSN catalog, with DLPs between 1980 and 2009 in blue, and those post 2009 in forest green).

The glaciers cladding Mt. Baker (and other Cascades volcanoes) also produce long period seismic events in the same band band as DLPs, and are primarily distinguised by source constraints near land-surface (Weaver and Malone, 1979; Caplain-Auerbach et al., 2009; Thelen et al., 2013; Allstadt & Malone, 2014). The relatively sparse seismic network around Mt. Baker results in only the largest LP events are sufficient to generate an event trigger in the PNSN's automated detection pipeline, meaning that smaller events are likely missed, or may be indistinguishable from glacier-sourced LPs (as discussed in Thelen et al., 2013). Moreover, surface events have only recently been retained in the PNSN catalog as standard practice, meaning that smaller DLPs may have been missed due to prior analyst policies.

Advances in seismic analysis tools and computing power give us the opportunity to rapidly re-analyze large volumes of data that was a substantial hurdle to earlier studies (e.g., Mousavi & Beroza, 2023). These include efficient template matching workflows (e.g,. Hotovec-Ellis and Jefferies, 2016; Chamberlain et al., 2017) and machine-learning enhanced models specifically trained for DLP earthquakes (Münchmeyer et al., 2024) and descriminating glacier-sourced events from earhquakes (Kharita et al., 2024). In this study we will conduct a survey of all LP events localized to Mt. Baker in the PNSN catalog to inform application of one or more of these methods on the entirety of digitally available waveform data for stations within a 30 km radius of Mt. Baker. Through this analysis we hope to address the following questions:

## Research Questions  
1) Does the PNSN catalog represent a comprehensive set of DLPs for Mt. Baker?  
2) Is the percieved uptick in cataloged DLPs in the past decade an observational artifact, as hypothesized?  
3) Can we identify dependable time-series-derived characteristics to differentiate DLP and glacier-sourced LP events at Mt. Baker?  

## Research Stages
1) Catalog Survey: review the [PNSN](https://pnsn.org/events?custom_search=true) and [REDpy](https://assets.pnsn.org/red/) catalogs to:
    - verify accurate classifications of 'LF' (for DLPs) and 'SU' (for shallow LPs)  
    - assess if DLPs and glacier LPs are triggering clusters in the repeating earthquake catalog.  

2) Catalog Classification Refinement: refine classifications of long period events using features identified in the literature and new features that arise from observations in step 1. Refined classes might be:  
    - deep long-period (DLP)  
    - shallow long-period (SLP) or glacier long-period (GLP)
    - under-constrained long-period (LP)   

3) Catalog Enhancement: enhance the PNSN long-period seismicity catalog at Mt. Baker using observations from steps 1-2 and one or more of the following analytic approaches:   
    - Waveform cross correlation / template matching using [EQcorrscan](https://eqcorrscan.readthedocs.io/en/latest/) (after Thelen et al., 2013). 
    - Detection/phase picking with the low-frequency earthquake detection model [LFEDetect](https://seisbench.readthedocs.io/en/stable/pages/documentation/models.html#seisbench.models.lfe_detect.LFEDetect) distributed with [SeisBench](https://seisbench.readthedocs.io).
    - Feature-driven detection/classification (e.g., Kharita et al., 2024). 

4) Interpretation: review the enhanced catalog to determine if meaningful patterns exist that better inform our understanding of magmatic processes at Mt. Baker.

## License  
This open-source repository is licensed under the attached GNUv3 license.

## Originating Author  
Nathan T. Stevens (PNSN Seismologist/Developer) - Research Mentor  

## Project Collaborators
Benz Poobua (ESS Undergraduate Researcher) - Research Mentee 
Renate Hartog (PNSN Network Manager) - Research Supervisor
Steve Malone (PNSN Director Emeritus)  
Alex Hutko (PNSN Research Seismologist)

## Collaborating Organizations  
 - Scientific Products Team, [Pacific Northwest Seismic Network](https://pnsn.org)  
 - ESS Researchers, [University of Washington Department of Earth and Space Sciences](https://ess.uw.edu)    
 - Research Staff, [USGS Cascade Volcano Observatory](https://www.usgs.gov/observatories/cvo)


## Works Cited
 - [Allstadt, K., and Malone, S.D. (2014) Swarms of repeating stick-slip ice-quakes triggered by snow loading at Mt. Rainier volcano. JGR-Earth Surface, 119, 1180-1203](https://doi.org/10.1002/2014JF003086).

 - [Caplan-Auerbach, J., Thelen, W.A., and Moran, S.C. (2009) An Unusual Cluster of Low-Frequency Earthquakes at Mount Baker, Washington, as Detected by a Local Broadband Network. EOS Trans. AGU 89, Fall Meeting Suppl., Paper number V23D-2111](https://ui.adsabs.harvard.edu/abs/2009AGUFM.V23D2111C/abstract)

 - [Chamberlain, C.J., Hopp, C.J., Boese, C.M., Warren-Smith, E., Chambers, D., Chu, S.X., Michailos, K., and Townsend, J. (2017) EQcorrscan: Repeating and Near-Repeating Earthquake Detection and Analysis in Python. SRL 89(1): 173-181](https://doi.org/10.1785/0220170151)

 - [Crider, J.G., Johnsen, K.H., and Williams-Jones, G. (2008) Thirty-year gravity change at Mount Baker volcano, Washington, USA: extracting the signal from under the ice. Geophys Res Lett 35:L20304– L20308](https://doi.org/10.1029/2008GL034921)

 - [Hotovec-Ellis, A.J., and Jeffries, C. (2016) Near Real-time Detection, Clustering, and Analysis of Repeating Earthquakes: Application to Mount St. Helens and Redoubt Volcanoes – Invited, presented at Seismological Society of America Annual Meeting, Reno, Nevada, 20 Apr.](https://code.usgs.gov/vsc/REDPy)

 - [Kharita, A., Denolle, M.A., West, M.E. (2024) Discrimination between icequakes and earthquakes in southern Alaska: an exploration of waveform features using Random Forest algorithm. GJI 237(2), 1189-11207](https://doi.org/10.1093/gji/ggae106)

 - [Malone, S.D. (1977) Summary of Seismicity and Gravity, pg. 19-25 in Frank, D., Meier, M.F., and Swanson, D., Assessmentof Increased THermal Activity at Mount Baker, Washington March 1975-March 1976. U.S. Geological Survey Professional Paper 1022-A](https://pubs.usgs.gov/pp/1022a/report.pdf)

 - [Mousavi, S.M., and Beroza, G.C. (2023) Machine Learning in Earthquake Seismology. Ann. Rev. Earth and Planet. Sci. 51: 105-129](https://doi.org/10.1145/annurev-earth-071822-100323)

 - [Münchmeyer, J., Giffard-Roisin, S., Malfante, M., Frank, W.B., Poli, P., Marsan, D., and Socquet, A. (2024) Deep learning detects uncataloged low-frequency earthquakes across regions. Seismica , 3(1)](https://doi.org/10.26443/seismica.v3i1.1185)

 - [Nichols, M.L., Malone, S.D., Moran, S.C., Thelen, W.A., and Vidale, J.E. (2011) Deep long-period earthquakes beneath Washington and Oregon volcanoes. J. Volc. and Geotherm. Res. 200(3–4), 116–128.](https://doi.org/10.1016/j.jvolgeores.2010.12.005)

 - [Thelen, W.A., Allstadt, K., De Angelis, S., Malone, S.D., Moran, S.C., and Vidale, J. (2013) Shallow repeating seismic events under an alpine glacier at Mount Rainier, Washington, USA. J. Glac., 59(214), 345-356.](https://doi.org/10.3189/2013Jog12J111)

 - [Weaver, C.S., and Malone, S.D. (1979) Seismic evidence for discrete glacier motion at the rock-ice interface. J. Glaciol 23:171–184](https://doi.org/10.3189/S0022143000029816)

