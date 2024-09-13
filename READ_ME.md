# DLP_Mt_Baker: Distinguishing Long-Period Seismic Sources at Mt. Baker Volcano
This repository hosts code and metadata for analyzing seismicity at Mt. Baker volcano in Washington state, USA,
focusing on detection and source-characterization of long-period earthquakes.

## Motivation  
Mt. Baker has been seismically monitored since 1972 by station UW.MBW and produces relatively few volcano-tectonic earthquakes, even during its period of unrest in 1975 (Crider et al., 2011; and references therein). However, it is the focus of most of the Deep Long Period seismicity (earthquakes with source solutions deeper than 10 km b.s.l. and dominant periods in the 0.2 - 1 second band) cataloged in the Cascades, accounting for 31 of 60 DLPs detected between 1980 and 2009 for all of Washington and Oregon (Nichols et al., 2011). Subsequent monitoring by the PNSN adds 21 further DLPs to the Baker catalog between 2010 and 2024. This percieved uptick in DLP occurrence may be an observational artifact, likely arising from a combination of improvements to instrumentation, increased station density around Mt. Baker (3 stations added in 2009, 2019, and 2023), and revised analyst practice motivated by findings from temporary deployments on Mt. Baker in 1975-76 (Weaver & Malone, 1976) and 2009 (Caplan-Auerbach et al., 2009) and detailed observation at other Cascades volcanoes (e.g., Nichols et al., 2011).

The glaciers cladding Mt. Baker (and other Cascades volcanoes) also produce long period seismic events in the same band band as DLPs, but are located near the land-surface (Weaver and Malone, 1979; Caplain-Auerbach et al., 2009; Thelen et al., 2013). In combination with the relatively sparse seismic network around Mt. Baker, only the largest LP events are sufficient to generate an event trigger in the PNSN's automated detection pipeline, meaning that smaller events are likely missed.

## Research Goals
1) Review the PNSN catalog to verify accurate classifications of 'LF' (for DLPs) and 'SU' (for shallow LPs).

2) Refine classifications of long period events using features identified in the literature and new features that arise from observations in step 1. Refined classes might be:  
  - deep long-period (DLP)  
  - shallow long-period (SLP)  
  - under-constrained long-period (LP)   

3) Refine/enhance the PNSN long-period seismicity catalog at Mt. Baker using observations from steps 1-2 and one or more of the following analytic approaches
  - Waveform cross correlation / template matching (e.g., Thelen et al., 2010; Allstadt & Malone, 2014)
  - Feature-driven detection/classification (e.g., Kharita et al., 2024)
  - Detection/phase picking with the low-frequency earthquake detection model [LFEDetect](https://seisbench.readthedocs.io/en/stable/pages/documentation/models.html#seisbench.models.lfe_detect.LFEDetect) distributed with [SeisBench](https://seisbench.readthedocs.io) 

4) Review the enhanced catalog to determine if meaningful patterns exist that better inform our understanding of magmatic processes at Mt. Baker.

## License  
This open-source repository is licensed under the attached GNUv3 license.

## Originating Author  
Nathan T. Stevens (PNSN)  
 

## Contributors  
 - [Pacific Northwest Seismic Network](https://pnsn.org) Scientific Products Team  
 - [University of Washington Department of Earth and Space Sciences](https://ess.uw.edu) Researchers  
 - [USGS Cascade Volcano Observatory](https://www.usgs.gov/observatories/cvo) Staff  
