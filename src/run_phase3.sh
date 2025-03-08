#!/bin/bash
: '
As the name impiles, runs everything in the current processing pipeline.

Starting point is having an up-to-date EventBank in `data/XML/QUAKE/BANK`

proceeds with:
 src/phase1_catalog_curation scripts starting with `step*`
 src/phase2_template_construction scripts starting with `step*`
 src/pahse3_grouping_analysis scripts starting with `step*`

'


python ./phase3_grouping_analysis/step1_assemble_distances.py
python ./phase3_grouping_analysis/step2_assess_ensemble_grouping.py
python ./phase3_grouping_analysis/step3_assess_single_station_grouping.py
