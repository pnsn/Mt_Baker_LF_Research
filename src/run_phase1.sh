#!/bin/bash
: '
As the name impiles, runs everything in the current processing pipeline.

Starting point is having an up-to-date EventBank in `data/XML/QUAKE/BANK`

proceeds with:
 src/phase1_catalog_curation scripts starting with `step*`
 src/phase2_template_construction scripts starting with `step*`
 src/pahse3_grouping_analysis scripts starting with `step*`

'

python ./phase1_catalog_curation/step1_select_subcatalog_events.py
python ./phase1_catalog_curation/step2_select_preferred_stations.py
python ./phase1_catalog_curation/step3_prep_template_metadata.py

