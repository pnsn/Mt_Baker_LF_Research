#!/bin/bash
: '
As the name impiles, runs everything in the current processing pipeline.

Starting point is having an up-to-date EventBank in `data/XML/QUAKE/BANK`

proceeds with:
 src/phase1_catalog_curation scripts starting with `step*`
 src/phase2_template_construction scripts starting with `step*`
 src/pahse3_grouping_analysis scripts starting with `step*`

'


python ./phase2_template_construction/step1_generate_single_station_templates.py
python ./phase2_template_construction/step2_cluster_single_station_tribes.py
