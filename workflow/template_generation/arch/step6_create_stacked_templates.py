# Imports
import logging, glob, os

from pathlib import Path

from eqcutil.core.clusteringtribe import ClusteringTribe

# Map absolute paths
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Processed Data Directory
PD_DIR = ROOT / "processed_data" / "workflow"
# Template Save Directory
CTSD = PD_DIR / 'templates' / 'single_station' / 'xcc_test'

cct = 0.9


# Get file names for templates
flist = glob.glob(str(CTSD/'*.tgz'))

# Iterate across stations
for file in flist:
    ctr = ClusteringTribe().read(file)
    # "Peel the onion" - apply high CC thresholds to get at sub-groups
