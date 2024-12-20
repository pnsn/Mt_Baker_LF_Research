import os, logging, glob
from pathlib import Path

import pandas as pd

from eqcutil.core.clusteringtribe import ClusteringTribe
from eqcutil.util.logging import setup_terminal_logger, CriticalExitHandler

Logger = setup_terminal_logger(os.path.split(__file__)[-1])
Logger.addHandler(CriticalExitHandler())
# Get absolute path of repo root directory
ROOT = Path(__file__).parent.parent.parent
# Get absolute path to pre-generated eventbank
EBBP = ROOT / "data" / "XML" / "QUAKE" / "BANK"
# Get absolute path to wavebank
WBBP = ROOT / "data" / "WF" / "BANK"
# Get RESP file directory
INVD = ROOT / 'data' / 'XML' / 'RESP'
# Output directory absolute path
DOUT = ROOT / "processed_data" / "workflow" / "templates"
# Status file save
STATUSFILE = DOUT / "template_construction_status.csv"
PFNAME = DOUT / "catalog_filtering_kwargs.csv"
# Parameter file to save
TCNAME = DOUT / "template_construction_kwargs.csv"


