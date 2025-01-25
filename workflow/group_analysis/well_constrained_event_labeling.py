import logging, pathlib, os

import pandas as pd
from sklearn.metrics import normalized_mutual_info_score

from eqcutil import ClusteringTribe, basic_logger_config

basic_logger_config(level=logging.INFO)
Logger = logging.getLogger(os.path.split(__file__)[-1])


ROOT = pathlib.Path(__file__).parent.parent.parent
Logger.warning(f'running from {ROOT}')

PD_DIR = ROOT / "processed_data" / "workflow" 
PREF_FILE = PD_DIR / "catalog" / "preferred_event_sta_picks.csv"
