from pathlib import Path
import pandas as pd


ROOT = Path(__file__).parent.parent.parent
DATA = ROOT / 'data' / 'Survey' / 'Classification_Responses.csv'

# Load Raw Responses
df = pd.read_csv(DATA)
# Pivot to place questions as index
df = df.T
hdr = []
opt_lines = []
class_lines = []
for idx, row in df.iterrows():
    if idx =='REVIEWER NUMBER']:
        hdr.append(row)
    elif 'Optional' in idx:
        opt_lines.append(row)
    else:
        class_lines.append(row)

# Parse anonymized Event Classifications
df_class = pd.DataFrame()
for _r in class_lines:
    df_class = pd.concat([df_class, _r], axis=1, ignore_index=False)

df_class = df_class.T
df_class.index = [f"uw{_ind.split('[')[-1][:-1]}" for _ind in df_class.index]

df_opt = pd.DataFrame()
