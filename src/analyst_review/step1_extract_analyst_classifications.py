import os
from random import sample
from pathlib import Path
import pandas as pd
# from sklearn.preprocessing import OneHotEncoder
# from sklearn.compose import ColumnTransformer
# from sklearn.metrics import silhouette_score, adjusted_mutual_info_score, normalized_mutual_info_score

# Absolute path for repository root directory
ROOT = Path(__file__).parent.parent.parent
# Absolute path to survey responses
DATA = ROOT / 'data' / 'Survey' / 'Classification_Responses.csv'
# Absolute path to AQMS etypes | uwEVID### pairs
AQMS = ROOT / 'data' / 'Events' / 'Mount_Baker_evids_etypes_10_JAN_2025.csv'
# Save file name
SAVENAME = ROOT / 'processed_data' / 'survey' / 'S1_extracted_reviewer_classes.csv'
# Load Raw Survey Responses
df = pd.read_csv(DATA)
# Pivot to place questions as index
df = df.T
# Split into header lines, optional response lines, and classification response lines
hdr = []
opt_lines = []
class_lines = []
for idx, row in df.iterrows():
    if idx =='REVIEWER NUMBER':
        hdr.append(row)
    elif 'Optional' in idx:
        opt_lines.append(row)
    else:
        class_lines.append(row)

# Parse anonymized Event Classifications
df_class = pd.DataFrame()
for _r in class_lines:
    df_class = pd.concat([df_class, _r], axis=1, ignore_index=False)

# Place reviewer numbers as columns
df_class = df_class.T
# Set index as full uwEVID#### 
df_class.index = [f"uw{_ind.split('[')[-1][:-1]}" for _ind in df_class.index]
df_class.index.name = 'EVID'

# Generate random sequence for reviewer shuffle
new_order = sample(range(4),4)
df_class = df_class.rename(columns={_k: new_order[_k] for _k in range(len(df_class.columns))})
# # Update reviewer numbers to "REVIEWER#" format
df_class = df_class.rename(columns={_k: f'REVIEWER{_k}' for _k in df_class.columns})

# Load AQMS etype/evid pairs
df_aqms = pd.read_csv(AQMS, index_col='uw_evid')
# Split down to series
ser_aqms = df_aqms.etype
# rename series
ser_aqms.name = 'AQMS'
# Leftjoin AQMS classifications to surveys
df_class = df_class.join(ser_aqms, how='left')

# pivot to access column names for sorting
df_class = df_class.T.sort_index().T



if not SAVENAME.parent.is_dir():
    os.makedirs(str(SAVENAME.parent))
df_class.to_csv(SAVENAME, index=True, header=True)


# df_cfac = df_class.copy()
# # Factorize Columns
# for col in df_class.columns:
#     df_cfac[col], _ = pd.factorize(df_cfac[col])

# # breakpoint()
# # Conduct pairwise score on reviewer pairs
# scores = []
# for _i, _ci in enumerate(df_class.columns):
#     for _j, _cj in enumerate(df_class.columns):
#         if _i < _j:
#             idf = df_cfac[(df_class[_ci].notna()) & (df_class[_cj].notna())]
#             print(len(idf))
            

#             line = [
#                 _ci,
#                 _cj,
#                 len(idf),
#                 adjusted_mutual_info_score(idf[_ci].values, idf[_cj].values)
#             ]

#             scores.append(line)


# breakpoint()
# rnotes = 
# # Iterate across lines
# for oline in opt_lines:
#     # Iterate across reviewer responses
#     for rno, row in oline.items():
#         # If non NaN row
#         if isinstance(row, str):
#             semi_split = row.split(';')
#             for ele in semi_split:
#                 nlsplit = ele.split('\n')
#             breakpoint()
        # if all(row.isna()):

