from psycopg2 import connect
import pandas as pd

def dbconnect(dbname='mountbaker',user='postgres',port='5252'):
    conn = connect(f'dbname={dbname} user={user} port={port}')
    return conn

def sbpick_insert(cursor, row, model_name, weight_name, weight_version=1, leapcorr=True)
    mapping = ''
    values = {}
    for _k, _v in row.items():
        if _k in ['start_time','peak_time','end_time']:
            if isinstance(_v, pd.Timestamp):
                values.update({_k: _v.timestamp()})
            elif isinstance(_v, (int, float)):
                values.update({_k: _v})
            elif _v == pd.NaT:
                continue
            else:
                raise TypeError
        

        mapping += f'%({_k})s'
            