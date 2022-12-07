"""Prepare data for Plotly Dash."""
from datetime import datetime

import ase.db.core
import numpy as np
import pandas as pd
from ase.db.table import Table


def create_dataframe(db: ase.db.core.Database, selection=None, columns=None):
    """Create Pandas DataFrame from database """
    if columns is None:
        columns = ['user', 'formula', 'osc_str','lambda_', 'ctime']
    query = list(db.select(selection=selection))
    table = []
    for i in query:
        row = []
        for column in columns:
            row+=[getattr(i,column )]
        row +=[f'<img src="data:image/png;;base64, {i.data.img}">']
        table +=[row]
    columns.append('Structure')
    df = pd.DataFrame(table,columns=columns)
    ini_time_for_now = datetime.now()
    df['ctime'] = df['ctime'].apply(lambda x: (-pd.Timedelta(x,"h").to_pytimedelta()+ini_time_for_now).strftime("%m/%d/%Y, %H:%M"))
    return df
