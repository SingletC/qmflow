"""Prepare data for Plotly Dash."""
from datetime import datetime

import ase.db.core
import numpy as np
import pandas as pd
from ase.db.table import Table


def create_dataframe(db: ase.db.core.Database, selection=None, columns=None):
    """Create Pandas DataFrame from database """
    if columns is None:
        columns = ['formula', 'name', 'osc_str','lambda_', 'ctime']
    query = list(db.select(selection=selection))
    table = []
    for i in query:
        row = []
        for column in columns:
            try:
                row+=[getattr(i,column )]
            except AttributeError:
                row+=[0]
        row +=[f'<img src="data:image/png;;base64, {i.data.get("img","")}">']
        table +=[row]
    columns.append('Structure')
    df = pd.DataFrame(table,columns=columns)

    df['ctime'] = df['ctime'].apply(lambda x: datetime.fromtimestamp(x * ase.db.core.YEAR + ase.db.core.T2000).strftime('%c'))
    return df
