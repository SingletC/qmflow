"""Prepare data for Plotly Dash."""
from datetime import datetime

import ase.db.core
import numpy as np
import pandas as pd
from ase.db.table import Table


def create_dataframe(db: ase.db.core.Database, selection=None, columns=None):
    """Create Pandas DataFrame from database """
    if columns is None:
        columns = ['formula', 'name', 'osc_str', 'lambda_', 'ctime', 'reaction', 'id', 'nto_type', 'note','delta_G',
                   'delta_G_TS','bn_index']
    query = list(db.select(selection=selection))
    table = []
    for i in query:
        row = []
        for column in columns:
            row += [getattr(i, column, 0)]

        row += [f'''<a href="details?smiles={(i.get("name", "")).replace('#', '$')}" target="_blank">
        <img src="data:image/png;;base64, {i.data.get("img", "")}">
        </a>''']
        table += [row]
    columns.append('Structure')
    df = pd.DataFrame(table, columns=columns)
    df['lifetime(ns)'] = round(
        1.5 / ((1e7 / df['lambda_']) ** 2 * df['osc_str']) * 1e9  # ns
        , 2)
    df['ctime'] = df['ctime'].apply(
        lambda x: datetime.fromtimestamp(x * ase.db.core.YEAR + ase.db.core.T2000))  # .strftime('%c'))
    return df
